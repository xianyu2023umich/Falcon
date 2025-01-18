Module ModCommunication

    use ModBlock,       only:   BlockType,GC_target
    use ModYinYang,     only:   ModYinYang_VecConv_1D
    use ModYinYangTree, only:   YYTree
    use ModGC,          only:   ModGC_FindGC_Single,&
                                ModGC_GetGC_Targets_single,&
                                ModGC_SetGC_Sources_Single
    use ModMath,        only:   ModMath_1D3D_interpolate_1D1D
    use ModControl,     only:   nMaxBlocksPerRank
    use ModParameters,  only:   MpiSize,MpiRank,&
                                ni,nj,nk,ng,nvar
    use MPI

    contains

    ! use Mpi_reduce to determine the global time step
    !
    subroutine ModCommunication_GlobalTimeStep(dt_local,dt_global)
        implicit none
        real,intent(in)                 ::  dt_local
        real,intent(out)                ::  dt_global
        integer                         ::  ierr
        call MPI_AllReduce(dt_local,dt_global,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ierr)
    end subroutine ModCommunication_GlobalTimeStep

    subroutine ModCommunication_SetGCAll(Tree)
        implicit none
        type(YYTree),target             ::  Tree                    ! the Tree
        type(BlockType),pointer         ::  Block1                  ! one Block
        integer                         ::  iBlock

        do iBlock=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iBlock)
            call ModGC_FindGC_Single(Tree,Block1)
        end do

        if(MpiRank==0) write(*,*)'Completed finding GC...'

        do iBlock=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iBlock)
            call ModGC_GetGC_Targets_single(Tree,Block1)
        end do
        
        if(MpiRank==0) write(*,*)'Completed setting GC_targets...'

        call ModCommunication_GetnGC_sources(Tree)

        if(MpiRank==0) write(*,*)'Completed counting nGCsource...'

        do iBlock=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iBlock)
            call ModGC_SetGC_Sources_Single(Tree,Block1)
        end do

        if(MpiRank==0) write(*,*)'Completed setting GC_sources...'
    end subroutine ModCommunication_SetGCAll

    subroutine ModCommunication_GetnGC_sources(Tree)
        implicit none
        type(YYTree),intent(in),target  ::  Tree
        type(BlockType),pointer         ::  Block1,Block2
        type(GC_target),pointer         ::  GC_target1,GC_source1
        integer                         ::  nGC_sources_table_local(Tree%NumLeafNodes)
        integer                         ::  nGC_sources_table_global(Tree%NumLeafNodes)
        integer                         ::  iLocalBlock,iGC_target
        integer                         ::  iRecv,nRecv,recv_message(3)
        integer                         ::  isend,nSend
        integer                         ::  ierr,request,status(MPI_STATUS_SIZE)
        
        ! Get the local table for GC_sources of each block
        nGC_sources_table_local=0
        do iLocalBlock=1,Tree%nLocalBlocks
            Block1=>Tree%LocalBlocks(iLocalBlock)
            do iGC_target=1,Block1%nGC_targets
                nGC_sources_table_local(Block1%GC_targets(iGC_target)%iBlock)=&
                    nGC_sources_table_local(Block1%GC_targets(iGC_target)%iBlock)+1
            end do
        end do


        ! Then reduce it to the global table
        call MPI_ALLReduce(nGC_sources_table_local,nGC_sources_table_global,&
            Tree%NumLeafNodes,MPI_integer,MPI_SUM,MPI_COMM_WORLD,ierr)

        ! Loop the local blocks and allocate GC_sources for them
        do iLocalBlock=1,Tree%nLocalBlocks
            Block1=>Tree%LocalBlocks(iLocalBlock)
            allocate(Block1%GC_sources(nGC_sources_table_global(Block1%iBlock)))
            Block1%nGC_sources=0
        end do

        ! See how many sends should be done for each block.
        ! And then allocate Block1%requests
        do iLocalBlock=1,Tree%nLocalBlocks
            nSend=0
            Block1=>Tree%LocalBlocks(iLocalBlock)
            do iGC_target=1,Block1%nGC_targets
                GC_target1=>Block1%GC_targets(iGC_target)
                if (GC_target1%iRank/=MpiRank) nSend=nSend+1
            end do
            if (nSend>0) allocate(Block1%requests(nSend))
        end do

        ! Then loop the local blocks and send the iBlock of GC_targets
        ! to the target ranks.

        do iLocalBlock=1,Tree%nLocalBlocks
            Block1=>Tree%LocalBlocks(iLocalBlock)
            isend=0
            do iGC_target=1,Block1%nGC_targets
                
                GC_target1=>Block1%GC_targets(iGC_target)
                if (GC_target1%iRank/=MpiRank) then
                    call MPI_ISEND([MpiRank,Block1%iBlock,GC_target1%iBlock],3,mpi_integer,&
                        GC_target1%iRank,1,MPI_COMM_WORLD,request,ierr)
                    
                    iSend=iSend+1
                    Block1%requests(iSend)=request
                else
                    Block2=>Tree%LocalBlocks(GC_target1%iBlock-Tree%iLeafNode_ranges(MpiRank,1)+1)
                    Block2%nGC_sources=Block2%nGC_sources+1
                    GC_source1=>Block2%GC_sources(Block2%nGC_sources)

                    GC_source1%iRank=MpiRank
                    GC_source1%iBlock=Block1%iBlock
                    GC_source1%if_yin=Block1%iBlock.le.Tree%NumLeafNodes_YinYang(1)
                    
                end if
            end do
        end do

        ! Get the n of messages to receive
        nRecv=0
        do iLocalBlock=1,Tree%nLocalBlocks
            Block1=>Tree%LocalBlocks(iLocalBlock)
            nRecv=nRecv+nGC_sources_table_global(Block1%iBlock)-Block1%nGC_sources
        end do

        ! Recv the iBlock messages for GC_sources
        do iRecv=1,nRecv
            call MPI_RECV(recv_message,3,mpi_integer,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,status,ierr)
            Block2=>Tree%LocalBlocks(recv_message(3)-Tree%iLeafNode_ranges(MpiRank,1)+1)
            Block2%nGC_sources=Block2%nGC_sources+1
            GC_source1=>Block2%GC_sources(Block2%nGC_sources)

            GC_source1%iRank=recv_message(1)
            GC_source1%iBlock=recv_message(2)
            GC_source1%if_yin=recv_message(2).le.Tree%NumLeafNodes_YinYang(1)
        end do

        ! At last, wait for all the requests to be completed.
        do iLocalBlock=1,Tree%nLocalBlocks
            Block1=>Tree%LocalBlocks(iLocalBlock)
            do iSend=1,size(Block1%requests)
                call MPI_WAIT(Block1%requests(iSend),status,ierr)
            end do
        end do
    end subroutine ModCommunication_GetnGC_sources

    subroutine ModCommunication_SendRecvGCAll(Tree,if_rk)
        implicit none
        type(YYTree),intent(in),target  ::  Tree
        type(BlockType),pointer         ::  Block1
        logical,intent(in)              ::  if_rk

        integer                         ::  iLocalBlock
        integer                         ::  iGC_target,iGC_source
        integer                         ::  iGC
        type(GC_target),pointer         ::  GC_source1,GC_target1
        real,pointer                    ::  primitive(:,:,:,:)
        real,allocatable                ::  recv_message(:,:)

        integer                         ::  tag,ierr,request,status(MPI_STATUS_SIZE),iSend

        ! Send all the GC_sources of all the 
        ! local blocks in the tree.
        !
        do iLocalBlock=1,Tree%nLocalBlocks
            Block1 => Tree%LocalBlocks(iLocalBlock)
            iSend=0

            ! See if rk primitive
            if (if_rk) then
                primitive => Block1%primitive_rk
            else
                primitive => Block1%primitive
            end if

            ! Then, loop all the GC_sources.
            ! Interpolate the primitive variables and mpi_send them.
            !
            do iGC_source=1,size(Block1%GC_sources)
                GC_source1 => Block1%GC_sources(iGC_source)
                
                ! Send GCs if it's at another MpiRank.
                if (GC_source1%iRank/=MpiRank) then
                    ! interpolate the primitive grid to the GCs
                    call ModMath_1D3D_interpolate_1D1D(primitive,nvar,ni,nj,nk,ng,&
                        Block1%xi,Block1%xj,Block1%xk,GC_source1%nGC,GC_source1%xijk_list,GC_source1%primitive_list)
                
                    ! If different if_yin then convert
                    if (Block1%if_yin .neqv. GC_source1%if_yin) then
                        GC_source1%primitive_list(:,2:4)=&
                            ModYinYang_VecConv_1D(GC_source1%xijk_list,GC_source1%primitive_list(:,2:4),GC_source1%nGC)
                    end if
                
                    ! Get the tag, which contains information for both the send block and receive block.
                    tag = ModCommunication_GetTag(Tree%iLeafNode_ranges,Block1%iBlock,GC_source1%iBlock,GC_source1%iRank)

                    call MPI_ISEND(GC_source1%primitive_list,&
                        nvar*GC_source1%nGC,mpi_real,GC_source1%iRank,tag,MPI_COMM_WORLD,request,ierr)
                    
                    ! Record the request
                    iSend=iSend+1
                    Block1%requests(iSend)=request
                end if 
            end do
        end do

        !print *,'Send done. Now start receiving.'

        ! Recv all the message and assign them to
        ! GC_targets of all the local blocks.
        !
        do iLocalBlock=1,Tree%nLocalBlocks
            Block1 => Tree%LocalBlocks(iLocalBlock)

            ! See if rk primitive
            if (if_rk) then
                primitive => Block1%primitive_rk
            else
                primitive => Block1%primitive
            end if
    
            do iGC_target=1,size(Block1%GC_targets)
                if (Block1%GC_targets(iGC_target)%iRank.ne.MpiRank) then
                    GC_target1 => Block1%GC_targets(iGC_target)
    
                    allocate(recv_message(GC_target1%nGC,nvar))
                    !allocate(recv_message_1(GC_target1%nGC))
    
                    tag=ModCommunication_GetTag(Tree%iLeafNode_ranges,GC_target1%iBlock,Block1%iBlock,MpiRank)

                    call MPI_RECV(recv_message,nvar*GC_target1%nGC,mpi_real,GC_target1%iRank,tag,MPI_COMM_WORLD,status,ierr)

                    do iGC = 1,GC_target1%nGC
                        primitive(  GC_target1%ijk_list(iGC,1),&
                                    GC_target1%ijk_list(iGC,2),&
                                    GC_target1%ijk_list(iGC,3),:)=recv_message(iGC,:)
                        !if(Block1%iBlock==6 .and. GC_target1%iBlock==14 .and. &
                        !    GC_target1%ijk_list(iGC,1)==32 .and. &
                        !    GC_target1%ijk_list(iGC,2)== 0 .and. &
                        !    GC_target1%ijk_list(iGC,3)== 1) print *,recv_message(iGC,3),iGC
                    end do
                    deallocate(recv_message)
                end if
            end do
        end do

        ! wait for all requests to be completed.

        do iLocalBlock=1,Tree%nLocalBlocks
            Block1=>Tree%LocalBlocks(iLocalBlock)
            do iSend=1,size(Block1%requests)
                call MPI_WAIT(Block1%requests(iSend),status,ierr)
            end do
        end do

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    end subroutine ModCommunication_SendRecvGCAll

    function ModCommunication_GetTag(iLeafNode_ranges,iBlockSend,iBlockRecv,iRankRecv) result(Tag)
        implicit none
        integer,intent(in)      :: iLeafNode_ranges(0:MpiSize-1,2)
        integer,intent(in)      :: iBlockSend,iBlockRecv
        integer,intent(in)      :: iRankRecv
        
        integer                 :: tag

        Tag=nMaxBlocksPerRank*iBlockSend+(iBlockRecv-iLeafNode_ranges(iRankRecv,1))
    end function ModCommunication_GetTag

    function ModCommunication_SolveTag(iLeafNode_ranges,Tag,iRankRecv) result(iBlockSendRecv)
        implicit none
        integer,intent(in)      :: iLeafNode_ranges(0:MpiSize-1,2)
        integer,intent(in)      :: Tag
        integer,intent(in)      :: iRankRecv

        integer                 :: iBlockSendRecv(2)

        iBlockSendRecv(1)=Tag/nMaxBlocksPerRank
        iBlockSendRecv(2)=mod(Tag,nMaxBlocksPerRank)+iLeafNode_ranges(iRankRecv,1)
    end function ModCommunication_SolveTag
end Module ModCommunication