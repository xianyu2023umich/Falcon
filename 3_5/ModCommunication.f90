module ModCommunication

    use ModMath,        only:   ModMath_1D3D_interpolate_1D1D,&
                                ModMath_3D_interpolate_1D
    use ModBlock,       only:   BlockType,Multigrid_level,GC_target
    use ModYinYang,     only:   ModYinYang_CoordConv_1D,&
                                ModYinYang_CoordConv_0D,&
                                ModYinYang_VecConv_1D
    use ModYinYangTree, only:   YYTree,GC_interface
    use ModGC,          only:   ModGC_FindGC_Single,&
                                ModGC_FindHC_Single,&
                                ModGC_GetGC_Targets_single,&
                                ModGC_GetHC_Targets_single
    use ModVariables,   only:   br_,bt_,bp_,vr_,vt_,vp_
    use ModParameters,  only:   MpiSize,MpiRank,&
                                ni,nj,nk,ng,nvar,Multigrid_nLevels
    use ModControl,     only:   nMaxBlocksPerRank
    use ModAllocation,  only:   ranges_of_ranks
    use ModMultigrid,   only:   ModMultigrid_InitAll
    use MPI

    implicit none

    contains

    subroutine ModCommunication_SetGCAll(Tree,if_HC)
        implicit none
        type(YYTree),target             ::  Tree                    ! the Tree
        logical,intent(in)              ::  if_HC                   ! if HC is needed
        type(BlockType),pointer         ::  Block1                  ! one Block
        integer                         ::  iBlock
        integer                         ::  iLevel
        type(Multigrid_level),pointer   ::  MGL1                    ! one Multigrid level
        integer                         ::  ierr

        ! First determine where each ghost cell belongs to
        ! by using the ModGC_FindGC_Single subroutine.
        do iBlock=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iBlock)
            call ModGC_FindGC_Single(Tree,Block1)
        end do
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (MpiRank==0) write(*,*)'Completed finding GC...'

        ! Then get the GC_targets for each block
        do iBlock=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iBlock)
            call ModGC_GetGC_Targets_single(Tree,Block1)
        end do
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (MpiRank==0) write(*,*)'Completed setting GC_targets...',MpiRank

        ! Then set the GC_sources for each block
        call ModCommunication_GetnGC_sources(Tree)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (MpiRank==0) write(*,*)'Completed setting nGC_sources...',MpiRank      

        ! Then set the GC_sources for each block
        call ModCommunication_SetGC_sources(Tree)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (MpiRank==0) write(*,*)'Completed setting GC_sources...',MpiRank

        ! Setting up GC senders and receivers
        call ModCommunication_SetGC_Sender(Tree)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (MpiRank==0) write(*,*)'Completed setting GC_senders...',MpiRank

        call ModCommunication_SetGC_Receiver(Tree)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (MpiRank==0) write(*,*)'Completed setting GC_receivers...',MpiRank

        if (if_HC) then
            ! Initialize multigrid
            call ModMultigrid_InitAll(Tree,10,.true.)
            if(MpiRank==0) write(*,*)'Completed initializing multigrid...'

            ! And now comes multigrid.
            do iLevel=1,Multigrid_nLevels
                if (MpiRank==0) write(*,*)'--------------------------------'
                if (MpiRank==0) write(*,*)'Starting multigrid for level ',iLevel,'...'

                do iBlock=1,size(Tree%LocalBlocks)
                    MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
                    call ModGC_FindHC_Single(Tree,MGL1)
                end do
                if (MpiRank==0) write(*,*)'Completed finding HC for level...'

                do iBlock=1,size(Tree%LocalBlocks)
                    MGL1=>Tree%LocalBlocks(iBlock)%Multigrid_levels(iLevel)
                    call ModGC_GetHC_Targets_single(Tree,MGL1)
                end do
                if (MpiRank==0) write(*,*)'Completed setting HC_targets...'

                call ModCommunication_GetnHC_sources(Tree,iLevel)
                if (MpiRank==0) write(*,*)'Completed setting nHC_sources...'

                call ModCommunication_SetHC_sources(Tree,iLevel)
                if (MpiRank==0) write(*,*)'Completed setting HC_sources...'

                if (MpiRank==0) write(*,*)'Completed multigrid for level ',iLevel,'...'
                if (MpiRank==0) write(*,*)'--------------------------------'
            end do

            call ModCommunication_SetHC_Sender(Tree)
            if (MpiRank==0) write(*,*)'Completed setting HC_senders...'

            call ModCommunication_SetHC_Receiver(Tree)
            if (MpiRank==0) write(*,*)'Completed setting HC_receivers...'        

        end if        
    end subroutine ModCommunication_SetGCAll

    subroutine ModCommunication_GetnGC_sources(Tree)
        implicit none
        type(YYTree),intent(in),target  ::  Tree
        type(BlockType),pointer         ::  Block_source,Block_target
        type(GC_target),pointer         ::  GC_target1,GC_source1
        integer                         ::  nGC_sources_table_local(Tree%NumLeafNodes)
        integer                         ::  nGC_sources_table_global(Tree%NumLeafNodes)
        integer                         ::  nGC_international_sources_table_local(Tree%NumLeafNodes)
        integer                         ::  nGC_international_sources_table_global(Tree%NumLeafNodes)
        integer                         ::  iLocalBlock,iGC_target,iGC_source
        integer                         ::  iRecv,nRecv
        integer,allocatable             ::  send_message_all(:),recv_message_all(:)
        integer                         ::  ierr,request,status(MPI_STATUS_SIZE)
        integer                         ::  nSend,iSend
        integer,allocatable             ::  requests_send(:),requests_recv(:)
        
        ! Get the local table for GC_sources of each block
        nGC_sources_table_local=0
        nGC_international_sources_table_local=0
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)

            do iGC_target=1,Block_source%nGC_targets
                GC_target1=>Block_source%GC_targets(iGC_target)

                ! This is all the GC_sources.
                nGC_sources_table_local(GC_target1%iBlock)=&
                        nGC_sources_table_local(GC_target1%iBlock)+1

                ! Only the international GC_sources.
                if (MpiRank/=GC_target1%iRank) then
                    nGC_international_sources_table_local(GC_target1%iBlock)=&
                        nGC_international_sources_table_local(GC_target1%iBlock)+1
                end if
            end do
        end do

        ! Then reduce it to the global table
        call MPI_ALLReduce(nGC_sources_table_local,nGC_sources_table_global,&
            Tree%NumLeafNodes,MPI_integer,MPI_SUM,MPI_COMM_WORLD,ierr)

        call MPI_ALLReduce(nGC_international_sources_table_local,nGC_international_sources_table_global,&
            Tree%NumLeafNodes,MPI_integer,MPI_SUM,MPI_COMM_WORLD,ierr)

        ! Loop the local blocks and allocate GC_sources for them
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)
            allocate(Block_source%GC_sources(nGC_sources_table_global(Block_source%iBlock)))
            Block_source%nGC_sources=0
        end do

        ! Before the loop, we need to get the n of messages to send
        nSend=0
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)
            do iGC_target=1,Block_source%nGC_targets
                GC_target1=>Block_source%GC_targets(iGC_target)

                if (GC_target1%iRank/=MpiRank) then
                    nSend=nSend+1
                end if
            end do
        end do

        ! Get the n of messages to receive
        nRecv=0
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)
            nRecv=nRecv+nGC_international_sources_table_global(Block_source%iBlock)
        end do

        allocate(send_message_all(nSend*4))
        allocate(recv_message_all(nRecv*4))
        allocate(requests_send(nSend))
        allocate(requests_recv(nRecv))

        ! Recv the iBlock messages for GC_sources
        do iRecv=1,nRecv
            call MPI_IRECV(recv_message_all(iRecv*4-3:iRecv*4),4,mpi_integer,&
                MPI_ANY_SOURCE,11,MPI_COMM_WORLD,requests_recv(iRecv),ierr)
        end do

        ! Then loop the local blocks and send the iBlock of GC_targets
        ! to the target ranks.

        iSend=0
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)
            do iGC_target=1,Block_source%nGC_targets
                GC_target1=>Block_source%GC_targets(iGC_target)

                if (GC_target1%iRank/=MpiRank) then
                    iSend=iSend+1
                    send_message_all(iSend*4-3:iSend*4)=[MpiRank,Block_source%iBlock,GC_target1%iBlock,GC_target1%nGC]

                    call MPI_ISEND(send_message_all(iSend*4-3:iSend*4),4,mpi_integer,&
                        GC_target1%iRank,11,MPI_COMM_WORLD,requests_send(iSend),ierr)
                else
                    Block_target=>Tree%LocalBlocks(GC_target1%iBlock-ranges_of_ranks(MpiRank,1)+1)
                    Block_target%nGC_sources=Block_target%nGC_sources+1
                    GC_source1=>Block_target%GC_sources(Block_target%nGC_sources)

                    GC_source1%iRank=MpiRank
                    GC_source1%iBlock=Block_source%iBlock
                    GC_source1%if_yin=Block_source%iBlock.le.Tree%NumLeafNodes_YinYang(1)
                    GC_source1%nGC=GC_target1%nGC
                    allocate(GC_source1%xijk_list(GC_source1%nGC,3))
                    allocate(GC_source1%primitive_list(GC_source1%nGC,nvar))
                end if
            end do
        end do
        ! Wait for all the sends to be done
        call MPI_WAITALL(nSend,requests_send,MPI_STATUSES_IGNORE,ierr)
        call MPI_WAITALL(nRecv,requests_recv,MPI_STATUSES_IGNORE,ierr)
        
        ! Assign and allocate GC_sources
        do iRecv=1,nRecv
            Block_target=>Tree%LocalBlocks(recv_message_all(iRecv*4-1)-ranges_of_ranks(MpiRank,1)+1)
            Block_target%nGC_sources=Block_target%nGC_sources+1

            GC_source1=>Block_target%GC_sources(Block_target%nGC_sources)

            GC_source1%iRank=recv_message_all(iRecv*4-3)
            GC_source1%iBlock=recv_message_all(iRecv*4-2)
            GC_source1%if_yin=recv_message_all(iRecv*4-2).le.Tree%NumLeafNodes_YinYang(1)
            GC_source1%nGC=recv_message_all(iRecv*4)
            allocate(GC_source1%xijk_list(GC_source1%nGC,3))
            allocate(GC_source1%primitive_list(GC_source1%nGC,nvar))
        end do

        deallocate(requests_recv)
        deallocate(requests_send)
        deallocate(send_message_all)
        deallocate(recv_message_all)

        ! At last we get the total n of requests

        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)
            nSend=0
            if (Block_source%nGC_sources>0) then
                do iGC_source=1,Block_source%nGC_sources
                    if (Block_source%GC_sources(iGC_source)%iRank/=MpiRank) nSend=nSend+1
                end do
                if (nSend>0) allocate(Block_source%GC_requests(nSend))
            endif 
        end do
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    end subroutine ModCommunication_GetnGC_sources

    subroutine ModCommunication_SetGC_sources(Tree)
        implicit none
        type(YYTree),intent(in),target  ::  Tree
        type(BlockType),pointer         ::  Block_source,Block_target

        integer                         ::  iLocalBlock
        integer                         ::  iGC_target,iGC_source
        type(GC_target),pointer         ::  GC_source1,GC_target1
        real,allocatable                ::  recv_message_real(:,:)
        integer,allocatable             ::  recv_message_int(:,:)

        integer                         ::  tag,ierr,request,status(MPI_STATUS_SIZE),count

        ! Send the GC_targets%xijk_list

        do iLocalBlock=1,Tree%nLocalBlocks
            Block_target => Tree%LocalBlocks(iLocalBlock)

            do iGC_target=1,size(Block_target%GC_targets)
                GC_target1 => Block_target%GC_targets(iGC_target)

                if (GC_target1%iRank/=MpiRank) then
                    tag = ModCommunication_GetTag(Block_target%iBlock,GC_target1%iBlock,MpiRank,GC_target1%iRank)
                    call MPI_ISEND(GC_target1%xijk_list,&
                        3*GC_target1%nGC,mpi_real,GC_target1%iRank,tag,MPI_COMM_WORLD,request,ierr)
                else
                    Block_source=>Tree%LocalBlocks(GC_target1%iBlock-ranges_of_ranks(MpiRank,1)+1)

                    ! Go over the GC_sources to see which one corresponds to the GC_target1.
                    do iGC_source=1,size(Block_source%GC_sources)
                        GC_source1=>Block_source%GC_sources(iGC_source)
                        if (GC_source1%iBlock==Block_target%iBlock) then
                            ! First record the ijk_list. It always means the ijk
                            ! at the target side, so it's the same for source.
                            GC_source1%ijk_list=GC_target1%ijk_list

                            ! Then record the xijk_list.
                            ! xijk_list is first calculated at the target side.
                            ! But it's later used at the source side for interpolation.
                            ! So we need to convert it to the source side with if_yin.
                            if (GC_source1%if_yin .eqv. Block_source%if_yin) then
                                GC_source1%xijk_list=GC_target1%xijk_list
                            else
                                GC_source1%xijk_list=ModYinYang_CoordConv_1D(GC_target1%xijk_list,GC_target1%nGC)
                            end if
                            exit
                        end if
                    end do
                end if ! GC_target1%iRank/=MpiRank
            end do ! iGC_target
        end do ! iLocalBlock
        

        ! Set GC_sources
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source => Tree%LocalBlocks(iLocalBlock)
    
            do iGC_source=1,size(Block_source%GC_sources)
                GC_source1 => Block_source%GC_sources(iGC_source)

                ! Local sends have been done before so only deal with the remote sends.
                
                if (GC_source1%iRank/=MpiRank) then
                    ! Allocate the buffer, get the tag and then receive
                    allocate(recv_message_real(GC_source1%nGC,3))
                    tag=ModCommunication_GetTag(GC_source1%iBlock,Block_source%iBlock,GC_source1%iRank,MpiRank)
                    call MPI_RECV(recv_message_real,3*GC_source1%nGC,mpi_real,GC_source1%iRank,tag,MPI_COMM_WORLD,status,ierr)

                    ! See whether same if_yin or not
                    if (GC_source1%if_yin .eqv. Block_source%if_yin) then
                        GC_source1%xijk_list=recv_message_real
                    else
                        GC_source1%xijk_list=ModYinYang_CoordConv_1D(recv_message_real,GC_source1%nGC)
                    endif
                    deallocate(recv_message_real)
                endif ! GC_source1%iRank/=MpiRank
            end do
        end do

        ! But -- these are just for xijk_list for remote sends.
        ! Do again for ijk_list.
        ! There's no need to care local sends, cuz local sends before include both xijk_list and ijk_list.

        do iLocalBlock=1,Tree%nLocalBlocks
            Block_target => Tree%LocalBlocks(iLocalBlock)
            do iGC_target=1,size(Block_target%GC_targets)
                GC_target1 => Block_target%GC_targets(iGC_target)
                if (GC_target1%iRank/=MpiRank) then
                    tag = ModCommunication_GetTag(Block_target%iBlock,GC_target1%iBlock,MpiRank,GC_target1%iRank)
                    call MPI_ISEND(GC_target1%ijk_list,&
                        3*GC_target1%nGC,mpi_integer,GC_target1%iRank,tag,MPI_COMM_WORLD,request,ierr)
                end if ! GC_target1%iRank/=MpiRank
            end do ! iGC_target
        end do ! iLocalBlock

        ! Set GC_sources
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source => Tree%LocalBlocks(iLocalBlock)
            do iGC_source=1,size(Block_source%GC_sources)
                GC_source1 => Block_source%GC_sources(iGC_source)

                ! Local sends have been done before so only deal with the remote sends.
                
                if (GC_source1%iRank/=MpiRank) then
                    ! Allocate the buffer, get the tag and then receive
                    allocate(recv_message_int(GC_source1%nGC,3))
                    tag=ModCommunication_GetTag(GC_source1%iBlock,Block_source%iBlock,GC_source1%iRank,MpiRank)
                    call MPI_RECV(recv_message_int,3*GC_source1%nGC,mpi_integer,GC_source1%iRank,tag,MPI_COMM_WORLD,status,ierr)

                    ! Assign the received message to the ijk_list. Then deallocate the buffer.
                    GC_source1%ijk_list=recv_message_int
                    deallocate(recv_message_int)
                endif ! GC_source1%iRank/=MpiRank
            end do
        end do
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    end subroutine ModCommunication_SetGC_sources

    subroutine ModCommunication_GetnHC_sources(Tree,iLevel)
        implicit none
        type(YYTree),intent(in),target  ::  Tree
        integer,intent(in)              ::  iLevel
        type(Multigrid_level),pointer   ::  MGL_source,MGL_target
        type(GC_target),pointer         ::  HC_target1,HC_source1
        integer                         ::  nHC_sources_table_local(Tree%NumLeafNodes)
        integer                         ::  nHC_sources_table_global(Tree%NumLeafNodes)
        integer                         ::  iLocalBlock,iHC_target
        integer                         ::  iRecv,nRecv,recv_message(4)
        integer                         ::  ierr,request,status(MPI_STATUS_SIZE)

        ! Get the local table for HC_sources of iLevelth level of each block
        nHC_sources_table_local=0
        do iLocalBlock=1,Tree%nLocalBlocks
            MGL_source=>Tree%LocalBlocks(iLocalBlock)%Multigrid_levels(iLevel)
            do iHC_target=1,MGL_source%nHC_targets
                nHC_sources_table_local(MGL_source%HC_targets(iHC_target)%iBlock)=&
                    nHC_sources_table_local(MGL_source%HC_targets(iHC_target)%iBlock)+1
            end do
        end do

        ! Then reduce it to the global table
        call MPI_ALLReduce(nHC_sources_table_local,nHC_sources_table_global,&
            Tree%NumLeafNodes,MPI_integer,MPI_SUM,MPI_COMM_WORLD,ierr)

        ! Loop the local blocks and allocate HC_sources for them
        do iLocalBlock=1,Tree%nLocalBlocks
            MGL_source=>Tree%LocalBlocks(iLocalBlock)%Multigrid_levels(iLevel)
            allocate(MGL_source%HC_sources(nHC_sources_table_global(MGL_source%iBlock)))
            MGL_source%nHC_sources=0
        end do

        ! Then loop the local blocks and send the iBlock of GC_targets
        ! to the target ranks.

        do iLocalBlock=1,Tree%nLocalBlocks
            MGL_source=>Tree%LocalBlocks(iLocalBlock)%Multigrid_levels(iLevel)
            do iHC_target=1,MGL_source%nHC_targets
                HC_target1=>MGL_source%HC_targets(iHC_target)
                
                if (HC_target1%iRank/=MpiRank) then
                    call MPI_ISEND([MpiRank,MGL_source%iBlock,HC_target1%iBlock,HC_target1%nGC],4,mpi_integer,&
                        HC_target1%iRank,1,MPI_COMM_WORLD,request,ierr)
                else
                    MGL_target=>Tree%LocalBlocks(HC_target1%iBlock-ranges_of_ranks(MpiRank,1)+1)%Multigrid_levels(iLevel)
                    MGL_target%nHC_sources=MGL_target%nHC_sources+1
                    HC_source1=>MGL_target%HC_sources(MGL_target%nHC_sources)

                    HC_source1%iRank=MpiRank
                    HC_source1%iBlock=MGL_source%iBlock
                    HC_source1%if_yin=MGL_source%iBlock.le.Tree%NumLeafNodes_YinYang(1)
                    HC_source1%nGC=HC_target1%nGC
                    allocate(HC_source1%xijk_list(HC_source1%nGC,3))
                    allocate(HC_source1%primitive_list(HC_source1%nGC,1))
                end if
            end do
        end do

        ! Get the n of messages to receive
        nRecv=0
        do iLocalBlock=1,Tree%nLocalBlocks
            MGL_source=>Tree%LocalBlocks(iLocalBlock)%Multigrid_levels(iLevel)
            nRecv=nRecv+nHC_sources_table_global(MGL_source%iBlock)-MGL_source%nHC_sources
        end do

        ! Recv the iBlock messages for GC_sources
        do iRecv=1,nRecv
            call MPI_RECV(recv_message,4,mpi_integer,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,status,ierr)
            MGL_target=>Tree%LocalBlocks(recv_message(3)-ranges_of_ranks(MpiRank,1)+1)%Multigrid_levels(iLevel)
            MGL_target%nHC_sources=MGL_target%nHC_sources+1
            HC_source1=>MGL_target%HC_sources(MGL_target%nHC_sources)

            HC_source1%iRank=recv_message(1)
            HC_source1%iBlock=recv_message(2)
            HC_source1%if_yin=recv_message(2).le.Tree%NumLeafNodes_YinYang(1)
            HC_source1%nGC=recv_message(4)
            allocate(HC_source1%xijk_list(HC_source1%nGC,3))
            allocate(HC_source1%primitive_list(HC_source1%nGC,1))
        end do

        

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    end subroutine ModCommunication_GetnHC_sources

    subroutine ModCommunication_SetHC_sources(Tree,iLevel)
        implicit none
        type(YYTree),intent(in),target  ::  Tree
        integer,intent(in)              ::  iLevel
        type(Multigrid_level),pointer   ::  MGL_source,MGL_target

        integer                         ::  iLocalBlock
        integer                         ::  iHC_target,iHC_source
        type(GC_target),pointer         ::  HC_source1,HC_target1
        real,allocatable                ::  recv_message_real(:,:)
        integer,allocatable             ::  recv_message_int(:,:)

        integer                         ::  tag,ierr,request,status(MPI_STATUS_SIZE)

        ! Send the HC_targets%xijk_list

        do iLocalBlock=1,Tree%nLocalBlocks
            MGL_target => Tree%LocalBlocks(iLocalBlock)%Multigrid_levels(iLevel)

            do iHC_target=1,size(MGL_target%HC_targets)
                HC_target1 => MGL_target%HC_targets(iHC_target)

                if (HC_target1%iRank/=MpiRank) then
                    tag = ModCommunication_GetTag(MGL_target%iBlock,HC_target1%iBlock,MpiRank,HC_target1%iRank)
                    call MPI_ISEND(HC_target1%xijk_list,&
                        3*HC_target1%nGC,mpi_real,HC_target1%iRank,tag,MPI_COMM_WORLD,request,ierr)
                else
                    MGL_source=>Tree%LocalBlocks(HC_target1%iBlock-ranges_of_ranks(MpiRank,1)+1)%Multigrid_levels(iLevel)

                    ! Go over the HC_sources to see which one corresponds to the HC_target1.
                    do iHC_source=1,size(MGL_source%HC_sources)
                        HC_source1=>MGL_source%HC_sources(iHC_source)
                        if (HC_source1%iBlock==MGL_target%iBlock) then
                            ! First record the ijk_list. It always means the ijk
                            ! at the target side, so it's the same for source.
                            HC_source1%ijk_list=HC_target1%ijk_list

                            ! Then record the xijk_list.
                            ! xijk_list is first calculated at the target side.
                            ! But it's later used at the source side for interpolation.
                            ! So we need to convert it to the source side with if_yin.
                            if (HC_source1%if_yin .eqv. MGL_source%if_yin) then
                                HC_source1%xijk_list=HC_target1%xijk_list
                            else
                                HC_source1%xijk_list=ModYinYang_CoordConv_1D(HC_target1%xijk_list,HC_target1%nGC)
                            end if
                            exit
                        end if
                    end do
                end if ! HC_target1%iRank/=MpiRank
            end do ! iHC_target
        end do ! iLocalBlock
        

        ! Set HC_sources
        do iLocalBlock=1,Tree%nLocalBlocks
            MGL_source => Tree%LocalBlocks(iLocalBlock)%Multigrid_levels(iLevel)
    
            do iHC_source=1,size(MGL_source%HC_sources)
                HC_source1 => MGL_source%HC_sources(iHC_source)

                ! Local sends have been done before so only deal with the remote sends.
                
                if (HC_source1%iRank/=MpiRank) then
                    ! Allocate the buffer, get the tag and then receive
                    allocate(recv_message_real(HC_source1%nGC,3))
                    tag=ModCommunication_GetTag(HC_source1%iBlock,MGL_source%iBlock,HC_source1%iRank,MpiRank)
                    call MPI_RECV(recv_message_real,3*HC_source1%nGC,mpi_real,HC_source1%iRank,tag,MPI_COMM_WORLD,status,ierr)

                    ! See whether same if_yin or not
                    if (HC_source1%if_yin .eqv. MGL_source%if_yin) then
                        HC_source1%xijk_list=recv_message_real
                    else
                        HC_source1%xijk_list=ModYinYang_CoordConv_1D(recv_message_real,HC_source1%nGC)
                    endif
                    deallocate(recv_message_real)
                endif ! HC_source1%iRank/=MpiRank
            end do
        end do

        ! But -- these are just for xijk_list for remote sends.
        ! Do again for ijk_list.
        ! There's no need to care local sends, cuz local sends before include both xijk_list and ijk_list.

        do iLocalBlock=1,Tree%nLocalBlocks
            MGL_target => Tree%LocalBlocks(iLocalBlock)%Multigrid_levels(iLevel)
            do iHC_target=1,size(MGL_target%HC_targets)
                HC_target1 => MGL_target%HC_targets(iHC_target)
                if (HC_target1%iRank/=MpiRank) then
                    tag = ModCommunication_GetTag(MGL_target%iBlock,HC_target1%iBlock,MpiRank,HC_target1%iRank)
                    call MPI_ISEND(HC_target1%ijk_list,&
                        3*HC_target1%nGC,mpi_integer,HC_target1%iRank,tag,MPI_COMM_WORLD,request,ierr)
                end if ! HC_target1%iRank/=MpiRank
            end do ! iHC_target
        end do ! iLocalBlock

        ! Set HC_sources
        do iLocalBlock=1,Tree%nLocalBlocks
            MGL_source => Tree%LocalBlocks(iLocalBlock)%Multigrid_levels(iLevel)
            do iHC_source=1,size(MGL_source%HC_sources)
                HC_source1 => MGL_source%HC_sources(iHC_source)

                ! Local sends have been done before so only deal with the remote sends.
                
                if (HC_source1%iRank/=MpiRank) then
                    ! Allocate the buffer, get the tag and then receive
                    allocate(recv_message_int(HC_source1%nGC,3))
                    tag=ModCommunication_GetTag(HC_source1%iBlock,MGL_source%iBlock,HC_source1%iRank,MpiRank)
                    call MPI_RECV(recv_message_int,3*HC_source1%nGC,mpi_integer,HC_source1%iRank,tag,MPI_COMM_WORLD,status,ierr)

                    ! Assign the received message to the ijk_list. Then deallocate the buffer.
                    HC_source1%ijk_list=recv_message_int
                    deallocate(recv_message_int)
                endif ! HC_source1%iRank/=MpiRank
            end do ! iHC_source
        end do ! iLocalBlock
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    end subroutine ModCommunication_SetHC_sources

    ! Set the GC sender for the tree. 
    ! First see all the distinct neighbour ranks. The number is the number of 
    ! senders to be allocated.
    ! Then see how many block pairs are there. There's NO NEED to see how many
    ! *distinct pairs* because different GC_sources always correspond to
    ! distinct pairs.

    subroutine ModCommunication_SetGC_Sender(Tree)
        implicit none
        type(YYTree),target             ::  Tree
        integer                         ::  iUniqueMpiRank
        integer                         ::  UniqueMpiRank_ptr
        integer                         ::  UniqueMpiRanks(MpiSize)

        ! iBlock_pair means different so choose another name
        integer                         ::  iofBlock_pair,iGC                   

        integer                         ::  iGC_Sender
        integer                         ::  iLocalBlock,iGC_source
        type(BlockType),pointer         ::  Block_source
        type(GC_target),pointer         ::  GC_source1
        type(GC_interface),pointer      ::  GC_Sender1

        ! Loop all local Blocks%GC_sources to find all neighbour ranks.
        ! First set UniqueMpiRanks to be -1. Loop, if find any new neighbour,
        ! record it there.

        UniqueMpiRank_ptr=1
        UniqueMpiRanks=-1

        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)
            do iGC_source=1,Block_source%nGC_sources
                GC_source1=>Block_source%GC_sources(iGC_source)

                if (GC_source1%iRank/=MpiRank) then
                    ! Loop all recorded unique ranks to see if it's already in there.
                    ! If not, record it.
                    do iUniqueMpiRank=1,UniqueMpiRank_ptr
                        if (UniqueMpiRanks(iUniqueMpiRank)==GC_source1%iRank) exit
                    end do
                    if (iUniqueMpiRank==UniqueMpiRank_ptr+1) then
                        UniqueMpiRanks(UniqueMpiRank_ptr)=GC_source1%iRank
                        UniqueMpiRank_ptr=UniqueMpiRank_ptr+1
                    end if
                end if
            end do
        end do

        ! Now the nPairs is UniqueMpiRank_ptr-1.
        ! Allocate the GC_Senders.
        Tree%nGC_Senders=UniqueMpiRank_ptr-1
        allocate(Tree%GC_Senders(Tree%nGC_Senders))
        Tree%GC_Senders(:)%iRank=UniqueMpiRanks(1:Tree%nGC_Senders)

        ! Get the iBlock_pairs for each GC_sender.
        ! Loop all the local blocks and GC_sources to get the unique pairs.
        ! The process is pretty similar to above, but we need a 2D array 
        ! for recording since it's a pair of blocks.
        do iGC_Sender=1,Tree%nGC_Senders

            GC_Sender1=>Tree%GC_Senders(iGC_Sender)
            GC_Sender1%nBlock_Pairs=0
            GC_Sender1%nGC=0

            do iLocalBlock=1,Tree%nLocalBlocks
                Block_source=>Tree%LocalBlocks(iLocalBlock)
                do iGC_source=1,Block_source%nGC_sources
                    GC_source1=>Block_source%GC_sources(iGC_source)

                    ! No need to record. Different GC_sources always correspond to
                    ! different block pairs.

                    if (GC_source1%iRank==GC_Sender1%iRank) then
                        GC_Sender1%nBlock_Pairs=GC_Sender1%nBlock_Pairs+1
                        GC_Sender1%nGC=GC_Sender1%nGC+GC_source1%nGC
                    end if

                end do
            end do

            allocate(GC_Sender1%message(GC_Sender1%nGC,nvar))
            allocate(GC_Sender1%iBlock_pairs(GC_Sender1%nBlock_Pairs,2))
            allocate(GC_Sender1%Block_Pairs_Ptrs(GC_Sender1%nBlock_Pairs,2))
            allocate(GC_Sender1%iLocalBlock_and_iGCtarget_list(GC_Sender1%nBlock_Pairs,2))

            !Loop again to find the iBlock_pairs and the Block_Pairs_Ptrs.
            !We need iGC and iBlock_pair to record the current position.

            iofBlock_pair=1
            iGC=1
            do iLocalBlock=1,Tree%nLocalBlocks
                Block_source=>Tree%LocalBlocks(iLocalBlock)
                do iGC_source=1,Block_source%nGC_sources
                    GC_source1=>Block_source%GC_sources(iGC_source)

                    if (GC_source1%iRank==GC_Sender1%iRank) then
                        ! Record the iBlock_pairs.
                        GC_Sender1%iBlock_pairs(iofBlock_pair,1)=Block_source%iBlock
                        GC_Sender1%iBlock_pairs(iofBlock_pair,2)=GC_source1%iBlock

                        ! Ptrs
                        GC_Sender1%Block_Pairs_Ptrs(iofBlock_pair,1)=iGC
                        iGC=iGC+GC_source1%nGC
                        GC_Sender1%Block_Pairs_Ptrs(iofBlock_pair,2)=iGC-1

                        ! iLocalBlock and iGCtarget
                        GC_Sender1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,1)=iLocalBlock
                        GC_Sender1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,2)=iGC_source

                        iofBlock_pair=iofBlock_pair+1
                    end if
                end do
            end do
        end do
    end subroutine ModCommunication_SetGC_Sender

    ! Then set the GC_receivers.
    ! First find neighbours like the GC_Senders so we know where to receive.
    ! Then find how many Block pairs in each Receiver is still the same.
    ! But the iBlock_pairs and ptrs are different -- we will get them by
    ! sending them from GC_Senders.

    subroutine ModCommunication_SetGC_Receiver(Tree)
        implicit none
        type(YYTree),target             ::  Tree
        integer                         ::  iUniqueMpiRank
        integer                         ::  UniqueMpiRank_ptr
        integer                         ::  UniqueMpiRanks(MpiSize)

        integer                         ::  iofBlock_pair

        integer                         ::  iLocalBlock,iGC_target,iGC_sender,iGC_receiver
        type(BlockType),pointer         ::  Block_target
        type(GC_target),pointer         ::  GC_target1
        type(GC_interface),pointer      ::  GC_Sender1,GC_Receiver1
        integer                         ::  ierr,request,status(MPI_STATUS_SIZE)
        integer,allocatable             ::  requests(:)

        ! Loop all local Blocks%GC_sources to find all neighbour ranks.
        ! First set UniqueMpiRanks to be -1. Loop, if find any new neighbour,
        ! record it there.

        UniqueMpiRank_ptr=1
        UniqueMpiRanks=-1

        do iLocalBlock=1,Tree%nLocalBlocks
            Block_target=>Tree%LocalBlocks(iLocalBlock)
            do iGC_target=1,Block_target%nGC_targets
                GC_target1=>Block_target%GC_targets(iGC_target)

                if (GC_target1%iRank/=MpiRank) then
                    ! Loop all recorded unique ranks to see if it's already in there.
                    ! If not, record it.
                    do iUniqueMpiRank=1,UniqueMpiRank_ptr
                        if (UniqueMpiRanks(iUniqueMpiRank)==GC_target1%iRank) exit
                    end do
                    if (iUniqueMpiRank==UniqueMpiRank_ptr+1) then
                        UniqueMpiRanks(UniqueMpiRank_ptr)=GC_target1%iRank
                        UniqueMpiRank_ptr=UniqueMpiRank_ptr+1
                    end if
                end if
            end do
        end do

        ! Now the nPairs is UniqueMpiRank_ptr-1.
        ! Allocate the GC_Receivers.
        Tree%nGC_Receivers=UniqueMpiRank_ptr-1
        allocate(Tree%GC_Receivers(Tree%nGC_Receivers))
        Tree%GC_Receivers(:)%iRank=UniqueMpiRanks(1:Tree%nGC_Receivers)

        ! Get the iBlock_pairs for each GC_receiver.
        ! Loop all the local blocks and GC_sources to get the unique pairs.
        ! The process is pretty similar to above, but we need a 2D array 
        ! for recording since it's a pair of blocks.

        do iGC_receiver=1,Tree%nGC_Receivers

            GC_Receiver1=>Tree%GC_Receivers(iGC_receiver)
            GC_Receiver1%nBlock_Pairs=0
            GC_Receiver1%nGC=0

            do iLocalBlock=1,Tree%nLocalBlocks
                Block_target=>Tree%LocalBlocks(iLocalBlock)
                do iGC_target=1,Block_target%nGC_targets
                    GC_target1=>Block_target%GC_targets(iGC_target)

                    ! No need to record. Different GC_sources always correspond to
                    ! different block pairs.

                    if (GC_target1%iRank==GC_Receiver1%iRank) then
                        GC_Receiver1%nBlock_Pairs=GC_Receiver1%nBlock_Pairs+1
                        GC_Receiver1%nGC=GC_Receiver1%nGC+GC_target1%nGC
                    end if

                end do
            end do

            allocate(GC_Receiver1%message(GC_Receiver1%nGC,nvar))
            allocate(GC_Receiver1%iBlock_pairs(GC_Receiver1%nBlock_Pairs,2))
            allocate(GC_Receiver1%Block_Pairs_Ptrs(GC_Receiver1%nBlock_Pairs,2))
            allocate(GC_Receiver1%iLocalBlock_and_iGCtarget_list(GC_Receiver1%nBlock_Pairs,2))
        end do

        ! Now we want to determine the iBlock_pairs and the Block_Pairs_Ptrs.
        ! We will send them from GC_Senders.

        allocate(requests(Tree%nGC_Senders))

        ! First send iBlock_pairs.
        do iGC_sender=1,Tree%nGC_Senders
            GC_Sender1=>Tree%GC_Senders(iGC_sender)
            call MPI_ISEND(GC_Sender1%iBlock_pairs,2*GC_Sender1%nBlock_Pairs,&
            mpi_integer,GC_Sender1%iRank,1,MPI_COMM_WORLD,request,ierr)
            requests(iGC_sender)=request
        end do

        do iGC_receiver=1,Tree%nGC_Receivers
            GC_Receiver1=>Tree%GC_Receivers(iGC_receiver)
            call MPI_RECV(GC_Receiver1%iBlock_pairs,2*GC_Receiver1%nBlock_Pairs,&
            mpi_integer,GC_Receiver1%iRank,1,MPI_COMM_WORLD,status,ierr)
            requests(iGC_receiver)=request
        end do

        call MPI_WAITALL(Tree%nGC_Senders,requests,MPI_STATUSES_IGNORE,ierr)

        ! Then send Block_Pairs_Ptrs.
        do iGC_sender=1,Tree%nGC_Senders
            GC_Sender1=>Tree%GC_Senders(iGC_sender)
            call MPI_ISEND(GC_Sender1%Block_Pairs_Ptrs,2*GC_Sender1%nBlock_Pairs,&
            mpi_integer,GC_Sender1%iRank,1,MPI_COMM_WORLD,request,ierr)
            requests(iGC_sender)=request
        end do

        do iGC_receiver=1,Tree%nGC_Receivers
            GC_Receiver1=>Tree%GC_Receivers(iGC_receiver)
            call MPI_RECV(GC_Receiver1%Block_Pairs_Ptrs,2*GC_Receiver1%nBlock_Pairs,&
            mpi_integer,GC_Receiver1%iRank,1,MPI_COMM_WORLD,status,ierr)
            requests(iGC_receiver)=request
        end do

        call MPI_WAITALL(Tree%nGC_Senders,requests,MPI_STATUSES_IGNORE,ierr)

        deallocate(requests)

        ! Ok. Now one more thing to do is to get the iLocalBlock_and_iGCtarget_list.
        ! This is pretty easy: since we have iBlock pairs, the second column tells
        ! which block pair it is -- which can be used to get the iLocalBlock. Then 
        ! we want to loop all the GC_targets of this local block to get the iGCtarget.

        do iGC_receiver=1,Tree%nGC_Receivers
            GC_Receiver1=>Tree%GC_Receivers(iGC_receiver)
            GC_Receiver1%iLocalBlock_and_iGCtarget_list(:,1)=&
                GC_Receiver1%iBlock_pairs(:,2)-ranges_of_ranks(MpiRank,1)+1
            
            do iofBlock_pair=1,GC_Receiver1%nBlock_Pairs
                iLocalBlock=GC_Receiver1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,1)
                Block_target=>Tree%LocalBlocks(iLocalBlock)
                do iGC_target=1,Block_target%nGC_targets
                    GC_target1=>Block_target%GC_targets(iGC_target)
                    if (GC_target1%iBlock==GC_Receiver1%iBlock_pairs(iofBlock_pair,1)) then
                        GC_Receiver1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,2)=iGC_target
                        exit
                    end if
                end do
            end do
        end do
    end subroutine ModCommunication_SetGC_Receiver

    ! Now let's do HC. It's much more complicated cuz we have MG levels.

    subroutine ModCommunication_SetHC_Sender(Tree)
        implicit none
        type(YYTree),target             ::  Tree
        integer                         ::  iUniqueMpiRank
        integer                         ::  UniqueMpiRank_ptr
        integer                         ::  UniqueMpiRanks(MpiSize,Multigrid_nLevels)

        ! iBlock_pair means different so choose another name
        integer                         ::  iofBlock_pair,iHC                   

        integer                         ::  iLevel
        integer                         ::  iHC_sender
        integer                         ::  iLocalBlock,iHC_source
        type(BlockType),pointer         ::  Block_source
        type(GC_target),pointer         ::  HC_source1
        type(GC_interface),pointer      ::  HC_Sender1

        ! Allocate the nHC_Senders_MGL and HC_Senders_MGL_ptrs.
        allocate(Tree%nHC_Senders_MGL(Multigrid_nLevels))
        allocate(Tree%HC_Senders_MGL_ptrs(Multigrid_nLevels,2))
        Tree%nHC_Senders_MGL=0

        ! Loop all the multigrid levels.
        do iLevel=1,Multigrid_nLevels
            ! Loop all local Blocks%GC_sources to find all neighbour ranks.
            ! First set UniqueMpiRanks to be -1. Loop, if find any new neighbour,
            ! record it there.

            UniqueMpiRank_ptr=1
            UniqueMpiRanks=-1

            do iLocalBlock=1,Tree%nLocalBlocks
                Block_source=>Tree%LocalBlocks(iLocalBlock)
                do iHC_source=1,Block_source%Multigrid_levels(iLevel)%nHC_sources
                    HC_source1=>Block_source%Multigrid_levels(iLevel)%HC_sources(iHC_source)

                    if (HC_source1%iRank/=MpiRank) then
                        ! Loop all recorded unique ranks to see if it's already in there.
                        ! If not, record it.
                        do iUniqueMpiRank=1,UniqueMpiRank_ptr
                            if (UniqueMpiRanks(iUniqueMpiRank,iLevel)==HC_source1%iRank) exit
                        end do
                        if (iUniqueMpiRank==UniqueMpiRank_ptr+1) then
                            UniqueMpiRanks(UniqueMpiRank_ptr,iLevel)=HC_source1%iRank
                            UniqueMpiRank_ptr=UniqueMpiRank_ptr+1
                        end if
                    end if
                end do
            end do

            ! Now we have the UniqueMPIRank-1 as the nHC_Senders for this level.
            ! However, we don't have HC_Senders' size until all levels are done.
            ! The Unique ranks of this level are stored in UniqueMpiRanks(:,iLevel).

            Tree%nHC_Senders_MGL(iLevel)=Tree%nHC_Senders_MGL(iLevel)+UniqueMpiRank_ptr-1
            Tree%HC_Senders_MGL_ptrs(iLevel,1)=Tree%nHC_Senders_MGL(iLevel)-UniqueMpiRank_ptr+1
            Tree%HC_Senders_MGL_ptrs(iLevel,2)=Tree%nHC_Senders_MGL(iLevel)
        end do

        ! Now we have the nHC_Senders_MGL and HC_Senders_MGL_ptrs.
        ! We can allocate the HC_Senders_MGL.

        allocate(Tree%HC_Senders_MGL(sum(Tree%nHC_Senders_MGL)))
        do iLevel=1,Multigrid_nLevels
            Tree%HC_Senders_MGL(&
                Tree%HC_Senders_MGL_ptrs(iLevel,1):Tree%HC_Senders_MGL_ptrs(iLevel,2)&
                )%iRank=UniqueMpiRanks(1:Tree%nHC_Senders_MGL(iLevel),iLevel)
        end do

        ! Now it's time to get the iBlock_pairs and Block_Pairs_Ptrs.
        ! Pretty similar to GC_Senders but just one more loop for levels.

        do iLevel=1,Multigrid_nLevels
            do iHC_sender=1,Tree%nHC_Senders_MGL(iLevel)

                HC_Sender1=>Tree%HC_Senders_MGL(Tree%HC_Senders_MGL_ptrs(iLevel,1)+iHC_sender-1)
                HC_Sender1%nBlock_Pairs=0
                HC_Sender1%nGC=0
    
                do iLocalBlock=1,Tree%nLocalBlocks
                    Block_source=>Tree%LocalBlocks(iLocalBlock)
                    do iHC_source=1,Block_source%Multigrid_levels(iLevel)%nHC_sources
                        HC_source1=>Block_source%Multigrid_levels(iLevel)%HC_sources(iHC_source)
    
                        ! No need to record. Different GC_sources always correspond to
                        ! different block pairs.
    
                        if (HC_source1%iRank/=MpiRank) then
                            HC_Sender1%nBlock_Pairs=HC_Sender1%nBlock_Pairs+1
                            HC_Sender1%nGC=HC_Sender1%nGC+HC_source1%nGC
                        end if
    
                    end do
                end do
    
                allocate(HC_Sender1%message(HC_Sender1%nGC,1))
                allocate(HC_Sender1%iBlock_pairs(HC_Sender1%nBlock_Pairs,2))
                allocate(HC_Sender1%Block_Pairs_Ptrs(HC_Sender1%nBlock_Pairs,2))
                allocate(HC_Sender1%iLocalBlock_and_iGCtarget_list(HC_Sender1%nBlock_Pairs,2))
    
                ! Loop again to find the iBlock_pairs and the Block_Pairs_Ptrs.
                ! We need iGC and iBlock_pair to record the current position.
    
                iofBlock_pair=1
                iHC=1
                do iLocalBlock=1,Tree%nLocalBlocks
                    Block_source=>Tree%LocalBlocks(iLocalBlock)
                    do iHC_source=1,Block_source%Multigrid_levels(iLevel)%nHC_sources
                        HC_source1=>Block_source%Multigrid_levels(iLevel)%HC_sources(iHC_source)
    
                        if (HC_source1%iRank/=MpiRank) then
                            ! Record the iBlock_pairs.
                            HC_Sender1%iBlock_pairs(iofBlock_pair,1)=Block_source%iBlock
                            HC_Sender1%iBlock_pairs(iofBlock_pair,2)=HC_source1%iBlock
                            iofBlock_pair=iofBlock_pair+1
    
                            ! Ptrs
                            HC_Sender1%Block_Pairs_Ptrs(iofBlock_pair,1)=iHC
                            iHC=iHC+HC_source1%nGC
                            HC_Sender1%Block_Pairs_Ptrs(iofBlock_pair,2)=iHC-1
    
                            ! iLocalBlock and iGCtarget
                            HC_Sender1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,1)=iLocalBlock
                            HC_Sender1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,2)=iHC_source
                        end if
                    end do
                end do
            end do 
        end do
    end subroutine ModCommunication_SetHC_Sender

    subroutine ModCommunication_SetHC_Receiver(Tree)
        implicit none
        type(YYTree),target             ::  Tree
        integer                         ::  iUniqueMpiRank
        integer                         ::  UniqueMpiRank_ptr
        integer                         ::  UniqueMpiRanks(MpiSize,Multigrid_nLevels)

        integer                         ::  iofBlock_pair

        integer                         ::  iLevel
        integer                         ::  iLocalBlock,iHC_target,iHC_sender,iHC_receiver
        type(BlockType),pointer         ::  Block_target
        type(GC_target),pointer         ::  HC_target1
        type(GC_interface),pointer      ::  HC_Sender1,HC_Receiver1
        integer                         ::  ierr,request,status(MPI_STATUS_SIZE)
        integer,allocatable             ::  requests(:)

        ! Loop all local Blocks%GC_sources to find all neighbour ranks.
        ! First set UniqueMpiRanks to be -1. Loop, if find any new neighbour,
        ! record it there.

        ! Allocate the nHC_Receivers_MGL and HC_Receivers_MGL_ptrs.
        allocate(Tree%nHC_Receivers_MGL(Multigrid_nLevels))
        allocate(Tree%HC_Receivers_MGL_ptrs(Multigrid_nLevels,2))
        Tree%nHC_Receivers_MGL=0

        ! Loop all the multigrid levels.
        do iLevel=1,Multigrid_nLevels
            ! Loop all local Blocks%GC_sources to find all neighbour ranks.
            ! First set UniqueMpiRanks to be -1. Loop, if find any new neighbour,
            ! record it there.

            UniqueMpiRank_ptr=1
            UniqueMpiRanks=-1

            do iLocalBlock=1,Tree%nLocalBlocks
                Block_target=>Tree%LocalBlocks(iLocalBlock)
                do iHC_target=1,Block_target%Multigrid_levels(iLevel)%nHC_targets
                    HC_target1=>Block_target%Multigrid_levels(iLevel)%HC_targets(iHC_target)

                    if (HC_target1%iRank/=MpiRank) then
                        ! Loop all recorded unique ranks to see if it's already in there.
                        ! If not, record it.
                        do iUniqueMpiRank=1,UniqueMpiRank_ptr
                            if (UniqueMpiRanks(iUniqueMpiRank,iLevel)==HC_target1%iRank) exit
                        end do
                        if (iUniqueMpiRank==UniqueMpiRank_ptr+1) then
                            UniqueMpiRanks(UniqueMpiRank_ptr,iLevel)=HC_target1%iRank
                            UniqueMpiRank_ptr=UniqueMpiRank_ptr+1
                        end if
                    end if
                end do
            end do

            ! Now we have the UniqueMPIRank-1 as the nHC_Senders for this level.
            ! However, we don't have HC_Senders' size until all levels are done.
            ! The Unique ranks of this level are stored in UniqueMpiRanks(:,iLevel).

            Tree%nHC_Receivers_MGL(iLevel)=Tree%nHC_Receivers_MGL(iLevel)+UniqueMpiRank_ptr-1
            Tree%HC_Receivers_MGL_ptrs(iLevel,1)=Tree%nHC_Receivers_MGL(iLevel)-UniqueMpiRank_ptr+1
            Tree%HC_Receivers_MGL_ptrs(iLevel,2)=Tree%nHC_Receivers_MGL(iLevel)
        end do

        ! Now we have the nHC_Receivers_MGL and HC_Receivers_MGL_ptrs.
        ! We can allocate the HC_Receivers_MGL.

        allocate(Tree%HC_Receivers_MGL(sum(Tree%nHC_Receivers_MGL)))
        do iLevel=1,Multigrid_nLevels
            Tree%HC_Receivers_MGL(&
                Tree%HC_Receivers_MGL_ptrs(iLevel,1):Tree%HC_Receivers_MGL_ptrs(iLevel,2)&
                )%iRank=UniqueMpiRanks(1:Tree%nHC_Receivers_MGL(iLevel),iLevel)
        end do

        ! Get the iBlock_pairs for each HC_receiver.
        ! Loop all the local blocks and HC_targets to get the unique pairs.
        ! The process is pretty similar to above, but we need a 2D array 
        ! for recording since it's a pair of blocks.
        ! One more loop for levels.

        do iLevel=1,Multigrid_nLevels
            do iHC_receiver=1,Tree%nHC_Receivers_MGL(iLevel)

                HC_Receiver1=>Tree%HC_Receivers_MGL(Tree%HC_Receivers_MGL_ptrs(iLevel,1)+iHC_receiver-1)
                HC_Receiver1%nBlock_Pairs=0
                HC_Receiver1%nGC=0
    
                do iLocalBlock=1,Tree%nLocalBlocks
                    Block_target=>Tree%LocalBlocks(iLocalBlock)
                    do iHC_target=1,Block_target%Multigrid_levels(iLevel)%nHC_targets
                        HC_target1=>Block_target%Multigrid_levels(iLevel)%HC_targets(iHC_target)
    
                        ! No need to record. Different HC_targets always correspond to
                        ! different block pairs.
    
                        if (HC_target1%iRank/=MpiRank) then
                            HC_Receiver1%nBlock_Pairs=HC_Receiver1%nBlock_Pairs+1
                            HC_Receiver1%nGC=HC_Receiver1%nGC+HC_target1%nGC
                        end if
    
                    end do
                end do
    
                allocate(HC_Receiver1%message(HC_Receiver1%nGC,1))
                allocate(HC_Receiver1%iBlock_pairs(HC_Receiver1%nBlock_Pairs,2))
                allocate(HC_Receiver1%Block_Pairs_Ptrs(HC_Receiver1%nBlock_Pairs,2))
                allocate(HC_Receiver1%iLocalBlock_and_iGCtarget_list(HC_Receiver1%nBlock_Pairs,2))
            end do
    
            ! Now we want to determine the iBlock_pairs and the Block_Pairs_Ptrs.
            ! We will send them from GC_Senders.
    
            allocate(requests(Tree%nHC_Senders_MGL(iLevel)))
    
            ! First send iBlock_pairs.
            do iHC_sender=1,Tree%nHC_Senders_MGL(iLevel)
                HC_Sender1=>Tree%HC_Senders_MGL(Tree%HC_Senders_MGL_ptrs(iLevel,1)+iHC_sender-1)
                call MPI_ISEND(HC_Sender1%iBlock_pairs,2*HC_Sender1%nBlock_Pairs,&
                mpi_integer,HC_Sender1%iRank,1,MPI_COMM_WORLD,request,ierr)
                requests(iHC_sender)=request
            end do
    
            do iHC_receiver=1,Tree%nHC_Receivers_MGL(iLevel)
                HC_Receiver1=>Tree%HC_Receivers_MGL(Tree%HC_Receivers_MGL_ptrs(iLevel,1)+iHC_receiver-1)
                call MPI_RECV(HC_Receiver1%iBlock_pairs,2*HC_Receiver1%nBlock_Pairs,&
                mpi_integer,HC_Receiver1%iRank,1,MPI_COMM_WORLD,status,ierr)
                requests(iHC_receiver)=request
            end do
    
            call MPI_WAITALL(Tree%nHC_Senders_MGL(iLevel),requests,MPI_STATUSES_IGNORE,ierr)
    
            ! Then send Block_Pairs_Ptrs.
            do iHC_sender=1,Tree%nHC_Senders_MGL(iLevel)
                HC_Sender1=>Tree%HC_Senders_MGL(Tree%HC_Senders_MGL_ptrs(iLevel,1)+iHC_sender-1)
                call MPI_ISEND(HC_Sender1%Block_Pairs_Ptrs,2*HC_Sender1%nBlock_Pairs,&
                mpi_integer,HC_Sender1%iRank,1,MPI_COMM_WORLD,request,ierr)
                requests(iHC_sender)=request
            end do
    
            do iHC_receiver=1,Tree%nHC_Receivers_MGL(iLevel)
                HC_Receiver1=>Tree%HC_Receivers_MGL(Tree%HC_Receivers_MGL_ptrs(iLevel,1)+iHC_receiver-1)
                call MPI_RECV(HC_Receiver1%Block_Pairs_Ptrs,2*HC_Receiver1%nBlock_Pairs,&
                mpi_integer,HC_Receiver1%iRank,1,MPI_COMM_WORLD,status,ierr)
                requests(iHC_receiver)=request
            end do
    
            call MPI_WAITALL(Tree%nHC_Senders_MGL(iLevel),requests,MPI_STATUSES_IGNORE,ierr)
    
            deallocate(requests)
    
            ! Ok. Now one more thing to do is to get the iLocalBlock_and_iHCtarget_list.
            ! This is pretty easy: since we have iBlock pairs, the second column tells
            ! which block pair it is -- which can be used to get the iLocalBlock. Then 
            ! we want to loop all the HC_targets of this local block to get the iHCtarget.
    
            do iHC_receiver=1,Tree%nHC_Receivers_MGL(iLevel)
                HC_Receiver1=>Tree%HC_Receivers_MGL(Tree%HC_Receivers_MGL_ptrs(iLevel,1)+iHC_receiver-1)
                HC_Receiver1%iLocalBlock_and_iGCtarget_list(:,1)=&
                    HC_Receiver1%iBlock_pairs(:,2)-ranges_of_ranks(MpiRank,1)+1
                
                do iofBlock_pair=1,HC_Receiver1%nBlock_Pairs
                    iLocalBlock=HC_Receiver1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,1)
                    Block_target=>Tree%LocalBlocks(iLocalBlock)
                    do iHC_target=1,Block_target%Multigrid_levels(iLevel)%nHC_targets
                        HC_target1=>Block_target%Multigrid_levels(iLevel)%HC_targets(iHC_target)
                        if (HC_target1%iBlock==HC_Receiver1%iBlock_pairs(iofBlock_pair,1)) then
                            HC_Receiver1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,2)=iHC_target
                            exit
                        end if
                    end do
                end do
            end do
        end do
    end subroutine ModCommunication_SetHC_Receiver

    ! Ok the new one works. I want a new new version that
    ! calculates EQN_update_R during sendrecv and waitall.
    ! But here is a tricky thing: in the two previous versions,
    ! the local sendrecvs are done BEFORE sendrecv waitall, while
    ! the global sendrecvs are done after those. So here I should 
    ! first move all local sendrecvs to the last loop.
    ! Then the second thing is to add a new input like if_rk to
    ! determine how to calculate EQN_update_R. if_rk should be
    ! renamed as if_rk_communicate, and the new one should be
    ! called if_rk_calc.
    ! And there's one thing to do before all these: define the
    ! EQN_update_R as a variable in BlockType. 

    

    ! In this new version i want to follow the suggestion
    ! from chatgpt that put the mpi_irecv before the mpi_isend.
    ! Hopefully this one will work and be faster.

    subroutine ModCommunication_SendRecvGC_new(Tree,if_rk)

        implicit none
        type(YYTree),target             ::  Tree
        logical,intent(in)              ::  if_rk

        type(BlockType),pointer         ::  Block_source,Block_target
        type(GC_target),pointer         ::  GC_source1,GC_target1
        integer                         ::  iLocalBlock
        integer                         ::  iGC_source,iGC_target
        real,pointer                    ::  primitive_source_IV(:,:,:,:)
        real,pointer                    ::  primitive_target_IV(:,:,:,:)

        integer                         ::  iGC
        integer                         ::  iGC_sender,iGC_receiver
        type(GC_interface),pointer      ::  GC_Sender1,GC_Receiver1
        integer                         ::  iofBlock_pair
        integer                         ::  ierr,status(MPI_STATUS_SIZE)
        integer,allocatable             ::  requests_send(:),requests_recv(:)

        ! Before getting the message in the GC_senders we need to interpolate the
        ! variables into each GC_sources.

        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)

            if (if_rk) then
                primitive_source_IV=>Block_source%primitive_rk_IV
            else
                primitive_source_IV=>Block_source%primitive_IV
            end if

            do iGC_source=1,Block_source%nGC_sources
                GC_source1=>Block_source%GC_sources(iGC_source)
                
                ! Interpolate the primitive variables.
                call ModMath_1D3D_interpolate_1D1D(primitive_source_IV,nvar,ni,nj,nk,ng,&
                    Block_source%xi_I,Block_source%xj_I,Block_source%xk_I,GC_source1%nGC,&
                    GC_source1%xijk_list,GC_source1%primitive_list)
                
                if (GC_source1%if_yin .neqv. Block_source%if_yin) then
                    if (br_>0) GC_source1%primitive_list(:,br_:bp_)=&
                        ModYinYang_VecConv_1D(GC_source1%xijk_list,&
                        GC_source1%primitive_list(:,br_:bp_),GC_source1%nGC)
                    if (vr_>0) GC_source1%primitive_list(:,vr_:vp_)=&
                        ModYinYang_VecConv_1D(GC_source1%xijk_list,&
                       GC_source1%primitive_list(:,vr_:vp_),GC_source1%nGC)
                end if

                if (GC_source1%iRank==MpiRank) then
                    Block_target=>Tree%LocalBlocks(GC_source1%iBlock-ranges_of_ranks(MpiRank,1)+1)

                    if (if_rk) then
                        primitive_target_IV=>Block_target%primitive_rk_IV
                    else
                        primitive_target_IV=>Block_target%primitive_IV
                    end if

                    do iGC=1,GC_source1%nGC
                        primitive_target_IV(GC_source1%ijk_list(iGC,1),&
                            GC_source1%ijk_list(iGC,2),&
                            GC_source1%ijk_list(iGC,3),:)=&
                            GC_source1%primitive_list(iGC,:)
                    end do
                end if
            end do
        end do

        ! Post the Irecv before the Isend.

        allocate(requests_recv(Tree%nGC_Receivers))

        do iGC_receiver=1,Tree%nGC_Receivers
            ! Receive the message.
            GC_Receiver1=>Tree%GC_Receivers(iGC_receiver)
            call MPI_IRECV(GC_Receiver1%message,nvar*GC_Receiver1%nGC,mpi_real,&
            GC_Receiver1%iRank,1,MPI_COMM_WORLD,requests_recv(iGC_receiver),ierr)

            do iofBlock_pair=1,GC_Receiver1%nBlock_Pairs
                GC_target1=>Tree%LocalBlocks(GC_Receiver1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,1))%&
                    GC_targets(GC_Receiver1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,2))
                GC_target1%primitive_list=&
                    GC_Receiver1%message(&
                    GC_Receiver1%Block_Pairs_Ptrs(iofBlock_pair,1):&
                    GC_Receiver1%Block_Pairs_Ptrs(iofBlock_pair,2),:)
            end do
        end do

        ! Now get the message and send them.
        allocate(requests_send(Tree%nGC_Senders))

        do iGC_sender=1,Tree%nGC_Senders
            GC_Sender1=>Tree%GC_Senders(iGC_sender)

            ! First get the message.
            do iofBlock_pair=1,GC_Sender1%nBlock_Pairs
                GC_source1=>Tree%LocalBlocks(GC_Sender1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,1))%&
                    GC_sources(GC_Sender1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,2))
                GC_Sender1%message(&
                    GC_Sender1%Block_Pairs_Ptrs(iofBlock_pair,1):&
                    GC_Sender1%Block_Pairs_Ptrs(iofBlock_pair,2),:)=&
                    GC_source1%primitive_list
            end do

            ! Send the message.
            call MPI_ISEND(GC_Sender1%message,nvar*GC_Sender1%nGC,mpi_real,&
                GC_Sender1%iRank,1,MPI_COMM_WORLD,requests_send(iGC_sender),ierr)
        end do

        call MPI_WAITALL(Tree%nGC_Receivers,requests_recv,MPI_STATUSES_IGNORE,ierr)
        call MPI_WAITALL(Tree%nGC_Senders,requests_send,MPI_STATUSES_IGNORE,ierr)

        ! Assign the GC_targets' primitive_list to primitive blocks.
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_target=>Tree%LocalBlocks(iLocalBlock)

            if (if_rk) then
                primitive_target_IV=>Block_target%primitive_rk_IV
            else
                primitive_target_IV=>Block_target%primitive_IV
            end if

            do iGC_target=1,Block_target%nGC_targets
                GC_target1=>Block_target%GC_targets(iGC_target)

                if (GC_target1%iRank/=MpiRank) then
                    do iGC=1,GC_target1%nGC
                        primitive_target_IV(GC_target1%ijk_list(iGC,1),&
                            GC_target1%ijk_list(iGC,2),&
                            GC_target1%ijk_list(iGC,3),:)=&
                            GC_target1%primitive_list(iGC,:)
                    end do
                end if
            end do
        end do

        deallocate(requests_send)
        deallocate(requests_recv)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    end subroutine ModCommunication_SendRecvGC_new

    ! With all above we can send and receive the messages
    ! using G(H)C_Senders and G(H)C_Receivers.
    ! This one for GC. GC has no MG levels.

    subroutine ModCommunication_SendRecvGC(Tree,if_rk)
        implicit none
        type(YYTree),target             ::  Tree
        logical,intent(in)              ::  if_rk

        type(BlockType),pointer         ::  Block_source,Block_target
        type(GC_target),pointer         ::  GC_source1,GC_target1
        integer                         ::  iLocalBlock
        integer                         ::  iGC_source,iGC_target
        real,pointer                    ::  primitive_source_IV(:,:,:,:)
        real,pointer                    ::  primitive_target_IV(:,:,:,:)

        integer                         ::  iGC
        integer                         ::  iGC_sender,iGC_receiver
        type(GC_interface),pointer      ::  GC_Sender1,GC_Receiver1
        integer                         ::  iofBlock_pair
        integer                         ::  ierr,request,status(MPI_STATUS_SIZE)
        integer,allocatable             ::  requests(:)

        ! Before getting the message in the GC_senders we need to interpolate the
        ! variables into each GC_sources.

        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)

            if (if_rk) then
                primitive_source_IV=>Block_source%primitive_rk_IV
            else
                primitive_source_IV=>Block_source%primitive_IV
            end if

            do iGC_source=1,Block_source%nGC_sources
                GC_source1=>Block_source%GC_sources(iGC_source)
                
                ! Interpolate the primitive variables.
                call ModMath_1D3D_interpolate_1D1D(primitive_source_IV,nvar,ni,nj,nk,ng,&
                    Block_source%xi_I,Block_source%xj_I,Block_source%xk_I,GC_source1%nGC,&
                    GC_source1%xijk_list,GC_source1%primitive_list)
                
                if (GC_source1%if_yin .neqv. Block_source%if_yin) then
                    if (br_>0) GC_source1%primitive_list(:,br_:bp_)=&
                        ModYinYang_VecConv_1D(GC_source1%xijk_list,&
                        GC_source1%primitive_list(:,br_:bp_),GC_source1%nGC)
                    if (vr_>0) GC_source1%primitive_list(:,vr_:vp_)=&
                        ModYinYang_VecConv_1D(GC_source1%xijk_list,&
                       GC_source1%primitive_list(:,vr_:vp_),GC_source1%nGC)
                end if

                if (GC_source1%iRank==MpiRank) then
                    Block_target=>Tree%LocalBlocks(GC_source1%iBlock-ranges_of_ranks(MpiRank,1)+1)

                    if (if_rk) then
                        primitive_target_IV=>Block_target%primitive_rk_IV
                    else
                        primitive_target_IV=>Block_target%primitive_IV
                    end if

                    do iGC=1,GC_source1%nGC
                        primitive_target_IV(GC_source1%ijk_list(iGC,1),&
                            GC_source1%ijk_list(iGC,2),&
                            GC_source1%ijk_list(iGC,3),:)=&
                            GC_source1%primitive_list(iGC,:)
                    end do
                end if
            end do
        end do

        ! Now get the message and send them.
        allocate(requests(Tree%nGC_Senders))

        do iGC_sender=1,Tree%nGC_Senders
            GC_Sender1=>Tree%GC_Senders(iGC_sender)

            ! First get the message.
            do iofBlock_pair=1,GC_Sender1%nBlock_Pairs
                GC_source1=>Tree%LocalBlocks(GC_Sender1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,1))%&
                    GC_sources(GC_Sender1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,2))
                GC_Sender1%message(&
                    GC_Sender1%Block_Pairs_Ptrs(iofBlock_pair,1):&
                    GC_Sender1%Block_Pairs_Ptrs(iofBlock_pair,2),:)=&
                    GC_source1%primitive_list
            end do

            ! Send the message.
            call MPI_ISEND(GC_Sender1%message,nvar*GC_Sender1%nGC,mpi_real,&
            GC_Sender1%iRank,1,MPI_COMM_WORLD,request,ierr)
            requests(iGC_sender)=request
        end do

        do iGC_receiver=1,Tree%nGC_Receivers
            ! Receive the message.
            GC_Receiver1=>Tree%GC_Receivers(iGC_receiver)
            call MPI_RECV(GC_Receiver1%message,nvar*GC_Receiver1%nGC,mpi_real,&
                GC_Receiver1%iRank,1,MPI_COMM_WORLD,status,ierr)

            do iofBlock_pair=1,GC_Receiver1%nBlock_Pairs
                GC_target1=>Tree%LocalBlocks(GC_Receiver1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,1))%&
                    GC_targets(GC_Receiver1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,2))
                GC_target1%primitive_list=&
                    GC_Receiver1%message(&
                    GC_Receiver1%Block_Pairs_Ptrs(iofBlock_pair,1):&
                    GC_Receiver1%Block_Pairs_Ptrs(iofBlock_pair,2),:)
            end do
        end do

        call MPI_WAITALL(Tree%nGC_Senders,requests,MPI_STATUSES_IGNORE,ierr)

        ! Assign the GC_targets' primitive_list to primitive blocks.
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_target=>Tree%LocalBlocks(iLocalBlock)

            if (if_rk) then
                primitive_target_IV=>Block_target%primitive_rk_IV
            else
                primitive_target_IV=>Block_target%primitive_IV
            end if

            do iGC_target=1,Block_target%nGC_targets
                GC_target1=>Block_target%GC_targets(iGC_target)

                if (GC_target1%iRank/=MpiRank) then
                    do iGC=1,GC_target1%nGC
                        primitive_target_IV(GC_target1%ijk_list(iGC,1),&
                            GC_target1%ijk_list(iGC,2),&
                            GC_target1%ijk_list(iGC,3),:)=&
                            GC_target1%primitive_list(iGC,:)
                    end do
                end if
            end do
        end do

        deallocate(requests)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
    end subroutine ModCommunication_SendRecvGC

    subroutine ModCommunication_SendRecvHC(Tree,iLevel,i_option)
        implicit none
        type(YYTree),target             ::  Tree
        integer,intent(in)              ::  iLevel
        integer,intent(in)              ::  i_option

        type(Multigrid_level),pointer   ::  MG_source,MG_target
        type(GC_target),pointer         ::  HC_source1,HC_target1
        integer                         ::  iLocalBlock
        integer                         ::  iHC_source,iHC_target
        real,pointer                    ::  primitive_source_III(:,:,:)
        real,pointer                    ::  primitive_target_III(:,:,:)

        integer                         ::  iHC
        integer                         ::  iHC_sender,iHC_receiver
        type(GC_interface),pointer      ::  HC_Sender1,HC_Receiver1
        integer                         ::  iofBlock_pair
        integer                         ::  ierr,request,status(MPI_STATUS_SIZE)
        integer,allocatable             ::  requests(:)

        ! Before getting the message in the HC_senders we need to interpolate the
        ! variables into each HC_sources.

        do iLocalBlock=1,Tree%nLocalBlocks
            MG_source=>Tree%LocalBlocks(iLocalBlock)%Multigrid_levels(iLevel)
            
            select case (i_option)
            case (1)
                primitive_source_III=>MG_source%p_III
            case (2)
                primitive_source_III=>MG_source%x_III
            case (3)
                primitive_source_III=>MG_source%z_III
            case default
                primitive_source_III=>MG_source%p_III
            end select

            do iHC_source=1,MG_source%nHC_sources
                HC_source1=>MG_source%HC_sources(iHC_source)

                ! Interpolate the primitive variables.
                call ModMath_3D_interpolate_1D(primitive_source_III,&
                    MG_source%mi,MG_source%mj,MG_source%mk,1,&
                    MG_source%xi_I,MG_source%xj_I,MG_source%xk_I,HC_source1%nGC,&
                    HC_source1%xijk_list,HC_source1%primitive_list)

                if (HC_source1%iRank==MpiRank) then
                    MG_target=>Tree%LocalBlocks(HC_source1%iBlock-ranges_of_ranks(MpiRank,1)+1)%Multigrid_levels(iLevel)

                    select case (i_option)
                    case (1)
                        primitive_target_III=>MG_target%p_III
                    case (2)
                        primitive_target_III=>MG_target%x_III
                    case (3)
                        primitive_target_III=>MG_target%z_III
                    case default
                        primitive_target_III=>MG_target%p_III
                    end select

                    do iHC=1,HC_source1%nGC
                        primitive_target_III(HC_source1%ijk_list(iHC,1),&
                            HC_source1%ijk_list(iHC,2),&
                            HC_source1%ijk_list(iHC,3))=&
                            HC_source1%primitive_list(iHC,1)
                    end do
                end if
            end do
        end do

        ! Now get the message and send them.
        allocate(requests(Tree%nHC_Senders_MGL(iLevel)))

        do iHC_sender=1,Tree%nHC_Senders_MGL(iLevel)
            HC_Sender1=>Tree%HC_Senders_MGL(Tree%HC_Senders_MGL_ptrs(iLevel,1)+iHC_sender-1)

            ! First get the message.
            do iofBlock_pair=1,HC_Sender1%nBlock_Pairs
                HC_source1=>Tree%LocalBlocks(HC_Sender1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,1))%&
                    Multigrid_levels(iLevel)%HC_sources(HC_Sender1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,2))
                HC_Sender1%message(&
                    HC_Sender1%Block_Pairs_Ptrs(iofBlock_pair,1):&
                    HC_Sender1%Block_Pairs_Ptrs(iofBlock_pair,2),:)=&
                    HC_source1%primitive_list
            end do

            ! Send the message.
            call MPI_ISEND(HC_Sender1%message,HC_Sender1%nGC,mpi_real,&
            HC_Sender1%iRank,1,MPI_COMM_WORLD,request,ierr)
            requests(iHC_sender)=request
        end do

        do iHC_receiver=1,Tree%nHC_Receivers_MGL(iLevel)
            ! Receive the message.
            HC_Receiver1=>Tree%HC_Receivers_MGL(Tree%HC_Receivers_MGL_ptrs(iLevel,1)+iHC_receiver-1)
            call MPI_RECV(HC_Receiver1%message,HC_Receiver1%nGC,mpi_real,HC_Receiver1%iRank,1,MPI_COMM_WORLD,status,ierr)

            do iofBlock_pair=1,HC_Receiver1%nBlock_Pairs
                HC_target1=>Tree%LocalBlocks(HC_Receiver1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,1))%&
                    Multigrid_levels(iLevel)%HC_targets(HC_Receiver1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,2))
                HC_target1%primitive_list=&
                    HC_Receiver1%message(&
                    HC_Receiver1%Block_Pairs_Ptrs(iofBlock_pair,1):&
                    HC_Receiver1%Block_Pairs_Ptrs(iofBlock_pair,2),:)
            end do
        end do

        call MPI_WAITALL(Tree%nHC_Senders_MGL(iLevel),requests,MPI_STATUSES_IGNORE,ierr)

        ! Assign the GC_targets' primitive_list to primitive blocks.

        do iLocalBlock=1,Tree%nLocalBlocks
            MG_target=>Tree%LocalBlocks(iLocalBlock)%Multigrid_levels(iLevel)

            select case (i_option)
            case (1)
                primitive_target_III=>MG_target%p_III
            case (2)
                primitive_target_III=>MG_target%x_III
            case (3)
                primitive_target_III=>MG_target%z_III
            case default
                primitive_target_III=>MG_target%p_III
            end select

            do iHC_target=1,MG_target%nHC_targets
                HC_target1=>MG_target%HC_targets(iHC_target)

                if (HC_target1%iRank/=MpiRank) then
                    do iHC=1,HC_target1%nGC
                        primitive_target_III(HC_target1%ijk_list(iHC,1),&
                            HC_target1%ijk_list(iHC,2),&
                            HC_target1%ijk_list(iHC,3))=&
                            HC_target1%primitive_list(iHC,1)
                    end do
                end if
            end do
        end do

        deallocate(requests)

        
    end subroutine ModCommunication_SendRecvHC

    ! About the tags.

    function ModCommunication_GetTag(iBlockSend,iBlockRecv,iRankSend,iRankRecv) result(Tag)
        implicit none
        integer,intent(in)      :: iBlockSend,iBlockRecv
        integer,intent(in)      :: iRankSend,iRankRecv 
        
        integer                 :: tag

        Tag=nMaxBlocksPerRank*(iBlockSend-ranges_of_ranks(iRankSend,1))+&
            (iBlockRecv-ranges_of_ranks(iRankRecv,1))
    end function ModCommunication_GetTag

    function ModCommunication_SolveTag(Tag,iRankSend,iRankRecv) result(iBlockSendRecv)
        implicit none
        integer,intent(in)      :: Tag
        integer,intent(in)      :: iRankRecv,iRankSend

        integer                 :: iBlockSendRecv(2)

        iBlockSendRecv(1)=Tag/nMaxBlocksPerRank+ranges_of_ranks(iRankSend,1)
        iBlockSendRecv(2)=mod(Tag,nMaxBlocksPerRank)+ranges_of_ranks(iRankRecv,1)
    end function ModCommunication_SolveTag
end module