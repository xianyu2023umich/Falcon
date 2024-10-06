Module ModCommunication

    use ModOcTree
    use ModTimeStep
    use ModMath,only : ModMath_1D3D_interpolate_1D1D
    use ModControl,only : nMaxBlocksPerRank
    use MPI

    contains

    !function ModCommunication_CountCells
    !    implicit none


    !end function ModCommunication_CountCells

    ! use Mpi_reduce to determine the global time step
    !
    subroutine ModCommunication_GlobalTimeStep(dt_local,dt_global)
        implicit none
        real,intent(in)         :: dt_local
        real,intent(out)        :: dt_global

        integer                 :: ierr

        call MPI_AllReduce(dt_local,dt_global,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ierr)
    end subroutine ModCommunication_GlobalTimeStep

    subroutine ModCommunication_SendRecvGCAll(Tree,rk_index,MpiSize,MpiRank)
        implicit none
        type(OcTree),intent(in),target  ::  Tree
        type(Block),pointer             ::  Block1
        integer,intent(in)              ::  rk_index
        integer,intent(in)              ::  MpiSize,MpiRank

        integer                         ::  iLocalBlock
        integer                         ::  iGC_target,iGC_source

        integer                         ::  iGC
        type(GC_target),pointer         ::  GC_source1,GC_target1
        real,pointer                    ::  primitive(:,:,:,:)
        real,allocatable                ::  recv_message(:,:)

        integer                         ::  tag,ierr,request,status(MPI_STATUS_SIZE)

        ! Send all the GC_sources of all the 
        ! local blocks in the tree.
        !
        do iLocalBlock=1,Tree%nLocalBlocks
            Block1 => Tree%LocalBlocks(iLocalBlock)

            ! first, select which rk primitive block
            !
            select case(rk_index)
            case(1)
                primitive => Block1%primitive
            case(2)
                primitive => Block1%primitive_k2
            case(3)
                primitive => Block1%primitive_k3
            case(4)
                primitive => Block1%primitive_k4
            end select

            ! Then, loop all the GC_sources.
            ! Interpolate the primitive variables and mpi_send them.
            !
            do iGC_source=1,size(Block1%GC_sources)
                
                ! assign GC_source pointer
                GC_source1 => Block1%GC_sources(iGC_source)

                !if(MpiRank==1 .and. GC_source1%iRank==7) then
                    !print *,GC_source1%primitive_list(:,53)
                    !print *,GC_source1%xijk_list(53,:),Block1%xi(1),Block1%xj(1),Block1%xk(3)
                    !print *,primitive(:,2,2,4)
                !end if

                ! interpolate the primitive grid to the GCs
                !
                call ModMath_1D3D_interpolate_1D1D(primitive,Block1%nvar,&
                    Block1%ni,Block1%nj,Block1%nk,Block1%ng,&
                    Block1%xi,Block1%xj,Block1%xk,GC_source1%nGC,GC_source1%xijk_list,GC_source1%primitive_list)
            
                ! Get the tag, which contains information for both the send block and receive block.
                !
                tag = ModCommunication_GetTag(MpiSize,Tree%iLeafNode_ranges,Block1%iBlock,GC_source1%iBlock,GC_source1%iRank)

                ! Send GCs.
                if (GC_source1%iRank/=MpiRank) call MPI_ISEND(GC_source1%primitive_list,&
                    Block1%nvar*GC_source1%nGC,mpi_real,GC_source1%iRank,tag,MPI_COMM_WORLD,request,ierr)
                
                !if(MpiRank==1 .and. GC_source1%iRank==7) then
                    !print *,GC_source1%primitive_list(:,53)
                    !print *,GC_source1%xijk_list(53,:),Block1%xi(1),Block1%xj(1),Block1%xk(3)
                    !print *,primitive(:,2,2,4)
                !end if
                
                !if (MpiRank==7 .and. GC_source1%iRank==0) print *,MpiRank,primitive(:,4,4,1)!GC_source1%primitive_list(4,:)
            end do
        end do


        do iLocalBlock=1,Tree%nLocalBlocks
            Block1 => Tree%LocalBlocks(iLocalBlock)

            select case(rk_index)
            case(1)
                primitive => Block1%primitive
            case(2)
                primitive => Block1%primitive_k2
            case(3)
                primitive => Block1%primitive_k3
            case(4)
                primitive => Block1%primitive_k4
            end select
    
            do iGC_target=1,size(Block1%GC_targets)
                if (Block1%GC_targets(iGC_target)%iRank.ne.MpiRank) then
    
                    GC_target1 => Block1%GC_targets(iGC_target)
    
                    allocate(recv_message(Block1%nvar,GC_target1%nGC))
                    !allocate(recv_message_1(GC_target1%nGC))
    
                    tag=ModCommunication_GetTag(MpiSize,Tree%iLeafNode_ranges,GC_target1%iBlock,Block1%iBlock,MpiRank)

                    call MPI_RECV(recv_message,Block1%nvar*GC_target1%nGC,mpi_real,GC_target1%iRank,tag,MPI_COMM_WORLD,status,ierr)
    
                    do iGC = 1,GC_target1%nGC
                        primitive(:,GC_target1%ijk_list(iGC,1),&
                                    GC_target1%ijk_list(iGC,2),&
                                    GC_target1%ijk_list(iGC,3))=recv_message(:,iGC)
                    end do
                    !if (MpiRank==0 .and. GC_target1%iRank==7) print *,MpiRank,recv_message(4,:)
                    deallocate(recv_message)
                end if
            end do
        end do

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    end subroutine ModCommunication_SendRecvGCAll

    ! Send all the GC_sources of all the 
    ! local blocks in the tree.
    !
    subroutine ModCommunication_SendGCAll(Tree,rk_index,MpiSize,MpiRank)
        implicit none
        type(OcTree),intent(in),target  ::  Tree
        type(Block),pointer             ::  Block1
        integer,intent(in)              ::  rk_index
        integer,intent(in)              ::  MpiSize,MpiRank

        integer                         ::  iLocalBlock

        do iLocalBlock=1,Tree%nLocalBlocks
            Block1 => Tree%LocalBlocks(iLocalBlock)

            call ModCommunication_SendGC_Single(Tree,Block1,rk_index,MpiSize,MpiRank)
        end do

    end subroutine ModCommunication_SendGCAll

    ! Send the GC_sources of one single block
    ! 
    subroutine ModCommunication_SendGC_Single(Tree,Block1,rk_index,MpiSize,MpiRank)
        implicit none
        type(OcTree),intent(in)         ::  Tree
        type(Block), intent(in),target  ::  Block1
        integer,     intent(in)         ::  rk_index
        integer,     intent(in)         ::  MpiSize,MpiRank

        integer                         ::  iGC_source
        type(GC_target),pointer         ::  GC_source1
        real,pointer                    ::  primitive(:,:,:,:)

        integer                         ::  tag,ierr,request

        ! first, select which rk primitive block
        !
        select case(rk_index)
        case(1)
            primitive => Block1%primitive
        case(2)
            primitive => Block1%primitive_k2
        case(3)
            primitive => Block1%primitive_k3
        case(4)
            primitive => Block1%primitive_k4
        end select

        ! Then, loop all the GC_sources.
        ! Interpolate the primitive variables and mpi_send them.
        !
        do iGC_source=1,size(Block1%GC_sources)
            
            GC_source1 => Block1%GC_sources(iGC_source)

            ! interpolate the primitive grid to the GCs
            !
            call ModMath_1D3D_interpolate_1D1D(primitive,Block1%nvar,&
                Block1%ni,Block1%nj,Block1%nk,Block1%ng,&
                Block1%xi,Block1%xj,Block1%xk,GC_source1%nGC,GC_source1%xijk_list,GC_source1%primitive_list)
            
            ! Get the tag, which contains information for both the send block and receive block.
            !
            tag = ModCommunication_GetTag(MpiSize,Tree%iLeafNode_ranges,Block1%iBlock,GC_source1%iBlock,GC_source1%iRank)

            ! Send GCs.
            if (GC_source1%iRank/=MpiRank) call MPI_ISEND(GC_source1%primitive_list,&
                    Block1%nvar*GC_source1%nGC,mpi_real,GC_source1%iRank,tag,MPI_COMM_WORLD,request,ierr)
        end do
    end subroutine ModCommunication_SendGC_Single

    subroutine ModCommunication_RecvGCAll(Tree,rk_index,MpiSize,MpiRank)
        implicit none
        type(OcTree),intent(in),target  :: Tree
        integer,intent(in)              :: rk_index
        integer,intent(in)              :: MpiSize,MpiRank

        type(Block),pointer             :: Block1
        integer                         :: iLocalBlock

        do iLocalBlock=1,Tree%nLocalBlocks
            Block1 => Tree%LocalBlocks(iLocalBlock)

            call ModCommunication_RecvGC_Single(Tree,Block1,rk_index,MpiSize,MpiRank)
        end do
    end subroutine ModCommunication_RecvGCAll

    subroutine ModCommunication_RecvGC_Single(Tree,Block1,rk_index,MpiSize,MpiRank)
        implicit none
        type(OcTree),intent(in)         :: Tree
        type(Block),intent(in),target   :: Block1
        integer,intent(in)              :: rk_index
        integer,intent(in)              :: MpiSize,MpiRank

        integer                         :: iGC
        type(GC_target),pointer         :: GC_target1
        integer                         :: iGC_target

        real,pointer                    :: primitive(:,:,:,:)
        real,allocatable                :: recv_message(:,:)
        integer                         :: status(MPI_STATUS_SIZE),tag,ierr

        ! first, select which rk primitive block
        !
        select case(rk_index)
        case(1)
            primitive => Block1%primitive
        case(2)
            primitive => Block1%primitive_k2
        case(3)
            primitive => Block1%primitive_k3
        case(4)
            primitive => Block1%primitive_k4
        end select

        ! Then, loop all the GC_targets.
        ! Receive the messages and assign them to primitives
        !
        do iGC_target=1,size(Block1%GC_targets)
            if (Block1%GC_targets(iGC_target)%iRank.ne.MpiRank) then

                ! assign GC_target1 pointer
                GC_target1 => Block1%GC_targets(iGC_target)

                ! allocate the receive buffer
                allocate(recv_message(Block1%nvar,GC_target1%nGC))

                ! get the tag
                tag=ModCommunication_GetTag(MpiSize,Tree%iLeafNode_ranges,GC_target1%iBlock,Block1%iBlock,MpiRank)

                ! receive the message
                call MPI_RECV(recv_message,Block1%nvar*GC_target1%nGC,mpi_real,GC_target1%iRank,tag,MPI_COMM_WORLD,status,ierr)

                ! assign received GC to primitives
                do iGC = 1,GC_target1%nGC
                    primitive(:,GC_target1%ijk_list(iGC,1),&
                                GC_target1%ijk_list(iGC,2),&
                                GC_target1%ijk_list(iGC,3))=recv_message(:,iGC)
                end do
                deallocate(recv_message)
            end if
        end do

    end subroutine ModCommunication_RecvGC_Single

    function ModCommunication_GetTag(MpiSize,iLeafNode_ranges,iBlockSend,iBlockRecv,iRankRecv) result(Tag)
        implicit none
        integer,intent(in)      :: MpiSize
        integer,intent(in)      :: iLeafNode_ranges(0:MpiSize-1,2)
        integer,intent(in)      :: iBlockSend,iBlockRecv
        integer,intent(in)      :: iRankRecv
        
        integer                 :: tag

        Tag=nMaxBlocksPerRank*iBlockSend+(iBlockRecv-iLeafNode_ranges(iRankRecv,1))
    end function ModCommunication_GetTag

    function ModCommunication_SolveTag(MpiSize,iLeafNode_ranges,Tag,iRankRecv) result(iBlockSendRecv)
        implicit none
        integer,intent(in)      :: MpiSize
        integer,intent(in)      :: iLeafNode_ranges(0:MpiSize-1,2)
        integer,intent(in)      :: Tag
        integer,intent(in)      :: iRankRecv

        integer                 :: iBlockSendRecv(2)

        iBlockSendRecv(1)=Tag/nMaxBlocksPerRank
        iBlockSendRecv(2)=mod(Tag,nMaxBlocksPerRank)+iLeafNode_ranges(iRankRecv,1)
    end function ModCommunication_SolveTag
end Module ModCommunication