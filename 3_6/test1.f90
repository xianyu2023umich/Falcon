program test1

    use ModBlock,           only: BlockType
    use ModParameters,      only:   MpiSize, MpiRank, r_range, nSteps, DoCheck,ni
    use ModReadParameters,  only:   ModReadParameters_read
    use ModYinYangTree,     only:   YYTree, YinYangTree_InitTree, YinYangTree_SetAll
    use ModCommunication,   only:   ModCommunication_SetGCAll, &
                                    ModCommunication_SendRecvGC_new
    use ModSavePlot,        only:   ModSave_DoAll
    use ModCheck,           only:   ModCheck_primitive
    use ModAdvance,         only:   ModAdvance_rk4
    use ModAMR,             only:   ModAMR_set_grid
    use ModInitiation,      only:   ModInitiation_rand_velocity
    use MPI

    implicit none

    type(YYTree),target     ::  Tree
    integer                 ::  ierr
    integer                 ::  iStep
    real(8)                 ::  dt
    integer :: iblock
    type(BlockType),pointer :: Block1

    call MPI_INIT(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, MpiSize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, MpiRank, ierr)

    call ModReadParameters_read('PARAM.in.orig',1)
    if (MpiRank==0) print *, 'read parameters'

    call test1_INITIATION(Tree)
    if (MpiRank==0) print *, 'initialized tree and state'

    if (DoCheck) call ModCheck_primitive(Tree,.false.)
    if (DoCheck) call ModCheck_primitive(Tree,.true.)

    call ModCommunication_SendRecvGC_new(Tree,.false.)
    call ModCommunication_SendRecvGC_new(Tree,.false.)
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if (MpiRank==0) print *,'Complete initial communication.'
    
    if (DoCheck) call ModCheck_primitive(Tree,.false.)
    if (DoCheck) call ModCheck_primitive(Tree,.true.)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if (MpiRank==0) print *,'Complete initial checks.'

    do iblock=1,2
        Block1 => Tree%LocalBlocks(iblock)
        print *,iBlock,Block1%gamma3_minus_1_III(1:ni,1,1)
    end do

    do iStep=1,nSteps
        call ModAdvance_rk4(Tree,.true.,dt)
        if (MpiRank==0) print *, 'Complete Advancing at iStep=',iStep,'dt=',dt

        if (DoCheck) call ModCheck_primitive(Tree,.false.)
        call ModSave_DoAll(Tree,iStep)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    end do

    !call ModCommunication_SendRecvGC_new(Tree,.false.)
    if (MpiRank==0) print *, 'refreshed ghost cells before final save'

    call ModSave_DoAll(Tree,nSteps)
    if (MpiRank==0) print *, 'saved final cube'

    call MPI_FINALIZE(ierr)

contains

    subroutine test1_INITIATION(Tree)
        type(YYTree),intent(inout),target :: Tree

        call YinYangTree_InitTree(Tree, r_range)
        if (MpiRank==0) print *, '  YinYangTree_InitTree'

        call ModAMR_set_grid(Tree)
        if (MpiRank==0) print *, '  ModAMR_set_grid'

        call YinYangTree_SetAll(Tree)
        if (MpiRank==0) print *, '  YinYangTree_SetAll'

        call ModCommunication_SetGCAll(Tree)
        if (MpiRank==0) print *, '  ModCommunication_SetGCAll'

        if (MpiRank==0 .or. MpiRank==MpiSize-1) then
            print *,'Tree%NumLeafNodes=',Tree%NumLeafNodes,'MpiRank=',MpiRank,' has ',Tree%nLocalBlocks,' blocks'
        end if

        call ModInitiation_rand_velocity(Tree, 1000.0d0)
        if (MpiRank==0) print *, '  ModInitiation_rand_velocity (v_rms=100)'
    end subroutine test1_INITIATION

end program test1
