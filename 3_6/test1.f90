program test1

    use ModBlock,           only:   BlockType
    use ModControl,         only:   t,dt,iStep,iStatus,if_param_file_opened
    use ModParameters,      only:   MpiSize, MpiRank, r_range, nSteps, DoCheck,&
                                    if_do_echo, nStepsEcho
    use ModReadParameters,  only:   ModReadParameters_read
    use ModYinYangTree,     only:   YYTree,&
                                    YinYangTree_InitTree,&
                                    YinYangTree_SetAll,&
                                    YinYangTree_UpdateAll
    use ModCommunication,   only:   ModCommunication_SetGCAll, &
                                    ModCommunication_SendRecvGC_new
    use ModSavePlot,        only:   ModSave_DoAll
    use ModCheck,           only:   ModCheck_primitive
    use ModAdvance,         only:   ModAdvance_rk4
    use ModAMR,             only:   ModAMR_set_grid
    use ModInitiation,      only:   ModInitiation_DoAll
    use MPI

    implicit none

    type(YYTree),target         ::  Tree
    integer                     ::  ierr
    character(len=*), parameter ::  param_file = 'PARAM.in'
    integer                     ::  Logical_unit_param_file = 42
    !integer                     :: iblock
    !type(BlockType),pointer     :: Block1

    iStatus=0
    call MPI_INIT(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, MpiSize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, MpiRank, ierr)

    ! First time to read parameters. This will read all parameters 
    ! until the checkpoint section. If there is no checkpoint section, 
    ! it will read all parameters.
    call ModReadParameters_read(param_file, Logical_unit_param_file)
    if (MpiRank==0) print *, 'read parameters'

    call test1_INITIATION(Tree)
    if (MpiRank==0) print *, 'initialized tree and state'
    if (DoCheck) call ModCheck_primitive(Tree,.true.)


    ! The main loop. Now I don't use iStep=1,nSteps
    ! since there might be checkpoints.
    iStatus=1

    ! Initialize iStep to 1 and t to 0.0
    iStep=1
    t=0.0

    ! The main loop.
    do
        ! Advance
        call ModAdvance_rk4(Tree,.true.,dt)

        ! Check
        if (DoCheck) call ModCheck_primitive(Tree,.false.)

        ! Echo
        if (if_do_echo .and. mod(iStep,nStepsEcho) == 0 .and. MpiRank == 0) then
            print *, 'Complete Advancing at iStep=',iStep,'dt=',dt
        end if

        ! Save
        call ModSave_DoAll(Tree,iStep)

        ! See if we have reached the total nSteps.
        ! If yes then read param again. 

        if (iStep >= nSteps) then
            ! Print checkpoint info.
            if (MpiRank==0) then
                print *,'Advancing stops at iStep=',iStep
                print *,'Attempting to read parameters for the next run...'
            end if

            ! If it is still opened. If yes then read it.
            ! If no then there should be no change to nStep,
            ! so the run will naturally stop.
            if (if_param_file_opened) then
                call ModReadParameters_read(param_file, Logical_unit_param_file)
                call YinYangTree_UpdateAll(Tree,1)
            end if
            
            ! If we have read parameters successfully, then it means we have more steps to do. 
            ! If not, it means we have reached the end of the parameter file and we can stop.
            if (iStep >= nSteps) then
                if (MpiRank==0) print *, 'No more steps to do. Stopping.'
                exit
            else
                if (MpiRank==0) print *, 'Successfully read parameters for the next run. Continuing.'
            end if
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        ! Update iStep and t. 
        iStep = iStep + 1
        t = t + dt
    end do

    call MPI_FINALIZE(ierr)
    iStatus=3

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

        call ModInitiation_DoAll(Tree)
        if (MpiRank==0) print *, '  ModInitiation_DoAll'

        call ModCommunication_SendRecvGC_new(Tree,.false.)
        call ModCommunication_SendRecvGC_new(Tree,.false.)
    
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (MpiRank==0) print *,'Complete initial communication.'
    end subroutine test1_INITIATION

end program test1
