program test1

    use ModBlock,           only:   BlockType
    use ModParameters,      only:   MpiSize, MpiRank, r_range,&
                                    CFL,nSteps,DoCheck
    use ModReadParameters,  only:   ModReadParameters_read
    use ModYinYangTree,     only:   YYTree, YinYangTree_InitTree, YinYangTree_SetAll, YinYangTree_DivideAll
    use ModCommunication,   only:   ModCommunication_SetGCAll
    use ModInitiation,      only:   ModInitiation_DoAll
    use ModSavePlot,        only:   ModSave_DoAll
    use ModCheck,           only:   ModCheck_primitive
    use ModAdvance,         only:   ModAdvance_rk4
    use ModAMR,             only:   ModAMR_set_grid
    use MPI

    implicit none

    type(YYTree),target     ::  Tree                        ! The Tree
    integer                 ::  ierr                        ! For MPI
    integer                 ::  iStep                       ! For step
    character(len=8)        ::  iStep_char
    real                    ::  dt                          ! Time step
    integer                 ::  i

    ! Initiate MPI and get MpiSize/Rank
    call MPI_INIT(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, MpiSize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, MpiRank, ierr)

    ! Read parameters. Do initiation
    call ModReadParameters_read('PARAM.in',1)
    call test1_INITIATION

    ! Main loop
    do iStep=1,nSteps

        ! Advance
        call ModAdvance_rk4(Tree,CFL,.True.,dt)
        if(MpiRank==0)print *,'Complete Advancing at iStep=',iStep,'dt=',dt

        ! Check if there is NaN
        if (DoCheck) call ModCheck_primitive(Tree,.false.)

        ! Save plot
        call ModSave_DoAll(Tree,iStep)
    end do

    call MPI_FINALIZE(ierr)

    contains

    subroutine test1_INITIATION
        ! Initialize the Yin-Yang tree
        call YinYangTree_InitTree(Tree, r_range)

        ! Divide all leaf nodes
        call ModAMR_set_grid(Tree)
        
        ! Complete setting up tree.
        call YinYangTree_SetAll(Tree, .true.,.false.)
        
        ! Set up ghost cells
        call ModCommunication_SetGCAll(Tree,.false.)

        ! Do initiation
        call ModInitiation_DoAll(Tree)
        
        if (MpiRank==0 .or. MpiRank==MpiSize-1) then
            print *,'Tree%NumLeafNodes=',Tree%NumLeafNodes,'MpiRank=',MpiRank,' has ',Tree%nLocalBlocks,' blocks'
        end if
    end subroutine test1_INITIATION
end program