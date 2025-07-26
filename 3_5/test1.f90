program test1

    ! This module tests if communications work.

    use ModBlock,           only:   BlockType
    use ModParameters,      only:   MpiSize, MpiRank, r_range, ni, nj, nk, ng,iGeometry,&
                                    nStepsSavePlot,CFL,rSave,r_range,nSteps
    use ModReadParameters,  only:   ModReadParameters_read
    use ModYinYangTree,     only:   YYTree, YinYangTree_InitTree, YinYangTree_SetAll, YinYangTree_DivideAll
    use ModCommunication,   only:   ModCommunication_SetGCAll
    !use ModPFSS,            only:   ModPFSS_setup, ModPFSS_solve
    !use ModMagnetogram,     only:   ModMagnetogram_dipole_magnetogram_ALL
    use ModInitiation,      only:   ModInitiation_harmonic
    !use ModVariables,       only:   rho1_,vr_,vt_,vp_,s1_
    use ModSavePlot,        only:   ModSave_Globe
    use ModCheck,           only:   ModCheck_primitive
    use ModAdvance,         only:   ModAdvance_rk4
    use ModAMR,             only:   ModAMR_set_grid
    use ModStratification,  only:   ModStratification_DoAll
    use MPI

    implicit none

    type(YYTree),target     ::  Tree                        ! The Tree
    integer                 ::  ierr                        ! For MPI
    integer                 ::  iStep                       ! For step
    character(len=8)        ::  iStep_char
    real                    ::  dt                          ! Time step
    integer                 ::  i
    logical                 ::  do_check=.true.
    !integer                 ::  iBlock
    !type(BlockType),pointer ::  Block1

    ! Initiate MPI and get MpiSize/Rank
    call MPI_INIT(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, MpiSize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, MpiRank, ierr)

    ! Read parameters. Do initiation
    call ModReadParameters_read('PARAM.in',1)
    write(*,*)'Parameters read.'
    call test1_INITIATION

    !call ModInitiation_harmonic(Tree)

    !print *,maxval(Tree%LocalBlocks(1)%Xi_rsst_III),minval(Tree%LocalBlocks(1)%Xi_rsst_III)

    ! Main loop
    do iStep=1,nSteps

        ! Check if there is NaN
        !if (do_check) call ModCheck_primitive(Tree,)

        ! Advance
        call ModAdvance_rk4(Tree,CFL,.True.,dt)
        if(MpiRank==0)print *,'Complete Advancing at iStep=',iStep,'dt=',dt

        ! Write plot
        write(iStep_char,'(I8)')iStep
        if (mod(iStep,nStepsSavePlot)==0 .and. iStep>-2400) then

            ! Get the filename and save the plot.
            do i=1,len_trim(iStep_char)
                if (iStep_char(i:i)==' ') iStep_char(i:i)='0'
            end do
            call ModSave_Globe(Tree,0.8,[180,360],'test1_'//iStep_char//'.dat',10)
        end if
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
        
        if (MpiRank==0) print *,'Tree%NumLeafNodes=',Tree%NumLeafNodes
        print *,'MpiRank=',MpiRank,' has ',Tree%nLocalBlocks,' blocks'
    end subroutine test1_INITIATION
end program