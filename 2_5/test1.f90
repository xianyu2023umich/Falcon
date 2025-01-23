program test1
    
    use ModCommunication,   only:   ModCommunication_SetGCAll
    use ModAMR,             only:   ModAMR_set_grid
    use ModCheck,           only:   ModCheck_primitive
    use ModSavePlot,        only:   ModSave_Globe
    use ModAdvance,         only:   ModAdvance_rk4,&
                                    ModAdvance_CommunicateAll
    use ModParameters,      only:   MpiSize,MpiRank,r_range,&
                                    ModelS_filename,nSteps,&
                                    nStepsSavePlot,nthSavePlot,nphSavePlot,rSave,CFL
    use ModReadParameters,  only:   ModReadParameters_read
    use ModYinYangTree,     only:   YYTree,&
                                    YinYangTree_InitTree,&
                                    YinYangTree_SetAll
    use ModStratification,  only:   ModStratification_DoAll
    use ModInitiation,      only:   ModInitiation_harmonic
    use ModAllocation,      only:   ModAllocation_GetRank
    use MPI

    implicit none

    type(YYTree),target     ::  Tree                        ! The Tree
    integer                 ::  i                           ! The grid

    integer                 ::  ierr                        ! For MPI
    integer                 ::  iStep                       ! For step
    character(len=8)        ::  iStep_char                  ! char of istep
    logical                 ::  if_advance
    logical                 ::  do_check
    real                    ::  dt

    if_advance=.true.
    do_check=.true.

    ! Initiate MPI and get MpiSize/Rank
    call MPI_INIT(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, MpiSize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, MpiRank, ierr)

    ! Read parameters
    call ModReadParameters_read("PARAM.in",1)
    call ModStratification_DoAll(ModelS_filename)

    ! Initiate the tree and set GC
    call test1_INITIATION
    !call test_harmonic
    !call ModSave_Globe(Tree,rSave,[nthSavePlot,nphSavePlot],filename='00000000.dat',logical_unit=2)

    ! Start the main loop !
    if (if_advance) then
        do iStep=1,nSteps

            !if(MpiRank==2) print *,Tree%LocalBlocks(2)%if_yin,Tree%LocalBlocks(2)%xijk_range
            

            ! Advance
            call ModAdvance_rk4(Tree,CFL,.True.,dt)
            if(MpiRank==0)print *,'Complete Advancing at iStep=',iStep,'dt=',dt

            ! Write plot
            write(iStep_char,'(I8)')iStep
            if (mod(iStep,nStepsSavePlot)==0) then
                

                do i=1,len_trim(iStep_char); if (iStep_char(i:i)==' ') iStep_char(i:i)='0'; end do
                call ModSave_Globe(Tree,rSave,[nthSavePlot,nphSavePlot],filename=iStep_char//'.dat',logical_unit=2)
            end if

            ! Check if there is NaN
            if (do_check) call ModCheck_primitive(Tree,MpiRank)

            
        end do
    end if

    ! Finalize the program.
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call mpi_finalize(ierr)

    contains

    subroutine test1_INITIATION
        call YinYangTree_InitTree(Tree,r_range)
        call ModAMR_set_grid(Tree)
        call YinYangTree_SetAll(Tree)
        print *,'I"m ',MpiRank,', I have',Tree%nLocalBlocks
        call ModCommunication_SetGCAll(Tree)
        call ModAdvance_CommunicateAll(Tree,.false.)
        call ModAdvance_CommunicateAll(Tree,.false.)
    end subroutine test1_INITIATION
end program test1