program test1
    
    use ModCommunication,   only:   ModCommunication_SetGCAll
    use ModAMR,             only:   ModAMR_set_grid
    use ModCheck,           only:   ModCheck_primitive
    use ModSavePlot,        only:   ModSave_Globe
    use ModAdvance,         only:   ModAdvance_rk4,&
                                    ModAdvance_CommunicateAll
    use ModParameters,      only:   r_range,&
                                    ModelS_filename,nSteps,&
                                    nStepsSavePlot,nthSavePlot,nphSavePlot,rSave,CFL
    use ModReadParameters,  only:   ModReadParameters_read
    use ModYinYangTree,     only:   YYTree,&
                                    YinYangTree_InitTree,&
                                    YinYangTree_SetAll
    use ModStratification,  only:   ModStratification_DoAll
    use MPI

    implicit none

    type(YYTree),target     ::  Tree                        ! The Tree
    integer                 ::  i                           ! The grid

    integer                 ::  ierr,MpiSize,MpiRank        ! For MPI
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
            call ModAdvance_rk4(Tree,CFL,.True.,MpiSize,MpiRank,dt)
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
        call YinYangTree_SetAll(Tree,MpiSize,MpiRank)
        print *,'I"m ',MpiRank,', I have',Tree%nLocalBlocks
        call ModCommunication_SetGCAll(Tree,MpiSize,MpiRank)
        call ModAdvance_CommunicateAll(Tree,.false.,MpiSize,MpiRank)
        call ModAdvance_CommunicateAll(Tree,.false.,MpiSize,MpiRank)
    end subroutine test1_INITIATION

    subroutine test_harmonic
        use ModParameters,  only:   ni,nj,nk,ng
        use ModYinYang,     only:   ModYinyang_CoordConv_0D,&
                                    ModYinYang_VecConv_0D
        use ModBlock
        use ModConst,       only:   dpi
        implicit none
        real                    ::  vec(3)
        real                    ::  coord(3)
        type(Block),pointer     ::  Block1
        integer                 ::  iLocalBlock,ir,it,ip

        do iLocalBlock=1,Tree%nLocalBlocks
            Block1=>Tree%LocalBlocks(iLocalBlock)
            do ip=-ng+1,ng+nk; do it=-ng+1,ng+nj; do ir=-ng+1,ng+ni
                Block1%primitive(ir,it,ip,:)=0.000000
                coord=[Block1%xi(ir),Block1%xj(it),Block1%xk(ip)]
                if (.not. Block1%if_yin) coord=ModYinyang_CoordConv_0D(coord)
                vec=sin(dpi*(coord(1)-r_range(1))/(r_range(2)-r_range(1)))*[1.,1.,1.]*1.e-2*&
                    sin(dpi*4.*coord(3))*sin(coord(2))**4
                Block1%primitive(ir,it,ip,2:4)=vec
                if (.not. Block1%if_yin) then
                    Block1%primitive(ir,it,ip,2:4)=&
                    ModYinYang_VecConv_0D(coord,vec)
                end if
            end do; end do; end do
            Block1%primitive_rk=Block1%primitive
        end do
        call ModSave_Globe(Tree,rSave,[nthSavePlot,nphSavePlot],filename='-00000000.dat',logical_unit=2)
        call ModAdvance_CommunicateAll(Tree,.false.,MpiSize,MpiRank)
        call ModSave_Globe(Tree,rSave,[nthSavePlot,nphSavePlot],filename='+00000000.dat',logical_unit=2)
    end subroutine test_harmonic

end program test1