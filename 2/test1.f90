program test1
    
    !use ModCommunication
    use ieee_arithmetic
    use ModTimeStep
    use ModSavePlot
    use ModAdvance

    implicit none

    type(YYTree),target     ::  Tree                        ! The Tree
    integer                 ::  nr,nt,np,ng=2,i             ! The grid
    real                    ::  r_range(2)                  ! r range
    real                    ::  r_save=2.05
    integer                 ::  nLevel=3,iLevel

    integer                 ::  ierr,MpiSize,MpiRank        ! For MPI
    integer                 ::  iStep                       ! For step
    character(len=8)        ::  iStep_char                  ! char of istep
    logical                 ::  if_advance
    real                    ::  dt

    integer                 ::  ir,it,ip,iLocalBlock        ! For test/debug
    type(Block),pointer     ::  Block1
    nr=30
    nt=24
    np=72

    if_advance=.true.

    ! Initiate MPI and get MpiSize/Rank
    call MPI_INIT(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, MpiSize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, MpiRank, ierr)

    ! Initiate the tree and set GC
    call test1_INITIATION
    !call test_harmonic
    !call test_yinyang

    if(.not.if_advance) then
        call ModSave_Globe(Tree,r_save,[180,360],filename='00000000.dat',logical_unit=2)
        call ModAdvance_CommunicateAll(Tree,.false.,MpiSize,MpiRank)
        !if(MpiRank==3)print *,Tree%LocalBlocks(1)%primitive(:,10,10,10)
        call ModAdvance_rk4(Tree,r_save,.True.,MpiSize,MpiRank,dt)
        call ModSave_Globe(Tree,1.51,[180,360],filename='00000001.dat',logical_unit=2)
        !if(MpiRank==3)print *,Tree%LocalBlocks(1)%primitive(:,10,10,10)
    end if

    !call test_yinyang
    !print *,Tree%LocalBlocks(32)%primitive(:,26,18,2)
    !print *,Tree%LocalBlocks(30)%primitive(2:4,10,10,10)
    !print *,ModYinYang_VecConv_0D([Tree%LocalBlocks(30)%xi(10),&
    !    Tree%LocalBlocks(30)%xj(10),&
    !    Tree%LocalBlocks(30)%xk(10)],Tree%LocalBlocks(30)%primitive(2:4,10,10,10))
    !call ModAdvance_CommunicateAll(Tree,.false.,MpiSize,MpiRank)
    !call ModSave_Globe(Tree,2.1,[180,360],filename='test.dat',logical_unit=2)
    !print *,'done.'

    !call ModSave_2D(Tree,'k',2.16,[200,200],filename='0.dat',logical_unit=2)

    !call ModSave_Globe(Tree,1.51,[180,360],filename='-00000000.dat',logical_unit=2)

    if (if_advance) then
        do iStep=0,8000
            if(MpiRank==0)print *,'Start Advancing at iStep=',iStep
            call ModAdvance_rk4(Tree,0.7,.True.,MpiSize,MpiRank,dt)
            if(MpiRank==0)print *,'End Advancing at iStep=',iStep
            if(MpiRank==0)print *,maxval(Tree%Localblocks(1)%primitive),dt

            write(iStep_char,'(I8)')iStep
            if (mod(iStep,50)==0) then
                    
                do i=1,len_trim(iStep_char)
                    if (iStep_char(i:i)==' ') iStep_char(i:i)='0'
                end do
                call ModSave_Globe(Tree,r_save,[180,360],filename=iStep_char//'.dat',logical_unit=2)
            end if
            if (MpiRank==0) print *,'iStep=',iStep,'Completed...'

            do iLocalBlock=1,Tree%nLocalBlocks
                Block1=>Tree%LocalBlocks(iLocalBlock)
                do ip=-ng+1,ng+np; do it=-ng+1,ng+nt; do ir=-ng+1,ng+nr
                    if(ieee_is_nan(Block1%primitive(ir,it,ip,1))) then
                        print *,Block1%iBlock,ir,it,ip
                        stop 1
                    end if
                end do; end do; end do
            end do
            
        end do
    end if

    ! Finalize the program.
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call mpi_finalize(ierr)

    contains

    subroutine test1_INITIATION
        !integer ::  iChild1
        r_range=[1.0,2.1]

        call YinYangTree_InitTree(Tree,r_range)
        do iLevel=1,nLevel-1
            print *,'Dividing iLevel=',iLevel
            call YinYangTree_DivideAll(Tree,[.False.,.True.,.True.])
            print *,'End Dividing iLevel=',iLevel
        end do
        call YinYangTree_SetAll(Tree,nr,nt,np,ng,MpiSize,MpiRank)
        print *,'Set done.'
        call ModCommunication_SetGCAll(Tree,MpiSize,MpiRank)
        print *,'Set GC done.'
        print *,'I"m ',MpiRank,', I have',Tree%nLocalBlocks
        call ModAdvance_CommunicateAll(Tree,.false.,MpiSize,MpiRank)
    end subroutine test1_INITIATION

    subroutine test_harmonic
        real    ::  vec(3)
        real    ::  coord(3)

        do iLocalBlock=1,Tree%nLocalBlocks
            Block1=>Tree%LocalBlocks(iLocalBlock)
            do ip=-ng+1,ng+np; do it=-ng+1,ng+nt; do ir=-ng+1,ng+nr
                Block1%primitive(:,ir,it,ip)=0.000000
                coord=[Block1%xi(ir),Block1%xj(it),Block1%xk(ip)]
                if (.not. Block1%if_yin) coord=ModYinyang_CoordConv_0D(coord)
                vec=sin(dpi*(coord(1)-r_range(1))/(r_range(2)-r_range(1)))*[1.,0.,1.]*0.0000001*&
                    sin(dpi*4.*coord(3))*sin(coord(2))**4
                Block1%primitive(2:4,ir,it,ip)=vec
                if (.not. Block1%if_yin) then
                    Block1%primitive(2:4,ir,it,ip)=&
                    ModYinYang_VecConv_0D(coord,vec)
                end if
            end do; end do; end do
        end do
        !call ModSave_Globe(Tree,r_save,[180,360],filename='-00000000.dat',logical_unit=2)
        call ModAdvance_CommunicateAll(Tree,.false.,MpiSize,MpiRank)
        call ModSave_Globe(Tree,r_save,[180,360],filename='+00000000.dat',logical_unit=2)
    end subroutine test_harmonic



    subroutine test_yinyang
        do iLocalBlock=1,Tree%nLocalBlocks
            Block1=>Tree%LocalBlocks(iLocalBlock)
            do ip=-ng+1,ng+np; do it=-ng+1,ng+nt; do ir=-ng+1,ng+nr
                Block1%primitive(1:5,ir,it,ip)=1.
                if (.not. Block1%if_yin) then
                    Block1%primitive(2:4,ir,it,ip)=&
                    ModYinYang_VecConv_0D(&
                    ModYinyang_CoordConv_0D([Block1%xi(ir),Block1%xj(it),Block1%xk(ip)]),&
                    [1.,1.,1.])
                end if
                
            end do; end do; end do
        end do
    end subroutine test_yinyang

end program test1