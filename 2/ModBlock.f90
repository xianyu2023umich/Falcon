module ModBlock

    type GC_target
        logical                     ::  if_yin                  ! if_yin of other
        integer                     ::  iBlock                  ! i of other Block
        integer                     ::  iRank                   ! i of other Rank
        integer                     ::  nGC                     ! n of GCs
        integer,allocatable         ::  ijk_list(:,:)           ! ijks of GCs
        real,allocatable            ::  xijk_list(:,:)          ! xijk of GCs
        real,allocatable            ::  primitive_list(:,:)     ! primitive of GCs
    end type GC_target

    type Block
        integer                     ::  iBlock                  ! i of block
        logical                     ::  if_yin                  ! if yin or yang
        logical                     ::  if_top,if_bottom        ! if top or bottom (or both)
        real                        ::  xijk_range(3,2)         ! range of the whole domain
        real                        ::  xijk_range_GC(3,2)      ! range including GC

        integer                     ::  nvar,ni,nj,nk,ng        ! variables & grid
        real,allocatable            ::  dxi,dxj,dxk             ! grid cell size
        real,allocatable            ::  xi(:),xj(:),xk(:)       ! grid

        real,allocatable            ::  primitive(:,:,:,:)      ! the primitive variables
        real,allocatable            ::  primitive_rk(:,:,:,:)   ! primitive tmp of rk
        integer,allocatable         ::  GC_iBlocks(:,:,:)       ! the iBlocks of EVERY grid (including non-GC)
        type(GC_target),allocatable ::  GC_targets(:)           ! GC targets 
        type(GC_target),allocatable ::  GC_sources(:)           ! GC targets 

        real,allocatable            ::  p0_list(:)
        real,allocatable            ::  rho0_list(:)
        real,allocatable            ::  te0_list(:)
        real,allocatable            ::  p0(:,:,:)
        real,allocatable            ::  rho0(:,:,:)
        real,allocatable            ::  te0(:,:,:)
        real,allocatable            ::  p1(:,:,:)

        ! This is important:
        !
        ! GC_sources only contain those in the other ranks
        ! for those in the same rank no need to store it.

        integer                     ::  nGC_targets,nGC_sources
    end type Block

    contains

    ! Initialize a single block

    subroutine ModBlock_Init(Block1,iBlock,xijk_range,nvar,ni,nj,nk,ng,if_yin,if_SSM)
        implicit none
        type(Block),target          ::  Block1                  ! the block to be initiated
        integer,intent(in)          ::  iBlock                  ! Global iBlock
        integer,intent(in)          ::  nvar,ni,nj,nk,ng        ! grid
        real,intent(in)             ::  xijk_range(3,2)         ! range
        logical,intent(in)          ::  if_yin                  ! if is yin or yang    
        logical,intent(in)          ::  if_SSM                  ! if initialize SSM     

        ! set iBlock and if_yin
        Block1%iBlock=iBlock
        Block1%if_yin=if_yin
        Block1%if_bottom=.false.
        Block1%if_top=.false.

        ! Initialize the grid
        call ModBlock_InitGrid(Block1,xijk_range,nvar,ni,nj,nk,ng)
        if (if_SSM) call ModBlock_Init_SSM(Block1)
    end subroutine ModBlock_Init

    subroutine ModBlock_InitGrid(Block1,xijk_range,nvar,ni,nj,nk,ng)
        implicit none
        type(Block),target          ::  Block1                  ! the block
        integer,intent(in)          ::  nvar,ni,nj,nk,ng        ! gird size
        real,intent(in)             ::  xijk_range(3,2)         ! range
        integer                     ::  i,j,k                   ! indices

        ! set range and sizes
        Block1%xijk_range=xijk_range
        Block1%ng=ng                                        
        Block1%ni=ni                                 
        Block1%nj=nj                                     
        Block1%nk=nk                                      
        Block1%nvar=nvar
        
        ! allocate the grid
        allocate(Block1%xi(-ng+1:ni+ng))
        allocate(Block1%xj(-ng+1:nj+ng))
        allocate(Block1%xk(-ng+1:nk+ng))

        ! allocate primitive and set it to 0
        allocate(Block1%primitive   (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,1:nvar))
        allocate(Block1%primitive_rk(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,1:nvar))
        Block1%primitive    =0.
        Block1%primitive_rk =0.
        
        ! initialize GC_iBlocks
        allocate(Block1%GC_iBlocks(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))  
        Block1%GC_iBlocks=-1                                           

        ! get xi,xj,xk arrays
        do i=1-ng,ni+ng; Block1%xi(i)=(xijk_range(1,1)*(ni-i+0.5)+xijk_range(1,2)*(i-0.5))/ni; end do
        do j=1-ng,nj+ng; Block1%xj(j)=(xijk_range(2,1)*(nj-j+0.5)+xijk_range(2,2)*(j-0.5))/nj; end do
        do k=1-ng,nk+ng; Block1%xk(k)=(xijk_range(3,1)*(nk-k+0.5)+xijk_range(3,2)*(k-0.5))/nk; end do
        
        ! get dxi, dxj, dxk
        Block1%dxi=(xijk_range(1,2)-xijk_range(1,1))/ni
        Block1%dxj=(xijk_range(2,2)-xijk_range(2,1))/nj
        Block1%dxk=(xijk_range(3,2)-xijk_range(3,1))/nk
        
        ! set xijk_range including GC
        Block1%xijk_range_GC(1,1)=Block1%xijk_range(1,1)-(ng-0.5)*Block1%dxi
        Block1%xijk_range_GC(2,1)=Block1%xijk_range(2,1)-(ng-0.5)*Block1%dxj
        Block1%xijk_range_GC(3,1)=Block1%xijk_range(3,1)-(ng-0.5)*Block1%dxk
        Block1%xijk_range_GC(1,2)=Block1%xijk_range(1,2)+(ng-0.5)*Block1%dxi
        Block1%xijk_range_GC(2,2)=Block1%xijk_range(2,2)+(ng-0.5)*Block1%dxj
        Block1%xijk_range_GC(3,2)=Block1%xijk_range(3,2)+(ng-0.5)*Block1%dxk

        Block1%nGC_targets=0
        Block1%nGC_sources=0
        
        call random_number(Block1%primitive)
        Block1%primitive=(Block1%primitive-0.5)*0.00001
    end subroutine ModBlock_InitGrid

    subroutine ModBlock_Init_SSM(Block1)
        use ModSSM_v0
        use ModParameter,only : paraGamma 

        implicit none
        type(Block),target          ::  Block1                  ! the block
        integer                     ::  j,k
        real                        ::  m

        m=1./(paraGamma-1.)

        allocate(Block1%p0_list     (-Block1%ng+1:Block1%ni+Block1%ng))
        allocate(Block1%rho0_list   (-Block1%ng+1:Block1%ni+Block1%ng))
        allocate(Block1%te0_list    (-Block1%ng+1:Block1%ni+Block1%ng))

        allocate(Block1%p0     (-Block1%ng+1:Block1%ni+Block1%ng,&
                                -Block1%ng+1:Block1%nj+Block1%ng,&
                                -Block1%ng+1:Block1%nk+Block1%ng))
        allocate(Block1%rho0   (-Block1%ng+1:Block1%ni+Block1%ng,&
                                -Block1%ng+1:Block1%nj+Block1%ng,&
                                -Block1%ng+1:Block1%nk+Block1%ng))
        allocate(Block1%te0    (-Block1%ng+1:Block1%ni+Block1%ng,&
                                -Block1%ng+1:Block1%nj+Block1%ng,&
                                -Block1%ng+1:Block1%nk+Block1%ng))
        allocate(Block1%p1     (-Block1%ng+1:Block1%ni+Block1%ng,&
                                -Block1%ng+1:Block1%nj+Block1%ng,&
                                -Block1%ng+1:Block1%nk+Block1%ng))

        Block1%rho0_list=ModSSM_v0_get_var0(Block1%ni+2*Block1%ng,Block1%xi,m,'rho0')
        Block1%p0_list  =ModSSM_v0_get_var0(Block1%ni+2*Block1%ng,Block1%xi,m,'p0'  )
        Block1%te0_list =ModSSM_v0_get_var0(Block1%ni+2*Block1%ng,Block1%xi,m,'te0' )

        do j=-Block1%ng+1,Block1%nj+Block1%ng; do k=-Block1%ng+1,Block1%nk+Block1%ng
            Block1%p0  (:,j,k)=Block1%p0_list
            Block1%rho0(:,j,k)=Block1%rho0_list
            Block1%te0 (:,j,k)=Block1%te0_list
        end do; end do
    end subroutine ModBlock_Init_SSM

    subroutine ModBlock_deallocate(Block1)
        implicit none
        type(Block),intent(inout)   ::  Block1
        integer                     ::  iGC_target,iGC_source
        if (allocated(Block1%xi)) deallocate(Block1%xi)
        if (allocated(Block1%xj)) deallocate(Block1%xj)
        if (allocated(Block1%xk)) deallocate(Block1%xk)

        if (allocated(Block1%primitive)) deallocate(Block1%primitive)
        if (allocated(Block1%primitive_rk)) deallocate(Block1%primitive_rk)

        if (allocated(Block1%GC_iBlocks)) deallocate(Block1%GC_iBlocks)
        if (allocated(Block1%GC_targets)) then
            do iGC_target=1,size(Block1%GC_targets)
                if (allocated(Block1%GC_targets(iGC_target)%ijk_list)) &
                    deallocate(Block1%GC_targets(iGC_target)%ijk_list)
                if (allocated(Block1%GC_targets(iGC_target)%xijk_list)) &
                    deallocate(Block1%GC_targets(iGC_target)%xijk_list)
            end do
            deallocate(Block1%GC_targets)
        end if
        if (allocated(Block1%GC_sources)) then
            do iGC_source=1,size(Block1%GC_sources)
                if (allocated(Block1%GC_sources(iGC_source)%ijk_list)) &
                    deallocate(Block1%GC_sources(iGC_source)%ijk_list)
                if (allocated(Block1%GC_sources(iGC_source)%xijk_list)) &
                    deallocate(Block1%GC_sources(iGC_source)%xijk_list)
            end do
            deallocate(Block1%GC_sources)
        end if

        if (allocated(Block1%p1))       deallocate(Block1%p1)
        if (allocated(Block1%p0))       deallocate(Block1%p0)
        if (allocated(Block1%rho0))     deallocate(Block1%rho0)
        if (allocated(Block1%te0))      deallocate(Block1%te0)
        if (allocated(Block1%p0_list))  deallocate(Block1%p0_list)
        if (allocated(Block1%rho0_list))deallocate(Block1%rho0_list)
        if (allocated(Block1%te0_list)) deallocate(Block1%te0_list)
    end subroutine ModBlock_deallocate

end module ModBlock