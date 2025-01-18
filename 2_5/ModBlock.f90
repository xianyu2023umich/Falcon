module ModBlock

    use ModParameters,  only:   ni,nj,nk,ng,nvar

    type GC_target
        logical                     ::  if_yin                  ! if_yin of other
        integer                     ::  iBlock                  ! i of other Block
        integer                     ::  iRank                   ! i of other Rank
        integer                     ::  nGC                     ! n of GCs
        integer,allocatable         ::  ijk_list(:,:)           ! ijks of GCs
        real,allocatable            ::  xijk_list(:,:)          ! xijk of GCs
        real,allocatable            ::  primitive_list(:,:)     ! primitive of GCs
    end type GC_target

    type BlockType
        integer                     ::  iBlock                  ! i of block
        logical                     ::  if_yin                  ! if yin or yang
        logical                     ::  if_top,if_bottom        ! if top or bottom (or both)
        real                        ::  xijk_range(3,2)         ! range of the whole domain
        real                        ::  xijk_range_GC(3,2)      ! range including GC

        real                        ::  dxi,dxj,dxk             ! grid cell size
        real,allocatable            ::  xi(:),xj(:),xk(:)       ! grid

        real,allocatable            ::  primitive(:,:,:,:)      ! the primitive variables
        real,allocatable            ::  primitive_rk(:,:,:,:)   ! primitive tmp of rk
        integer,allocatable         ::  GC_iBlocks(:,:,:)       ! the iBlocks of EVERY grid (including non-GC)
        type(GC_target),allocatable ::  GC_targets(:)           ! GC targets 
        type(GC_target),allocatable ::  GC_sources(:)           ! GC targets 
        integer,allocatable         ::  requests(:)

        real,allocatable            ::  g_list(:)               ! background gravity
        real,allocatable            ::  p0_list(:)              ! background pressure
        real,allocatable            ::  rho0_list(:)            ! background density
        real,allocatable            ::  te0_list(:)             ! background temperature
        real,allocatable            ::  gamma1_list(:)          ! background gamma1
        real,allocatable            ::  gamma3_list(:)          ! background gamma3
        real,allocatable            ::  diffusion_list(:)       ! background diffusion
        real,allocatable            ::  cooling_list(:)         ! artificial cooling
        real,allocatable            ::  Xi_rsst_list(:)         ! For rsst

        real,allocatable            ::  gamma1(:,:,:)
        real,allocatable            ::  gamma3_minus_1(:,:,:)
        real,allocatable            ::  rho0(:,:,:)
        real,allocatable            ::  g_over_rho0(:,:,:)
        real,allocatable            ::  p0_over_rho0(:,:,:)
        real,allocatable            ::  rho0T0(:,:,:)
        real,allocatable            ::  total_heat(:,:,:)
        real,allocatable            ::  p1(:,:,:)
        real,allocatable            ::  Xi_rsst(:,:,:)

        ! This is important:
        ! GC_sources only contain those in the other ranks
        ! for those in the same rank no need to store it.

        integer                     ::  nGC_targets,nGC_sources
    end type BlockType

    contains

    ! Initialize a single block

    subroutine ModBlock_Init(Block1,iBlock,xijk_range,if_yin,if_SSM,if_use_actual_nvar)
        implicit none
        type(BlockType),target      ::  Block1                  ! the block to be initiated
        integer,intent(in)          ::  iBlock                  ! Global iBlock
        real,intent(in)             ::  xijk_range(3,2)         ! range
        logical,intent(in)          ::  if_yin                  ! if is yin or yang    
        logical,intent(in)          ::  if_SSM                  ! if initialize SSM     
        logical,intent(in)          ::  if_use_actual_nvar      ! if use nvar=5 or 1

        ! set iBlock and if_yin
        Block1%iBlock=iBlock
        Block1%if_yin=if_yin
        Block1%if_bottom=.false.
        Block1%if_top=.false.

        ! Initialize the grid
        call ModBlock_InitGrid(Block1,xijk_range)

        ! Initialize vars
        ! allocate primitive and set it to 0
        if (if_use_actual_nvar) then
            allocate(Block1%primitive   (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,1:nvar))
            allocate(Block1%primitive_rk(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,1:nvar))
        else
            allocate(Block1%primitive   (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,1:1))
            allocate(Block1%primitive_rk(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,1:1))
        end if
        Block1%primitive    =0.
        Block1%primitive_rk =0.
        call random_number(Block1%primitive)
        Block1%primitive=(Block1%primitive-0.5)*1.e-2
        Block1%primitive_rk=Block1%primitive
        ! ModelS initiation
        if (if_SSM) call ModBlock_Init_ModelS(Block1)
        !if (if_SSM) print *,Block1%Xi_rsst_list
        !if (if_SSM .and. Block1%xijk_range(1,2)>0.98) print *,Block1%cooling_list,Block1%xijk_range(1,:)
    end subroutine ModBlock_Init

    subroutine ModBlock_InitGrid(Block1,xijk_range)
        implicit none
        type(BlockType),target      ::  Block1                  ! the block
        real,intent(in)             ::  xijk_range(3,2)         ! range
        integer                     ::  i,j,k                   ! indices

        ! set range and sizes
        Block1%xijk_range=xijk_range                                    
        
        ! allocate the grid
        allocate(Block1%xi(-ng+1:ni+ng))
        allocate(Block1%xj(-ng+1:nj+ng))
        allocate(Block1%xk(-ng+1:nk+ng))
        
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

        Block1%ngC_targets=0
        Block1%ngC_sources=0
    end subroutine ModBlock_InitGrid

    subroutine ModBlock_Init_ModelS(Block1)

        ! Use the scales from ModelS.
        use ModStratification,  only:   ModStratification_get_vars,&
                                        rho0__bar,p0__bar,T0__bar,&
                                        g__bar,heat__bar

        implicit none
        type(BlockType),target      ::  Block1                  ! the block
        real                        ::  g__CGS,rho0__CGS,&      ! ModelS vars of one height
                                        p0__CGS,T0__CGS,&
                                        gamma1,gamma3,&
                                        kap__CGS,&
                                        diffusion__CGS,&
                                        cooling__CGS,&
                                        Xi_rsst
        integer                     ::  i,j,k

        ! Allocate the 1D arrays.
        allocate(Block1%g_list          (-ng+1:ni+ng))
        allocate(Block1%p0_list         (-ng+1:ni+ng))
        allocate(Block1%rho0_list       (-ng+1:ni+ng))
        allocate(Block1%te0_list        (-ng+1:ni+ng))
        allocate(Block1%gamma1_list     (-ng+1:ni+ng))
        allocate(Block1%gamma3_list     (-ng+1:ni+ng))
        allocate(Block1%diffusion_list  (-ng+1:ni+ng))
        allocate(Block1%cooling_list    (-ng+1:ni+ng))
        allocate(Block1%Xi_rsst_list    (-ng+1:ni+ng))

        ! Allocate 3D arrays
        allocate(Block1%gamma1         (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        allocate(Block1%gamma3_minus_1 (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        allocate(Block1%rho0           (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        allocate(Block1%g_over_rho0    (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        allocate(Block1%p0_over_rho0   (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        allocate(Block1%rho0T0         (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        allocate(Block1%total_heat     (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        allocate(Block1%p1             (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        allocate(Block1%Xi_rsst        (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))

        ! Get the 1D profiles
        do i=-ng+1,ni+ng
            call ModStratification_get_vars(Block1%xi(i),&
                g__CGS,rho0__CGS,p0__CGS,T0__CGS,gamma1,gamma3,kap__CGS,&
                diffusion__CGS,cooling__CGS,Xi_rsst)
            
            ! Set the scales
            Block1%g_list(i)        =   g__CGS/g__bar
            Block1%p0_list(i)       =   p0__CGS/p0__bar
            Block1%rho0_list(i)     =   rho0__CGS/rho0__bar
            Block1%te0_list(i)      =   T0__CGS/T0__bar
            Block1%gamma1_list(i)   =   gamma1
            Block1%gamma3_list(i)   =   gamma3
            Block1%diffusion_list(i)=   diffusion__CGS/heat__bar
            Block1%cooling_list(i)  =   cooling__CGS/heat__bar   
            Block1%Xi_rsst_list(i)  =   Xi_rsst
        end do

        ! Fulfill the 3D arrays with 1D profiles
        do j=-ng+1,nj+ng; do k=-ng+1,nk+ng

            Block1%gamma1(:,j,k)        =Block1%gamma1_list
            Block1%gamma3_minus_1(:,j,k)=Block1%gamma3_list-1.0
            Block1%rho0(:,j,k)          =Block1%rho0_list
            Block1%g_over_rho0(:,j,k)   =Block1%g_list/Block1%rho0_list
            Block1%p0_over_rho0(:,j,k)  =Block1%p0_list/Block1%rho0_list
            Block1%rho0T0(:,j,k)        =Block1%rho0_list*Block1%te0_list
            Block1%total_heat(:,j,k)    =Block1%diffusion_list+Block1%cooling_list
            Block1%Xi_rsst(:,j,k)       =Block1%Xi_rsst_list
        end do; end do
    end subroutine ModBlock_Init_ModelS

    subroutine ModBlock_deallocate(Block1)
        implicit none
        type(BlockType),intent(inout)   ::  Block1
        integer                         ::  iGC_target,iGC_source

        if (allocated(Block1%xi))               deallocate(Block1%xi)
        if (allocated(Block1%xj))               deallocate(Block1%xj)
        if (allocated(Block1%xk))               deallocate(Block1%xk)

        if (allocated(Block1%primitive))        deallocate(Block1%primitive)
        if (allocated(Block1%primitive_rk))     deallocate(Block1%primitive_rk)

        if (allocated(Block1%GC_iBlocks))       deallocate(Block1%GC_iBlocks)
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

        if (allocated(Block1%gamma1))           deallocate(Block1%gamma1)
        if (allocated(Block1%gamma3_minus_1))   deallocate(Block1%gamma3_minus_1)
        if (allocated(Block1%rho0))             deallocate(Block1%rho0)
        if (allocated(Block1%g_over_rho0))      deallocate(Block1%g_over_rho0)
        if (allocated(Block1%p0_over_rho0))     deallocate(Block1%p0_over_rho0)
        if (allocated(Block1%rho0T0))           deallocate(Block1%rho0T0)
        if (allocated(Block1%total_heat))       deallocate(Block1%total_heat)
        if (allocated(Block1%Xi_rsst))          deallocate(Block1%Xi_rsst)
        if (allocated(Block1%p1))               deallocate(Block1%p1)

        if (allocated(Block1%g_list))           deallocate(Block1%g_list)
        if (allocated(Block1%p0_list))          deallocate(Block1%p0_list)
        if (allocated(Block1%rho0_list))        deallocate(Block1%rho0_list)
        if (allocated(Block1%te0_list))         deallocate(Block1%te0_list)
        if (allocated(Block1%gamma1_list))      deallocate(Block1%gamma1_list)
        if (allocated(Block1%gamma3_list))      deallocate(Block1%gamma3_list)
        if (allocated(Block1%diffusion_list))   deallocate(Block1%diffusion_list)
        if (allocated(Block1%cooling_list))     deallocate(Block1%cooling_list)
        if (allocated(Block1%p1))               deallocate(Block1%p1)
    end subroutine ModBlock_deallocate

end module ModBlock