module ModBlock
    
    use ModConst,               only:   M_sun__CGS,G_CGS,R_sun__CGS
    use ModParameters,          only:   ni,nj,nk,ng,iGeometry
    use ModStratification_new,  only:   ModStratification_new_get_vars

    type GC_target
        logical                     ::  if_yin                  ! if_yin of other
        integer                     ::  iBlock                  ! i of other Block
        integer                     ::  iRank                   ! i of other Rank
        integer                     ::  nGC                     ! n of GCs
        integer,allocatable         ::  ijk_list(:,:)           ! ijks of GCs
        real(8),allocatable         ::  xijk_list(:,:)          ! xijk of GCs
        real(8),allocatable         ::  primitive_list(:,:)     ! primitive of GCs
    end type GC_target

    ! All the units are in CGS

    type BlockType

        ! This is about the block itself.
        integer                     ::  iBlock                  ! i of block
        logical                     ::  if_initialized=.false.  ! if initialized

        ! About the geometry        
        logical                     ::  if_yin                  ! if yin or yang in spherical
        logical                     ::  if_top,if_bottom        ! if top or bottom (or both)
        real(8)                     ::  xijk_range(3,2)         ! range of the whole domain
        real(8)                     ::  xijk_range_GC(3,2)      ! range including GC

        ! About the grid
        real(8)                     ::  dxi,dxj,dxk             ! Grid cell size in CGS
        real(8),allocatable         ::  xi_I(:),xj_I(:),xk_I(:) ! Cell center array
        real(8),allocatable         ::  xi_F(:),xj_F(:),xk_F(:) ! Cell faces array
        real(8),allocatable         ::  V_LLL(:,:,:)            ! Cell volume
        real(8),allocatable         ::  Si_FLL(:,:,:),&         ! Face areas
                                        Sj_LFL(:,:,:),&
                                        Sk_LLF(:,:,:)
        real(8),allocatable         ::  Di_FLL(:,:,:),&         ! Distance between adjacent cell centers.
                                        Dj_LFL(:,:,:),&
                                        Dk_LLF(:,:,:)

        ! About the variables
        integer                     ::  itype_equation          ! type of equation.
        integer                     ::  nvar
        integer                     ::  rho1_=-1,&
                                        rho_=-1,&
                                        vr_=-1,&
                                        vt_=-1,&
                                        vp_=-1,&
                                        s1_=-1,&
                                        br_=-1,&
                                        bt_=-1,&
                                        bp_=-1,&
                                        psi_=-1,&
                                        te_=-1,&
                                        te1_=-1,&
                                        w_plus_=-1,&
                                        w_minus_=-1

        real(8),pointer             ::  primitive(:,:,:,:)      ! the pointer. To be used in the equations.
        real(8),allocatable         ::  primitive_IV(:,:,:,:)   ! the primitive variables
        real(8),allocatable         ::  primitive_rk_IV(:,:,:,:)! primitive tmp of rk
        real(8),allocatable         ::  primitive_corner_IV(:,:,:,:) ! primitive at cell corners (for output)
        real(8),allocatable         ::  EQN_update_R_IV(:,:,:,:)! update of R during RK
        integer,allocatable         ::  GC_iBlocks_III(:,:,:)   ! the iBlocks of EVERY grid (including non-GC)

        ! About GC
        type(GC_target),allocatable ::  GC_targets(:)           ! GC targets 
        type(GC_target),allocatable ::  GC_sources(:)           ! GC targets 
        integer,allocatable         ::  GC_requests(:)
        integer                     ::  nSend,nRecv

        ! About the convection zone
        logical                     ::  if_SSM
        logical                     ::  if_involve_B = .false.  ! if magnetic field is active for this block
        real(8),allocatable         ::  g_I(:)               ! background gravity
        real(8),allocatable         ::  p0_I(:)              ! background pressure
        real(8),allocatable         ::  rho0_I(:)            ! background density
        real(8),allocatable         ::  te0_I(:)             ! background temperature
        real(8),allocatable         ::  gamma1_I(:)          ! background gamma1
        real(8),allocatable         ::  gamma3_I(:)          ! background gamma3
        real(8),allocatable         ::  diffusion_I(:)       ! background diffusion
        real(8),allocatable         ::  cooling_I(:)         ! artificial cooling
        real(8),allocatable         ::  diffusion_flux_I(:)  ! background diffusion flux
        real(8),allocatable         ::  cooling_flux_I(:)    ! artificial cooling flux
        real(8),allocatable         ::  Xi_rsst_I(:)         ! For rsst

        ! 3D arrays for convection zone
        real(8),allocatable         ::  gamma1_III(:,:,:)
        real(8),allocatable         ::  gamma3_minus_1_III(:,:,:)
        real(8),allocatable         ::  rho0_III(:,:,:)
        real(8),allocatable         ::  g_over_rho0_III(:,:,:)
        real(8),allocatable         ::  p0_over_rho0_III(:,:,:)
        real(8),allocatable         ::  rho0T0_III(:,:,:)
        real(8),allocatable         ::  total_heat_III(:,:,:)
        real(8),allocatable         ::  p1_III(:,:,:)
        real(8),allocatable         ::  Xi_rsst_III(:,:,:)
        real(8),allocatable         ::  v_wave_III(:,:,:)
        real(8),allocatable         ::  h_LLL(:,:,:)            ! local min physical cell size

        ! 3D arrays for corona
        real(8),allocatable         ::  p_III(:,:,:)
        real(8),allocatable         ::  g_III(:,:,:)

        ! This is important:
        ! GC_sources only contain those in the other ranks
        ! for those in the same rank no need to store it.

        integer                     ::  nGC_targets,nGC_sources
    end type BlockType

    contains

    subroutine ModBlock_Init(Block1,iBlock,xijk_range,if_yin,itype_equation,if_use_actual_nvar,&
        if_involve_B_here)
        implicit none
        type(BlockType),target      ::  Block1                  ! the block
        integer,intent(in)          ::  iBlock                  ! Global iBlock
        real(8),intent(in)          ::  xijk_range(3,2)         ! range
        logical,intent(in)          ::  if_yin                  ! if is yin or yang   
        integer,intent(in)          ::  itype_equation          ! type of equation
        logical,intent(in)          ::  if_use_actual_nvar
        logical,intent(in)          ::  if_involve_B_here       ! if involve B for this block

        ! Set iBlock, if_yin,itype_equation
        Block1%iBlock=iBlock
        Block1%if_yin=if_yin
        Block1%itype_equation=itype_equation
        
        ! If not HD, then set if_involve_B. Other wise always false.
        if (itype_equation==0) then
            Block1%if_involve_B=.false.
        else
            Block1%if_involve_B=if_involve_B_here
        end if

        ! Initialize the boundaries.
        Block1%if_bottom=.false.
        Block1%if_top=.false.

        call ModBlock_InitGrid(Block1,xijk_range)
        select case(itype_equation)
        case(0) ! Dynamo HD
            Block1%nvar=5
            Block1%rho1_=1
            Block1%vr_=2
            Block1%vt_=3
            Block1%vp_=4
            Block1%s1_=5
            call ModBlock_InitPrimitives(Block1,if_use_actual_nvar)
            call ModBlock_InitSSM(Block1)
        case(1) ! Dynamo MHD
            Block1%nvar=9
            Block1%rho1_=1
            Block1%vr_=2
            Block1%vt_=3
            Block1%vp_=4
            Block1%br_=5
            Block1%bt_=6
            Block1%bp_=7
            Block1%s1_=8
            Block1%psi_=9
            call ModBlock_InitPrimitives(Block1,if_use_actual_nvar)
            call ModBlock_InitSSM(Block1)
        case(2) ! Corona MHD
            Block1%nvar=10
            Block1%rho_=1
            Block1%vr_=2
            Block1%vt_=3
            Block1%vp_=4
            Block1%br_=5
            Block1%bt_=6
            Block1%bp_=7
            Block1%te_=8
            Block1%w_plus_=9
            Block1%w_minus_=10
            call ModBlock_InitPrimitives(Block1,if_use_actual_nvar)
            call ModBlock_InitCorona(Block1)
        case default
            stop 'ModBlock_Init: Invalid itype_equation'
        end select        

        ! It's been initialized.
        Block1%if_initialized=.true.
    end subroutine ModBlock_Init

    subroutine ModBlock_InitSSM(Block1)
        implicit none
        type(BlockType),target      ::  Block1                  ! the block
        real(8)                     ::  g__CGS,rho0__CGS,&      ! ModelS vars of one height
                                        p0__CGS,T0__CGS,&
                                        gamma1,gamma3,&
                                        kap__CGS,&
                                        diffusion__CGS,&
                                        cooling__CGS,&
                                        diffusion_flux__CGS,&
                                        cooling_flux__CGS,&
                                        Xi_rsst
        integer                     ::  i,j,k

        Block1%if_SSM=.true.

        ! Allocate the 1D arrays.
        allocate(Block1%g_I             (-ng+1:ni+ng))
        allocate(Block1%p0_I            (-ng+1:ni+ng))
        allocate(Block1%rho0_I          (-ng+1:ni+ng))
        allocate(Block1%te0_I           (-ng+1:ni+ng))
        allocate(Block1%gamma1_I        (-ng+1:ni+ng))
        allocate(Block1%gamma3_I        (-ng+1:ni+ng))
        allocate(Block1%diffusion_I     (-ng+1:ni+ng))
        allocate(Block1%cooling_I       (-ng+1:ni+ng))
        allocate(Block1%diffusion_flux_I(-ng+1:ni+ng))
        allocate(Block1%cooling_flux_I  (-ng+1:ni+ng))
        allocate(Block1%Xi_rsst_I       (-ng+1:ni+ng))

        ! Allocate 3D arrays
        allocate(Block1%gamma1_III          (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        allocate(Block1%gamma3_minus_1_III  (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        allocate(Block1%rho0_III            (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        allocate(Block1%g_over_rho0_III     (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        allocate(Block1%p0_over_rho0_III    (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        allocate(Block1%rho0T0_III          (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        allocate(Block1%total_heat_III      (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        allocate(Block1%p1_III              (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        allocate(Block1%Xi_rsst_III         (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))

        ! Allocate wave speed array
        allocate(Block1%v_wave_III          (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))

        ! Get the 1D profiles
        do i=-ng+1,ni+ng
           call ModStratification_new_get_vars(Block1%xi_I(i),&
               g__CGS,rho0__CGS,p0__CGS,T0__CGS,gamma1,gamma3,kap__CGS,&
               diffusion__CGS,cooling__CGS,diffusion_flux__CGS,cooling_flux__CGS,Xi_rsst)
            
           ! Set the scales
           Block1%g_I(i)        =   g__CGS
           Block1%p0_I(i)       =   p0__CGS
           Block1%rho0_I(i)     =   rho0__CGS
           Block1%te0_I(i)      =   T0__CGS
           Block1%gamma1_I(i)   =   gamma1
           Block1%gamma3_I(i)   =   gamma3
           Block1%diffusion_I(i)=   diffusion__CGS
           Block1%cooling_I(i)  =   cooling__CGS   
           Block1%diffusion_flux_I(i)=diffusion_flux__CGS
           Block1%cooling_flux_I(i)  =cooling_flux__CGS
           Block1%Xi_rsst_I(i)  =   Xi_rsst
        end do

        ! Fulfill the 3D arrays with 1D profiles
        do j=-ng+1,nj+ng
           do k=-ng+1,nk+ng

                Block1%gamma1_III(:,j,k)        =Block1%gamma1_I
                Block1%gamma3_minus_1_III(:,j,k)=Block1%gamma3_I-1.0
                Block1%rho0_III(:,j,k)          =Block1%rho0_I
                Block1%g_over_rho0_III(:,j,k)   =Block1%g_I/Block1%rho0_I
                Block1%p0_over_rho0_III(:,j,k)  =Block1%p0_I/Block1%rho0_I
                Block1%rho0T0_III(:,j,k)        =Block1%rho0_I*Block1%te0_I
                Block1%total_heat_III(:,j,k)    =Block1%diffusion_I+Block1%cooling_I
                Block1%Xi_rsst_III(:,j,k)       =Block1%Xi_rsst_I
           end do
        end do
    end subroutine ModBlock_InitSSM
    
    subroutine ModBlock_InitPrimitives(Block1,if_use_actual_nvar)
        implicit none
        type(BlockType),target      ::  Block1                  ! the block
        logical,intent(in)          ::  if_use_actual_nvar      ! if no then nvar_here=1

        ! Initialize vars
        ! allocate primitive and set it to 0
        if (if_use_actual_nvar) then
            allocate(Block1%primitive_IV    (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,1:Block1%nvar))
            allocate(Block1%primitive_rk_IV (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,1:Block1%nvar))
            allocate(Block1%EQN_update_R_IV (1:ni,1:nj,1:nk,1:Block1%nvar))
        else
            allocate(Block1%primitive_IV    (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,1:1))
            allocate(Block1%primitive_rk_IV (-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,1:1))
            allocate(Block1%EQN_update_R_IV (1:ni,1:nj,1:nk,1:1))
        end if
        Block1%primitive_IV     = 0.
        Block1%primitive_rk_IV  = 0.
        Block1%EQN_update_R_IV  = 0.
    end subroutine ModBlock_InitPrimitives
    
    subroutine ModBlock_InitGrid(Block1,xijk_range)
        implicit none
        type(BlockType),target      ::  Block1                  ! the block
        real(8),intent(in)          ::  xijk_range(3,2)         ! range
        integer                     ::  i,j,k                   ! indices

        ! set range and sizes and geometry and if_SSM
        Block1%xijk_range=xijk_range                              
        
        ! allocate the grid
        allocate(Block1%xi_I(-ng+1:ni+ng))
        allocate(Block1%xj_I(-ng+1:nj+ng))
        allocate(Block1%xk_I(-ng+1:nk+ng))
        allocate(Block1%xi_F(1:ni+1))
        allocate(Block1%xj_F(1:nj+1))
        allocate(Block1%xk_F(1:nk+1))
        
        ! initialize GC_iBlocks
        allocate(Block1%GC_iBlocks_III(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))  
        Block1%GC_iBlocks_III=-1                                           

        ! Get the positions of the cell centers
        do i=1-ng,ni+ng; Block1%xi_I(i)=(xijk_range(1,1)*(ni-i+0.5)+xijk_range(1,2)*(i-0.5))/ni; end do
        do j=1-ng,nj+ng; Block1%xj_I(j)=(xijk_range(2,1)*(nj-j+0.5)+xijk_range(2,2)*(j-0.5))/nj; end do
        do k=1-ng,nk+ng; Block1%xk_I(k)=(xijk_range(3,1)*(nk-k+0.5)+xijk_range(3,2)*(k-0.5))/nk; end do
        
        ! Get dxi, dxj, dxk
        Block1%dxi=(xijk_range(1,2)-xijk_range(1,1))/ni
        Block1%dxj=(xijk_range(2,2)-xijk_range(2,1))/nj
        Block1%dxk=(xijk_range(3,2)-xijk_range(3,1))/nk

        ! Get the positions of the cell faces
        do i=1,ni+1; Block1%xi_F(i)=(xijk_range(1,1)*(ni+1-i)+xijk_range(1,2)*(i-1))/ni; end do
        do j=1,nj+1; Block1%xj_F(j)=(xijk_range(2,1)*(nj+1-j)+xijk_range(2,2)*(j-1))/nj; end do
        do k=1,nk+1; Block1%xk_F(k)=(xijk_range(3,1)*(nk+1-k)+xijk_range(3,2)*(k-1))/nk; end do
        
        ! Set the face areas and cell volumes based on the geometry
        allocate(   Block1%Si_FLL(1:ni+1,1:nj,1:nk),&
                    Block1%Sj_LFL(1:ni,1:nj+1,1:nk),&
                    Block1%Sk_LLF(1:ni,1:nj,1:nk+1))
        allocate(   Block1%Di_FLL(1:ni+1,1:nj,1:nk),&
                    Block1%Dj_LFL(1:ni,1:nj+1,1:nk),&
                    Block1%Dk_LLF(1:ni,1:nj,1:nk+1))
        allocate(   Block1%V_LLL(1:ni,1:nj,1:nk))

        select case(iGeometry)
        case(0)
            do k=1,nk
                do j=1,nj
                    do i=1,ni+1
                        Block1%Si_FLL(i,j,k)=(Block1%xj_F(j+1)-Block1%xj_F(j))*(Block1%xk_F(k+1)-Block1%xk_F(k))
                    end do
                end do 
            end do
            do k=1,nk
                do j=1,nj+1
                    do i=1,ni
                        Block1%Sj_LFL(i,j,k)=(Block1%xi_F(i+1)-Block1%xi_F(i))*(Block1%xk_F(k+1)-Block1%xk_F(k))
                    end do
                end do
            end do
            do k=1,nk+1
                do j=1,nj
                    do i=1,ni
                        Block1%Sk_LLF(i,j,k)=(Block1%xj_F(j+1)-Block1%xj_F(j))*(Block1%xi_F(i+1)-Block1%xi_F(i))
                    end do
                end do
            end do
            do k=1,nk
                do j=1,nj
                    do i=1,ni
                        Block1%V_LLL(i,j,k)=(Block1%xi_F(i+1)-Block1%xi_F(i))*&
                                    (Block1%xj_F(j+1)-Block1%xj_F(j))*(Block1%xk_F(k+1)-Block1%xk_F(k))
                    end do
                end do
            end do
            Block1%Di_FLL=Block1%dxi
            Block1%Dj_LFL=Block1%dxj
            Block1%Dk_LLF=Block1%dxk
        case(1)
            do k=1,nk
                do j=1,nj
                    do i=1,ni+1
                        Block1%Si_FLL(i,j,k)=Block1%xi_F(i)**2*&
                                    (cos(Block1%xj_F(j))-cos(Block1%xj_F(j+1)))*&
                                    (Block1%xk_F(k+1)-Block1%xk_F(k))
                    end do
                end do
            end do
            do k=1,nk
                do j=1,nj+1
                    do i=1,ni
                        Block1%Sj_LFL(i,j,k)=0.5*(Block1%xi_F(i+1)**2-Block1%xi_F(i)**2)*(Block1%xk_F(k+1)-Block1%xk_F(k))*&
                                        sin(Block1%xj_F(j))
                    end do
                end do
            end do
            do k=1,nk+1
                do j=1,nj
                    do i=1,ni
                        Block1%Sk_LLF(i,j,k)=0.5*(Block1%xi_F(i+1)**2-Block1%xi_F(i)**2)*(Block1%xj_F(j+1)-Block1%xj_F(j))
                    end do
                end do
            end do
            do k=1,nk
                do j=1,nj
                    do i=1,ni
                        Block1%V_LLL(i,j,k)=1.0/3.0*(Block1%xi_F(i+1)**3-Block1%xi_F(i)**3)*&
                                    (cos(Block1%xj_F(j))-cos(Block1%xj_F(j+1)))*&
                                    (Block1%xk_F(k+1)-Block1%xk_F(k))
                    end do
                end do
            end do
            Block1%Di_FLL=Block1%dxi
            do k=1,nk
                do j=1,nj+1
                    do i=1,ni
                        Block1%Dj_LFL(i,j,k)=Block1%xi_I(i)*Block1%dxj
                    end do
                end do
            end do
            do k=1,nk+1
                do j=1,nj
                    do i=1,ni
                        Block1%Dk_LLF(i,j,k)=Block1%xi_I(i)*Block1%dxk*sin(Block1%xj_I(j))
                    end do
                end do
            end do
            allocate(Block1%h_LLL(1:ni,1:nj,1:nk))
            do k=1,nk
                do j=1,nj
                    do i=1,ni
                        Block1%h_LLL(i,j,k)=min(Block1%dxi, &
                            Block1%xi_I(i)*Block1%dxj, &
                            Block1%xi_I(i)*sin(Block1%xj_I(j))*Block1%dxk)
                    end do
                end do
            end do
        end select
        
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

    subroutine ModBlock_InitCorona(Block1)
        implicit none
        type(BlockType),target      ::  Block1                  ! the block
        integer                     ::  j,k                     ! jk

        ! Allocate 3D arrays
        allocate(Block1%p_III(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        allocate(Block1%g_III(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
        Block1%p_III=0.

        ! Get the gravity based on solar mass and radius
        do k=1,nk
            do j=1,nj
                Block1%g_III(:,j,k)=G_CGS*M_sun__CGS/(Block1%xi_I(1:ni)**2)
            end do
        end do
    end subroutine ModBlock_InitCorona

    subroutine ModBlock_deallocate(Block1)
        implicit none
        type(BlockType)             ::  Block1

        if (Block1%if_initialized) then
            deallocate(Block1%xi_I)
            deallocate(Block1%xj_I)
            deallocate(Block1%xk_I)
            deallocate(Block1%xi_F)
            deallocate(Block1%xj_F)
            deallocate(Block1%xk_F)
            deallocate(Block1%GC_iBlocks_III)
            deallocate(Block1%Si_FLL)
            deallocate(Block1%Sj_LFL)
            deallocate(Block1%Sk_LLF)
            deallocate(Block1%Di_FLL)
            deallocate(Block1%Dj_LFL)
            deallocate(Block1%Dk_LLF)
            deallocate(Block1%V_LLL)
            deallocate(Block1%primitive_IV)
            deallocate(Block1%primitive_rk_IV)
            deallocate(Block1%EQN_update_R_IV)
        end if

        if (Block1%if_SSM) then
            deallocate(Block1%g_I         )
            deallocate(Block1%p0_I        )
            deallocate(Block1%rho0_I      )
            deallocate(Block1%te0_I       )
            deallocate(Block1%gamma1_I    )
            deallocate(Block1%gamma3_I    )
            deallocate(Block1%diffusion_I )
            deallocate(Block1%cooling_I   )
            deallocate(Block1%diffusion_flux_I)
            deallocate(Block1%cooling_flux_I)
            deallocate(Block1%Xi_rsst_I   )
    
            ! DeAllocate 3D arrays
            deallocate(Block1%gamma1_III         )
            deallocate(Block1%gamma3_minus_1_III )
            deallocate(Block1%rho0_III           )
            deallocate(Block1%g_over_rho0_III    )
            deallocate(Block1%p0_over_rho0_III   )
            deallocate(Block1%rho0T0_III         )
            deallocate(Block1%total_heat_III     )
            deallocate(Block1%p1_III             )
            deallocate(Block1%Xi_rsst_III        )
    
            ! DeAllocate wave speed array
            deallocate(Block1%v_wave_III         )
        end if
    end subroutine ModBlock_deallocate
end module ModBlock