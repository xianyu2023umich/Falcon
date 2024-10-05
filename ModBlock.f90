module ModBlock

    use ModVariable

    type GC_target
        integer             ::  iBlock
        integer             ::  iRank
        integer             ::  nGC
        integer,allocatable ::  ijk_list(:,:)
        real,allocatable    ::  xijk_list(:,:)
        real,allocatable    ::  primitive_list(:,:)
    end type GC_target

    type Block
        integer             :: iBlock
        real                :: xijk_range(3,2)            ! range of the whole domain
        real                :: xijk_range_GC(3,2)         ! range including GC

        integer             :: nvar,ni,nj,nk,ng           ! variables & grid
        real                :: dxi,dxj,dxk                ! grid cell size
        real,allocatable    :: xi(:),xj(:),xk(:)          ! grid

        integer             :: rho1_,vi1_,vj1_,vk1_,s1_   ! positions of the vars

        real,allocatable    :: primitive(:,:,:,:)   ,&    ! the primitive variables
                            primitive_k2(:,:,:,:),&    ! |
                            primitive_k3(:,:,:,:),&    ! | used for runge kutta
                            primitive_k4(:,:,:,:)      ! |

        integer,allocatable         :: GC_iBlocks(:,:,:)  ! the iBlocks of EVERY grid (including non-GC)
        type(GC_target),allocatable :: GC_targets(:)      ! GC targets 
        type(GC_target),allocatable :: GC_sources(:)      ! GC targets 

        ! This is important:
        !
        ! GC_sources only contain those in the other ranks
        ! for those in the same rank no need to store it.

        integer                     :: nGC_targets,nGC_sources
    end type Block

    contains

    ! initiate a single block
    !
    ! initiate the grid, variables sort of things

    subroutine ModBlock_Init(Block1,iBlock,xijk_range,ni,nj,nk,ng,NameEqn,rk_order)
        implicit none
        type(Block),target          :: Block1           ! the block to be initiated
        integer,intent(in)          :: iBlock           ! Global iBlock
        integer,intent(in)          :: ni,nj,nk,ng      ! grid
        real,intent(in)             :: xijk_range(3,2)   ! range
        character(len=*),intent(in) :: NameEqn          ! which equation to be used
        integer,intent(in)          :: rk_order         ! runge kutta order

        integer                     :: nvar             ! how many variables

        nvar = ModVariable_SetNvar(NameEqn)             ! get how many variables

        Block1%iBlock=iBlock

        call ModBlock_InitGrid(Block1,xijk_range,nvar,ni,nj,nk,ng,rk_order)
        call ModBlock_InitVarInfo(Block1,NameEqn)
    end subroutine ModBlock_Init

    ! initiage the grid
    !
    !
    
    subroutine ModBlock_InitGrid(Block1,xijk_range,nvar,ni,nj,nk,ng,rk_order)
        implicit none
        type(Block),target :: Block1                        ! the block
        integer,intent(in) :: nvar,ni,nj,nk,ng,rk_order     ! gird size
        real,intent(in)    :: xijk_range(3,2)               ! range

        integer            :: i,j,k                         ! indices

        Block1%xijk_range=xijk_range                        ! set range
        Block1%ng=ng                                        ! |
        Block1%ni=ni                                        ! | set sizes
        Block1%nj=nj                                        ! |
        Block1%nk=nk                                        ! |
        Block1%nvar=nvar

        allocate(Block1%xi(-ng+1:ni+ng))                    ! |
        allocate(Block1%xj(-ng+1:nj+ng))                    ! | allocate grids
        allocate(Block1%xk(-ng+1:nk+ng))                    ! |

        select case(rk_order)     ! set the grid according to the runge kutta order
        case(1)
            allocate(Block1%primitive(1:nvar,-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
            Block1%primitive=0.
        case(2)
            allocate(Block1%primitive(1:nvar,-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
            allocate(Block1%primitive_k2(1:nvar,-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
            Block1%primitive=0.
            Block1%primitive_k2=0.
        case(4)
            allocate(Block1%primitive(1:nvar,-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
            allocate(Block1%primitive_k2(1:nvar,-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
            allocate(Block1%primitive_k3(1:nvar,-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
            allocate(Block1%primitive_k4(1:nvar,-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))
            Block1%primitive=0.
            Block1%primitive_k2=0.
            Block1%primitive_k3=0.
            Block1%primitive_k4=0.
        end select
        
        allocate(Block1%GC_iBlocks(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng))  ! initiate GC_iBlocks
        Block1%GC_iBlocks=-1                                              !

        do i=1-ng,ni+ng; Block1%xi(i)=(xijk_range(1,1)*(ni-i+0.5)+xijk_range(1,2)*(i-0.5))/ni; end do
        do j=1-ng,nj+ng; Block1%xj(j)=(xijk_range(2,1)*(nj-j+0.5)+xijk_range(2,2)*(j-0.5))/nj; end do
        do k=1-ng,nk+ng; Block1%xk(k)=(xijk_range(3,1)*(nk-k+0.5)+xijk_range(3,2)*(k-0.5))/nk; end do
        !do k=1-ng,nk+ng; Block1%xk(k)=(xijk_range(3,1)*(nk-k)+xijk_range(3,2)*(k-1))/(nk-1); end do

        Block1%dxi=(xijk_range(1,2)-xijk_range(1,1))/ni
        Block1%dxj=(xijk_range(2,2)-xijk_range(2,1))/nj
        Block1%dxk=(xijk_range(3,2)-xijk_range(3,1))/nk

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

    subroutine ModBlock_InitVarInfo(Block1,NameEqn)
        implicit none

        type(Block),target :: Block1
        character(len=*) :: NameEqn

        Block1%rho1_=ModVariable_SetVarIndex(NameEqn,'rho1')
        Block1%vi1_=ModVariable_SetVarIndex(NameEqn,'vi')
        Block1%vj1_=ModVariable_SetVarIndex(NameEqn,'vj')
        Block1%vk1_=ModVariable_SetVarIndex(NameEqn,'vk')
        Block1%s1_=ModVariable_SetVarIndex(NameEqn,'s1')
    end subroutine ModBlock_InitVarInfo

    subroutine ModBlock_deallocate(Block1)
        implicit none
        type(Block),intent(inout)   :: Block1
        integer                     :: iGC_target,iGC_source
        if (allocated(Block1%xi)) deallocate(Block1%xi)
        if (allocated(Block1%xj)) deallocate(Block1%xj)
        if (allocated(Block1%xk)) deallocate(Block1%xk)

        if (allocated(Block1%primitive)) deallocate(Block1%primitive)
        if (allocated(Block1%primitive_k2)) deallocate(Block1%primitive_k2)
        if (allocated(Block1%primitive_k3)) deallocate(Block1%primitive_k3)
        if (allocated(Block1%primitive_k4)) deallocate(Block1%primitive_k4)

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
    end subroutine ModBlock_deallocate
end module ModBlock