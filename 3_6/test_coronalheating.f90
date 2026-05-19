program test_coronalheating

    ! Use a fixed potential field with fixed density to see the evolution of the w_plus and w_minus

    use ModConst,           only:   dpi,m_p__CGS,R_sun__CGS
    use ModBlock,           only:   BlockType,ModBlock_Init
    use ModCoronalHeating,  only:   ModCoronalHeating_TurbulentCacade
    use ModParameters,      only:   ni,nj,nk,ng,CFL
    use ModMHD,             only:   ModMHD_AlfvenVelocityVector_3D
    use ModSpherical,       only:   ModSpherical_abs
    implicit none

    type(BlockType),target     ::  Block1
    real(8)             ::  xijk_range(3,2)
    real(8),allocatable ::  AlfvenVelocityVector_IV(:,:,:,:),Va_LLL(:,:,:)
    integer             ::  i,j,k
    real(8)             ::  monopole_charge,monopoles_xyz(3,2)
    real(8)             ::  r,theta,phi,x,y,z,dx,dy,dz
    real(8)             ::  r_monopoles(2)
    real(8)             ::  bxyz(3)
    real(8)             ::  dt
    integer             ::  istep
    logical             ::  if_rk_input,if_rk_output
    integer             ::  rk_index

    ni=100
    nj=1
    nk=1
    ng=2
    CFL=0.75

    allocate(AlfvenVelocityVector_IV(1:ni,1:nj,1:nk,3))
    allocate(Va_LLL(1:ni,1:nj,1:nk))

    xijk_range(1,:) = [1.0d0,2.0d0]*R_sun__CGS
    xijk_range(2,:) = [45.0d0,135.0d0]*dpi/180.0d0
    xijk_range(3,:) = [-45.0d0,45.0d0]*dpi/180.0d0

    call ModBlock_Init(Block1,1,xijk_range,.true.,2,.true.)
    Block1%primitive_IV=0.0d0

    ! assume coronal number density is uniformly 10^10 cm^-3

    Block1%primitive_IV(:,:,:,Block1%rho_)=1.0d10*m_p__CGS
    

    ! Initialize dipole field. The two monopoles are at [0.95,0,\pm 0.1] Rs.
    ! B=monopole_charge/r^3 * hat{r}
    
    !call Initialize_dipole
    call Initialize_tube
    print *,Block1%xi_I(-ng+1:ni+ng)

    ! Set boundary condition for w_plus and w_minus. Nothing else is updated so no need for BC.
    ! If Br>0 then w_plus=1.0e5, w_minus mirror, otherwise w_plus mirror, w_minus=1.0e5.

    do k=-ng+1,nk+ng
        do j=-ng+1,nj+ng
            do i=-ng+1,0
                if (Block1%primitive_IV(i,j,k,Block1%br_)>0.0d0) then
                    Block1%primitive_IV(i,j,k,Block1%w_plus_)=1.0e5
                    Block1%primitive_IV(i,j,k,Block1%w_minus_)=Block1%primitive_IV(1-i,j,k,Block1%w_minus_)
                else
                    Block1%primitive_IV(i,j,k,Block1%w_plus_)=Block1%primitive_IV(1-i,j,k,Block1%w_plus_)
                    Block1%primitive_IV(i,j,k,Block1%w_minus_)=1.0e5
                end if

            end do
        end do
    end do

    Block1%primitive_rk_IV=Block1%primitive_IV

    AlfvenVelocityVector_IV=&
            ModMHD_AlfvenVelocityVector_3D(ni,nj,nk,ng,&
            Block1%primitive_IV(:,:,:,Block1%rho_),&
            Block1%primitive_IV(:,:,:,Block1%br_:Block1%bp_))

    Va_LLL=ModSpherical_abs(ni,nj,nk,0,AlfvenVelocityVector_IV)
    dt=CFL*min(Block1%dxi,Block1%dxj*minval(Block1%xi_I),Block1%dxk*minval(Block1%xi_I))/maxval(Va_LLL)

    !call ModCoronalHeating_TurbulentCacade(Block1,.false.)

    !stop 1

    do istep=1,10000
        print *,istep
        do rk_index=1,4
            ! Decide if rk or ot for the input and output primitives
            if_rk_input=(rk_index>1)
            if_rk_output=(rk_index<4)

            call Boundary_Condition

            ! Get the EQN_update_R
            call ModCoronalHeating_TurbulentCacade(Block1,if_rk_input)

            ! Get the next RK
            if (if_rk_output) then
                Block1%primitive=>Block1%primitive_rk_IV
            else
                Block1%primitive=>Block1%primitive_IV
            end if  

            Block1%primitive(1:ni,1:nj,1:nk,Block1%w_plus_:Block1%w_minus_)=&
                Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%w_plus_:Block1%w_minus_)+dt*Block1%EQN_update_R_IV(1:ni,1:nj,1:nk,Block1%w_plus_:Block1%w_minus_)/(5.0d0-rk_index)
     
            !print *,Block1%primitive_rk_IV(:,1,1,Block1%w_plus_)  
        end do

        print *,Block1%primitive_IV(:,1,1,Block1%w_plus_)

    end do

    contains

    subroutine Boundary_Condition

        if (if_rk_input) then
            Block1%primitive=>Block1%primitive_rk_IV
        else
            Block1%primitive=>Block1%primitive_IV
        end if

        do k=-ng+1,nk+ng
            do j=-ng+1,nj+ng

                ! Bottom
                do i=-ng+1,0
                    if (Block1%primitive(i,j,k,Block1%br_)>0.0d0) then
                        Block1%primitive(i,j,k,Block1%w_plus_)=1.0e5
                        Block1%primitive(i,j,k,Block1%w_minus_)=Block1%primitive(1-i,j,k,Block1%w_minus_)
                    else
                        Block1%primitive(i,j,k,Block1%w_plus_)=Block1%primitive(1-i,j,k,Block1%w_plus_)
                        Block1%primitive(i,j,k,Block1%w_minus_)=1.0e5
                    end if
    
                end do

                ! Top

                do i=ni+1,ni+ng
                    if (Block1%primitive(i,j,k,Block1%br_)>0.0d0) then
                        Block1%primitive(i,j,k,Block1%w_plus_)=Block1%primitive(2*ni+1-i,j,k,Block1%w_plus_)
                        Block1%primitive(i,j,k,Block1%w_minus_)=Block1%primitive(2*ni+1-i,j,k,Block1%w_minus_)
                    else
                        Block1%primitive(i,j,k,Block1%w_plus_)=Block1%primitive(2*ni+1-i,j,k,Block1%w_plus_)
                        Block1%primitive(i,j,k,Block1%w_minus_)=Block1%primitive(2*ni+1-i,j,k,Block1%w_minus_)
                    end if
    
                end do
            end do
        end do

    end subroutine Boundary_Condition
    
    ! Straight tube
    ! div B=0 -> B=const/r^2

    subroutine Initialize_tube
        do k=-ng+1,nk+ng
            do j=-ng+1,nj+ng
                Block1%primitive_IV(:,j,k,Block1%br_)=10.0/(Block1%xi_I(-ng+1:ni+ng)/R_sun__CGS)**2
            end do
        end do
        
        Block1%primitive_IV(:,:,:,Block1%bt_)=0.0d0
        Block1%primitive_IV(:,:,:,Block1%bp_)=0.0d0
    end subroutine Initialize_tube

    subroutine Initialize_dipole
        monopole_charge=1.0d35
        monopoles_xyz(:,1)=[0.95,0.0,0.1]
        monopoles_xyz(:,2)=[0.95,0.0,-0.1]

        do k=-ng+1,nk+ng
            do j=-ng+1,nj+ng
                do i=-ng+1,ni+ng
                    r=Block1%xi_I(i)
                    theta=Block1%xj_I(j)
                    phi=Block1%xk_I(k)
                    x=r*sin(theta)*cos(phi)
                    y=r*sin(theta)*sin(phi)
                    z=r*cos(theta)
                    r=sqrt(x**2+y**2+z**2)
                    dx=x-monopoles_xyz(1,1)
                    dy=y-monopoles_xyz(2,1)
                    dz=z-monopoles_xyz(3,1)
                    r_monopoles(1)=sqrt(dx**2+dy**2+dz**2)
                    bxyz(1)=monopole_charge/r_monopoles(1)**3*dx
                    bxyz(2)=monopole_charge/r_monopoles(1)**3*dy
                    bxyz(3)=monopole_charge/r_monopoles(1)**3*dz


                    dx=x-monopoles_xyz(1,2)
                    dy=y-monopoles_xyz(2,2)
                    dz=z-monopoles_xyz(3,2)
                    r_monopoles(2)=sqrt(dx**2+dy**2+dz**2)
                    bxyz(1)=bxyz(1)-monopole_charge/r_monopoles(2)**3*dx
                    bxyz(2)=bxyz(2)-monopole_charge/r_monopoles(2)**3*dy
                    bxyz(3)=bxyz(3)-monopole_charge/r_monopoles(2)**3*dz

                    Block1%primitive_IV(i,j,k,Block1%br_)=sin(theta)*cos(phi)*bxyz(1)+sin(theta)*sin(phi)*bxyz(2)+cos(theta)*bxyz(3)
                    Block1%primitive_IV(i,j,k,Block1%bt_)=cos(theta)*cos(phi)*bxyz(1)+cos(theta)*sin(phi)*bxyz(2)-sin(theta)*bxyz(3)
                    Block1%primitive_IV(i,j,k,Block1%bp_)=-sin(phi)*bxyz(1)+cos(phi)*bxyz(2)
                end do
            end do
        end do
    end subroutine Initialize_dipole

    !print *,Block1%primitive_IV(:,1,1,Block1%w_plus_)


end program test_coronalheating