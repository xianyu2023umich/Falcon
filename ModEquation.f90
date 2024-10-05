Module ModEquation

    use ieee_arithmetic
    use ModBlock
    use ModSSM_v0
    use ModDeviation

    contains 

    ! The dynamo model from Fan 2003.

    subroutine ModEquation_Dynamo_HD_test
        
        integer :: ni,nj,nk,ng
        real :: dxi,dxj,dxk
        real :: primitive(5,-1:12,-1:12,-1:13),EQN_update_R(5,10,10,11)
        real :: xk(-1:12)
        integer :: i

        ni=10; nj=10; nk=11; ng=2
        dxi=0.1; dxj=0.1; dxk=0.1

        do i=-1,12
            xk(i)=(i-1)*dxk
            primitive(:,:,:,i)=1.
        end do

        call ModEquation_Dynamo_HD(primitive,&
            ni,nj,nk,ng,dxi,dxj,dxk,xk,'cartesian',EQN_update_R,.True.,.True.)

        print *,'last'
        !print *,EQN_update_R(1,5,5,1:nk)
        !print *,primitive(1,5,5,-1:12)
    end subroutine 

    subroutine ModEquation_Dynamo_HD(primitive,ni,nj,nk,ng,dxi,dxj,dxk,xk,geometry,EQN_update_R,&
            if_top,if_bottom)

        use ModParameter, only : paraXi,paraGamma,paraDelta_r,paraRe,paraPr

        implicit none

        character(len=*),intent(in) :: geometry
        integer,intent(in) :: ni,nj,nk,ng
        real,intent(in) :: dxi,dxj,dxk
        real,intent(in) :: xk(-ng+1:nk+ng)
        real,intent(inout) :: primitive(5,-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        logical,intent(in) :: if_top,if_bottom

        real :: p0(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        real :: p1(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        real :: te0(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        real :: rho0(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        real :: divergence_v(1:ni,1:nj,1:nk)
        real :: shear(1:3,1:3,1:ni,1:nj,1:nk)
        real :: rho0_list(-ng+1:nk+ng),p0_list(-ng+1:nk+ng),te0_list(-ng+1:nk+ng)
        
        integer :: i,j,k,direction1,direction2
        real :: m
        
        real,intent(out) :: EQN_update_R(1:5,1:ni,1:nj,1:nk)

        ! set the ghost cells for upper and lower boundaries

        if (if_top) then
            do i=nk+1,nk+ng
                primitive([1,2,3],:,:,i)=(primitive([1,2,3],:,:,nk-1)*8.-primitive([1,2,3],:,:,nk-2))/7.
                primitive(4,:,:,i )=0.!-primitive(4,:,:,2*nk-i)
                primitive(4,:,:,nk)=0.
                primitive(5,:,:,i )=0.!-primitive(5,:,:,2*nk-i)
                primitive(5,:,:,nk)=0.
            end do
        end if

        if (if_bottom) then
            do i=-ng+1,0
                primitive([1,2,3],:,:,i)=(primitive([1,2,3],:,:,nk-2)*8.-primitive([1,2,3],:,:,nk-3))/7.
                primitive(4,:,:,i)=0.!-primitive(4,:,:,2-i)
                primitive(4,:,:,1)=0.
                primitive(5,:,:,i)=0.!-primitive(5,:,:,2-i)
                primitive(5,:,:,1)=0.
            end do
        end if

        ! Perparations

        m=1./(paraGamma-1.)

        EQN_update_R(1:5,1:ni,1:nj,1:nk)=0.


        select case(geometry)
        case('cartesian')

            ! preparations

            ! get the rho p te profiles from the standard solar model

            rho0_list=ModSSM_v0_get_var0(nk+2*ng,xk,m,'rho0')
            p0_list=ModSSM_v0_get_var0(nk+2*ng,xk,m,'p0')
            te0_list=ModSSM_v0_get_var0(nk+2*ng,xk,m,'te0')

            ! and copy them to 3D domain

            do i=-ng+1,ni+ng; do j=-ng+1,nj+ng
                rho0(i,j,:)=rho0_list; p0(i,j,:)=p0_list; te0(i,j,:)=te0_list
            end do; end do

            ! p1 = gamma * p0 * (rho1/rho0 + s1) (non-dimensionalized)
            
            p1=paraGamma*p0*(primitive(1,:,:,:)/rho0+primitive(5,:,:,:))

            ! shear = partial vi / partial xj as a tensor.

            do direction1=1,3
                do direction2=1,3
                    shear(direction1,direction2,:,:,:)=&
                        ModDeviation_1st_O4(primitive(direction1+1,:,:,:),ni,nj,nk,ng,dxi,dxj,dxk,direction2,geometry)
                end do
            end do

            ! div v = shear_ii

            divergence_v=shear(1,1,:,:,:)+shear(2,2,:,:,:)+shear(3,3,:,:,:)
            
            ! end preparations.
            ! now starting equations.

            ! =============== Eqn 1 =================
            
            EQN_update_R(1,:,:,:)=(-1./paraXi**2)*(1./paraDelta_r/8.)*(&
                ModDeviation_1st_O4(rho0*primitive(2,:,:,:),ni,nj,nk,ng,dxi,dxj,dxk,1,geometry)+&
                ModDeviation_1st_O4(rho0*primitive(3,:,:,:),ni,nj,nk,ng,dxi,dxj,dxk,2,geometry)+&
                ModDeviation_1st_O4(rho0*primitive(4,:,:,:),ni,nj,nk,ng,dxi,dxj,dxk,3,geometry))
            
            ! ============== Eqn 2-5 ================  
            
            do direction1=1,3

                

                ! first term: inertial force
                
                EQN_update_R(direction1+1,:,:,:)=-&
                    primitive(2,1:ni,1:nj,1:nk)*&
                    ModDeviation_1st_O4(primitive(direction1+1,:,:,:),ni,nj,nk,ng,dxi,dxj,dxk,1,geometry)-&
                    primitive(3,1:ni,1:nj,1:nk)*&
                    ModDeviation_1st_O4(primitive(direction1+1,:,:,:),ni,nj,nk,ng,dxi,dxj,dxk,2,geometry)-&
                    primitive(4,1:ni,1:nj,1:nk)*&
                    ModDeviation_1st_O4(primitive(direction1+1,:,:,:),ni,nj,nk,ng,dxi,dxj,dxk,3,geometry)

                ! third term: gravity
                
                if (direction1 .eq. 3) then
                    EQN_update_R(direction1+1,:,:,:)=EQN_update_R(direction1+1,:,:,:)-primitive(1,1:ni,1:nj,1:nk)/&
                    rho0(1:ni,1:nj,1:nk)
                end if

                if (direction2==3) then
                    EQN_update_R(direction1+1,:,:,:)=EQN_update_R(direction1+1,:,:,:)+&
                        1./(rho0(1:ni,1:nj,1:nk)*ParaRe)*&
                        ModDeviation_1st_O4(rho0,ni,nj,nk,ng,dxi,dxj,dxk,direction2,geometry)*&
                        (shear(direction1,direction2,:,:,:)+shear(direction2,direction1,:,:,:))
                end if

                ! forth term: viscos force

                do direction2=1,3
                    
                    if (direction1==direction2) then
                        EQN_update_R(direction1+1,:,:,:)=EQN_update_R(direction1+1,:,:,:)+1./(3.*ParaRe)*&
                            ModDeviation_2nd_O3(primitive(direction1+1,:,:,:),&
                            ni,nj,nk,ng,dxi,dxj,dxk,direction1,geometry)
                    else
                        EQN_update_R(direction1+1,:,:,:)=EQN_update_R(direction1+1,:,:,:)+1./(3.*ParaRe)*&
                            ModDeviation_2D_O4(primitive(direction1+1,:,:,:),&
                            ni,nj,nk,ng,dxi,dxj,dxk,direction1,direction2,geometry)
                    end if
                    
                    

                    EQN_update_R(direction1+1,:,:,:)=EQN_update_R(direction1+1,:,:,:)+1./ParaRe*&
                        ModDeviation_2nd_O3(primitive(direction1+1,:,:,:),&
                        ni,nj,nk,ng,dxi,dxj,dxk,direction2,geometry)

                    if (direction2==3) then
                        EQN_update_R(direction1+1,:,:,:)=EQN_update_R(direction1+1,:,:,:)+&
                            1./(rho0(1:ni,1:nj,1:nk)*ParaRe)*&
                            ModDeviation_1st_O4(rho0,ni,nj,nk,ng,dxi,dxj,dxk,direction2,geometry)*&
                            (shear(direction1,direction2,:,:,:)+shear(direction2,direction1,:,:,:))
                    end if
                    
                end do

                if (direction1==3) then
                    EQN_update_R(direction1+1,:,:,:)=EQN_update_R(direction1+1,:,:,:)-&
                        2./(3.*rho0(1:ni,1:nj,1:nk)*ParaRe)*divergence_v*&
                        ModDeviation_1st_O4(rho0,ni,nj,nk,ng,dxi,dxj,dxk,direction1,geometry)
                end if

                ! viscous term ends
            end do

            
            
            ! boundary. only affect p1

            if (if_top) then
                do k=nk+1,nk+ng
                    p1(1:ni,1:nj,k)=&
                        (EQN_update_R(4,1:ni,1:nj,nk)*12.*dxk*rho0(1:ni,1:nj,nk)+8.*p1(1:ni,1:nj,nk-1)-p1(1:ni,1:nj,nk-2))/7.
                end do
            end if

            if (if_bottom) then
                do k=-ng+1,0
                    p1(1:ni,1:nj,k)=&
                        (-EQN_update_R(4,1:ni,1:nj,1)*12.*dxk*rho0(1:ni,1:nj,1)+8.*p1(1:ni,1:nj,2)-p1(1:ni,1:nj,3))/7.
                end do
            end if

            do direction1=1,3 ! pressure gradient in the momentum equations.
                EQN_update_R(direction1+1,:,:,:)=EQN_update_R(direction1+1,:,:,:)-&
                    ModDeviation_1st_O4(p1,ni,nj,nk,ng,dxi,dxj,dxk,direction1,geometry)/&
                    rho0(1:ni,1:nj,1:nk)
            end do

            ! ============== Eqn 5 ================


            ! first advection term

            EQN_update_R(5,:,:,:)=-&
                primitive(2,1:ni,1:nj,1:nk)*ModDeviation_1st_O4(primitive(5,:,:,:),ni,nj,nk,ng,dxi,dxj,dxk,1,geometry)-&
                primitive(3,1:ni,1:nj,1:nk)*ModDeviation_1st_O4(primitive(5,:,:,:),ni,nj,nk,ng,dxi,dxj,dxk,2,geometry)-&
                primitive(4,1:ni,1:nj,1:nk)*ModDeviation_1st_O4(primitive(5,:,:,:),ni,nj,nk,ng,dxi,dxj,dxk,3,geometry)

            ! second advection term

            EQN_update_R(5,:,:,:)=EQN_update_R(5,:,:,:)+primitive(4,1:ni,1:nj,1:nk)/&
                p0(1:ni,1:nj,1:nk)/8.

            ! heat conduction
            
            do direction1=1,2
                EQN_update_R(5,:,:,:)=EQN_update_R(5,:,:,:)+1./paraRe/paraPr*&
                    ModDeviation_2nd_O3(primitive(5,:,:,:),ni,nj,nk,ng,dxi,dxj,dxk,direction1,geometry)
            end do
            EQN_update_R(5,:,:,:)=EQN_update_R(5,:,:,:)+1./paraRe/paraPr*&
                ModDeviation_2nd_O3(rho0*te0*primitive(5,:,:,:),ni,nj,nk,ng,dxi,dxj,dxk,3,geometry)/&
                (te0(1:ni,1:nj,1:nk)*rho0(1:ni,1:nj,1:nk))

            ! viscous heating
            
            

            do direction1=1,3
                do direction2=1,3
                    EQN_update_R(5,:,:,:)=EQN_update_R(5,:,:,:)+1./(paraRe*paraGamma)*(paraGamma-1)*&
                        rho0(1:ni,1:nj,1:nk)/p0(1:ni,1:nj,1:nk)*&
                        shear(direction1,direction2,:,:,:)*&
                        (shear(direction1,direction2,:,:,:)+shear(direction2,direction1,:,:,:))
                end do
            end do

            EQN_update_R(5,:,:,:)=EQN_update_R(5,:,:,:)-2./(3.*paraRe*paraGamma)*(paraGamma-1)*&
                rho0(1:ni,1:nj,1:nk)/p0(1:ni,1:nj,1:nk)*divergence_v**2
            
        end select
        
    end subroutine ModEquation_Dynamo_HD



end Module ModEquation