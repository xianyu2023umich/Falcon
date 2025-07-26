module ModEquation

    use ModBlock,       only:   BlockType
    use ModSpherical,   only:   ModSpherical_Grad_f,&
                                ModSpherical_A_dot_Grad_f,&
                                ModSpherical_div,&
                                ModSpherical_A_dot_nabla_B,&
                                ModSpherical_cross,&
                                ModSpherical_curl
    use ModDiffusion,    only:   ModDiffusion_Aritificial_1
    use ModBoundary,     only:   ModBoundary_Dynamo_HD_primitives

    use ModParameters,  only:   ni,nj,nk,ng,nvar,ModelS_delta,ModelS_heating_ratio
    use ModVariables,   only:   rho1_,vr_,vp_,s1_

    implicit none

    ! Coeff1 is 2*\bar{\rho} * k_{B} * \bar{T} / m_{p}. p=2nkBT.
    real                            ::  Coeff1=1.0

    contains

    subroutine ModEquation_Dynamo_HD(Block1,if_rk,EQN_update_R)
        implicit none
        type(BlockType),target      ::  Block1                      
        logical,intent(in)          ::  if_rk
        real,pointer                ::  primitive(:,:,:,:)
        integer                     ::  ivar
        real,intent(out)            ::  EQN_update_R(1:ni,1:nj,1:nk,1:nvar)
        real                        ::  tmp(-ng+1:ng+ni,-ng+1:ng+nj,-ng+1:ng+nk,vr_:vp_)
        
        ! If_rk
        if (if_rk) then
            primitive=>Block1%primitive_rk_IV
        else
            primitive=>Block1%primitive_IV
        end if
        
        ! Perparations
        ! Get m; Set R=0.; Set p1; Set Boundary.
        EQN_update_R=0.

        call ModBoundary_Dynamo_HD_primitives(Block1,if_rk)

        Block1%p1_III=Block1%gamma1_III*Block1%p0_over_rho0_III*primitive(:,:,:,rho1_)+&
            Block1%gamma3_minus_1_III*Block1%rho0T0_III*primitive(:,:,:,s1_)
        !call ModBoundary_Dynamo_HD_p1(Block1)
        
        ! EQN rho1_
        do ivar=vr_,vp_
            tmp(:,:,:,ivar)=Block1%rho0_III*primitive(:,:,:,ivar)
        end do
        EQN_update_R(:,:,:,rho1_)=-1.0/(Block1%Xi_rsst_III(1:ni,1:nj,1:nk)**2*ModelS_delta)*&
            ModSpherical_div(ni,nj,nk,ng,Block1%xi_I,Block1%xj_I,Block1%dxi,Block1%dxj,Block1%dxk,tmp)
        
        ! EQN vr_:vp_ Inertial Force
        EQN_update_R(:,:,:,vr_:vp_)=-&
            ModSpherical_A_dot_nabla_B(ni,nj,nk,ng,&
            Block1%xi_I,Block1%xj_I,Block1%dxi,Block1%dxj,Block1%dxk,&
            primitive(:,:,:,vr_:vp_),primitive(:,:,:,vr_:vp_))
        
        ! EQN vr_:vp_ pressure gradient
        tmp(1:ni,1:nj,1:nk,:)=&
            ModSpherical_Grad_f(ni,nj,nk,ng,&
            Block1%xi_I,Block1%xj_I,Block1%dxi,Block1%dxj,Block1%dxk,Block1%p1_III)
        do ivar=vr_,vp_
            EQN_update_R(:,:,:,ivar)=EQN_update_R(:,:,:,ivar)-&
                tmp(1:ni,1:nj,1:nk,ivar)/&
                Block1%rho0_III(1:ni,1:nj,1:nk)
        end do
        
        ! EQN vr Gravity
        EQN_update_R(:,:,:,vr_)=EQN_update_R(:,:,:,vr_)-&
            Block1%g_over_rho0_III*primitive(1:ni,1:nj,1:nk,rho1_)

        ! EQN s1_ advection term of s1
        EQN_update_R(:,:,:,s1_)=-&
            ModSpherical_A_dot_Grad_f(ni,nj,nk,ng,&
            Block1%xi_I,Block1%xj_I,Block1%dxi,Block1%dxj,Block1%dxk,&
            primitive(:,:,:,vr_:vp_),primitive(:,:,:,s1_))
        
        ! EQN s1_ heating
        EQN_update_R(:,:,:,s1_)=EQN_update_R(:,:,:,s1_)+&
            ModelS_heating_ratio*(Block1%total_heat_III(1:ni,1:nj,1:nk))/Block1%rho0T0_III(1:ni,1:nj,1:nk)

        call ModDiffusion_Aritificial_1(Block1,EQN_update_R,2,if_rk)
    end subroutine ModEquation_Dynamo_HD
    
    subroutine ModEquation_CORONA_V1(Block1,if_rk,EQN_update_R)
        implicit none
        type(BlockType),target      ::  Block1                      
        logical,intent(in)          ::  if_rk
        real,pointer                ::  primitive(:,:,:,:)

        real                        ::  gamma=5./3.
        integer                     ::  ivar
        real,intent(out)            ::  EQN_update_R(1:ni,1:nj,1:nk,1:nvar)
        real                        ::  tmp(-ng+1:ng+ni,-ng+1:ng+nj,-ng+1:ng+nk,2:4)
        real                        ::  p(-ng+1:ng+ni,-ng+1:ng+nj,-ng+1:ng+nk)

        ! If_rk
        if (if_rk) then
            primitive=>Block1%primitive_rk_IV
        else
            primitive=>Block1%primitive_IV
        end if

        ! Here should be a boundary.

        ! Preparation
        p=Coeff1*primitive(:,:,:,1)*primitive(:,:,:,8)
        EQN_update_R=0.

        ! The rho equation
        do ivar=2,4
            tmp(:,:,:,ivar)=primitive(:,:,:,1)*primitive(:,:,:,ivar)
        end do

        EQN_update_R(:,:,:,1)=-&
            ModSpherical_div(ni,nj,nk,ng,Block1%xi_I,Block1%xj_I,Block1%dxi,Block1%dxj,Block1%dxk,tmp)

        ! The v_i equation
        ! Inertial force
        EQN_update_R(:,:,:,2:4)=-&
            ModSpherical_A_dot_nabla_B(ni,nj,nk,ng,&
            Block1%xi_I,Block1%xj_I,Block1%dxi,Block1%dxj,Block1%dxk,&
            primitive(:,:,:,2:4),primitive(:,:,:,2:4))
        
        ! Gas pressure gradient + Lorentz force
        tmp(1:ni,1:nj,1:nk,:)=&
            ModSpherical_Grad_f(ni,nj,nk,ng,&
            Block1%xi_I,Block1%xj_I,Block1%dxi,Block1%dxj,Block1%dxk,p)
        tmp(1:ni,1:nj,1:nk,:)=tmp(1:ni,1:nj,1:nk,:)+                        &
            ModSpherical_cross(ni,nj,nk,0,                                  &
                ModSpherical_curl(ni,nj,nk,ng,                              &
                    Block1%xi_I,Block1%xj_I,Block1%dxi,Block1%dxj,Block1%dxk,   &
                    primitive(:,:,:,5:7)),                              &
                primitive(1:ni,1:nj,1:nk,5:7))

        do ivar=2,4
            EQN_update_R(:,:,:,ivar)=EQN_update_R(:,:,:,ivar)-&
                tmp(1:ni,1:nj,1:nk,ivar)/&
                primitive(1:ni,1:nj,1:nk,1)
        end do

        ! Gravity
        EQN_update_R(:,:,:,2)=EQN_update_R(:,:,:,2)-&
            primitive(1:ni,1:nj,1:nk,1)*Block1%g_over_rho0_III

        ! Magnetic induction equation

        ! Magnetic field
        EQN_update_R(:,:,:,5:7)=&
            ModSpherical_curl(ni,nj,nk,ng,&
            Block1%xi_I,Block1%xj_I,Block1%dxi,Block1%dxj,Block1%dxk,&
            ModSpherical_cross(ni,nj,nk,ng,primitive(:,:,:,2:4),primitive(:,:,:,5:7)))

        ! Temperature
        do ivar=2,4
            tmp(:,:,:,ivar)=primitive(:,:,:,ivar)*primitive(:,:,:,8)
        end do
        EQN_update_R(:,:,:,8)=-&
            ModSpherical_div(ni,nj,nk,ng,Block1%xi_I,Block1%xj_I,&
            Block1%dxi,Block1%dxj,Block1%dxk,tmp)
        EQN_update_R(:,:,:,8)=EQN_update_R(:,:,:,8)-&
            (gamma-2.0)*primitive(1:ni,1:nj,1:nk,8)*&
            ModSpherical_div(ni,nj,nk,ng,Block1%xi_I,Block1%xj_I,&
            Block1%dxi,Block1%dxj,Block1%dxk,primitive(:,:,:,2:4))

    end subroutine ModEquation_CORONA_V1
end module 