module ModEquation

    use ieee_arithmetic
    use ModBlock,       only:   BlockType
    use ModSpherical,   only:   ModSpherical_Grad_f,&
                                ModSpherical_A_dot_Grad_f,&
                                ModSpherical_div,&
                                ModSpherical_A_dot_nabla_B,&
                                ModSpherical_cross,&
                                ModSpherical_curl
    use ModBoundary,    only:   ModBoundary_Dynamo_HD_primitives,&
                                ModBoundary_Dynamo_HD_p1
    use ModDiffusion,   only:   ModDiffusion_Aritificial_1
    use ModParameters,  only:   ni,nj,nk,ng,nvar,ModelS_delta,ModelS_heating_ratio
    use ModVariables,   only:   rho1_,vr_,vt_,vp_,br_,bt_,bp_,s1_

    contains 

    ! The dynamo model from Hotta 2014.
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
            primitive=>Block1%primitive_rk
        else
            primitive=>Block1%primitive
        end if
        
        ! Perparations
        ! Get m; Set R=0.; Set p1; Set Boundary.
        EQN_update_R=0.
        
        Block1%p1=Block1%Gamma1*Block1%p0_over_rho0*primitive(:,:,:,rho1_)+&
            Block1%Gamma3_minus_1*Block1%rho0T0*primitive(:,:,:,s1_)

        call ModBoundary_Dynamo_HD_primitives(Block1,if_rk)
        !call ModBoundary_Dynamo_HD_p1(Block1)
        
        ! EQN rho1_
        do ivar=vr_,vp_
            tmp(:,:,:,ivar)=Block1%rho0*primitive(:,:,:,ivar)
        end do
        EQN_update_R(:,:,:,rho1_)=-1.0/(Block1%Xi_rsst(1:ni,1:nj,1:nk)**2*ModelS_delta)*&
            ModSpherical_div(ni,nj,nk,ng,Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,tmp)
        
        ! EQN vr_:vp_ Inertial Force
        EQN_update_R(:,:,:,vr_:vp_)=-&
            ModSpherical_A_dot_nabla_B(ni,nj,nk,ng,&
            Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,&
            primitive(:,:,:,vr_:vp_),primitive(:,:,:,vr_:vp_))
        
        ! EQN vr_:vp_ pressure gradient
        tmp(1:ni,1:nj,1:nk,:)=&
            ModSpherical_Grad_f(ni,nj,nk,ng,&
            Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,Block1%p1)
        do ivar=vr_,vp_
            EQN_update_R(:,:,:,ivar)=EQN_update_R(:,:,:,ivar)-&
                tmp(1:ni,1:nj,1:nk,ivar)/&
                Block1%rho0(1:ni,1:nj,1:nk)
        end do
        
        ! EQN vr Gravity
        EQN_update_R(:,:,:,vr_)=EQN_update_R(:,:,:,vr_)-&
            Block1%g_over_rho0*primitive(1:ni,1:nj,1:nk,rho1_)

        ! EQN s1_ advection term of s1
        EQN_update_R(:,:,:,s1_)=-&
            ModSpherical_A_dot_Grad_f(ni,nj,nk,ng,&
            Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,&
            primitive(:,:,:,vr_:vp_),primitive(:,:,:,s1_))
        
        ! EQN s1_ heating
        EQN_update_R(:,:,:,s1_)=EQN_update_R(:,:,:,s1_)+&
            ModelS_heating_ratio*(Block1%total_heat(1:ni,1:nj,1:nk))/Block1%rho0T0(1:ni,1:nj,1:nk)

        
    end subroutine ModEquation_Dynamo_HD

    subroutine ModEquation_Dynamo_MHD(Block1,if_rk,EQN_update_R)
        implicit none
        type(BlockType),target      ::  Block1                      
        logical,intent(in)          ::  if_rk
        real,pointer                ::  primitive(:,:,:,:)
        integer                     ::  ivar
        real,intent(out)            ::  EQN_update_R(1:ni,1:nj,1:nk,1:nvar)
        real                        ::  tmp(-ng+1:ng+ni,-ng+1:ng+nj,-ng+1:ng+nk,vr_:vp_)
        real                        ::  DivB(1:ni,1:nj,1:nk)
        
        ! If_rk
        if (if_rk) then
            primitive=>Block1%primitive_rk
        else
            primitive=>Block1%primitive
        end if
        
        ! Perparations
        ! Get m; Set R=0.; Set p1; Set Boundary.
        EQN_update_R=0.
        
        Block1%p1=Block1%Gamma1*Block1%p0_over_rho0*primitive(:,:,:,rho1_)+&
            Block1%Gamma3_minus_1*Block1%rho0T0*primitive(:,:,:,s1_)

        call ModBoundary_Dynamo_HD_primitives(Block1,if_rk)
        !call ModBoundary_Dynamo_HD_p1(Block1)
        
        ! EQN rho1_
        do ivar=vr_,vp_
            tmp(:,:,:,ivar)=Block1%rho0*primitive(:,:,:,ivar)
        end do
        EQN_update_R(:,:,:,rho1_)=-1.0/(Block1%Xi_rsst(1:ni,1:nj,1:nk)**2*ModelS_delta)*&
            ModSpherical_div(ni,nj,nk,ng,Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,tmp)
        
        ! EQN vr_:vp_ Inertial Force
        EQN_update_R(:,:,:,vr_:vp_)=-&
            ModSpherical_A_dot_nabla_B(ni,nj,nk,ng,&
            Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,&
            primitive(:,:,:,vr_:vp_),primitive(:,:,:,vr_:vp_))
        
        ! EQN vr_:vp_ pressure gradient & Lorentz force
        tmp(1:ni,1:nj,1:nk,:)=&
            -ModSpherical_Grad_f(ni,nj,nk,ng,&
            Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,Block1%p1)
        tmp(1:ni,1:nj,1:nk,:)=tmp(1:ni,1:nj,1:nk,:)+                        &
            ModSpherical_cross(ni,nj,nk,0,                                  &
                ModSpherical_curl(ni,nj,nk,ng,                              &
                    Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,   &
                    primitive(:,:,:,Br_:Bp_)),                              &
                primitive(1:ni,1:nj,1:nk,Br_:Bp_))
        do ivar=vr_,vp_
            EQN_update_R(:,:,:,ivar)=EQN_update_R(:,:,:,ivar)+&
                tmp(1:ni,1:nj,1:nk,ivar)/&
                Block1%rho0(1:ni,1:nj,1:nk)
        end do
        
        ! EQN vr Gravity
        EQN_update_R(:,:,:,vr_)=EQN_update_R(:,:,:,vr_)-&
            Block1%g_over_rho0*primitive(1:ni,1:nj,1:nk,rho1_)
            
        ! EQN s1_ advection term of s1
        EQN_update_R(:,:,:,s1_)=-&
            ModSpherical_A_dot_Grad_f(ni,nj,nk,ng,&
            Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,&
            primitive(:,:,:,vr_:vp_),primitive(:,:,:,s1_))
        
        ! EQN s1_ heating
        EQN_update_R(:,:,:,s1_)=EQN_update_R(:,:,:,s1_)+&
            ModelS_heating_ratio*(Block1%total_heat(1:ni,1:nj,1:nk))/Block1%rho0T0(1:ni,1:nj,1:nk)

        ! Magnetic field
        EQN_update_R(:,:,:,Br_:Bp_)=&
            ModSpherical_curl(ni,nj,nk,ng,&
            Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,&
            ModSpherical_cross(ni,nj,nk,ng,primitive(:,:,:,vr_:vp_),primitive(:,:,:,Br_:Bp_)))

        ! Div B correction
        DivB=ModSpherical_div(ni,nj,nk,ng,Block1%xi,Block1%xj,&
            Block1%dxi,Block1%dxj,Block1%dxk,primitive(:,:,:,br_:bp_))
        
        do ivar=vr_,vp_
            EQN_update_R(:,:,:,ivar)=EQN_update_R(:,:,:,ivar)-DivB*primitive(1:ni,1:nj,1:nk,br_+ivar-vr_)/Block1%rho0(1:ni,1:nj,1:nk)
            EQN_update_R(:,:,:,br_+ivar-vr_)=EQN_update_R(:,:,:,br_+ivar-vr_)-DivB*primitive(1:ni,1:nj,1:nk,ivar)
        end do
    end subroutine ModEquation_Dynamo_MHD
end module ModEquation