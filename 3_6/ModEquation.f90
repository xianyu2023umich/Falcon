module ModEquation

    use ModBlock,       only:   BlockType
    use ModSpherical,   only:   ModSpherical_Grad_f,&
                                ModSpherical_A_dot_Grad_f,&
                                ModSpherical_div,&
                                ModSpherical_A_dot_nabla_B,&
                                ModSpherical_cross,&
                                ModSpherical_curl
    use ModDiffusion,   only:   ModDiffusion_Aritificial_1
    use ModBoundary,    only:   ModBoundary_Dynamo_HD_primitives,&
                                ModBoundary_Dynamo_MHD_primitives

    use ModParameters,  only:   ni,nj,nk,ng,nvar,ModelS_heating_ratio
    use ModConst,       only:   gamma_ideal_gas

    implicit none

    contains
    

    ! Dynamo: first order fluctuations (e.g. rho1)
    ! Corona: Zeroth order variables (e.g. total rho)

    subroutine ModEquation_Dynamo_HD(Blc1,if_rk,EQN_update_R)
        implicit none
        type(BlockType),target      ::  Blc1                      
        logical,intent(in)          ::  if_rk
        real(8),intent(out),optional ::  EQN_update_R(:,:,:,:)
        
        ! If_rk
        if (if_rk) then
            Blc1%primitive=>Blc1%primitive_rk_IV
        else
            Blc1%primitive=>Blc1%primitive_IV
        end if
        
        ! Perparations
        ! Set R=0.; Set p1; Set Boundary.
        EQN_update_R=0.
        call ModBoundary_Dynamo_HD_primitives(Blc1,if_rk)
        call ModEquation_Dynamo_Get_p1(Blc1)
        
        ! Main equations
        call ModEquation_Dynamo_Mass_Conservation(Blc1)
        call ModEquation_Dynamo_Inertial_Force(Blc1)
        call ModEquation_Dynamo_Pressure_Gradient(Blc1)
        call ModEquation_Dynamo_Gravity(Blc1)
        call ModEquation_Dynamo_Entropy_Advection(Blc1)
        call ModEquation_Dynamo_Entropy_Heating(Blc1)
        if (present(EQN_update_R)) EQN_update_R = Blc1%EQN_update_R_IV(1:ni,1:nj,1:nk,:)
    end subroutine ModEquation_Dynamo_HD

    ! Linearized EOS: p1 = γ1*(p0/ρ0)*ρ1 + (γ3−1)*ρ0*T0*s1
    subroutine ModEquation_Dynamo_Get_p1(Blc1)
        implicit none
        type(BlockType),target      ::  Blc1

        Blc1%p1_III=Blc1%gamma1_III*Blc1%p0_over_rho0_III*Blc1%primitive(:,:,:,Blc1%rho1_)+&
            Blc1%gamma3_minus_1_III*Blc1%rho0T0_III*Blc1%primitive(:,:,:,Blc1%s1_)
    end subroutine ModEquation_Dynamo_Get_p1



    ! RSST continuity: ∂ρ1/∂t = −∇·(ρ0 v) / (ξ²_rsst * δ)
    ! The ξ² factor reduces the effective sound speed (Reduced Sound Speed Technique).
    subroutine ModEquation_Dynamo_Mass_Conservation(Blc1)
        implicit none
        type(BlockType),target      ::  Blc1
        integer                     ::  ivar
        real(8)                     ::  tmp(-ng+1:ng+ni,-ng+1:ng+nj,-ng+1:ng+nk,Blc1%vr_:Blc1%vp_)

        do ivar=Blc1%vr_,Blc1%vp_
            tmp(:,:,:,ivar)=Blc1%rho0_III*Blc1%primitive(:,:,:,ivar)
        end do
        Blc1%EQN_update_R_IV(:,:,:,Blc1%rho1_)=-1.0/Blc1%Xi_rsst_III(1:ni,1:nj,1:nk)**2*&
            ModSpherical_div(ni,nj,nk,ng,Blc1%xi_I,Blc1%xj_I,Blc1%dxi,Blc1%dxj,Blc1%dxk,tmp)
    end subroutine ModEquation_Dynamo_Mass_Conservation



    ! Advective acceleration: ∂v/∂t = −(v·∇)v
    subroutine ModEquation_Dynamo_Inertial_Force(Blc1)
        implicit none
        type(BlockType),target      ::  Blc1

        Blc1%EQN_update_R_IV(:,:,:,Blc1%vr_:Blc1%vp_)=-&
            ModSpherical_A_dot_nabla_B(ni,nj,nk,ng,&
            Blc1%xi_I,Blc1%xj_I,Blc1%dxi,Blc1%dxj,Blc1%dxk,&
            Blc1%primitive(:,:,:,Blc1%vr_:Blc1%vp_),Blc1%primitive(:,:,:,Blc1%vr_:Blc1%vp_))
    end subroutine ModEquation_Dynamo_Inertial_Force



    ! Pressure gradient: ∂v/∂t += −∇p1 / ρ0
    subroutine ModEquation_Dynamo_Pressure_Gradient(Blc1)
        implicit none
        type(BlockType),target      ::  Blc1
        integer                     ::  ivar
        real(8)                     ::  tmp(1:ni,1:nj,1:nk,Blc1%vr_:Blc1%vp_)

        tmp=ModSpherical_Grad_f(ni,nj,nk,ng,&
            Blc1%xi_I,Blc1%xj_I,Blc1%dxi,Blc1%dxj,Blc1%dxk,Blc1%p1_III)
        
        do ivar=Blc1%vr_,Blc1%vp_
            Blc1%EQN_update_R_IV(:,:,:,ivar)=Blc1%EQN_update_R_IV(:,:,:,ivar)-&
                tmp(:,:,:,ivar)/Blc1%rho0_III(1:ni,1:nj,1:nk)
        end do
    end subroutine ModEquation_Dynamo_Pressure_Gradient



    ! Lorentz force: ∂v/∂t += (∇×B)×B / ρ0  (i.e. J×B/ρ0 with J = ∇×B)
    subroutine ModEquation_Dynamo_Lorentz_Force(Blc1)
        implicit none
        type(BlockType),target      ::  Blc1
        integer                     ::  ivar
        real(8)                     ::  tmp(1:ni,1:nj,1:nk,Blc1%vr_:Blc1%vp_)
        real(8)                     ::  curl_B(1:ni,1:nj,1:nk,Blc1%br_:Blc1%bp_)

        curl_B=ModSpherical_curl(                                   &
            ni,nj,nk,ng,Blc1%xi_I,Blc1%xj_I,                    &
            Blc1%dxi,Blc1%dxj,Blc1%dxk,                       &
            Blc1%primitive(:,:,:,Blc1%br_:Blc1%bp_))

        tmp=ModSpherical_cross(ni,nj,nk,0,curl_B,                   &
                Blc1%primitive(1:ni,1:nj,1:nk,Blc1%br_:Blc1%bp_))
        
        do ivar=Blc1%vr_,Blc1%vp_
            Blc1%EQN_update_R_IV(:,:,:,ivar)=Blc1%EQN_update_R_IV(:,:,:,ivar)+&
                tmp(:,:,:,ivar)/Blc1%rho0_III(1:ni,1:nj,1:nk)
        end do
    end subroutine ModEquation_Dynamo_Lorentz_Force



    ! Buoyancy: ∂vr/∂t += −(g/ρ0) ρ1  (radial component only)
    subroutine ModEquation_Dynamo_Gravity(Blc1)
        implicit none
        type(BlockType),target      ::  Blc1

        Blc1%EQN_update_R_IV(:,:,:,Blc1%vr_)=Blc1%EQN_update_R_IV(:,:,:,Blc1%vr_)-&
            Blc1%g_over_rho0_III(1:ni,1:nj,1:nk)*Blc1%primitive(1:ni,1:nj,1:nk,Blc1%rho1_)
    end subroutine ModEquation_Dynamo_Gravity



    ! Entropy advection: ∂s1/∂t = −(v·∇)s1
    subroutine ModEquation_Dynamo_Entropy_Advection(Blc1)
        implicit none
        type(BlockType),target      ::  Blc1

        Blc1%EQN_update_R_IV(:,:,:,Blc1%s1_)=-&
            ModSpherical_A_dot_Grad_f(ni,nj,nk,ng,&
            Blc1%xi_I,Blc1%xj_I,Blc1%dxi,Blc1%dxj,Blc1%dxk,&
            Blc1%primitive(:,:,:,Blc1%vr_:Blc1%vp_),Blc1%primitive(:,:,:,Blc1%s1_))
    end subroutine ModEquation_Dynamo_Entropy_Advection



    ! Entropy source from heating: ∂s1/∂t += Q / (ρ0 T0)
    subroutine ModEquation_Dynamo_Entropy_Heating(Blc1)
        implicit none
        type(BlockType),target      ::  Blc1

        Blc1%EQN_update_R_IV(:,:,:,Blc1%s1_)=Blc1%EQN_update_R_IV(:,:,:,Blc1%s1_)+&
            ModelS_heating_ratio*(Blc1%total_heat_III(1:ni,1:nj,1:nk))/Blc1%rho0T0_III(1:ni,1:nj,1:nk)
    end subroutine ModEquation_Dynamo_Entropy_Heating

    ! Ideal induction: ∂B/∂t = ∇×(v×B)
    subroutine ModEquation_Dynamo_Induction_Equation(Blc1)
        implicit none
        type(BlockType),target      ::  Blc1

        Blc1%EQN_update_R_IV(:,:,:,Blc1%br_:Blc1%bp_)=&
            ModSpherical_curl(ni,nj,nk,ng,&
            Blc1%xi_I,Blc1%xj_I,Blc1%dxi,Blc1%dxj,Blc1%dxk,&
            ModSpherical_cross(ni,nj,nk,ng,Blc1%primitive(:,:,:,Blc1%vr_:Blc1%vp_),Blc1%primitive(:,:,:,Blc1%br_:Blc1%bp_)))
    end subroutine ModEquation_Dynamo_Induction_Equation

    subroutine ModEquation_Dynamo_MHD(Blc1,if_rk,EQN_update_R)
        implicit none
        type(BlockType),target      ::  Blc1                      
        logical,intent(in)          ::  if_rk
        real(8),intent(out),optional ::  EQN_update_R(:,:,:,:)
        integer                     ::  ivar
        real(8)                     ::  DivB(1:ni,1:nj,1:nk)
        
        ! If_rk
        if (if_rk) then
            Blc1%primitive=>Blc1%primitive_rk_IV
        else
            Blc1%primitive=>Blc1%primitive_IV
        end if
        
        ! Perparations
        ! Set R=0.; Set p1; Set Boundary.
        Blc1%EQN_update_R_IV=0.

        call ModBoundary_Dynamo_MHD_primitives(Blc1,if_rk)
        call ModEquation_Dynamo_Get_p1(Blc1)
        
        ! Main equations
        call ModEquation_Dynamo_Mass_Conservation(Blc1)
        call ModEquation_Dynamo_Inertial_Force(Blc1)
        call ModEquation_Dynamo_Pressure_Gradient(Blc1)
        call ModEquation_Dynamo_Lorentz_Force(Blc1)
        call ModEquation_Dynamo_Gravity(Blc1)
        call ModEquation_Dynamo_Entropy_Advection(Blc1)
        call ModEquation_Dynamo_Entropy_Heating(Blc1)
        call ModEquation_Dynamo_Induction_Equation(Blc1)
        
        ! Powell div-B source term: ∂B/∂t −= (∇·B) v  (removes divergence errors)
        DivB=ModSpherical_div(ni,nj,nk,ng,Blc1%xi_I,Blc1%xj_I,&
            Blc1%dxi,Blc1%dxj,Blc1%dxk,Blc1%primitive(:,:,:,Blc1%br_:Blc1%bp_))

        do ivar=Blc1%vr_,Blc1%vp_
            Blc1%EQN_update_R_IV(:,:,:,Blc1%br_+ivar-Blc1%vr_)=Blc1%EQN_update_R_IV(:,:,:,Blc1%br_+ivar-Blc1%vr_)-&
                DivB*Blc1%primitive(1:ni,1:nj,1:nk,ivar)
        end do
        if (present(EQN_update_R)) EQN_update_R = Blc1%EQN_update_R_IV(1:ni,1:nj,1:nk,:)
    end subroutine ModEquation_Dynamo_MHD






    subroutine ModEquation_Corona_Mass_Conservation(Blc1)
        implicit none
        type(BlockType),target      ::  Blc1
        integer                     ::  ivar
        real(8)                     ::  tmp(-ng+1:ng+ni,-ng+1:ng+nj,-ng+1:ng+nk,Blc1%vr_:Blc1%vp_)

        do ivar=Blc1%vr_,Blc1%vp_
            tmp(:,:,:,ivar)=Blc1%primitive(:,:,:,Blc1%rho_)*Blc1%primitive(:,:,:,ivar)
        end do
        Blc1%EQN_update_R_IV(:,:,:,Blc1%rho_)=-ModSpherical_div(ni,nj,nk,ng,&
            Blc1%xi_I,Blc1%xj_I,Blc1%dxi,Blc1%dxj,Blc1%dxk,tmp)
    end subroutine ModEquation_Corona_Mass_Conservation

    ! Inertial force is the same as the dynamo.

    subroutine ModEquation_Corona_Inertial_Force(Blc1)
        implicit none
        type(BlockType),target      ::  Blc1

        Blc1%EQN_update_R_IV(:,:,:,Blc1%vr_:Blc1%vp_)=-&
            ModSpherical_A_dot_nabla_B(ni,nj,nk,ng,&
            Blc1%xi_I,Blc1%xj_I,Blc1%dxi,Blc1%dxj,Blc1%dxk,&
            Blc1%primitive(:,:,:,Blc1%vr_:Blc1%vp_),Blc1%primitive(:,:,:,Blc1%vr_:Blc1%vp_))
    end subroutine ModEquation_Corona_Inertial_Force

    ! In the corona p is not defined. p is from T and rho.

    subroutine ModEquation_Corona_Pressure_Gradient(Blc1)
        implicit none
        type(BlockType),target      ::  Blc1

        Blc1%EQN_update_R_IV(:,:,:,Blc1%vr_:Blc1%vp_)=&
        Blc1%EQN_update_R_IV(:,:,:,Blc1%vr_:Blc1%vp_)-&
            ModSpherical_Grad_f(ni,nj,nk,ng,&
            Blc1%xi_I,Blc1%xj_I,Blc1%dxi,Blc1%dxj,Blc1%dxk,Blc1%p_III)
    end subroutine ModEquation_Corona_Pressure_Gradient

    subroutine ModEquation_Corona_Gravity(Blc1)
        implicit none
        type(BlockType),target      ::  Blc1

        Blc1%EQN_update_R_IV(:,:,:,Blc1%vr_)=&
        Blc1%EQN_update_R_IV(:,:,:,Blc1%vr_)-&
            Blc1%g_III(1:ni,1:nj,1:nk)
    end subroutine ModEquation_Corona_Gravity


    subroutine ModEquation_Corona_Lorentz_Force(Blc1)
        implicit none
        type(BlockType),target      ::  Blc1
        integer                     ::  ivar
        real(8)                     ::  tmp(1:ni,1:nj,1:nk,Blc1%vr_:Blc1%vp_)
        real(8)                     ::  curl_B(1:ni,1:nj,1:nk,Blc1%br_:Blc1%bp_)

        curl_B=ModSpherical_curl(                                   &
            ni,nj,nk,ng,Blc1%xi_I,Blc1%xj_I,                    &
            Blc1%dxi,Blc1%dxj,Blc1%dxk,                       &
            Blc1%primitive(:,:,:,Blc1%br_:Blc1%bp_))

        tmp=ModSpherical_cross(ni,nj,nk,0,curl_B,                   &
                Blc1%primitive(1:ni,1:nj,1:nk,Blc1%br_:Blc1%bp_))
        
        do ivar=Blc1%vr_,Blc1%vp_
            Blc1%EQN_update_R_IV(:,:,:,ivar)=Blc1%EQN_update_R_IV(:,:,:,ivar)+&
                tmp(:,:,:,ivar)/Blc1%primitive(1:ni,1:nj,1:nk,Blc1%rho_)
        end do
    end subroutine ModEquation_Corona_Lorentz_Force

    subroutine ModEquation_Corona_Induction_Equation(Blc1)
        implicit none
        type(BlockType),target      ::  Blc1

        Blc1%EQN_update_R_IV(:,:,:,Blc1%br_:Blc1%bp_)=&
            ModSpherical_curl(ni,nj,nk,ng,&
            Blc1%xi_I,Blc1%xj_I,Blc1%dxi,Blc1%dxj,Blc1%dxk,&
            ModSpherical_cross(ni,nj,nk,ng,&
            Blc1%primitive(:,:,:,Blc1%vr_:Blc1%vp_),Blc1%primitive(:,:,:,Blc1%br_:Blc1%bp_)))
    end subroutine ModEquation_Corona_Induction_Equation

    subroutine ModEquation_Corona_Temperature_Advection(Blc1)
        implicit none
        type(BlockType),target      ::  Blc1
        integer                     ::  ivar
        real(8)                     ::  tmp(-ng+1:ng+ni,-ng+1:ng+nj,-ng+1:ng+nk,Blc1%vr_:Blc1%vp_)

        do ivar=Blc1%vr_,Blc1%vp_
            tmp(:,:,:,ivar)=Blc1%primitive(:,:,:,Blc1%Te_)*Blc1%primitive(:,:,:,ivar)
        end do
        Blc1%EQN_update_R_IV(:,:,:,Blc1%Te_)=-ModSpherical_div(ni,nj,nk,ng,&
            Blc1%xi_I,Blc1%xj_I,Blc1%dxi,Blc1%dxj,Blc1%dxk,tmp)
    end subroutine ModEquation_Corona_Temperature_Advection

    subroutine ModEquation_Corona_Temperature_adiabatic(Blc1)
        implicit none
        type(BlockType),target      ::  Blc1

        Blc1%EQN_update_R_IV(:,:,:,Blc1%Te_)=&
        Blc1%EQN_update_R_IV(:,:,:,Blc1%Te_)-&
            (gamma_ideal_gas-2.0)*Blc1%primitive(1:ni,1:nj,1:nk,Blc1%Te_)*&
            ModSpherical_div(ni,nj,nk,ng,&
            Blc1%xi_I,Blc1%xj_I,Blc1%dxi,Blc1%dxj,Blc1%dxk,&
            Blc1%primitive(:,:,:,Blc1%vr_:Blc1%vp_))
    end subroutine ModEquation_Corona_Temperature_adiabatic

    subroutine ModEquation_Corona_MHD(Blc1,if_rk)
        implicit none
        type(BlockType),target      ::  Blc1
        logical,intent(in)          ::  if_rk
        integer                     ::  ivar
        real(8)                     ::  DivB(1:ni,1:nj,1:nk)

        ! If_rk
        if (if_rk) then
            Blc1%primitive=>Blc1%primitive_rk_IV
        else
            Blc1%primitive=>Blc1%primitive_IV
        end if

        Blc1%EQN_update_R_IV=0.

        call ModBoundary_Dynamo_MHD_primitives(Blc1,if_rk)
        call ModEquation_Corona_Mass_Conservation(Blc1)
        call ModEquation_Corona_Inertial_Force(Blc1)
        call ModEquation_Corona_Pressure_Gradient(Blc1)
        call ModEquation_Corona_Gravity(Blc1)
        call ModEquation_Corona_Lorentz_Force(Blc1)
        call ModEquation_Corona_Induction_Equation(Blc1)
        call ModEquation_Corona_Temperature_Advection(Blc1)
        call ModEquation_Corona_Temperature_adiabatic(Blc1)

        ! Div B correction
        DivB=ModSpherical_div(ni,nj,nk,ng,Blc1%xi_I,Blc1%xj_I,&
            Blc1%dxi,Blc1%dxj,Blc1%dxk,Blc1%primitive(:,:,:,Blc1%br_:Blc1%bp_))
        
        do ivar=Blc1%vr_,Blc1%vp_
            Blc1%EQN_update_R_IV(:,:,:,Blc1%br_+ivar-Blc1%vr_)=Blc1%EQN_update_R_IV(:,:,:,Blc1%br_+ivar-Blc1%vr_)-&
                DivB*Blc1%primitive(1:ni,1:nj,1:nk,ivar)
        end do
    end subroutine ModEquation_Corona_MHD
end module 