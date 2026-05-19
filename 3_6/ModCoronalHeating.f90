module ModCoronalHeating

    use ModBlock,       only:   BlockType
    use ModMHD,         only:   ModMHD_AlfvenVelocityVector_3D
    use ModParameters,  only:   ni,nj,nk,ng,LperpSqrtB
    use ModSpherical,   only:   ModSpherical_div, ModSpherical_curl, ModSpherical_dot, &
                                ModSpherical_abs, ModSpherical_Grad_f

    implicit none

    contains

    subroutine ModCoronalHeating_TurbulentCacade(Block1,if_rk)
        Implicit none
        type(BlockType),target      ::  Block1                  ! the block
        logical,intent(in)          ::  if_rk
        real(8)                     ::  AlfvenVelocityVector_IV(-ng+1:ng+ni,-ng+1:ng+nj,-ng+1:ng+nk,3)
        real(8)                     ::  tmp(-ng+1:ng+ni,-ng+1:ng+nj,-ng+1:ng+nk,Block1%vr_:Block1%vp_)
        real(8)                     ::  pm_sign
        integer                     ::  ivar,pm,w_here_

        if (if_rk) then
            Block1%primitive=>Block1%primitive_rk_IV
        else
            Block1%primitive=>Block1%primitive_IV
        end if

        AlfvenVelocityVector_IV=&
            ModMHD_AlfvenVelocityVector_3D(ni,nj,nk,ng,&
            Block1%primitive(:,:,:,Block1%rho_),&
            Block1%primitive(:,:,:,Block1%br_:Block1%bp_))
        
        do pm=1,2
            ! Get the sign of the pm and w_here_
            if (pm==1) then
                pm_sign=1.0d0
                w_here_=Block1%w_plus_
            else
                pm_sign=-1.0d0
                w_here_=Block1%w_minus_
            end if

            ! First term: - div (w_plus/minus * (U +- V_A))
            ! U +- V_A
            tmp(:,:,:,Block1%vr_:Block1%vp_)=&
                Block1%primitive(:,:,:,Block1%vr_:Block1%vp_)+&
                pm_sign*AlfvenVelocityVector_IV(:,:,:,1:3)

            ! Then get w_plus (or w_minus) * (U +- V_A)
            do ivar=Block1%vr_,Block1%vp_
                tmp(:,:,:,ivar)=tmp(:,:,:,ivar)*Block1%primitive(:,:,:,w_here_)
            end do


            ! Then get the divergence of the result
            Block1%EQN_update_R_IV(:,:,:,w_here_)=-ModSpherical_div(ni,nj,nk,ng,&
                Block1%xi_I,Block1%xj_I,Block1%dxi,Block1%dxj,Block1%dxk,tmp)

            !print *,111
            !print *,tmp(:,1,1,Block1%vr_)
            !print *,tmp(:,1,0,Block1%vr_)
            !print *,222


            ! Second term: - 0.5 * div(u) * w_plus/minus
            Block1%EQN_update_R_IV(:,:,:,w_here_)=&
                Block1%EQN_update_R_IV(:,:,:,w_here_)-&
                0.5d0*ModSpherical_div(ni,nj,nk,ng,&
                Block1%xi_I,Block1%xj_I,Block1%dxi,Block1%dxj,Block1%dxk,&
                Block1%primitive(:,:,:,Block1%vr_:Block1%vp_))*&
                Block1%primitive(:,:,:,w_here_)
        end do

        !call ModCoronalHeating_TurbulentReflectionDissipation(Block1,AlfvenVelocityVector_IV)
    end subroutine ModCoronalHeating_TurbulentCacade

    ! Source term for w_pm is \mp R sqrt(w_plus * w_minus) - Gamma_pm * w_pm

    subroutine ModCoronalHeating_TurbulentReflectionDissipation(Block1,AlfvenVelocityVector_IV)
        implicit none
        type(BlockType),target      ::  Block1                  ! the block
        real(8),intent(in)          ::  AlfvenVelocityVector_IV(-ng+1:ng+ni,-ng+1:ng+nj,-ng+1:ng+nk,3)
        real(8)                     ::  b_tot_LLL(1:ni,1:nj,1:nk)
        real(8)                     ::  R_imb_LLL(1:ni,1:nj,1:nk),R_LLL(1:ni,1:nj,1:nk)
        real(8)                     ::  GAMMA_pm_LLL(1:ni,1:nj,1:nk,2)

        b_tot_LLL=max(ModSpherical_abs(ni,nj,nk,0,&
            Block1%primitive(1:ni,1:nj,1:nk,Block1%br_:Block1%bp_)),1.0e-30)

        R_imb_LLL=ModCoronalHeating_TurbulentReflectionRate(Block1,AlfvenVelocityVector_IV,b_tot_LLL)
        GAMMA_pm_LLL=ModCoronalHeating_TurbulentDissipationRate(Block1,b_tot_LLL)

        ! R = min(R_imb,max(Gamma_plus,Gamma_minus))*(max(1-2sqrt(w_plus/w_minus),0)-max(1-2sqrt(w_minus/w_plus),0))

        R_LLL=min(R_imb_LLL,max(GAMMA_pm_LLL(:,:,:,Block1%w_plus_),GAMMA_pm_LLL(:,:,:,Block1%w_minus_)))
        R_LLL=R_LLL*(&
            max(1.0d0-2.0d0*sqrt(Block1%primitive(1:ni,1:nj,1:nk,Block1%w_plus_)/Block1%primitive(1:ni,1:nj,1:nk,Block1%w_minus_)),0.0d0)-&
            max(1.0d0-2.0d0*sqrt(Block1%primitive(1:ni,1:nj,1:nk,Block1%w_minus_)/Block1%primitive(1:ni,1:nj,1:nk,Block1%w_plus_)),0.0d0) &
            )

        Block1%EQN_update_R_IV(:,:,:,Block1%w_plus_)=&
            Block1%EQN_update_R_IV(:,:,:,Block1%w_plus_)-&
            R_LLL*sqrt(Block1%primitive(1:ni,1:nj,1:nk,Block1%w_plus_)*Block1%primitive(1:ni,1:nj,1:nk,Block1%w_minus_))-&
            GAMMA_pm_LLL(:,:,:,Block1%w_plus_)*Block1%primitive(1:ni,1:nj,1:nk,Block1%w_plus_)

        Block1%EQN_update_R_IV(:,:,:,Block1%w_minus_)=&
            Block1%EQN_update_R_IV(:,:,:,Block1%w_minus_)+&
            R_LLL*sqrt(Block1%primitive(1:ni,1:nj,1:nk,Block1%w_plus_)*Block1%primitive(1:ni,1:nj,1:nk,Block1%w_minus_))-&
            GAMMA_pm_LLL(:,:,:,Block1%w_minus_)*Block1%primitive(1:ni,1:nj,1:nk,Block1%w_minus_)

    end subroutine ModCoronalHeating_TurbulentReflectionDissipation

    ! R_imb = sqrt( (b dot (nabla cross u))^2 + (Va dot (nabla log |Va|))^2)
    ! See Van der Holst et al. 2014 Eq. (38)

    function ModCoronalHeating_TurbulentReflectionRate(Block1,AlfvenVelocityVector_IV,b_tot_LLL) result(R_imb_LLL)
        Implicit none
        type(BlockType),target      ::  Block1                  ! the block
        real(8),intent(in)          ::  AlfvenVelocityVector_IV(-ng+1:ng+ni,-ng+1:ng+nj,-ng+1:ng+nk,3)
        real(8)                     ::  vorticity_LLL(1:ni,1:nj,1:nk,3)
        real(8),intent(in)          ::  b_tot_LLL(1:ni,1:nj,1:nk)
        real(8)                     ::  b_hat_IV(1:ni,1:nj,1:nk,3)
        real(8)                     ::  Va_III(-ng+1:ng+ni,-ng+1:ng+nj,-ng+1:ng+nk)
        real(8)                     ::  log_Va_gradient_LLL(1:ni,1:nj,1:nk,3)
        real(8)                     ::  R_imb_LLL(1:ni,1:nj,1:nk)

        vorticity_LLL=ModSpherical_curl(ni,nj,nk,ng,&
            Block1%xi_I,Block1%xj_I,Block1%dxi,Block1%dxj,Block1%dxk,&
            Block1%primitive(:,:,:,Block1%vr_:Block1%vp_))

        b_hat_IV(:,:,:,1)=Block1%primitive(1:ni,1:nj,1:nk,Block1%br_)/b_tot_LLL
        b_hat_IV(:,:,:,2)=Block1%primitive(1:ni,1:nj,1:nk,Block1%bt_)/b_tot_LLL
        b_hat_IV(:,:,:,3)=Block1%primitive(1:ni,1:nj,1:nk,Block1%bp_)/b_tot_LLL

        Va_III=ModSpherical_abs(ni,nj,nk,ng,AlfvenVelocityVector_IV)

        log_Va_gradient_LLL=ModSpherical_Grad_f(ni,nj,nk,ng,&
            Block1%xi_I,Block1%xj_I,Block1%dxi,Block1%dxj,Block1%dxk,&
            log(Va_III))

        R_imb_LLL=sqrt((ModSpherical_dot(ni,nj,nk,0,vorticity_LLL,b_hat_IV))**2+&
            (ModSpherical_dot(ni,nj,nk,0,log_Va_gradient_LLL,Va_III(1:ni,1:nj,1:nk)))**2)
    end function ModCoronalHeating_TurbulentReflectionRate

    ! Gamma_pm = 2 / L_perp * sqrt(w_mp / rho)
    ! L_perp = (L_perp sqrt(B))/sqrt(B)

    function ModCoronalHeating_TurbulentDissipationRate(Block1,btot_LLL) result(GAMMA_pm_LLL)
        implicit none
        type(BlockType),target      ::  Block1                  ! the block
        real(8),intent(in)          ::  btot_LLL(1:ni,1:nj,1:nk)
        real(8)                     ::  GAMMA_pm_LLL(1:ni,1:nj,1:nk,Block1%w_plus_:Block1%w_minus_)
        real(8)                     ::  L_perp_LLL(1:ni,1:nj,1:nk)

        L_perp_LLL=LperpSqrtB*sqrt(btot_LLL)

        GAMMA_pm_LLL(:,:,:,Block1%w_plus_)=2.0d0/L_perp_LLL*sqrt(Block1%primitive(:,:,:,Block1%w_minus_)/&
            Block1%primitive(:,:,:,Block1%rho_))

        GAMMA_pm_LLL(:,:,:,Block1%w_minus_)=2.0d0/L_perp_LLL*sqrt(Block1%primitive(:,:,:,Block1%w_plus_)/&
            Block1%primitive(:,:,:,Block1%rho_))
    end function ModCoronalHeating_TurbulentDissipationRate
end module









