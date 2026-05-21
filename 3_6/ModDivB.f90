module ModDivB

    use ModBlock,      only: BlockType
    use ModSpherical,  only: ModSpherical_div, ModSpherical_Grad_f
    use ModParameters, only: ni, nj, nk, ng

    contains

    ! GLM div-B cleaning (Dedner 2002).
    ! Adds to EQN_update_R_IV:
    !   dB/dt  += -∇ψ
    !   dψ/dt   = -c_h² ∇·B  - ψ/τ,   τ = h_LLL/c_h,  c_h = v_wave (local)
    subroutine ModDivB_GLM(Block1)
        implicit none
        type(BlockType), intent(inout) :: Block1

        real(8) :: divB(1:ni, 1:nj, 1:nk)
        real(8) :: grad_psi(1:ni, 1:nj, 1:nk, 1:3)

        divB = ModSpherical_div(ni, nj, nk, ng, &
                    Block1%xi_I, Block1%xj_I, Block1%dxi, Block1%dxj, Block1%dxk, &
                    Block1%primitive(:,:,:, Block1%br_:Block1%bp_))

        grad_psi = ModSpherical_Grad_f(ni, nj, nk, ng, &
                    Block1%xi_I, Block1%xj_I, Block1%dxi, Block1%dxj, Block1%dxk, &
                    Block1%primitive(:,:,:, Block1%psi_))

        ! dB/dt += -∇ψ
        Block1%EQN_update_R_IV(:,:,:, Block1%br_:Block1%bp_) = &
            Block1%EQN_update_R_IV(:,:,:, Block1%br_:Block1%bp_) - grad_psi

        ! dψ/dt = -c_h² ∇·B - ψ/τ,  τ = h_LLL/c_h,  c_h = v_wave
        Block1%EQN_update_R_IV(:,:,:, Block1%psi_) = &
            -Block1%v_wave_III(1:ni,1:nj,1:nk)**2 * divB &
            - Block1%v_wave_III(1:ni,1:nj,1:nk) / Block1%h_LLL * &
              Block1%primitive(1:ni,1:nj,1:nk, Block1%psi_)

    end subroutine ModDivB_GLM

end module ModDivB
