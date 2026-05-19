program test_divB_sph
    implicit none

    ! ---- Grid --------------------------------------------------------
    integer, parameter :: nr = 32, nt = 32, np = 64
    real(8), parameter :: pi    = acos(-1.0d0)
    real(8), parameter :: r_min = 0.7d0,    r_max = 0.96d0
    real(8), parameter :: t_min = pi/4.0d0, t_max = 3.0d0*pi/4.0d0
    real(8), parameter :: p_min = 0.0d0,    p_max = 2.0d0*pi
    real(8), parameter :: dr    = (r_max-r_min)/nr
    real(8), parameter :: dt    = (t_max-t_min)/nt
    real(8), parameter :: dp    = (p_max-p_min)/np

    ! ---- Method parameters -------------------------------------------
    real(8), parameter :: noise_fraction = 0.01d0
    integer, parameter :: n_iter  = 30
    integer, parameter :: n_print = 5
    real(8), parameter :: mu    = 0.4d0   ! Rempel stability: mu < 2/3
    real(8), parameter :: c_h   = 1.0d0
    real(8), parameter :: CFL   = 0.8d0
    real(8), parameter :: alpha = 1.0d0   ! tau = alpha * dh/c_h

    ! ---- Coordinate arrays (interior cells) --------------------------
    real(8) :: r_arr(nr), t_arr(nt)

    ! ---- Field arrays ------------------------------------------------
    real(8) :: B0(nr, nt, np, 3)
    real(8) :: nB(nr, nt, np, 3)
    real(8) :: v(nr, nt, np, 3)
    real(8) :: EMF(nr, nt, np, 3)

    integer :: ir, it, ip, irp, irm, itp, itm, ipp, ipm, seed_size
    integer, allocatable :: seed(:)
    real(8) :: E_B0, E_dB, scale, dh_min

    ! ---- Coordinates -------------------------------------------------
    do ir = 1, nr
        r_arr(ir) = r_min + (ir - 0.5d0) * dr
    end do
    do it = 1, nt
        t_arr(it) = t_min + (it - 0.5d0) * dt
    end do

    ! minimum physical cell spacing: dr vs r*dθ vs r*sinθ*dφ
    dh_min = min(dr, r_min*dt, r_min*sin(t_min)*dp)

    ! ---- Analytically and discretely div-free background B0 ----------
    ! ∇·B = (1/r²)∂(r²Br)/∂r + (1/r sinθ)[∂(sinθ Bθ)/∂θ + ∂Bφ/∂φ]
    ! Each term independently zero (for any stencil order):
    !   Br = sin(θ)sin(φ)/r²        → r²Br = f(θ,φ)    → ∂(r²Br)/∂r  = 0
    !   Bθ = sin(2π r_n)cos(φ)/sinθ → sinθ Bθ = g(r,φ) → ∂(sinθ Bθ)/∂θ = 0
    !   Bφ = cos(θ)sin(π r_n)       → no φ-dep          → ∂Bφ/∂φ      = 0
    do ip = 1, np
        do it = 1, nt
            do ir = 1, nr
                B0(ir,it,ip,1) = sin(t_arr(it)) * sin(p_min+(ip-0.5d0)*dp) / r_arr(ir)**2
                B0(ir,it,ip,2) = sin(2*pi*(r_arr(ir)-r_min)/(r_max-r_min)) * &
                                  cos(p_min+(ip-0.5d0)*dp) / sin(t_arr(it))
                B0(ir,it,ip,3) = cos(t_arr(it)) * sin(pi*(r_arr(ir)-r_min)/(r_max-r_min))
            end do
        end do
    end do

    ! ---- Induction-equation perturbation: nB = ∇×(v×B0) ---------------
    ! Random v → EMF = v×B0 → nB = ∇×EMF.
    ! Analytically div-free but non-zero discrete divergence, mimicking
    ! truncation-error divB in cell-centred MHD codes.
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 12345
    call random_seed(put=seed)
    call random_number(v)
    v = v - 0.5d0

    do ip = 1, np
        do it = 1, nt
            do ir = 1, nr
                EMF(ir,it,ip,1) = v(ir,it,ip,2)*B0(ir,it,ip,3) - v(ir,it,ip,3)*B0(ir,it,ip,2)
                EMF(ir,it,ip,2) = v(ir,it,ip,3)*B0(ir,it,ip,1) - v(ir,it,ip,1)*B0(ir,it,ip,3)
                EMF(ir,it,ip,3) = v(ir,it,ip,1)*B0(ir,it,ip,2) - v(ir,it,ip,2)*B0(ir,it,ip,1)
            end do
        end do
    end do

    ! 2nd-order spherical curl, same BCs as div (periodic φ, clamped r/θ)
    ! (∇×F)_r = [∂(sinθ Fφ)/∂θ - ∂Fθ/∂φ] / (r sinθ)
    ! (∇×F)_θ = [∂Fr/∂φ/sinθ  - ∂(r Fφ)/∂r] / r
    ! (∇×F)_φ = [∂(r Fθ)/∂r   - ∂Fr/∂θ]     / r
    do ip = 1, np
        ipp = modulo(ip,   np) + 1;  ipm = modulo(ip-2, np) + 1
        do it = 1, nt
            itp = min(it+1, nt);  itm = max(it-1, 1)
            do ir = 1, nr
                irp = min(ir+1, nr);  irm = max(ir-1, 1)
                nB(ir,it,ip,1) = &
                    (sin(t_arr(itp))*EMF(ir,itp,ip,3) - sin(t_arr(itm))*EMF(ir,itm,ip,3)) / &
                        (2*dt * r_arr(ir) * sin(t_arr(it))) - &
                    (EMF(ir,it,ipp,2) - EMF(ir,it,ipm,2)) / &
                        (2*dp * r_arr(ir) * sin(t_arr(it)))
                nB(ir,it,ip,2) = &
                    (EMF(ir,it,ipp,1) - EMF(ir,it,ipm,1)) / &
                        (2*dp * r_arr(ir) * sin(t_arr(it))) - &
                    (r_arr(irp)*EMF(irp,it,ip,3) - r_arr(irm)*EMF(irm,it,ip,3)) / &
                        (2*dr * r_arr(ir))
                nB(ir,it,ip,3) = &
                    (r_arr(irp)*EMF(irp,it,ip,2) - r_arr(irm)*EMF(irm,it,ip,2)) / &
                        (2*dr * r_arr(ir)) - &
                    (EMF(ir,itp,ip,1) - EMF(ir,itm,ip,1)) / &
                        (2*dt * r_arr(ir))
            end do
        end do
    end do

    E_B0  = mag_energy_sph(B0)
    E_dB  = mag_energy_sph(nB)
    scale = noise_fraction * sqrt(E_B0 / E_dB)
    nB    = nB * scale
    E_dB  = mag_energy_sph(nB)

    print '(A)', '# Div-B cleaning: spherical wedge, 2nd-order div & grad'
    print '(A,F5.3,A)', '# perturbation = curl(v x B0), v random, scaled to', noise_fraction*100, '% energy'
    print '(A)', '#   r:[0.7,0.96]  th:[pi/4,3pi/4]  ph:[0,2pi]'
    print '(A,3F10.6)', '#   dr / r_min*dt / r_min*sin(t_min)*dp =', &
        dr, r_min*dt, r_min*sin(t_min)*dp
    print '(A,F10.6)', '#   dh_min =', dh_min
    print '(A,ES10.3,A,ES10.3)', '#   E_B0 =', E_B0, '   E_dB =', E_dB
    print '(A)', '#'
    print '(A)', '# Columns: iter | divB_rms/divB0 | (E-E_B0)/E_dB | max(divB)/(Brms/dh)'

    print '(/,A)', '#  -- Rempel per-direction (h = local physical spacing)'
    call run_rempel(B0, nB, .true.)

    print '(/,A)', '#  -- Rempel uniform min-h'
    call run_rempel(B0, nB, .false.)

    print '(/,A,3F7.3)', '#  -- GLM  c_h/CFL/alpha =', c_h, CFL, alpha
    call run_glm(B0, nB)

contains

    ! Rempel (2009) iterative cleaning.
    ! per_dir=T: h = local physical spacing (hr=dr, ht=r*dθ, hp=r*sinθ*dφ)
    ! per_dir=F: h = dh_min everywhere
    subroutine run_rempel(B0, nB, per_dir)
        real(8), intent(in) :: B0(nr,nt,np,3), nB(nr,nt,np,3)
        logical, intent(in) :: per_dir

        real(8) :: B(nr,nt,np,3), divB(nr,nt,np)
        real(8) :: dBr(nr,nt,np), dBt(nr,nt,np), dBp(nr,nt,np)
        real(8) :: E_ref, E_pert0, divB_rms_0, divB_rms, E_now, Brms, max_divB
        real(8) :: hr, ht, hp, cr, ct, cp
        integer :: ir, it, ip, iter, irp, irm, itp, itm, ipp, ipm

        B = B0 + nB
        E_ref   = mag_energy_sph(B0)
        E_pert0 = mag_energy_sph(nB)

        call compute_divB_sph(B, divB)
        divB_rms_0 = rms3d(divB)
        E_now      = mag_energy_sph(B)
        Brms       = sqrt(2.0d0 * E_now)
        max_divB   = maxval(abs(divB))
        print '(I9, 3ES16.5)', 0, 1.0d0, (E_now-E_ref)/E_pert0, max_divB/(Brms/dh_min)

        do iter = 1, n_iter
            do ip = 1, np
                ipp = modulo(ip,   np) + 1;  ipm = modulo(ip-2, np) + 1
                do it = 1, nt
                    itp = min(it+1, nt);  itm = max(it-1, 1)
                    do ir = 1, nr
                        irp = min(ir+1, nr);  irm = max(ir-1, 1)

                        if (per_dir) then
                            hr = dr
                            ht = r_arr(ir) * dt
                            hp = r_arr(ir) * sin(t_arr(it)) * dp
                        else
                            hr = dh_min;  ht = dh_min;  hp = dh_min
                        end if
                        cr = mu * hr**2 / dr                              * 0.5d0
                        ct = mu * ht**2 / (r_arr(ir)*dt)                  * 0.5d0
                        cp = mu * hp**2 / (r_arr(ir)*sin(t_arr(it))*dp)   * 0.5d0

                        dBr(ir,it,ip) = cr * (divB(irp,it,ip) - divB(irm,it,ip))
                        dBt(ir,it,ip) = ct * (divB(ir,itp,ip) - divB(ir,itm,ip))
                        dBp(ir,it,ip) = cp * (divB(ir,it,ipp) - divB(ir,it,ipm))
                    end do
                end do
            end do

            B(:,:,:,1) = B(:,:,:,1) + dBr
            B(:,:,:,2) = B(:,:,:,2) + dBt
            B(:,:,:,3) = B(:,:,:,3) + dBp

            call compute_divB_sph(B, divB)

            if (mod(iter, n_print) == 0) then
                divB_rms = rms3d(divB)
                E_now    = mag_energy_sph(B)
                Brms     = sqrt(2.0d0 * E_now)
                max_divB = maxval(abs(divB))
                print '(I9, 3ES16.5)', iter, &
                    divB_rms/divB_rms_0, (E_now-E_ref)/E_pert0, max_divB/(Brms/dh_min)
            end if
        end do

        divB_rms = rms3d(divB)
        E_now    = mag_energy_sph(B)
        Brms     = sqrt(2.0d0 * E_now)
        max_divB = maxval(abs(divB))
        print '(A,3ES11.4)', '#     final =', &
            divB_rms/divB_rms_0, (E_now-E_ref)/E_pert0, max_divB/(Brms/dh_min)
    end subroutine run_rempel


    ! GLM (Dedner 2002): dB/dt = -∇ψ,  dψ/dt = -c_h² ∇·B - ψ/τ
    ! dt = dh/c_h * CFL,  τ = dh/c_h * α
    ! ψ^{n+1} = exp(-dt/τ)*ψ^n - dt*c_h²*(∇·B)^n
    subroutine run_glm(B0, nB)
        real(8), intent(in) :: B0(nr,nt,np,3), nB(nr,nt,np,3)

        real(8) :: B(nr,nt,np,3), psi(nr,nt,np), divB(nr,nt,np)
        real(8) :: gpr(nr,nt,np), gpt(nr,nt,np), gpp(nr,nt,np)
        real(8) :: E_ref, E_pert0, divB_rms_0, divB_rms, E_now, Brms, max_divB
        real(8) :: dt_glm, tau, decay
        integer :: ir, it, ip, iter, irp, irm, itp, itm, ipp, ipm

        dt_glm = dh_min / c_h * CFL
        tau    = dh_min / c_h * alpha
        decay  = exp(-dt_glm / tau)
        print '(A,3ES11.3)', '#     dh_min/dt/tau =', dh_min, dt_glm, tau

        B   = B0 + nB
        psi = 0.0d0
        E_ref   = mag_energy_sph(B0)
        E_pert0 = mag_energy_sph(nB)

        call compute_divB_sph(B, divB)
        divB_rms_0 = rms3d(divB)
        E_now      = mag_energy_sph(B)
        Brms       = sqrt(2.0d0 * E_now)
        max_divB   = maxval(abs(divB))
        print '(I9, 3ES16.5)', 0, 1.0d0, (E_now-E_ref)/E_pert0, max_divB/(Brms/dh_min)

        do iter = 1, n_iter
            ! 2nd-order spherical gradient of ψ
            ! BCs: periodic φ, zero-gradient r and θ
            do ip = 1, np
                ipp = modulo(ip,   np) + 1;  ipm = modulo(ip-2, np) + 1
                do it = 1, nt
                    itp = min(it+1, nt);  itm = max(it-1, 1)
                    do ir = 1, nr
                        irp = min(ir+1, nr);  irm = max(ir-1, 1)
                        gpr(ir,it,ip) = (psi(irp,it,ip) - psi(irm,it,ip)) / (2*dr)
                        gpt(ir,it,ip) = (psi(ir,itp,ip) - psi(ir,itm,ip)) / (2*dt*r_arr(ir))
                        gpp(ir,it,ip) = (psi(ir,it,ipp) - psi(ir,it,ipm)) / &
                                         (2*dp*r_arr(ir)*sin(t_arr(it)))
                    end do
                end do
            end do

            B(:,:,:,1) = B(:,:,:,1) - dt_glm * gpr
            B(:,:,:,2) = B(:,:,:,2) - dt_glm * gpt
            B(:,:,:,3) = B(:,:,:,3) - dt_glm * gpp

            call compute_divB_sph(B, divB)
            psi = decay * psi - dt_glm * c_h**2 * divB

            if (mod(iter, n_print) == 0) then
                divB_rms = rms3d(divB)
                E_now    = mag_energy_sph(B)
                Brms     = sqrt(2.0d0 * E_now)
                max_divB = maxval(abs(divB))
                print '(I9, 3ES16.5)', iter, &
                    divB_rms/divB_rms_0, (E_now-E_ref)/E_pert0, max_divB/(Brms/dh_min)
            end if
        end do

        divB_rms = rms3d(divB)
        E_now    = mag_energy_sph(B)
        Brms     = sqrt(2.0d0 * E_now)
        max_divB = maxval(abs(divB))
        print '(A,3ES11.4)', '#     final =', &
            divB_rms/divB_rms_0, (E_now-E_ref)/E_pert0, max_divB/(Brms/dh_min)
    end subroutine run_glm


    ! 2nd-order spherical divergence.
    ! ∇·B = (1/r²)∂(r²Br)/∂r + (1/r sinθ)[∂(sinθ Bθ)/∂θ + ∂Bφ/∂φ]
    ! BCs: periodic φ, zero-gradient r and θ (clamp indices).
    subroutine compute_divB_sph(B, divB)
        real(8), intent(in)  :: B(nr,nt,np,3)
        real(8), intent(out) :: divB(nr,nt,np)
        integer :: ir, it, ip, irp, irm, itp, itm, ipp, ipm
        real(8) :: r2p, r2m

        do ip = 1, np
            ipp = modulo(ip,   np) + 1;  ipm = modulo(ip-2, np) + 1
            do it = 1, nt
                itp = min(it+1, nt);  itm = max(it-1, 1)
                do ir = 1, nr
                    irp = min(ir+1, nr);  irm = max(ir-1, 1)
                    r2p = r_arr(irp)**2;  r2m = r_arr(irm)**2
                    divB(ir,it,ip) = &
                        (r2p*B(irp,it,ip,1) - r2m*B(irm,it,ip,1)) / (2*dr * r_arr(ir)**2) + &
                        (sin(t_arr(itp))*B(ir,itp,ip,2) - sin(t_arr(itm))*B(ir,itm,ip,2)) / &
                            (2*dt * r_arr(ir) * sin(t_arr(it))) + &
                        (B(ir,it,ipp,3) - B(ir,it,ipm,3)) / (2*dp * r_arr(ir) * sin(t_arr(it)))
                end do
            end do
        end do
    end subroutine compute_divB_sph


    real(8) function rms3d(f)
        real(8), intent(in) :: f(nr,nt,np)
        rms3d = sqrt(sum(f**2) / real(nr*nt*np, 8))
    end function rms3d

    real(8) function mag_energy_sph(B)
        real(8), intent(in) :: B(nr,nt,np,3)
        real(8) :: wt, vol_tot
        integer :: ir, it, ip
        mag_energy_sph = 0.0d0
        vol_tot        = 0.0d0
        do ip = 1, np
            do it = 1, nt
                do ir = 1, nr
                    wt             = r_arr(ir)**2 * sin(t_arr(it))
                    mag_energy_sph = mag_energy_sph + &
                        0.5d0*wt*(B(ir,it,ip,1)**2 + B(ir,it,ip,2)**2 + B(ir,it,ip,3)**2)
                    vol_tot        = vol_tot + wt
                end do
            end do
        end do
        mag_energy_sph = mag_energy_sph / vol_tot
    end function mag_energy_sph

end program test_divB_sph
