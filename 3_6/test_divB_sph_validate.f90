! Validation: thin equatorial shell, all-periodic BCs, 2nd-order spherical div & grad.
! No ModSpherical — operators written directly like the Cartesian test.
! Should match Cartesian isotropic 1:1:1 results closely.
program test_divB_sph_validate
    implicit none

    ! ---- Grid --------------------------------------------------------
    integer, parameter :: nr = 32, nt = 32, np = 64
    real(8), parameter :: pi    = acos(-1.0d0)
    ! Thin equatorial shell, isotropic: dr ≈ r_mid*dt ≈ r_mid*dp
    real(8), parameter :: r_min = 0.9d0, r_max = 1.0d0
    real(8), parameter :: r_mid = 0.5d0*(r_min+r_max)
    real(8), parameter :: dr    = (r_max-r_min)/nr
    real(8), parameter :: t_rng = real(nt,8)*dr/r_mid   ! ≈ 0.1053 rad
    real(8), parameter :: t_min = pi/2.0d0 - t_rng/2.0d0
    real(8), parameter :: t_max = pi/2.0d0 + t_rng/2.0d0
    real(8), parameter :: dt    = t_rng/nt
    real(8), parameter :: p_rng = real(np,8)*dr/r_mid   ! ≈ 0.2105 rad
    real(8), parameter :: dp    = p_rng/np

    ! ---- Method parameters -------------------------------------------
    real(8), parameter :: noise_fraction = 0.01d0
    integer, parameter :: n_iter  = 30
    integer, parameter :: n_print = 5
    real(8), parameter :: mu    = 0.4d0
    real(8), parameter :: c_h   = 1.0d0
    real(8), parameter :: CFL   = 0.8d0
    real(8), parameter :: alpha = 1.0d0

    ! ---- Coordinate arrays (interior cells only) ---------------------
    real(8) :: r_arr(nr), t_arr(nt)

    ! ---- Field arrays ------------------------------------------------
    real(8) :: B0(nr, nt, np, 3)
    real(8) :: nB(nr, nt, np, 3)

    integer :: ir, it, ip, seed_size
    integer, allocatable :: seed(:)
    real(8) :: E_B0, E_dB, scale, dh_min
    real(8) :: r_n, t_n, p_n

    ! ---- Coordinates -------------------------------------------------
    do ir = 1, nr
        r_arr(ir) = r_min + (ir - 0.5d0) * dr
    end do
    do it = 1, nt
        t_arr(it) = t_min + (it - 0.5d0) * dt
    end do

    dh_min = min(dr, r_min*dt, r_min*sin(t_min)*dp)

    ! ---- Analytically & discretely div-free B0 -----------------------
    ! r_n,t_n,p_n ∈ [0,1] normalized over each range.
    !   Br = sin(2π t_n)*sin(2π p_n)/r²  →  r²Br independent of r
    !   Bθ = sin(2π r_n)*sin(2π p_n)/sinθ → sinθ Bθ independent of θ
    !   Bφ = sin(2π t_n)*sin(2π r_n)      → no φ-dep
    do ip = 1, np
        p_n = (ip - 0.5d0) / np
        do it = 1, nt
            t_n = (t_arr(it) - t_min) / t_rng
            do ir = 1, nr
                r_n = (r_arr(ir) - r_min) / (r_max - r_min)
                B0(ir,it,ip,1) = sin(2*pi*t_n) * sin(2*pi*p_n) / r_arr(ir)**2
                B0(ir,it,ip,2) = sin(2*pi*r_n) * sin(2*pi*p_n) / sin(t_arr(it))
                B0(ir,it,ip,3) = sin(2*pi*t_n) * sin(2*pi*r_n)
            end do
        end do
    end do

    ! ---- Random perturbation -----------------------------------------
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 12345
    call random_seed(put=seed)
    call random_number(nB)
    nB    = nB - 0.5d0
    E_B0  = mag_energy_sph(B0)
    E_dB  = mag_energy_sph(nB)
    scale = noise_fraction * sqrt(E_B0 / E_dB)
    nB    = nB * scale
    E_dB  = mag_energy_sph(nB)

    print '(A)', '# VALIDATION: near-Cartesian spherical, 2nd-order div & grad, all-periodic'
    print '(A,3F10.6)', '#   dr / r_mid*dt / r_mid*dp =', dr, r_mid*dt, r_mid*dp
    print '(A,F10.6)', '#   dh_min =', dh_min
    print '(A,ES10.3,A,ES10.3)', '#   E_B0 =', E_B0, '   E_dB =', E_dB
    print '(A)', '#'
    print '(A)', '# Columns: iter | divB_rms/divB0 | (E-E_B0)/E_dB | max(divB)/(Brms/dh)'

    print '(/,A)', '#  -- Rempel per-direction'
    call run_rempel(B0, nB, .true.)

    print '(/,A)', '#  -- GLM'
    call run_glm(B0, nB)

contains

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
                    itp = modulo(it,   nt) + 1;  itm = modulo(it-2, nt) + 1
                    do ir = 1, nr
                        irp = modulo(ir,   nr) + 1;  irm = modulo(ir-2, nr) + 1

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
            ! 2nd-order spherical gradient of psi (all-periodic)
            do ip = 1, np
                ipp = modulo(ip,   np) + 1;  ipm = modulo(ip-2, np) + 1
                do it = 1, nt
                    itp = modulo(it,   nt) + 1;  itm = modulo(it-2, nt) + 1
                    do ir = 1, nr
                        irp = modulo(ir,   nr) + 1;  irm = modulo(ir-2, nr) + 1
                        gpr(ir,it,ip) = (psi(irp,it,ip) - psi(irm,it,ip)) / (2*dr)
                        gpt(ir,it,ip) = (psi(ir,itp,ip) - psi(ir,itm,ip)) / (2*dt*r_arr(ir))
                        gpp(ir,it,ip) = (psi(ir,it,ipp) - psi(ir,it,ipm)) / (2*dp*r_arr(ir)*sin(t_arr(it)))
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


    ! 2nd-order spherical divergence (all-periodic)
    ! ∇·B = (1/r²)∂(r²Br)/∂r + (1/r sinθ)[∂(sinθ Bθ)/∂θ + ∂Bφ/∂φ]
    subroutine compute_divB_sph(B, divB)
        real(8), intent(in)  :: B(nr,nt,np,3)
        real(8), intent(out) :: divB(nr,nt,np)
        integer :: ir, it, ip, irp, irm, itp, itm, ipp, ipm
        real(8) :: r2p, r2m

        do ip = 1, np
            ipp = modulo(ip,   np) + 1;  ipm = modulo(ip-2, np) + 1
            do it = 1, nt
                itp = modulo(it,   nt) + 1;  itm = modulo(it-2, nt) + 1
                do ir = 1, nr
                    irp = modulo(ir,   nr) + 1;  irm = modulo(ir-2, nr) + 1
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

end program test_divB_sph_validate
