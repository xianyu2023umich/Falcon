program test_divB_cleaning
    implicit none

    ! ---- Grid -------------------------------------------------------
    integer,  parameter :: nx = 32, ny = 32, nz = 32
    real(8),  parameter :: pi = acos(-1.0d0)
    ! ---- Common parameters ------------------------------------------
    real(8),  parameter :: noise_fraction = 0.01d0   ! |δB|_rms / |B0|_rms
    integer,  parameter :: n_iter         = 30
    integer,  parameter :: n_print        = 5
    ! ---- Rempel parameters ------------------------------------------
    ! Stability: mu < 2/3 (Fourier analysis)
    real(8),  parameter :: mu    = 0.4d0
    ! ---- GLM parameters ---------------------------------------------
    ! dt   = dh/c_h * CFL           (CFL condition on the ψ-wave speed)
    ! tau  = dh/c_h * alpha         (damping timescale)
    ! ψ update: exponential decay + source (avoids instability for large dt/tau)
    !   ψ^{n+1} = exp(-dt/tau) * ψ^n - dt * c_h^2 * divB^n
    real(8),  parameter :: c_h   = 1.0d0   ! wave speed (arbitrary here)
    real(8),  parameter :: CFL   = 0.8d0
    real(8),  parameter :: alpha = 1.0d0   ! tau = alpha * dh/c_h
    ! -----------------------------------------------------------------

    ! Background potential field (analytically div-free):
    !   B0x = sin(2π j/ny) — depends only on y  → ∂B0x/∂x = 0
    !   B0y = sin(2π k/nz) — depends only on z  → ∂B0y/∂y = 0
    !   B0z = sin(2π i/nx) — depends only on x  → ∂B0z/∂z = 0
    ! Discrete ∇·B0 = 0 exactly, for any dx/dy/dz.
    real(8) :: B0x(nx,ny,nz), B0y(nx,ny,nz), B0z(nx,ny,nz)
    real(8) :: nBx(nx,ny,nz), nBy(nx,ny,nz), nBz(nx,ny,nz)

    integer :: i, j, k, seed_size
    integer, allocatable :: seed(:)
    real(8) :: E_B0, E_dB, scale

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 12345
    call random_seed(put=seed)

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                B0x(i,j,k) = sin(2*pi*(j-0.5d0)/ny)
                B0y(i,j,k) = sin(2*pi*(k-0.5d0)/nz)
                B0z(i,j,k) = sin(2*pi*(i-0.5d0)/nx)
            end do
        end do
    end do

    call random_number(nBx);  nBx = nBx - 0.5d0
    call random_number(nBy);  nBy = nBy - 0.5d0
    call random_number(nBz);  nBz = nBz - 0.5d0
    E_B0  = mag_energy(B0x, B0y, B0z)
    E_dB  = mag_energy(nBx, nBy, nBz)
    scale = noise_fraction * sqrt(E_B0 / E_dB)
    nBx = nBx * scale;  nBy = nBy * scale;  nBz = nBz * scale
    E_dB  = mag_energy(nBx, nBy, nBz)

    print '(A)', '# Div-B cleaning comparison: Rempel vs GLM'
    print '(A,F6.3,A)', '# noise_fraction =', noise_fraction*100, '%  (of B0 rms)'
    print '(A,I3,A,I3)', '# n_iter =', n_iter, '   n_print =', n_print
    print '(A,F5.3,A,F5.3,A,F5.3,A,F5.3)', &
        '# Rempel: mu =', mu, '   GLM: c_h =', c_h, '  CFL =', CFL, '  alpha =', alpha
    print '(A,ES10.3,A,ES10.3)', '#   E_B0 =', E_B0, '   E_dB =', E_dB
    print '(A)', '#'
    print '(A)', '# Columns: iter | divB_rms/divB0 | (E-E_B0)/E_dB | max(divB)/(Brms/dh)'

    call run_test(B0x, B0y, B0z, nBx, nBy, nBz, 1.0d0,  1.0d0,  1.0d0,  'isotropic   1:1:1')
    call run_test(B0x, B0y, B0z, nBx, nBy, nBz, 1.0d0,  5.0d0, 25.0d0,  'anisotropic 1:5:25')

contains

    subroutine run_test(B0x, B0y, B0z, nBx, nBy, nBz, dx, dy, dz, label)
        real(8), intent(in)          :: B0x(nx,ny,nz), B0y(nx,ny,nz), B0z(nx,ny,nz)
        real(8), intent(in)          :: nBx(nx,ny,nz), nBy(nx,ny,nz), nBz(nx,ny,nz)
        real(8), intent(in)          :: dx, dy, dz
        character(len=*), intent(in) :: label
        real(8) :: dh

        dh = min(dx, dy, dz)

        print '(/,A,A)', '# === ', label
        print '(A,3ES10.3,A,ES10.3)', '#     dx/dy/dz =', dx, dy, dz, '   dh =', dh

        print '(/,A)', '#  -- Rempel per-direction h'
        call run_rempel(B0x, B0y, B0z, nBx, nBy, nBz, dx, dy, dz, dx, dy, dz)

        if (abs(dh - dx) > 1d-10 .or. abs(dh - dy) > 1d-10 .or. abs(dh - dz) > 1d-10) then
            print '(/,A)', '#  -- Rempel uniform min-h'
            call run_rempel(B0x, B0y, B0z, nBx, nBy, nBz, dx, dy, dz, dh, dh, dh)
        end if

        print '(/,A,F5.3,A,F5.3,A,F5.3)', &
            '#  -- GLM  c_h =', c_h, '  CFL =', CFL, '  alpha =', alpha
        call run_glm(B0x, B0y, B0z, nBx, nBy, nBz, dx, dy, dz)
    end subroutine run_test


    ! Rempel (2009) iterative method.
    ! hx,hy,hz are the h values used in the correction:
    !   per-direction: hx=dx, hy=dy, hz=dz  →  dBi = mu*(hi^2/di)*(divB+ - divB-)/2
    !   uniform min-h: hx=hy=hz=dh
    subroutine run_rempel(B0x, B0y, B0z, nBx, nBy, nBz, dx, dy, dz, hx, hy, hz)
        real(8), intent(in) :: B0x(nx,ny,nz), B0y(nx,ny,nz), B0z(nx,ny,nz)
        real(8), intent(in) :: nBx(nx,ny,nz), nBy(nx,ny,nz), nBz(nx,ny,nz)
        real(8), intent(in) :: dx, dy, dz, hx, hy, hz

        real(8) :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz)
        real(8) :: divB(nx,ny,nz), dBx(nx,ny,nz), dBy(nx,ny,nz), dBz(nx,ny,nz)
        real(8) :: divB_rms_0, divB_rms, E_now, E_ref, E_pert0
        real(8) :: Brms, max_divB, dh, cx, cy, cz
        integer :: i, j, k, iter, ip, im, jp, jm, kp, km

        dh = min(dx, dy, dz)
        cx = mu * hx**2 / dx * 0.5d0
        cy = mu * hy**2 / dy * 0.5d0
        cz = mu * hz**2 / dz * 0.5d0

        Bx = B0x + nBx;  By = B0y + nBy;  Bz = B0z + nBz
        E_ref   = mag_energy(B0x, B0y, B0z)
        E_pert0 = mag_energy(nBx, nBy, nBz)

        call compute_divB(Bx, By, Bz, divB, dx, dy, dz)
        divB_rms_0 = rms3d(divB)
        E_now      = mag_energy(Bx, By, Bz)
        Brms       = sqrt(2.0d0 * E_now)
        max_divB   = maxval(abs(divB))
        print '(I9, 3ES16.5)', 0, 1.0d0, (E_now-E_ref)/E_pert0, max_divB/(Brms/dh)

        do iter = 1, n_iter
            do k = 1, nz
                kp = modulo(k,   nz) + 1;  km = modulo(k-2, nz) + 1
                do j = 1, ny
                    jp = modulo(j,   ny) + 1;  jm = modulo(j-2, ny) + 1
                    do i = 1, nx
                        ip = modulo(i,   nx) + 1;  im = modulo(i-2, nx) + 1
                        dBx(i,j,k) = cx * (divB(ip,j,k) - divB(im,j,k))
                        dBy(i,j,k) = cy * (divB(i,jp,k) - divB(i,jm,k))
                        dBz(i,j,k) = cz * (divB(i,j,kp) - divB(i,j,km))
                    end do
                end do
            end do
            Bx = Bx + dBx;  By = By + dBy;  Bz = Bz + dBz
            call compute_divB(Bx, By, Bz, divB, dx, dy, dz)
            if (mod(iter, n_print) == 0) then
                divB_rms = rms3d(divB)
                E_now    = mag_energy(Bx, By, Bz)
                Brms     = sqrt(2.0d0 * E_now)
                max_divB = maxval(abs(divB))
                print '(I9, 3ES16.5)', iter, &
                    divB_rms/divB_rms_0, (E_now-E_ref)/E_pert0, max_divB/(Brms/dh)
            end if
        end do

        divB_rms = rms3d(divB)
        E_now    = mag_energy(Bx, By, Bz)
        Brms     = sqrt(2.0d0 * E_now)
        max_divB = maxval(abs(divB))
        print '(A,3ES11.4)', '#     final =', &
            divB_rms/divB_rms_0, (E_now-E_ref)/E_pert0, max_divB/(Brms/dh)
    end subroutine run_rempel


    ! GLM (Dedner 2002) mixed hyperbolic-parabolic cleaning.
    ! Equations:  dB/dt = -grad(psi)
    !             dpsi/dt = -c_h^2 * div(B) - psi/tau
    ! dt  = dh/c_h * CFL
    ! tau = dh/c_h * alpha
    ! psi update uses exponential decay to handle large dt/tau stably:
    !   psi^{n+1} = exp(-dt/tau) * psi^n - dt * c_h^2 * divB^n
    subroutine run_glm(B0x, B0y, B0z, nBx, nBy, nBz, dx, dy, dz)
        real(8), intent(in) :: B0x(nx,ny,nz), B0y(nx,ny,nz), B0z(nx,ny,nz)
        real(8), intent(in) :: nBx(nx,ny,nz), nBy(nx,ny,nz), nBz(nx,ny,nz)
        real(8), intent(in) :: dx, dy, dz

        real(8) :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz)
        real(8) :: psi(nx,ny,nz), divB(nx,ny,nz)
        real(8) :: gpx(nx,ny,nz), gpy(nx,ny,nz), gpz(nx,ny,nz)
        real(8) :: divB_rms_0, divB_rms, E_now, E_ref, E_pert0
        real(8) :: Brms, max_divB, dh, dt, tau, decay
        integer :: i, j, k, iter, ip, im, jp, jm, kp, km

        dh    = min(dx, dy, dz)
        dt    = dh / c_h * CFL
        tau   = dh / c_h * alpha
        decay = exp(-dt / tau)

        print '(A,3ES11.3)', '#     dh/dt/tau =', dh, dt, tau

        Bx = B0x + nBx;  By = B0y + nBy;  Bz = B0z + nBz
        psi = 0.0d0
        E_ref   = mag_energy(B0x, B0y, B0z)
        E_pert0 = mag_energy(nBx, nBy, nBz)

        call compute_divB(Bx, By, Bz, divB, dx, dy, dz)
        divB_rms_0 = rms3d(divB)
        E_now      = mag_energy(Bx, By, Bz)
        Brms       = sqrt(2.0d0 * E_now)
        max_divB   = maxval(abs(divB))
        print '(I9, 3ES16.5)', 0, 1.0d0, (E_now-E_ref)/E_pert0, max_divB/(Brms/dh)

        do iter = 1, n_iter
            ! Gradient of psi  (central differences)
            do k = 1, nz
                kp = modulo(k,   nz) + 1;  km = modulo(k-2, nz) + 1
                do j = 1, ny
                    jp = modulo(j,   ny) + 1;  jm = modulo(j-2, ny) + 1
                    do i = 1, nx
                        ip = modulo(i,   nx) + 1;  im = modulo(i-2, nx) + 1
                        gpx(i,j,k) = (psi(ip,j,k) - psi(im,j,k)) / (2*dx)
                        gpy(i,j,k) = (psi(i,jp,k) - psi(i,jm,k)) / (2*dy)
                        gpz(i,j,k) = (psi(i,j,kp) - psi(i,j,km)) / (2*dz)
                    end do
                end do
            end do

            ! Update B: dB/dt = -grad(psi)
            Bx = Bx - dt * gpx
            By = By - dt * gpy
            Bz = Bz - dt * gpz

            ! Recompute divB on updated B
            call compute_divB(Bx, By, Bz, divB, dx, dy, dz)

            ! Update psi: exponential decay + source
            ! dpsi/dt = -c_h^2 * divB - psi/tau
            psi = decay * psi - dt * c_h**2 * divB

            if (mod(iter, n_print) == 0) then
                divB_rms = rms3d(divB)
                E_now    = mag_energy(Bx, By, Bz)
                Brms     = sqrt(2.0d0 * E_now)
                max_divB = maxval(abs(divB))
                print '(I9, 3ES16.5)', iter, &
                    divB_rms/divB_rms_0, (E_now-E_ref)/E_pert0, max_divB/(Brms/dh)
            end if
        end do

        divB_rms = rms3d(divB)
        E_now    = mag_energy(Bx, By, Bz)
        Brms     = sqrt(2.0d0 * E_now)
        max_divB = maxval(abs(divB))
        print '(A,3ES11.4)', '#     final =', &
            divB_rms/divB_rms_0, (E_now-E_ref)/E_pert0, max_divB/(Brms/dh)
    end subroutine run_glm


    subroutine compute_divB(Bx, By, Bz, divB, dx, dy, dz)
        real(8), intent(in)  :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz)
        real(8), intent(out) :: divB(nx,ny,nz)
        real(8), intent(in)  :: dx, dy, dz
        integer :: i, j, k, ip, im, jp, jm, kp, km

        do k = 1, nz
            kp = modulo(k,   nz) + 1;  km = modulo(k-2, nz) + 1
            do j = 1, ny
                jp = modulo(j,   ny) + 1;  jm = modulo(j-2, ny) + 1
                do i = 1, nx
                    ip = modulo(i,   nx) + 1;  im = modulo(i-2, nx) + 1
                    divB(i,j,k) = (Bx(ip,j,k) - Bx(im,j,k)) / (2*dx) + &
                                  (By(i,jp,k) - By(i,jm,k)) / (2*dy) + &
                                  (Bz(i,j,kp) - Bz(i,j,km)) / (2*dz)
                end do
            end do
        end do
    end subroutine compute_divB

    real(8) function rms3d(f)
        real(8), intent(in) :: f(nx,ny,nz)
        rms3d = sqrt(sum(f**2) / real(nx*ny*nz, 8))
    end function rms3d

    real(8) function mag_energy(Bx, By, Bz)
        real(8), intent(in) :: Bx(nx,ny,nz), By(nx,ny,nz), Bz(nx,ny,nz)
        mag_energy = 0.5d0 * sum(Bx**2 + By**2 + Bz**2) / real(nx*ny*nz, 8)
    end function mag_energy

end program test_divB_cleaning
