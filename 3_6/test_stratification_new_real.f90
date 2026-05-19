program test_stratification_new_real

    use ModLookUpTable,       only: ModLookUpTable_Read, nLookUpTables
    use ModConst,             only: dpi,R_sun__CGS
    use ModEOS,               only: ModEOS_init
    use ModOpacity,           only: ModOpacity_init, Opacity_table, nlogT_opacity,nlogR_opacity,&
                                    logT_opacity,logR_opacity
    use ModMath,              only: ModMath_1D_interpol_0D, ModMath_2D_interpol_1D
    use ModStratification_new, only: ModStratification_new_Init, nPoints_stratification,&
                                     r_Rsun_max_simulation,&
                                     r_Rsun_stratification, P_stratification, T_stratification,&
                                     Rho_stratification, g_stratification, gamma1_stratification,&
                                     c_sound_stratification, Xi_stratification,&
                                     Diffusion_stratification, Cooling_stratification

    implicit none

    character(len=*), parameter :: eos_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/FreeEOS_v51_x0.75_z0.02.dat'
    character(len=*), parameter :: entropy_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/FreeEOS_v51_x0.75_z0.02_entropy_contour_9.252752.dat'
    character(len=*), parameter :: opacity_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/opacity_x70_y28_z2.dat'
    integer :: i

    call ModLookUpTable_Read(eos_file, logical_unit=1)
    call ModLookUpTable_Read(entropy_file, logical_unit=1)
    call ModLookUpTable_Read(opacity_file, logical_unit=1)

    call ModEOS_init
    call ModOpacity_init
    call ModStratification_new_Init

    print *, 'nLookUpTables = ', nLookUpTables
    print *, 'nPoints_stratification = ', nPoints_stratification
    print *, 'gamma1 range = ', minval(gamma1_stratification), maxval(gamma1_stratification)
    print *, 'first 5 points: r  P  T  rho  g  c_s  Xi'
    do i = 1, min(5, nPoints_stratification)
        print '(1x,7es14.5)', r_Rsun_stratification(i), P_stratification(i), T_stratification(i), &
            Rho_stratification(i), g_stratification(i), c_sound_stratification(i), Xi_stratification(i)
    end do
    print *, 'last 5 points: r  P  T  rho  g  c_s  Xi'
    do i = max(1, nPoints_stratification-4), nPoints_stratification
        print '(1x,7es14.5)', r_Rsun_stratification(i), P_stratification(i), T_stratification(i), &
            Rho_stratification(i), g_stratification(i), c_sound_stratification(i), Xi_stratification(i)
    end do

    print *, 'at r_Rsun_max_simulation = ', r_Rsun_max_simulation
    print '(1x,a,es14.5)', '  Xi(r_Rsun_max_simulation) = ', ModMath_1D_interpol_0D(Xi_stratification,&
        nPoints_stratification,0,r_Rsun_stratification,r_Rsun_max_simulation)
    print '(1x,a,es14.5)', '  log10(rho) = ', ModMath_1D_interpol_0D(log10(Rho_stratification),&
        nPoints_stratification,0,r_Rsun_stratification,r_Rsun_max_simulation)
    print '(1x,a,es14.5)', '  log10(T)   = ', ModMath_1D_interpol_0D(log10(T_stratification),&
        nPoints_stratification,0,r_Rsun_stratification,r_Rsun_max_simulation)

    call print_opacity_probe(0.8d0)

    print *, 'total power from r=0.7 to 0.99 (erg/s)'
    print '(1x,a,es24.16)', '  diffusion = ', integrate_shell_power(Diffusion_stratification, 0.7d0, 0.99d0)
    print '(1x,a,es24.16)', '  cooling   = ', integrate_shell_power(Cooling_stratification, 0.7d0, 0.99d0)
    print '(1x,a,es24.16)', '  net       = ', integrate_shell_power(Diffusion_stratification, 0.7d0, 0.99d0) + &
        integrate_shell_power(Cooling_stratification, 0.7d0, 0.99d0)

    call print_log_profile_samples
    call print_diffusion_profile_samples
    call print_cooling_profile_samples

contains

    subroutine print_log_profile_samples
        real(8), parameter :: sample_radii(7) = [0.7d0, 0.75d0, 0.8d0, 0.85d0, 0.9d0, 0.95d0, 1.0d0]
        integer :: iSample, iPoint

        print *, 'sampled log profile points: r  log10(rho)  log10(T)'
        do iSample = 1, size(sample_radii)
            iPoint = minloc(abs(r_Rsun_stratification - sample_radii(iSample)), dim=1)
            print '(1x,3es14.5)', r_Rsun_stratification(iPoint), log10(Rho_stratification(iPoint)), log10(T_stratification(iPoint))
        end do
    end subroutine print_log_profile_samples

    subroutine print_diffusion_profile_samples
        real(8), parameter :: sample_radii(7) = [0.7d0, 0.75d0, 0.8d0, 0.85d0, 0.9d0, 0.95d0, 1.0d0]
        integer :: iSample, iPoint

        print *, 'sampled diffusion points: r  diffusion'
        do iSample = 1, size(sample_radii)
            iPoint = minloc(abs(r_Rsun_stratification - sample_radii(iSample)), dim=1)
            print '(1x,2es14.5)', r_Rsun_stratification(iPoint), Diffusion_stratification(iPoint)
        end do
    end subroutine print_diffusion_profile_samples

    subroutine print_cooling_profile_samples
        real(8), parameter :: sample_radii(5) = [0.97d0, 0.98d0, 0.985d0, 0.99d0, 1.0d0]
        integer :: iSample, iPoint

        print *, 'sampled cooling points near the boundary: r  cooling'
        do iSample = 1, size(sample_radii)
            iPoint = minloc(abs(r_Rsun_stratification - sample_radii(iSample)), dim=1)
            print '(1x,2es14.5)', r_Rsun_stratification(iPoint), Cooling_stratification(iPoint)
        end do
    end subroutine print_cooling_profile_samples

    subroutine print_opacity_probe(rProbe)
        real(8), intent(in) :: rProbe
        integer :: iPoint
        real(8) :: logTR(1,2)
        real(8) :: logKappa(1)
        real(8) :: logT_local, logR_local

        iPoint = minloc(abs(r_Rsun_stratification - rProbe), dim=1)
        logT_local = log10(T_stratification(iPoint))
        logR_local = log10(Rho_stratification(iPoint)) - 3.0d0*logT_local + 18.0d0
        logTR(1,1) = logT_local
        logTR(1,2) = logR_local
        logKappa = ModMath_2D_interpol_1D(Opacity_table%data_3D(:,:,3),&
            nlogT_opacity,nlogR_opacity,0,logT_opacity,logR_opacity,1,logTR)

        print *, 'opacity probe near r = ', rProbe, ' index = ', iPoint
        print '(1x,a,es14.5)', '  r      = ', r_Rsun_stratification(iPoint)
        print '(1x,a,es14.5)', '  logT   = ', logT_local
        print '(1x,a,es14.5)', '  logR   = ', logR_local
        print '(1x,a,es14.5)', '  log10(kappa) = ', logKappa(1)
        print '(1x,a,es14.5)', '  kappa        = ', 10.0d0**logKappa(1)
    end subroutine print_opacity_probe

    function integrate_shell_power(profile, r_min, r_max) result(total_power)
        real(8), intent(in) :: profile(:)
        real(8), intent(in) :: r_min, r_max
        real(8) :: total_power
        integer :: i, i_min, i_max
        real(8) :: shell_average

        i_min = minloc(abs(r_Rsun_stratification - r_min), dim=1)
        i_max = minloc(abs(r_Rsun_stratification - r_max), dim=1)
        total_power = 0.0d0

        do i = i_min, i_max - 1
            shell_average = 0.5d0 * (profile(i) * r_Rsun_stratification(i)**2 + &
                                    profile(i+1) * r_Rsun_stratification(i+1)**2)
            total_power = total_power + shell_average * (r_Rsun_stratification(i+1) - r_Rsun_stratification(i))
        end do

        total_power = total_power * 4.0d0 * dpi * R_sun__CGS**3
    end function integrate_shell_power

end program test_stratification_new_real