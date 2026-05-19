module ModStratification_new

    use ModConst,       only: dpi,G_CGS,rad_a__CGS,speed_c__CGS,R_sun__CGS
    use ModLookUpTable, only: LookUpTables,LookUpTable,nLookUpTables
    use ModEOS,         only: EOS_table,nlogQ_EOS_table,nlogT_EOS_table,logQ_EOS_table,logT_EOS_table
    use ModOpacity,     only: Opacity_table,nlogT_opacity,nlogR_opacity,logT_opacity,logR_opacity
    use ModMath,        only: ModMath_1D_interpol_0D,ModMath_1D_interpol_1D,&
                              ModMath_2D_interpolate_1D,ModMath_2D_interpol_1D,&
                              ModMath_value_locate_1D

    implicit none

    type(LookUpTable),pointer ::  entropy_contour
    integer             ::  nEntropyContour
    

    real(8)             ::  logS
    real(8),allocatable ::  logRho_of_entropy_contour(:)
    real(8),allocatable ::  logQ_of_entropy_contour(:)
    real(8),allocatable ::  logT_of_entropy_contour(:)
    real(8),allocatable ::  logQ_T_of_entropy_contour(:,:)
    real(8),allocatable ::  logT_R_of_entropy_contour(:,:)
    real(8),allocatable ::  logP_of_entropy_contour(:)
    

    ! The values at the base of the convection zone.
    real(8)             ::  r_Rsun_base=0.7
    real(8)             ::  Rho_base=10.0**(-0.6841807359334023)
    real(8)             ::  g_base=54340.03688603464
    real(8)             ::  P_base

    ! stratification profiles
    real(8)             ::  dr_Rsun=0.0001
    real(8)             ::  r_Rsun_max_profile=1.0025
    real(8)             ::  r_Rsun_max_cooling=0.99d0
    real(8)             ::  r_Rsun_max_simulation=0.99
    integer             ::  nPoints_stratification
    real(8),allocatable ::  r_Rsun_stratification(:)
    real(8),allocatable ::  Rho_stratification(:)
    real(8),allocatable ::  P_stratification(:)
    real(8),allocatable ::  T_stratification(:)
    real(8),allocatable ::  g_stratification(:)
    real(8),allocatable ::  gamma1_stratification(:)
    real(8),allocatable ::  gamma3_stratification(:)
    real(8),allocatable ::  c_sound_stratification(:)
    real(8),allocatable ::  Xi_stratification(:)
    real(8),allocatable ::  Kappa_stratification(:)
    real(8),allocatable ::  BlackBodyFlux_stratification(:)
    real(8),allocatable ::  CoolingFlux_stratification(:)
    real(8),allocatable ::  Diffusion_stratification(:)
    real(8),allocatable ::  Cooling_stratification(:)

    contains

    ! This subroutine initializes the stratification module.
    ! It first determines which are the entropy contour and EOS table.
    ! Then it gets p(rho), which is used to calculate the stratification profiles.

    subroutine ModStratification_new_Init
        implicit none
        integer :: iLookUpTable

        ! Determine which are the entropy contour and EOS table.

        do iLookUpTable=1,nLookUpTables
            if (LookUpTables(iLookUpTable)%name == 'Constant_entropy_contour') then
                entropy_contour => LookUpTables(iLookUpTable)
                nEntropyContour = entropy_contour%sizes(1)

                call ModStratification_new_nudge_entropy_contour_off_upper_edge

                ! Get logQ, logT, and logRho of the entropy contour.
                allocate(logRho_of_entropy_contour(nEntropyContour))
                allocate(logP_of_entropy_contour(nEntropyContour))
                allocate(logQ_of_entropy_contour(nEntropyContour))
                allocate(logT_of_entropy_contour(nEntropyContour))
                allocate(logQ_T_of_entropy_contour(nEntropyContour,2))
                allocate(logT_R_of_entropy_contour(nEntropyContour,2))

                logQ_of_entropy_contour = entropy_contour%data_2D(:,1)
                logT_of_entropy_contour = entropy_contour%data_2D(:,2)
                logQ_T_of_entropy_contour = entropy_contour%data_2D(:,:)
                logRho_of_entropy_contour = logQ_of_entropy_contour + 2*logT_of_entropy_contour - 12

                ! R is for opacity table. logR = logRho - 3*logT + 18.
                logT_R_of_entropy_contour(:,1) = logT_of_entropy_contour
                logT_R_of_entropy_contour(:,2) = logRho_of_entropy_contour - 3.0*logT_of_entropy_contour + 18.0
                exit
            end if
        end do

        call ModStratification_new_get_pressure_density
        call ModStratification_new_allocate
        call ModStratification_new_calc_profiles
        call ModStratification_new_calc_opacity
        call ModStratification_new_calc_Xi
        call ModStratification_new_calc_gamma3
        call ModStratification_new_calc_heating
        call ModStratification_new_artificial_cooling
    end subroutine ModStratification_new_Init

    subroutine ModStratification_new_nudge_entropy_contour_off_upper_edge
        implicit none

        if (allocated(entropy_contour%data_2D)) then
            entropy_contour%data_2D(:,2) = min(entropy_contour%data_2D(:,2), 8.199999d0)
        end if
    end subroutine ModStratification_new_nudge_entropy_contour_off_upper_edge

    ! This subroutine gets pressure as a function of density for the entropy contour.
    ! The first column of entropy contour is logQ, second logT. Entropy is constant.
    ! logQ=logRho - 2 * logT + 12.
    ! So we should first get logRho directly from logQ and logT. Then we should
    ! interoplate EOS table with logQ and logT to get pressure.

    subroutine ModStratification_new_get_pressure_density
        implicit none

        ! Interpolate EOS table with logQ and logT to get pressure.
        logP_of_entropy_contour = ModMath_2D_interpol_1D(EOS_table%data_3D(:,:,3),&
            nlogQ_EOS_table,nlogT_EOS_table,0,logQ_EOS_table,logT_EOS_table,&
            nEntropyContour,logQ_T_of_entropy_contour)
    end subroutine ModStratification_new_get_pressure_density

    subroutine ModStratification_new_allocate
        implicit none
        ! Get n_points
        nPoints_stratification = floor((r_Rsun_max_profile - r_Rsun_base) / dr_Rsun) + 1

        ! Allocate arrays
        allocate(r_Rsun_stratification(nPoints_stratification))
        allocate(Rho_stratification(nPoints_stratification  ))
        allocate(P_stratification(nPoints_stratification))
        allocate(T_stratification(nPoints_stratification))
        allocate(g_stratification(nPoints_stratification))
        allocate(gamma1_stratification(nPoints_stratification))
        allocate(gamma3_stratification(nPoints_stratification))
        allocate(c_sound_stratification(nPoints_stratification))
        allocate(Xi_stratification(nPoints_stratification))
        allocate(BlackBodyFlux_stratification(nPoints_stratification))
        allocate(CoolingFlux_stratification(nPoints_stratification))
        allocate(Diffusion_stratification(nPoints_stratification))
        allocate(Cooling_stratification(nPoints_stratification))
        allocate(Kappa_stratification(nPoints_stratification))
    end subroutine ModStratification_new_allocate

    ! This subroutine calculates the stratification profiles.
    ! Here r is dimensionless (in units of R_sun__CGS), so the derivatives
    ! need the solar-radius conversion factor.
    ! d(logP)/dr = -R_sun__CGS * rho * g / P
    ! dg/dr = -2 * g / r + 4 * pi * G * rho * R_sun__CGS

    subroutine ModStratification_new_calc_profiles
        implicit none
        integer :: iPoint
        real(8) :: k1(2),k2(2),k3(2),k4(2)
        real(8) :: r_Rsun_current,P_current,g_current,rho_current

        ! Set the base pressure
        P_base = 10.0d0**ModMath_1D_interpol_0D(logP_of_entropy_contour,&
            nEntropyContour,0,logRho_of_entropy_contour,log10(rho_base))

        ! Iterate P, g and r. So initialize them.

        r_Rsun_current = r_Rsun_base
        g_current = g_base
        P_current = P_base
        iPoint = 1

        ! Calculate the stratification profiles using rk4.
        do while (r_Rsun_current < r_Rsun_max_profile)
            r_Rsun_stratification(iPoint) = r_Rsun_current
            P_stratification(iPoint) = P_current
            g_stratification(iPoint) = g_current

            rho_current = 10.0d0**ModMath_1D_interpol_0D(logRho_of_entropy_contour,&
                nEntropyContour,0,logP_of_entropy_contour,log10(P_current))

            k1 = ModStratification_new_rk4_derivative(P_current,g_current,r_Rsun_current)
            k2 = ModStratification_new_rk4_derivative(P_current + 0.5*dr_Rsun*k1(1),&
                g_current + 0.5*dr_Rsun*k1(2), r_Rsun_current + 0.5*dr_Rsun)
            k3 = ModStratification_new_rk4_derivative(P_current + 0.5*dr_Rsun*k2(1),&
                g_current + 0.5*dr_Rsun*k2(2), r_Rsun_current + 0.5*dr_Rsun)
            k4 = ModStratification_new_rk4_derivative(P_current + dr_Rsun*k3(1),&
                g_current + dr_Rsun*k3(2), r_Rsun_current + dr_Rsun)
            
            P_current = P_current + (dr_Rsun/6.0)*(k1(1) + 2.0*k2(1) + 2.0*k3(1) + k4(1))
            g_current = g_current + (dr_Rsun/6.0)*(k1(2) + 2.0*k2(2) + 2.0*k3(2) + k4(2))
            r_Rsun_current = r_Rsun_current + dr_Rsun
            
            ! Store the profiles
            iPoint = iPoint + 1
        end do

        ! Use interpol to get logRho and logT (thus Rho and T) profiles.

        Rho_stratification=10.0**&
            ModMath_1D_interpol_1D(logRho_of_entropy_contour,&
            nEntropyContour,0,logP_of_entropy_contour,&
            nPoints_stratification,log10(P_stratification))

        T_stratification=10.0**&
            ModMath_1D_interpol_1D(logT_of_entropy_contour,&
            nEntropyContour,0,logP_of_entropy_contour,&
            nPoints_stratification,log10(P_stratification))
    end subroutine ModStratification_new_calc_profiles

    subroutine ModStratification_new_calc_opacity
        implicit none
        real(8) :: logT_R_stratification(nPoints_stratification,2)
        
        ! Interpolate opacity table with logT and logR to get Kappa.
        logT_R_stratification(:,1) = log10(T_stratification)
        logT_R_stratification(:,2) = log10(Rho_stratification) - 3.0d0*logT_R_stratification(:,1) + 18.0d0

        Kappa_stratification = ModMath_2D_interpol_1D(Opacity_table%data_3D(:,:,3),&
            nlogT_opacity,nlogR_opacity,0,logT_opacity,logR_opacity,&
            nPoints_stratification,logT_R_stratification)
        Kappa_stratification=10.0**Kappa_stratification
    end subroutine ModStratification_new_calc_opacity

    ! To run rk4 I need a function to get the RHS of following equations
    ! based on current P, g, and r.
    ! dP/dr = -(R_sun__CGS * rho * g)
    ! dg/dr = -2 * g / r + 4 * pi * G * rho * R_sun__CGS

    function ModStratification_new_rk4_derivative(p,g,r) result(derivatives)
        implicit none
        real(8),intent(in)  :: p,g,r
        real(8)             :: rho
        real(8)             :: derivatives(2)
        
        rho = 10.d0**ModMath_1D_interpol_0D(logRho_of_entropy_contour,&
            nEntropyContour,0,logP_of_entropy_contour,log10(p))
        
        derivatives(1) = -(R_sun__CGS * rho * g)
        derivatives(2) = -2.0 * g / r + 4.0 * dpi * G_CGS * rho * R_sun__CGS
    end function ModStratification_new_rk4_derivative

    subroutine ModStratification_new_calc_heating
        implicit none
        real(8)                     ::  flux_r2(nPoints_stratification),&
                                        dT_dr(nPoints_stratification)

        ! First get the temperature gradient
        dT_dr(2:nPoints_stratification-1)=&
            (T_stratification(3:nPoints_stratification)-T_stratification(1:nPoints_stratification-2))/&
            (r_Rsun_stratification(3:nPoints_stratification)-r_Rsun_stratification(1:nPoints_stratification-2))/R_sun__CGS

        ! Linear extrapolation to the two ends
        dT_dr(1)=dT_dr(2)*2.0-dT_dr(3)
        dT_dr(nPoints_stratification)=dT_dr(nPoints_stratification-1)*2.0-dT_dr(nPoints_stratification-2)

        ! Then get the flux
        BlackBodyFlux_stratification=&
            -1.333333333333*speed_c__CGS*rad_a__CGS*&
            (T_stratification(:))**3*&
            dT_dr/(Kappa_stratification*Rho_stratification)
        flux_r2=BlackBodyFlux_stratification*r_Rsun_stratification**2

        ! Then get the diffusion term
        Diffusion_stratification(2:nPoints_stratification-1)=-&
            (flux_r2(3:nPoints_stratification)-flux_r2(1:nPoints_stratification-2))/&
            (r_Rsun_stratification(3:nPoints_stratification)-r_Rsun_stratification(1:nPoints_stratification-2))/R_sun__CGS/&
            r_Rsun_stratification(2:nPoints_stratification-1)**2

        ! Linear extrapolation to the two ends
        Diffusion_stratification(1)=Diffusion_stratification(2)*2.0-Diffusion_stratification(3)
        Diffusion_stratification(nPoints_stratification)=Diffusion_stratification(nPoints_stratification-1)*2.0-&
            Diffusion_stratification(nPoints_stratification-2)
    end subroutine ModStratification_new_calc_heating

    ! This subroutine calculates Xi for reduced speed of sound technique.
    ! See Hotta et al. 2012 or 2014 for more.

    subroutine ModStratification_new_calc_Xi
        implicit none
        real(8) :: c_sound_r_max_simulation

        ! Gamma1 profile using interpolation of EOS table.
        Gamma1_stratification=&
            ModMath_2D_interpol_1D(EOS_table%data_3D(:,:,16),&
            nlogQ_EOS_table,nlogT_EOS_table,0,logQ_EOS_table,logT_EOS_table,&
            nPoints_stratification,logQ_T_of_entropy_contour)

        ! Get the speed of sound profile.
        ! c_s = sqrt(gamma1 * P / rho)
        c_sound_stratification=sqrt(gamma1_stratification*&
            P_stratification/Rho_stratification)

        ! Get c_sound at r_max_simulation
        c_sound_r_max_simulation = &
            ModMath_1D_interpol_0D(c_sound_stratification,&
            nPoints_stratification,0,r_Rsun_stratification,r_Rsun_max_simulation)

        Xi_stratification=c_sound_stratification/c_sound_r_max_simulation
    end subroutine ModStratification_new_calc_Xi

    subroutine ModStratification_new_calc_gamma3
        implicit none
        real(8) :: logQ_T_stratification(nPoints_stratification,2)

        logQ_T_stratification(:,1) = log10(Rho_stratification) - 2.0d0*log10(T_stratification) + 12.0d0
        logQ_T_stratification(:,2) = log10(T_stratification)

        gamma3_stratification = ModMath_2D_interpol_1D(EOS_table%data_3D(:,:,17),&
            nlogQ_EOS_table,nlogT_EOS_table,0,logQ_EOS_table,logT_EOS_table,&
            nPoints_stratification,logQ_T_stratification)
    end subroutine ModStratification_new_calc_gamma3

    subroutine ModStratification_new_artificial_cooling
        implicit none
        real(8) :: cooling_center
        real(8) :: cooling_width
        real(8) :: p_center
        real(8) :: rho_center
        real(8) :: g_center
        real(8) :: cooling_flux_0
        real(8) :: cooling_flux_r2(nPoints_stratification)

        cooling_center = r_Rsun_max_cooling
        if (cooling_center < 0.0d0) cooling_center = r_Rsun_max_simulation

        p_center = ModMath_1D_interpol_0D(P_stratification,&
            nPoints_stratification,0,r_Rsun_stratification,cooling_center)
        rho_center = ModMath_1D_interpol_0D(Rho_stratification,&
            nPoints_stratification,0,r_Rsun_stratification,cooling_center)
        g_center = ModMath_1D_interpol_0D(g_stratification,&
            nPoints_stratification,0,r_Rsun_stratification,cooling_center)

        ! Use the same scale-height style width as the old stratification module.
        cooling_width = 2.0d0 * p_center / (g_center * rho_center) / R_sun__CGS

        ! Anchor the cooling profile to the surface diffusion flux.
        cooling_flux_0 = BlackBodyFlux_stratification(1)

        cooling_flux_r2 = cooling_flux_0 * &
            exp(-(r_Rsun_stratification-cooling_center)**2 / cooling_width**2) * &
            r_Rsun_stratification(1)**2 / r_Rsun_stratification**2

        Cooling_stratification = cooling_flux_r2 * &
            (r_Rsun_stratification-cooling_center) / cooling_width**2 / R_sun__CGS * 2.0d0
        CoolingFlux_stratification = cooling_flux_r2

    end subroutine ModStratification_new_artificial_cooling

    subroutine ModStratification_new_get_vars(r,&
            g__CGS,rho0__CGS,p0__CGS,T0__CGS,gamma1,gamma3,kap__CGS,&
            diffusion__CGS,cooling__CGS,diffusion_flux__CGS,cooling_flux__CGS,Xi_rsst)
        implicit none
        real(8),intent(in)           :: r
        real(8),intent(out)          :: g__CGS,rho0__CGS,p0__CGS,T0__CGS,&
                                        gamma1,gamma3,kap__CGS
        real(8),intent(out),optional :: diffusion__CGS,cooling__CGS,diffusion_flux__CGS,cooling_flux__CGS,Xi_rsst
        real(8)                      :: r_rsun

        ! The stratification profiles are tabulated in solar-radius units,
        ! while callers pass CGS radii.
        r_rsun = r / R_sun__CGS

        g__CGS    = ModMath_1D_interpol_0D(g_stratification,     nPoints_stratification,0,r_Rsun_stratification,r_rsun)
        rho0__CGS = ModMath_1D_interpol_0D(Rho_stratification,   nPoints_stratification,0,r_Rsun_stratification,r_rsun)
        p0__CGS   = ModMath_1D_interpol_0D(P_stratification,     nPoints_stratification,0,r_Rsun_stratification,r_rsun)
        T0__CGS   = ModMath_1D_interpol_0D(T_stratification,     nPoints_stratification,0,r_Rsun_stratification,r_rsun)
        gamma1    = ModMath_1D_interpol_0D(gamma1_stratification, nPoints_stratification,0,r_Rsun_stratification,r_rsun)
        gamma3    = ModMath_1D_interpol_0D(gamma3_stratification, nPoints_stratification,0,r_Rsun_stratification,r_rsun)
        kap__CGS  = ModMath_1D_interpol_0D(Kappa_stratification, nPoints_stratification,0,r_Rsun_stratification,r_rsun)

        if (present(diffusion__CGS))      diffusion__CGS      = ModMath_1D_interpol_0D(Diffusion_stratification,      nPoints_stratification,0,r_Rsun_stratification,r_rsun)
        if (present(cooling__CGS))        cooling__CGS        = ModMath_1D_interpol_0D(Cooling_stratification,        nPoints_stratification,0,r_Rsun_stratification,r_rsun)
        if (present(diffusion_flux__CGS)) diffusion_flux__CGS = ModMath_1D_interpol_0D(BlackBodyFlux_stratification, nPoints_stratification,0,r_Rsun_stratification,r_rsun)
        if (present(cooling_flux__CGS))   cooling_flux__CGS   = ModMath_1D_interpol_0D(CoolingFlux_stratification,   nPoints_stratification,0,r_Rsun_stratification,r_rsun)
        if (present(Xi_rsst))             Xi_rsst             = ModMath_1D_interpol_0D(Xi_stratification,             nPoints_stratification,0,r_Rsun_stratification,r_rsun)
    end subroutine ModStratification_new_get_vars

    ! Find the r that bisects ∫|d(ln ρ₀)/dr| dr over r_range (both in CGS).
    ! This places AMR splits at acoustically equal intervals rather than geometric ones.
    function ModStratification_new_get_middle_r(r_range) result(r_middle)
        implicit none
        real(8), intent(in) :: r_range(2)
        real(8)             :: r_middle

        integer, parameter  :: n = 200
        real(8)             :: dr, r
        real(8)             :: rho0_arr(n+1), f(n+1), cumul(n+1)
        real(8)             :: g_d, p_d, T_d, g1_d, g3_d, k_d   ! dummy outputs
        real(8)             :: half_total
        integer             :: i

        dr = (r_range(2) - r_range(1)) / dble(n)

        ! Sample rho0 on uniform grid
        do i = 1, n+1
            r = r_range(1) + (i-1)*dr
            call ModStratification_new_get_vars(r, g_d, rho0_arr(i), p_d, T_d, g1_d, g3_d, k_d)
        end do

        ! f = |d(ln rho0)/dr| via finite differences
        f(1)   = abs(rho0_arr(2)   - rho0_arr(1)  ) / (dr       * rho0_arr(1)  )
        f(n+1) = abs(rho0_arr(n+1) - rho0_arr(n)  ) / (dr       * rho0_arr(n+1))
        do i = 2, n
            f(i) = abs(rho0_arr(i+1) - rho0_arr(i-1)) / (2.0d0*dr * rho0_arr(i)  )
        end do

        ! Cumulative trapezoid integral
        cumul(1) = 0.0d0
        do i = 2, n+1
            cumul(i) = cumul(i-1) + 0.5d0*(f(i)+f(i-1))*dr
        end do
        half_total = 0.5d0 * cumul(n+1)

        ! Find index where cumulative crosses half_total, then interpolate
        r_middle = r_range(2)
        do i = 2, n+1
            if (cumul(i) >= half_total) then
                r_middle = (r_range(1) + (i-2)*dr) + &
                           (half_total - cumul(i-1)) / (cumul(i) - cumul(i-1)) * dr
                exit
            end if
        end do

    end function ModStratification_new_get_middle_r

end module ModStratification_new