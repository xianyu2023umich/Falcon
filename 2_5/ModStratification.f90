Module ModStratification

    use ModConst,       only:   dpi,speed_c__CGS,rad_a__CGS,R_sun__CGS
    use ModParameters,  only:   ModelS_delta,ModelS_c_sound__CGS,ModelS_rmax,&
                                ModelS_dc_type,ModelS_dc_rmax
    use ModVariables,   only:   x__bar,t__bar,v__bar,&
                                rho1__bar,rho0__bar,g__bar,&
                                p1__bar,p0__bar,T0__bar,&
                                s1__bar,heat__bar,B__bar

    implicit none

    ! Declare the 1D arrays to store the profiles.
    ! All in CGS unit.

    integer             ::      ModelS_npoints
    real(8),allocatable ::      ModelS_r_list(:)                    ,&
                                ModelS_g_list__CGS(:)               ,&
                                ModelS_log_rho0_list__CGS(:)        ,&
                                ModelS_log_p0_list__CGS(:)          ,&
                                ModelS_log_T0_list(:)               ,&
                                ModelS_gamma1_list(:)               ,&
                                ModelS_gamma3_list(:)               ,&
                                ModelS_log_kap_list__CGS(:)         ,&
                                ModelS_Diffusion_list__CGS(:)       ,&
                                ModelS_cooling_list__CGS(:)         ,&
                                ModelS_Xi_list(:)                   ,& 
                                ModelS_rho0Hc_inv(:)                ,&
                                ModelS_rho0Hc_inv_int(:) 
    contains

    subroutine ModStratification_DoAll(filename)
        implicit none
        character(len=*),intent(in) ::  filename        !   filename
        
        call ModStratification_read_lookuptable(filename,logical_unit=2,head_size=1)
        call ModStratification_calc_heating
        call ModStratification_calc_Xi
        call ModStratification_set_scales
        call ModStratification_set_scale_height
        print *,ModelS_dc_rmax
    end subroutine ModStratification_DoAll

    ! Read the stratification lookuptable for the profiles.

    subroutine ModStratification_read_lookuptable(filename,logical_unit,head_size)
        implicit none
        character(len=*),intent(in) ::  filename        !   filename
        integer,intent(in)          ::  logical_unit    !   logical unit
        integer,intent(in)          ::  head_size       !   n lines of head

        character                   ::  header
        integer                     ::  ihead
        integer                     ::  iline
        integer                     ::  nvars           !   n of vars and samples

        open(unit=logical_unit, file=filename, status='old', action='read')

        ! Read the head
        do ihead=1,head_size
            read(logical_unit, '(A)') header
        end do

        ! Read data size and allocate 
        read(logical_unit,*) nvars,ModelS_npoints
        allocate(ModelS_r_list(ModelS_npoints),&
            ModelS_g_list__CGS          (ModelS_npoints),&
            ModelS_log_rho0_list__CGS   (ModelS_npoints),&
            ModelS_log_p0_list__CGS     (ModelS_npoints),&
            ModelS_log_T0_list          (ModelS_npoints),&
            ModelS_gamma1_list          (ModelS_npoints),&
            ModelS_gamma3_list          (ModelS_npoints),&
            ModelS_log_kap_list__CGS    (ModelS_npoints),&
            ModelS_Diffusion_list__CGS  (ModelS_npoints),&
            ModelS_cooling_list__CGS    (ModelS_npoints),&
            ModelS_Xi_list              (ModelS_npoints),&
            ModelS_rho0Hc_inv           (ModelS_npoints),&
            ModelS_rho0Hc_inv_int       (ModelS_npoints))

        ! Then read the file
        do iline=1,ModelS_npoints
            read(logical_unit,*)ModelS_r_list(iline),&
                                ModelS_g_list__CGS(iline),&
                                ModelS_log_rho0_list__CGS(iline),&
                                ModelS_log_p0_list__CGS(iline),&
                                ModelS_log_T0_list(iline),&
                                ModelS_gamma1_list(iline),&
                                ModelS_gamma3_list(iline),&
                                ModelS_log_kap_list__CGS(iline)
        end do

        close(logical_unit)
    end subroutine ModStratification_read_lookuptable

    ! Interface to get the variables at one height.
    ! Heatings are optional.

    subroutine ModStratification_get_vars(r,&
            g__CGS,rho0__CGS,p0__CGS,T0__CGS,gamma1,gamma3,kap__CGS,&
            diffusion__CGS,cooling__CGS,Xi_rsst)
        implicit none
        real,intent(in)             ::  r
        real,intent(out)            ::  g__CGS,rho0__CGS,p0__CGS,T0__CGS,&
                                        gamma1,gamma3,kap__CGS
        real,intent(out),optional   ::  diffusion__CGS,cooling__CGS,Xi_rsst
        
        real                        ::  posi_r,weight_r
        integer                     ::  posi_r_int

        ! Get the position of the input r value
        posi_r=(r-ModelS_r_list(1))/(ModelS_r_list(ModelS_npoints)-ModelS_r_list(1))*(ModelS_npoints-1)+1.0
        posi_r_int=floor(posi_r)
        weight_r=posi_r-posi_r_int

        ! Interpolate the table to posi_r
        ! Constant extrapolation outside two ends

        if (posi_r<1) then

            g__CGS          =       ModelS_g_list__CGS          (1)
            rho0__CGS       =10.0** ModelS_log_rho0_list__CGS   (1)
            p0__CGS         =10.0** ModelS_log_p0_list__CGS     (1)
            T0__CGS         =10.0** ModelS_log_T0_list          (1)
            gamma1          =       ModelS_gamma1_list          (1)
            gamma3          =       ModelS_gamma3_list          (1)
            kap__CGS        =10.0** ModelS_log_kap_list__CGS    (1)
            if (present(diffusion__CGS)) &
            diffusion__CGS  =       ModelS_Diffusion_list__CGS  (1)
            if (present(cooling__CGS)) &
            cooling__CGS    =       ModelS_cooling_list__CGS    (1)
            if (present(Xi_rsst)) &
            Xi_rsst         =       ModelS_Xi_list              (1)

        else if (posi_r>=ModelS_npoints) then

            g__CGS          =       ModelS_g_list__CGS          (ModelS_npoints)
            rho0__CGS       =10.0** ModelS_log_rho0_list__CGS   (ModelS_npoints)
            p0__CGS         =10.0** ModelS_log_p0_list__CGS     (ModelS_npoints)
            T0__CGS         =10.0** ModelS_log_T0_list          (ModelS_npoints)
            gamma1          =       ModelS_gamma1_list          (ModelS_npoints)
            gamma3          =       ModelS_gamma3_list          (ModelS_npoints)
            kap__CGS        =10.0** ModelS_log_kap_list__CGS    (ModelS_npoints)
            if (present(diffusion__CGS)) &
            diffusion__CGS  =       ModelS_Diffusion_list__CGS  (ModelS_npoints)
            if (present(cooling__CGS)) &
            cooling__CGS    =       ModelS_cooling_list__CGS    (ModelS_npoints)
            if (present(Xi_rsst)) &
            Xi_rsst         =       ModelS_Xi_list              (ModelS_npoints)

        else

            g__CGS          =      (ModelS_g_list__CGS          (posi_r_int)*(1.0-weight_r) + &
                                    ModelS_g_list__CGS          (posi_r_int+1)*weight_r)
            rho0__CGS       =10.0**(ModelS_log_rho0_list__CGS   (posi_r_int)*(1.0-weight_r) + &
                                    ModelS_log_rho0_list__CGS   (posi_r_int+1)*weight_r)
            p0__CGS         =10.0**(ModelS_log_p0_list__CGS     (posi_r_int)*(1.0-weight_r) + &
                                    ModelS_log_p0_list__CGS     (posi_r_int+1)*weight_r)
            T0__CGS         =10.0**(ModelS_log_T0_list          (posi_r_int)*(1.0-weight_r) + &
                                    ModelS_log_T0_list          (posi_r_int+1)*weight_r)
            gamma1          =      (ModelS_gamma1_list          (posi_r_int)*(1.0-weight_r) + &
                                    ModelS_gamma1_list          (posi_r_int+1)*weight_r)
            gamma3          =      (ModelS_gamma3_list          (posi_r_int)*(1.0-weight_r) + &
                                    ModelS_gamma3_list          (posi_r_int+1)*weight_r)
            kap__CGS        =10.0**(ModelS_log_kap_list__CGS    (posi_r_int)*(1.0-weight_r) + &
                                    ModelS_log_kap_list__CGS    (posi_r_int+1)*weight_r)

            if (present(diffusion__CGS)) &
            diffusion__CGS  =      (ModelS_diffusion_list__CGS  (posi_r_int)*(1.0-weight_r) + &
                                    ModelS_diffusion_list__CGS  (posi_r_int+1)*weight_r)
            
            if (present(cooling__CGS)) &
            cooling__CGS    =      (ModelS_cooling_list__CGS    (posi_r_int)*(1.0-weight_r) + &
                                    ModelS_cooling_list__CGS    (posi_r_int+1)*weight_r)
                            
            if (present(Xi_rsst)) &
            Xi_rsst         =      (ModelS_Xi_list              (posi_r_int)*(1.0-weight_r) + &
                                    ModelS_Xi_list              (posi_r_int+1)*weight_r)
        end if
        
    end subroutine ModStratification_get_vars

    ! Calculate the diffusion and cooling
    subroutine ModStratification_calc_heating
        implicit none
        real                        ::  g_rmax__CGS,rho0_rmax__CGS,&
                                        p0_rmax__CGS,T0_rmax__CGS,&
                                        gamma1_rmax,gamma3_rmax,&
                                        kap_rmax__CGS
        real(8)                     ::  flux(ModelS_npoints),&
                                        flux__CGS_r2(ModelS_npoints),&
                                        flux__CGS_r2_s(ModelS_npoints),&
                                        dT_dr__CGS(ModelS_npoints)

        ! As preparation get the rmax values                                            
        call ModStratification_get_vars(ModelS_rmax,g_rmax__CGS,rho0_rmax__CGS,&
            p0_rmax__CGS,T0_rmax__CGS,gamma1_rmax,gamma3_rmax,kap_rmax__CGS)
        
        ! First get the temperature gradient
        dT_dr__CGS(2:ModelS_npoints-1)=&
            (10.0**ModelS_log_T0_list(3:ModelS_npoints)-10.0**ModelS_log_T0_list(1:ModelS_npoints-2))/&
            (ModelS_r_list(3:ModelS_npoints)-ModelS_r_list(1:ModelS_npoints-2))/R_sun__CGS
        
        ! Linear extrapolation to the two ends
        dT_dr__CGS(1)=dT_dr__CGS(2)*2.0-dT_dr__CGS(3)
        dT_dr__CGS(ModelS_npoints)=dT_dr__CGS(ModelS_npoints-1)*2.0-dT_dr__CGS(ModelS_npoints-2)

        ! Then get the flux
        flux=-4./3.*speed_c__CGS*rad_a__CGS*&
            (10.0**ModelS_log_T0_list)**3*&
            dT_dr__CGS/(10.0**ModelS_log_kap_list__CGS)/&
            (10.0**ModelS_log_rho0_list__CGS)
        flux__CGS_r2=flux*ModelS_r_list**2

        ! Then get the diffusion term
        ModelS_Diffusion_list__CGS(2:ModelS_npoints-1)=-&
            (flux__CGS_r2(3:ModelS_npoints)-flux__CGS_r2(1:ModelS_npoints-2))/&
            (ModelS_r_list(3:ModelS_npoints)-ModelS_r_list(1:ModelS_npoints-2))/R_sun__CGS/&
            ModelS_r_list**2
        
        ! Linear extrapolation to the two ends
        ModelS_Diffusion_list__CGS(1)=ModelS_Diffusion_list__CGS(2)*2.0-ModelS_Diffusion_list__CGS(3)
        ModelS_Diffusion_list__CGS(ModelS_npoints)=ModelS_Diffusion_list__CGS(ModelS_npoints-1)*2.0-&
            ModelS_Diffusion_list__CGS(ModelS_npoints-2)
        
        ! Now, get the artificial cooling term
        ! Get the dc first
        if (ModelS_dc_rmax<0.0 .and. ModelS_dc_type=='scale') then
            ModelS_dc_rmax=2.0*p0_rmax__CGS/(g_rmax__CGS*rho0_rmax__CGS)/R_sun__CGS
        end if
        
        flux__CGS_r2_s=flux__CGS_r2(1)*exp(-(ModelS_r_list-ModelS_rmax)**2/ModelS_dc_rmax**2)
        ModelS_cooling_list__CGS=flux__CGS_r2(1)*&
            exp(-(ModelS_r_list-ModelS_rmax)**2/ModelS_dc_rmax**2)/&
            ModelS_r_list**2*&
            (ModelS_r_list-ModelS_rmax)/ &
            ModelS_dc_rmax**2/ &
            R_sun__CGS*2.0

        !ModelS_cooling_list__CGS(2:ModelS_npoints-1)=-&
        !    (flux__CGS_r2_s(3:ModelS_npoints)-flux__CGS_r2_s(1:ModelS_npoints-2))/&
        !    (ModelS_r_list(3:ModelS_npoints)-ModelS_r_list(1:ModelS_npoints-2))/R_sun__CGS/&
        !    ModelS_r_list**2
        
        ! Linear extrapolation to the two ends
        !ModelS_cooling_list__CGS(1)=ModelS_cooling_list__CGS(2)*2.0-ModelS_cooling_list__CGS(3)
        !ModelS_cooling_list__CGS(ModelS_npoints)=ModelS_cooling_list__CGS(ModelS_npoints-1)*2.0-&
        !    ModelS_cooling_list__CGS(ModelS_npoints-2)

        !print *,1
        !print *,ModelS_cooling_list__CGS(1:2900:10)
        !print *,ModelS_diffusion_list__CGS(1:2900:10)
        !print *,ModelS_dc_rmax
        !print *,ModelS_rmax,ModelS_cooling_list__CGS(29001),flux_r2(1)
        !print *,sum(ModelS_r_list(1:2001)**2*ModelS_cooling_list__CGS(1:2001)),sum(ModelS_r_list(1:2001)**2*ModelS_diffusion_list__CGS(1:2001))
    end subroutine ModStratification_calc_heating

    ! Calculate Xi(r)
    subroutine ModStratification_calc_Xi
        implicit none
        ModelS_Xi_list=sqrt(ModelS_gamma1_list*&
            10.0**ModelS_log_p0_list__CGS/&
            10.0**ModelS_log_rho0_list__CGS)/ModelS_c_sound__CGS
    end subroutine ModStratification_calc_Xi

    ! Set the scale for variables
    ! See my overleaf for explanation
    subroutine ModStratification_set_scales
        implicit none

        x__bar      =R_sun__CGS
        rho0__bar   =10.0**ModelS_log_rho0_list__CGS(1)
        rho1__bar   =rho0__bar*ModelS_delta
        g__bar      =ModelS_g_list__CGS(1)
        v__bar      =sqrt(ModelS_delta*x__bar*g__bar)
        t__bar      =x__bar/v__bar
        p1__bar     =rho0__bar*v__bar**2
        p0__bar     =p1__bar/ModelS_delta
        T0__bar     =10.0**ModelS_log_T0_list(1)
        s1__bar     =p1__bar/(rho0__bar*T0__bar)
        heat__bar   =rho0__bar*T0__bar*s1__bar/t__bar
        B__bar       =sqrt(4.0*dpi*rho0__bar*v__bar)
        !print *,ModelS_delta,x__bar,v__bar,t__bar,rho0__bar,rho1__bar,g__bar,p0__bar,p1__bar,T0__bar,s1__bar,heat__bar
        
    end subroutine ModStratification_set_scales

    subroutine ModStratification_set_scale_height
        implicit none

        integer                     ::  ipoint

        ModelS_rho0Hc_inv(2:ModelS_npoints-1)=&
            (   10.0**ModelS_log_rho0_list__CGS(3:ModelS_npoints)       &
            -   10.0**ModelS_log_rho0_list__CGS(1:ModelS_npoints-2) )/  &
            (   ModelS_r_list(3:ModelS_npoints)                         &
            -   ModelS_r_list(1:ModelS_npoints-2)                   )/R_sun__CGS

        ! Linear extrapolation to the two ends
        ModelS_rho0Hc_inv(1) =&
            ModelS_rho0Hc_inv(2)*2.0-ModelS_rho0Hc_inv(3)
        ModelS_rho0Hc_inv(ModelS_npoints)=&
            ModelS_rho0Hc_inv(ModelS_npoints-1)*2.0-ModelS_rho0Hc_inv(ModelS_npoints-2)

        ModelS_rho0Hc_inv=abs(ModelS_rho0Hc_inv)/10.0**ModelS_log_rho0_list__CGS

        ! Get the integral of it
        ModelS_rho0Hc_inv(1)=0

        do ipoint=2,ModelS_npoints
            ModelS_rho0Hc_inv_int(ipoint)=ModelS_rho0Hc_inv_int(ipoint-1)+&
                (ModelS_rho0Hc_inv(ipoint)+ModelS_rho0Hc_inv(ipoint-1))/2.0*&
                (ModelS_r_list(ipoint)-ModelS_r_list(ipoint-1))
        end do
    end subroutine ModStratification_set_scale_height

    function ModStratification_get_rho0HC_inv_int(r) result(Hc_inv_int)
        implicit none

        real,intent(in)             ::  r
        real(8)                     ::  Hc_inv_int
        real                        ::  posi_r,weight_r
        integer                     ::  posi_r_int

        ! Interpolate the table to posi_r
        ! Constant extrapolation outside two ends

        posi_r=(r-ModelS_r_list(1))/(ModelS_r_list(ModelS_npoints)-ModelS_r_list(1))*(ModelS_npoints-1)+1.0

        if (posi_r<1) then

            Hc_inv_int  =   ModelS_rho0HC_inv_int(1)

        else if (posi_r>=ModelS_npoints) then

            Hc_inv_int  =   ModelS_rho0HC_inv_int(ModelS_npoints)
        else
            ! Get the position of the input r value
            
            posi_r_int=floor(posi_r)
            weight_r=posi_r-posi_r_int

            Hc_inv_int  =   (ModelS_rho0HC_inv_int(posi_r_int)*(1.0-weight_r) + &
                    ModelS_rho0HC_inv_int(posi_r_int+1)*weight_r)
        end if
    end function ModStratification_get_rho0HC_inv_int

    function ModStratification_get_r_from_rho0HC_inv_int(rho0HC_inv_int) result(r)
        implicit none
        real(8),intent(in)          ::  rho0HC_inv_int
        integer                     ::  ipoint
        real(8)                     ::  r,weight_r

        r=-1.0

        if (rho0HC_inv_int < ModelS_rho0HC_inv_int(1)) then
            r=ModelS_r_list(1)
        else if (rho0HC_inv_int > ModelS_rho0HC_inv_int(ModelS_npoints)) then
            r=ModelS_r_list(ModelS_npoints)
        else
            do ipoint=1,ModelS_npoints-1
                if ((rho0HC_inv_int- ModelS_rho0HC_inv_int(ipoint))*&
                    (rho0HC_inv_int- ModelS_rho0HC_inv_int(ipoint+1))<0.0) then
                    
                    weight_r=abs(rho0HC_inv_int-ModelS_rho0HC_inv_int(ipoint))/&
                        abs(ModelS_rho0HC_inv_int(ipoint+1)-ModelS_rho0HC_inv_int(ipoint))

                    r=  ModelS_r_list(ipoint)*(1.0-weight_r)+&
                        ModelS_r_list(ipoint+1)*weight_r
                end if
            end do
        end if
    end function ModStratification_get_r_from_rho0HC_inv_int


    function ModStratification_get_middle_r(r_range) result(r_middle)
        implicit none
        real                        ::  r_middle
        real,intent(in)             ::  r_range(2)

        r_middle = ModStratification_get_r_from_rho0HC_inv_int(&
            (ModStratification_get_rho0HC_inv_int(r_range(1))+&
            ModStratification_get_rho0HC_inv_int(r_range(2)))/2.0)
    end function ModStratification_get_middle_r

end module ModStratification