Module ModSSM_v0

    !real :: ModSSM_rho_inner,ModSSM_p_inner,ModSSM_T_inner,ModSSM_delta_inner

    contains

    function ModSSM_v0_get_var0(nz,z,m,NameVar) result(var0)

        integer,intent(in) :: nz
        real,intent(in) :: z(nz)
        real,intent(in) :: m
        character(len=*),intent(in) :: NameVar
        
        real :: poly(nz), var0(nz)

        poly=1.-z/(m+1)
        select case(NameVar)
        case('rho0')
            var0=poly**m
        case('p0')
            var0=poly**(m+1)
        case('te0')
            var0=poly
        end select
    end function ModSSM_v0_get_var0

    !function ModSSM_v0_get_s0
    !end function ModSSM_v0_get_s0

end Module ModSSM_v0