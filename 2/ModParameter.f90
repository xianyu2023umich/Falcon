Module ModParameter

    ! Reynolds number

    real :: paraRe=300.

    ! Prandtl number

    real :: paraPr=1.

    ! delta_r, i.e., the delta at the bottom boundary

    real :: paraDelta_r=1.e-4

    ! xi, used for RSST (reduced speed of sound technique)

    real :: paraXi=10.

    ! Gamma

    real :: paraGamma=5./3.

    contains

    subroutine ModParameter_SetRealConst(NameConst,RealConst)

        implicit none
        
        character(len=*),intent(in) :: NameConst
        real :: RealConst

        select case(NameConst)
        case('Re')
            paraRe=RealConst
        case('Pr')
            paraPr=RealConst
        case('delta_r')
            paraDelta_r=RealConst
        end select
    end subroutine ModParameter_SetRealConst
    

end Module ModParameter