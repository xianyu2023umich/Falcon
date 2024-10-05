Module ModVariable

    contains

    function ModVariable_SetNvar(NameEqn) result(nvar)

        implicit none

        character(len=*),intent(in) :: NameEqn
        integer :: nvar

        select case(NameEqn)
        case('Dynamo_HD')
            nvar=5
        end select
    end function ModVariable_SetNvar

    function ModVariable_SetVarIndex(NameEqn,NameVar) result(index)

        implicit none

        character(len=*),intent(in) :: NameEqn,NameVar
        integer :: index

        index=-1

        select case(NameEqn)
        case('dynamo_HD')
            select case(NameVar)
            case('rho1')
                index=1
            case('vi')
                index=2
            case('vj')
                index=3
            case('vk')
                index=4
            case('s1')
                index=5
            case default
            end select
        case default
        end select

    end function ModVariable_SetVarIndex

    function ModVariable_SetVarScale(NameEqn,NameVar) result(scale)

        use ModParameter, only : paraDelta_r

        implicit none

        character(len=*),intent(in) :: NameEqn,NameVar
        real :: scale

        scale=1.

        select case(NameEqn)
        case('dynamo_HD')
            select case(NameVar)
            case('x')
                scale=1.
            case('v')
                scale=1.
            case('t')
                scale=1.
            case('rho0')
                scale=1.
            case('p0')
                scale=1.
            case('s0')
                scale=1.
            case('rho1')
                scale=8.*paraDelta_r
            case('p1')
                scale=8.*paraDelta_r
            case('s1')
                scale=8.*paraDelta_r
            case default
            end select
        case default
        end select

    end function ModVariable_SetVarScale


end Module ModVariable