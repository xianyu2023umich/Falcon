module ModLinReconstruct
    use ModParameters, only: ni, nj, nk, ng, nvar
    implicit none

contains

    pure elemental function minmod3(a,b,c) result(mm)
        implicit none
        real(8), intent(in) :: a,b,c
        real(8) :: mm

        if (a > 0.0d0 .and. b > 0.0d0 .and. c > 0.0d0) then
            mm = min(a,b,c)
        else if (a < 0.0d0 .and. b < 0.0d0 .and. c < 0.0d0) then
            mm = max(a,b,c)
        else
            mm = 0.0d0
        endif
    end function minmod3


    pure function minmod_slope_x(u,i,j,k,ivar) result(du)
        implicit none

        real(8), intent(in) :: u(1-ng:ni+ng, 1-ng:nj+ng, 1-ng:nk+ng, 1:nvar)
        integer, intent(in) :: i,j,k,ivar
        real(8) :: du
        real(8) :: a,b,c

        a = 0.5d0*(u(i+1,j,k,ivar) - u(i-1,j,k,ivar))
        b = 2.0d0*(u(i+1,j,k,ivar) - u(i  ,j,k,ivar))
        c = 2.0d0*(u(i  ,j,k,ivar) - u(i-1,j,k,ivar))

        du = minmod3(a,b,c)
    end function minmod_slope_x


    pure function minmod_slope_y(u,i,j,k,ivar) result(du)
        implicit none

        real(8), intent(in) :: u(1-ng:ni+ng, 1-ng:nj+ng, 1-ng:nk+ng, 1:nvar)
        integer, intent(in) :: i,j,k,ivar
        real(8) :: du
        real(8) :: a,b,c

        a = 0.5d0*(u(i,j+1,k,ivar) - u(i,j-1,k,ivar))
        b = 2.0d0*(u(i,j+1,k,ivar) - u(i,j  ,k,ivar))
        c = 2.0d0*(u(i,j  ,k,ivar) - u(i,j-1,k,ivar))

        du = minmod3(a,b,c)
    end function minmod_slope_y


    pure function minmod_slope_z(u,i,j,k,ivar) result(du)
        implicit none

        real(8), intent(in) :: u(1-ng:ni+ng, 1-ng:nj+ng, 1-ng:nk+ng, 1:nvar)
        integer, intent(in) :: i,j,k,ivar
        real(8) :: du
        real(8) :: a,b,c

        a = 0.5d0*(u(i,j,k+1,ivar) - u(i,j,k-1,ivar))
        b = 2.0d0*(u(i,j,k+1,ivar) - u(i,j,k  ,ivar))
        c = 2.0d0*(u(i,j,k  ,ivar) - u(i,j,k-1,ivar))

        du = minmod3(a,b,c)
    end function minmod_slope_z

end module ModLinReconstruct