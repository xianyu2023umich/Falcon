program test_value_locate_1D

    real(8) :: x(5)
    real(8) :: fout(3)
    x = -[1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0]
    print *, ModMath_value_locate_1D(x, 3.6d0)

    call ModMath_1D_interpol_1D(x**2,5,0,x,3,-[0.5d0,2.5d0,5.5d0],fout)

    print *, fout

    contains

    subroutine ModMath_1D_interpol_1D(f,ni,ng,xi,nout,xi_out,fout)

        implicit none

        integer,intent(in) :: ni,ng,nout
        real(8),intent(in) :: f(-ng+1:ni+ng),xi(-ng+1:ni+ng),xi_out(nout)
        real(8),intent(out):: fout(nout)

        real(8) :: w_i
        integer :: posi_i_integer
        integer :: i

        do i=1,nout
            posi_i_integer = ModMath_value_locate_1D(xi,xi_out(i))
            w_i = (xi_out(i) - xi(posi_i_integer)) / (xi(posi_i_integer+1) - xi(posi_i_integer))
            fout(i) = f(posi_i_integer) * (1.0d0 - w_i) + f(posi_i_integer+1) * w_i
        end do
    end subroutine ModMath_1D_interpol_1D

    function ModMath_value_locate_1D(x, xv) result(i)
        real(8), intent(in) :: x(:)
        real(8), intent(in) :: xv
        integer :: lo, hi, mid, n
        logical :: ascending

        n = size(x)
        if (n < 2) then
            i = 0         ! or some error code
            return
        end if

        ascending = (x(n) > x(1))

        ! Handle outside range
        if (ascending) then
            if (xv <= x(1)) then
            i = 1; return
            else if (xv >= x(n)) then
            i = n-1; return
            end if
        else
            if (xv >= x(1)) then
            i = 1; return
            else if (xv <= x(n)) then
            i = n-1; return
            end if
        end if

        lo = 1
        hi = n

        do while (hi - lo > 1)
            mid = (lo + hi)/2
            if (ascending) then
            if (xv >= x(mid)) then
                lo = mid
            else
                hi = mid
            end if
            else
            if (xv <= x(mid)) then
                lo = mid
            else
                hi = mid
            end if
            end if
        end do

        i = lo
    end function ModMath_value_locate_1D
end program test_value_locate_1D