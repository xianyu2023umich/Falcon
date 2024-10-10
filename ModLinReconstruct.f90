module ModLinReconstruct

    implicit none

    contains

    subroutine ModLinReconstruct_minmod(nvars,ni,nj,nk,ng,direction,u,d_u)

        integer,intent(in)      ::  nvars,ni,nj,nk,ng,direction
        real,intent(in)         ::  u(nvars,-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        real,intent(out)        ::  d_u(nvars,-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1)
        real                    ::  u_l(nvars,-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1),&
                                    u_r(nvars,-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1),&
                                    a(nvars,-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1),&
                                    b(nvars,-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1),&
                                    c(nvars,-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1)

        ! first step: get the u at left and right
        ! for every grid we care 

        select case(direction)
        case(1)
            u_l(:,-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1)=&
                u(:,-ng+1:ni+ng-2,-ng+2:nj+ng-1,-ng+2:nk+ng-1)
            u_r(:,-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1)=&
                u(:,-ng+3:ni+ng,-ng+2:nj+ng-1,-ng+2:nk+ng-1)
        case(2)
            u_l(:,-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1)=&
                u(:,-ng+2:ni+ng-1,-ng+1:nj+ng-2,-ng+2:nk+ng-1)
            u_r(:,-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1)=&
                u(:,-ng+2:ni+ng-1,-ng+3:nj+ng,-ng+2:nk+ng-1)
        case(3)
            u_l(:,-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1)=&
                u(:,-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+1:nk+ng-2)
            u_r(:,-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1)=&
                u(:,-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+3:nk+ng)
        end select

        ! then get the three elements in minmod:
        ! \Delta_u = minmod( (u_r-u_l)/2 , 2(u_r-u) , 2(u-u_l) )
            
        a=(u_r-u_l)*0.5
        b=2.0*(u_r-u(:,-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1))
        c=2.0*(u(:,-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1)-u_l)

        ! initialize d_u and get it using minmod

        d_u=0.0
        
        where (a*b>0.0 .and. b*c>0.0)
            d_u=merge(a,b,abs(a)<abs(b))
            d_u=merge(d_u,c,abs(d_u)<abs(c))
        end where
    end subroutine ModLinReconstruct_minmod



end module ModLinReconstruct