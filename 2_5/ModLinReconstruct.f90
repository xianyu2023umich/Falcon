module ModLinReconstruct

    contains

    subroutine ModLinReconstruct_minmod(nvars,ni,nj,nk,ng,direction,u,d_u)
        implicit none
        integer,intent(in)      ::  nvars,ni,nj,nk,ng,direction
        real,intent(in)         ::  u(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,nvars)
        real,intent(out)        ::  d_u(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,nvars)
        real                    ::  u_l(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,nvars),&
                                    u_r(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,nvars),&
                                    a(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,nvars),&
                                    b(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,nvars),&
                                    c(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,nvars)

        ! first step: get the u at left and right
        ! for every grid we care 

        select case(direction)
        case(1)
            u_l(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)=&
                u(-ng+1:ni+ng-2,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)
            u_r(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)=&
                u(-ng+3:ni+ng,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)
        case(2)
            u_l(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)=&
                u(-ng+2:ni+ng-1,-ng+1:nj+ng-2,-ng+2:nk+ng-1,:)
            u_r(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)=&
                u(-ng+2:ni+ng-1,-ng+3:nj+ng,-ng+2:nk+ng-1,:)
        case(3)
            u_l(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)=&
                u(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+1:nk+ng-2,:)
            u_r(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)=&
                u(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+3:nk+ng,:)
        end select

        ! then get the three elements in minmod:
        ! \Delta_u = minmod( (u_r-u_l)/2 , 2(u_r-u) , 2(u-u_l) )
            
        a=(u_r-u_l)*0.5
        b=2.0*(u_r-u(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:))
        c=2.0*(u(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)-u_l)

        ! initialize d_u and get it using minmod

        d_u=0.5
        d_u=(sign(d_u,max(0.0,a*b)*max(0.0,b*c)*1.e10-1.e-30)+0.5)*&
            sign(min(abs(a),abs(b),abs(c)),a)
    end subroutine ModLinReconstruct_minmod

    subroutine ModLinReconstruct_minmod_new(nvars,ni,nj,nk,ng,direction,u,d_u)
        implicit none
        integer,intent(in)      ::  nvars,ni,nj,nk,ng,direction
        real,intent(in)         ::  u(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,nvars)
        real,intent(out)        ::  d_u(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,nvars)
        real                    ::  a1,b1,c1
        integer                 ::  ivar,i,j,k

        select case(direction)
        case(1)
            do ivar=1,nvars; do k=-ng+2,nk+ng-1; do j=-ng+2,nj+ng-1; do i=-ng+2,ni+ng-1
                a1=(u(i+1,j,k,ivar)-u(i-1,j,k,ivar))*0.5
                b1=2.0*(u(i+1,j,k,ivar)-u(i,j,k,ivar))
                c1=2.0*(u(i,j,k,ivar)-u(i-1,j,k,ivar))
                if (a1*b1>0.0 .and. c1*b1>0.0) then
                    d_u(i,j,k,ivar)=sign(minval(abs([a1,b1,c1])),a1)
                else
                    d_u(i,j,k,ivar)=0.
                end if
            end do; end do; end do; end do
        case(2)
            do ivar=1,nvars; do k=-ng+2,nk+ng-1; do j=-ng+2,nj+ng-1; do i=-ng+2,ni+ng-1
                a1=(u(i,j+1,k,ivar)-u(i,j-1,k,ivar))*0.5
                b1=2.0*(u(i,j+1,k,ivar)-u(i,j,k,ivar))
                c1=2.0*(u(i,j,k,ivar)-u(i,j-1,k,ivar))
                if (a1*b1>0.0 .and. c1*b1>0.0) then
                    d_u(i,j,k,ivar)=sign(minval(abs([a1,b1,c1])),a1)
                else
                    d_u(i,j,k,ivar)=0.
                end if
            end do; end do; end do; end do
        case(3)
            do ivar=1,nvars; do k=-ng+2,nk+ng-1; do j=-ng+2,nj+ng-1; do i=-ng+2,ni+ng-1
                a1=(u(i,j,k+1,ivar)-u(i,j,k-1,ivar))*0.5
                b1=2.0*(u(i,j,k+1,ivar)-u(i,j,k,ivar))
                c1=2.0*(u(i,j,k,ivar)-u(i,j,k-1,ivar))
                if (a1*b1>0.0 .and. c1*b1>0.0) then
                    d_u(i,j,k,ivar)=sign(minval(abs([a1,b1,c1])),a1)
                else
                    d_u(i,j,k,ivar)=0.
                end if
            end do; end do; end do; end do
        end select
    end subroutine ModLinReconstruct_minmod_new
end module ModLinReconstruct