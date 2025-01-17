module ModMath

    contains

    function ModMath_IfLinesInterSect(range1,range2) result(if_intersect)
        implicit none
        real,intent(in)     :: range1(2),range2(2)
        logical             :: if_intersect
        
        if_intersect=.false.

        if((range1(1)-range2(1))*(range1(1)-range2(2))<0.0 .or. &
           (range1(2)-range2(1))*(range1(2)-range2(2))<0.0 .or. &
           (range2(1)-range1(1))*(range2(1)-range1(2))<0.0 .or. &
           (range2(2)-range1(1))*(range2(2)-range1(2))<0.0) if_intersect=.true.
    end function

    function ModMath_IfBlocksInterSect(xijk_range1,xijk_range2) result(if_intersect)
        implicit none
        real,intent(in)     ::  xijk_range1(3,2),xijk_range2(3,2)   ! the xijk_ranges of two blocks
        logical             ::  if_intersect                        ! the output
        
        integer             ::  direction

        do direction=1,3
            if_intersect = ModMath_IfLinesInterSect(xijk_range1(direction,:),&
                xijk_range2(direction,:))
            if(.not.if_intersect)exit
        end do
    end function ModMath_IfBlocksInterSect

    function ModMath_if_belong_1D(x,f,posi) result(if_belong)

        integer,intent(in) :: f(:),x

        integer,intent(out) :: posi
        integer :: i
        logical :: if_belong

        if_belong=.false.
        posi=-1

        do i=1,size(f)
            if (f(i)==x) then
                if_belong=.true.
                posi=i
                exit
            else
            end if
        end do
    end function ModMath_if_belong_1D

    subroutine ModMath_3D_interpolate_1D(f,ni,nj,nk,ng,xi,xj,xk,nout,xijk_out,fout,debug)

        implicit none

        integer,intent(in) :: ni,nj,nk,ng,nout
        real,intent(in) :: xi(-ng+1:ni+ng),xj(-ng+1:nj+ng),xk(-ng+1:nk+ng),f(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        real,intent(in) :: xijk_out(nout,3)

        real :: fout(nout)

        real :: dxi,dxj,dxk,posi_i,posi_j,posi_k,w_i,w_j,w_k
        integer :: posi_i_integer,posi_j_integer,posi_k_integer
        integer :: i
        logical,optional :: debug

        dxi=xi(2)-xi(1); dxj=xj(2)-xj(1); dxk=xk(2)-xk(1)

        do i=1,nout

            posi_i=-ng+1.+(xijk_out(i,1)-xi(-ng+1))/dxi
            posi_j=-ng+1.+(xijk_out(i,2)-xj(-ng+1))/dxj
            posi_k=-ng+1.+(xijk_out(i,3)-xk(-ng+1))/dxk

            posi_i_integer=floor(posi_i)
            posi_j_integer=floor(posi_j)
            posi_k_integer=floor(posi_k)

            w_i=posi_i-posi_i_integer
            w_j=posi_j-posi_j_integer
            w_k=posi_k-posi_k_integer

            fout(i)=&
                f(posi_i_integer  ,posi_j_integer  ,posi_k_integer  )*(1.-w_i)*(1.-w_j)*(1.-w_k)+&
                f(posi_i_integer+1,posi_j_integer  ,posi_k_integer  )*(   w_i)*(1.-w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer+1,posi_k_integer  )*(1.-w_i)*(   w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer  ,posi_k_integer+1)*(1.-w_i)*(1.-w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer+1,posi_k_integer  )*(   w_i)*(   w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer+1,posi_k_integer+1)*(1.-w_i)*(   w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer  ,posi_k_integer+1)*(   w_i)*(1.-w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer+1,posi_k_integer+1)*(   w_i)*(   w_j)*(   w_k)
            
            if (present(debug)) then
                !print *,-ng+1.+(xijk_out(i,3)-xk(-ng+1))/dxk
                !print *,xk
            end if
        end do
    end subroutine ModMath_3D_interpolate_1D

    subroutine ModMath_1D3D_interpolate_1D1D(f,nvar,ni,nj,nk,ng,xi,xj,xk,nout,xijk_out,fout)

        implicit none

        integer,intent(in) :: nvar,ni,nj,nk,ng,nout
        real,intent(in) :: xi(-ng+1:ni+ng),xj(-ng+1:nj+ng),xk(-ng+1:nk+ng),f(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,nvar)
        real,intent(in) :: xijk_out(nout,3)

        real :: fout(nout,nvar)

        real :: dxi,dxj,dxk,posi_i,posi_j,posi_k,w_i,w_j,w_k
        integer :: posi_i_integer,posi_j_integer,posi_k_integer
        integer :: i

        dxi=xi(2)-xi(1); dxj=xj(2)-xj(1); dxk=xk(2)-xk(1)

        do i=1,nout

            posi_i=-ng+1.+(xijk_out(i,1)-xi(-ng+1))/dxi
            posi_j=-ng+1.+(xijk_out(i,2)-xj(-ng+1))/dxj
            posi_k=-ng+1.+(xijk_out(i,3)-xk(-ng+1))/dxk
            posi_i=max(min(posi_i,ni+0.0),1.0)
            posi_j=max(min(posi_j,nj+0.0),1.0)
            posi_k=max(min(posi_k,nk+0.0),1.0)

            posi_i_integer=floor(posi_i)
            posi_j_integer=floor(posi_j)
            posi_k_integer=floor(posi_k)

            w_i=posi_i-posi_i_integer
            w_j=posi_j-posi_j_integer
            w_k=posi_k-posi_k_integer

            fout(i,:)=&
                f(posi_i_integer  ,posi_j_integer  ,posi_k_integer  ,:)*(1.-w_i)*(1.-w_j)*(1.-w_k)+&
                f(posi_i_integer+1,posi_j_integer  ,posi_k_integer  ,:)*(   w_i)*(1.-w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer+1,posi_k_integer  ,:)*(1.-w_i)*(   w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer  ,posi_k_integer+1,:)*(1.-w_i)*(1.-w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer+1,posi_k_integer  ,:)*(   w_i)*(   w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer+1,posi_k_integer+1,:)*(1.-w_i)*(   w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer  ,posi_k_integer+1,:)*(   w_i)*(1.-w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer+1,posi_k_integer+1,:)*(   w_i)*(   w_j)*(   w_k)
        end do
    end subroutine ModMath_1D3D_interpolate_1D1D

    subroutine ModMath_1D3D_interpolate_1D2D(f,nvar,ni,nj,nk,ng,xi,xj,xk,nij_out,xi_out,xj_out,xk_out,fout)

        implicit none

        integer,intent(in) :: nvar,ni,nj,nk,ng
        integer,intent(in) :: nij_out(2)
        real,intent(in) :: xi(-ng+1:ni+ng),xj(-ng+1:nj+ng),xk(-ng+1:nk+ng),f(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,nvar)
        real,intent(in) :: xi_out(nij_out(1),nij_out(2)),xj_out(nij_out(1),nij_out(2)),xk_out(nij_out(1),nij_out(2))

        real :: fout(nvar,nij_out(1),nij_out(2))

        real :: dxi,dxj,dxk,posi_i,posi_j,posi_k,w_i,w_j,w_k
        integer :: posi_i_integer,posi_j_integer,posi_k_integer
        integer :: i_out,j_out

        dxi=xi(2)-xi(1); dxj=xj(2)-xj(1); dxk=xk(2)-xk(1)

        do i_out=1,nij_out(1); do j_out=1,nij_out(2)

            posi_i=-ng+1.+(xi_out(i_out,j_out)-xi(-ng+1))/dxi
            posi_j=-ng+1.+(xj_out(i_out,j_out)-xj(-ng+1))/dxj
            posi_k=-ng+1.+(xk_out(i_out,j_out)-xk(-ng+1))/dxk

            posi_i_integer=floor(posi_i)
            posi_j_integer=floor(posi_j)
            posi_k_integer=floor(posi_k)

            w_i=posi_i-posi_i_integer
            w_j=posi_j-posi_j_integer
            w_k=posi_k-posi_k_integer

            fout(:,i_out,j_out)=&
                f(posi_i_integer  ,posi_j_integer  ,posi_k_integer  ,:)*(1.-w_i)*(1.-w_j)*(1.-w_k)+&
                f(posi_i_integer+1,posi_j_integer  ,posi_k_integer  ,:)*(   w_i)*(1.-w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer+1,posi_k_integer  ,:)*(1.-w_i)*(   w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer  ,posi_k_integer+1,:)*(1.-w_i)*(1.-w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer+1,posi_k_integer  ,:)*(   w_i)*(   w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer+1,posi_k_integer+1,:)*(1.-w_i)*(   w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer  ,posi_k_integer+1,:)*(   w_i)*(1.-w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer+1,posi_k_integer+1,:)*(   w_i)*(   w_j)*(   w_k)

        end do; end do

    end subroutine ModMath_1D3D_interpolate_1D2D
end module ModMath