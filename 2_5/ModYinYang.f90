module ModYinYang

    contains

    ! Functions to convert Yin rtp to Yang rtp
    ! or change Yang rtp to Yin rtp

    function ModYinYang_CoordConv_0D(rtp_in) result(rtp_out)
        implicit none
        real,intent(in)     ::  rtp_in(3)       ! input
        real                ::  rtp_out(3)      ! output

        rtp_out(1)=rtp_in(1)
        rtp_out(2)=acos((sin(rtp_in(2))*sin(rtp_in(3))))
        rtp_out(3)=sign(acos(-sin(rtp_in(2))*cos(rtp_in(3))/&
            sqrt(sin(rtp_in(2))**2*cos(rtp_in(3))**2+cos(rtp_in(2))**2)),cos(rtp_in(2)))
    end function ModYinYang_CoordConv_0D

    function ModYinYang_CoordConv_1D(rtp_in,n) result(rtp_out)
        implicit none
        real,intent(in)     ::  rtp_in(n,3)     ! input
        integer,intent(in)  ::  n               ! num of points
        real                ::  rtp_out(n,3)    ! output

        rtp_out(:,1)=rtp_in(:,1)
        rtp_out(:,2)=acos((sin(rtp_in(:,2))*sin(rtp_in(:,3))))
        rtp_out(:,3)=sign(acos(-sin(rtp_in(:,2))*cos(rtp_in(:,3))/&
            sqrt(sin(rtp_in(:,2))**2*cos(rtp_in(:,3))**2+cos(rtp_in(:,2))**2)),cos(rtp_in(:,2)))
    end function ModYinYang_CoordConv_1D


    ! Functions to convert vectors

    function ModYinYang_VecConv_0D(rtp_in,vec_in) result(vec_out)
        implicit none
        real,intent(in)     ::  rtp_in(3)       ! input coord
        real,intent(in)     ::  vec_in(3)       ! input vector
        real                ::  rtp_out(3)      ! output coord
        real                ::  vec_out(3)      ! output vector

        rtp_out=ModYinYang_CoordConv_0D(rtp_in)

        vec_out(1)=vec_in(1)
        vec_out(2)=-sin(rtp_in(3))*sin(rtp_out(3))*vec_in(2)-cos(rtp_in(3))/sin(rtp_out(2))*vec_in(3)
        vec_out(3)=-sin(rtp_in(3))*sin(rtp_out(3))*vec_in(3)+cos(rtp_in(3))/sin(rtp_out(2))*vec_in(2)
    end function ModYinYang_VecConv_0D

    function ModYinYang_VecConv_1D(rtp_in,vec_in,n) result(vec_out)
        implicit none
        real,intent(in)     ::  rtp_in(n,3)   ! input coord
        real,intent(in)     ::  vec_in(n,3)   ! input vector
        integer,intent(in)  ::  n
        real                ::  rtp_out(n,3)  ! output coord
        real                ::  vec_out(n,3)  ! output vector

        rtp_out=ModYinYang_CoordConv_1D(rtp_in,n)

        vec_out(:,1)=vec_in(:,1)
        vec_out(:,2)=-sin(rtp_in(:,3))*sin(rtp_out(:,3))*vec_in(:,2)-cos(rtp_in(:,3))/sin(rtp_out(:,2))*vec_in(:,3)
        vec_out(:,3)=-sin(rtp_in(:,3))*sin(rtp_out(:,3))*vec_in(:,3)+cos(rtp_in(:,3))/sin(rtp_out(:,2))*vec_in(:,2)
    end function ModYinYang_VecConv_1D


    ! If I have an rtp_range in one branch, but I want to 
    ! get the minmum rectangular that covers this rtp_range
    ! IN THE OTHER BRANCH, then use this function.

    function ModYinYang_GetOtherRange(rtp_range_in) result(rtp_range_out)
        implicit none
        real                ::  rtp_range_in(3,2)
        real                ::  rtp_corner_out(3),tp_corners_out(2,2,2)
        real                ::  rtp_range_out(3,2)
        integer             ::  it,ip

        rtp_range_out(1,:)=rtp_range_in(1,:)

        do it=1,2; do ip=1,2
            rtp_corner_out=ModYinYang_CoordConv_0D(&
                [rtp_range_in(1,1),rtp_range_in(2,it),rtp_range_in(3,ip)])
            tp_corners_out(1:2,it,ip)=rtp_corner_out(2:3)
        end do; end do

        rtp_range_out(2,1)=minval(tp_corners_out(1,:,:))
        rtp_range_out(2,2)=maxval(tp_corners_out(1,:,:))
        rtp_range_out(3,1)=minval(tp_corners_out(2,:,:))
        rtp_range_out(3,2)=maxval(tp_corners_out(2,:,:))
    end function ModYinYang_GetOtherRange
end module ModYinYang