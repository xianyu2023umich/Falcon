Module ModDerivative

    contains

    function ModDerivative_1st_O4_3D(f,ni,nj,nk,ng,dxi,dxj,dxk,direction) result(df_dx)

        implicit none
        integer,intent(in) :: ni,nj,nk,ng,direction
        real,intent(in) :: f(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        real,intent(in) :: dxi,dxj,dxk

        real :: df_dx(1:ni,1:nj,1:nk)
        
        select case(direction)
        case(1)
            df_dx=(1./12./dxi)*( &
                -     f(( 3):(ni+2),1:nj,1:nk) &
                + 8.* f(( 2):(ni+1),1:nj,1:nk) &
                - 8.* f(( 0):(ni-1),1:nj,1:nk) &
                +     f((-1):(ni-2),1:nj,1:nk))
        case(2)
            df_dx=(1./12./dxj)*( &
                -     f(1:ni,( 3):(nj+2),1:nk) &
                + 8.* f(1:ni,( 2):(nj+1),1:nk) &
                - 8.* f(1:ni,( 0):(nj-1),1:nk) &
                +     f(1:ni,(-1):(nj-2),1:nk))
        case(3)
            df_dx=(1./12./dxk)*( &
                -     f(1:ni,1:nj,( 3):(nk+2)) &
                + 8.* f(1:ni,1:nj,( 2):(nk+1)) &
                - 8.* f(1:ni,1:nj,( 0):(nk-1)) &
                +     f(1:ni,1:nj,(-1):(nk-2)))
        end select
    end function ModDerivative_1st_O4_3D

    function ModDerivative_1st_O2(f,ni,nj,nk,ng,dxi,dxj,dxk,direction,geometry) result(df_dx)
        implicit none
        integer,intent(in) :: ni,nj,nk,ng,direction
        real,intent(in) :: f(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        real,intent(in) :: dxi,dxj,dxk
        character(len=*),intent(in) :: geometry

        real :: df_dx(1:ni,1:nj,1:nk)

        select case(geometry)
        case('cartesian')
            select case(direction)
            case(1)
                df_dx=(1./2./dxi)*( &
                    + f((2):(ni+1),1:nj,1:nk) &
                    - f((0):(ni-1),1:nj,1:nk))
            case(2)
                df_dx=(1./2./dxj)*( &
                    + f(1:ni,(2):(nj+1),1:nk) &
                    - f(1:ni,(0):(nj-1),1:nk))
            case(3)
                df_dx=(1./2./dxk)*( &
                    + f(1:ni,1:nj,(2):(nk+1)) &
                    - f(1:ni,1:nj,(0):(nk-1)))
            end select
        end select
    end function ModDerivative_1st_O2

    function ModDerivative_pxpy_O4(f,ni,nj,nk,ng,dxi,dxj,dxk,direction1,direction2,geometry) result(d2f_dx_dy)
        implicit none
        integer,intent(in) :: ni,nj,nk,ng,direction1,direction2
        real,intent(in) :: f(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        real,intent(in) :: dxi,dxj,dxk
        character(len=*),intent(in) :: geometry

        real :: df_dx(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        real :: d2f_dx_dy(1:ni,1:nj,1:nk)

        select case(geometry)
        case('cartesian')
            select case(direction1)
            case(1)
                df_dx(1:ni,:,:)=(1./12./dxi)*( &
                    -     f(( 3):(ni+2),:,:) &
                    + 8.* f(( 2):(ni+1),:,:) &
                    - 8.* f(( 0):(ni-1),:,:) &
                    +     f((-1):(ni-2),:,:))
            case(2)
                df_dx(:,1:nj,:)=(1./12./dxj)*( &
                    -     f(:,( 3):(nj+2),:) &
                    + 8.* f(:,( 2):(nj+1),:) &
                    - 8.* f(:,( 0):(nj-1),:) &
                    +     f(:,(-1):(nj-2),:))
            case(3)
                df_dx(:,:,1:nk)=(1./12./dxk)*( &
                    -     f(:,:,( 3):(nk+2)) &
                    + 8.* f(:,:,( 2):(nk+1)) &
                    - 8.* f(:,:,( 0):(nk-1)) &
                    +     f(:,:,(-1):(nk-2)))
            end select
            select case(direction2)
            case(1)
                d2f_dx_dy=(1./12./dxi)*( &
                    -     df_dx(( 3):(ni+2),1:nj,1:nk) &
                    + 8.* df_dx(( 2):(ni+1),1:nj,1:nk) &
                    - 8.* df_dx(( 0):(ni-1),1:nj,1:nk) &
                    +     df_dx((-1):(ni-2),1:nj,1:nk))
            case(2)
                d2f_dx_dy=(1./12./dxj)*( &
                    -     df_dx(1:ni,( 3):(nj+2),1:nk) &
                    + 8.* df_dx(1:ni,( 2):(nj+1),1:nk) &
                    - 8.* df_dx(1:ni,( 0):(nj-1),1:nk) &
                    +     df_dx(1:ni,(-1):(nj-2),1:nk))
            case(3)
                d2f_dx_dy=(1./12./dxk)*( &
                    -     df_dx(1:ni,1:nj,( 3):(nk+2)) &
                    + 8.* df_dx(1:ni,1:nj,( 2):(nk+1)) &
                    - 8.* df_dx(1:ni,1:nj,( 0):(nk-1)) &
                    +     df_dx(1:ni,1:nj,(-1):(nk-2)))
            end select
        end select

        
    end function ModDerivative_pxpy_O4

    function ModDerivative_2nd_O3(f,ni,nj,nk,ng,dxi,dxj,dxk,direction,geometry) result(d2f_dx2)
        implicit none
        integer,intent(in) :: ni,nj,nk,ng,direction
        real,intent(in) :: f(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        real,intent(in) :: dxi,dxj,dxk
        character(len=*),intent(in) :: geometry

        real :: d2f_dx2(1:ni,1:nj,1:nk)

        select case(geometry)
        case('cartesian')
            select case(direction)
            case(1)
                d2f_dx2=(1./12./dxi**2)*( &
                    -      f(( 3):(ni+2),1:nj,1:nk) &
                    + 16.* f(( 2):(ni+1),1:nj,1:nk) &
                    - 30.* f(( 1):(ni  ),1:nj,1:nk) &
                    + 16.* f(( 0):(ni-1),1:nj,1:nk) &
                    -      f((-1):(ni-2),1:nj,1:nk))
            case(2)
                d2f_dx2=(1./12./dxj**2)*( &
                    -      f(1:ni,( 3):(nj+2),1:nk) &
                    + 16.* f(1:ni,( 2):(nj+1),1:nk) &
                    - 30.* f(1:ni,( 1):(nj  ),1:nk) &
                    + 16.* f(1:ni,( 0):(nj-1),1:nk) &
                    -      f(1:ni,(-1):(nj-2),1:nk))
            case(3)
                d2f_dx2=(1./12./dxk**2)*( &
                    -      f(1:ni,1:nj,( 3):(nk+2)) &
                    + 16.* f(1:ni,1:nj,( 2):(nk+1)) &
                    - 30.* f(1:ni,1:nj,( 1):(nk  )) &
                    + 16.* f(1:ni,1:nj,( 0):(nk-1)) &
                    -      f(1:ni,1:nj,(-1):(nk-2)))
            end select
        end select
    end function ModDerivative_2nd_O3

end Module ModDerivative