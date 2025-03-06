module ModSpherical

    use ModDerivative

    contains

    function ModSpherical_cross(nr,nt,np,ng,A,B) result(C)
        implicit none
        integer,intent(in)          ::  nr,nt,np,ng
        real,intent(in)             ::  A(-ng+1:nr+ng,-ng+1:nt+ng,-ng+1:np+ng,3),&
                                        B(-ng+1:nr+ng,-ng+1:nt+ng,-ng+1:np+ng,3)
        real                        ::  C(-ng+1:nr+ng,-ng+1:nt+ng,-ng+1:np+ng,3)

        C(:,:,:,1)=A(:,:,:,2)*B(:,:,:,3)-A(:,:,:,3)*B(:,:,:,2)
        C(:,:,:,2)=A(:,:,:,3)*B(:,:,:,1)-A(:,:,:,1)*B(:,:,:,3)
        C(:,:,:,3)=A(:,:,:,1)*B(:,:,:,2)-A(:,:,:,2)*B(:,:,:,1)
    end function ModSpherical_cross

    function ModSpherical_curl(nr,nt,np,ng,r,t,dr,dt,dp,A) result(curl_A)
        implicit none
        integer,intent(in)          ::  nr,nt,np,ng
        real,intent(in)             ::  r(-ng+1:ng+nr),t(-ng+1:ng+nt)
        real,intent(in)             ::  dr,dt,dp
        real,intent(in)             ::  A(-ng+1:ng+nr,-ng+1:ng+nt,-ng+1:ng+np,1:3)
        real                        ::  At_r(-ng+1:ng+nr,-ng+1:ng+nt,-ng+1:ng+np)
        real                        ::  Ap_r(-ng+1:ng+nr,-ng+1:ng+nt,-ng+1:ng+np)
        real                        ::  Ar_sint(-ng+1:ng+nr,-ng+1:ng+nt,-ng+1:ng+np)
        real                        ::  Ap_sint(-ng+1:ng+nr,-ng+1:ng+nt,-ng+1:ng+np)
        real                        ::  sint(-ng+1:ng+nt),r_inverse(1:nr),sint_inverse(1:nt)
        real                        ::  curl_A(1:nr,1:nt,1:np,1:3)
        integer                     ::  ir,it,ip,idirection

        ! Preparations
        sint=sin(t)
        r_inverse=1./r(1:nr)
        sint_inverse=1./sint(1:nt)

        ! Get scaled Ap_sint, Ap_r, At_r
        do ip=-ng+1,ng+np; do it=-ng+1,ng+nt
            At_r(:,it,ip)=A(:,it,ip,2)*r
            Ap_r(:,it,ip)=A(:,it,ip,3)*r
        end do; end do

        do ip=-ng+1,ng+np; do ir=-ng+1,ng+nr
            Ar_sint(ir,:,ip)=A(ir,:,ip,1)/sint
            Ap_sint(ir,:,ip)=A(ir,:,ip,3)*sint
        end do; end do

        ! Get the derivatives
        ! r component
        curl_A(:,:,:,1)                                                 &
            =   ModDerivative_1st_O4_3D(Ap_sint,nr,nt,np,ng,dr,dt,dp,2) &
            -   ModDerivative_1st_O4_3D(A(:,:,:,2),nr,nt,np,ng,dr,dt,dp,3)
        
        ! th component
        curl_A(:,:,:,2)                                                 &
            =   ModDerivative_1st_O4_3D(Ar_sint,nr,nt,np,ng,dr,dt,dp,3) &
            -   ModDerivative_1st_O4_3D(Ap_r,nr,nt,np,ng,dr,dt,dp,1)
        
        ! ph component
        curl_A(:,:,:,3)                                                 &
            =   ModDerivative_1st_O4_3D(At_r,nr,nt,np,ng,dr,dt,dp,1)    &
            -   ModDerivative_1st_O4_3D(A(:,:,:,1),nr,nt,np,ng,dr,dt,dp,2)

        ! Rescale
        do idirection=1,3; do ip=1,np; do it=1,nt
            curl_A(:,it,ip,idirection)=curl_A(:,it,ip,idirection)*r_inverse
        end do; end do; end do

        do ip=1,np; do ir=1,nr
            curl_A(ir,:,ip,1)=curl_A(ir,:,ip,1)*sint_inverse
        end do; end do
    end function ModSpherical_curl

    function ModSpherical_div(nr,nt,np,ng,r,t,dr,dt,dp,A) result(div_A)
        implicit none
        integer,intent(in)          ::  nr,nt,np,ng
        real,intent(in)             ::  r(-ng+1:ng+nr),t(-ng+1:ng+nt)
        real,intent(in)             ::  dr,dt,dp
        real,intent(in)             ::  A(-ng+1:ng+nr,-ng+1:ng+nt,-ng+1:ng+np,1:3)
        real                        ::  A_scaled(-ng+1:ng+nr,-ng+1:ng+nt,-ng+1:ng+np,1:2)
        real                        ::  r2(-ng+1:ng+nr),sint(-ng+1:ng+nt),&
                                        r2_inverse(1:nr),r_inverse(1:nr),sint_inverse(1:nt)
        real                        ::  div_A(1:nr,1:nt,1:np),div_A23(1:nr,1:nt,1:np)
        integer                     ::  ir,it,ip

        ! Preparations
        r2=r**2
        sint=sin(t)
        r2_inverse=1./r2(1:nr)
        r_inverse=1./r(1:nr)
        sint_inverse=1./sint(1:nt)

        ! Get the scaled f
        do ip=-ng+1,ng+np; do it=-ng+1,ng+nt
            A_scaled(:,it,ip,1)=A(:,it,ip,1)*r2
        end do; end do

        do ip=-ng+1,ng+np; do ir=-ng+1,ng+nr
            A_scaled(ir,:,ip,2)=A(ir,:,ip,2)*sint
        end do; end do

        ! r term
        div_A=ModDerivative_1st_O4_3D(A_scaled(:,:,:,1),nr,nt,np,ng,dr,dt,dp,1)
        do ip=1,np; do it=1,nt
            div_A(:,it,ip)=div_A(:,it,ip)*r2_inverse
        end do; end do

        ! th & ph terms
        div_A23=ModDerivative_1st_O4_3D(A_scaled(:,:,:,2),nr,nt,np,ng,dr,dt,dp,2)
        div_A23=div_A23+ModDerivative_1st_O4_3D(A(:,:,:,3),nr,nt,np,ng,dr,dt,dp,3)
        do ip=1,np; do ir=1,nr
            div_A23(ir,:,ip)=div_A23(ir,:,ip)*sint_inverse
        end do; end do
        do ip=1,np; do it=1,nt
            div_A23(:,it,ip)=div_A23(:,it,ip)*r_inverse
        end do; end do

        div_A=div_A+div_A23
    end function ModSpherical_div

    function ModSpherical_Grad_f(nr,nt,np,ng,r,t,dr,dt,dp,f) result(Grad_f)
        implicit none
        integer,intent(in)          ::  nr,nt,np,ng
        real,intent(in)             ::  r(-ng+1:ng+nr),&
                                        t(-ng+1:ng+nt)
        real,intent(in)             ::  dr,dt,dp
        real,intent(in)             ::  f(-ng+1:ng+nr,-ng+1:ng+nt,-ng+1:ng+np)
        real                        ::  Grad_f(1:nr,1:nt,1:np,1:3)
        real                        ::  r_inverse(1:nr),sint_inverse(1:nt)
        integer                     ::  ir,it,ip

        ! Preparations
        r_inverse=1./r(1:nr)
        sint_inverse=1./sin(t(1:nt))

        ! r term
        Grad_f(:,:,:,1)=ModDerivative_1st_O4_3D(f,nr,nt,np,ng,dr,dt,dp,1)

        ! th term
        Grad_f(:,:,:,2)=ModDerivative_1st_O4_3D(f,nr,nt,np,ng,dr,dt,dp,2)
        do ip=1,np; do it=1,nt
            Grad_f(:,it,ip,2)=Grad_f(:,it,ip,2)*r_inverse
        end do; end do

        ! ph term
        Grad_f(:,:,:,3)=ModDerivative_1st_O4_3D(f,nr,nt,np,ng,dr,dt,dp,3)
        do ip=1,np; do it=1,nt
            Grad_f(:,it,ip,3)=Grad_f(:,it,ip,3)*r_inverse
        end do; end do
        do ip=1,np; do ir=1,nr
            Grad_f(ir,:,ip,3)=Grad_f(ir,:,ip,3)*sint_inverse
        end do; end do
    end function ModSpherical_Grad_f

    function ModSpherical_A_dot_Grad_f(nr,nt,np,ng,r,t,dr,dt,dp,A,f) result(A_dot_Grad_f)
        implicit none
        integer,intent(in)          ::  nr,nt,np,ng
        real,intent(in)             ::  r(-ng+1:ng+nr),&
                                        t(-ng+1:ng+nt)
        real,intent(in)             ::  dr,dt,dp
        real,intent(in)             ::  A(-ng+1:ng+nr,-ng+1:ng+nt,-ng+1:ng+np,1:3)
        real,intent(in)             ::  f(-ng+1:ng+nr,-ng+1:ng+nt,-ng+1:ng+np)
        real                        ::  A_dot_Grad_f(1:nr,1:nt,1:np),tmp(1:nr,1:nt,1:np)
        real                        ::  r_inverse(1:nr),sint_inverse(1:nt)
        integer                     ::  ir,it,ip

        ! Preparations
        r_inverse=1./r(1:nr)
        sint_inverse=1./sin(t(1:nt))
        ! r term
        A_dot_Grad_f=A(1:nr,1:nt,1:np,1)*ModDerivative_1st_O4_3D(f,nr,nt,np,ng,dr,dt,dp,1)

        ! th term
        tmp=A(1:nr,1:nt,1:np,2)*ModDerivative_1st_O4_3D(f,nr,nt,np,ng,dr,dt,dp,2)
        do ip=1,np; do it=1,nt
            A_dot_Grad_f(:,it,ip)=A_dot_Grad_f(:,it,ip)+tmp(:,it,ip)*r_inverse
        end do; end do

        ! ph term
        tmp=A(1:nr,1:nt,1:np,3)*ModDerivative_1st_O4_3D(f,nr,nt,np,ng,dr,dt,dp,3)
        do ip=1,np; do it=1,nt
            tmp(:,it,ip)=tmp(:,it,ip)*r_inverse
        end do; end do
        do ip=1,np; do ir=1,nr
            A_dot_Grad_f(ir,:,ip)=A_dot_Grad_f(ir,:,ip)+sint_inverse*tmp(ir,:,ip)
        end do; end do
    end function ModSpherical_A_dot_Grad_f

    function ModSpherical_A_dot_nabla_B(nr,nt,np,ng,r,t,dr,dt,dp,A,B) result(A_dot_nabla_B)
        implicit none
        integer,intent(in)          ::  nr,nt,np,ng
        real,intent(in)             ::  r(-ng+1:ng+nr),&
                                        t(-ng+1:ng+nt)
        real,intent(in)             ::  dr,dt,dp
        real,intent(in)             ::  A(-ng+1:ng+nr,-ng+1:ng+nt,-ng+1:ng+np,1:3)
        real,intent(in)             ::  B(-ng+1:ng+nr,-ng+1:ng+nt,-ng+1:ng+np,1:3)
        real                        ::  A_dot_nabla_B(1:nr,1:nt,1:np,1:3)
        real                        ::  tmp(1:nr,1:nt,1:np)
        real                        ::  cotant(1:nt),r_inverse(1:nr)
        integer                     ::  direction,ir,it,ip

        ! Preparations
        r_inverse=1./r(1:nr)
        cotant=cos(t(1:nt))/sin(t(1:nt))
        ! First three terms for direction i is simply A * nabla B_{i}
        do direction=1,3
            A_dot_nabla_B(:,:,:,direction)=&
                ModSpherical_A_dot_Grad_f(nr,nt,np,ng,r,t,dr,dt,dp,A,B(:,:,:,direction))
        end do

        ! However, since \hat{e_{i}} changes in spherical coordinate,
        ! There are several crosstalk terms due to nabla \hat{e_{i}}.
        ! First do the additional terms for r:

        tmp=-A(1:nr,1:nt,1:np,2)*B(1:nr,1:nt,1:np,2)-A(1:nr,1:nt,1:np,3)*B(1:nr,1:nt,1:np,3)
        do ip=1,np; do it=1,nt
            A_dot_nabla_B(:,it,ip,1)=A_dot_nabla_B(:,it,ip,1)+tmp(:,it,ip)*r_inverse
        end do; end do

        ! Then for th:

        tmp=-A(1:nr,1:nt,1:np,3)*B(1:nr,1:nt,1:np,3)
        do ip=1,np; do ir=1,nr
            tmp(ir,:,ip)=tmp(ir,:,ip)*cotant
        end do; end do
        tmp=tmp+A(1:nr,1:nt,1:np,2)*B(1:nr,1:nt,1:np,1)
        do ip=1,np; do it=1,nt
            A_dot_nabla_B(:,it,ip,2)=A_dot_nabla_B(:,it,ip,2)+tmp(:,it,ip)*r_inverse
        end do; end do

        ! Then for ph:

        tmp=A(1:nr,1:nt,1:np,3)*B(1:nr,1:nt,1:np,2)
        do ip=1,np; do ir=1,nr
            tmp(ir,:,ip)=tmp(ir,:,ip)*cotant
        end do; end do
        tmp=tmp+A(1:nr,1:nt,1:np,3)*B(1:nr,1:nt,1:np,1)
        do ip=1,np; do it=1,nt
            A_dot_nabla_B(:,it,ip,3)=A_dot_nabla_B(:,it,ip,3)+tmp(:,it,ip)*r_inverse
        end do; end do
    end function ModSpherical_A_dot_nabla_B


end module ModSpherical