module ModMHD

    use ModConst,      only:   mu0__CGS

    implicit none

    contains

    function ModMHD_AlfvenVelocityVector_3D(ni,nj,nk,ng,rho_III,B_IV) result(AlfvenVelocityVector_IV)
        implicit none
        integer,intent(in)          ::  ni,nj,nk,ng
        real(8),intent(in)          ::  rho_III(-ng+1:ng+ni,-ng+1:ng+nj,-ng+1:ng+nk)
        real(8),intent(in)          ::  B_IV(-ng+1:ng+ni,-ng+1:ng+nj,-ng+1:ng+nk,3)
        real(8)                     ::  AlfvenVelocityVector_IV(-ng+1:ng+ni,-ng+1:ng+nj,-ng+1:ng+nk,3)
        integer                     ::  ivar

        do ivar=1,3
            AlfvenVelocityVector_IV(:,:,:,ivar)=B_IV(:,:,:,ivar)/sqrt(mu0__CGS*rho_III)
        end do

    end function ModMHD_AlfvenVelocityVector_3D

end module ModMHD