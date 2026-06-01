module ModDiffusion

    use ModParameters,      only: iEquation
    use ModBlock,           only: BlockType
    use ModParameters,      only: ni,nj,nk,ng,nvar
    use ModLinReconstruct,  only: minmod_slope_x, minmod_slope_y, minmod_slope_z

    implicit none

    contains

    subroutine ModDiffusion_Artificial_1(Block1,h)
        implicit none

        type(BlockType), intent(inout), target :: Block1
        integer, intent(in) :: h

        integer :: i,j,k,ivar,nvar_here
        real(8) :: hh
        real(8), allocatable :: eps_den(:)
        real(8), allocatable :: flux_x(:,:,:,:)
        real(8), allocatable :: flux_y(:,:,:,:)
        real(8), allocatable :: flux_z(:,:,:,:)

        nvar_here = nvar
        if (iEquation == 1 .and. .not. Block1%if_involve_B) nvar_here = nvar - 4

        hh = real(h,8)

        allocate(eps_den(nvar_here))

        do ivar = 1,nvar_here
            eps_den(ivar) = 100.0d0 * epsilon(1.0d0) * &
                max(maxval(abs(Block1%primitive(:,:,:,ivar))), tiny(1.0d0))
        enddo

        ! x-direction flux: interfaces between i-1 and i, indexed 1:ni+1
        allocate(flux_x(1:ni+1,1:nj,1:nk,1:nvar_here))
        call get_flux_x(Block1,hh,eps_den,flux_x)

        do ivar = 1,nvar_here
            do k = 1,nk
                do j = 1,nj
                    do i = 1,ni
                        Block1%EQN_update_R_IV(i,j,k,ivar) = Block1%EQN_update_R_IV(i,j,k,ivar) &
                            + (flux_x(i,j,k,ivar) - flux_x(i+1,j,k,ivar)) / Block1%dxi
                    enddo
                enddo
            enddo
        enddo

        deallocate(flux_x)

        ! y-direction flux: interfaces between j-1 and j, indexed 1:nj+1
        allocate(flux_y(1:ni,1:nj+1,1:nk,1:nvar_here))
        call get_flux_y(Block1,hh,eps_den,flux_y)

        do ivar = 1,nvar_here
            do k = 1,nk
                do j = 1,nj
                    do i = 1,ni
                        Block1%EQN_update_R_IV(i,j,k,ivar) = Block1%EQN_update_R_IV(i,j,k,ivar) &
                            + (flux_y(i,j,k,ivar) - flux_y(i,j+1,k,ivar)) &
                            / (Block1%dxj * Block1%xi_I(i))
                    enddo
                enddo
            enddo
        enddo

        deallocate(flux_y)

        ! z-direction flux: interfaces between k-1 and k, indexed 1:nk+1
        allocate(flux_z(1:ni,1:nj,1:nk+1,1:nvar_here))
        call get_flux_z(Block1,hh,eps_den,flux_z)

        do ivar = 1,nvar_here
            do k = 1,nk
                do j = 1,nj
                    do i = 1,ni
                        Block1%EQN_update_R_IV(i,j,k,ivar) = Block1%EQN_update_R_IV(i,j,k,ivar) &
                            + (flux_z(i,j,k,ivar) - flux_z(i,j,k+1,ivar)) &
                            / (Block1%dxk * Block1%xi_I(i) * sin(Block1%xj_I(j)))
                    enddo
                enddo
            enddo
        enddo

        deallocate(flux_z)
        deallocate(eps_den)

    end subroutine ModDiffusion_Artificial_1


    subroutine get_flux_x(Block1,hh,eps_den,flux_x)
        implicit none

        type(BlockType), intent(in) :: Block1
        real(8), intent(in) :: hh
        real(8), intent(in) :: eps_den(:)
        real(8), intent(out) :: flux_x(1:ni+1,1:nj,1:nk,1:size(eps_den))

        integer :: i,j,k,ivar
        integer :: iL,iR
        real(8) :: slope_l, slope_r
        real(8) :: ul, ur
        real(8) :: du, dface, r, phi, cface

        do ivar = 1,size(eps_den)
            do k = 1,nk
                do j = 1,nj
                    do i = 1,ni+1

                        iL = i - 1
                        iR = i

                        slope_l = minmod_slope_x(Block1%primitive,iL,j,k,ivar)
                        slope_r = minmod_slope_x(Block1%primitive,iR,j,k,ivar)

                        ul = Block1%primitive(iL,j,k,ivar) + 0.5d0*slope_l
                        ur = Block1%primitive(iR,j,k,ivar) - 0.5d0*slope_r

                        du    = Block1%primitive(iR,j,k,ivar) - Block1%primitive(iL,j,k,ivar)
                        dface = ur - ul

                        phi = 0.0d0
                        if (dface*du > 0.0d0 .and. abs(du) > eps_den(ivar)) then
                            r = dface / du
                            r = max(0.0d0, min(1.0d0, r))
                            phi = max(0.0d0, 1.0d0 + hh*(r - 1.0d0))
                        endif

                        cface = 0.5d0*(Block1%v_wave_III(iL,j,k) + Block1%v_wave_III(iR,j,k))

                        flux_x(i,j,k,ivar) = -0.5d0*cface*phi*dface

                    enddo
                enddo
            enddo
        enddo

    end subroutine get_flux_x


    subroutine get_flux_y(Block1,hh,eps_den,flux_y)
        implicit none

        type(BlockType), intent(in) :: Block1
        real(8), intent(in) :: hh
        real(8), intent(in) :: eps_den(:)
        real(8), intent(out) :: flux_y(1:ni,1:nj+1,1:nk,1:size(eps_den))

        integer :: i,j,k,ivar
        integer :: jL,jR
        real(8) :: slope_l, slope_r
        real(8) :: ul, ur
        real(8) :: du, dface, r, phi, cface

        do ivar = 1,size(eps_den)
            do k = 1,nk
                do j = 1,nj+1
                    do i = 1,ni

                        jL = j - 1
                        jR = j

                        slope_l = minmod_slope_y(Block1%primitive,i,jL,k,ivar)
                        slope_r = minmod_slope_y(Block1%primitive,i,jR,k,ivar)

                        ul = Block1%primitive(i,jL,k,ivar) + 0.5d0*slope_l
                        ur = Block1%primitive(i,jR,k,ivar) - 0.5d0*slope_r

                        du    = Block1%primitive(i,jR,k,ivar) - Block1%primitive(i,jL,k,ivar)
                        dface = ur - ul

                        phi = 0.0d0
                        if (dface*du > 0.0d0 .and. abs(du) > eps_den(ivar)) then
                            r = dface / du
                            r = max(0.0d0, min(1.0d0, r))
                            phi = max(0.0d0, 1.0d0 + hh*(r - 1.0d0))
                        endif

                        cface = 0.5d0*(Block1%v_wave_III(i,jL,k) + Block1%v_wave_III(i,jR,k))

                        flux_y(i,j,k,ivar) = -0.5d0*cface*phi*dface

                    enddo
                enddo
            enddo
        enddo

    end subroutine get_flux_y


    subroutine get_flux_z(Block1,hh,eps_den,flux_z)
        implicit none

        type(BlockType), intent(in) :: Block1
        real(8), intent(in) :: hh
        real(8), intent(in) :: eps_den(:)
        real(8), intent(out) :: flux_z(1:ni,1:nj,1:nk+1,1:size(eps_den))

        integer :: i,j,k,ivar
        integer :: kL,kR
        real(8) :: slope_l, slope_r
        real(8) :: ul, ur
        real(8) :: du, dface, r, phi, cface

        do ivar = 1,size(eps_den)
            do k = 1,nk+1
                do j = 1,nj
                    do i = 1,ni

                        kL = k - 1
                        kR = k

                        slope_l = minmod_slope_z(Block1%primitive,i,j,kL,ivar)
                        slope_r = minmod_slope_z(Block1%primitive,i,j,kR,ivar)

                        ul = Block1%primitive(i,j,kL,ivar) + 0.5d0*slope_l
                        ur = Block1%primitive(i,j,kR,ivar) - 0.5d0*slope_r

                        du    = Block1%primitive(i,j,kR,ivar) - Block1%primitive(i,j,kL,ivar)
                        dface = ur - ul

                        phi = 0.0d0
                        if (dface*du > 0.0d0 .and. abs(du) > eps_den(ivar)) then
                            r = dface / du
                            r = max(0.0d0, min(1.0d0, r))
                            phi = max(0.0d0, 1.0d0 + hh*(r - 1.0d0))
                        endif

                        cface = 0.5d0*(Block1%v_wave_III(i,j,kL) + Block1%v_wave_III(i,j,kR))

                        flux_z(i,j,k,ivar) = -0.5d0*cface*phi*dface

                    enddo
                enddo
            enddo
        enddo

    end subroutine get_flux_z

end module ModDiffusion