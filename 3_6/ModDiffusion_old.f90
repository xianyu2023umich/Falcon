module ModDiffusion

    use ModParameters,      only:   iEquation
    use ModBlock,           only:   BlockType
    use ModParameters,      only:   ni,nj,nk,ng,nvar
    use ModLinReconstruct,  only:   ModLinReconstruct_minmod

    contains

    ! Definitions: l & r: left and right values at the face.
    ! plus & minus: cell center values for left and right cells of the face.

    subroutine ModDiffusion_Artificial_1(Block1,h)
        implicit none
        type(BlockType),target          ::  Block1
        integer,intent(in)              ::  h

        integer                         ::  i,j,ivar,nvar_here,direction1
        real(8),allocatable             ::  d_primitive(:,:,:,:)
        real(8),allocatable             ::  flux(:,:,:,:),&
                                            phi(:,:,:,:),&
                                            u_r_minus_u_l(:,:,:,:),&
                                            u_plus_minus_u_minus(:,:,:,:),&
                                            scales(:),&
                                            epsilon_denominator(:)

        ! Determine nvar_here. Allocate d_primitive based on nvar_here.
        nvar_here = nvar
        if (iEquation==1 .and. .not. Block1%if_involve_B) nvar_here = nvar-4
        allocate(d_primitive(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,nvar_here))

        ! Get the scales of the vars, which will be used
        ! to set the epsilon in the denominator when calculating phi.
        allocate(scales(nvar_here),&
                 epsilon_denominator(nvar_here))

        ! The scale should be no smaller than 1.e-30.
        do ivar=1,nvar_here
            scales(ivar)=maxval(abs(Block1%primitive(:,:,:,ivar)))
            epsilon_denominator(ivar)=1.0d2*epsilon(1.0d0)*max(1.e-30,scales(ivar))
        end do

        ! Get the diffusion flux and update EQN_update_R_IV for every direction.
        do direction1=1,3

            ! use minmod to find \Delta u
            call ModLinReconstruct_minmod(nvar_here,ni,nj,nk,ng,direction1,Block1%primitive,d_primitive)

            ! get phi & flux
            select case(direction1)
            case(1)
                ! Allocate all needed arrays
                allocate(   flux(ni+1,nj,nk,nvar_here),         &
                            phi(ni+1,nj,nk,nvar_here),          &
                            u_r_minus_u_l(ni+1,nj,nk,nvar_here),&
                            u_plus_minus_u_minus(ni+1,nj,nk,nvar_here) )

                ! Get u_plus-u_minus & u_r-u_l for the i+1/2 face.

                u_plus_minus_u_minus=&
                    Block1%primitive(1:ni+1,1:nj,1:nk,1:nvar_here)  &
                    -Block1%primitive(0:ni,1:nj,1:nk,1:nvar_here)
                
                u_r_minus_u_l=u_plus_minus_u_minus-0.5*(        &
                    d_primitive(1:ni+1,1:nj,1:nk,1:nvar_here)   &
                   +d_primitive(0:ni,1:nj,1:nk,1:nvar_here))

                
                ! Correct too small u_plus_minus_u_minus to avoid too large phi.
                do ivar=1,nvar_here
                    u_plus_minus_u_minus(:,:,:,ivar)=sign(&
                        max(abs(u_plus_minus_u_minus(:,:,:,ivar)),epsilon_denominator(ivar)),&
                        u_plus_minus_u_minus(:,:,:,ivar))
                end do

                ! Get \Phi_{h}
                ! phi = max[0, 1+((u_r-u_l)/(u_i+1-u_i)-1)],
                ! for the i+1/2 face.
                phi=max(0.0,1.0+h*(-1.0+u_r_minus_u_l/u_plus_minus_u_minus))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux=-0.5*phi*u_r_minus_u_l
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,nvar_here
                    flux(:,:,:,ivar)=flux(:,:,:,ivar)*&
                    (Block1%v_wave_III(1:ni+1,:,:)+Block1%v_wave_III(0:ni,:,:))*0.5
                end do

                ! update EQN_update_R
                do ivar=1,nvar_here
                    Block1%EQN_update_R_IV(:,:,:,ivar)=Block1%EQN_update_R_IV(:,:,:,ivar)+&
                        (flux(1:ni,:,:,ivar)-flux(2:ni+1,:,:,ivar))/Block1%dxi
                end do
            case(2)
                allocate(   flux(ni,nj+1,nk,nvar_here),&
                            phi(ni,nj+1,nk,nvar_here),          &
                            u_r_minus_u_l(ni,nj+1,nk,nvar_here),&
                            u_plus_minus_u_minus(ni,nj+1,nk,nvar_here) )

                ! Get u_plus-u_minus & u_r-u_l for the i+1/2 face.

                u_plus_minus_u_minus=&
                    Block1%primitive(1:ni,1:nj+1,1:nk,1:nvar_here)  &
                    -Block1%primitive(1:ni,0:nj,1:nk,1:nvar_here)
                
                u_r_minus_u_l=u_plus_minus_u_minus-0.5*(        &
                    d_primitive(1:ni,1:nj+1,1:nk,1:nvar_here)   &
                   +d_primitive(1:ni,0:nj,1:nk,1:nvar_here))

                
                ! Correct too small u_plus_minus_u_minus to avoid too large phi.
                do ivar=1,nvar_here
                    u_plus_minus_u_minus(:,:,:,ivar)=sign(&
                        max(abs(u_plus_minus_u_minus(:,:,:,ivar)),epsilon_denominator(ivar)),&
                        u_plus_minus_u_minus(:,:,:,ivar))
                end do

                ! Get \Phi_{h}
                ! phi = max[0, 1+((u_r-u_l)/(u_i+1-u_i)-1)],
                ! for the i+1/2 face.
                phi=max(0.0,1.0+h*(-1.0+u_r_minus_u_l/u_plus_minus_u_minus))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux=-0.5*phi*u_r_minus_u_l
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,nvar_here
                    flux(:,:,:,ivar)=flux(:,:,:,ivar)*&
                    (Block1%v_wave_III(:,1:nj+1,:)+Block1%v_wave_III(:,0:nj,:))*0.5
                end do

                do ivar=1,nvar_here
                    do i=1,ni
                        Block1%EQN_update_R_IV(i,:,:,ivar)=Block1%EQN_update_R_IV(i,:,:,ivar)+&
                            (flux(i,1:nj,:,ivar)-flux(i,2:nj+1,:,ivar))/(Block1%dxj*Block1%xi_I(i))
                    end do
                end do
            case(3)
                allocate(   flux(ni,nj,nk+1,nvar_here),&
                            phi(ni,nj,nk+1,nvar_here),          &
                            u_r_minus_u_l(ni,nj,nk+1,nvar_here),&
                            u_plus_minus_u_minus(ni,nj,nk+1,nvar_here) )

                ! Get u_plus-u_minus & u_r-u_l for the i+1/2 face.

                u_plus_minus_u_minus=&
                    Block1%primitive(1:ni,1:nj,1:nk+1,1:nvar_here)  &
                    -Block1%primitive(1:ni,1:nj,0:nk,1:nvar_here)
                
                u_r_minus_u_l=u_plus_minus_u_minus-0.5*(        &
                    d_primitive(1:ni,1:nj,1:nk+1,1:nvar_here)   &
                   +d_primitive(1:ni,1:nj,0:nk,1:nvar_here))

                
                ! Correct too small u_plus_minus_u_minus to avoid too large phi.
                do ivar=1,nvar_here
                    u_plus_minus_u_minus(:,:,:,ivar)=sign(&
                        max(abs(u_plus_minus_u_minus(:,:,:,ivar)),epsilon_denominator(ivar)),&
                        u_plus_minus_u_minus(:,:,:,ivar))
                end do

                ! Get \Phi_{h}
                ! phi = max[0, 1+((u_r-u_l)/(u_i+1-u_i)-1)],
                ! for the i+1/2 face.
                phi=max(0.0,1.0+h*(-1.0+u_r_minus_u_l/u_plus_minus_u_minus))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux=-0.5*phi*u_r_minus_u_l
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,nvar_here
                    flux(:,:,:,ivar)=flux(:,:,:,ivar)*&
                    (Block1%v_wave_III(:,:,1:nk+1)+Block1%v_wave_III(:,:,0:nk))*0.5
                end do

                do ivar=1,nvar_here
                    do j=1,nj
                        do i=1,ni
                            Block1%EQN_update_R_IV(i,j,:,ivar)=Block1%EQN_update_R_IV(i,j,:,ivar)+&
                                (flux(i,j,1:nk,ivar)-flux(i,j,2:nk+1,ivar))/(Block1%dxk*Block1%xi_I(i)*sin(Block1%xj_I(j)))
                        end do
                    end do
                end do
            end select

            ! Deallocate arrays because they will be reallocated.
            deallocate(flux,phi,u_r_minus_u_l,u_plus_minus_u_minus)
        end do
    end subroutine ModDiffusion_Artificial_1

end module ModDiffusion