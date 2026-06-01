module ModDiffusion

    use ModBlock,           only:   BlockType
    use ModParameters,      only:   ni,nj,nk,ng,nvar
    use ModLinReconstruct,  only:   ModLinReconstruct_minmod

    contains

    subroutine ModDiffusion_Artificial_1(Block1,h)
        implicit none
        type(BlockType),target          ::  Block1
        integer,intent(in)              ::  h

        integer                         ::  direction1
        real(8)                         ::  c(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        real(8)                         ::  d_primitive(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,nvar)
        real(8),allocatable             ::  flux(:,:,:,:),phi(:,:,:,:)
        integer                         ::  i,j,ivar

        ! Get the sound speed
        c=Block1%v_wave_III

        ! Get the primitive pointer
        Block1%primitive=>Block1%primitive_IV

        do direction1=1,3

            ! get total speed c
            !c=abs(Block1%primitive(:,:,:,direction1+1))+c_s

            ! use minmod to find \Delta u
            call ModLinReconstruct_minmod(nvar,ni,nj,nk,ng,direction1,Block1%primitive,d_primitive)

            ! get phi & flux
            select case(direction1)
            case(1)
                allocate(flux(ni+1,nj,nk,nvar),phi(ni+1,nj,nk,nvar))

                ! first get \Phi_{h}
                ! phi = max[0, 1+((u_r-u_l)/(u_i+1-u_i)-1)],
                ! for the i+1/2 face.
                phi(1:ni+1,:,:,:)=&
                    max(0.0,1.0+h*((Block1%primitive(1:ni+1,1:nj,1:nk,:)-d_primitive(1:ni+1,1:nj,1:nk,:)*0.5-&
                    Block1%primitive(0:ni,1:nj,1:nk,:)-d_primitive(0:ni,1:nj,1:nk,:)*0.5)/&
                    (Block1%primitive(1:ni+1,1:nj,1:nk,:)-Block1%primitive(0:ni,1:nj,1:nk,:))-1))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux(1:ni+1,:,:,:)=&
                    -0.5*phi*(Block1%primitive(1:ni+1,1:nj,1:nk,:)-d_primitive(1:ni+1,1:nj,1:nk,:)*0.5-&
                    Block1%primitive(0:ni,1:nj,1:nk,:)-d_primitive(0:ni,1:nj,1:nk,:)*0.5)
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,nvar
                    flux(:,:,:,ivar)=flux(:,:,:,ivar)*(c(1:ni+1,:,:)+c(0:ni,:,:))*0.5
                end do

                !do ivar=vr_,vp_
                !    flux(:,:,:,ivar)=flux(:,:,:,ivar)+&
                !        (Block1%primitive(1:ni+1,:,:,ivar)+Block1%primitive(0:ni,:,:,ivar))*&
                !        0.5*flux(:,:,:,rho1_)
                !end do

                ! update EQN_update_R
                do ivar=1,nvar
                    Block1%EQN_update_R_IV(:,:,:,ivar)=Block1%EQN_update_R_IV(:,:,:,ivar)+&
                        (flux(1:ni,:,:,ivar)-flux(2:ni+1,:,:,ivar))/Block1%dxi
                end do
            case(2)
                allocate(flux(ni,nj+1,nk,nvar),phi(ni,nj+1,nk,nvar))

                ! first get \Phi_{h}
                ! phi = max[0, 1+((u_r-u_l)/(u_i+1-u_i)-1)],
                ! for the i+1/2 face.
                phi(:,1:nj+1,:,:)=&
                    max(0.0,1.0+h*((Block1%primitive(1:ni,1:nj+1,1:nk,:)-d_primitive(1:ni,1:nj+1,1:nk,:)*0.5-&
                    Block1%primitive(1:ni,0:nj,1:nk,:)-d_primitive(1:ni,0:nj,1:nk,:)*0.5)/&
                    (Block1%primitive(1:ni,1:nj+1,1:nk,:)-Block1%primitive(1:ni,0:nj,1:nk,:))-1))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux(:,1:nj+1,:,:)=&
                    -0.5*phi*(Block1%primitive(1:ni,1:nj+1,1:nk,:)-d_primitive(1:ni,1:nj+1,1:nk,:)*0.5-&
                    Block1%primitive(1:ni,0:nj,1:nk,:)-d_primitive(1:ni,0:nj,1:nk,:)*0.5)
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,nvar
                    flux(:,:,:,ivar)=flux(:,:,:,ivar)*(c(:,1:nj+1,:)+c(:,0:nj,:))*0.5
                end do

                !do ivar=vr_,vp_
                !    flux(:,:,:,ivar)=flux(:,:,:,ivar)+&
                !        (Block1%primitive(:,1:nj+1,:,ivar)+Block1%primitive(:,0:nj,:,ivar))*&
                !        0.5*flux(:,:,:,rho1_)
                !end do

                do ivar=1,nvar
                    do i=1,ni
                        Block1%EQN_update_R_IV(i,:,:,ivar)=Block1%EQN_update_R_IV(i,:,:,ivar)+&
                            (flux(i,1:nj,:,ivar)-flux(i,2:nj+1,:,ivar))/(Block1%dxj*Block1%xi_I(i))
                    end do
                end do
            case(3)
                allocate(flux(ni,nj,nk+1,nvar),phi(ni,nj,nk+1,nvar))

                ! first get \Phi_{h}
                ! phi = max[0, 1+((u_r-u_l)/(u_i+1-u_i)-1)],
                ! for the i+1/2 face.
                phi(:,:,1:nk+1,:)=&
                    max(0.0,1.0+h*((Block1%primitive(1:ni,1:nj,1:nk+1,:)-d_primitive(1:ni,1:nj,1:nk+1,:)*0.5-&
                    Block1%primitive(1:ni,1:nj,0:nk,:)-d_primitive(1:ni,1:nj,0:nk,:)*0.5)/&
                    (Block1%primitive(1:ni,1:nj,1:nk+1,:)-Block1%primitive(1:ni,1:nj,0:nk,:))-1))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux(:,:,1:nk+1,:)=&
                    -0.5*phi*(Block1%primitive(1:ni,1:nj,1:nk+1,:)-d_primitive(1:ni,1:nj,1:nk+1,:)*0.5-&
                    Block1%primitive(1:ni,1:nj,0:nk,:)-d_primitive(1:ni,1:nj,0:nk,:)*0.5)
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,nvar
                    flux(:,:,:,ivar)=flux(:,:,:,ivar)*(c(:,:,1:nk+1)+c(:,:,0:nk))*0.5
                end do

                !do ivar=vr_,vp_
                !    flux(:,:,:,ivar)=flux(:,:,:,ivar)+&
                !        (Block1%primitive(:,:,1:nk+1,ivar)+Block1%primitive(:,:,0:nk,ivar))*&
                !        0.5*flux(:,:,:,rho1_)
                !end do

                do ivar=1,nvar
                    do j=1,nj
                        do i=1,ni
                            Block1%EQN_update_R_IV(i,j,:,ivar)=Block1%EQN_update_R_IV(i,j,:,ivar)+&
                                (flux(i,j,1:nk,ivar)-flux(i,j,2:nk+1,ivar))/(Block1%dxk*Block1%xi_I(i)*sin(Block1%xj_I(j)))
                        end do
                    end do
                end do
            end select
            deallocate(flux,phi)
        end do
    end subroutine ModDiffusion_Artificial_1

end module ModDiffusion