module ModDiffusion

    use ModBlock,           only:   BlockType
    use ModParameters,      only:   ni,nj,nk,ng,ModelS_delta,nvar
    use ModLinReconstruct,  only:   ModLinReconstruct_minmod
    use ModVariables,       only:   rho1_,vr_,vt_,vp_,br_,bt_,bp_,s1_

    contains

    subroutine ModDiffusion_Aritificial_1(Block1,EQN_update_R,h,if_rk)
        implicit none
        type(BlockType),target          ::  Block1
        real,intent(inout)              ::  EQN_update_R(1:ni,1:nj,1:nk,1:nvar)
        integer,intent(in)              ::  h
        logical,intent(in)              ::  if_rk

        integer                         ::  direction1
        real                            ::  c_s(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng),&
                                            c(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        real,pointer                    ::  primitive(:,:,:,:)
        real                            ::  d_primitive(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,nvar)
        real,allocatable                ::  flux(:,:,:,:),phi(:,:,:,:)
        integer                         ::  i,j,ivar

        ! Get the sound speed
        c_s=1./Block1%Xi_rsst*sqrt(Block1%gamma1/ModelS_delta*Block1%p0_over_rho0)

        ! Get the primitive pointer based on if_rk
        if (if_rk) then
            primitive=>Block1%primitive_rk
        else
            primitive=>Block1%primitive
        end if

        do direction1=1,3

            ! get total speed c
            c=abs(primitive(:,:,:,direction1+1))+c_s

            ! use minmod to find \Delta u
            call ModLinReconstruct_minmod(nvar,ni,nj,nk,ng,direction1,primitive,d_primitive)

            ! get phi & flux
            select case(direction1)
            case(1)
                allocate(flux(ni+1,nj,nk,nvar),phi(ni+1,nj,nk,nvar))

                ! first get \Phi_{h}
                ! phi = max[0, 1+((u_r-u_l)/(u_i+1-u_i)-1)],
                ! for the i+1/2 face.
                phi(1:ni+1,:,:,:)=&
                    max(0.0,1.0+h*((primitive(1:ni+1,1:nj,1:nk,:)-d_primitive(1:ni+1,1:nj,1:nk,:)*0.5-&
                    primitive(0:ni,1:nj,1:nk,:)-d_primitive(0:ni,1:nj,1:nk,:)*0.5)/&
                    (primitive(1:ni+1,1:nj,1:nk,:)-primitive(0:ni,1:nj,1:nk,:))-1))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux(1:ni+1,:,:,:)=&
                    -0.5*phi*(primitive(1:ni+1,1:nj,1:nk,:)-d_primitive(1:ni+1,1:nj,1:nk,:)*0.5-&
                    primitive(0:ni,1:nj,1:nk,:)-d_primitive(0:ni,1:nj,1:nk,:)*0.5)
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,nvar
                    flux(:,:,:,ivar)=flux(:,:,:,ivar)*(c(1:ni+1,:,:)+c(0:ni,:,:))*0.5
                end do

                !do ivar=vr_,vp_
                !    flux(:,:,:,ivar)=flux(:,:,:,ivar)+&
                !        (primitive(1:ni+1,:,:,ivar)+primitive(0:ni,:,:,ivar))*&
                !        0.5*flux(:,:,:,rho1_)
                !end do

                ! update EQN_update_R
                do ivar=1,nvar
                    EQN_update_R(:,:,:,ivar)=EQN_update_R(:,:,:,ivar)+&
                        (flux(1:ni,:,:,ivar)-flux(2:ni+1,:,:,ivar))/Block1%dxi
                end do
            case(2)
                allocate(flux(ni,nj+1,nk,nvar),phi(ni,nj+1,nk,nvar))

                ! first get \Phi_{h}
                ! phi = max[0, 1+((u_r-u_l)/(u_i+1-u_i)-1)],
                ! for the i+1/2 face.
                phi(:,1:nj+1,:,:)=&
                    max(0.0,1.0+h*((primitive(1:ni,1:nj+1,1:nk,:)-d_primitive(1:ni,1:nj+1,1:nk,:)*0.5-&
                    primitive(1:ni,0:nj,1:nk,:)-d_primitive(1:ni,0:nj,1:nk,:)*0.5)/&
                    (primitive(1:ni,1:nj+1,1:nk,:)-primitive(1:ni,0:nj,1:nk,:))-1))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux(:,1:nj+1,:,:)=&
                    -0.5*phi*(primitive(1:ni,1:nj+1,1:nk,:)-d_primitive(1:ni,1:nj+1,1:nk,:)*0.5-&
                    primitive(1:ni,0:nj,1:nk,:)-d_primitive(1:ni,0:nj,1:nk,:)*0.5)
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,nvar
                    flux(:,:,:,ivar)=flux(:,:,:,ivar)*(c(:,1:nj+1,:)+c(:,0:nj,:))*0.5
                end do

                !do ivar=vr_,vp_
                !    flux(:,:,:,ivar)=flux(:,:,:,ivar)+&
                !        (primitive(:,1:nj+1,:,ivar)+primitive(:,0:nj,:,ivar))*&
                !        0.5*flux(:,:,:,rho1_)
                !end do

                do ivar=1,nvar
                    do i=1,ni
                        EQN_update_R(i,:,:,ivar)=EQN_update_R(i,:,:,ivar)+&
                            (flux(i,1:nj,:,ivar)-flux(i,2:nj+1,:,ivar))/(Block1%dxj*Block1%xi(i))
                    end do
                end do
            case(3)
                allocate(flux(ni,nj,nk+1,nvar),phi(ni,nj,nk+1,nvar))

                ! first get \Phi_{h}
                ! phi = max[0, 1+((u_r-u_l)/(u_i+1-u_i)-1)],
                ! for the i+1/2 face.
                phi(:,:,1:nk+1,:)=&
                    max(0.0,1.0+h*((primitive(1:ni,1:nj,1:nk+1,:)-d_primitive(1:ni,1:nj,1:nk+1,:)*0.5-&
                    primitive(1:ni,1:nj,0:nk,:)-d_primitive(1:ni,1:nj,0:nk,:)*0.5)/&
                    (primitive(1:ni,1:nj,1:nk+1,:)-primitive(1:ni,1:nj,0:nk,:))-1))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux(:,:,1:nk+1,:)=&
                    -0.5*phi*(primitive(1:ni,1:nj,1:nk+1,:)-d_primitive(1:ni,1:nj,1:nk+1,:)*0.5-&
                    primitive(1:ni,1:nj,0:nk,:)-d_primitive(1:ni,1:nj,0:nk,:)*0.5)
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,nvar
                    flux(:,:,:,ivar)=flux(:,:,:,ivar)*(c(:,:,1:nk+1)+c(:,:,0:nk))*0.5
                end do

                !do ivar=vr_,vp_
                !    flux(:,:,:,ivar)=flux(:,:,:,ivar)+&
                !        (primitive(:,:,1:nk+1,ivar)+primitive(:,:,0:nk,ivar))*&
                !        0.5*flux(:,:,:,rho1_)
                !end do

                do ivar=1,nvar
                    do j=1,nj; do i=1,ni
                        EQN_update_R(i,j,:,ivar)=EQN_update_R(i,j,:,ivar)+&
                            (flux(i,j,1:nk,ivar)-flux(i,j,2:nk+1,ivar))/(Block1%dxk*Block1%xi(i)*sin(Block1%xj(j)))
                    end do; end do
                end do
            end select
            deallocate(flux,phi)
        end do
    end subroutine ModDiffusion_Aritificial_1

end module ModDiffusion