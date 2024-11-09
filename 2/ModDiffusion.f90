module ModDiffusion

    use ModLinReconstruct

    contains

    subroutine ModDiffusion_Aritificial_1(primitive,ni,nj,nk,ng,xi,xj,dxi,dxj,dxk,EQN_update_R,&
        p0,rho0,h)

        use ModParameter, only : paraXi,paraGamma,paraDelta_r
        implicit none
        integer,intent(in)              ::  ni,nj,nk,ng
        integer,intent(in)              ::  h
        real,intent(in)                 ::  dxi,dxj,dxk
        real,intent(in)                 ::  primitive(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,5) 
        real,intent(in)                 ::  p0(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng),&
                                            rho0(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        real,intent(in)                 ::  xi(-ng+1:ng+ni),xj(-ng+1:ng+nj)

        integer                         ::  direction1,direction2
        real                            ::  c_s(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng),&
                                            c(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        real                            ::  d_primitive(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,5)
        real,allocatable                ::  flux(:,:,:,:),phi(:,:,:,:)
        real,intent(inout)              ::  EQN_update_R(1:ni,1:nj,1:nk,1:5)
        integer                         ::  i,j

        integer ::  ivar

        c_s=1./paraXi*sqrt(paraGamma/(8.*paraDelta_r))*p0/rho0

        do direction1=1,3

            ! get total speed c
            c=abs(primitive(:,:,:,direction1+1))+c_s

            ! use minmod to find \Delta u
            call ModLinReconstruct_minmod(5,ni,nj,nk,ng,direction1,primitive,d_primitive)

            ! get phi & flux
            select case(direction1)
            case(1)
                allocate(flux(ni+1,nj,nk,5),phi(ni+1,nj,nk,5))

                ! first get \Phi_{h}
                ! phi = max[0, 1+((u_r-u_l)/(u_i+1-u_i)-1)],
                ! for the i+1/2 face.
                phi(1:ni+1,:,:,:)=&
                    max(0.0,1.0+h*((primitive(1:ni+1,1:nj,1:nk,:)-d_primitive(1:ni+1,1:nj,1:nk,:)*0.5-&
                    primitive(0:ni,1:nj,1:nk,:)-d_primitive(0:ni,1:nj,1:nk,:))/&
                    (primitive(1:ni+1,1:nj,1:nk,:)-primitive(0:ni,1:nj,1:nk,:))-1))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux(1:ni+1,:,:,:)=&
                    -0.5*phi*(primitive(1:ni+1,1:nj,1:nk,:)-d_primitive(1:ni+1,1:nj,1:nk,:)*0.5-&
                    primitive(0:ni,1:nj,1:nk,:)-d_primitive(0:ni,1:nj,1:nk,:)*0.5)
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,5
                    flux(:,:,:,ivar)=flux(:,:,:,ivar)*(c(1:ni+1,:,:)+c(0:ni,:,:))*0.5
                end do

                do direction2=1,3
                    flux(:,:,:,direction2+1)=flux(:,:,:,direction2+1)+&
                        (primitive(1:ni+1,:,:,direction2+1)+primitive(0:ni,:,:,direction2+1))*&
                        0.5*flux(:,:,:,1)
                end do

                ! update EQN_update_R
                do ivar=1,5
                    EQN_update_R(:,:,:,ivar)=EQN_update_R(:,:,:,ivar)+&
                        (flux(1:ni,:,:,ivar)-flux(2:ni+1,:,:,ivar))/dxi
                end do
            case(2)
                allocate(flux(ni,nj+1,nk,5),phi(ni,nj+1,nk,5))

                ! first get \Phi_{h}
                ! phi = max[0, 1+((u_r-u_l)/(u_i+1-u_i)-1)],
                ! for the i+1/2 face.
                phi(:,1:nj+1,:,:)=&
                    max(0.0,1.0+h*((primitive(1:ni,1:nj+1,1:nk,:)-d_primitive(1:ni,1:nj+1,1:nk,:)*0.5-&
                    primitive(1:ni,0:nj,1:nk,:)-d_primitive(1:ni,0:nj,1:nk,:))/&
                    (primitive(1:ni,1:nj+1,1:nk,:)-primitive(1:ni,0:nj,1:nk,:))-1))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux(:,1:nj+1,:,:)=&
                    -0.5*phi*(primitive(1:ni,1:nj+1,1:nk,:)-d_primitive(1:ni,1:nj+1,1:nk,:)*0.5-&
                    primitive(1:ni,0:nj,1:nk,:)-d_primitive(1:ni,0:nj,1:nk,:)*0.5)
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,5
                    flux(:,:,:,ivar)=flux(:,:,:,ivar)*(c(:,1:nj+1,:)+c(:,0:nj,:))*0.5
                end do

                do direction2=1,3
                    flux(:,:,:,direction2+1)=flux(:,:,:,direction2+1)+&
                        (primitive(:,1:nj+1,:,direction2+1)+primitive(:,0:nj,:,direction2+1))*&
                        0.5*flux(:,:,:,1)
                end do

                do ivar=1,5
                    do i=1,ni
                        EQN_update_R(i,:,:,ivar)=EQN_update_R(i,:,:,ivar)+&
                            (flux(i,1:nj,:,ivar)-flux(i,2:nj+1,:,ivar))/(dxj*xi(i))
                    end do
                end do
            case(3)
                allocate(flux(ni,nj,nk+1,5),phi(ni,nj,nk+1,5))

                ! first get \Phi_{h}
                ! phi = max[0, 1+((u_r-u_l)/(u_i+1-u_i)-1)],
                ! for the i+1/2 face.
                phi(:,:,1:nk+1,:)=&
                    max(0.0,1.0+h*((primitive(1:ni,1:nj,1:nk+1,:)-d_primitive(1:ni,1:nj,1:nk+1,:)*0.5-&
                    primitive(1:ni,1:nj,0:nk,:)-d_primitive(1:ni,1:nj,0:nk,:))/&
                    (primitive(1:ni,1:nj,1:nk+1,:)-primitive(1:ni,1:nj,0:nk,:))-1))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux(:,:,1:nk+1,:)=&
                    -0.5*phi*(primitive(1:ni,1:nj,1:nk+1,:)-d_primitive(1:ni,1:nj,1:nk+1,:)*0.5-&
                    primitive(1:ni,1:nj,0:nk,:)-d_primitive(1:ni,1:nj,0:nk,:)*0.5)
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,5
                    flux(:,:,:,ivar)=flux(:,:,:,ivar)*(c(:,:,1:nk+1)+c(:,:,0:nk))*0.5
                end do

                do direction2=1,3
                    flux(:,:,:,direction2+1)=flux(:,:,:,direction2+1)+&
                        (primitive(:,:,1:nk+1,direction2+1)+primitive(:,:,0:nk,direction2+1))*&
                        0.5*flux(:,:,:,1)
                end do

                do ivar=1,5
                    do j=1,nj; do i=1,ni
                        EQN_update_R(i,j,:,ivar)=EQN_update_R(i,j,:,ivar)+&
                            (flux(i,j,1:nk,ivar)-flux(i,j,2:nk+1,ivar))/(dxk*xi(i)*sin(xj(j)))
                    end do; end do
                end do
            end select
            deallocate(flux,phi)
        end do
    end subroutine ModDiffusion_Aritificial_1

end module ModDiffusion