module ModDiffusion

    use ModBlock,           only:   BlockType
    use ModParameters,      only:   ni,nj,nk,ng,nvar,iEquation,if_involve_B
    use ModLinReconstruct,  only:   ModLinReconstruct_minmod

    implicit none

    private

    public :: ModDiffusion_DoAll, &
              ModDiffusion_Artificial_1

    ! The n of vars for diffusion
    integer :: nvar_diffusion

    ! Allocate several arrays from the start so no need to allocate/deallocate.
    ! Probably no need for phi.
    real(8), allocatable                ::  d_primitive(:,:,:,:)
    real(8), allocatable                ::  flux_x(:,:,:,:), flux_y(:,:,:,:), flux_z(:,:,:,:)
    real(8), allocatable                ::  u_r_minus_u_l_x(:,:,:,:),&
                                            u_i_plus_1_minus_u_i_x(:,:,:,:),&
                                            u_r_minus_u_l_y(:,:,:,:),&
                                            u_i_plus_1_minus_u_i_y(:,:,:,:),&
                                            u_r_minus_u_l_z(:,:,:,:),&
                                            u_i_plus_1_minus_u_i_z(:,:,:,:)

    ! The flag of whether above have been allocated.
    logical :: if_diffusion_allocated = .false.

    contains

    subroutine ModDiffusion_DoAll
        ! Get the n of var here
        if (iEquation == 1 .and. .not. if_involve_B) then
            nvar_diffusion = nvar - 4
        else
            nvar_diffusion = nvar
        end if

        ! Deallocate if already allocated, then allocate.
        ! This is to make sure if there's a change in nvar_diffusion.
        if (if_diffusion_allocated) call ModDiffusion_DeAllocate
        call ModDiffusion_Allocate
    end subroutine ModDiffusion_DoAll

    subroutine ModDiffusion_Allocate
        allocate(d_primitive(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,nvar_diffusion))
        allocate(flux_x(ni+1,nj,nk,nvar_diffusion))
        allocate(flux_y(ni,nj+1,nk,nvar_diffusion))
        allocate(flux_z(ni,nj,nk+1,nvar_diffusion))
        allocate(u_r_minus_u_l_x(ni+1,nj,nk,nvar_diffusion))
        allocate(u_r_minus_u_l_y(ni,nj+1,nk,nvar_diffusion))
        allocate(u_r_minus_u_l_z(ni,nj,nk+1,nvar_diffusion))
        allocate(u_i_plus_1_minus_u_i_x(ni+1,nj,nk,nvar_diffusion))
        allocate(u_i_plus_1_minus_u_i_y(ni,nj+1,nk,nvar_diffusion))
        allocate(u_i_plus_1_minus_u_i_z(ni,nj,nk+1,nvar_diffusion))
        if_diffusion_allocated = .true.
    end subroutine ModDiffusion_Allocate

    subroutine ModDiffusion_DeAllocate
        deallocate(d_primitive)
        deallocate(flux_x,flux_y,flux_z)
        deallocate(u_r_minus_u_l_x,u_i_plus_1_minus_u_i_x)
        deallocate(u_r_minus_u_l_y,u_i_plus_1_minus_u_i_y)
        deallocate(u_r_minus_u_l_z,u_i_plus_1_minus_u_i_z)
        if_diffusion_allocated = .false.
    end subroutine ModDiffusion_DeAllocate

    subroutine ModDiffusion_Artificial_1(Block1,h)
        implicit none
        type(BlockType),target          ::  Block1
        integer,intent(in)              ::  h

        integer                         ::  direction1
        real(8)                         ::  c(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        integer                         ::  ivar

        ! Get the primitive pointer
        Block1%primitive=>Block1%primitive_IV

        do direction1=1,3

            ! c=c_s+|v_i|
            c=abs(Block1%primitive(:,:,:,direction1+1))+Block1%v_sound_III

            ! use minmod to find \Delta u
            call ModLinReconstruct_minmod(nvar_diffusion,ni,nj,nk,ng,direction1,Block1%primitive,d_primitive)

            ! get phi & flux
            select case(direction1)
            case(1)
                ! First get u differences.
                u_i_plus_1_minus_u_i_x=Block1%primitive(1:ni+1,1:nj,1:nk,:)-&
                    Block1%primitive(0:ni,1:nj,1:nk,:)
                
                u_r_minus_u_l_x=u_i_plus_1_minus_u_i_x-0.5*&
                    (d_primitive(1:ni+1,1:nj,1:nk,:)+d_primitive(0:ni,1:nj,1:nk,:))

                flux_x=max(0.0,1.0+h*(u_r_minus_u_l_x/u_i_plus_1_minus_u_i_x-1))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux_x=-0.5*flux_x*u_r_minus_u_l_x
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,nvar_diffusion
                    flux_x(:,:,:,ivar)=flux_x(:,:,:,ivar)*(c(1:ni+1,:,:)+c(0:ni,:,:))*0.5

                    ! face area
                    flux_x(:,:,:,ivar)=flux_x(:,:,:,ivar)*Block1%Si_FLL
                end do

                ! update EQN_update_R
                do ivar=1,nvar_diffusion
                    Block1%EQN_update_R_IV(:,:,:,ivar)=Block1%EQN_update_R_IV(:,:,:,ivar)+&
                        (flux_x(1:ni,:,:,ivar)-flux_x(2:ni+1,:,:,ivar))/Block1%V_LLL
                end do
            case(2)
                ! First get u differences.
                u_i_plus_1_minus_u_i_y=Block1%primitive(1:ni,1:nj+1,1:nk,:)-&
                    Block1%primitive(1:ni,0:nj,1:nk,:)
                
                u_r_minus_u_l_y=u_i_plus_1_minus_u_i_y-0.5*&
                    (d_primitive(1:ni,1:nj+1,1:nk,:)+d_primitive(1:ni,0:nj,1:nk,:))

                flux_y=max(0.0,1.0+h*(u_r_minus_u_l_y/u_i_plus_1_minus_u_i_y-1))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux_y=-0.5*flux_y*u_r_minus_u_l_y

                ! multiply flux by c_{i+1/2} and face area
                do ivar=1,nvar_diffusion
                    flux_y(:,:,:,ivar)=flux_y(:,:,:,ivar)*(c(:,1:nj+1,:)+c(:,0:nj,:))*0.5

                    ! face area
                    flux_y(:,:,:,ivar)=flux_y(:,:,:,ivar)*Block1%Sj_LFL
                end do

                do ivar=1,nvar_diffusion
                    Block1%EQN_update_R_IV(:,:,:,ivar)=Block1%EQN_update_R_IV(:,:,:,ivar)+&
                            (flux_y(:,1:nj,:,ivar)-flux_y(:,2:nj+1,:,ivar))/Block1%V_LLL
                end do
            case(3)
                ! First get u differences.
                u_i_plus_1_minus_u_i_z=Block1%primitive(1:ni,1:nj,1:nk+1,:)-&
                    Block1%primitive(1:ni,1:nj,0:nk,:)
                
                u_r_minus_u_l_z=u_i_plus_1_minus_u_i_z-0.5*&
                    (d_primitive(1:ni,1:nj,1:nk+1,:)+d_primitive(1:ni,1:nj,0:nk,:))

                flux_z=max(0.0,1.0+h*(u_r_minus_u_l_z/u_i_plus_1_minus_u_i_z-1))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux_z=-0.5*flux_z*u_r_minus_u_l_z
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,nvar_diffusion
                    flux_z(:,:,:,ivar)=flux_z(:,:,:,ivar)*(c(:,:,1:nk+1)+c(:,:,0:nk))*0.5

                    ! face area
                    flux_z(:,:,:,ivar)=flux_z(:,:,:,ivar)*Block1%Sk_LLF
                end do

                do ivar=1,nvar_diffusion
                    Block1%EQN_update_R_IV(:,:,:,ivar)=Block1%EQN_update_R_IV(:,:,:,ivar)+&
                        (flux_z(:,:,1:nk,ivar)-flux_z(:,:,2:nk+1,ivar))/Block1%V_LLL
                end do
            end select
        end do
    end subroutine ModDiffusion_Artificial_1

end module ModDiffusion