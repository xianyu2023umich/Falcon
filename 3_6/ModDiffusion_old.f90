module ModDiffusion

    use ModBlock,           only:   BlockType
    use ModParameters,      only:   ni,nj,nk,ng,nvar,iEquation,&
                                    if_involve_B,if_account_diffused_energy,diffusion_h
    use ModLinReconstruct,  only:   ModLinReconstruct_minmod

    implicit none

    private

    public :: ModDiffusion_DoAll, &
              ModDiffusion_HyperArtificial

    ! The n of vars for diffusion
    integer :: nvar_diffusion

    ! Allocate several arrays from the start so no need to allocate/deallocate.
    ! Probably no need for phi.
    real(8), allocatable                ::  SpecDiffHeat(:,:,:)    
    real(8), allocatable                ::  d_primitive(:,:,:,:)
    real(8), allocatable                ::  flux_r(:,:,:,:), flux_t(:,:,:,:), flux_p(:,:,:,:)
    real(8), allocatable                ::  u_r_minus_u_l_r(:,:,:,:),&
                                            u_i_plus_1_minus_u_i_r(:,:,:,:),&
                                            u_r_minus_u_l_t(:,:,:,:),&
                                            u_i_plus_1_minus_u_i_t(:,:,:,:),&
                                            u_r_minus_u_l_p(:,:,:,:),&
                                            u_i_plus_1_minus_u_i_p(:,:,:,:)

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
        allocate(flux_r(ni+1,nj,nk,nvar_diffusion))
        allocate(flux_t(ni,nj+1,nk,nvar_diffusion))
        allocate(flux_p(ni,nj,nk+1,nvar_diffusion))
        allocate(u_r_minus_u_l_r(ni+1,nj,nk,nvar_diffusion))
        allocate(u_r_minus_u_l_t(ni,nj+1,nk,nvar_diffusion))
        allocate(u_r_minus_u_l_p(ni,nj,nk+1,nvar_diffusion))
        allocate(u_i_plus_1_minus_u_i_r(ni+1,nj,nk,nvar_diffusion))
        allocate(u_i_plus_1_minus_u_i_t(ni,nj+1,nk,nvar_diffusion))
        allocate(u_i_plus_1_minus_u_i_p(ni,nj,nk+1,nvar_diffusion))
        if_diffusion_allocated = .true.
        if (if_account_diffused_energy) then
            allocate(SpecDiffHeat(1:ni,1:nj,1:nk))
        end if
    end subroutine ModDiffusion_Allocate

    subroutine ModDiffusion_DeAllocate
        deallocate(d_primitive)
        deallocate(flux_r,flux_t,flux_p)
        deallocate(u_r_minus_u_l_r,u_i_plus_1_minus_u_i_r)
        deallocate(u_r_minus_u_l_t,u_i_plus_1_minus_u_i_t)
        deallocate(u_r_minus_u_l_p,u_i_plus_1_minus_u_i_p)
        if (allocated(SpecDiffHeat)) deallocate(SpecDiffHeat)
        if_diffusion_allocated = .false.
    end subroutine ModDiffusion_DeAllocate

    subroutine ModDiffusion_HyperArtificial(Block1)
        implicit none
        type(BlockType),target          ::  Block1

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
                u_i_plus_1_minus_u_i_r=Block1%primitive(1:ni+1,1:nj,1:nk,:)-&
                    Block1%primitive(0:ni,1:nj,1:nk,:)
                
                u_r_minus_u_l_r=u_i_plus_1_minus_u_i_r-0.5*&
                    (d_primitive(1:ni+1,1:nj,1:nk,:)+d_primitive(0:ni,1:nj,1:nk,:))

                flux_r=max(0.0,1.0+diffusion_h*(u_r_minus_u_l_r/u_i_plus_1_minus_u_i_r-1))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux_r=-0.5*flux_r*u_r_minus_u_l_r
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,nvar_diffusion
                    flux_r(:,:,:,ivar)=flux_r(:,:,:,ivar)*(c(1:ni+1,:,:)+c(0:ni,:,:))*0.5

                    ! face area
                    flux_r(:,:,:,ivar)=flux_r(:,:,:,ivar)*Block1%Si_FLL
                end do

                ! update EQN_update_R
                do ivar=1,nvar_diffusion
                    Block1%EQN_update_R_IV(:,:,:,ivar)=Block1%EQN_update_R_IV(:,:,:,ivar)+&
                        (flux_r(1:ni,:,:,ivar)-flux_r(2:ni+1,:,:,ivar))/Block1%V_LLL
                end do
            case(2)
                ! First get u differences.
                u_i_plus_1_minus_u_i_t=Block1%primitive(1:ni,1:nj+1,1:nk,:)-&
                    Block1%primitive(1:ni,0:nj,1:nk,:)
                
                u_r_minus_u_l_t=u_i_plus_1_minus_u_i_t-0.5*&
                    (d_primitive(1:ni,1:nj+1,1:nk,:)+d_primitive(1:ni,0:nj,1:nk,:))

                flux_t=max(0.0,1.0+diffusion_h*(u_r_minus_u_l_t/u_i_plus_1_minus_u_i_t-1))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux_t=-0.5*flux_t*u_r_minus_u_l_t

                ! multiply flux by c_{i+1/2} and face area
                do ivar=1,nvar_diffusion
                    flux_t(:,:,:,ivar)=flux_t(:,:,:,ivar)*(c(:,1:nj+1,:)+c(:,0:nj,:))*0.5

                    ! face area
                    flux_t(:,:,:,ivar)=flux_t(:,:,:,ivar)*Block1%Sj_LFL
                end do

                do ivar=1,nvar_diffusion
                    Block1%EQN_update_R_IV(:,:,:,ivar)=Block1%EQN_update_R_IV(:,:,:,ivar)+&
                            (flux_t(:,1:nj,:,ivar)-flux_t(:,2:nj+1,:,ivar))/Block1%V_LLL
                end do
            case(3)
                ! First get u differences.
                u_i_plus_1_minus_u_i_p=Block1%primitive(1:ni,1:nj,1:nk+1,:)-&
                    Block1%primitive(1:ni,1:nj,0:nk,:)
                
                u_r_minus_u_l_p=u_i_plus_1_minus_u_i_p-0.5*&
                    (d_primitive(1:ni,1:nj,1:nk+1,:)+d_primitive(1:ni,1:nj,0:nk,:))

                flux_p=max(0.0,1.0+diffusion_h*(u_r_minus_u_l_p/u_i_plus_1_minus_u_i_p-1))
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux_p=-0.5*flux_p*u_r_minus_u_l_p
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,nvar_diffusion
                    flux_p(:,:,:,ivar)=flux_p(:,:,:,ivar)*(c(:,:,1:nk+1)+c(:,:,0:nk))*0.5

                    ! face area
                    flux_p(:,:,:,ivar)=flux_p(:,:,:,ivar)*Block1%Sk_LLF
                end do

                do ivar=1,nvar_diffusion
                    Block1%EQN_update_R_IV(:,:,:,ivar)=Block1%EQN_update_R_IV(:,:,:,ivar)+&
                        (flux_p(:,:,1:nk,ivar)-flux_p(:,:,2:nk+1,ivar))/Block1%V_LLL
                end do
            end select
        end do

        if (if_account_diffused_energy) then
            call ModDiffusion_AccountDiffusedEnergy(Block1)
        end if
    end subroutine ModDiffusion_HyperArtificial

    ! The artificial diffusoin dissipates energy.
    ! This subroutine estimate this dissipated energy and
    ! account it as a heating, thus appears in the perturbed
    ! entropy equation.
    !
    ! For a velocity flux f_i (i for direction), the
    ! kinetic energy change is f_i*(u_r-u_l)*rho0. We should
    ! minus this change in the entropy equation, which appears as
    ! -f_i*(u_r-u_l)*rho0/(rho0*T0)=-f_i*(u_r-u_l)/T0.
    !
    ! Also we need to partition this energy into two cells at two
    ! sides of the face, so there's a 0.5. Also divide the volume.
    ! So the final formula is -0.5*f_i*(u_r-u_l)/(T0*V_LLL).
    subroutine ModDiffusion_AccountDiffusedEnergy(Block1)
        implicit none
        type(BlockType),target          ::  Block1

        ! The left face of all cells is from 1 to ni.
        ! The right is from 2 to ni+1.

        SpecDiffHeat=0.d0

        SpecDiffHeat=SpecDiffHeat-flux_r(1:ni,:,:,Block1%vr_)*&
            (u_r_minus_u_l_r(1:ni,:,:,Block1%vr_))
        SpecDiffHeat=SpecDiffHeat-flux_r(2:ni+1,:,:,Block1%vr_)*&
            (u_r_minus_u_l_r(2:ni+1,:,:,Block1%vr_))
        SpecDiffHeat=SpecDiffHeat-flux_t(:,1:nj,:,Block1%vr_)*&
            (u_r_minus_u_l_t(:,1:nj,:,Block1%vr_))
        SpecDiffHeat=SpecDiffHeat-flux_t(:,2:nj+1,:,Block1%vr_)*&
            (u_r_minus_u_l_t(:,2:nj+1,:,Block1%vr_))
        SpecDiffHeat=SpecDiffHeat-flux_p(:,:,1:nk,Block1%vr_)*&
            (u_r_minus_u_l_p(:,:,1:nk,Block1%vr_))
        SpecDiffHeat=SpecDiffHeat-flux_p(:,:,2:nk+1,Block1%vr_)*&
            (u_r_minus_u_l_p(:,:,2:nk+1,Block1%vr_))

        SpecDiffHeat=SpecDiffHeat*0.5d0
        Block1%EQN_update_R_IV(:,:,:,Block1%s1_)=&
            Block1%EQN_update_R_IV(:,:,:,Block1%s1_)+&
            SpecDiffHeat/(Block1%T0_LLL*Block1%V_LLL)
    end subroutine  ModDiffusion_AccountDiffusedEnergy

end module ModDiffusion
