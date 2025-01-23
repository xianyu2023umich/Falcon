module ModEquation

    use ieee_arithmetic
    use ModBlock,       only:   BlockType
    use ModSpherical,   only:   ModSpherical_Grad_f,&
                                ModSpherical_A_dot_Grad_f,&
                                ModSpherical_div,&
                                ModSpherical_A_dot_nabla_B
    use ModBoundary,    only:   ModBoundary_Dynamo_HD_primitives,&
                                ModBoundary_Dynamo_HD_p1
    use ModDiffusion,   only:   ModDiffusion_Aritificial_1
    use ModParameters,  only:   ni,nj,nk,ng,ModelS_delta,ModelS_heating_ratio

    contains 

    ! The dynamo model from Hotta 2014.
    subroutine ModEquation_Dynamo_HD(Block1,if_rk,EQN_update_R)
        implicit none
        type(BlockType),target      ::  Block1                      
        logical,intent(in)          ::  if_rk
        real,pointer                ::  primitive(:,:,:,:)
        integer                     ::  direction
        real,intent(out)            ::  EQN_update_R(1:ni,1:nj,1:nk,1:5)
        real                        ::  tmp(-ng+1:ng+ni,-ng+1:ng+nj,-ng+1:ng+nk,3)
        
        ! If_rk
        if (if_rk) then
            primitive=>Block1%primitive_rk
        else
            primitive=>Block1%primitive
        end if
        
        ! Perparations
        ! Get m; Set R=0.; Set p1; Set Boundary.
        EQN_update_R=0.
        
        Block1%p1=Block1%Gamma1*Block1%p0_over_rho0*primitive(:,:,:,1)+&
            Block1%Gamma3_minus_1*Block1%rho0T0*primitive(:,:,:,5)

        call ModBoundary_Dynamo_HD_primitives(Block1,if_rk)
        !call ModBoundary_Dynamo_HD_p1(Block1)
        
        ! EQN 1
        do direction=1,3
            tmp(:,:,:,direction)=Block1%rho0*primitive(:,:,:,direction+1)
        end do
        EQN_update_R(:,:,:,1)=-1.0/(Block1%Xi_rsst(1:ni,1:nj,1:nk)**2*ModelS_delta)*&
            ModSpherical_div(ni,nj,nk,ng,Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,tmp)
        
        ! EQN 2-4 Inertial Force
        EQN_update_R(:,:,:,2:4)=-&
            ModSpherical_A_dot_nabla_B(ni,nj,nk,ng,&
            Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,&
            primitive(:,:,:,2:4),primitive(:,:,:,2:4))
        
        ! EQN 2-4 pressure gradient
        tmp(1:ni,1:nj,1:nk,:)=&
            ModSpherical_Grad_f(ni,nj,nk,ng,&
            Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,Block1%p1)
        do direction=1,3
            EQN_update_R(:,:,:,direction+1)=EQN_update_R(:,:,:,direction+1)-&
                tmp(1:ni,1:nj,1:nk,direction)/&
                Block1%rho0(1:ni,1:nj,1:nk)
        end do
        
        ! EQN 2 Gravity
        EQN_update_R(:,:,:,2)=EQN_update_R(:,:,:,2)-&
            Block1%g_over_rho0*primitive(1:ni,1:nj,1:nk,1)

        ! EQN5 advection term of s1
        EQN_update_R(:,:,:,5)=-&
            ModSpherical_A_dot_Grad_f(ni,nj,nk,ng,&
            Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,&
            primitive(:,:,:,2:4),primitive(:,:,:,5))
        
        ! EQN5 heating
        EQN_update_R(:,:,:,5)=EQN_update_R(:,:,:,5)+&
            ModelS_heating_ratio*(Block1%total_heat(1:ni,1:nj,1:nk))/Block1%rho0T0(1:ni,1:nj,1:nk)

        ! Aritificial diffusion
        call ModDiffusion_Aritificial_1(Block1,EQN_update_R,2,if_rk)
    end subroutine ModEquation_Dynamo_HD
end module ModEquation