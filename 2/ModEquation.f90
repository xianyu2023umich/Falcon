module ModEquation

    use ieee_arithmetic
    use ModBlock
    use ModSSM_v0
    use ModSpherical
    use ModBoundary
    use ModDiffusion

    contains 

    ! The dynamo model from Hotta 2019.
    subroutine ModEquation_Dynamo_HD(Block1,if_rk,EQN_update_R)

        use ModParameter, only : paraXi,paraGamma,paraDelta_r

        implicit none
        type(Block),target          ::  Block1                      
        logical,intent(in)          ::  if_rk
        real,pointer                ::  primitive(:,:,:,:)
        integer                     ::  direction
        real                        ::  m
        real,intent(out)            ::  EQN_update_R(1:Block1%ni,1:Block1%nj,1:Block1%nk,1:5)
        real                        ::  tmp(-Block1%ng+1:Block1%ng+Block1%ni,&
                                            -Block1%ng+1:Block1%ng+Block1%nj,&
                                            -Block1%ng+1:Block1%ng+Block1%nk,3)
        
        ! If_rk
        if (if_rk) then
            primitive=>Block1%primitive_rk
        else
            primitive=>Block1%primitive
        end if
        
        ! Perparations
        ! Get m; Set R=0.; Set p1; Set Boundary.
        m=1./(paraGamma-1.)
        EQN_update_R=0.
        Block1%p1=paraGamma*Block1%p0*(primitive(:,:,:,1)/Block1%rho0+primitive(:,:,:,5))
        call ModBoundary_Dynamo_HD_primitives(Block1)
        
        ! EQN 1
        do direction=1,3
            tmp(:,:,:,direction)=Block1%rho0*primitive(:,:,:,direction+1)
        end do
        EQN_update_R(:,:,:,1)=(-1./paraXi**2)*(1./paraDelta_r/8.)*&
            ModSpherical_div(Block1%ni,Block1%nj,Block1%nk,Block1%ng,&
            Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,tmp)

        ! EQN 2-4 Inertial Force
        EQN_update_R(:,:,:,2:4)=-&
            ModSpherical_A_dot_nabla_B(Block1%ni,Block1%nj,Block1%nk,Block1%ng,&
            Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,&
            primitive(:,:,:,2:4),primitive(:,:,:,2:4))
        
        ! EQN 2-4 pressure gradient
        tmp(1:Block1%ni,1:Block1%nj,1:Block1%nk,:)=&
            ModSpherical_Grad_f(Block1%ni,Block1%nj,Block1%nk,Block1%ng,&
            Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,Block1%p1)
        do direction=1,3
            EQN_update_R(:,:,:,direction+1)=EQN_update_R(:,:,:,direction+1)-&
                tmp(1:Block1%ni,1:Block1%nj,1:Block1%nk,direction)/&
                Block1%rho0(1:Block1%ni,1:Block1%nj,1:Block1%nk)
        end do
        
        ! EQN 2 Gravity
        EQN_update_R(:,:,:,2)=EQN_update_R(:,:,:,2)-&
            primitive(1:Block1%ni,1:Block1%nj,1:Block1%nk,1)/&
            Block1%rho0(1:Block1%ni,1:Block1%nj,1:Block1%nk)

        ! EQN5 advection term of s1
        EQN_update_R(:,:,:,5)=-&
            ModSpherical_A_dot_Grad_f(Block1%ni,Block1%nj,Block1%nk,Block1%ng,&
            Block1%xi,Block1%xj,Block1%dxi,Block1%dxj,Block1%dxk,&
            primitive(:,:,:,2:4),primitive(:,:,:,5))
        
        ! EQN5 advection term of s0
        EQN_update_R(:,:,:,5)=EQN_update_R(:,:,:,5)+&
            primitive(1:Block1%ni,1:Block1%nj,1:Block1%nk,2)/&
            Block1%p0(1:Block1%ni,1:Block1%nj,1:Block1%nk)/8.
        
        ! Aritificial diffusion
        call ModDiffusion_Aritificial_1(primitive,&
            Block1%ni,Block1%nj,Block1%nk,Block1%ng,Block1%xi,Block1%xj,&
            Block1%dxi,Block1%dxj,Block1%dxk,EQN_update_R,Block1%p0,Block1%rho0,2)
        
    end subroutine ModEquation_Dynamo_HD
end module ModEquation