module ModTimeStep

    use ieee_arithmetic
    use ModOcTree
    use ModSSM_v0

    contains

    subroutine ModTimeStep_Dynamo_HD(Tree,CFL_ad,CFL_df,dt)
        implicit none

        type(OcTree),intent(in),target  :: Tree             ! Tree
        real,intent(in)                 :: CFL_ad,CFL_df    ! CFL

        type(Block),pointer             :: Block1           ! Block pointer
        integer                         :: iBlock           ! i of block

        real                            :: dt1              ! dt of 1 Block
        real,intent(out)                :: dt               ! dt considering all local blocks

        dt=1.e20

        do iBlock=1,size(Tree%LocalBlocks) 
            Block1 => Tree%LocalBlocks(iBlock)

            call ModTimeStep_Dynamo_HD_Block(Block1%primitive,Block1%ni,Block1%nj,Block1%nk,&
                Block1%ng,Block1%dxi,Block1%dxj,Block1%dxk,Block1%xk,'cartesian',CFL_ad,CFL_df,dt1)

            dt=min(dt,dt1)
        end do
    end subroutine ModTimeStep_Dynamo_HD

    subroutine ModTimeStep_Dynamo_HD_Block(primitive,ni,nj,nk,ng,dxi,dxj,dxk,xk,geometry,&
        CFL_ad,CFL_df,dt)

        use ModParameter, only : paraXi,paraGamma,paraDelta_r,paraRe,paraPr

        implicit none

        character(len=*),intent(in) :: geometry
        integer,intent(in) :: ni,nj,nk,ng
        real,intent(in) :: dxi,dxj,dxk
        real,intent(in) :: xk(-ng+1:nk+ng)
        real,intent(in) :: primitive(5,-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        real,intent(in) :: CFL_ad,CFL_df

        real :: rho0_list(-ng+1:nk+ng),p0_list(-ng+1:nk+ng)
        real :: m
        real :: dt_adv,dt_diff

        real,intent(out) :: dt

        m=1./(paraGamma-1.)

        rho0_list=ModSSM_v0_get_var0(nk+2*ng,xk,m,'rho0')
        p0_list=ModSSM_v0_get_var0(nk+2*ng,xk,m,'p0')

        select case(geometry)
        case('cartesian')
            dt_adv = min(dxi,dxj,dxk)/&
                (maxval(abs(primitive(2:4,:,:,:)))+&
                1./paraXi*sqrt(paraGamma/(8.*paraDelta_r))*maxval(p0_list/rho0_list))

            dt_diff = min(dxi**2,dxj**2,dxk**2) * paraRe * min(1.,paraPr,1./paraPr)

            dt=min(CFL_ad*dt_adv,CFL_df*dt_diff)
        end select

    end subroutine ModTimeStep_Dynamo_HD_Block
end module ModTimeStep