module ModTimeStep

    use ModBlock,       only:   BlockType
    use ModYinYangTree, only:   YYTree
    use ModParameters,  only:   ModelS_delta,iGeometry
    use ModVariables,   only:   vr_,vt_,vp_,br_,bt_,bp_

    contains

    subroutine ModTimeStep_Dynamo_HD(Tree,CFL_ad,dt)
        implicit none
        type(YYTree),intent(in),target  ::  Tree                ! Tree
        real,intent(in)                 ::  CFL_ad              ! CFL
        type(BlockType),pointer         ::  Block1              ! Block pointer
        integer                         ::  iBlock              ! i of block   
        real                            ::  dt1                 ! dt of 1 Block
        real,intent(out)                ::  dt                  ! dt considering all local blocks

        ! Initialize dt to be very big
        dt=1.e20
        dt1=1.e20

        ! Loop all the blocks to get minimum dt.
        do iBlock=1,size(Tree%LocalBlocks) 
            Block1 => Tree%LocalBlocks(iBlock)

            call ModTimeStep_Dynamo_HD_single(Block1,CFL_ad,dt1)
            dt=min(dt,dt1)
        end do
    end subroutine ModTimeStep_Dynamo_HD

    subroutine ModTimeStep_Dynamo_HD_single(Block1,CFL_ad,dt)
        implicit none
        type(BlockType),intent(inout)   ::  Block1              ! The block
        real,intent(in)                 ::  CFL_ad
        real                            ::  dl_min
        real,intent(inout)              ::  dt

        select case(iGeometry)
        case(1)
            ! Wave speed
            Block1%v_wave_III=&
                    1./Block1%Xi_rsst_III*sqrt(Block1%gamma1_III/ModelS_delta*Block1%p0_over_rho0_III) + &
                    sqrt(   Block1%primitive_IV(:,:,:,vr_)**2  + &
                            Block1%primitive_IV(:,:,:,vt_)**2  + &
                            Block1%primitive_IV(:,:,:,vp_)**2)

            dl_min=min(Block1%dxi,Block1%dxj*minval(Block1%xi_I),Block1%dxk*minval(Block1%xi_I))
            dt = CFL_ad * dl_min / maxval(Block1%v_wave_III)
            !print *,maxval(Block1%Xi_rsst_III),minval(Block1%Xi_rsst_III)
        end select
    end subroutine ModTimeStep_Dynamo_HD_single

    subroutine ModTimeStep_Dynamo_MHD(Tree,CFL_ad,dt)
        implicit none
        type(YYTree),intent(in),target  ::  Tree                ! Tree
        real,intent(in)                 ::  CFL_ad              ! CFL
        type(BlockType),pointer         ::  Block1              ! Block pointer
        integer                         ::  iBlock              ! i of block   
        real                            ::  dt1                 ! dt of 1 Block
        real,intent(out)                ::  dt                  ! dt considering all local blocks

        ! Initialize dt to be very big
        dt=1.e20
        dt1=1.e20

        ! Loop all the blocks to get minimum dt.
        do iBlock=1,size(Tree%LocalBlocks) 
            Block1 => Tree%LocalBlocks(iBlock)

            call ModTimeStep_Dynamo_MHD_single(Block1,'spherical',CFL_ad,dt1)
            dt=min(dt,dt1)
        end do
    end subroutine ModTimeStep_Dynamo_MHD

    subroutine ModTimeStep_Dynamo_MHD_single(Block1,geometry,CFL_ad,dt)
        implicit none
        type(BlockType),intent(inout)   ::  Block1              ! The block
        character(len=*),intent(in)     ::  geometry
        real,intent(in)                 ::  CFL_ad
        real                            ::  dl_min
        real,intent(inout)              ::  dt

        select case(geometry)
        case('spherical')
            ! Wave speed
            Block1%v_wave_III=&
                    1./Block1%Xi_rsst_III*sqrt(Block1%gamma1_III/ModelS_delta*Block1%p0_over_rho0_III) + &
                    sqrt(   Block1%primitive_IV(:,:,:,vr_)**2  + &
                            Block1%primitive_IV(:,:,:,vt_)**2  + &
                            Block1%primitive_IV(:,:,:,vp_)**2) + &
                    sqrt(   Block1%primitive_IV(:,:,:,br_)**2  + &
                            Block1%primitive_IV(:,:,:,bt_)**2  + &
                            Block1%primitive_IV(:,:,:,bp_)**2)/sqrt(Block1%rho0_III)

            dl_min=min(Block1%dxi,Block1%dxj*minval(Block1%xi_I),Block1%dxk*minval(Block1%xi_I))
            dt = CFL_ad * dl_min / maxval(Block1%v_wave_III)
        end select
    end subroutine ModTimeStep_Dynamo_MHD_single
end module ModTimeStep