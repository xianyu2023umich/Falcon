module ModTimeStep

    use ModBlock,       only:   BlockType
    use ModYinYangTree, only:   YYTree
    use ModParameters,  only:   ModelS_delta
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

        ! Loop all the blocks to get minimum dt.
        do iBlock=1,size(Tree%LocalBlocks) 
            Block1 => Tree%LocalBlocks(iBlock)

            call ModTimeStep_Dynamo_HD_single(Block1,'spherical',CFL_ad,dt1)
            dt=min(dt,dt1)
        end do
    end subroutine ModTimeStep_Dynamo_HD

    subroutine ModTimeStep_Dynamo_HD_single(Block1,geometry,CFL_ad,dt)
        implicit none
        type(BlockType),intent(in)      ::  Block1              ! The block
        character(len=*),intent(in)     ::  geometry
        real,intent(in)                 ::  CFL_ad
        real                            ::  dt_adv
        real,intent(out)                ::  dt

        select case(geometry)
        case('spherical')
            dt_adv = min(Block1%dxi,Block1%dxj*minval(Block1%xi),Block1%dxk*minval(Block1%xi))/&
                (maxval(abs(Block1%primitive(:,:,:,vr_:vp_)))+&
                 maxval(1./Block1%Xi_rsst(:,1,1)*sqrt(Block1%Gamma1_list/ModelS_delta*Block1%p0_list/Block1%rho0_list)))

            dt = CFL_ad * dt_adv
        end select
    end subroutine ModTimeStep_Dynamo_HD_single
end module ModTimeStep