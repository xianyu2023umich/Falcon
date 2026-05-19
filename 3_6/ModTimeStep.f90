module ModTimeStep

    use ModBlock,       only:   BlockType
    use ModYinYangTree, only:   YYTree
    use ModParameters,  only:   iGeometry,ni,nj,nk,MpiRank
    use ModConst,       only:   gamma_ideal_gas,k_B__CGS,m_p__CGS,mu0__CGS,R_sun__CGS

    contains

    subroutine ModTimeStep_Dynamo_HD(Tree,CFL_ad,dt)
        implicit none
        type(YYTree),intent(in),target  ::  Tree                ! Tree
        real(8),intent(in)              ::  CFL_ad              ! CFL
        type(BlockType),pointer         ::  Block1              ! Block pointer
        integer                         ::  iBlock              ! i of block   
        real(8)                         ::  dt1                 ! dt of 1 Block
        real(8),intent(out)             ::  dt                  ! dt considering all local blocks

        ! Initialize dt to be very big
        dt=1.e20
        dt1=1.e20

        ! Loop all the blocks to get minimum dt.
        do iBlock=1,size(Tree%LocalBlocks) 
            Block1 => Tree%LocalBlocks(iBlock)

            call ModTimeStep_Dynamo_HD_single(Block1,CFL_ad,dt1)
            dt=min(dt,dt1)

            if (dt .gt. 1.e10) then
                print *,'Detected big dt.'
                print *,MpiRank,iBlock,maxval(Block1%v_wave_III)
                print *,MpiRank,iBlock,maxval(Block1%Xi_rsst_III)
                print *,MpiRank,iBlock,minval(Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%vr_)**2  + &
                            Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%vt_)**2  + &
                            Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%vp_)**2)

                stop 1
            end if
        end do
    end subroutine ModTimeStep_Dynamo_HD

    subroutine ModTimeStep_Dynamo_HD_single(Block1,CFL_ad,dt)
        implicit none
        type(BlockType),intent(inout)   ::  Block1              ! The block
        real(8),intent(in)              ::  CFL_ad
        real(8)                         ::  dl_min
        real(8),intent(inout)           ::  dt

        select case(iGeometry)
        case(1)
            ! Wave speed
            Block1%v_wave_III(1:ni,1:nj,1:nk)=1./Block1%Xi_rsst_III(1:ni,1:nj,1:nk)*&
                    sqrt(Block1%gamma1_III(1:ni,1:nj,1:nk)*Block1%p0_over_rho0_III(1:ni,1:nj,1:nk)) + &
                    sqrt(   Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%vr_)**2  + &
                            Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%vt_)**2  + &
                            Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%vp_)**2)

            dl_min=min(Block1%dxi,Block1%dxj*minval(Block1%xi_I),Block1%dxk*minval(Block1%xi_I))
            dt = CFL_ad * dl_min / maxval(Block1%v_wave_III(1:ni,1:nj,1:nk))
            !print *,maxval(Block1%Xi_rsst_III),minval(Block1%Xi_rsst_III)
        end select
    end subroutine ModTimeStep_Dynamo_HD_single

    subroutine ModTimeStep_Dynamo_MHD(Tree,CFL_ad,dt)
        implicit none
        type(YYTree),intent(in),target  ::  Tree                ! Tree
        real(8),intent(in)              ::  CFL_ad              ! CFL
        type(BlockType),pointer         ::  Block1              ! Block pointer
        integer                         ::  iBlock              ! i of block   
        real(8)                         ::  dt1                 ! dt of 1 Block
        real(8),intent(out)             ::  dt                  ! dt considering all local blocks

        ! Initialize dt to be very big
        dt=1.e20
        dt1=1.e20

        ! Loop all the blocks to get minimum dt.
        do iBlock=1,size(Tree%LocalBlocks) 
            Block1 => Tree%LocalBlocks(iBlock)

            call ModTimeStep_Dynamo_MHD_single(Block1,1,CFL_ad,dt1)
            dt=min(dt,dt1)
        end do
    end subroutine ModTimeStep_Dynamo_MHD

    subroutine ModTimeStep_Dynamo_MHD_single(Block1,iGeometry,CFL_ad,dt)
        implicit none
        type(BlockType),intent(inout)   ::  Block1              ! The block
        integer,intent(in)              ::  iGeometry
        real(8),intent(in)              ::  CFL_ad
        real(8)                         ::  dl_min
        real(8),intent(inout)           ::  dt

        select case(iGeometry)
        case(1)
            ! Wave speed
            Block1%v_wave_III(1:ni,1:nj,1:nk)=1./Block1%Xi_rsst_III(1:ni,1:nj,1:nk)*&
                    sqrt(Block1%gamma1_III(1:ni,1:nj,1:nk)*Block1%p0_over_rho0_III(1:ni,1:nj,1:nk)) + &
                    sqrt(   Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%vr_)**2  + &
                            Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%vt_)**2  + &
                            Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%vp_)**2) + &
                    sqrt(   Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%br_)**2  + &
                            Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%bt_)**2  + &
                            Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%bp_)**2) / &
                    sqrt(Block1%rho0_III(1:ni,1:nj,1:nk))

            dl_min=min(Block1%dxi,Block1%dxj*minval(Block1%xi_I),Block1%dxk*minval(Block1%xi_I))
            dt = CFL_ad * dl_min / maxval(Block1%v_wave_III(1:ni,1:nj,1:nk))
        end select
    end subroutine ModTimeStep_Dynamo_MHD_single

    subroutine ModTimeStep_Corona_single(Block1,CFL_ad,dt)
        implicit none
        type(BlockType),intent(inout)   ::  Block1              ! The block
        real(8),intent(in)              ::  CFL_ad
        real(8)                         ::  dl_min
        real(8),intent(inout)           ::  dt

        select case(iGeometry)
        case(1)
            ! Wave speed
            Block1%v_wave_III(1:ni,1:nj,1:nk)=&
                    sqrt(gamma_ideal_gas*k_B__CGS/m_p__CGS*Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%te_)) + &
                    sqrt(   Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%vr_)**2  + &
                            Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%vt_)**2  + &
                            Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%vp_)**2) + &
                    sqrt(   Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%br_)**2  + &
                            Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%bt_)**2  + &
                            Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%bp_)**2) / &
                    sqrt(mu0__CGS*Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%rho_))

            dl_min=min(Block1%dxi,Block1%dxj*minval(Block1%xi_I),Block1%dxk*minval(Block1%xi_I))
            dt = CFL_ad * dl_min / maxval(Block1%v_wave_III(1:ni,1:nj,1:nk))
        end select
    end subroutine ModTimeStep_Corona_single
end module ModTimeStep
