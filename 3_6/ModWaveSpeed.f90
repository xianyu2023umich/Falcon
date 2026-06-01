module ModWaveSpeed

    use ModConst,      only: dpi
    use ModBlock,      only: BlockType
    use ModParameters, only: ni, nj, nk, iEquation

    contains

    subroutine ModWaveSpeed_v_Alfven_Dynamo(Block1)
        implicit none
        type(BlockType), intent(inout) :: Block1

        Block1%v_wave_III = Block1%v_wave_III + &
                sqrt(Block1%primitive_IV(:,:,:,Block1%br_)**2 + &
                     Block1%primitive_IV(:,:,:,Block1%bt_)**2 + &
                     Block1%primitive_IV(:,:,:,Block1%bp_)**2) / &
                sqrt(Block1%rho0_III(:,:,:)*4.d0*dpi)
    end subroutine ModWaveSpeed_v_Alfven_Dynamo

    subroutine ModWaveSpeed_Dynamo(Block1)
        implicit none
        type(BlockType), intent(inout) :: Block1

        Block1%v_wave_III = Block1%v_sound_III + &
            sqrt(Block1%primitive_IV(:,:,:,Block1%vr_)**2 + &
                 Block1%primitive_IV(:,:,:,Block1%vt_)**2 + &
                 Block1%primitive_IV(:,:,:,Block1%vp_)**2)
        
        if (iEquation==1 .and. Block1%if_involve_B) then
            call ModWaveSpeed_v_Alfven_Dynamo(Block1)
            Block1%v_wave_III = &
                Block1%v_wave_III + &
                Block1%v_alfven_III
        end if

    end subroutine ModWaveSpeed_Dynamo
end module ModWaveSpeed
