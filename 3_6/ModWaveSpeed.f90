module ModWaveSpeed

    use ModBlock,      only: BlockType
    use ModParameters, only: ni, nj, nk

    contains

    subroutine ModWaveSpeed_Dynamo(Block1)
        implicit none
        type(BlockType), intent(inout) :: Block1

        Block1%v_wave_III(1:ni,1:nj,1:nk) = &
            1./Block1%Xi_rsst_III(1:ni,1:nj,1:nk) * &
            sqrt(Block1%gamma1_III(1:ni,1:nj,1:nk) * Block1%p0_over_rho0_III(1:ni,1:nj,1:nk)) + &
            sqrt(Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%vr_)**2 + &
                 Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%vt_)**2 + &
                 Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%vp_)**2)
        
        if (Block1%if_involve_B) then
            Block1%v_wave_III(1:ni,1:nj,1:nk) = Block1%v_wave_III(1:ni,1:nj,1:nk) + &
                sqrt(Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%br_)**2 + &
                     Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%bt_)**2 + &
                     Block1%primitive_IV(1:ni,1:nj,1:nk,Block1%bp_)**2) / &
                sqrt(Block1%rho0_III(1:ni,1:nj,1:nk))
        end if

    end subroutine ModWaveSpeed_Dynamo
end module ModWaveSpeed
