module ModInitiation

    use ModBlock,       only:   BlockType
    use ModYinYang,     only:   ModYinyang_CoordConv_0D,&
                                ModYinYang_VecConv_0D
    use ModYinYangTree, only:   YYTree
    use ModParameters,  only:   r_range,ni,nj,nk,ng,Initiation_type_index
    use ModConst,       only:   dpi
    contains

    subroutine ModInitiation_DoAll(Tree)
        implicit none
        type(YYTree),target     ::  Tree

        select case (Initiation_type_index)
        case (1)
            call ModInitiation_harmonic(Tree)
        case (2)
            call ModInitiation_all_ones(Tree)
        case default
    end select
    end subroutine

    subroutine ModInitiation_harmonic(Tree)
        
        implicit none
        type(YYTree),target     ::  Tree
        real(8)                 ::  vec(3)
        real(8)                 ::  coord(3)
        type(BlockType),pointer ::  Block1
        integer                 ::  iLocalBlock,ir,it,ip

        
        do iLocalBlock=1,Tree%nLocalBlocks
            
            Block1=>Tree%LocalBlocks(iLocalBlock)

            do ip=-ng+1,ng+nk
                do it=-ng+1,ng+nj
                    do ir=-ng+1,ng+ni
                        Block1%primitive_IV(ir,it,ip,:)=0.000000
                        coord=[Block1%xi_I(ir),Block1%xj_I(it),Block1%xk_I(ip)]
                        if (.not. Block1%if_yin) coord=ModYinyang_CoordConv_0D(coord)
                        vec=sin(dpi*(coord(1)-r_range(1))/(r_range(2)-r_range(1)))*[1.,1.,1.]*1.e-2*&
                            sin(12.*coord(3))*sin(coord(2))**4
                        Block1%primitive_IV(ir,it,ip,Block1%vr_:Block1%vp_)=vec
                        if (.not. Block1%if_yin) then
                            Block1%primitive_IV(ir,it,ip,Block1%vr_:Block1%vp_)=&
                            ModYinYang_VecConv_0D(coord,vec)
                        end if
                        
                    end do
                    
                end do
                
            end do
            Block1%primitive_rk_IV=Block1%primitive_IV
        end do
    end subroutine ModInitiation_harmonic

    subroutine ModInitiation_all_ones(Tree)
        implicit none
        type(YYTree),target     ::  Tree
        type(BlockType),pointer ::  Block1
        integer                 ::  iLocalBlock,ir,it,ip

        do iLocalBlock=1,Tree%nLocalBlocks
            Block1=>Tree%LocalBlocks(iLocalBlock)
            do ip=-ng+1,ng+nk
                do it=-ng+1,ng+nj
                    do ir=-ng+1,ng+ni
                        Block1%primitive_IV(ir,it,ip,:)=1.0d0
                    end do
                end do
            end do
            Block1%primitive_rk_IV=Block1%primitive_IV
        end do
    end subroutine ModInitiation_all_ones

    subroutine ModInitiation_rand_velocity(Tree, v_rms)
        implicit none
        type(YYTree), target   :: Tree
        real(8), intent(in)    :: v_rms

        type(BlockType), pointer :: Block1
        integer :: iLocalBlock, ir, it, ip
        real(8) :: coord(3), vec(3)

        do iLocalBlock = 1, Tree%nLocalBlocks
            Block1 => Tree%LocalBlocks(iLocalBlock)
            do ip = 1, nk
                do it = 1, nj
                    do ir = 1, ni
                        call random_number(vec)
                        vec = (vec - 0.5d0) * 2.0d0 * v_rms

                        coord = [Block1%xi_I(ir), Block1%xj_I(it), Block1%xk_I(ip)]
                        if (.not. Block1%if_yin) coord = ModYinyang_CoordConv_0D(coord)

                        Block1%primitive_IV(ir, it, ip, Block1%vr_:Block1%vp_) = vec
                        if (.not. Block1%if_yin) &
                            Block1%primitive_IV(ir, it, ip, Block1%vr_:Block1%vp_) = &
                                ModYinYang_VecConv_0D(coord, vec)
                    end do
                end do
            end do
            Block1%primitive_rk_IV = Block1%primitive_IV
        end do
    end subroutine ModInitiation_rand_velocity

end module ModInitiation
