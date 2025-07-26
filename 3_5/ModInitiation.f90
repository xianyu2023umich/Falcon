module ModInitiation

    use ModBlock,       only:   BlockType
    use ModYinYang,     only:   ModYinyang_CoordConv_0D,&
                                ModYinYang_VecConv_0D
    use ModYinYangTree, only:   YYTree
    use ModParameters,  only:   r_range,ni,nj,nk,ng,&
                                rSave,nthSavePlot,nphSavePlot
    use ModConst,       only:   dpi
    use ModVariables,   only:   vr_,vt_,vp_
    contains

    subroutine ModInitiation_harmonic(Tree)
        
        implicit none
        type(YYTree),target     ::  Tree
        real                    ::  vec(3)
        real                    ::  coord(3)
        type(BlockType),pointer ::  Block1
        integer                 ::  iLocalBlock,ir,it,ip

        do iLocalBlock=1,Tree%nLocalBlocks
            Block1=>Tree%LocalBlocks(iLocalBlock)
            do ip=-ng+1,ng+nk
                do it=-ng+1,ng+nj
                    do ir=-ng+1,ng+ni
                        Block1%primitive_IV(ir,it,ip,:)=1.000000
                        coord=[Block1%xi_I(ir),Block1%xj_I(it),Block1%xk_I(ip)]
                        if (.not. Block1%if_yin) coord=ModYinyang_CoordConv_0D(coord)
                        vec=sin(dpi*(coord(1)-r_range(1))/(r_range(2)-r_range(1)))*[1.,1.,1.]*1.e-2*&
                            sin(12.*coord(3))*sin(coord(2))**4
                        Block1%primitive_IV(ir,it,ip,vr_:vp_)=vec
                        if (.not. Block1%if_yin) then
                            Block1%primitive_IV(ir,it,ip,vr_:vp_)=&
                            ModYinYang_VecConv_0D(coord,vec)
                        end if
                    end do
                end do
            end do
            Block1%primitive_rk_IV=Block1%primitive_IV
        end do
    end subroutine ModInitiation_harmonic

end module ModInitiation
