module ModCheck

    use ieee_arithmetic
    use ModBlock,       only:   BlockType
    use ModYinYangTree, only:   YYTree
    use ModParameters,  only:   ni,nj,nk,ng

    contains

    subroutine ModCheck_primitive(Tree,MpiRank)
        implicit none
        type(YYTree),target     ::  Tree
        type(BlockType),pointer ::  Block1
        integer,intent(in)      ::  MpiRank
        integer                 ::  iLocalBlock
        integer                 ::  ir,it,ip

        ! Loop all the local blocks
        do iLocalBlock=1,Tree%nLocalBlocks
            Block1=>Tree%LocalBlocks(iLocalBlock)

            ! Loop the primitives
            do ip=-ng+1,ng+nk; do it=-ng+1,ng+nj; do ir=-ng+1,ng+ni
                if(ieee_is_nan(Block1%primitive(ir,it,ip,1))) then
                    write(*,*)'Error: Detecting NAN at MpiRank=',MpiRank,&
                        'iBlock=',Block1%iBlock,'ir,it,ip=',ir,it,ip,&
                        'r,t,p=',Block1%xi(ir),Block1%xj(it),Block1%xk(ip)
                    stop 1
                end if
            end do; end do; end do
        end do
    end subroutine ModCheck_primitive
end module ModCheck