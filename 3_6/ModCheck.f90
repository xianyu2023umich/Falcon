module ModCheck

    use ieee_arithmetic
    use ModBlock,       only:   BlockType
    use ModYinYangTree, only:   YYTree
    use ModParameters,  only:   ni,nj,nk,ng,nvar,MpiRank
    use ModYinYang,     only:   ModYinYang_CoordConv_0D
    use MPI

    contains

    subroutine ModCheck_primitive(Tree,if_rk)
        implicit none
        type(YYTree),target     ::  Tree
        logical,intent(in)      ::  if_rk
        type(BlockType),pointer ::  Block1
        integer                 ::  iLocalBlock
        integer                 ::  ir,it,ip,ivar
        real(8)                 ::  coord(3)
        integer                 ::  ierr

        if (if_rk) then
            ! Loop all the local blocks
            do iLocalBlock=1,Tree%nLocalBlocks
                Block1=>Tree%LocalBlocks(iLocalBlock)

                ! Loop the primitives
                do ip=-ng+1,ng+nk
                    do it=-ng+1,ng+nj
                        do ir=-ng+1,ng+ni
                            do ivar=1,nvar
                                if(ieee_is_nan(Block1%primitive_rk_IV(ir,it,ip,ivar))) then
                                    coord=[Block1%xi_I(ir),Block1%xj_I(it),Block1%xk_I(ip)]
                                    if (.not. Block1%if_yin) coord=ModYinYang_CoordConv_0D(coord)
                                    write(*,*)'Error: Detecting NAN at MpiRank=',MpiRank,&
                                        'iBlock=',Block1%iBlock,'ir,it,ip,ivar=',ir,it,ip,ivar,&
                                        'r,t,p=',coord
                                    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        else
            ! Loop all the local blocks
            do iLocalBlock=1,Tree%nLocalBlocks
                Block1=>Tree%LocalBlocks(iLocalBlock)

                ! Loop the primitives
                do ip=-ng+1,ng+nk
                    do it=-ng+1,ng+nj
                        do ir=-ng+1,ng+ni
                            do ivar=1,nvar
                                if(ieee_is_nan(Block1%primitive_IV(ir,it,ip,ivar))) then
                                    coord=[Block1%xi_I(ir),Block1%xj_I(it),Block1%xk_I(ip)]
                                    if (.not. Block1%if_yin) coord=ModYinYang_CoordConv_0D(coord)
                                    write(*,*)'Error: Detecting NAN at MpiRank=',MpiRank,&
                                        'iBlock=',Block1%iBlock,'ir,it,ip,ivar=',ir,it,ip,ivar,&
                                        'r,t,p=',coord
                                    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end if


    end subroutine ModCheck_primitive

    subroutine ModCheck_DiffCool_Power(Tree)
        implicit none
        type(YYTree),target     ::  Tree
        type(BlockType),pointer ::  Block1
        integer                 ::  iLocalBlock,ir,it,ip
        real(8)                 ::  P_diff_local,P_cool_local
        real(8)                 ::  P_diff_total,P_cool_total
        integer                 ::  ierr

        P_diff_local = 0.0d0
        P_cool_local = 0.0d0

        do iLocalBlock = 1, Tree%nLocalBlocks
            Block1 => Tree%LocalBlocks(iLocalBlock)
            do ip = 1, nk
                do it = 1, nj
                    do ir = 1, ni
                        P_diff_local = P_diff_local + Block1%diffusion_I(ir) * Block1%V_LLL(ir,it,ip)
                        P_cool_local = P_cool_local + Block1%cooling_I(ir)   * Block1%V_LLL(ir,it,ip)
                    end do
                end do
            end do
        end do

        call MPI_Reduce(P_diff_local, P_diff_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_Reduce(P_cool_local, P_cool_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        if (MpiRank==0) write(*,'(a,es13.5,a,es13.5)') &
            'Diffusion power [erg/s] =', P_diff_total, &
            '   Cooling power [erg/s] =', P_cool_total
    end subroutine ModCheck_DiffCool_Power

end module ModCheck