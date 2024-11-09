module ModBoundary

    use ModOcTree

    contains

    ! Determine whether the block is in the top or bottom (or both)
    !
    subroutine ModBundary_if_top_bottom(Block1,xijk_range,if_top,if_bottom)

        implicit none

        type(Block),intent(in)  :: Block1               ! The input Block
        real,intent(in)         :: xijk_range(3,2)      ! xijk_range of whole tree
        logical,intent(out)     :: if_top,if_bottom     ! if top or bottom
        

        if (Block1%xk(Block1%nk+Block1%ng) .gt. xijk_range(3,2)) then
            if_top=.True.
        else
            if_top=.False.
        end if
        if (Block1%xk(-Block1%ng+1) .lt. xijk_range(3,1)) then
            if_bottom=.True.
        else
            if_bottom=.False.
        end if
    end subroutine ModBundary_if_top_bottom

end module ModBoundary