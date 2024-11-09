module ModBoundary

    use ModBlock

    contains

    ! Determine whether the block is in the top or bottom (or both)
    !
    subroutine ModBoundary_if_top_bottom(Block1,r_range)
        implicit none
        type(Block),intent(inout)   ::  Block1              ! The input Block
        real,intent(in)             ::  r_range(2)          ! r_range of whole tree

        if (Block1%xi(Block1%ni+Block1%ng) .gt. r_range(2)) then
            Block1%if_top=.True.
        else
            Block1%if_top=.False.
        end if
        if (Block1%xi(-Block1%ng+1) .lt. r_range(1)) then
            Block1%if_bottom=.True.
        else
            Block1%if_bottom=.False.
        end if
    end subroutine ModBoundary_if_top_bottom

    ! Set the primitives in ghost cell for Dynamo HD

    subroutine ModBoundary_Dynamo_HD_primitives(Block1)
        implicit none
        type(Block),intent(inout)   ::  Block1
        integer                     ::  i

        if (Block1%if_top) then
            do i=Block1%ni+1,Block1%ni+Block1%ng
                Block1%primitive(i,:,:,[1,3,4])=(Block1%primitive(2*Block1%ni+1-i,:,:,[1,3,4]))
                Block1%primitive(i,:,:,2)=-Block1%primitive(2*Block1%ni+1-i,:,:,2)
                Block1%primitive(i,:,:,5)=-Block1%primitive(2*Block1%ni+1-i,:,:,5)
            end do
        end if
        if (Block1%if_bottom) then
            do i=-Block1%ng+1,0
                Block1%primitive(i,:,:,[1,3,4])=(Block1%primitive(1-i,:,:,[1,3,4]))
                Block1%primitive(i,:,:,2)=-Block1%primitive(1-i,:,:,2)
                Block1%primitive(i,:,:,5)=-Block1%primitive(1-i,:,:,5)
            end do
        end if
    end subroutine ModBoundary_Dynamo_HD_primitives
    
    subroutine ModBoundary_Dynamo_HD_p1(Block1,EQN_update_R)
        implicit none
        type(Block),intent(inout)   ::  Block1
        real,intent(in)             ::  EQN_update_R(:,:,:,:)
        integer                     ::  i

        if (Block1%if_top) then
            do i=Block1%ni+1,Block1%ni+Block1%ng
                Block1%p1(i,1:Block1%nj,1:Block1%nk)=&
                    (EQN_update_R(Block1%ni,:,:,4)*24.*Block1%dxi*&
                    Block1%rho0(Block1%ni,1:Block1%nj,1:Block1%nk)+&
                    27.*Block1%p1(Block1%ni-1,1:Block1%nj,1:Block1%nk)-&
                    Block1%p1(Block1%ni-2,1:Block1%nj,1:Block1%nk))/26.
            end do
        end if

        if (Block1%if_bottom) then
            do i=-Block1%ng+1,0
                Block1%p1(i,1:Block1%nj,1:Block1%nk)=&
                    (EQN_update_R(1,:,:,4)*24.*Block1%dxi*&
                    Block1%rho0(1,1:Block1%nj,1:Block1%nk)+&
                    27.*Block1%p1(2,1:Block1%nj,1:Block1%nk)-&
                    Block1%p1(3,1:Block1%nj,1:Block1%nk))/26.
            end do
        end if

    end subroutine ModBoundary_Dynamo_HD_p1

end module ModBoundary