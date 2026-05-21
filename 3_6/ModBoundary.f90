module ModBoundary

    use ModBlock,       only:   BlockType
    use ModParameters,  only:   ni,nj,nk,ng

    contains

    ! Set the primitives in ghost cell for Dynamo HD

    subroutine ModBoundary_Dynamo_HD_primitives(Block1)
        implicit none
        type(BlockType),target          ::  Block1
        integer                         ::  i

        ! Set boundary condition
        if (Block1%if_top) then
            do i=ni+1,ni+ng
                ! mirror boundary condition for rho1, vt, vp, s1
                ! zero boundary condition for vr
                Block1%primitive(i,:,:,Block1%rho1_)=Block1%primitive(2*ni+1-i,:,:,Block1%rho1_)
                Block1%primitive(i,:,:,Block1%vt_)=Block1%primitive(2*ni+1-i,:,:,Block1%vt_)
                Block1%primitive(i,:,:,Block1%vp_)=Block1%primitive(2*ni+1-i,:,:,Block1%vp_)
                Block1%primitive(i,:,:,Block1%s1_)=Block1%primitive(2*ni+1-i,:,:,Block1%s1_)
                Block1%primitive(i,:,:,Block1%vr_)=-Block1%primitive(2*ni+1-i,:,:,Block1%vr_)
            end do
        end if
        if (Block1%if_bottom) then
            do i=-ng+1,0
                Block1%primitive(i,:,:,Block1%rho1_)=Block1%primitive(1-i,:,:,Block1%rho1_)
                Block1%primitive(i,:,:,Block1%vt_)=Block1%primitive(1-i,:,:,Block1%vt_)
                Block1%primitive(i,:,:,Block1%vp_)=Block1%primitive(1-i,:,:,Block1%vp_)
                Block1%primitive(i,:,:,Block1%s1_)=Block1%primitive(1-i,:,:,Block1%s1_)
                Block1%primitive(i,:,:,Block1%vr_)=-Block1%primitive(1-i,:,:,Block1%vr_)
            end do
        end if
    end subroutine ModBoundary_Dynamo_HD_primitives

    subroutine ModBoundary_Dynamo_MHD_primitives(Block1)
        implicit none
        type(BlockType),target          ::  Block1
        integer                         ::  i

        ! Set boundary condition
        if (Block1%if_top) then

            ! For the upper boundary, let's start from
            ! 1. rho1 s1 and vrtp the same with HD case
            ! 2. br is mirror, while bt and bp are negative mirror.
            do i=ni+1,ni+ng
                Block1%primitive(i,:,:,Block1%rho1_)=Block1%primitive(2*ni+1-i,:,:,Block1%rho1_)
                Block1%primitive(i,:,:,Block1%vt_)=Block1%primitive(2*ni+1-i,:,:,Block1%vt_)
                Block1%primitive(i,:,:,Block1%vp_)=Block1%primitive(2*ni+1-i,:,:,Block1%vp_)
                Block1%primitive(i,:,:,Block1%s1_)=Block1%primitive(2*ni+1-i,:,:,Block1%s1_)
                Block1%primitive(i,:,:,Block1%vr_)=-Block1%primitive(2*ni+1-i,:,:,Block1%vr_)
                Block1%primitive(i,:,:,Block1%br_)=Block1%primitive(2*ni+1-i,:,:,Block1%br_)
                Block1%primitive(i,:,:,Block1%bt_)=-Block1%primitive(2*ni+1-i,:,:,Block1%bt_)
                Block1%primitive(i,:,:,Block1%bp_)=-Block1%primitive(2*ni+1-i,:,:,Block1%bp_)
            end do
        end if
        if (Block1%if_bottom) then

            ! For the inner boundary, the difference is
            ! that br is negative mirror while bt and bp are positive mirror.
            ! which is like a perfect conductor boundary condition.
            do i=-ng+1,0
                Block1%primitive(i,:,:,Block1%rho1_)=Block1%primitive(1-i,:,:,Block1%rho1_)
                Block1%primitive(i,:,:,Block1%vt_)=Block1%primitive(1-i,:,:,Block1%vt_)
                Block1%primitive(i,:,:,Block1%vp_)=Block1%primitive(1-i,:,:,Block1%vp_)
                Block1%primitive(i,:,:,Block1%s1_)=Block1%primitive(1-i,:,:,Block1%s1_)
                Block1%primitive(i,:,:,Block1%vr_)=-Block1%primitive(1-i,:,:,Block1%vr_)
                Block1%primitive(i,:,:,Block1%br_)=-Block1%primitive(1-i,:,:,Block1%br_)
                Block1%primitive(i,:,:,Block1%bt_)=Block1%primitive(1-i,:,:,Block1%bt_)
                Block1%primitive(i,:,:,Block1%bp_)=Block1%primitive(1-i,:,:,Block1%bp_)
            end do
        end if
    end subroutine ModBoundary_Dynamo_MHD_primitives

end module ModBoundary