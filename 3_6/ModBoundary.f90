module ModBoundary

    use ModBlock,       only:   BlockType
    use ModParameters,  only:   ni,nj,nk,ng

    contains

    ! Set the primitives in ghost cell for Dynamo HD

    subroutine ModBoundary_Dynamo_HD_primitives(Block1,if_rk)
        implicit none
        type(BlockType),target          ::  Block1
        logical,intent(in)              ::  if_rk
        real(8),pointer                 ::  primitive(:,:,:,:)
        integer                         ::  i

        ! Assign the pointer
        if (if_rk) then
            primitive=>Block1%primitive_rk_IV
        else
            primitive=>Block1%primitive_IV
        end if

        ! Set boundary condition
        if (Block1%if_top) then
            do i=ni+1,ni+ng
                ! mirror boundary condition for rho1, vt, vp, s1
                ! zero boundary condition for vr
                primitive(i,:,:,Block1%rho1_)=primitive(2*ni+1-i,:,:,Block1%rho1_)
                primitive(i,:,:,Block1%vt_)=primitive(2*ni+1-i,:,:,Block1%vt_)
                primitive(i,:,:,Block1%vp_)=primitive(2*ni+1-i,:,:,Block1%vp_)
                primitive(i,:,:,Block1%s1_)=primitive(2*ni+1-i,:,:,Block1%s1_)
                primitive(i,:,:,Block1%vr_)=-primitive(2*ni+1-i,:,:,Block1%vr_)
            end do
        end if
        if (Block1%if_bottom) then
            do i=-ng+1,0
                primitive(i,:,:,Block1%rho1_)=primitive(1-i,:,:,Block1%rho1_)
                primitive(i,:,:,Block1%vt_)=primitive(1-i,:,:,Block1%vt_)
                primitive(i,:,:,Block1%vp_)=primitive(1-i,:,:,Block1%vp_)
                primitive(i,:,:,Block1%s1_)=primitive(1-i,:,:,Block1%s1_)
                primitive(i,:,:,Block1%vr_)=-primitive(1-i,:,:,Block1%vr_)
            end do
        end if
    end subroutine ModBoundary_Dynamo_HD_primitives

    subroutine ModBoundary_Dynamo_MHD_primitives(Block1,if_rk)
        implicit none
        type(BlockType),target          ::  Block1
        logical,intent(in)              ::  if_rk
        real(8),pointer                 ::  primitive(:,:,:,:)
        integer                         ::  i

        ! Assign the pointer
        if (if_rk) then
            primitive=>Block1%primitive_rk_IV
        else
            primitive=>Block1%primitive_IV
        end if

        ! Set boundary condition
        if (Block1%if_top) then
            do i=ni+1,ni+ng
                primitive(i,:,:,Block1%rho1_)=primitive(2*ni+1-i,:,:,Block1%rho1_)
                primitive(i,:,:,Block1%vt_)=primitive(2*ni+1-i,:,:,Block1%vt_)
                primitive(i,:,:,Block1%vp_)=primitive(2*ni+1-i,:,:,Block1%vp_)
                primitive(i,:,:,Block1%s1_)=primitive(2*ni+1-i,:,:,Block1%s1_)
                primitive(i,:,:,Block1%vr_)=-primitive(2*ni+1-i,:,:,Block1%vr_)
                primitive(i,:,:,Block1%br_)=primitive(2*ni+1-i,:,:,Block1%br_)
                primitive(i,:,:,Block1%bt_)=primitive(2*ni+1-i,:,:,Block1%bt_)
                primitive(i,:,:,Block1%bp_)=primitive(2*ni+1-i,:,:,Block1%bp_)
            end do
        end if
        if (Block1%if_bottom) then
            do i=-ng+1,0
                primitive(i,:,:,Block1%rho1_)=primitive(1-i,:,:,Block1%rho1_)
                primitive(i,:,:,Block1%vt_)=primitive(1-i,:,:,Block1%vt_)
                primitive(i,:,:,Block1%vp_)=primitive(1-i,:,:,Block1%vp_)
                primitive(i,:,:,Block1%s1_)=primitive(1-i,:,:,Block1%s1_)
                primitive(i,:,:,Block1%vr_)=-primitive(1-i,:,:,Block1%vr_)
                primitive(i,:,:,Block1%br_)=primitive(1-i,:,:,Block1%br_)
                primitive(i,:,:,Block1%bt_)=primitive(1-i,:,:,Block1%bt_)
                primitive(i,:,:,Block1%bp_)=primitive(1-i,:,:,Block1%bp_)
            end do
        end if
    end subroutine ModBoundary_Dynamo_MHD_primitives

end module ModBoundary