module ModBoundary

    use ModBlock,       only:   BlockType
    use ModParameters,  only:   ni,nj,nk,ng
    use ModVariables,   only:   rho1_,vr_,vt_,vp_,br_,bt_,bp_,s1_

    contains

    ! Set the primitives in ghost cell for Dynamo HD

    subroutine ModBoundary_Dynamo_HD_primitives(Block1,if_rk)
        implicit none
        type(BlockType),target          ::  Block1
        logical,intent(in)              ::  if_rk
        real,pointer                    ::  primitive(:,:,:,:)
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
                primitive(i,:,:,[rho1_,vt_,vp_,s1_])=primitive(2*ni+1-i,:,:,[rho1_,vt_,vp_,s1_])
                primitive(i,:,:,vr_)=-primitive(2*ni+1-i,:,:,vr_)
            end do
        end if
        if (Block1%if_bottom) then
            do i=-ng+1,0
                primitive(i,:,:,[rho1_,vt_,vp_,s1_])=primitive(1-i,:,:,[rho1_,vt_,vp_,s1_])
                primitive(i,:,:,vr_)=-primitive(1-i,:,:,vr_)
            end do
        end if
    end subroutine ModBoundary_Dynamo_HD_primitives

    subroutine ModBoundary_Dynamo_MHD_primitives(Block1,if_rk)
        implicit none
        type(BlockType),target          ::  Block1
        logical,intent(in)              ::  if_rk
        real,pointer                    ::  primitive(:,:,:,:)
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
                primitive(i,:,:,[rho1_,vt_,vp_,s1_,bt_,bp_])=primitive(2*ni+1-i,:,:,[rho1_,vt_,vp_,s1_,bt_,bp_])
                primitive(i,:,:,[vr_,br_])=-primitive(2*ni+1-i,:,:,[vr_,br_])
            end do
        end if
        if (Block1%if_bottom) then
            do i=-ng+1,0
                primitive(i,:,:,[rho1_,vt_,vp_,s1_,bt_,bp_])=primitive(1-i,:,:,[rho1_,vt_,vp_,s1_,bt_,bp_])
                primitive(i,:,:,[vr_,br_])=-primitive(1-i,:,:,[vr_,br_])
            end do
        end if
    end subroutine ModBoundary_Dynamo_MHD_primitives

end module ModBoundary