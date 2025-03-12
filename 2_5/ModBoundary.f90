module ModBoundary

    use ModBlock,       only:   BlockType
    use ModParameters,  only:   ni,nj,nk,ng
    use ModVariables,   only:   rho1_,vr_,vt_,vp_,br_,bt_,bp_,s1_

    contains

    ! Determine whether the block is in the top or bottom (or both)
    !
    subroutine ModBoundary_if_top_bottom(Block1,r_range)
        implicit none
        type(BlockType),intent(inout)   ::  Block1              ! The input Block
        real,intent(in)                 ::  r_range(2)          ! r_range of whole tree

        if (Block1%xi(ni+ng) .gt. r_range(2)) then
            Block1%if_top=.True.
        else
            Block1%if_top=.False.
        end if
        if (Block1%xi(-ng+1) .lt. r_range(1)) then
            Block1%if_bottom=.True.
        else
            Block1%if_bottom=.False.
        end if
    end subroutine ModBoundary_if_top_bottom

    ! Set the primitives in ghost cell for Dynamo HD

    subroutine ModBoundary_Dynamo_HD_primitives(Block1,if_rk)
        implicit none
        type(BlockType),target          ::  Block1
        logical,intent(in)              ::  if_rk
        real,pointer                    ::  primitive(:,:,:,:)
        integer                         ::  i

        ! Assign the pointer
        if (if_rk) then
            primitive=>Block1%primitive_rk
        else
            primitive=>Block1%primitive
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
            primitive=>Block1%primitive_rk
        else
            primitive=>Block1%primitive
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
    
    subroutine ModBoundary_Dynamo_HD_p1(Block1)
        implicit none
        type(BlockType),intent(inout)   ::  Block1
        integer                         ::  i

        if (Block1%if_top) then
            do i=ng+1,ng+ng
                Block1%p1(i,1:ng,1:ng)=&

                    (27.*   Block1%p1(ng-1,1:ng,1:ng)              &   ! Derivative term
                    -       Block1%p1(ng-2,1:ng,1:ng))/26.         &   

                    -       12.0/26.0*Block1%dxi                                        &   ! Gravity term 12.0=24.0/2.0
                    *       Block1%primitive  (rho1_,ng,  1:ng,1:ng)*  &   ! rho1 at the boundary
                           (Block1%g_over_rho0  (ng,  1:ng,1:ng)*  &   ! g at the boundary
                            Block1%rho0         (ng,  1:ng,1:ng)+  &
                            Block1%g_over_rho0  (ng+1,1:ng,1:ng)*  &
                            Block1%rho0         (ng+1,1:ng,1:ng))
            end do
        end if

        if (Block1%if_bottom) then
            do i=-ng+1,0
                Block1%p1(i,1:ng,1:ng)=&

                    (27.*   Block1%p1(2,1:ng,1:ng)                        &   ! Derivative term
                    -       Block1%p1(3,1:ng,1:ng))/26.                   &   

                    +       12.0/26.0*Block1%dxi                                        &   ! Gravity term 12.0=24.0/2.0
                    *       Block1%primitive  (rho1_,1,1:ng,1:ng)*            &   ! rho1 at the boundary
                           (Block1%g_over_rho0  (1,1:ng,1:ng)*            &   ! g at the boundary
                            Block1%rho0         (1,1:ng,1:ng)+            &
                            Block1%g_over_rho0  (0,1:ng,1:ng)*            &
                            Block1%rho0         (0,1:ng,1:ng))
            end do
        end if

    end subroutine ModBoundary_Dynamo_HD_p1

end module ModBoundary