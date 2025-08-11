module ModAMR

    use ModMath,        only:   ModMath_IfLinesInterSect
    use ModParameters,  only:   rLevelInitial,nAMRs,AMRs,AMRType
    use ModYinYangTree, only:   YYTree,TreeNode,&
                                YinYangTree_Divide,&
                                YinYangTree_Divide_r

    contains

    ! The main subroutine to set the grid for the yinyang Tree.
    subroutine ModAMR_set_grid(Tree)
        implicit none
        type(YYTree),target     ::  Tree
        integer                 ::  iAMR

        ! the first step is to set r grid.
        call ModAMR_Divide_r_grid(Tree)
        
        ! Loop nAMRs times
        ! for appropriate nAMRs values
        if (nAMRs>0 .and. nAMRs<=16) then
            
            do iAMR=1,nAMRs
                call ModAMR_Divide_OneBranch(Tree,Tree%Yin,iAMR)
                call ModAMR_Divide_OneBranch(Tree,Tree%Yang,iAMR)
            end do
        end if
    end subroutine ModAMR_set_grid

    subroutine ModAMR_Divide_r_grid(Tree)
        implicit none
        type(YYTree)                ::  Tree
        integer                     ::  rLevel

        if (rLevelInitial==0) then
        else if (rLevelInitial>10) then
            write(*,*)'R level too high.'
            stop 1
        else
            do rLevel=1,rLevelInitial
                call ModAMR_Divide_r_grid_OneBranch(Tree,Tree%Yin)
                call ModAMR_Divide_r_grid_OneBranch(Tree,Tree%Yang)
            end do
        end if
    end subroutine ModAMR_Divide_r_grid

    recursive subroutine ModAMR_Divide_r_grid_OneBranch(Tree,Node)
        implicit none
        type(YYTree)                ::  Tree
        type(TreeNode),target       ::  Node
        integer                     ::  iChild

        if (Node%if_leaf) then
            call YinYangTree_Divide_r(Tree,Node,i_option=1)
        else
            do iChild=1,size(Node%Children)
                call ModAMR_Divide_r_grid_OneBranch(Tree,Node%Children(iChild))
            end do
        end if
    end subroutine ModAMR_Divide_r_grid_OneBranch

    ! Perform division for one AMR_iLevel on one branch
    recursive subroutine ModAMR_Divide_OneBranch(Tree,Node,iAMR)
        implicit none
        type(YYTree)                ::  Tree
        type(TreeNode),target       ::  Node
        integer,intent(in)          ::  iAMR
        type(AMRType),pointer       ::  AMR1

        integer                     ::  iChild

        AMR1=>AMRs(iAMR)

        if (Node%if_leaf) then
            ! If it's global, then we only need to check if the 
            ! block's r_range intersects with the AMR range.
            ! If it's not global, then the AMR is only done for blocks that:
            ! 1. If_yin=.true.
            ! 2. rtp_range intersects with the AMR rtp_range.
            if (AMR1%if_global) then
                if (ModMath_IfLinesInterSect(Node%rtp_range(1,:),AMR1%rtp_range(1,:))) &
                    call YinYangTree_Divide(Tree,Node,AMR1%rtp_if_divide(1:3))
            else
                if (Node%if_yin) then
                    if (ModMath_IfLinesInterSect(Node%rtp_range(1,:),AMR1%rtp_range(1,:)) .and. &
                        ModMath_IfLinesInterSect(Node%rtp_range(2,:),AMR1%rtp_range(2,:)) .and. &
                        ModMath_IfLinesInterSect(Node%rtp_range(3,:),AMR1%rtp_range(3,:))) &
                        call YinYangTree_Divide(Tree,Node,AMR1%rtp_if_divide(1:3))
                end if
            end if
        else
            do iChild=1,size(Node%Children)
                call ModAMR_Divide_OneBranch(Tree,&
                    Node%Children(iChild),iAMR)
            end do
        end if
    end subroutine ModAMR_Divide_OneBranch
end module ModAMR