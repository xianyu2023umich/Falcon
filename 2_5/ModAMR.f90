module ModAMR

    use ModMath,        only:   ModMath_IfLinesInterSect
    use ModYinYangTree, only:   YYTree,TreeNode,&
                                YinYangTree_Divide

    integer             ::      AMR_nLevels,AMR_iLevel
    real,allocatable    ::      AMR_r_ranges(:,:)
    logical,allocatable ::      AMR_rtp_if_divide(:,:)

    contains

    ! The main subroutine to set the grid for the yinyang Tree.
    subroutine ModAMR_set_grid(Tree)
        implicit none
        type(YYTree),target     ::  Tree
        
        ! Loop AMR_nLevels times
        ! for appropriate AMR_nLevels values
        if (AMR_nLevels>0 .and. AMR_nLevels<=16) then
            
            do AMR_iLevel=1,AMR_nLevels
                call ModAMR_DivideAll(Tree,[AMR_rtp_if_divide(AMR_iLevel,:)])
            end do
        end if
    end subroutine ModAMR_set_grid

    ! Perform division for one AMR_iLevel
    subroutine ModAMR_DivideAll(Tree,rtp_if_divide)
        implicit none
        type(YYTree)                ::  Tree
        logical,intent(in)          ::  rtp_if_divide(3)

        call ModAMR_Divide_OneBranch(Tree,Tree%Yin,rtp_if_divide)
        call ModAMR_Divide_OneBranch(Tree,Tree%Yang,rtp_if_divide)
    end subroutine ModAMR_DivideAll

    ! Perform division for one AMR_iLevel on one branch
    recursive subroutine ModAMR_Divide_OneBranch(Tree,Node,rtp_if_divide)
        implicit none
        type(YYTree)                ::  Tree
        type(TreeNode),target       ::  Node
        logical,intent(in)          ::  rtp_if_divide(3)

        integer                     ::  iChild

        if (Node%if_leaf) then
            if (ModMath_IfLinesInterSect(Node%rtp_range(1,:),AMR_r_ranges(:,AMR_iLevel))) &
                call YinYangTree_Divide(Tree,Node,rtp_if_divide)
        else
            do iChild=1,size(Node%Children)
                call ModAMR_Divide_OneBranch(Tree,&
                    Node%Children(iChild),rtp_if_divide)
            end do
        end if
    end subroutine ModAMR_Divide_OneBranch
end module ModAMR