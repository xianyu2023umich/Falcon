module ModAMR

    use ModMath,        only:   ModMath_IfLinesInterSect
    use ModParameters,  only:   rLevelInitial
    use ModYinYangTree, only:   YYTree,TreeNode,&
                                YinYangTree_Divide,&
                                YinYangTree_Divide_r

    integer             ::      AMR_nLevels,AMR_iLevel
    real,allocatable    ::      AMR_r_ranges(:,:)
    logical,allocatable ::      AMR_rtp_if_divide(:,:)

    contains

    ! The main subroutine to set the grid for the yinyang Tree.
    subroutine ModAMR_set_grid(Tree)
        implicit none
        type(YYTree),target     ::  Tree

        ! the first step is to set r grid.
        call ModAMR_Divide_r_grid(Tree)
        
        ! Loop AMR_nLevels times
        ! for appropriate AMR_nLevels values
        if (AMR_nLevels>0 .and. AMR_nLevels<=16) then
            
            do AMR_iLevel=1,AMR_nLevels
                call ModAMR_DivideAll(Tree,[AMR_rtp_if_divide(AMR_iLevel,:)])
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
            call YinYangTree_Divide_r(Tree,Node)
        else
            do iChild=1,size(Node%Children)
                call ModAMR_Divide_r_grid_OneBranch(Tree,Node%Children(iChild))
            end do
        end if
    end subroutine ModAMR_Divide_r_grid_OneBranch

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