module ModYinYangTree

    use ModConst, only : dpi
    use ModBoundary
    use ModYinYang
    use ModAllocation

    ! one node
    type TreeNode
        logical                     ::  if_leaf                 ! if this is the leaf node.
        logical                     ::  if_yin                  ! If yin
        real                        ::  rtp_range(3,2)          ! the coordinate ranges for three directions.
        integer                     ::  iLeafNode               ! which leaf node it is.
        type(TreeNode),pointer      ::  children(:)             ! the children of the node.
    end type TreeNode

    type YYTree
        type(TreeNode)              ::  Yin,Yang                ! Yin & Yang Root

        integer                     ::  nLocalBlocks            ! number of local blocks
        integer,allocatable         ::  iBlocksGlobal(:)        ! list for local iblocks
        type(Block),allocatable     ::  LocalBlocks(:)          ! the local blocks of this rank
        integer,allocatable         ::  iLeafNode_ranges(:,:)   ! the iNode table for each rank

        integer                     ::  NumLeafNodes            ! total number of leaf nodes     
        integer                     ::  NumLeafNodes_YinYang(2) ! n of leaf nodes in YinYang respectively 
        logical                     ::  if_labelled             ! if labelled

        real                        ::  r_range(2)              ! r range of tree
    end type YYTree

    contains

    ! Initiate the YYtree
    ! Compared with the common trees, there's no need to give
    ! rtp_range since YYTree covers the whole globe.
    ! However, one DOES need to tell the r range.
    ! Set unlabel, and initialize Yin & Yang as two roots

    subroutine YinYangTree_InitTree(Tree,r_range)
        implicit none
        type(YYTree)                ::  Tree                ! the generated tree
        real,intent(in)             ::  r_range(2)          ! radial range
        real                        ::  rtp_range_root(3,2) ! rtp_range for yin/yang root

        ! The Tree is initially unlabelld
        Tree%if_labelled=.False.

        ! Set r_range for the tree
        ! Get the rtp_range for Yin/Yang root
        Tree%r_range=r_range
        rtp_range_root(1,:)=r_range
        rtp_range_root(2,:)=[dpi*0.25,dpi*0.75]
        rtp_range_root(3,:)=[-dpi*0.75,dpi*0.75]

        ! Initialize the Yin/Yang root
        call YinYangTree_InitNode(Tree%Yin, rtp_range_root,.true. )
        call YinYangTree_InitNode(Tree%Yang,rtp_range_root,.false.)
    end subroutine YinYangTree_InitTree

    ! initiate the node
    ! which simply set the range and if_leaf for the node.

    subroutine YinYangTree_InitNode(Node,rtp_range,if_yin)
        implicit none
        type(TreeNode),target       ::  Node                ! The node to initialize
        real,intent(in)             ::  rtp_range(3,2)      ! r th ph range
        logical,intent(in)          ::  if_yin

        ! assign everything
        Node%if_leaf    =   .True.
        Node%rtp_range  =   rtp_range
        Node%if_yin     =   if_yin
    end subroutine YinYangTree_InitNode

    ! Divide all leaf nodes

    subroutine YinYangTree_DivideAll(Tree,rtp_if_divide)
        implicit none
        type(YYTree)                ::  Tree
        logical,intent(in)          ::  rtp_if_divide(3)

        call YinYangTree_Divide_OneBranch(Tree,Tree%Yin,rtp_if_divide)
        call YinYangTree_Divide_OneBranch(Tree,Tree%Yang,rtp_if_divide)
    end subroutine YinYangTree_DivideAll

    recursive subroutine YinYangTree_Divide_OneBranch(Tree,Node,rtp_if_divide)
        implicit none
        type(YYTree)                ::  Tree
        type(TreeNode),target       ::  Node
        logical,intent(in)          ::  rtp_if_divide(3)

        integer                     ::  iChild

        if (Node%if_leaf) then
            call YinYangTree_Divide(Tree,Node,rtp_if_divide)
        else
            do iChild=1,size(Node%Children)
                call YinYangTree_Divide_OneBranch(Tree,&
                    Node%Children(iChild),rtp_if_divide)
            end do
        end if
    end subroutine YinYangTree_Divide_OneBranch

    ! Divide one node into 2/4/8 children

    subroutine YinYangTree_Divide(Tree,Node,rtp_if_divide)
        implicit none
        type(YYTree)                ::  Tree                ! Tree
        type(TreeNode),target       ::  Node                ! The node to divide
        logical,intent(in)          ::  rtp_if_divide(3)    ! if to divide for rtp

        integer                     ::  iChild              ! i of children
        integer                     ::  ir,it,ip            ! indice
        integer                     ::  N_rtp(3)            ! N of childern for rtp
        real                        ::  rtp_range(3,2)      ! rtp range of Node
        real                        ::  rtp_range1(3,2)     ! rtp range of child
        
        ! get N_rtp
        N_rtp=merge(2,1,rtp_if_divide)
        
        ! Only do for the leaf node
        if (Node%if_leaf) then

            ! Every time a new divide is performed,
            ! reset if_labelled to false.
            Tree%if_labelled=.false.

            ! 8 children.
            allocate(Node%children(N_rtp(1)*N_rtp(2)*N_rtp(3)))
            iChild=1

            ! use rtp+range=Node%rtp_range to simply expressions
            rtp_range=Node%rtp_range

            do ir=1,N_rtp(1); do it=1,N_rtp(2); do ip=1,N_rtp(3)

                ! get the rtp_range1 for the child node

                rtp_range1(1,:)=&
                    rtp_range(1,1)  ** ([N_rtp(1)-ir+1.0, N_rtp(1)-ir+0.0]/N_rtp(1))/ &
                    rtp_range(1,2)  ** ([        -ir+1.0,         -ir+0.0]/N_rtp(1))

                rtp_range1(2,:)=&
                    rtp_range(2,1)  *   [N_rtp(2)-it+1.0, N_rtp(2)-it+0.0]/N_rtp(2) - &
                    rtp_range(2,2)  *   [        -it+1.0,         -it+0.0]/N_rtp(2)

                rtp_range1(3,:)=&
                    rtp_range(3,1)  *   [N_rtp(3)-ip+1.0, N_rtp(3)-ip+0.0]/N_rtp(3) - &
                    rtp_range(3,2)  *   [        -ip+1.0,         -ip+0.0]/N_rtp(3)

                ! initialize this Node
                ! the children should have the same if_yin as their mom
                call YinYangTree_InitNode(Node%children(iChild),rtp_range1,Node%if_yin)

                iChild=iChild+1
        end do; end do; end do; end if

        ! No longer leaf
        Node%if_leaf=.False.
    end subroutine YinYangTree_Divide

    ! Set the YYTree
    ! First Label the nodes and get total number. Do yin then Yang.
    ! Then get the local blocks and set & allocate them

    subroutine YinYangTree_SetAll(Tree,MpiSize,MpiRank)
        implicit none
        type(YYTree)                ::  Tree                ! Tree
        integer,intent(in)          ::  MpiSize,MpiRank     ! MPI

        integer,allocatable         ::  iLocalLeafNodes(:)  ! local nodes
        real,allocatable            ::  rtp_ranges(:,:,:)   ! of local blocks
        logical                     ::  if_yin
        integer                     ::  iLocalBlock,&       ! Local Block index
                                        iLocalBlock_Current
        
        ! initialize NumLeafNodes to be 0  
        Tree%NumLeafNodes=0
        Tree%NumLeafNodes_YinYang=0

        ! First label Yin and then Yang
        call YinYangTree_LabelNodes(Tree,Tree%Yin)
        Tree%NumLeafNodes_YinYang(1)=Tree%NumLeafNodes
        call YinYangTree_LabelNodes(Tree,Tree%Yang)
        Tree%NumLeafNodes_YinYang(2)=&
            Tree%NumLeafNodes-Tree%NumLeafNodes_YinYang(1)
        ! Get nodes table for each rank
        ! and local nodes
        !
        call ModAllocation_GetTable(Tree%NumLeafNodes,Tree%iLeafNode_ranges,MpiSize)
        call ModAllocation_GetNodes(Tree%NumLeafNodes,iLocalLeafNodes,MpiSize,MpiRank)
        ! Allocate blocks and block info
        !
        allocate(rtp_ranges(3,2,size(iLocalLeafNodes)))
        allocate(Tree%LocalBlocks(size(iLocalLeafNodes)))
        allocate(Tree%iBlocksGlobal(size(iLocalLeafNodes)))
        Tree%iBlocksGlobal=iLocalLeafNodes
        Tree%nLocalBlocks=size(iLocalLeafNodes)
        ! Find Blocks info
        !
        iLocalBlock_Current=1
        call YinYangTree_GetLocalBlocksInfo(Tree,Tree%Yin ,rtp_ranges,iLocalBlock_Current)
        call YinYangTree_GetLocalBlocksInfo(Tree,Tree%Yang,rtp_ranges,iLocalBlock_Current)
        ! At last initiate the blocks
        !
        do iLocalBlock=1,size(iLocalLeafNodes)
            ! decide if this block is yin or yang
            if_yin=(Tree%iBlocksGlobal(iLocalBlock).le.Tree%NumLeafNodes_YinYang(1))
            ! Initiate the Block.
            ! See if it's at the boundary (top or bottom or both)
            call ModBlock_Init(Tree%LocalBlocks(iLocalBlock),Tree%iBlocksGlobal(iLocalBlock),&
            rtp_ranges(:,:,iLocalBlock),if_yin=if_yin,if_ssm=.true.,if_use_actual_nvar=.true.)
            call ModBoundary_if_top_bottom(Tree%LocalBlocks(iLocalBlock),Tree%r_range)
        end do    
    end subroutine YinYangTree_SetAll

    ! Get the blocks info ( the xijk_ranges, specifically)
    ! from the Tree
    ! Get either from Yin or Yang

    recursive subroutine YinYangTree_GetLocalBlocksInfo(Tree,Node,xijk_ranges,iBlockLocal_Current)
        implicit none
        type(YYTree)                :: Tree
        type(TreeNode),target       :: Node
        real,intent(inout)          :: xijk_ranges(3,2,Tree%nLocalBlocks)
        integer,intent(inout)       :: iBlockLocal_Current

        integer :: iChild

        if (Node%if_leaf) then
            ! if it is leaf node then label it.
            ! num plus 1.

            if (Node%iLeafNode==Tree%iBlocksGlobal(iBlockLocal_Current)) then
                xijk_ranges(:,:,iBlockLocal_Current)=Node%rtp_range
                if (iBlockLocal_Current<size(Tree%iBlocksGlobal)) &
                    iBlockLocal_Current=iBlockLocal_Current+1
            end if
        else
            ! if not leaf node then go into children.

            do ichild=1,size(Node%children)
                call YinYangTree_GetLocalBlocksInfo(Tree,Node%children(ichild),xijk_ranges,iBlockLocal_Current)
            end do
        end if
    end subroutine YinYangTree_GetLocalBlocksInfo

    ! Label all the nodes from 1 to num

    recursive subroutine YInYangTree_LabelNodes(Tree,Node)
        implicit none
        type(YYTree)                :: Tree   ! the tree
        type(TreeNode),target       :: Node   ! the current node to deal with
        integer                     :: ichild
        
        if (Node%if_leaf) then
            ! if it is leaf node then label it.
            ! num plus 1.

            Tree%NumLeafNodes=Tree%NumLeafNodes+1
            Node%iLeafNode=Tree%NumLeafNodes

        else
            do ichild=1,size(Node%children)
                call YinYangTree_LabelNodes(Tree,Node%children(ichild))
            end do
        end if
    end subroutine YinYangTree_LabelNodes

    ! see if the 

    function YinYangTree_IfInsideNode(Node,rtp,if_yin) result(if_inside)
        implicit none
        type(TreeNode),intent(in)   ::  Node                ! the node
        real,intent(in)             ::  rtp(3)              ! the point
        real                        ::  rtp_here(3)         ! rtp for yy of node
        logical,intent(in)          ::  if_yin              ! if the point is yin
        logical                     ::  if_inside           ! output

        ! If the if_yin of point is different from that
        ! of the Node, convert it.
        rtp_here=merge(rtp,ModYinYang_CoordConv_0D(rtp),if_yin .eqv. Node%if_yin)

        ! See if inside
        if ((rtp_here(1)-Node%rtp_range(1,1))*(rtp_here(1)-Node%rtp_range(1,2)) .le. 0. .and. &
            (rtp_here(2)-Node%rtp_range(2,1))*(rtp_here(2)-Node%rtp_range(2,2)) .le. 0. .and. &
            (rtp_here(3)-Node%rtp_range(3,1))*(rtp_here(3)-Node%rtp_range(3,2)) .le. 0.) then
            
            if_inside=.True.
        else
            if_inside=.False.
        end if
    end function YinYangTree_IfInsideNode

    ! Find where the point is

    function YinYangTree_FindNode(Tree,rtp,if_yin) result(Node1)
        implicit none
        type(YYTree)                ::  Tree                ! The tree
        real,intent(in)             ::  rtp(3)              ! rtp of point
        logical,intent(in)          ::  if_yin              ! if_yin of point
        type(TreeNode),pointer      ::  Node1               ! result

        ! Initialize it to be Null
        Node1 => Null()

        if (YinYangTree_IfInsideNode(Tree%Yin,rtp,if_yin)) then
            Node1=>YinYangTree_FindNode_OneBranch(Tree%Yin,rtp,if_yin)
        else
            Node1=>YinYangTree_FindNode_OneBranch(Tree%Yang,rtp,if_yin)
        end if
    end function YinYangTree_FindNode

    ! Find where the point is for a single Yin OR Yang branch.

    recursive function YinYangTree_FindNode_OneBranch(Node0,rtp,if_yin) result(Node1)
        implicit none
        type(TreeNode),target       ::  Node0               ! Current Node
        real,intent(in)             ::  rtp(3)              ! Point coord
        real                        ::  rtp_here(3)         ! rtp here
        logical,intent(in)          ::  if_yin              ! if yin of point
        integer                     ::  NumChildren,iChild  ! about Children
        type(TreeNode),pointer      ::  Node1,Child         ! Pointers

        ! Initialize
        Node1 => Null()

        ! If the if_yin of point is different from that
        ! of the Node, convert it.
        rtp_here=merge(rtp,ModYinYang_CoordConv_0D(rtp),if_yin .eqv. Node0%if_yin)

        ! First, see if the point is inside this node

        if (YinYangTree_IfInsideNode(Node0,rtp,if_yin)) then

            ! If this node is leaf then assign it
            ! Else, loop the children

            if (Node0%if_leaf .eqv. .True.) then
                Node1 => Node0
            else
                NumChildren=size(Node0%Children)
                do iChild=1,NumChildren
                    Child=>Node0%Children(iChild)
                    if (YinYangTree_IfInsideNode(Child,rtp,if_yin)) Node1 => &
                        YinYangTree_FindNode_OneBranch(Child,rtp,if_yin)
                    
                    ! If find it then exit
                    if (associated(Node1)) exit
                end do
            end if
        end if
    end function YinYangTree_FindNode_OneBranch

    ! Get the rtp_range of a given iNode

    function YinYangTree_Get_rtp_range(Tree,iLeafNode) result(rtp_range)
        implicit none
        type(YYTree)                ::  Tree                ! The tree
        integer,intent(in)          ::  iLeafNode           ! The iLeafNode to Get
        real                        ::  rtp_range(3,2)      ! The output

        rtp_range=0.
        if (iLeafNode .le. Tree%NumLeafNodes_YinYang(1)) then
            call YinYangTree_Get_rtp_range_OneBranch(Tree%Yin,iLeafNode,rtp_range)
        else
            call YinYangTree_Get_rtp_range_OneBranch(Tree%Yang,iLeafNode,rtp_range)
        end if
    end function YinYangTree_Get_rtp_range

    recursive subroutine YinYangTree_Get_rtp_range_OneBranch(Node0,iLeafNode,rtp_range)
        implicit none
        type(TreeNode),intent(in)   ::  Node0
        integer,intent(in)          ::  iLeafNode
        real,intent(inout)          ::  rtp_range(3,2)
        integer                     ::  iChild

        if (Node0%if_leaf) then
            if (Node0%iLeafNode==iLeafNode)rtp_range=Node0%rtp_range
        else
            do iChild=1,size(Node0%Children)
                call YinYangTree_Get_rtp_range_OneBranch(Node0%Children(iChild),iLeafNode,rtp_range)
                if (maxval(rtp_range)>0.0 .or. minval(rtp_range)<0.0) exit
            end do
        end if
    end subroutine YinYangTree_Get_rtp_range_OneBranch


end module ModYinYangTree