Module ModOcTree

    use ModAllocation
    use ModBlock

    ! one node
    type OcTreeNode
        logical :: if_leaf          ! if this is the leaf node.
        real    :: xijk_range(3,2)  ! the coordinate ranges for three directions.
        integer :: iLeafNode        ! which leaf node it is.
        integer :: iNode            ! which node it is, including non-leaf nodes.

        type(OcTreeNode),pointer :: mom
        type(OcTreeNode),pointer :: children(:)! the children of the node.
    end type OcTreeNode

    ! the whole tree.
    ! starts from the root node and contains some blocks
    type OcTree
        type(OcTreeNode)        :: RootNode         ! the Root Node
        
        integer                 :: nLocalBlocks     ! number of local blocks
        integer,allocatable     :: iBlocksGlobal(:) ! list for local iblocks
        type(Block),allocatable :: LocalBlocks(:)   ! the local blocks of this rank
        integer,allocatable     :: iLeafNode_ranges(:,:)! the iNode table for each rank

        integer :: NumLeafNodes            ! total number of leaf nodes     
        integer :: NumAllNodes             ! total number of nodes including non-leaf ones. 
        real    :: xijk_range(3,2)         ! the coordinate ranges for three directions
        logical :: if_labelled             ! if labelled
        logical :: ijk_if_periodic(3) = &  ! if each direction is periodic
            [.True.,.True.,.False.]
    end type OcTree

    contains

    ! initiate the tree.

    subroutine OcTree_InitTree(Tree,xijk_range,ijk_if_periodic)
        implicit none
        type(OcTree)                :: Tree
        real,intent(in)             :: xijk_range(3,2)
        logical,intent(in),optional :: ijk_if_periodic(3)

        ! the tree is not labelled at first.
        Tree%if_labelled=.False.

        ! periodic
        if (present(ijk_if_periodic)) Tree%ijk_if_periodic=ijk_if_periodic

        ! ranges
        Tree%xijk_range=xijk_range

        ! initiate the rot node of the tree
        call OcTree_InitNode(Tree%RootNode,xijk_range)
    end subroutine OcTree_InitTree

    ! initiate the node
    ! which simply set the range and if_leaf for the node.

    subroutine OcTree_InitNode(Node,xijk_range)
        implicit none
        type(OcTreeNode),target :: Node
        real,intent(in)         :: xijk_range(3,2)

        Node%if_leaf=.True.
        Node%xijk_range=xijk_range
    end subroutine OcTree_InitNode

    ! divide one node into 2/4/8 children

    subroutine OcTree_Divide(Tree,Node,ijk_if_divide)
        implicit none
        type(OcTree)             :: Tree
        type(OcTreeNode),target  :: Node
        type(OcTreeNode),pointer :: Child
        logical,intent(in)       :: ijk_if_divide(3)

        integer :: i,j,k,iChild
        integer :: N_ijk(3)

        real :: xijk_range(3,2)
        real :: xijk_range_child(3,2)
        real :: dxi,dxj,dxk

        Tree%if_labelled=.false.

        xijk_range=Node%xijk_range

        do i=1,3
            if (ijk_if_divide(i)) then; N_ijk(i)=2; else; N_ijk(i)=1; end if;
        end do

        iChild=1
        allocate(Node%children(N_ijk(1)*N_ijk(2)*N_ijk(3)))
        
        if (Node%if_leaf .eqv. .True. .and. N_ijk(1)*N_ijk(2)*N_ijk(3) .gt. 1) then
            do i=1,N_ijk(1)
                dxi=(xijk_range(1,2)-xijk_range(1,1))/N_ijk(1)
                xijk_range_child(1,:)=xijk_range(1,1)+[0.,dxi]+dxi*(i-1)
                do j=1,N_ijk(2)
                    dxj=(xijk_range(2,2)-xijk_range(2,1))/N_ijk(2)
                    xijk_range_child(2,:)=xijk_range(2,1)+[0.,dxj]+dxj*(j-1)
                    do k=1,N_ijk(3)
                        dxk=(xijk_range(3,2)-xijk_range(3,1))/N_ijk(3)
                        xijk_range_child(3,:)=xijk_range(3,1)+[0.,dxk]+dxk*(k-1)

                        child => Node%children(iChild)
                        call OcTree_InitNode(child,xijk_range_child)
                        child%mom=>Node

                        iChild=iChild+1
        end do; end do; end do; end if

        Node%if_leaf=.False.

    end subroutine OcTree_Divide

    ! this subroutine should be run after dividing all the blocks
    ! first label the nodes determine the total number
    ! then determine the blocks for the local rank and allocate blocks
    ! then find the info of these blocks and initiate them

    subroutine OcTree_SetAll(Tree,ni,nj,nk,ng,NameEqn,rk_order,MpiSize,MpiRank)
        implicit none
        type(OcTree),intent(inout)  :: Tree                  ! Tree
        integer,intent(in)          :: ni,nj,nk,ng           ! grid
        integer,intent(in)          :: rk_order              ! runge kutta order
        character(len=*)            :: NameEqn               ! Name of equation
        integer,intent(in)          :: MpiSize,MpiRank       ! mpi info

        integer,allocatable         :: iLocalLeafNodes(:)             ! local nodes
        real,allocatable            :: xijk_ranges(:,:,:)

        integer                     :: iLocalBlock_Current,iLocalBlock

        iLocalBlock_Current=1

        ! Label nodes and get total number
        !
        Tree%NumLeafNodes=0
        Tree%NumAllNodes =0
        call OcTree_LabelNodes(Tree,Tree%RootNode)
        ! Get nodes table for each rank
        ! and local nodes
        !
        call ModAllocation_GetTable(Tree%NumLeafNodes,Tree%iLeafNode_ranges,MpiSize)
        call ModAllocation_GetNodes(Tree%NumLeafNodes,iLocalLeafNodes,MpiSize,MpiRank)
        ! Allocate blocks and block info
        !
        allocate(xijk_ranges(3,2,size(iLocalLeafNodes)))
        allocate(Tree%LocalBlocks(size(iLocalLeafNodes)))
        allocate(Tree%iBlocksGlobal(size(iLocalLeafNodes)))
        Tree%iBlocksGlobal=iLocalLeafNodes
        Tree%nLocalBlocks=size(iLocalLeafNodes)
        ! Find Blocks info
        !
        call OcTree_GetLocalBlocksInfo(Tree,Tree%RootNode,xijk_ranges,iLocalBlock_Current)
        ! At last initiate the blocks
        !
        do iLocalBlock=1,size(iLocalLeafNodes)
            call ModBlock_Init(Tree%LocalBlocks(iLocalBlock),Tree%iBlocksGlobal(iLocalBlock),&
                xijk_ranges(:,:,iLocalBlock),&
                ni,nj,nk,ng,NameEqn,rk_order)
        end do
    end subroutine OcTree_SetAll

    ! get the blocks info ( the xijk_ranges, specifically)
    ! from the Tree

    recursive subroutine OcTree_GetLocalBlocksInfo(Tree,Node,xijk_ranges,iBlockLocal_Current)
        implicit none
        type(OcTree) :: Tree
        type(OcTreeNode),intent(in),target :: Node
        real,intent(inout) :: xijk_ranges(3,2,Tree%nLocalBlocks)
        integer,intent(inout) :: iBlockLocal_Current

        integer :: iChild

        if (Node%if_leaf) then
            ! if it is leaf node then label it.
            ! num plus 1.

            if (Node%iLeafNode==Tree%iBlocksGlobal(iBlockLocal_Current)) then
                xijk_ranges(:,:,iBlockLocal_Current)=Node%xijk_range
                iBlockLocal_Current=iBlockLocal_Current+1
            end if
        else
            ! if not leaf node then go into children.

            do ichild=1,size(Node%children)
                call OcTree_GetLocalBlocksInfo(Tree,Node%children(ichild),xijk_ranges,iBlockLocal_Current)
            end do
        end if
    end subroutine OcTree_GetLocalBlocksInfo

    ! Label all the nodes from 1 to num

    recursive subroutine OcTree_LabelNodes(Tree,Node)
        implicit none
        type(OcTree)                          :: Tree   ! the tree
        type(OcTreeNode),intent(inout),target :: Node   ! the current node to deal with

        integer :: ichild
        
        if (Node%if_leaf) then
            ! if it is leaf node then label it.
            ! num plus 1.

            Tree%NumLeafNodes=Tree%NumLeafNodes+1
            Node%iLeafNode=Tree%NumLeafNodes

            ! for all node do the same thing

            Tree%NumAllNodes=Tree%NumAllNodes+1
            Node%iNode=Tree%NumAllNodes
        else
            ! if not leaf node then do nothing to Num(i)LeafNode
            ! but still need to do all node.

            Tree%NumAllNodes=Tree%NumAllNodes+1
            Node%iNode=Tree%NumAllNodes

            do ichild=1,size(Node%children)
                call OcTree_LabelNodes(Tree,Node%children(ichild))
            end do
        end if
    end subroutine OcTree_LabelNodes

    function OcTree_IfInsideNode(Node,xijk) result(if_inside)
        implicit none
        type(OcTreeNode),intent(in) :: Node
        real,intent(in) :: xijk(3)

        logical :: if_inside

        if ((xijk(1)-Node%xijk_range(1,1))*(xijk(1)-Node%xijk_range(1,2)) .le. 0. .and. &
            (xijk(2)-Node%xijk_range(2,1))*(xijk(2)-Node%xijk_range(2,2)) .le. 0. .and. &
            (xijk(3)-Node%xijk_range(3,1))*(xijk(3)-Node%xijk_range(3,2)) .le. 0.) then
            if_inside=.True.
        else
            if_inside=.False.
        end if
    end function OcTree_IfInsideNode

    recursive function OcTree_FindNode(Node0,ijk) result(Node1)
        implicit none
        type(OcTreeNode),intent(in),target :: Node0
        real,intent(in) :: ijk(3)

        integer :: NumChildren,iChild
        type(OcTreeNode),pointer :: child
        type(OcTreeNode),pointer :: Node1

        Node1 => Null()

        if (OcTree_IfInsideNode(Node0,ijk)) then
            if (Node0%if_leaf .eqv. .True.) then
                Node1 => Node0
            else
                NumChildren=size(Node0%children)
                do iChild=1,NumChildren
                    child => Node0%children(iChild)
                    if (OcTree_IfInsideNode(Child,ijk)) Node1 => OcTree_FindNode(Child,ijk)
                    if (associated(Node1)) exit
                end do
            end if
        end if
    end function OcTree_FindNode

end Module ModOcTree