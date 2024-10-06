Module ModGC

    use ModOcTree
    use ModMath,only : ModMath_1D3D_interpolate_1D1D,ModMath_count_unique_elements_3D,ModMath_IfBlocksInterSect
    use ModAllocation,only : ModAllocation_GetNodes,ModAllocation_GetRank

    contains

    ! communicate local GC_targets

    subroutine ModGC_CommunicateGCLocal(Tree,MpiRank,rk_index)
        implicit none
        type(OcTree),target             ::  Tree
        integer,intent(in)              ::  MpiRank,rk_index

        integer                         ::  iLocalBlock
        integer                         ::  iGC_Target,iGC
        type(Block),pointer             ::  Block1,Block2
        type(GC_target),pointer         ::  GC_target1
        real,pointer                    ::  primitive_send(:,:,:,:),&
                                            primitive_recv(:,:,:,:)

        real,allocatable :: primitive_GC(:,:)

        ! loop all the local blocks

        do iLocalBlock=1,size(Tree%LocalBlocks)

            Block1=>Tree%LocalBlocks(iLocalBlock)

            ! loop all the GC_targets
            !
            do iGC_Target=1,size(Block1%GC_targets)

                ! check if the current GC_target belongs to 
                ! the same rank

                if (Block1%GC_targets(iGC_Target)%iRank==MpiRank) then

                    ! assign the pointers
                    ! allocate primitive_GC

                    GC_target1=>Block1%GC_targets(iGC_Target)
                    Block2=>Tree%LocalBlocks(GC_target1%iBlock-Tree%iLeafNode_ranges(MpiRank,1)+1)
                    allocate(primitive_GC(Block1%nvar,GC_target1%nGC))

                    ! rk_index

                    select case(rk_index)
                    case(1) 
                        primitive_recv=>Block1%primitive
                        primitive_send=>Block2%primitive
                    case(2)
                        primitive_recv=>Block1%primitive_k2
                        primitive_send=>Block2%primitive_k2
                    case(3)
                        primitive_recv=>Block1%primitive_k3
                        primitive_send=>Block2%primitive_k3
                    case(4)
                        primitive_recv=>Block1%primitive_k4
                        primitive_send=>Block2%primitive_k4
                    end select

                    ! interpolate the block to primitive_GC

                    call ModMath_1D3D_interpolate_1D1D(primitive_send,Block2%nvar,&
                        Block2%ni,Block2%nj,Block2%nk,Block2%ng,&
                        Block2%xi,Block2%xj,Block2%xk,GC_target1%nGC,GC_target1%xijk_list,&
                        primitive_GC)

                    ! give it to the current block
    
                    do iGC=1,GC_target1%nGC
                        primitive_recv(:, GC_target1%ijk_list(iGC,1),&
                                        GC_target1%ijk_list(iGC,2),&
                                        GC_target1%ijk_list(iGC,3))=primitive_GC(:,iGC)
                    end do

                    deallocate(primitive_GC)
                end if
            end do
        end do

    end subroutine ModGC_CommunicateGCLocal

    ! Get the GC_sources for every local blocks of
    ! the tree.
    !
    subroutine ModGC_GetGC_SourcesAll(Tree,MpiSize)
        implicit none
        type(OcTree),intent(in),target  :: Tree
        type(Block),pointer             :: Block1
        integer,intent(in)              :: MpiSize

        integer                         :: iLocalBlock

        do iLocalBlock=1,Tree%nLocalBlocks
            Block1=>Tree%LocalBlocks(iLocalBlock)
            call ModGC_GetGC_Sources_Single(Tree,Block1,Tree%iBlocksGlobal(iLocalBlock),MpiSize)
        end do
    end subroutine ModGC_GetGC_SourcesAll

    ! Get GC_sources for one single block
    ! 
    ! In this subroutine you will see two definitions:
    !
    ! 1. neighbour: any block of which the xijk_range including its GCs
    !               intersects with the xijk_range of our particular block1
    !
    ! 2. GC neighbour: any block that has at least one GC in our block1
    !
    ! Obviously that {GC neighbour} is a subset of {neighbour}. In most of
    ! the cases they should be identical, but MIGHT be different when I 
    ! later go into the YinYang grid. So just a preparation for that.
    ! 
    subroutine ModGC_GetGC_Sources_Single(Tree,Block1,iBlockGlobal,MpiSize)
        implicit none
        type(OcTree),intent(in)         :: Tree
        type(Block),intent(inout),target:: Block1
        integer,intent(in)              :: iBlockGlobal
        integer,intent(in)              :: MpiSize
        
        integer                         :: NumNeighbours,NumGCNeighbours
        integer                         :: iNeighbour,iGCNeighbour,iGC_target
        type(GC_target),pointer         :: GC_target1,GC_source1
        real,allocatable                :: xijk_ranges_Neighbours(:,:,:)
        integer,allocatable             :: iBlocks_Neighbours(:)
        logical,allocatable             :: if_GC_neighbours(:)
        type(Block),target              :: Block2

        ! first count the number of neighbours.
        ! default it to be 0
        !
        
        NumNeighbours=0
        call ModGC_CountNeighbours(Tree,Tree%RootNode,Block1,NumNeighbours)

        ! allocate:
        ! 1. iBlocks_Neighbours: the iBlocksGlobal of all neighbours
        ! 2. xijk_ranges_neighbours: xijk_ranges of them
        ! 3. if_GC_neighbours: if each of them is GC neighbour
        !
        allocate(iBlocks_Neighbours(NumNeighbours))
        allocate(xijk_ranges_Neighbours(NumNeighbours,3,2))
        allocate(if_GC_neighbours(NumNeighbours)); if_GC_neighbours=.false.

        ! then, find the neighbours
        ! initialize iNeighbour to be 0
        !
        iNeighbour=0
        call ModGC_FindNeighbours(Tree,Tree%RootNode,Block1,NumNeighbours,iNeighbour,&
            iBlocks_Neighbours,xijk_ranges_Neighbours)

        ! count the number of GC_neighbours
        !
        NumGCNeighbours=0
        do iNeighbour=1,NumNeighbours
            ! initiate each block in turn to get its GC positions
            ! and then locate the target blocks of the GCs
            !
            call ModBlock_InitGrid(Block2,xijk_ranges_Neighbours(iNeighbour,:,:),&
                1,Block1%ni,Block1%nj,Block1%nk,Block1%ng,1)
            call ModGC_FindGC_Single(Tree,Block2)
            call ModGC_SetGC_Single(Tree,Block2,MpiSize)

            ! loop by all the unique GC target blocks to see whether
            ! at least one GC falls in block1.
            !
            do iGC_target=1,Block2%nGC_targets
                if (Block2%GC_targets(iGC_target)%iBlock==iBlockGlobal) then
                    NumGCNeighbours=NumGCNeighbours+1
                    exit
                end if
            end do
            
            ! at last deallocate
            call ModBlock_deallocate(Block2)
        end do

        ! GC_sources should have NumGCNeighbours elements
        !
        Block1%nGC_sources=NumGCNeighbours
        allocate(Block1%GC_sources(NumGCNeighbours))

        ! The GC_targets of the GC_neighbours targeting at block1
        ! are copied to block1's GC_sources
        !
        ! Initialize iGCNeighbour to 0
        !
        iGCNeighbour=0
        do iNeighbour=1,NumNeighbours
            ! Initiate the grid. Find GC positions and targets
            !
            call ModBlock_InitGrid(Block2,xijk_ranges_Neighbours(iNeighbour,:,:),&
                1,Block1%ni,Block1%nj,Block1%nk,Block1%ng,1)
            call ModGC_FindGC_Single(Tree,Block2)
            call ModGC_SetGC_Single(Tree,Block2,MpiSize)

            ! loop the GC targets to find those fall into block1
            ! and copy them to GC_sources of block1
            !
            do iGC_target=1,Block2%nGC_targets
                if (Block2%GC_targets(iGC_target)%iBlock==iBlockGlobal) then
                    ! the current index += 1
                    ! Set two pointers to simplify the expressions
                    !
                    iGCNeighbour=iGCNeighbour+1
                    GC_source1=>Block1%GC_sources(iGCNeighbour)
                    GC_target1=>Block2%GC_targets(iGC_target)

                    ! copy the iBlock of GC_target1 into GC_source1
                    ! then find iRank
                    !
                    GC_source1%iBlock=iBlocks_Neighbours(iNeighbour)
                    call ModAllocation_GetRank(Tree%NumLeafNodes,&
                        GC_source1%iBlock,MpiSize,GC_source1%iRank)
                    
                    
                    ! copy nGC then allocate and copy (x)ijk_list
                    !
                    GC_source1%nGC=GC_target1%nGC
                    allocate(GC_source1%ijk_list(GC_source1%nGC,3))
                    allocate(GC_source1%xijk_list(GC_source1%nGC,3))
                    allocate(GC_source1%primitive_list(Block1%nvar,GC_source1%nGC))
                    GC_source1%ijk_list=GC_target1%ijk_list
                    GC_source1%xijk_list=GC_target1%xijk_list
                    exit
                end if
            end do
            call ModBlock_deallocate(Block2)
        end do

        ! at last, deallocate the two lists/
        !
        deallocate(iBlocks_Neighbours)
        deallocate(xijk_ranges_Neighbours)
        deallocate(if_GC_neighbours)
    end subroutine ModGC_GetGC_Sources_Single

    ! Count the number of neighbours
    !
    ! The definition of neighbour:
    ! 
    ! If block1%xijk_range intersects with 
    ! block2%xijk_range_GC, then block2 is a
    ! neighbour of block1. 
    !
    recursive subroutine ModGC_CountNeighbours(Tree,CurrentNode,Block1,NumNeighbours)
        implicit none
        type(OcTree),intent(in)         :: Tree
        type(Block),intent(in)          :: Block1
        type(OcTreeNode),intent(in)     :: CurrentNode  ! this is NOT the Node of the Block1
                                                        ! BUT the Current Node
        integer,intent(inout)           :: NumNeighbours

        real                            :: xijk_range_GC_CurrentNode(3,2)
        integer                         :: iChild

        real :: dxi_CurrentNode,dxj_CurrentNode,dxk_CurrentNode

        if (CurrentNode%if_leaf) then
            dxi_CurrentNode=(CurrentNode%xijk_range(1,2)-CurrentNode%xijk_range(1,1))/Block1%ni
            dxj_CurrentNode=(CurrentNode%xijk_range(2,2)-CurrentNode%xijk_range(2,1))/Block1%nj
            dxk_CurrentNode=(CurrentNode%xijk_range(3,2)-CurrentNode%xijk_range(3,1))/Block1%nk

            xijk_range_GC_CurrentNode(1,1)=CurrentNode%xijk_range(1,1)-(Block1%ng-0.5)*dxi_CurrentNode
            xijk_range_GC_CurrentNode(2,1)=CurrentNode%xijk_range(2,1)-(Block1%ng-0.5)*dxj_CurrentNode
            xijk_range_GC_CurrentNode(3,1)=CurrentNode%xijk_range(3,1)-(Block1%ng-0.5)*dxk_CurrentNode

            xijk_range_GC_CurrentNode(1,2)=CurrentNode%xijk_range(1,2)+(Block1%ng-0.5)*dxi_CurrentNode
            xijk_range_GC_CurrentNode(2,2)=CurrentNode%xijk_range(2,2)+(Block1%ng-0.5)*dxj_CurrentNode
            xijk_range_GC_CurrentNode(3,2)=CurrentNode%xijk_range(3,2)+(Block1%ng-0.5)*dxk_CurrentNode
            if (ModMath_IfBlocksInterSect(Block1%xijk_range,xijk_range_GC_CurrentNode,&
                Tree%ijk_if_periodic,Tree%xijk_range)) NumNeighbours=NumNeighbours+1
        else
            ! if not leaf node then go into children.

            do ichild=1,size(CurrentNode%children)
                call ModGC_CountNeighbours(Tree,CurrentNode%Children(iChild),Block1,NumNeighbours)
            end do
        end if
    end subroutine ModGC_CountNeighbours

    ! find the neighbour blocks
    !
    recursive subroutine ModGC_FindNeighbours(Tree,CurrentNode,Block1,&
        NumNeighbours,iNeighbour,iBlocks_Neighbours,xijk_ranges_Neighbours)
        implicit none
        type(OcTree),intent(in)         :: Tree
        type(Block),intent(in)          :: Block1
        type(OcTreeNode),intent(in)     :: CurrentNode  ! this is NOT the Node of the Block1
                                                        ! BUT the Current Node
        integer,intent(in)              :: NumNeighbours
        integer,intent(inout)           :: iNeighbour
        integer,intent(inout)           :: iBlocks_Neighbours(NumNeighbours)
        real,intent(inout)              :: xijk_ranges_Neighbours(NumNeighbours,3,2)

        real                            :: xijk_range_GC_CurrentNode(3,2)
        integer                         :: iChild

        real :: dxi_CurrentNode,dxj_CurrentNode,dxk_CurrentNode

        if (CurrentNode%if_leaf) then
            dxi_CurrentNode=(CurrentNode%xijk_range(1,2)-CurrentNode%xijk_range(1,1))/Block1%ni
            dxj_CurrentNode=(CurrentNode%xijk_range(2,2)-CurrentNode%xijk_range(2,1))/Block1%nj
            dxk_CurrentNode=(CurrentNode%xijk_range(3,2)-CurrentNode%xijk_range(3,1))/Block1%nk

            xijk_range_GC_CurrentNode(1,1)=CurrentNode%xijk_range(1,1)-(Block1%ng-0.5)*dxi_CurrentNode
            xijk_range_GC_CurrentNode(2,1)=CurrentNode%xijk_range(2,1)-(Block1%ng-0.5)*dxj_CurrentNode
            xijk_range_GC_CurrentNode(3,1)=CurrentNode%xijk_range(3,1)-(Block1%ng-0.5)*dxk_CurrentNode

            xijk_range_GC_CurrentNode(1,2)=CurrentNode%xijk_range(1,2)+(Block1%ng-0.5)*dxi_CurrentNode
            xijk_range_GC_CurrentNode(2,2)=CurrentNode%xijk_range(2,2)+(Block1%ng-0.5)*dxj_CurrentNode
            xijk_range_GC_CurrentNode(3,2)=CurrentNode%xijk_range(3,2)+(Block1%ng-0.5)*dxk_CurrentNode
            if (ModMath_IfBlocksInterSect(Block1%xijk_range,xijk_range_GC_CurrentNode,&
                Tree%ijk_if_periodic,Tree%xijk_range)) then
                iNeighbour=iNeighbour+1
                
                iBlocks_Neighbours(iNeighbour)=CurrentNode%iLeafNode
                xijk_ranges_Neighbours(iNeighbour,:,:)=CurrentNode%xijk_range
            end if
        else
            ! if not leaf node then go into children.

            do ichild=1,size(CurrentNode%children)
                call ModGC_FindNeighbours(Tree,CurrentNode%Children(iChild),&
                    Block1,NumNeighbours,iNeighbour,iBlocks_Neighbours,&
                    xijk_ranges_Neighbours)
            end do
        end if
    end subroutine ModGC_FindNeighbours

    subroutine ModGC_InitGCAll(Tree,MpiSize)

        implicit none

        type(OcTree),intent(inout),target   :: Tree        ! the Tree
        integer     ,intent(in)             :: MpiSize

        call ModGC_FindGCAll(Tree)
        call ModGC_SetGCAll(Tree,MpiSize)
    end subroutine ModGC_InitGCAll

    ! find where each GC belongs to
    !
    subroutine ModGC_FindGCAll(Tree)

        implicit none

        type(OcTree),intent(inout),target:: Tree        ! the Tree
        type(Block),pointer              :: Block1      ! one Block

        integer :: iBlock

        do iBlock=1,size(Tree%LocalBlocks)

            Block1=>Tree%LocalBlocks(iBlock)
            call ModGC_FindGC_Single(Tree,Block1)
        end do
    end subroutine ModGC_FindGCAll

    subroutine ModGC_FindGC_Single(Tree,Block1)
        implicit none
        type(OcTree),intent(in)         :: Tree
        type(Block),intent(inout)       :: Block1      ! one Block
        type(OcTreeNode),pointer        :: Node1       ! one Node

        integer :: i,j,k,ni,nj,nk,ng
        real :: xijk(3)
        real :: xijk_range(3,2)
        logical :: ijk_if_periodic(3)

        Block1%GC_iBlocks=-1                        ! initial value -1.
        ijk_if_periodic=Tree%ijk_if_periodic
        xijk_range=Tree%RootNode%xijk_range

        ni=Block1%ni
        nj=Block1%nj
        nk=Block1%nk
        ng=Block1%ng
            
        do i=-ng+1,ni+ng; do j=-ng+1,nj+ng; do k=-ng+1,nk+ng

            xijk=[Block1%xi(i),Block1%xj(j),Block1%xk(k)]

            if (ijk_if_periodic(1)) then
                if (xijk(1) .lt. xijk_range(1,1)) xijk(1)=xijk(1)+(xijk_range(1,2)-xijk_range(1,1))
                if (xijk(1) .gt. xijk_range(1,2)) xijk(1)=xijk(1)-(xijk_range(1,2)-xijk_range(1,1))
            end if 

            if (ijk_if_periodic(2)) then
                if (xijk(2) .lt. xijk_range(2,1)) xijk(2)=xijk(2)+(xijk_range(2,2)-xijk_range(2,1))
                if (xijk(2) .gt. xijk_range(2,2)) xijk(2)=xijk(2)-(xijk_range(2,2)-xijk_range(2,1))
            end if 

            if (ijk_if_periodic(3)) then
                if (xijk(3) .lt. xijk_range(2,1)) xijk(3)=xijk(3)+(xijk_range(3,2)-xijk_range(3,1))
                if (xijk(3) .gt. xijk_range(3,2)) xijk(3)=xijk(3)-(xijk_range(3,2)-xijk_range(3,1))
            end if 

            if (.not. (i .ge. 1 .and. j .ge. 1 .and. k .ge. 1 .and. &
                i .le. ni .and. j .le. nj .and. k .le. nk)) then
                    
                if (OcTree_IfInsideNode(Tree%RootNode,xijk)) then
                    Node1=>OcTree_FindNode(Tree%RootNode,xijk)
                    Block1%GC_iBlocks(i,j,k)=Node1%iLeafNode
                end if

            end if
        end do; end do; end do
    end subroutine ModGC_FindGC_Single
    
    subroutine ModGC_SetGCAll(Tree,MpiSize)
        implicit none
        type(OcTree),target :: Tree
        type(Block),pointer :: Block1
        integer,intent(in)  :: MpiSize

        integer :: iBlockLocal

        do iBlockLocal=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iBlockLocal)
            call ModGC_SetGC_Single(Tree,Block1,MpiSize)
        end do
    end subroutine ModGC_SetGCAll

    ! set the GC_targets for a single block

    subroutine ModGC_SetGC_Single(Tree,Block1,MpiSize)
        implicit none
        type(OcTree)            ::  Tree                                 ! Tree
        type(Block),target      ::  Block1                               ! The block to set
        integer,intent(in)      ::  MpiSize                              ! MpiSize

        integer                 ::  i,j,k,iUniqueGCiBlock,&
                                    iUniqueGCiBlockGlobal_current,&
                                    ni,nj,nk,ng,iGC_current

        type(GC_target),pointer ::  GC_target1
        integer,allocatable     ::  UniqueGCiBlocksGlobal(:),nGC_list(:)
        integer                 ::  nUniqueGCiBlocksGlobal
        logical                 ::  ijk_if_periodic(3)
        real                    ::  xijk(3)
        real                    ::  xijk_range(3,2)

        ijk_if_periodic=Tree%ijk_if_periodic
        xijk_range=Tree%RootNode%xijk_range

        ni=Block1%ni
        nj=Block1%nj
        nk=Block1%nk
        ng=Block1%ng

        ! count how many unique blocks are there within the whole block containing GCs.
        ! obviously this should be the number of GC_targets + 1, since it includes the block itself.
        !
        UniqueGCiBlocksGlobal=ModMath_count_unique_elements_3D(Block1%GC_iBlocks,ni+2*ng,nj+2*ng,nk+2*ng,&
            nUniqueGCiBlocksGlobal,nGC_list)

        Block1%nGC_targets=nUniqueGCiBlocksGlobal-1
        allocate(Block1%GC_targets(Block1%nGC_targets))

        ! loop all the GC_targets
        ! initiate the current index to be 0
        !
        iUniqueGCiBlockGlobal_current=0

        do iUniqueGCiBlock=1,nUniqueGCiBlocksGlobal

            ! if this is not this block itself then begin
            !
            if (UniqueGCiBlocksGlobal(iUniqueGCiBlock) .ne. -1) then

                iUniqueGCiBlockGlobal_current=iUniqueGCiBlockGlobal_current+1

                GC_target1=>Block1%GC_targets(iUniqueGCiBlockGlobal_current)
                GC_target1%iBlock=UniqueGCiBlocksGlobal(iUniqueGCiBlock)

                ! see which rank the block of this GC_target belongs to
                !
                call ModAllocation_GetRank(Tree%NumLeafNodes,GC_target1%iBlock,&
                    MpiSize,GC_target1%iRank)
                
                ! allocate
                !
                GC_target1%nGC=nGC_list(iUniqueGCiBlock)
                allocate(GC_target1%xijk_list(nGC_list(iUniqueGCiBlock),3))
                allocate(GC_target1%ijk_list(nGC_list(iUniqueGCiBlock),3))

                ! loop all the points to fill in the GC_target1
                !
                iGC_current=0
                do i=-ng+1,ni+ng; do j=-ng+1,nj+ng; do k=-ng+1,nk+ng
                    if (Block1%GC_iBlocks(i,j,k) .eq. UniqueGCiBlocksGlobal(iUniqueGCiBlock)) then

                        iGC_current=iGC_current+1
                        GC_target1%ijk_list(iGC_current,:)=[i,j,k]
                        xijk=[Block1%xi(i),Block1%xj(j),Block1%xk(k)]
                        
                        if (ijk_if_periodic(1)) then
                            if (xijk(1) .lt. xijk_range(1,1)) xijk(1)=xijk(1)+(xijk_range(1,2)-xijk_range(1,1))
                            if (xijk(1) .gt. xijk_range(1,2)) xijk(1)=xijk(1)-(xijk_range(1,2)-xijk_range(1,1))
                        end if 
        
                        if (ijk_if_periodic(2)) then
                            if (xijk(2) .lt. xijk_range(2,1)) xijk(2)=xijk(2)+(xijk_range(2,2)-xijk_range(2,1))
                            if (xijk(2) .gt. xijk_range(2,2)) xijk(2)=xijk(2)-(xijk_range(2,2)-xijk_range(2,1))
                        end if 
        
                        if (ijk_if_periodic(3)) then
                            if (xijk(3) .lt. xijk_range(2,1)) xijk(3)=xijk(3)+(xijk_range(3,2)-xijk_range(3,1))
                            if (xijk(3) .gt. xijk_range(3,2)) xijk(3)=xijk(3)-(xijk_range(3,2)-xijk_range(3,1))
                        end if 

                        GC_target1%xijk_list(iGC_current,:)=xijk
                    end if
                end do; end do; end do
                
            end if
            
        end do
        
        deallocate(UniqueGCiBlocksGlobal)
        deallocate(nGC_list)

    end subroutine ModGC_SetGC_Single

end module ModGC