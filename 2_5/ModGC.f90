Module ModGC

    use ModYinYangTree
    use ModMath
    use ModAllocation
    use ModParameters,  only:   ni,nj,nk,ng,nvar

    contains

    ! communicate local GC_targets

    subroutine ModGC_CommunicateGCLocal(Tree,MpiRank,if_rk)
        implicit none
        type(YYTree),target             ::  Tree
        integer,intent(in)              ::  MpiRank
        logical,intent(in)              ::  if_rk

        integer                         ::  iLocalBlock
        integer                         ::  iGC_Target,iGC
        type(Block),pointer             ::  Block_target,Block_source
        type(GC_target),pointer         ::  GC_target1
        real,pointer                    ::  primitive_send(:,:,:,:),&
                                            primitive_recv(:,:,:,:)

        real,allocatable :: primitive_GC(:,:)

        ! loop all the local blocks

        do iLocalBlock=1,size(Tree%LocalBlocks)
            Block_target=>Tree%LocalBlocks(iLocalBlock)

            ! loop all the GC_targets and check if
            ! the current GC_target belongs to the same rank

            do iGC_Target=1,size(Block_target%GC_targets)
                if (Block_target%GC_targets(iGC_Target)%iRank==MpiRank) then

                    ! assign the pointers
                    ! allocate primitive_GC

                    GC_target1=>Block_target%GC_targets(iGC_Target)
                    Block_source=>Tree%LocalBlocks(GC_target1%iBlock-Tree%iLeafNode_ranges(MpiRank,1)+1)
                    allocate(primitive_GC(GC_target1%nGC,nvar))

                    ! if_rk

                    if (if_rk) then
                        primitive_recv=>Block_target%primitive_rk
                        primitive_send=>Block_source%primitive_rk
                    else
                        primitive_recv=>Block_target%primitive
                        primitive_send=>Block_source%primitive
                    end if

                    ! Interpolate the block to primitive_GC
                    ! If different if_yin then convert coord
                    ! and then convert v vector

                    if (Block_target%if_yin .eqv. Block_source%if_yin) then
                        call ModMath_1D3D_interpolate_1D1D(primitive_send,nvar,ni,nj,nk,ng,&
                            Block_source%xi,Block_source%xj,Block_source%xk,GC_target1%nGC,&
                            GC_target1%xijk_list,primitive_GC)
                    else
                        call ModMath_1D3D_interpolate_1D1D(primitive_send,nvar,ni,nj,nk,ng,&
                            Block_source%xi,Block_source%xj,Block_source%xk,GC_target1%nGC,&
                            ModYinYang_CoordConv_1D(GC_target1%xijk_list,GC_target1%nGC),&
                            primitive_GC)
                        !print *,primitive_GC(1,2:4)
                        primitive_GC(:,2:4)=&
                            ModYinYang_VecConv_1D(ModYinYang_CoordConv_1D(GC_target1%xijk_list,GC_target1%nGC),&
                                primitive_GC(:,2:4),GC_target1%nGC)
                        !print *,primitive_GC(1,2:4)
                        !print *,11,Block_target%if_yin
                    end if

                    ! assign it to the current block
    
                    do iGC=1,GC_target1%nGC
                        primitive_recv( GC_target1%ijk_list(iGC,1),&
                                        GC_target1%ijk_list(iGC,2),&
                                        GC_target1%ijk_list(iGC,3),:)=primitive_GC(iGC,:)
                    end do

                    deallocate(primitive_GC)
                end if
            end do
        end do

    end subroutine ModGC_CommunicateGCLocal

    ! Find where each GC belongs to
    ! for all the local blocks
    subroutine ModGC_SetGCAll(Tree,MpiSize)
        implicit none
        type(YYTree),target             ::  Tree                    ! the Tree
        type(Block),pointer             ::  Block1                  ! one Block
        integer                         ::  MpiSize
        integer                         ::  iBlock

        ! Loop each block and call 
        ! ModGC_FindGC_Single
        do iBlock=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iBlock)
            call ModGC_FindGC_Single(Tree,Block1)
        end do

        print *,'Find GC done.'

        do iBlock=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iBlock)
            call ModGC_GetGC_Targets_single(Tree,Block1,MpiSize)
        end do

        print *,'GCTarget done.'

        do iBlock=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iBlock)
            call ModGC_GetGC_Sources_Single(Tree,Block1,MpiSize)
        end do

        print *,'GCSource done.'
    end subroutine ModGC_SetGCAll

    ! Do for only one block
    subroutine ModGC_FindGC_Single(Tree,Block1)
        implicit none
        type(YYTree),intent(in)         ::  Tree
        type(Block),intent(inout)       ::  Block1                  ! one Block
        type(TreeNode),pointer          ::  Node1                   ! one Node
        integer                         ::  i,j,k                   ! grid
        real                            ::  rtp(3)                  ! rtp of one point

        ! Initialize to be -1
        Block1%GC_iBlocks=-1                        
        
        ! Find the node of each GC
        do i=-ng+1,ni+ng; do j=-ng+1,nj+ng; do k=-ng+1,nk+ng
            rtp=[Block1%xi(i),Block1%xj(j),Block1%xk(k)]

            if (i<1 .or. j<1 .or. k<1 .or. i>ni .or. j>nj .or. k>nk) then
                Node1=>YinYangTree_FindNode(Tree,rtp,Block1%if_yin)
                if (associated(Node1))then
                    Block1%GC_iBlocks(i,j,k)=Node1%iLeafNode
                end if
                
            end if
        end do; end do; end do
    end subroutine ModGC_FindGC_Single

    ! set the GC_targets for a single block

    subroutine ModGC_GetGC_Targets_single(Tree,Block1,MpiSize)
        implicit none
        type(YYTree)                    ::  Tree                                 ! Tree
        type(Block),target              ::  Block1                               ! The block to set
        integer,intent(in)              ::  MpiSize                              ! MpiSize

        integer                         ::  i,j,k,iUniqueGCiBlock,&
                                            iUniqueGCiBlockGlobal_current,&
                                            iGC_current

        type(GC_target),pointer         ::  GC_target1
        integer,allocatable             ::  UniqueGCiBlocksGlobal(:),nGC_list(:)
        integer                         ::  nUniqueGCiBlocksGlobal

        ! Initializations to avoid warning message.
        nUniqueGCiBlocksGlobal=0

        ! count how many unique blocks are there within the whole block containing GCs.
        ! obviously this should be the number of GC_targets + 1, since it includes the block itself.
        UniqueGCiBlocksGlobal=ModMath_count_unique_elements_3D(Block1%GC_iBlocks,ni+2*ng,nj+2*ng,nk+2*ng,&
            nUniqueGCiBlocksGlobal,nGC_list)

        Block1%nGC_targets=nUniqueGCiBlocksGlobal-1
        allocate(Block1%GC_targets(Block1%nGC_targets))

        ! loop all the GC_targets
        ! initiate the current index to be 0
        iUniqueGCiBlockGlobal_current=0

        do iUniqueGCiBlock=1,nUniqueGCiBlocksGlobal

            ! if this is not this block itself then begin
            if (UniqueGCiBlocksGlobal(iUniqueGCiBlock) .ne. -1) then

                iUniqueGCiBlockGlobal_current=iUniqueGCiBlockGlobal_current+1

                GC_target1=>Block1%GC_targets(iUniqueGCiBlockGlobal_current)
                GC_target1%iBlock=UniqueGCiBlocksGlobal(iUniqueGCiBlock)
                GC_target1%if_yin=GC_target1%iBlock.le.Tree%NumLeafNodes_YinYang(1)

                ! see which rank the block of this GC_target belongs to
                call ModAllocation_GetRank(Tree%NumLeafNodes,GC_target1%iBlock,&
                    MpiSize,GC_target1%iRank)
                
                ! allocate
                GC_target1%nGC=nGC_list(iUniqueGCiBlock)
                allocate(GC_target1%xijk_list(nGC_list(iUniqueGCiBlock),3))
                allocate(GC_target1%ijk_list(nGC_list(iUniqueGCiBlock),3))

                ! loop all the points to fill in the GC_target1
                iGC_current=0
                do i=-ng+1,ni+ng; do j=-ng+1,nj+ng; do k=-ng+1,nk+ng
                    if (Block1%GC_iBlocks(i,j,k) .eq. UniqueGCiBlocksGlobal(iUniqueGCiBlock)) then

                        iGC_current=iGC_current+1
                        GC_target1%ijk_list (iGC_current,:)=[i,j,k]
                        GC_target1%xijk_list(iGC_current,:)=[Block1%xi(i),Block1%xj(j),Block1%xk(k)]
                    end if
                end do; end do; end do
            end if
        end do
        
        deallocate(UniqueGCiBlocksGlobal)
        deallocate(nGC_list)
    end subroutine ModGC_GetGC_Targets_single

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
        type(YYTree)    ,intent(in)     ::  Tree                            ! The Tree
        type(Block)     ,intent(in)     ::  Block1                          ! The Block  
        type(TreeNode)  ,intent(in)     ::  CurrentNode                     ! Current Node
        integer         ,intent(inout)  ::  NumNeighbours                   ! Num of neighbours

        real                            ::  rtp_range_GC_CurrentNode(3,2)
        integer                         ::  iChild
        real                            ::  dxi_CurrentNode,&
                                            dxj_CurrentNode,&
                                            dxk_CurrentNode

        if (CurrentNode%if_leaf) then

            dxi_CurrentNode=(CurrentNode%rtp_range(1,2)-CurrentNode%rtp_range(1,1))/ni
            dxj_CurrentNode=(CurrentNode%rtp_range(2,2)-CurrentNode%rtp_range(2,1))/nj
            dxk_CurrentNode=(CurrentNode%rtp_range(3,2)-CurrentNode%rtp_range(3,1))/nk

            ! The ranges CONSIDERING GCs
            rtp_range_GC_CurrentNode(1,1)=CurrentNode%rtp_range(1,1)-(ng-0.5)*dxi_CurrentNode
            rtp_range_GC_CurrentNode(1,2)=CurrentNode%rtp_range(1,2)+(ng-0.5)*dxi_CurrentNode
            rtp_range_GC_CurrentNode(2,1)=CurrentNode%rtp_range(2,1)-(ng-0.5)*dxj_CurrentNode
            rtp_range_GC_CurrentNode(2,2)=CurrentNode%rtp_range(2,2)+(ng-0.5)*dxj_CurrentNode
            rtp_range_GC_CurrentNode(3,1)=CurrentNode%rtp_range(3,1)-(ng-0.5)*dxk_CurrentNode 
            rtp_range_GC_CurrentNode(3,2)=CurrentNode%rtp_range(3,2)+(ng-0.5)*dxk_CurrentNode

            if (CurrentNode%if_yin.eqv.Block1%if_yin) then   
            else
                rtp_range_GC_CurrentNode=ModYinYang_GetOtherRange(rtp_range_GC_CurrentNode)
            end if

            ! See if they possibily intersect
            if (ModMath_IfBlocksInterSect(Block1%xijk_range,rtp_range_GC_CurrentNode)) &
                    NumNeighbours=NumNeighbours+1
            
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
        NumNeighbours,iNeighbour,iBlocks_Neighbours,rtp_ranges_Neighbours)
        implicit none
        type(YYTree),intent(in)         ::  Tree
        type(Block),intent(in)          ::  Block1
        type(TreeNode),intent(in)       ::  CurrentNode
        integer,intent(in)              ::  NumNeighbours
        integer,intent(inout)           ::  iNeighbour
        integer,intent(inout)           ::  iBlocks_Neighbours(NumNeighbours)
        real,intent(inout)              ::  rtp_ranges_Neighbours(NumNeighbours,3,2)

        real                            ::  rtp_range_GC_CurrentNode(3,2)
        integer                         ::  iChild

        real :: dxi_CurrentNode,dxj_CurrentNode,dxk_CurrentNode

        if (CurrentNode%if_leaf) then
            dxi_CurrentNode=(CurrentNode%rtp_range(1,2)-CurrentNode%rtp_range(1,1))/ni
            dxj_CurrentNode=(CurrentNode%rtp_range(2,2)-CurrentNode%rtp_range(2,1))/nj
            dxk_CurrentNode=(CurrentNode%rtp_range(3,2)-CurrentNode%rtp_range(3,1))/nk

            ! The ranges CONSIDERING GCs
            rtp_range_GC_CurrentNode(1,1)=CurrentNode%rtp_range(1,1)-(ng-0.5)*dxi_CurrentNode
            rtp_range_GC_CurrentNode(1,2)=CurrentNode%rtp_range(1,2)+(ng-0.5)*dxi_CurrentNode
            rtp_range_GC_CurrentNode(2,1)=CurrentNode%rtp_range(2,1)-(ng-0.5)*dxj_CurrentNode
            rtp_range_GC_CurrentNode(2,2)=CurrentNode%rtp_range(2,2)+(ng-0.5)*dxj_CurrentNode
            rtp_range_GC_CurrentNode(3,1)=CurrentNode%rtp_range(3,1)-(ng-0.5)*dxk_CurrentNode 
            rtp_range_GC_CurrentNode(3,2)=CurrentNode%rtp_range(3,2)+(ng-0.5)*dxk_CurrentNode

            if (CurrentNode%if_yin.eqv.Block1%if_yin) then
            else
                rtp_range_GC_CurrentNode=ModYinYang_GetOtherRange(rtp_range_GC_CurrentNode)
            end if
            if (ModMath_IfBlocksInterSect(Block1%xijk_range,rtp_range_GC_CurrentNode)) then
                iNeighbour=iNeighbour+1
                
                iBlocks_Neighbours(iNeighbour)=CurrentNode%iLeafNode
                rtp_ranges_Neighbours(iNeighbour,:,:)=CurrentNode%rtp_range
            end if
        else
            ! if not leaf node then go into children.

            do ichild=1,size(CurrentNode%children)
                call ModGC_FindNeighbours(Tree,CurrentNode%Children(iChild),&
                    Block1,NumNeighbours,iNeighbour,iBlocks_Neighbours,&
                    rtp_ranges_Neighbours)
            end do
        end if
    end subroutine ModGC_FindNeighbours

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
    subroutine ModGC_GetGC_Sources_Single(Tree,Block1,MpiSize)
        implicit none
        type(YYTree),intent(in)         :: Tree
        type(Block),target              :: Block1
        integer,intent(in)              :: MpiSize
        
        integer                         :: NumNeighbours,NumGCNeighbours
        integer                         :: iNeighbour,iGCNeighbour,iGC_target
        type(GC_target),pointer         :: GC_target1,GC_source1
        real,allocatable                :: rtp_ranges_Neighbours(:,:,:)
        integer,allocatable             :: iBlocks_Neighbours(:)
        logical,allocatable             :: if_GC_neighbours(:)
        type(Block),target              :: Block2

        ! First count the number of neighbours. 
        ! Loop Yin and Yang.
        ! Default it to be 0
        NumNeighbours=0
        call ModGC_CountNeighbours(Tree,Tree%Yin,Block1,NumNeighbours)
        call ModGC_CountNeighbours(Tree,Tree%Yang,Block1,NumNeighbours)

        ! allocate:
        ! 1. iBlocks_Neighbours: the iBlocksGlobal of all neighbours
        ! 2. xijk_ranges_neighbours: xijk_ranges of them
        ! 3. if_GC_neighbours: if each of them is GC neighbour
        !
        allocate(iBlocks_Neighbours(NumNeighbours))
        allocate(rtp_ranges_Neighbours(NumNeighbours,3,2))
        allocate(if_GC_neighbours(NumNeighbours)); if_GC_neighbours=.false.

        ! Then, find the neighbours. Call Yin and Yang.
        ! initialize iNeighbour to be 0
        !
        iNeighbour=0
        call ModGC_FindNeighbours(Tree,Tree%Yin,Block1,NumNeighbours,iNeighbour,&
            iBlocks_Neighbours,rtp_ranges_Neighbours)
        call ModGC_FindNeighbours(Tree,Tree%Yang,Block1,NumNeighbours,iNeighbour,&
            iBlocks_Neighbours,rtp_ranges_Neighbours)

        ! count the number of GC_neighbours
        !
        NumGCNeighbours=0
        do iNeighbour=1,NumNeighbours
            ! initiate each block in turn to get its GC positions
            ! and then locate the target blocks of the GCs
            !
            call ModBlock_InitGrid(Block2,rtp_ranges_Neighbours(iNeighbour,:,:))
            call ModGC_FindGC_Single(Tree,Block2)
            call ModGC_GetGC_targets_single(Tree,Block2,MpiSize)

            ! loop by all the unique GC target blocks to see whether
            ! at least one GC falls in block1.
            !
            do iGC_target=1,Block2%nGC_targets
                if (Block2%GC_targets(iGC_target)%iBlock==Block1%iBlock) then
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
            call ModBlock_InitGrid(Block2,rtp_ranges_Neighbours(iNeighbour,:,:))
            call ModGC_FindGC_Single(Tree,Block2)
            call ModGC_GetGC_targets_single(Tree,Block2,MpiSize)

            ! loop the GC targets to find those fall into block1
            ! and copy them to GC_sources of block1
            !
            do iGC_target=1,Block2%nGC_targets
                if (Block2%GC_targets(iGC_target)%iBlock==Block1%iBlock) then
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
                    allocate(GC_source1%primitive_list(GC_source1%nGC,nvar))
                    GC_source1%ijk_list=GC_target1%ijk_list

                    if (Block2%if_yin .eqv. Block1%if_yin) then
                        GC_source1%xijk_list=GC_target1%xijk_list
                    else
                        GC_source1%xijk_list=ModYinYang_CoordConv_1D(GC_target1%xijk_list,GC_target1%nGC)
                    end if
                    exit
                end if
            end do
            call ModBlock_deallocate(Block2)
        end do

        ! at last, deallocate the two lists/
        !
        deallocate(iBlocks_Neighbours)
        deallocate(rtp_ranges_Neighbours)
        deallocate(if_GC_neighbours)
    end subroutine ModGC_GetGC_Sources_Single

    ! This subroutine is called after ModCommunication_GetnGC_sources
    ! completes getting the n of GC_sources for each Block.
    ! This subroutine aims to allocate the GC_sources.
    subroutine ModGC_SetGC_Sources_Single(Tree,Block1,MpiSize)
        implicit none
        type(YYTree),intent(in)         ::  Tree
        type(Block),target,intent(in)   ::  Block1
        integer,intent(in)              ::  MpiSize

        integer                         ::  iGC_source,iGC_target
        type(GC_target),pointer         ::  GC_target1,GC_source1
        type(Block),target              ::  Block2

        ! Loop all the GC_sources
        do iGC_source=1,Block1%nGC_sources
            GC_source1=>Block1%GC_sources(iGC_source)
            !print *,GC_source1%iBlock
            ! First find the rtp_range of the target Block
            ! based on the iBlock of the GC_source.
            ! Then initiate the grid. Find GC positions and targets

            call ModBlock_Init(Block2,GC_source1%iBlock,&
                YinYangTree_Get_rtp_range(Tree,GC_source1%iBlock),&
                if_yin=GC_source1%if_yin,if_SSM=.false.,if_use_actual_nvar=.false.)
            call ModGC_FindGC_Single(Tree,Block2)
            call ModGC_GetGC_targets_single(Tree,Block2,MpiSize)

            ! loop the GC targets to find those fall into block1
            ! and copy them to GC_sources of block1
            !
            do iGC_target=1,Block2%nGC_targets

                if (Block2%GC_targets(iGC_target)%iBlock==Block1%iBlock) then
                    GC_target1=>Block2%GC_targets(iGC_target)
                    GC_source1%nGC=GC_target1%nGC
                    allocate(GC_source1%xijk_list(GC_source1%nGC,3))
                    allocate(GC_source1%primitive_list(GC_source1%nGC,nvar))

                    if (Block2%if_yin .eqv. Block1%if_yin) then
                        GC_source1%xijk_list=GC_target1%xijk_list
                    else
                        GC_source1%xijk_list=ModYinYang_CoordConv_1D(GC_target1%xijk_list,GC_target1%nGC)
                    end if
                    exit
                end if
            end do
            call ModBlock_deallocate(Block2)
        end do

    end subroutine ModGC_SetGC_Sources_Single
    
end module ModGC