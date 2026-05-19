module ModGC
    ! Couple of things to do this time:
    ! 1. Now we have multigrid. Each level should be done separately.
    ! 2. Previously, there's a sendrecv between each two neighbouring blocks.
    !   There might be several sends between two mpi ranks -- which is not efficient.
    !   We need a GC_Send and GC_Recv owned by each rank to gather all the message
    !   of the local blocks.
    ! 3. Delete several previous unused subroutines. e.g. ModGC_CommunicateGCLocal --
    !   I've moved it to ModCommunication.f90.

    use ModYinYangTree,     only:   YYTree,TreeNode,&
                                    YinYangTree_FindNode
    use ModBlock,           only:   BlockType,GC_target
    use ModParameters,      only:   MpiRank,MpiSize,&
                                    ni,nj,nk,ng,nvar,&
                                    Multigrid_nLevels
    use ModAllocation,      only:   ModAllocation_GetRank,&
                                    ranges_of_ranks

    contains

    ! There is no need to write set all.
    ! That should be done in mod communication.
    ! For one block, loop each cell to see where it belongs to.
    subroutine ModGC_FindGC_Single(Tree,Block1)
        implicit none
        type(YYTree),intent(in)         ::  Tree
        type(BlockType),intent(inout)   ::  Block1                  ! one Block
        type(TreeNode),pointer          ::  Node1                   ! one Node
        integer                         ::  i,j,k                   ! grid
        real(8)                         ::  rtp(3)                  ! rtp of one point

        ! Initialize to be -777
        Block1%GC_iBlocks_III=-777                     
        
        ! Find the node of each GC
        do i=-ng+1,ni+ng
            do j=-ng+1,nj+ng
                do k=-ng+1,nk+ng
                    rtp=[Block1%xi_I(i),Block1%xj_I(j),Block1%xk_I(k)]

                    if (i<1 .or. j<1 .or. k<1 .or. i>ni .or. j>nj .or. k>nk) then
                        Node1=>YinYangTree_FindNode(Tree,rtp,Block1%if_yin)
                        if (associated(Node1))then
                            Block1%GC_iBlocks_III(i,j,k)=Node1%iLeafNode
                        end if
                    end if
                end do
            end do
        end do
    end subroutine ModGC_FindGC_Single

    subroutine ModGC_GetGC_Targets_single(Tree,Block1)
        implicit none
        type(YYTree)                    ::  Tree                                 ! Tree
        type(BlockType),target          ::  Block1                               ! The block to set

        integer                         ::  i,j,k,l
        integer                         ::  iUniqueGCiBlock,iGC_current

        type(GC_target),pointer         ::  GC_target1
        integer,allocatable             ::  UniqueGCiBlocksGlobal(:)
        integer,allocatable             ::  nGC_list(:)
        integer                         ::  nUniqueGCiBlocksGlobal
        integer                         ::  GC_tmp1((ni+2*ng)*(nj+2*ng)*(nk+2*ng))
        integer                         ::  GC_tmp2((ni+2*ng)*(nj+2*ng)*(nk+2*ng))
        logical                         ::  if_new_GCiBlock

        ! Initialize
        nUniqueGCiBlocksGlobal=0
        GC_tmp2=0

        ! Count how many unique elements are there
        ! before allocating the arrays
        ! At the same time count how many GCs for 
        ! each unique block.

        do i=-ng+1,ni+ng
            do j=-ng+1,nj+ng
                do k=-ng+1,nk+ng
                    if ((i<1 .or. i>ni .or. j<1 .or. j>nj .or. k<1 .or. k>nk)  .and. Block1%GC_iBlocks_III(i,j,k) /= -777) then
                        if_new_GCiBlock=.True.
                        if (nUniqueGCiBlocksGlobal .gt. 0) then
                            do l=1,nUniqueGCiBlocksGlobal
                                if (GC_tmp1(l) == Block1%GC_iBlocks_III(i,j,k)) then 
                                    if_new_GCiBlock=.false.
                                    GC_tmp2(l)=GC_tmp2(l)+1
                                end if
                            end do
                        end if ! if (nUniqueGCiBlocksGlobal .gt. 0)
                        if (if_new_GCiBlock) then 
                            nUniqueGCiBlocksGlobal=nUniqueGCiBlocksGlobal+1
                            GC_tmp1(nUniqueGCiBlocksGlobal)=Block1%GC_iBlocks_III(i,j,k)
                            GC_tmp2(nUniqueGCiBlocksGlobal)=1
                        end if ! if (if_new_GCiBlock)
                    end if
                end do
            end do
        end do ! do i=-ng+1,ni+ng; do j=-ng+1,nj+ng; do k=-ng+1,nk+ng

        ! Allocate based on what we got above
        allocate(UniqueGCiBlocksGlobal  (   nUniqueGCiBlocksGlobal))
        allocate(nGC_list               (   nUniqueGCiBlocksGlobal))
        UniqueGCiBlocksGlobal   =   GC_tmp1(1: nUniqueGCiBlocksGlobal)
        nGC_list                =   GC_tmp2(1: nUniqueGCiBlocksGlobal) 
        Block1%nGC_targets      =           nUniqueGCiBlocksGlobal
        allocate(Block1%GC_targets(         nUniqueGCiBlocksGlobal))

        ! loop all the GC_targets

        do iUniqueGCiBlock=1,nUniqueGCiBlocksGlobal

            GC_target1=>Block1%GC_targets(iUniqueGCiBlock)
            GC_target1%iBlock=UniqueGCiBlocksGlobal(iUniqueGCiBlock)
            GC_target1%if_yin=GC_target1%iBlock.le.Tree%NumLeafNodes_YinYang(1)

            ! see which rank the block of this GC_target belongs to
            GC_target1%iRank=ModAllocation_GetRank(GC_target1%iBlock)
            
            ! allocate
            GC_target1%nGC=nGC_list(iUniqueGCiBlock)
            allocate(GC_target1%xijk_list(nGC_list(iUniqueGCiBlock),3))
            allocate(GC_target1%ijk_list(nGC_list(iUniqueGCiBlock),3))
            allocate(GC_target1%primitive_list(nGC_list(iUniqueGCiBlock),nvar))
            GC_target1%primitive_list(:,:)=0.0

            ! loop all the points to fill in the GC_target1
            iGC_current=0
            do i=-ng+1,ni+ng
                do j=-ng+1,nj+ng
                    do k=-ng+1,nk+ng
                        if (i<1 .or. i>ni .or. j<1 .or. j>nj .or. k<1 .or. k>nk) then
                            if (Block1%GC_iBlocks_III(i,j,k) .eq. UniqueGCiBlocksGlobal(iUniqueGCiBlock)) then

                                iGC_current=iGC_current+1
                                GC_target1%ijk_list (iGC_current,:)=[i,j,k]
                                GC_target1%xijk_list(iGC_current,:)=[Block1%xi_I(i),Block1%xj_I(j),Block1%xk_I(k)]
                            end if
                        end if
                    end do
                end do
            end do
        end do
        
        deallocate(UniqueGCiBlocksGlobal)
        deallocate(nGC_list)
    end subroutine ModGC_GetGC_Targets_single
end module ModGC