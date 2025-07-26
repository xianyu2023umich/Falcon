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
    use ModBlock,           only:   BlockType,Multigrid_level,GC_target
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
        real                            ::  rtp(3)                  ! rtp of one point

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

    subroutine ModGC_FindHC_Single(Tree,MGL1)
        implicit none
        type(YYTree),intent(in)         ::  Tree
        type(Multigrid_level)           ::  MGL1
        type(TreeNode),pointer          ::  Node1                   ! one Node
        integer                         ::  i,j,k                   ! grid
        real                            ::  rtp(3)                  ! rtp of one point

        ! Find the node of each GC
        do i=0,MGL1%mi+1
            do j=0,MGL1%mj+1
                do k=0,MGL1%mk+1
                    rtp=[MGL1%xi_I(i),MGL1%xj_I(j),MGL1%xk_I(k)]

                    if (i<1 .or. j<1 .or. k<1 .or. i>MGL1%mi .or. j>MGL1%mj .or. k>MGL1%mk) then
                        Node1=>YinYangTree_FindNode(Tree,rtp,MGL1%if_yin)
                        if (associated(Node1))then
                            MGL1%HC_iBlock_III(i,j,k)=Node1%iLeafNode
                        end if
                    end if
                end do
            end do
        end do
    end subroutine ModGC_FindHC_Single

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

    subroutine ModGC_GetHC_Targets_single(Tree,MGL1)
        implicit none
        type(YYTree)                    ::  Tree                                 ! Tree
        type(Multigrid_level),target    ::  MGL1                                 ! The multigrid level to set

        integer                         ::  i,j,k,l
        integer                         ::  iUniqueHCiBlock,iHC_current

        type(GC_target),pointer         ::  HC_target1
        integer,allocatable             ::  UniqueHCiBlocksGlobal(:)
        integer,allocatable             ::  nHC_list(:)
        integer                         ::  nUniqueHCiBlocksGlobal
        integer                         ::  HC_tmp1((MGL1%mi+2)*(MGL1%mj+2)*(MGL1%mk+2))
        integer                         ::  HC_tmp2((MGL1%mi+2)*(MGL1%mj+2)*(MGL1%mk+2))
        integer                         ::  mi,mj,mk
        integer                         ::  ijk_if_HC(1:3)
        logical                         ::  if_new_HCiBlock

        ! Initialize
        nUniqueHCiBlocksGlobal=0
        HC_tmp2=0

        ! Count how many unique elements are there
        ! before allocating the arrays
        ! At the same time count how many GCs for 
        ! each unique block.
        ! But -- only for those ghost cells involved
        ! in the heat conduction.

        mi=MGL1%mi
        mj=MGL1%mj
        mk=MGL1%mk

        do i=0,mi+1
            do j=0,mj+1
                do k=0,mk+1
                    ijk_if_HC=0
                    if (i==0 .or. i==mi+1) ijk_if_HC(1)=1
                    if (j==0 .or. j==mj+1) ijk_if_HC(2)=1
                    if (k==0 .or. k==mk+1) ijk_if_HC(3)=1
                    if (sum(ijk_if_HC)==1 .and. MGL1%HC_iBlock_III(i,j,k) /= -777) then
                        if_new_HCiBlock=.True.
                        if (nUniqueHCiBlocksGlobal .gt. 0) then
                            do l=1,nUniqueHCiBlocksGlobal
                                if (HC_tmp1(l) == MGL1%HC_iBlock_III(i,j,k)) then 
                                    if_new_HCiBlock=.false.
                                    HC_tmp2(l)=HC_tmp2(l)+1
                                end if
                            end do
                        end if
                        if (if_new_HCiBlock) then 
                            nUniqueHCiBlocksGlobal=nUniqueHCiBlocksGlobal+1
                            HC_tmp1(nUniqueHCiBlocksGlobal)=MGL1%HC_iBlock_III(i,j,k) 
                            HC_tmp2(nUniqueHCiBlocksGlobal)=1
                        end if
                    end if
                end do
            end do
        end do

        ! Allocate based on what we got above
        allocate(UniqueHCiBlocksGlobal  (   nUniqueHCiBlocksGlobal))
        allocate(nHC_list               (   nUniqueHCiBlocksGlobal))
        UniqueHCiBlocksGlobal   =   HC_tmp1(1: nUniqueHCiBlocksGlobal)
        nHC_list                =   HC_tmp2(1: nUniqueHCiBlocksGlobal) 
        MGL1%nHC_targets      =           nUniqueHCiBlocksGlobal
        allocate(MGL1%HC_targets(         nUniqueHCiBlocksGlobal))

        ! loop all the GC_targets

        do iUniqueHCiBlock=1,nUniqueHCiBlocksGlobal

            HC_target1=>MGL1%HC_targets(iUniqueHCiBlock)
            HC_target1%iBlock=UniqueHCiBlocksGlobal(iUniqueHCiBlock)
            HC_target1%if_yin=HC_target1%iBlock.le.Tree%NumLeafNodes_YinYang(1)

            ! see which rank the block of this GC_target belongs to
            HC_target1%iRank=ModAllocation_GetRank(HC_target1%iBlock)
            
            ! allocate
            HC_target1%nGC=nHC_list(iUniqueHCiBlock)
            allocate(HC_target1%xijk_list(nHC_list(iUniqueHCiBlock),3))
            allocate(HC_target1%ijk_list(nHC_list(iUniqueHCiBlock),3))
            allocate(HC_target1%primitive_list(nHC_list(iUniqueHCiBlock),1))
            HC_target1%primitive_list(:,:)=0.0

            ! loop all the points to fill in the GC_target1
            iHC_current=0
            do i=0,mi+1
                do j=0,mj+1
                    do k=0,mk+1
                        ijk_if_HC=0
                        if (i==0 .or. i==mi+1) ijk_if_HC(1)=1
                        if (j==0 .or. j==mj+1) ijk_if_HC(2)=1
                        if (k==0 .or. k==mk+1) ijk_if_HC(3)=1
                        if (sum(ijk_if_HC)==1) then
                            if (MGL1%HC_iBlock_III(i,j,k) .eq. UniqueHCiBlocksGlobal(iUniqueHCiBlock)) then

                                iHC_current=iHC_current+1
                                HC_target1%ijk_list (iHC_current,:)=[i,j,k]
                                HC_target1%xijk_list(iHC_current,:)=[MGL1%xi_I(i),MGL1%xj_I(j),MGL1%xk_I(k)]
                            end if
                        end if
                    end do
                end do
            end do
        end do
        
        deallocate(UniqueHCiBlocksGlobal)
        deallocate(nHC_list)
    end subroutine ModGC_GetHC_Targets_single

    ! Ok that's it temporarily no more.
    ! GC(HC) sources are set in ModCommunication.f90
    ! Previous versions reserve some very old methods, which worked very slowly.
end module ModGC