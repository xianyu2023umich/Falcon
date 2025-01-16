Module ModAllocation

    contains

    ! Get the table including the start and end node for each rank

    subroutine ModAllocation_GetTable(NumLeafNodes,iLeafNode_ranges,MpiSize)
        implicit none
        integer,intent(in)          ::  NumLeafNodes            ! Num of leaf nodes   
        integer,intent(in)          ::  MpiSize                 ! Mpi Size
        integer,allocatable         ::  iLeafNode_ranges(:,:)   ! output

        integer                     ::  nNodesPerRank           ! num Nodes per rank
        integer                     ::  iLeafNodeStart,&        ! start i node
                                        iLeafNodeEnd            ! end i node
        integer                     ::  MpiRank                 ! MpiRank
        
        ! First allocate the table and initialize it to from -1 to -2
        allocate(iLeafNode_ranges(0:MpiSize-1,2))
        iLeafNode_ranges(0:MpiSize-1,1)=-1
        iLeafNode_ranges(0:MpiSize-1,2)=-2

        ! Num of nodes per rank
        nNodesPerRank = NumLeafNodes / MpiSize
        if (mod(NumLeafNodes,MpiSize).ne.0) nNodesPerRank = nNodesPerRank + 1

        do MpiRank=0,MpiSize-1
            iLeafNodeStart = 1 + MpiRank * nNodesPerRank
            iLeafNodeEnd = min(nNodesPerRank + MpiRank * nNodesPerRank,NumLeafNodes)
            iLeafNode_ranges(MpiRank,:)=[iLeafNodeStart,iLeafNodeEnd]
        end do

    end subroutine ModAllocation_GetTable

    ! Get the nodes for this rank

    subroutine ModAllocation_GetNodes(NumLeafNodes,iLeafNodes,MpiSize,MpiRank)
        implicit none
        integer,intent(in)          ::  NumLeafNodes            
        integer,intent(in)          ::  MpiSize,MpiRank
        integer,allocatable         ::  iLeafNodes(:)
        integer                     ::  nNodesPerRank
        integer                     ::  iLeafNodeStart,&
                                        iLeafNodeEnd
        integer                     ::  iLeafNodeLocal

        ! determine how many nodes are there in one rank at most
        nNodesPerRank = NumLeafNodes / MpiSize
        if (mod(NumLeafNodes,MpiSize).ne.0) nNodesPerRank = nNodesPerRank + 1

        ! get the start and end iNodes of this rank
        iLeafNodeStart = 1 + MpiRank * nNodesPerRank
        iLeafNodeEnd = min(nNodesPerRank + MpiRank * nNodesPerRank,NumLeafNodes)

        ! then the total length is iNodeEnd - iNodeStart + 1
        ! so we can allocate iNodes
        allocate(iLeafNodes(iLeafNodeEnd - iLeafNodeStart + 1))

        ! fulfill iNodes
        do iLeafNodeLocal = 1,iLeafNodeEnd - iLeafNodeStart + 1
            iLeafNodes(iLeafNodeLocal) = iLeafNodeStart + (iLeafNodeLocal - 1)
        end do
    end subroutine ModAllocation_GetNodes

    subroutine ModAllocation_GetRank(NumLeafNodes,iLeafNode,MpiSize,MpiRank)

        implicit none

        integer,intent(in)          :: NumLeafNodes
        integer,intent(in)          :: iLeafNode,MpiSize
        integer,intent(out)         :: MpiRank
        integer                     :: nNodesPerRank

        ! determine how many nodes are there in one rank at most
        nNodesPerRank = NumLeafNodes / MpiSize
        if (mod(NumLeafNodes,MpiSize).ne.0) nNodesPerRank = nNodesPerRank + 1
        MpiRank = (iLeafNode - 1) / nNodesPerRank
    end subroutine ModAllocation_GetRank

end module ModAllocation