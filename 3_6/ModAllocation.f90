Module ModAllocation

    use ModParameters,only  :   MpiSize,MpiRank

    integer,allocatable     ::  ranges_of_ranks(:,:)
    integer,allocatable     ::  ranks_of_iNodes(:)
    integer,allocatable     ::  iLocalLeafNodes(:)
    integer                 ::  nLocalLeafNodes
    integer                 ::  ModAllocation_status=0

    contains

    subroutine ModAllocation_Init(NumLeafNodes)
        implicit none
        integer,intent(in)  ::  NumLeafNodes
        integer             ::  iLeafNode,iRank

        allocate(ranges_of_ranks(0:MpiSize-1,2))
        allocate(ranks_of_iNodes(NumLeafNodes))

        ! Get the ranges as function of ranks
        ! First initialize ileafnode to 1 waiting 
        ! to be added up in the loop.

        iLeafNode=1
        do iRank=0,MpiSize-1
            ranges_of_ranks(iRank,1)=iLeafNode

            ! If iRank<mod(NumLeafNodes,MpiSize),
            ! there should be NumLeafNodes/MpiSize+1
            ! nodes. Otherwise just NumLeafNodes/MpiSize.

            iLeafNode=iLeafNode+merge(&
                NumLeafNodes/MpiSize+1,&
                NumLeafNodes/MpiSize,&
                iRank<mod(NumLeafNodes,MpiSize))

            ranges_of_ranks(iRank,2)=iLeafNode-1
        end do

        ! THen get the corresponding rank of each Node.
        ! Initialize Rank=0 then loop.

        iRank=0
        do iLeafNode=1,NumLeafNodes

            ! If current iLeafNode exceed top of
            ! current iRank's range then move
            ! to the next one.

            if (iLeafNode>ranges_of_ranks(iRank,2)) iRank=iRank+1
            ranks_of_iNodes(iLeafNode)=iRank
        end do

        ! Then get the local iNode ranges
        nLocalLeafNodes=ranges_of_ranks(MpiRank,2)-ranges_of_ranks(MpiRank,1)+1
        allocate(iLocalLeafNodes(nLocalLeafNodes))
        
        do iLeafNode=1,nLocalLeafNodes
            iLocalLeafNodes(iLeafNode)=ranges_of_ranks(MpiRank,1)+iLeafNode-1
        end do

        ! The module has been initialized.
        ModAllocation_status=1
    end subroutine ModAllocation_Init

    function ModAllocation_GetRank(iLeafNode) result(iRank)
        implicit none
        integer,intent(in)  ::  iLeafNode
        integer             ::  iRank

        iRank=ranks_of_iNodes(iLeafNode)
    end function ModAllocation_GetRank

end module ModAllocation