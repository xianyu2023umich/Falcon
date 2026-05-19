module ModCommunication

    use ModMath,        only:   ModMath_1D3D_interpolate_1D1D,&
                                ModMath_3D_interpolate_1D
    use ModBlock,       only:   BlockType,GC_target
    use ModYinYang,     only:   ModYinYang_CoordConv_1D,&
                                ModYinYang_CoordConv_0D,&
                                ModYinYang_VecConv_1D
    use ModYinYangTree, only:   YYTree,GC_interface
    use ModGC,          only:   ModGC_FindGC_Single,&
                                ModGC_GetGC_Targets_single
    use ModParameters,  only:   MpiSize,MpiRank,&
                                ni,nj,nk,ng,nvar
    use ModControl,     only:   nMaxBlocksPerRank
    use ModAllocation,  only:   ranges_of_ranks
    use MPI

    implicit none

    contains

    subroutine ModCommunication_SetGCAll(Tree)
        implicit none
        type(YYTree),target             ::  Tree                    ! the Tree
        type(BlockType),pointer         ::  Block1                  ! one Block
        integer                         ::  iBlock
        integer                         ::  ierr

        ! First determine where each ghost cell belongs to
        ! by using the ModGC_FindGC_Single subroutine.
        do iBlock=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iBlock)
            call ModGC_FindGC_Single(Tree,Block1)
        end do
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (MpiRank==0) write(*,*)'Completed finding GC...'

        ! Then get the GC_targets for each block
        do iBlock=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iBlock)
            call ModGC_GetGC_Targets_single(Tree,Block1)
        end do
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (MpiRank==0) write(*,*)'Completed setting GC_targets...'

        ! Then set the GC_sources for each block
        call ModCommunication_GetnGC_sources(Tree)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (MpiRank==0) write(*,*)'Completed setting nGC_sources...'

        ! Then set the GC_sources for each block
        call ModCommunication_SetGC_sources(Tree)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (MpiRank==0) write(*,*)'Completed setting GC_sources...'

        ! Setting up GC senders and receivers
        call ModCommunication_SetGC_Sender(Tree)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (MpiRank==0) write(*,*)'Completed setting GC_senders...'

        call ModCommunication_SetGC_Receiver(Tree)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if (MpiRank==0) write(*,*)'Completed setting GC_receivers...'
    end subroutine ModCommunication_SetGCAll

    subroutine ModCommunication_GetnGC_sources(Tree)
        implicit none
        type(YYTree),target             ::  Tree
        type(BlockType),pointer         ::  Block_source,Block_target
        type(GC_target),pointer         ::  GC_target1,GC_source1
        integer                         ::  nGC_sources_table_local(Tree%NumLeafNodes)
        integer                         ::  nGC_sources_table_global(Tree%NumLeafNodes)
        integer                         ::  nGC_international_sources_table_local(Tree%NumLeafNodes)
        integer                         ::  nGC_international_sources_table_global(Tree%NumLeafNodes)
        integer                         ::  iLocalBlock,iGC_target,iGC_source
        integer                         ::  iRecv,nRecv
        integer,allocatable             ::  send_message_all(:),recv_message_all(:)
        integer                         ::  ierr
        integer                         ::  nSend,iSend
        integer,allocatable             ::  requests_send(:),requests_recv(:)
        
        ! Get the local table for GC_sources of each block
        nGC_sources_table_local=0
        nGC_international_sources_table_local=0
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)

            do iGC_target=1,Block_source%nGC_targets
                GC_target1=>Block_source%GC_targets(iGC_target)

                ! This is all the GC_sources.
                nGC_sources_table_local(GC_target1%iBlock)=&
                        nGC_sources_table_local(GC_target1%iBlock)+1

                ! Only the international GC_sources.
                if (MpiRank/=GC_target1%iRank) then
                    nGC_international_sources_table_local(GC_target1%iBlock)=&
                        nGC_international_sources_table_local(GC_target1%iBlock)+1
                end if
            end do
        end do

        ! Then reduce it to the global table
        call MPI_ALLReduce(nGC_sources_table_local,nGC_sources_table_global,&
            Tree%NumLeafNodes,MPI_integer,MPI_SUM,MPI_COMM_WORLD,ierr)

        call MPI_ALLReduce(nGC_international_sources_table_local,nGC_international_sources_table_global,&
            Tree%NumLeafNodes,MPI_integer,MPI_SUM,MPI_COMM_WORLD,ierr)

        ! Loop the local blocks and allocate GC_sources for them
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)
            allocate(Block_source%GC_sources(nGC_sources_table_global(Block_source%iBlock)))
            Block_source%nGC_sources=0
        end do

        ! Before the loop, we need to get the n of messages to send
        Tree%LocalBlocks(:)%nSend=0
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)
            do iGC_target=1,Block_source%nGC_targets
                GC_target1=>Block_source%GC_targets(iGC_target)

                if (GC_target1%iRank/=MpiRank) then
                    Block_source%nSend=Block_source%nSend+1
                end if
            end do
        end do
        nSend=sum(Tree%LocalBlocks(:)%nSend)

        ! Get the n of messages to receive
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)
            Block_source%nRecv=nGC_international_sources_table_global(Block_source%iBlock)
        end do
        nRecv=sum(Tree%LocalBlocks(:)%nRecv)

        allocate(send_message_all(nSend*4))
        allocate(recv_message_all(nRecv*4))
        allocate(requests_send(nSend))
        allocate(requests_recv(nRecv))

        ! Recv the iBlock messages for GC_sources
        do iRecv=1,nRecv
            call MPI_IRECV(recv_message_all(iRecv*4-3:iRecv*4),4,mpi_integer,&
                MPI_ANY_SOURCE,11,MPI_COMM_WORLD,requests_recv(iRecv),ierr)
        end do

        ! Then loop the local blocks and send the iBlock of GC_targets
        ! to the target ranks.

        iSend=0
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)
            do iGC_target=1,Block_source%nGC_targets
                GC_target1=>Block_source%GC_targets(iGC_target)

                if (GC_target1%iRank/=MpiRank) then
                    iSend=iSend+1
                    send_message_all(iSend*4-3:iSend*4)=[MpiRank,Block_source%iBlock,GC_target1%iBlock,GC_target1%nGC]

                    call MPI_ISEND(send_message_all(iSend*4-3:iSend*4),4,mpi_integer,&
                        GC_target1%iRank,11,MPI_COMM_WORLD,requests_send(iSend),ierr)
                else
                    Block_target=>Tree%LocalBlocks(GC_target1%iBlock-ranges_of_ranks(MpiRank,1)+1)
                    Block_target%nGC_sources=Block_target%nGC_sources+1
                    GC_source1=>Block_target%GC_sources(Block_target%nGC_sources)

                    GC_source1%iRank=MpiRank
                    GC_source1%iBlock=Block_source%iBlock
                    GC_source1%if_yin=Block_source%iBlock.le.Tree%NumLeafNodes_YinYang(1)
                    GC_source1%nGC=GC_target1%nGC
                    allocate(GC_source1%xijk_list(GC_source1%nGC,3))
                    allocate(GC_source1%ijk_list(GC_source1%nGC,3))
                    allocate(GC_source1%primitive_list(GC_source1%nGC,nvar))
                end if
            end do
        end do
        ! Wait for all the sends to be done
        call MPI_WAITALL(nSend,requests_send,MPI_STATUSES_IGNORE,ierr)
        call MPI_WAITALL(nRecv,requests_recv,MPI_STATUSES_IGNORE,ierr)
        
        ! Assign and allocate GC_sources
        do iRecv=1,nRecv
            Block_target=>Tree%LocalBlocks(recv_message_all(iRecv*4-1)-ranges_of_ranks(MpiRank,1)+1)
            Block_target%nGC_sources=Block_target%nGC_sources+1

            GC_source1=>Block_target%GC_sources(Block_target%nGC_sources)

            GC_source1%iRank=recv_message_all(iRecv*4-3)
            GC_source1%iBlock=recv_message_all(iRecv*4-2)
            GC_source1%if_yin=recv_message_all(iRecv*4-2).le.Tree%NumLeafNodes_YinYang(1)
            GC_source1%nGC=recv_message_all(iRecv*4)
            allocate(GC_source1%xijk_list(GC_source1%nGC,3))
            allocate(GC_source1%ijk_list(GC_source1%nGC,3))
            allocate(GC_source1%primitive_list(GC_source1%nGC,nvar))
        end do

        deallocate(requests_recv)
        deallocate(requests_send)
        deallocate(send_message_all)
        deallocate(recv_message_all)

        ! At last we get the total n of requests

        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)
            nSend=0
            if (Block_source%nGC_sources>0) then
                do iGC_source=1,Block_source%nGC_sources
                    if (Block_source%GC_sources(iGC_source)%iRank/=MpiRank) nSend=nSend+1
                end do
                if (nSend>0) allocate(Block_source%GC_requests(nSend))
            endif 
        end do
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    end subroutine ModCommunication_GetnGC_sources

    subroutine ModCommunication_SetGC_sources(Tree)
        implicit none
        type(YYTree),intent(in),target  ::  Tree
        type(BlockType),pointer         ::  Block_source,Block_target

        integer                         ::  iLocalBlock
        integer                         ::  iGC_target,iGC_source
        type(GC_target),pointer         ::  GC_source1,GC_target1
        integer                         ::  iSend,iRecv,nSend,nRecv
        integer,allocatable             ::  requests_send(:),requests_recv(:)

        integer                         ::  tag,ierr

        ! xijk_list

        ! Prepare the requests_recv
        iRecv=0
        nRecv=sum(Tree%LocalBlocks(:)%nRecv)
        allocate(requests_recv(nRecv))

        ! Post the irecvs at first
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source => Tree%LocalBlocks(iLocalBlock)
    
            do iGC_source=1,size(Block_source%GC_sources)
                GC_source1 => Block_source%GC_sources(iGC_source)

                ! Local sends have been done before so only deal with the remote sends.
                
                if (GC_source1%iRank/=MpiRank) then
                    iRecv=iRecv+1
                    ! Allocate the buffer, get the tag and then receive
                    tag=ModCommunication_GetTag(GC_source1%iBlock,Block_source%iBlock,GC_source1%iRank,MpiRank)
                    call MPI_IRECV(GC_source1%xijk_list,3*GC_source1%nGC,mpi_double_precision,GC_source1%iRank,tag,MPI_COMM_WORLD,&
                    requests_recv(iRecv),ierr)
                endif ! GC_source1%iRank/=MpiRank
            end do
        end do

        ! Send the GC_targets%xijk_list
        iSend=0
        nSend=sum(Tree%LocalBlocks(:)%nSend)
        allocate(requests_send(nSend))

        ! Post the isends at first
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_target => Tree%LocalBlocks(iLocalBlock)

            do iGC_target=1,size(Block_target%GC_targets)
                GC_target1 => Block_target%GC_targets(iGC_target)

                if (GC_target1%iRank/=MpiRank) then
                    iSend=iSend+1
                    tag = ModCommunication_GetTag(Block_target%iBlock,GC_target1%iBlock,MpiRank,GC_target1%iRank)
                    call MPI_ISEND(GC_target1%xijk_list,&
                        3*GC_target1%nGC,mpi_double_precision,GC_target1%iRank,tag,MPI_COMM_WORLD,requests_send(iSend),ierr)
                else
                    Block_source=>Tree%LocalBlocks(GC_target1%iBlock-ranges_of_ranks(MpiRank,1)+1)

                    ! Go over the GC_sources to see which one corresponds to the GC_target1.
                    do iGC_source=1,size(Block_source%GC_sources)
                        GC_source1=>Block_source%GC_sources(iGC_source)
                        if (GC_source1%iBlock==Block_target%iBlock) then
                            ! First record the ijk_list. It always means the ijk
                            ! at the target side, so it's the same for source.
                            GC_source1%ijk_list=GC_target1%ijk_list

                            ! Then record the xijk_list.
                            ! xijk_list is first calculated at the target side.
                            ! But it's later used at the source side for interpolation.
                            ! So we need to convert it to the source side with if_yin.
                            if (GC_source1%if_yin .eqv. Block_source%if_yin) then
                                GC_source1%xijk_list=GC_target1%xijk_list
                            else
                                GC_source1%xijk_list=ModYinYang_CoordConv_1D(GC_target1%xijk_list,GC_target1%nGC)
                            end if
                            exit
                        end if
                    end do
                end if ! GC_target1%iRank/=MpiRank
            end do ! iGC_target
        end do ! iLocalBlock

        ! Wait for all the recvs to be done
        call MPI_WAITALL(nRecv,requests_recv,MPI_STATUSES_IGNORE,ierr)
        call MPI_WAITALL(nSend,requests_send,MPI_STATUSES_IGNORE,ierr)

        ! YinYang conversion
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source => Tree%LocalBlocks(iLocalBlock)
    
            do iGC_source=1,size(Block_source%GC_sources)
                GC_source1 => Block_source%GC_sources(iGC_source)
                
                if (GC_source1%iRank/=MpiRank) then
                    ! See whether same if_yin or not
                    if (GC_source1%if_yin .neqv. Block_source%if_yin) then
                        GC_source1%xijk_list=ModYinYang_CoordConv_1D(GC_source1%xijk_list,GC_source1%nGC)
                    endif
                endif ! GC_source1%iRank/=MpiRank
            end do
        end do

        ! ijk_list
        iRecv=0
        requests_recv=MPI_REQUEST_NULL

        ! Post the irecvs
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source => Tree%LocalBlocks(iLocalBlock)
            do iGC_source=1,size(Block_source%GC_sources)
                GC_source1 => Block_source%GC_sources(iGC_source)

                ! Local sends have been done before so only deal with the remote sends.
                
                if (GC_source1%iRank/=MpiRank) then
                    iRecv=iRecv+1
                    tag=ModCommunication_GetTag(GC_source1%iBlock,Block_source%iBlock,GC_source1%iRank,MpiRank)
                    call MPI_IRECV(GC_source1%ijk_list,3*GC_source1%nGC,mpi_integer,GC_source1%iRank,tag,MPI_COMM_WORLD,&
                    requests_recv(iRecv),ierr)
                endif ! GC_source1%iRank/=MpiRank
            end do
        end do

        ! Send the GC_targets%ijk_list
        iSend=0
        requests_send=MPI_REQUEST_NULL

        ! Post the isends then

        do iLocalBlock=1,Tree%nLocalBlocks
            Block_target => Tree%LocalBlocks(iLocalBlock)
            do iGC_target=1,size(Block_target%GC_targets)
                GC_target1 => Block_target%GC_targets(iGC_target)
                if (GC_target1%iRank/=MpiRank) then
                    iSend=iSend+1
                    tag = ModCommunication_GetTag(Block_target%iBlock,GC_target1%iBlock,MpiRank,GC_target1%iRank)
                    call MPI_ISEND(GC_target1%ijk_list,&
                        3*GC_target1%nGC,mpi_integer,GC_target1%iRank,tag,MPI_COMM_WORLD,requests_send(iSend),ierr)
                end if ! GC_target1%iRank/=MpiRank
            end do ! iGC_target
        end do ! iLocalBlock

        ! Wait for all the recvs to be done
        call MPI_WAITALL(nRecv,requests_recv,MPI_STATUSES_IGNORE,ierr)
        call MPI_WAITALL(nSend,requests_send,MPI_STATUSES_IGNORE,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    end subroutine ModCommunication_SetGC_sources

    ! Set the GC sender for the tree. 
    ! First see all the distinct neighbour ranks. The number is the number of 
    ! senders to be allocated.
    ! Then see how many block pairs are there. There's NO NEED to see how many
    ! *distinct pairs* because different GC_sources always correspond to
    ! distinct pairs.

    subroutine ModCommunication_SetGC_Sender(Tree)
        implicit none
        type(YYTree),target             ::  Tree
        integer                         ::  iUniqueMpiRank
        integer                         ::  UniqueMpiRank_ptr
        integer                         ::  UniqueMpiRanks(MpiSize)

        ! iBlock_pair means different so choose another name
        integer                         ::  iofBlock_pair,iGC                   

        integer                         ::  iGC_Sender
        integer                         ::  iLocalBlock,iGC_source
        type(BlockType),pointer         ::  Block_source
        type(GC_target),pointer         ::  GC_source1
        type(GC_interface),pointer      ::  GC_Sender1

        ! Loop all local Blocks%GC_sources to find all neighbour ranks.
        ! First set UniqueMpiRanks to be -1. Loop, if find any new neighbour,
        ! record it there.

        UniqueMpiRank_ptr=1
        UniqueMpiRanks=-1

        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)
            do iGC_source=1,Block_source%nGC_sources
                GC_source1=>Block_source%GC_sources(iGC_source)

                if (GC_source1%iRank/=MpiRank) then
                    ! Loop all recorded unique ranks to see if it's already in there.
                    ! If not, record it.
                    do iUniqueMpiRank=1,UniqueMpiRank_ptr
                        if (UniqueMpiRanks(iUniqueMpiRank)==GC_source1%iRank) exit
                    end do
                    if (iUniqueMpiRank==UniqueMpiRank_ptr+1) then
                        UniqueMpiRanks(UniqueMpiRank_ptr)=GC_source1%iRank
                        UniqueMpiRank_ptr=UniqueMpiRank_ptr+1
                    end if
                end if
            end do
        end do

        ! Now the nPairs is UniqueMpiRank_ptr-1.
        ! Allocate the GC_Senders.
        Tree%nGC_Senders=UniqueMpiRank_ptr-1
        allocate(Tree%GC_Senders(Tree%nGC_Senders))
        Tree%GC_Senders(:)%iRank=UniqueMpiRanks(1:Tree%nGC_Senders)

        ! Get the iBlock_pairs for each GC_sender.
        ! Loop all the local blocks and GC_sources to get the unique pairs.
        ! The process is pretty similar to above, but we need a 2D array 
        ! for recording since it's a pair of blocks.
        do iGC_Sender=1,Tree%nGC_Senders

            GC_Sender1=>Tree%GC_Senders(iGC_Sender)
            GC_Sender1%nBlock_Pairs=0
            GC_Sender1%nGC=0

            do iLocalBlock=1,Tree%nLocalBlocks
                Block_source=>Tree%LocalBlocks(iLocalBlock)
                do iGC_source=1,Block_source%nGC_sources
                    GC_source1=>Block_source%GC_sources(iGC_source)

                    ! No need to record. Different GC_sources always correspond to
                    ! different block pairs.

                    if (GC_source1%iRank==GC_Sender1%iRank) then
                        GC_Sender1%nBlock_Pairs=GC_Sender1%nBlock_Pairs+1
                        GC_Sender1%nGC=GC_Sender1%nGC+GC_source1%nGC
                    end if

                end do
            end do

            allocate(GC_Sender1%message(GC_Sender1%nGC,nvar))
            allocate(GC_Sender1%iBlock_pairs(GC_Sender1%nBlock_Pairs,2))
            allocate(GC_Sender1%Block_Pairs_Ptrs(GC_Sender1%nBlock_Pairs,2))
            allocate(GC_Sender1%iLocalBlock_and_iGCtarget_list(GC_Sender1%nBlock_Pairs,2))

            !Loop again to find the iBlock_pairs and the Block_Pairs_Ptrs.
            !We need iGC and iBlock_pair to record the current position.

            iofBlock_pair=1
            iGC=1
            do iLocalBlock=1,Tree%nLocalBlocks
                Block_source=>Tree%LocalBlocks(iLocalBlock)
                do iGC_source=1,Block_source%nGC_sources
                    GC_source1=>Block_source%GC_sources(iGC_source)

                    if (GC_source1%iRank==GC_Sender1%iRank) then
                        ! Record the iBlock_pairs.
                        GC_Sender1%iBlock_pairs(iofBlock_pair,1)=Block_source%iBlock
                        GC_Sender1%iBlock_pairs(iofBlock_pair,2)=GC_source1%iBlock

                        ! Ptrs
                        GC_Sender1%Block_Pairs_Ptrs(iofBlock_pair,1)=iGC
                        iGC=iGC+GC_source1%nGC
                        GC_Sender1%Block_Pairs_Ptrs(iofBlock_pair,2)=iGC-1

                        ! iLocalBlock and iGCtarget
                        GC_Sender1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,1)=iLocalBlock
                        GC_Sender1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,2)=iGC_source

                        iofBlock_pair=iofBlock_pair+1
                    end if
                end do
            end do
        end do
    end subroutine ModCommunication_SetGC_Sender

    ! Then set the GC_receivers.
    ! First find neighbours like the GC_Senders so we know where to receive.
    ! Then find how many Block pairs in each Receiver is still the same.
    ! But the iBlock_pairs and ptrs are different -- we will get them by
    ! sending them from GC_Senders.

    subroutine ModCommunication_SetGC_Receiver(Tree)
        implicit none
        type(YYTree),target             ::  Tree
        integer                         ::  iUniqueMpiRank
        integer                         ::  UniqueMpiRank_ptr
        integer                         ::  UniqueMpiRanks(MpiSize)

        integer                         ::  iofBlock_pair

        integer                         ::  iLocalBlock,iGC_target,iGC_sender,iGC_receiver
        type(BlockType),pointer         ::  Block_target
        type(GC_target),pointer         ::  GC_target1
        type(GC_interface),pointer      ::  GC_Sender1,GC_Receiver1
        integer                         ::  ierr,request,status(MPI_STATUS_SIZE)
        integer,allocatable             ::  requests(:)

        ! Loop all local Blocks%GC_sources to find all neighbour ranks.
        ! First set UniqueMpiRanks to be -1. Loop, if find any new neighbour,
        ! record it there.

        UniqueMpiRank_ptr=1
        UniqueMpiRanks=-1

        do iLocalBlock=1,Tree%nLocalBlocks
            Block_target=>Tree%LocalBlocks(iLocalBlock)
            do iGC_target=1,Block_target%nGC_targets
                GC_target1=>Block_target%GC_targets(iGC_target)

                if (GC_target1%iRank/=MpiRank) then
                    ! Loop all recorded unique ranks to see if it's already in there.
                    ! If not, record it.
                    do iUniqueMpiRank=1,UniqueMpiRank_ptr
                        if (UniqueMpiRanks(iUniqueMpiRank)==GC_target1%iRank) exit
                    end do
                    if (iUniqueMpiRank==UniqueMpiRank_ptr+1) then
                        UniqueMpiRanks(UniqueMpiRank_ptr)=GC_target1%iRank
                        UniqueMpiRank_ptr=UniqueMpiRank_ptr+1
                    end if
                end if
            end do
        end do

        ! Now the nPairs is UniqueMpiRank_ptr-1.
        ! Allocate the GC_Receivers.
        Tree%nGC_Receivers=UniqueMpiRank_ptr-1
        allocate(Tree%GC_Receivers(Tree%nGC_Receivers))
        Tree%GC_Receivers(:)%iRank=UniqueMpiRanks(1:Tree%nGC_Receivers)

        ! Get the iBlock_pairs for each GC_receiver.
        ! Loop all the local blocks and GC_sources to get the unique pairs.
        ! The process is pretty similar to above, but we need a 2D array 
        ! for recording since it's a pair of blocks.

        do iGC_receiver=1,Tree%nGC_Receivers

            GC_Receiver1=>Tree%GC_Receivers(iGC_receiver)
            GC_Receiver1%nBlock_Pairs=0
            GC_Receiver1%nGC=0

            do iLocalBlock=1,Tree%nLocalBlocks
                Block_target=>Tree%LocalBlocks(iLocalBlock)
                do iGC_target=1,Block_target%nGC_targets
                    GC_target1=>Block_target%GC_targets(iGC_target)

                    ! No need to record. Different GC_sources always correspond to
                    ! different block pairs.

                    if (GC_target1%iRank==GC_Receiver1%iRank) then
                        GC_Receiver1%nBlock_Pairs=GC_Receiver1%nBlock_Pairs+1
                        GC_Receiver1%nGC=GC_Receiver1%nGC+GC_target1%nGC
                    end if

                end do
            end do

            allocate(GC_Receiver1%message(GC_Receiver1%nGC,nvar))
            allocate(GC_Receiver1%iBlock_pairs(GC_Receiver1%nBlock_Pairs,2))
            allocate(GC_Receiver1%Block_Pairs_Ptrs(GC_Receiver1%nBlock_Pairs,2))
            allocate(GC_Receiver1%iLocalBlock_and_iGCtarget_list(GC_Receiver1%nBlock_Pairs,2))
        end do

        ! Now we want to determine the iBlock_pairs and the Block_Pairs_Ptrs.
        ! We will send them from GC_Senders.

        allocate(requests(Tree%nGC_Senders))

        ! First send iBlock_pairs.
        do iGC_sender=1,Tree%nGC_Senders
            GC_Sender1=>Tree%GC_Senders(iGC_sender)
            call MPI_ISEND(GC_Sender1%iBlock_pairs,2*GC_Sender1%nBlock_Pairs,&
            mpi_integer,GC_Sender1%iRank,1,MPI_COMM_WORLD,request,ierr)
            requests(iGC_sender)=request
        end do

        do iGC_receiver=1,Tree%nGC_Receivers
            GC_Receiver1=>Tree%GC_Receivers(iGC_receiver)
            call MPI_RECV(GC_Receiver1%iBlock_pairs,2*GC_Receiver1%nBlock_Pairs,&
            mpi_integer,GC_Receiver1%iRank,1,MPI_COMM_WORLD,status,ierr)
        end do

        call MPI_WAITALL(Tree%nGC_Senders,requests,MPI_STATUSES_IGNORE,ierr)

        ! Then send Block_Pairs_Ptrs.
        do iGC_sender=1,Tree%nGC_Senders
            GC_Sender1=>Tree%GC_Senders(iGC_sender)
            call MPI_ISEND(GC_Sender1%Block_Pairs_Ptrs,2*GC_Sender1%nBlock_Pairs,&
            mpi_integer,GC_Sender1%iRank,1,MPI_COMM_WORLD,request,ierr)
            requests(iGC_sender)=request
        end do

        do iGC_receiver=1,Tree%nGC_Receivers
            GC_Receiver1=>Tree%GC_Receivers(iGC_receiver)
            call MPI_RECV(GC_Receiver1%Block_Pairs_Ptrs,2*GC_Receiver1%nBlock_Pairs,&
            mpi_integer,GC_Receiver1%iRank,1,MPI_COMM_WORLD,status,ierr)
        end do

        call MPI_WAITALL(Tree%nGC_Senders,requests,MPI_STATUSES_IGNORE,ierr)


        deallocate(requests)

        ! Ok. Now one more thing to do is to get the iLocalBlock_and_iGCtarget_list.
        ! This is pretty easy: since we have iBlock pairs, the second column tells
        ! which block pair it is -- which can be used to get the iLocalBlock. Then 
        ! we want to loop all the GC_targets of this local block to get the iGCtarget.

        do iGC_receiver=1,Tree%nGC_Receivers
            GC_Receiver1=>Tree%GC_Receivers(iGC_receiver)
            GC_Receiver1%iLocalBlock_and_iGCtarget_list(:,1)=&
                GC_Receiver1%iBlock_pairs(:,2)-ranges_of_ranks(MpiRank,1)+1
            
            do iofBlock_pair=1,GC_Receiver1%nBlock_Pairs
                iLocalBlock=GC_Receiver1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,1)
                Block_target=>Tree%LocalBlocks(iLocalBlock)
                do iGC_target=1,Block_target%nGC_targets
                    GC_target1=>Block_target%GC_targets(iGC_target)
                    if (GC_target1%iBlock==GC_Receiver1%iBlock_pairs(iofBlock_pair,1)) then
                        GC_Receiver1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,2)=iGC_target
                        exit
                    end if
                end do
            end do
        end do
    end subroutine ModCommunication_SetGC_Receiver

    ! Ok the new one works. I want a new new version that
    ! calculates EQN_update_R during sendrecv and waitall.
    ! But here is a tricky thing: in the two previous versions,
    ! the local sendrecvs are done BEFORE sendrecv waitall, while
    ! the global sendrecvs are done after those. So here I should 
    ! first move all local sendrecvs to the last loop.
    ! Then the second thing is to add a new input like if_rk to
    ! determine how to calculate EQN_update_R. if_rk should be
    ! renamed as if_rk_communicate, and the new one should be
    ! called if_rk_calc.
    ! And there's one thing to do before all these: define the
    ! EQN_update_R as a variable in BlockType. 

    

    ! In this new version i want to follow the suggestion
    ! from chatgpt that put the mpi_irecv before the mpi_isend.
    ! Hopefully this one will work and be faster.

    subroutine ModCommunication_SendRecvGC_new(Tree,if_rk)

        implicit none
        type(YYTree),target             ::  Tree
        logical,intent(in)              ::  if_rk

        type(BlockType),pointer         ::  Block_source,Block_target
        type(GC_target),pointer         ::  GC_source1,GC_target1
        integer                         ::  iLocalBlock
        integer                         ::  iGC_source,iGC_target
        real(8),pointer                 ::  primitive_source_IV(:,:,:,:)
        real(8),pointer                 ::  primitive_target_IV(:,:,:,:)

        integer                         ::  iGC
        integer                         ::  iGC_sender,iGC_receiver
        type(GC_interface),pointer      ::  GC_Sender1,GC_Receiver1
        integer                         ::  iofBlock_pair
        integer                         ::  ierr
        integer,allocatable             ::  requests_send(:),requests_recv(:)

        ! Before getting the message in the GC_senders we need to interpolate the
        ! variables into each GC_sources.

        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)

            if (if_rk) then
                primitive_source_IV=>Block_source%primitive_rk_IV
            else
                primitive_source_IV=>Block_source%primitive_IV
            end if

            do iGC_source=1,Block_source%nGC_sources
                GC_source1=>Block_source%GC_sources(iGC_source)
                
                ! Interpolate the primitive variables.
                call ModMath_1D3D_interpolate_1D1D(primitive_source_IV,nvar,ni,nj,nk,ng,&
                    Block_source%xi_I,Block_source%xj_I,Block_source%xk_I,GC_source1%nGC,&
                    GC_source1%xijk_list,GC_source1%primitive_list)
                
                if (GC_source1%if_yin .neqv. Block_source%if_yin) then
                    if (Block_source%br_>0) GC_source1%primitive_list(:,Block_source%br_:Block_source%bp_)=&
                        ModYinYang_VecConv_1D(GC_source1%xijk_list,&
                        GC_source1%primitive_list(:,Block_source%br_:Block_source%bp_),GC_source1%nGC)
                    if (Block_source%vr_>0) GC_source1%primitive_list(:,Block_source%vr_:Block_source%vp_)=&
                        ModYinYang_VecConv_1D(GC_source1%xijk_list,&
                       GC_source1%primitive_list(:,Block_source%vr_:Block_source%vp_),GC_source1%nGC)
                end if

                if (GC_source1%iRank==MpiRank) then
                    Block_target=>Tree%LocalBlocks(GC_source1%iBlock-ranges_of_ranks(MpiRank,1)+1)

                    if (if_rk) then
                        primitive_target_IV=>Block_target%primitive_rk_IV
                    else
                        primitive_target_IV=>Block_target%primitive_IV
                    end if

                    do iGC=1,GC_source1%nGC
                        primitive_target_IV(GC_source1%ijk_list(iGC,1),&
                            GC_source1%ijk_list(iGC,2),&
                            GC_source1%ijk_list(iGC,3),:)=&
                            GC_source1%primitive_list(iGC,:)
                    end do
                end if
            end do
        end do

        ! Post the Irecv before the Isend.

        allocate(requests_recv(Tree%nGC_Receivers))

        do iGC_receiver=1,Tree%nGC_Receivers
            ! Receive the message.
            GC_Receiver1=>Tree%GC_Receivers(iGC_receiver)
            call MPI_IRECV(GC_Receiver1%message,nvar*GC_Receiver1%nGC,mpi_double_precision,&
            GC_Receiver1%iRank,1,MPI_COMM_WORLD,requests_recv(iGC_receiver),ierr)

            do iofBlock_pair=1,GC_Receiver1%nBlock_Pairs
                GC_target1=>Tree%LocalBlocks(GC_Receiver1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,1))%&
                    GC_targets(GC_Receiver1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,2))
                GC_target1%primitive_list=&
                    GC_Receiver1%message(&
                    GC_Receiver1%Block_Pairs_Ptrs(iofBlock_pair,1):&
                    GC_Receiver1%Block_Pairs_Ptrs(iofBlock_pair,2),:)
            end do
        end do

        ! Now get the message and send them.
        allocate(requests_send(Tree%nGC_Senders))

        do iGC_sender=1,Tree%nGC_Senders
            GC_Sender1=>Tree%GC_Senders(iGC_sender)

            ! First get the message.
            do iofBlock_pair=1,GC_Sender1%nBlock_Pairs
                GC_source1=>Tree%LocalBlocks(GC_Sender1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,1))%&
                    GC_sources(GC_Sender1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,2))
                GC_Sender1%message(&
                    GC_Sender1%Block_Pairs_Ptrs(iofBlock_pair,1):&
                    GC_Sender1%Block_Pairs_Ptrs(iofBlock_pair,2),:)=&
                    GC_source1%primitive_list
            end do

            ! Send the message.
            call MPI_ISEND(GC_Sender1%message,nvar*GC_Sender1%nGC,mpi_double_precision,&
                GC_Sender1%iRank,1,MPI_COMM_WORLD,requests_send(iGC_sender),ierr)
        end do

        call MPI_WAITALL(Tree%nGC_Receivers,requests_recv,MPI_STATUSES_IGNORE,ierr)
        call MPI_WAITALL(Tree%nGC_Senders,requests_send,MPI_STATUSES_IGNORE,ierr)

        ! Assign the GC_targets' primitive_list to primitive blocks.
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_target=>Tree%LocalBlocks(iLocalBlock)

            if (if_rk) then
                primitive_target_IV=>Block_target%primitive_rk_IV
            else
                primitive_target_IV=>Block_target%primitive_IV
            end if

            do iGC_target=1,Block_target%nGC_targets
                GC_target1=>Block_target%GC_targets(iGC_target)

                if (GC_target1%iRank/=MpiRank) then
                    do iGC=1,GC_target1%nGC
                        primitive_target_IV(GC_target1%ijk_list(iGC,1),&
                            GC_target1%ijk_list(iGC,2),&
                            GC_target1%ijk_list(iGC,3),:)=&
                            GC_target1%primitive_list(iGC,:)
                    end do
                end if
            end do
        end do

        deallocate(requests_send)
        deallocate(requests_recv)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    end subroutine ModCommunication_SendRecvGC_new

    ! With all above we can send and receive the messages
    ! using G(H)C_Senders and G(H)C_Receivers.
    ! This one for GC. GC has no MG levels.

    subroutine ModCommunication_SendRecvGC(Tree,if_rk)
        implicit none
        type(YYTree),target             ::  Tree
        logical,intent(in)              ::  if_rk

        type(BlockType),pointer         ::  Block_source,Block_target
        type(GC_target),pointer         ::  GC_source1,GC_target1
        integer                         ::  iLocalBlock
        integer                         ::  iGC_source,iGC_target
        real(8),pointer                 ::  primitive_source_IV(:,:,:,:)
        real(8),pointer                 ::  primitive_target_IV(:,:,:,:)

        integer                         ::  iGC
        integer                         ::  iGC_sender,iGC_receiver
        type(GC_interface),pointer      ::  GC_Sender1,GC_Receiver1
        integer                         ::  iofBlock_pair
        integer                         ::  ierr,request,status(MPI_STATUS_SIZE)
        integer,allocatable             ::  requests(:)

        ! Before getting the message in the GC_senders we need to interpolate the
        ! variables into each GC_sources.

        do iLocalBlock=1,Tree%nLocalBlocks
            Block_source=>Tree%LocalBlocks(iLocalBlock)

            if (if_rk) then
                primitive_source_IV=>Block_source%primitive_rk_IV
            else
                primitive_source_IV=>Block_source%primitive_IV
            end if

            do iGC_source=1,Block_source%nGC_sources
                GC_source1=>Block_source%GC_sources(iGC_source)
                
                ! Interpolate the primitive variables.
                call ModMath_1D3D_interpolate_1D1D(primitive_source_IV,nvar,ni,nj,nk,ng,&
                    Block_source%xi_I,Block_source%xj_I,Block_source%xk_I,GC_source1%nGC,&
                    GC_source1%xijk_list,GC_source1%primitive_list)
                
                if (GC_source1%if_yin .neqv. Block_source%if_yin) then
                    if (Block_source%br_>0) GC_source1%primitive_list(:,Block_source%br_:Block_source%bp_)=&
                        ModYinYang_VecConv_1D(GC_source1%xijk_list,&
                        GC_source1%primitive_list(:,Block_source%br_:Block_source%bp_),GC_source1%nGC)
                    if (Block_source%vr_>0) GC_source1%primitive_list(:,Block_source%vr_:Block_source%vp_)=&
                        ModYinYang_VecConv_1D(GC_source1%xijk_list,&
                       GC_source1%primitive_list(:,Block_source%vr_:Block_source%vp_),GC_source1%nGC)
                end if

                if (GC_source1%iRank==MpiRank) then
                    Block_target=>Tree%LocalBlocks(GC_source1%iBlock-ranges_of_ranks(MpiRank,1)+1)

                    if (if_rk) then
                        primitive_target_IV=>Block_target%primitive_rk_IV
                    else
                        primitive_target_IV=>Block_target%primitive_IV
                    end if

                    do iGC=1,GC_source1%nGC
                        primitive_target_IV(GC_source1%ijk_list(iGC,1),&
                            GC_source1%ijk_list(iGC,2),&
                            GC_source1%ijk_list(iGC,3),:)=&
                            GC_source1%primitive_list(iGC,:)
                    end do
                end if
            end do
        end do

        ! Now get the message and send them.
        allocate(requests(Tree%nGC_Senders))

        do iGC_sender=1,Tree%nGC_Senders
            GC_Sender1=>Tree%GC_Senders(iGC_sender)

            ! First get the message.
            do iofBlock_pair=1,GC_Sender1%nBlock_Pairs
                GC_source1=>Tree%LocalBlocks(GC_Sender1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,1))%&
                    GC_sources(GC_Sender1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,2))
                GC_Sender1%message(&
                    GC_Sender1%Block_Pairs_Ptrs(iofBlock_pair,1):&
                    GC_Sender1%Block_Pairs_Ptrs(iofBlock_pair,2),:)=&
                    GC_source1%primitive_list
            end do

            ! Send the message.
            call MPI_ISEND(GC_Sender1%message,nvar*GC_Sender1%nGC,mpi_double_precision,&
            GC_Sender1%iRank,1,MPI_COMM_WORLD,request,ierr)
            requests(iGC_sender)=request
        end do

        do iGC_receiver=1,Tree%nGC_Receivers
            ! Receive the message.
            GC_Receiver1=>Tree%GC_Receivers(iGC_receiver)
            call MPI_RECV(GC_Receiver1%message,nvar*GC_Receiver1%nGC,mpi_double_precision,&
                GC_Receiver1%iRank,1,MPI_COMM_WORLD,status,ierr)

            do iofBlock_pair=1,GC_Receiver1%nBlock_Pairs
                GC_target1=>Tree%LocalBlocks(GC_Receiver1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,1))%&
                    GC_targets(GC_Receiver1%iLocalBlock_and_iGCtarget_list(iofBlock_pair,2))
                GC_target1%primitive_list=&
                    GC_Receiver1%message(&
                    GC_Receiver1%Block_Pairs_Ptrs(iofBlock_pair,1):&
                    GC_Receiver1%Block_Pairs_Ptrs(iofBlock_pair,2),:)
            end do
        end do

        call MPI_WAITALL(Tree%nGC_Senders,requests,MPI_STATUSES_IGNORE,ierr)

        ! Assign the GC_targets' primitive_list to primitive blocks.
        do iLocalBlock=1,Tree%nLocalBlocks
            Block_target=>Tree%LocalBlocks(iLocalBlock)

            if (if_rk) then
                primitive_target_IV=>Block_target%primitive_rk_IV
            else
                primitive_target_IV=>Block_target%primitive_IV
            end if

            do iGC_target=1,Block_target%nGC_targets
                GC_target1=>Block_target%GC_targets(iGC_target)

                if (GC_target1%iRank/=MpiRank) then
                    do iGC=1,GC_target1%nGC
                        primitive_target_IV(GC_target1%ijk_list(iGC,1),&
                            GC_target1%ijk_list(iGC,2),&
                            GC_target1%ijk_list(iGC,3),:)=&
                            GC_target1%primitive_list(iGC,:)
                    end do
                end if
            end do
        end do

        deallocate(requests)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
    end subroutine ModCommunication_SendRecvGC

    ! About the tags.

    function ModCommunication_GetTag(iBlockSend,iBlockRecv,iRankSend,iRankRecv) result(Tag)
        implicit none
        integer,intent(in)      :: iBlockSend,iBlockRecv
        integer,intent(in)      :: iRankSend,iRankRecv 
        
        integer                 :: tag

        Tag=nMaxBlocksPerRank*(iBlockSend-ranges_of_ranks(iRankSend,1))+&
            (iBlockRecv-ranges_of_ranks(iRankRecv,1))
    end function ModCommunication_GetTag

    function ModCommunication_SolveTag(Tag,iRankSend,iRankRecv) result(iBlockSendRecv)
        implicit none
        integer,intent(in)      :: Tag
        integer,intent(in)      :: iRankRecv,iRankSend

        integer                 :: iBlockSendRecv(2)

        iBlockSendRecv(1)=Tag/nMaxBlocksPerRank+ranges_of_ranks(iRankSend,1)
        iBlockSendRecv(2)=mod(Tag,nMaxBlocksPerRank)+ranges_of_ranks(iRankRecv,1)
    end function ModCommunication_SolveTag
end module