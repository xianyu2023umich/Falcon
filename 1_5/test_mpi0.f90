program test_mpi0

    !use ModAllocation
    !use mpi

    integer :: ierr,size,rank
    integer :: data,destination,tag
    integer :: status(MPI_STATUS_SIZE)
    integer :: i,j
    real :: k

    integer,allocatable :: iNodes(:)

    call mpi_init(ierr)

    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    ! Get the rank of the process
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    !print *,'Hello World!'
    !print *,size,rank,ierr

    if (rank==0) then

        do destination=1,size-1
            tag=destination
            data=destination ** 2
            call mpi_send(data,1,mpi_integer,destination,tag,MPI_COMM_WORLD, ierr)
        end do
        print *,'Sent'

    else
        tag=rank
        call mpi_recv(data,1,mpi_integer,0,tag,MPI_COMM_WORLD,status,ierr)

        print *,'Received, data=',data,"I'm ",rank
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    call sleep(0.001)

    call ModAllocation_GetNodes(11,iNodes,Size,Rank)

    print *,'I"m ',rank,', my nodes are',iNodes

    

    call mpi_finalize(ierr)



end program test_mpi0