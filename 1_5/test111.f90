program test111

    use mpi
    implicit none

    integer :: ierr, flag, source, tag, comm
    integer :: status(MPI_STATUS_SIZE)

    call MPI_INIT(ierr)
    comm = MPI_COMM_WORLD

    ! Probe for messages
    call MPI_IPROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, flag, status, ierr)

    if (flag .ne. 0) then
        print *, "Message is available"
    else
        print *, "No message available"
    end if

    call MPI_FINALIZE(ierr)


end program test111