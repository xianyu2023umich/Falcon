program test0

    ! This program is used to test the speed of the two versions of the rebinning function.
    ! Receive an input from terminal as n_test, which is the number of tests.
    use ModMath
    implicit none

    real,allocatable :: f(:,:,:)
    real,allocatable :: f_prolongate_version1(:,:,:),f_prolongate_version2(:,:,:)

    integer :: ni,nj,nk
    integer :: i,j,k
    integer :: n_test,i_test
    real    :: t1,t2
    character(len=32) :: arg
    integer :: status

    ! Read command line argument
    if (command_argument_count() < 1    ) then
        print *, "Usage: ./test0.exe <n_test>"
        print *, "Example: ./test0.exe 10"
        stop
    end if
    
    call get_command_argument(1, arg, status=status)
    if (status /= 0) then
        print *, "Error reading command line argument"
        stop
    end if
    
    read(arg, *, iostat=status) n_test
    if (status /= 0) then
        print *, "Error: n_test must be an integer"
        stop
    end if

    print *, "n_test = ", n_test
    ! This is a realistic case.
    ni=30; nj=30; nk=90

    allocate(f(ni,nj,nk))
    allocate(f_prolongate_version1(ni*2,nj*2,nk*2))
    allocate(f_prolongate_version2(ni*2,nj*2,nk*2))

    do k=1,nk
        do j=1,nj
            do i=1,ni
                f(i,j,k)=i+j+k
            end do
        end do
    end do

    f_prolongate_version1=ModMath_prolongate_3D_factor2_version1(f,ni,nj,nk)
    f_prolongate_version2=ModMath_prolongate_3D_factor2_version2(f,ni,nj,nk)

    print *,f_prolongate_version1(1,:,:)-f_prolongate_version2(1,:,:) 

    ! Now test speed
    call cpu_time(t1)
    do i_test=1,n_test
        f_prolongate_version1=ModMath_prolongate_3D_factor2_version1(f,ni,nj,nk)
    end do
    call cpu_time(t2)
    print *, "Time for version 1: ", t2-t1



    call cpu_time(t1)
    do i_test=1,n_test
        f_prolongate_version2=ModMath_prolongate_3D_factor2_version2(f,ni,nj,nk)
    end do
    call cpu_time(t2)
    print *, "Time for version 2: ", t2-t1

    deallocate(f)
    deallocate(f_prolongate_version1)
    deallocate(f_prolongate_version2)
end program test0