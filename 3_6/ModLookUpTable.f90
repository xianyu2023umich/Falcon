module ModLookUpTable

    implicit none
    
    type LookUpTable
        character(len=256)  ::  name
        integer             ::  iLookUpTable
        integer             ::  dimension
        integer,allocatable ::  sizes(:)

        ! One of the following is used depending on the dimension:  
        real(8),allocatable ::  data_1D(:)
        real(8),allocatable ::  data_2D(:,:)
        real(8),allocatable ::  data_3D(:,:,:)
    end type LookUpTable

    integer                         ::  nLookUpTables=0
    integer                         ::  iLookUpTable=0
    type(LookUpTable),target        ::  LookUpTables(8)     ! Maximum 8 dimensions
    type(LookUpTable),pointer       ::  LookUpTable_ptr

    contains

    ! This subroutine reads any lookuptable with the standard format:
    ! first line: name
    ! second line: dimension
    ! third line: sizes(dimension)
    ! fourth line: the var names. I don't read them, just use the order.
    ! from fifth line: the data.

    subroutine ModLookUpTable_Read(filename,logical_unit)
        implicit none
        character(len=*),intent(in) ::  filename
        integer,intent(in)          ::  logical_unit

        nLookUpTables=nLookUpTables+1
        iLookUpTable=iLookUpTable+1
        LookUpTable_ptr=>LookUpTables(iLookUpTable)

        open(unit=logical_unit, file=filename, status='old', action='read')

        ! Read the name
        read(logical_unit,*) LookUpTable_ptr%name

        ! Read the dimension
        read(logical_unit,*) LookUpTable_ptr%dimension

        ! Read the sizes
        allocate(LookUpTable_ptr%sizes(LookUpTable_ptr%dimension))
        read(logical_unit,*) LookUpTable_ptr%sizes

        ! Then skip the var names
        read(logical_unit,*)

        ! Read the data
        select case (LookUpTable_ptr%dimension)
        case (1)
            call ModLookUpTable_Read_1D(logical_unit)
        case (2)
            call ModLookUpTable_Read_2D(logical_unit)
        case (3)
            call ModLookUpTable_Read_3D(logical_unit)
        end select

        close(logical_unit)
    end subroutine ModLookUpTable_Read

    ! 1D case only has lines each with one value.

    subroutine ModLookUpTable_Read_1D(logical_unit)
        implicit none
        integer,intent(in)              ::  logical_unit
        integer                         ::  iLine

        allocate(LookUpTable_ptr%data_1D(LookUpTable_ptr%sizes(1)))
        do iLine=1,LookUpTable_ptr%sizes(1)
            read(logical_unit,*) LookUpTable_ptr%data_1D(iLine)
        end do
    end subroutine ModLookUpTable_Read_1D

    ! 2D bit different: each line has multiple values given by sizes(2).

    subroutine ModLookUpTable_Read_2D(logical_unit)
        implicit none
        integer,intent(in)              ::  logical_unit
        integer                         ::  iLine

        allocate(LookUpTable_ptr%data_2D(LookUpTable_ptr%sizes(1),LookUpTable_ptr%sizes(2)))
        do iLine=1,LookUpTable_ptr%sizes(1)
            read(logical_unit,*) LookUpTable_ptr%data_2D(iLine,:)
        end do
    end subroutine ModLookUpTable_Read_2D

    ! 3D: size(1)*size(2) number of lines each with sizes(3) values.
    ! size(1) is the outer loop and size(2) is the inner loop.

    subroutine ModLookUpTable_Read_3D(logical_unit)
        implicit none
        integer,intent(in)              ::  logical_unit
        integer                         ::  iLine,jLine

        allocate(LookUpTable_ptr%data_3D(LookUpTable_ptr%sizes(1),LookUpTable_ptr%sizes(2),LookUpTable_ptr%sizes(3)))
        do iLine=1,LookUpTable_ptr%sizes(1)
            do jLine=1,LookUpTable_ptr%sizes(2)
                read(logical_unit,*) LookUpTable_ptr%data_3D(iLine,jLine,:)
            end do
        end do
    end subroutine ModLookUpTable_Read_3D

end module ModLookUpTable