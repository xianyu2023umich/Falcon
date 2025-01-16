Module ModLookUpTable



    contains

    subroutine ModLookUpTable_readfile(filename,logical_unit,head_size,data)
        implicit none
        character(len=*),intent(in)     ::          filename        !   filename
        integer,intent(in)              ::          logical_unit    !   logical unit
        integer,intent(in)              ::          head_size       !   n lines of head
        real,pointer,intent(inout)      ::          data(:,:)       !   store the data

        character                       ::          header
        integer                         ::          ihead
        integer                         ::          iline
        integer                         ::          nvars,nsamples  !   n of vars and samples

        open(unit=logical_unit, file=filename, status='old', action='read')

        ! Read the head
        do ihead=1,head_size
            read(logical_unit, '(A)') header
        end do

        ! Read data size
        read(logical_unit,*) nvars,nsamples
        allocate(data(nvars,nsamples))

        do iline=1,nsamples
            read(logical_unit,*) data(:,iline)
        end do




    end subroutine ModLookUpTable_readfile



end module ModLookUpTable