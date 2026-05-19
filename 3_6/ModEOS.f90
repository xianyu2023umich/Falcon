module ModEOS

    use ModLookUpTable, only: LookUpTables,LookUpTable,nLookUpTables

    implicit none

    type(LookUpTable),pointer   ::  EOS_table
    integer                     ::  nlogQ_EOS_table,nlogT_EOS_table
    real(8),allocatable         ::  logQ_EOS_table(:)
    real(8),allocatable         ::  logT_EOS_table(:)

    contains

    subroutine ModEOS_init
        implicit none
        integer ::  iLookUpTable

        do iLookUpTable=1,nLookUpTables
            
            if (LookUpTables(iLookUpTable)%name == 'EOS') then
                EOS_table => LookUpTables(iLookUpTable)
                nlogQ_EOS_table = EOS_table%sizes(1)
                nlogT_EOS_table = EOS_table%sizes(2)

                ! Get logQ and logT grids of the EOS table.
                allocate(logQ_EOS_table(nlogQ_EOS_table))
                allocate(logT_EOS_table(nlogT_EOS_table))
                logQ_EOS_table = EOS_table%data_3D(:,1,1)
                logT_EOS_table = EOS_table%data_3D(1,:,2)
            end if
        end do

    end subroutine ModEOS_init
end module