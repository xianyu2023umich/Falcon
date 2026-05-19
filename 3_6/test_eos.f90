program test_eos

    use ModLookUpTable, only: ModLookUpTable_Read, LookUpTables, nLookUpTables
    use ModEOS,         only: ModEOS_init, EOS_table, nlogQ_EOS_table, nlogT_EOS_table,&
                              logQ_EOS_table, logT_EOS_table

    implicit none

    character(len=*), parameter :: eos_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/FreeEOS_v51_x0.75_z0.02.dat'

    call ModLookUpTable_Read(eos_file, logical_unit=1)
    call ModEOS_init

    call assert_true(nLookUpTables == 1, 'expected one lookup table')
    call assert_true(trim(LookUpTables(1)%name) == 'EOS', 'unexpected table name')
    call assert_true(associated(EOS_table, LookUpTables(1)), 'EOS_table is not associated with the EOS lookup table')
    call assert_true(LookUpTables(1)%dimension == 3, 'EOS lookup table is not 3D')
    call assert_true(all(LookUpTables(1)%sizes == [527, 306, 23]), 'unexpected EOS table sizes')
    call assert_true(nlogQ_EOS_table == 527, 'nlogQ_EOS_table mismatch')
    call assert_true(nlogT_EOS_table == 306, 'nlogT_EOS_table mismatch')
    call assert_true(maxval(abs(logQ_EOS_table - LookUpTables(1)%data_3D(:,1,1))) <= 0d0, 'logQ_EOS_table copy mismatch')
    call assert_true(maxval(abs(logT_EOS_table - LookUpTables(1)%data_3D(1,:,2))) <= 0d0, 'logT_EOS_table copy mismatch')
    call assert_true(all(logQ_EOS_table(2:) > logQ_EOS_table(:size(logQ_EOS_table)-1)), 'logQ grid is not strictly increasing')
    call assert_true(all(logT_EOS_table(2:) > logT_EOS_table(:size(logT_EOS_table)-1)), 'logT grid is not strictly increasing')

    print *, 'ModEOS regression test passed.'

contains

    subroutine assert_true(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            print *, 'ModEOS regression test failed: ', trim(message)
            stop 1
        end if
    end subroutine assert_true

end program test_eos