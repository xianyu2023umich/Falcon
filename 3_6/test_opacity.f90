program test_opacity

    use ModLookUpTable, only: ModLookUpTable_Read, LookUpTables, nLookUpTables
    use ModOpacity,     only: ModOpacity_init, Opacity_table, nlogT_opacity, nlogR_opacity,&
                              logT_opacity, logR_opacity

    implicit none

    character(len=*), parameter :: opacity_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/opacity_x70_y28_z2.dat'
    real(8), parameter :: tol = 0d0

    call ModLookUpTable_Read(opacity_file, logical_unit=1)
    call ModOpacity_init

    call assert_true(nLookUpTables == 1, 'expected one lookup table')
    call assert_true(trim(LookUpTables(1)%name) == 'OPACITY', 'unexpected table name')
    call assert_true(LookUpTables(1)%dimension == 3, 'opacity table is not 3D')
    call assert_true(all(LookUpTables(1)%sizes == [70, 19, 3]), 'unexpected opacity table sizes')
    call assert_true(associated(Opacity_table, LookUpTables(1)), 'Opacity_table is not associated with the loaded table')
    call assert_true(nlogT_opacity == 70, 'nlogT_opacity mismatch')
    call assert_true(nlogR_opacity == 19, 'nlogR_opacity mismatch')
    call assert_true(maxval(abs(logT_opacity - LookUpTables(1)%data_3D(:,1,1))) <= tol, 'logT_opacity copy mismatch')
    call assert_true(maxval(abs(logR_opacity - LookUpTables(1)%data_3D(1,:,2))) <= tol, 'logR_opacity copy mismatch')
    call assert_true(all(logT_opacity(2:) >= logT_opacity(:size(logT_opacity)-1)), 'logT grid is not monotonic')
    call assert_true(all(logR_opacity(2:) >= logR_opacity(:size(logR_opacity)-1)), 'logR grid is not monotonic')

    print *, 'ModOpacity regression test passed.'

contains

    subroutine assert_true(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            print *, 'ModOpacity regression test failed: ', trim(message)
            stop 1
        end if
    end subroutine assert_true

end program test_opacity