program test_stratification_new

    use ModLookUpTable,       only: LookUpTables, nLookUpTables
    use ModEOS,               only: ModEOS_init
    use ModStratification_new, only: ModStratification_new_get_pressure_density,&
                                     nEntropyContour, logQ_T_of_entropy_contour,&
                                     logP_of_entropy_contour

    implicit none

    real(8), parameter :: tol = 1d-12
    real(8) :: expected(2)

    call build_synthetic_eos_table
    call prepare_queries

    call ModEOS_init

    call ModStratification_new_get_pressure_density

    expected = [0.50d0, 1.00d0]

    call assert_true(nLookUpTables == 1, 'expected one synthetic lookup table')
    call assert_true(nEntropyContour == 2, 'expected two synthetic query points')
    call assert_true(allocated(logP_of_entropy_contour), 'logP_of_entropy_contour not allocated')
    call assert_true(all(abs(logP_of_entropy_contour - expected) <= tol), 'interpolated values do not match the synthetic expectation')

    print *, 'ModStratification_new regression test passed.'

contains

    subroutine build_synthetic_eos_table
        integer :: iPoint, jPoint

        nLookUpTables = 1
        LookUpTables(1)%name = 'EOS'
        LookUpTables(1)%dimension = 3
        allocate(LookUpTables(1)%sizes(3))
        LookUpTables(1)%sizes = [2, 2, 3]
        allocate(LookUpTables(1)%data_3D(2,2,3))

        do jPoint = 1, 2
            do iPoint = 1, 2
                LookUpTables(1)%data_3D(iPoint,jPoint,1) = real(iPoint - 1, 8)
                LookUpTables(1)%data_3D(iPoint,jPoint,2) = real(jPoint - 1, 8)
                LookUpTables(1)%data_3D(iPoint,jPoint,3) = LookUpTables(1)%data_3D(iPoint,jPoint,1) + &
                                                         0.5d0 * LookUpTables(1)%data_3D(iPoint,jPoint,2)
            end do
        end do
    end subroutine build_synthetic_eos_table

    subroutine prepare_queries
        nEntropyContour = 2
        allocate(logQ_T_of_entropy_contour(nEntropyContour,2))
        allocate(logP_of_entropy_contour(nEntropyContour))

        logQ_T_of_entropy_contour(1,:) = [0.25d0, 0.5d0]
        logQ_T_of_entropy_contour(2,:) = [0.75d0, 0.5d0]
    end subroutine prepare_queries

    subroutine assert_true(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            print *, 'ModStratification_new regression test failed: ', trim(message)
            stop 1
        end if
    end subroutine assert_true

end program test_stratification_new