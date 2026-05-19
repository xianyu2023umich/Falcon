program test_ModBlock

    use ModConst,              only: dpi
    use ModParameters,         only: ni, nj, nk, ng, iGeometry, iEquation
    use ModLookUpTable,        only: ModLookUpTable_Read
    use ModEOS,                only: ModEOS_init
    use ModOpacity,            only: ModOpacity_init
    use ModStratification_new, only: ModStratification_new_Init, ModStratification_new_get_vars
    use ModBlock,              only: BlockType, ModBlock_Init

    implicit none

    character(len=*), parameter :: eos_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/FreeEOS_v51_x0.75_z0.02.dat'
    character(len=*), parameter :: entropy_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/FreeEOS_v51_x0.75_z0.02_entropy_contour_9.252752.dat'
    character(len=*), parameter :: opacity_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/opacity_x70_y28_z2.dat'

    type(BlockType), target :: block
    real(8) :: xijk_range(3,2)
    ni = 12
    nj = 1
    nk = 1
    ng = 2
    iGeometry = 0
    iEquation = 0

    call ModLookUpTable_Read(eos_file, logical_unit=1)
    call ModLookUpTable_Read(entropy_file, logical_unit=1)
    call ModLookUpTable_Read(opacity_file, logical_unit=1)

    call ModEOS_init
    call ModOpacity_init
    call ModStratification_new_Init

    xijk_range(1,:) = [0.70d0, 1.00d0]
    xijk_range(2,:) = [0.0d0, 1.0d0]
    xijk_range(3,:) = [0.0d0, 1.0d0]

    call ModBlock_Init(block, 1, xijk_range, .true., iEquation, .true.)

    call assert_equal_int(block%nvar, 5, 'nvar')
    call assert_close(block%xi_I(1), xijk_range(1,1) + 0.5d0 * (xijk_range(1,2) - xijk_range(1,1)) / real(ni,8), 1.0d-12, 'first xi cell center')

    call check_block_point(1)
    call check_block_point(max(1, ni/2))
    call check_block_point(ni)

    call check_3d_replication(1)
    call check_3d_replication(max(1, ni/2))
    call check_3d_replication(ni)

    print *, 'ModBlock regression test passed.'

contains

    subroutine check_block_point(iPoint)
        integer, intent(in) :: iPoint
        real(8) :: g2, rho2, p2, T2, gamma12, gamma32, kap2
        real(8) :: diffusion2, cooling2, diffusion_flux2, cooling_flux2, Xi2

        call ModStratification_new_get_vars(block%xi_I(iPoint), g2, rho2, p2, T2, gamma12, gamma32, kap2, &
            diffusion2, cooling2, diffusion_flux2, cooling_flux2, Xi2)

        print *, 'Sample point index =', iPoint
        print '(1x,a,es14.5)', '  r = ', block%xi_I(iPoint)
        print '(1x,a,es14.5,a,es14.5)', '  g   = ', block%g_I(iPoint), '  direct = ', g2
        print '(1x,a,es14.5,a,es14.5)', '  p0  = ', block%p0_I(iPoint), '  direct = ', p2
        print '(1x,a,es14.5,a,es14.5)', '  rho = ', block%rho0_I(iPoint), '  direct = ', rho2
        print '(1x,a,es14.5,a,es14.5)', '  T   = ', block%te0_I(iPoint), '  direct = ', T2
        print '(1x,a,es14.5,a,es14.5)', '  gamma1 = ', block%gamma1_I(iPoint), '  direct = ', gamma12
        print '(1x,a,es14.5,a,es14.5)', '  gamma3 = ', block%gamma3_I(iPoint), '  direct = ', gamma32
        print '(1x,a,es14.5,a,es14.5)', '  Xi  = ', block%Xi_rsst_I(iPoint), '  direct = ', Xi2
        print '(1x,a,es14.5,a,es14.5)', '  diffusion = ', block%diffusion_I(iPoint), '  direct = ', diffusion2
        print '(1x,a,es14.5,a,es14.5)', '  cooling = ', block%cooling_I(iPoint), '  direct = ', cooling2
        print '(1x,a,es14.5,a,es14.5)', '  flux diff = ', block%diffusion_flux_I(iPoint), '  direct = ', diffusion_flux2
        print '(1x,a,es14.5,a,es14.5)', '  flux cool = ', block%cooling_flux_I(iPoint), '  direct = ', cooling_flux2

        call assert_close(block%g_I(iPoint), g2, 1.0d-12, 'gravity match')
        call assert_close(block%p0_I(iPoint), p2, 1.0d-12, 'pressure match')
        call assert_close(block%rho0_I(iPoint), rho2, 1.0d-12, 'density match')
        call assert_close(block%te0_I(iPoint), T2, 1.0d-12, 'temperature match')
        call assert_close(block%gamma1_I(iPoint), gamma12, 1.0d-12, 'gamma1 match')
        call assert_close(block%gamma3_I(iPoint), gamma32, 1.0d-12, 'gamma3 match')
        call assert_close(block%diffusion_I(iPoint), diffusion2, 1.0d-12, 'diffusion match')
        call assert_close(block%cooling_I(iPoint), cooling2, 1.0d-12, 'cooling match')
        call assert_close(block%diffusion_flux_I(iPoint), diffusion_flux2, 1.0d-12, 'diffusion flux match')
        call assert_close(block%cooling_flux_I(iPoint), cooling_flux2, 1.0d-12, 'cooling flux match')
        call assert_close(block%Xi_rsst_I(iPoint), Xi2, 1.0d-12, 'Xi match')
    end subroutine check_block_point

    subroutine check_3d_replication(iPoint)
        integer, intent(in) :: iPoint
        integer :: jPoint, kPoint

        do jPoint = 1, nj
            do kPoint = 1, nk
                call assert_close(block%gamma1_III(iPoint,jPoint,kPoint), block%gamma1_I(iPoint), 1.0d-12, 'gamma1 3D replication')
                call assert_close(block%rho0_III(iPoint,jPoint,kPoint), block%rho0_I(iPoint), 1.0d-12, 'rho0 3D replication')
                call assert_close(block%g_over_rho0_III(iPoint,jPoint,kPoint), block%g_I(iPoint)/block%rho0_I(iPoint), 1.0d-12, 'g/rho0 3D replication')
                call assert_close(block%p0_over_rho0_III(iPoint,jPoint,kPoint), block%p0_I(iPoint)/block%rho0_I(iPoint), 1.0d-12, 'p0/rho0 3D replication')
                call assert_close(block%rho0T0_III(iPoint,jPoint,kPoint), block%rho0_I(iPoint)*block%te0_I(iPoint), 1.0d-12, 'rho0T0 3D replication')
                call assert_close(block%total_heat_III(iPoint,jPoint,kPoint), block%diffusion_I(iPoint)+block%cooling_I(iPoint), 1.0d-12, 'total heat 3D replication')
                call assert_close(block%Xi_rsst_III(iPoint,jPoint,kPoint), block%Xi_rsst_I(iPoint), 1.0d-12, 'Xi 3D replication')
            end do
        end do
    end subroutine check_3d_replication

    subroutine assert_equal_int(actual, expected, label)
        integer, intent(in) :: actual, expected
        character(len=*), intent(in) :: label

        if (actual /= expected) then
            print *, 'ASSERTION FAILED:', trim(label), ' expected=', expected, ' actual=', actual
            stop 1
        end if
    end subroutine assert_equal_int

    subroutine assert_close(actual, expected, tolerance, label)
        real(8), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tolerance * max(1.0d0, abs(expected))) then
            print *, 'ASSERTION FAILED:', trim(label)
            print *, '  actual  =', actual
            print *, '  expected=', expected
            print *, '  abs diff=', abs(actual - expected)
            stop 1
        end if
    end subroutine assert_close

end program test_ModBlock
