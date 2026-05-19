program test_ModBlock_integral

    use ModConst,              only: dpi, R_sun__CGS
    use ModLookUpTable,        only: ModLookUpTable_Read
    use ModEOS,                only: ModEOS_init
    use ModOpacity,            only: ModOpacity_init
    use ModStratification_new, only: ModStratification_new_Init
    use ModParameters,         only: ni, nj, nk, ng, iGeometry, iEquation
    use ModBlock,              only: BlockType, ModBlock_Init

    implicit none

    character(len=*), parameter :: eos_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/FreeEOS_v51_x0.75_z0.02.dat'
    character(len=*), parameter :: entropy_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/FreeEOS_v51_x0.75_z0.02_entropy_contour_9.252752.dat'
    character(len=*), parameter :: opacity_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/opacity_x70_y28_z2.dat'

    type(BlockType), target :: block
    real(8) :: xijk_range(3,2)
    real(8) :: diffusion_total, cooling_total, net_total
    real(8) :: ratio

    ni = 2900
    nj = 1
    nk = 1
    ng = 2
    iGeometry = 1
    iEquation = 0

    call ModLookUpTable_Read(eos_file, logical_unit=1)
    call ModLookUpTable_Read(entropy_file, logical_unit=1)
    call ModLookUpTable_Read(opacity_file, logical_unit=1)

    call ModEOS_init
    call ModOpacity_init
    call ModStratification_new_Init

    xijk_range(1,:) = [0.70d0, 0.99d0]
    xijk_range(2,:) = [0.90d0, 0.96d0]
    xijk_range(3,:) = [0.00d0, 0.06d0]

    call ModBlock_Init(block, 1, xijk_range, .true., iEquation, .true.)

    diffusion_total = integrate_source(block%diffusion_I, block%V_LLL(:,1,1))
    cooling_total = integrate_source(block%cooling_I, block%V_LLL(:,1,1))
    net_total = diffusion_total + cooling_total
    ratio = abs(net_total) / max(abs(diffusion_total), abs(cooling_total))

    print *, 'ModBlock integral test'
    print *, '  radial range =', xijk_range(1,1), 'to', xijk_range(1,2)
    print *, '  grid size    =', ni, nj, nk
    print '(1x,a,es24.16)', '  diffusion total = ', diffusion_total
    print '(1x,a,es24.16)', '  cooling total   = ', cooling_total
    print '(1x,a,es24.16)', '  net total       = ', net_total
    print '(1x,a,es24.16)', '  cancellation ratio = ', ratio

    if (ratio > 1.0d-2) then
        print *, 'ERROR: heating and cooling do not cancel closely enough.'
        stop 1
    end if

    print *, 'ModBlock integral test passed.'

contains

    function integrate_source(profile, volume) result(total)
        real(8), intent(in) :: profile(:)
        real(8), intent(in) :: volume(:)
        real(8) :: total

        total = sum(profile(1:ni) * volume(1:ni)) * R_sun__CGS**3
    end function integrate_source

end program test_ModBlock_integral
