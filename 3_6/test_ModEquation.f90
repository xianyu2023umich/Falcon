program test_ModEquation

    use ModConst,              only: dpi
    use ModLookUpTable,        only: ModLookUpTable_Read
    use ModEOS,                only: ModEOS_init
    use ModOpacity,            only: ModOpacity_init
    use ModStratification_new, only: ModStratification_new_Init
    use ModParameters,         only: ni, nj, nk, ng, iGeometry, iEquation
    use ModBlock,              only: BlockType, ModBlock_Init
    use ModBoundary,           only: ModBoundary_Dynamo_HD_primitives, ModBoundary_Dynamo_MHD_primitives
    use ModSpherical,          only: ModSpherical_div
    use ModEquation,           only: ModEquation_Dynamo_HD, ModEquation_Dynamo_MHD, &
                                     ModEquation_Dynamo_Get_p1, ModEquation_Dynamo_Mass_Conservation, &
                                     ModEquation_Dynamo_Inertial_Force, ModEquation_Dynamo_Pressure_Gradient, &
                                     ModEquation_Dynamo_Gravity, ModEquation_Dynamo_Entropy_Advection, &
                                     ModEquation_Dynamo_Entropy_Heating, ModEquation_Dynamo_Lorentz_Force, &
                                     ModEquation_Dynamo_Induction_Equation

    implicit none

    character(len=*), parameter :: eos_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/FreeEOS_v51_x0.75_z0.02.dat'
    character(len=*), parameter :: entropy_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/FreeEOS_v51_x0.75_z0.02_entropy_contour_9.252752.dat'
    character(len=*), parameter :: opacity_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/opacity_x70_y28_z2.dat'

    type(BlockType), target :: block_hd_driver, block_hd_manual
    type(BlockType), target :: block_mhd_driver, block_mhd_manual
    real(8) :: xijk_range(3,2)
    real(8), allocatable :: eqn_manual(:,:,:,:)
    integer, allocatable :: seed(:)
    integer :: seed_size
    real(8) :: max_diff_hd, max_diff_mhd

    ni = 40
    nj = 2
    nk = 2
    ng = 2
    iGeometry = 1

    call ModLookUpTable_Read(eos_file, logical_unit=1)
    call ModLookUpTable_Read(entropy_file, logical_unit=1)
    call ModLookUpTable_Read(opacity_file, logical_unit=1)

    call ModEOS_init
    call ModOpacity_init
    call ModStratification_new_Init

    xijk_range(1,:) = [0.70d0, 0.99d0]
    xijk_range(2,:) = [0.75d0, 0.85d0]
    xijk_range(3,:) = [0.15d0, 0.25d0]

    allocate(eqn_manual(1:ni,1:nj,1:nk,1:8))

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    call init_seed(seed)

    iEquation = 0
    call random_seed(put=seed)
    call ModBlock_Init(block_hd_driver, 1, xijk_range, .true., iEquation, .true.)
    call random_seed(put=seed)
    call ModBlock_Init(block_hd_manual, 2, xijk_range, .true., iEquation, .true.)

    call ModEquation_Dynamo_HD(block_hd_driver, .false., eqn_manual(:,:,:,1:5))
    call manual_hd_rhs(block_hd_manual)
    eqn_manual(:,:,:,1:5) = block_hd_manual%EQN_update_R_IV(1:ni,1:nj,1:nk,1:5)
    max_diff_hd = maxval(abs(block_hd_driver%EQN_update_R_IV(1:ni,1:nj,1:nk,1:5) - eqn_manual(:,:,:,1:5)))
    call assert_close(max_diff_hd, 0.0d0, 1.0d-12, 'HD driver/manual max diff')

    iEquation = 1
    call random_seed(put=seed)
    call ModBlock_Init(block_mhd_driver, 3, xijk_range, .true., iEquation, .true.)
    call random_seed(put=seed)
    call ModBlock_Init(block_mhd_manual, 4, xijk_range, .true., iEquation, .true.)

    call ModEquation_Dynamo_MHD(block_mhd_driver, .false., eqn_manual)
    call manual_mhd_rhs(block_mhd_manual)
    eqn_manual = block_mhd_manual%EQN_update_R_IV(1:ni,1:nj,1:nk,1:8)
    max_diff_mhd = maxval(abs(block_mhd_driver%EQN_update_R_IV(1:ni,1:nj,1:nk,1:8) - eqn_manual))
    call assert_close(max_diff_mhd, 0.0d0, 1.0d-12, 'MHD driver/manual max diff')

    print *, 'ModEquation regression test passed.'
    print '(1x,a,es14.5)', '  HD max diff  = ', max_diff_hd
    print '(1x,a,es14.5)', '  MHD max diff = ', max_diff_mhd

contains

    subroutine init_seed(seed)
        integer, intent(out) :: seed(:)
        integer :: i

        do i = 1, size(seed)
            seed(i) = 12345 + 37 * i
        end do
    end subroutine init_seed

    subroutine manual_hd_rhs(block)
        type(BlockType), intent(inout), target :: block

        block%primitive => block%primitive_IV
        block%EQN_update_R_IV = 0.0d0
        call ModBoundary_Dynamo_HD_primitives(block, .false.)
        call ModEquation_Dynamo_Get_p1(block)
        call ModEquation_Dynamo_Mass_Conservation(block)
        call ModEquation_Dynamo_Inertial_Force(block)
        call ModEquation_Dynamo_Pressure_Gradient(block)
        call ModEquation_Dynamo_Gravity(block)
        call ModEquation_Dynamo_Entropy_Advection(block)
        call ModEquation_Dynamo_Entropy_Heating(block)
    end subroutine manual_hd_rhs

    subroutine manual_mhd_rhs(block)
        type(BlockType), intent(inout), target :: block
        integer :: ivar
        real(8) :: DivB(1:ni,1:nj,1:nk)

        block%primitive => block%primitive_IV
        block%EQN_update_R_IV = 0.0d0
        call ModBoundary_Dynamo_MHD_primitives(block, .false.)
        call ModEquation_Dynamo_Get_p1(block)
        call ModEquation_Dynamo_Mass_Conservation(block)
        call ModEquation_Dynamo_Inertial_Force(block)
        call ModEquation_Dynamo_Pressure_Gradient(block)
        call ModEquation_Dynamo_Lorentz_Force(block)
        call ModEquation_Dynamo_Gravity(block)
        call ModEquation_Dynamo_Entropy_Advection(block)
        call ModEquation_Dynamo_Entropy_Heating(block)
        call ModEquation_Dynamo_Induction_Equation(block)

        DivB = ModSpherical_div(ni,nj,nk,ng,block%xi_I,block%xj_I,&
            block%dxi,block%dxj,block%dxk,block%primitive(:,:,:,block%br_:block%bp_))
        do ivar = block%vr_, block%vp_
            block%EQN_update_R_IV(:,:,:,block%br_+ivar-block%vr_) = block%EQN_update_R_IV(:,:,:,block%br_+ivar-block%vr_)-&
                DivB*block%primitive(1:ni,1:nj,1:nk,ivar)
        end do
    end subroutine manual_mhd_rhs

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

end program test_ModEquation
