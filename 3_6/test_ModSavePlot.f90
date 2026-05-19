program test_ModSavePlot

    use ModBlock,              only: BlockType
    use ModLookUpTable,        only: ModLookUpTable_Read
    use ModEOS,                only: ModEOS_init
    use ModOpacity,            only: ModOpacity_init
    use ModStratification_new, only: ModStratification_new_Init, ModStratification_new_get_middle_r
    use ModParameters,         only: MpiSize, MpiRank, r_range, ni, nj, nk, ng
    use ModConst,              only: R_sun__CGS
    use ModReadParameters,     only: ModReadParameters_read
    use ModYinYangTree,        only: YYTree, YinYangTree_InitTree, YinYangTree_SetAll
    use ModCommunication,      only: ModCommunication_SetGCAll, ModCommunication_SendRecvGC
    use ModAMR,                only: ModAMR_set_grid
    use ModInitiation,         only: ModInitiation_harmonic
    use ModSavePlot,           only: ModSave_AllCells
    use MPI

    implicit none

    character(len=*), parameter :: eos_file     = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/FreeEOS_v51_x0.75_z0.02.dat'
    character(len=*), parameter :: entropy_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/FreeEOS_v51_x0.75_z0.02_entropy_contour_9.252752.dat'
    character(len=*), parameter :: opacity_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/opacity_x70_y28_z2.dat'

    type(YYTree),target :: Tree
    integer             :: ierr
    integer,parameter   :: logical_unit = 10

    call MPI_INIT(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, MpiSize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, MpiRank, ierr)

    if (MpiRank==0) print *, 'start ModSave_AllCells smoke test'
    call ModReadParameters_read('PARAM.in.allcells', 1)
    if (MpiRank==0) print *, 'read parameters'

    call ModLookUpTable_Read(eos_file,     logical_unit=1)
    call ModLookUpTable_Read(entropy_file, logical_unit=1)
    call ModLookUpTable_Read(opacity_file, logical_unit=1)
    call ModEOS_init
    call ModOpacity_init
    call ModStratification_new_Init
    if (MpiRank==0) print *, 'initialized stratification tables'
    if (MpiRank==0) print '(a,f8.5,a)', &
        '  get_middle_r([0.7,0.99]) = ', &
        ModStratification_new_get_middle_r([0.7d0*R_sun__CGS, 0.99d0*R_sun__CGS]) / R_sun__CGS, &
        ' Rs  (expect close to 0.99)'

    ! --- Tree (no AMR) ---
    call YinYangTree_InitTree(Tree, r_range)
    if (MpiRank==0) print *, 'YinYangTree_InitTree'
    call ModAMR_set_grid(Tree)
    if (MpiRank==0) print *, 'ModAMR_set_grid'
    call YinYangTree_SetAll(Tree)
    if (MpiRank==0) print *, 'YinYangTree_SetAll'
    call ModCommunication_SetGCAll(Tree)
    if (MpiRank==0) print *, 'ModCommunication_SetGCAll'

    if (MpiRank==0 .or. MpiRank==MpiSize-1) &
        print '(a,i0,a,i0,a)', '  rank ', MpiRank, ' has ', Tree%nLocalBlocks, ' blocks'

    ! --- Initiation ---
    call ModInitiation_harmonic(Tree)
    if (MpiRank==0) print *, 'ModInitiation_harmonic'

    call ModCommunication_SendRecvGC(Tree, .false.)
    if (MpiRank==0) print *, 'ModCommunication_SendRecvGC'

    ! --- Save ---
    call ModSave_AllCells(Tree, 0, logical_unit)
    if (MpiRank==0) print *, 'ModSave_AllCells done'

    call MPI_FINALIZE(ierr)

end program test_ModSavePlot
