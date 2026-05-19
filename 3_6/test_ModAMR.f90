program test_ModAMR

    use ModBlock,              only: BlockType
    use ModLookUpTable,        only: ModLookUpTable_Read
    use ModEOS,                only: ModEOS_init
    use ModOpacity,            only: ModOpacity_init
    use ModStratification_new, only: ModStratification_new_Init
    use ModParameters,         only: MpiSize, MpiRank, r_range, CFL, ni, nj, nk, ng
    use ModReadParameters,     only: ModReadParameters_read
    use ModYinYangTree,        only: YYTree, YinYangTree_InitTree, YinYangTree_SetAll
    use ModCommunication,      only: ModCommunication_SetGCAll
    use ModAMR,                only: ModAMR_set_grid
    use MPI

    implicit none

    character(len=*), parameter :: eos_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/FreeEOS_v51_x0.75_z0.02.dat'
    character(len=*), parameter :: entropy_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/FreeEOS_v51_x0.75_z0.02_entropy_contour_9.252752.dat'
    character(len=*), parameter :: opacity_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/opacity_x70_y28_z2.dat'

    type(YYTree),target     ::  Tree
    integer                 ::  ierr

    call MPI_INIT(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, MpiSize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, MpiRank, ierr)

    if (MpiRank==0) print *, 'start ModAMR smoke test'
    call ModReadParameters_read('PARAM.in.LOCAL_AMR',1)
    if (MpiRank==0) print *, 'read parameters'
    call ModLookUpTable_Read(eos_file, logical_unit=1)
    call ModLookUpTable_Read(entropy_file, logical_unit=1)
    call ModLookUpTable_Read(opacity_file, logical_unit=1)
    call ModEOS_init
    call ModOpacity_init
    call ModStratification_new_Init
    if (MpiRank==0) print *, 'initialized stratification tables'
    call test_ModAMR_INITIATION(Tree)
    if (MpiRank==0) print *, 'initialized AMR tree and state'

    call MPI_FINALIZE(ierr)

contains

    subroutine test_ModAMR_INITIATION(Tree)
        type(YYTree),intent(inout),target :: Tree
        type(BlockType),pointer :: Block1
        integer :: iLocalBlock

        call YinYangTree_InitTree(Tree, r_range)
        if (MpiRank==0) print *, '  YinYangTree_InitTree'

        call ModAMR_set_grid(Tree)
        if (MpiRank==0) print *, '  ModAMR_set_grid'

        call YinYangTree_SetAll(Tree, .true.)
        if (MpiRank==0) print *, '  YinYangTree_SetAll'

        call ModCommunication_SetGCAll(Tree)
        if (MpiRank==0) print *, '  ModCommunication_SetGCAll'

        do iLocalBlock=1,Tree%nLocalBlocks
            Block1=>Tree%LocalBlocks(iLocalBlock)
            Block1%primitive_IV=0.0d0
            Block1%primitive_rk_IV=0.0d0
        end do
        if (MpiRank==0) print *, '  zeroed primitives'

        if (MpiRank==0 .or. MpiRank==MpiSize-1) then
            print *,'Tree%NumLeafNodes=',Tree%NumLeafNodes,'MpiRank=',MpiRank,' has ',Tree%nLocalBlocks,' blocks'
        end if
    end subroutine test_ModAMR_INITIATION

end program test_ModAMR
