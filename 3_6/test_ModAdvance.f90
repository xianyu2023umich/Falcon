program test_ModAdvance

    use ModBlock,           only: BlockType
    use ModLookUpTable,     only: ModLookUpTable_Read
    use ModEOS,             only: ModEOS_init
    use ModOpacity,         only: ModOpacity_init
    use ModStratification_new, only: ModStratification_new_Init
    use ModParameters,      only: MpiSize, MpiRank, r_range, CFL, ni, nj, nk, ng
    use ModReadParameters,  only: ModReadParameters_read
    use ModYinYangTree,     only: YYTree, YinYangTree_InitTree, YinYangTree_SetAll
    use ModCommunication,   only: ModCommunication_SetGCAll
    use ModAdvance,         only: ModAdvance_rk4
    use MPI

    implicit none

    character(len=*), parameter :: eos_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/FreeEOS_v51_x0.75_z0.02.dat'
    character(len=*), parameter :: entropy_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/FreeEOS_v51_x0.75_z0.02_entropy_contour_9.252752.dat'
    character(len=*), parameter :: opacity_file = '/Users/xianyu/Documents/research/fortran/Falcon/3_6/opacity_x70_y28_z2.dat'

    type(YYTree),target     ::  Tree
    integer                 ::  ierr
    real(8)                 ::  dt

    call MPI_INIT(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, MpiSize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, MpiRank, ierr)

    if (MpiRank==0) print *, 'start ModAdvance smoke test'
    call ModReadParameters_read('PARAM.in.LOCAL_AMR',1)
    if (MpiRank==0) print *, 'read parameters'
    call ModLookUpTable_Read(eos_file, logical_unit=1)
    call ModLookUpTable_Read(entropy_file, logical_unit=1)
    call ModLookUpTable_Read(opacity_file, logical_unit=1)
    call ModEOS_init
    call ModOpacity_init
    call ModStratification_new_Init
    if (MpiRank==0) print *, 'initialized stratification tables'
    call test_ModAdvance_INITIATION
    if (MpiRank==0) print *, 'initialized tree and state'
    call test_ModAdvance_PrintFirstBlockInfo(Tree)

    call ModAdvance_rk4(Tree,CFL,.true.,dt)
    if (MpiRank==0) print *, 'finished ModAdvance_rk4'

    if (MpiRank==0) then
        print *, 'ModAdvance smoke test passed.'
        print '(1x,a,es14.5)', '  dt = ', dt
    end if

    call MPI_FINALIZE(ierr)

contains

    subroutine test_ModAdvance_INITIATION
        type(BlockType),pointer :: Block1
        integer :: iLocalBlock

        call YinYangTree_InitTree(Tree, r_range)
        if (MpiRank==0) print *, '  YinYangTree_InitTree'
        call YinYangTree_SetAll(Tree, .false.)
        if (MpiRank==0) print *, '  YinYangTree_SetAll'
        call ModCommunication_SetGCAll(Tree)
        if (MpiRank==0) print *, '  ModCommunication_SetGCAll'

        do iLocalBlock=1,Tree%nLocalBlocks
            Block1=>Tree%LocalBlocks(iLocalBlock)
            Block1%primitive_IV=0.0d0
            Block1%primitive_rk_IV=0.0d0
        end do
        if (MpiRank==0) print *, '  zeroed primitives'
    end subroutine test_ModAdvance_INITIATION

    subroutine test_ModAdvance_PrintFirstBlockInfo(Tree)
        type(YYTree), intent(in), target :: Tree
        type(BlockType), pointer :: Block1

        if (Tree%nLocalBlocks <= 0) then
            if (MpiRank==0) print *, '  no local blocks on rank 0'
            return
        end if

        Block1 => Tree%LocalBlocks(1)

        if (MpiRank==0) then
            print *, '  first local block summary'
            print '(1x,a,i0)', '    rank = ', MpiRank
            print '(1x,a,i0)', '    iBlock = ', Block1%iBlock
            print '(1x,a,l1)', '    if_yin = ', Block1%if_yin
            print '(1x,a,3(i0,1x))', '    ni nj nk = ', ni, nj, nk
            print '(1x,a,i0)', '    ng = ', ng
            print '(1x,a,3(es14.5,1x))', '    dxi dxj dxk = ', Block1%dxi, Block1%dxj, Block1%dxk
            print '(1x,a,2(es14.5,1x))', '    x-range = ', Block1%xijk_range(1,1), Block1%xijk_range(1,2)
            print '(1x,a,2(es14.5,1x))', '    y-range = ', Block1%xijk_range(2,1), Block1%xijk_range(2,2)
            print '(1x,a,2(es14.5,1x))', '    z-range = ', Block1%xijk_range(3,1), Block1%xijk_range(3,2)
            print '(1x,a,2(es14.5,1x))', '    GC x-range = ', Block1%xijk_range_GC(1,1), Block1%xijk_range_GC(1,2)
            print '(1x,a,2(es14.5,1x))', '    GC y-range = ', Block1%xijk_range_GC(2,1), Block1%xijk_range_GC(2,2)
            print '(1x,a,2(es14.5,1x))', '    GC z-range = ', Block1%xijk_range_GC(3,1), Block1%xijk_range_GC(3,2)
        end if
    end subroutine test_ModAdvance_PrintFirstBlockInfo

end program test_ModAdvance
