program test1
    
    !use ModCommunication
    use ieee_arithmetic
    use ModTimeStep
    use ModSavePlot
    use ModAdvance

    implicit none

    type(OcTree),target     ::  Tree
    integer                 ::  ni,nj,nk,ng=2,i,j,k
    integer                 ::  iStep,ierr,MpiSize,MpiRank
    real                    ::  xijk_range(3,2)
    character(len=8)        ::  iStep_char

    ni=25
    nj=25
    nk=25

    ! Initiate MPI and get MpiSize/Rank
    call MPI_INIT(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, MpiSize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, MpiRank, ierr)

    ! Initiate the tree and set GC
    call test1_initial_tree_1
    print *,'I"m ',MpiRank,', I have',Tree%nLocalBlocks

    !call ModSave_2D(Tree,'k',2.16,[200,200],filename='0.dat',logical_unit=2)

    do iStep=1,7000
        call ModAdvance_rk4(Tree,0.7,0.7,.True.,MpiSize,MpiRank)
        if(MpiRank==0)print *,maxval(Tree%Localblocks(1)%primitive)

        write(iStep_char,'(I8)')iStep
        if (mod(iStep,50)==0) then
                
            do i=1,len_trim(iStep_char)
                if (iStep_char(i:i)==' ') iStep_char(i:i)='0'
            end do
            call ModSave_2D(Tree,'i',0.1,[100,100],filename=iStep_char//'.dat',logical_unit=2)
        end if

        !do i=-ng+1,ni+ng; do j=-ng+1,nj+ng; do k=-ng+1,nk+ng
        !    if (ieee_is_nan(Tree%localblocks(1)%primitive_k2(1,i,j,k))) then
        !        print *,i,j,k
        !        print *,Tree%localblocks(1)%primitive_k2(:,i,j,k)
        !        stop 1
        !    end if            
        !end do; end do; end do

        !if (MpiRank==0)print *,Tree%LocalBlocks(1)%primitive(1,:,12,12)


        if (MpiRank==0) print *,'iStep=',iStep,'Completed...'
        !if (MpiRank==0) print *,maxval(Tree%LocalBlocks(1)%primitive(4,:,:,:))
    end do

    
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call mpi_finalize(ierr)

    contains

    subroutine test1_initial_tree_1
        integer ::  iChild1
        xijk_range(1,:)=[-1.09,1.09]
        xijk_range(2,:)=[-1.09,1.09]
        xijk_range(3,:)=[0.,2.16]

        call OcTree_InitTree(Tree,xijk_range,ijk_if_periodic=[.True.,.True.,.False.])
        call OcTree_Divide(Tree,Tree%RootNode,[.True.,.True.,.True.])
        do iChild1=1,8
            call OcTree_Divide(Tree,Tree%RootNode%children(iChild1),[.True.,.True.,.True.])
        end do
        call OcTree_SetAll(Tree,ni,nj,nk,ng,'Dynamo_HD',4,MpiSize,MpiRank)
        call ModGC_InitGCAll(Tree,MpiSize)
        call ModGC_GetGC_SourcesAll(Tree,MpiSize)
    end subroutine test1_initial_tree_1

end program test1