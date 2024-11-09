module ModSavePlot

    use ModOcTree
    use ModMath
    use MPI

    contains

    subroutine ModSave_2D(Tree,direction,position,nij_out,xij_out_range,filename,logical_unit)
        implicit none

        type(OcTree),intent(in),target  :: Tree
        character(len=*),intent(in)     :: direction
        real,intent(in)                 :: position
        integer,intent(in)              :: nij_out(2)
        real,intent(in),optional        :: xij_out_range(2,2)
        real                            :: xij_out_range_here(2,2)
        character(len=*),intent(in)     :: filename
        integer,intent(in)              :: logical_unit
        integer                         :: MpiRank,MpiSubRank,MpiSubSize

        integer                         :: iLocalBlock,i,j
        type(Block),pointer             :: Block1
        logical                         :: if_inside_single
        integer                         :: if_join,MPI_COMM_SAVE_SUBSET,ierr

        real                            :: save_primitive_local(Tree%LocalBlocks(1)%nvar,nij_out(1),nij_out(2))
        real                            :: save_primitive_global(Tree%LocalBlocks(1)%nvar,nij_out(1),nij_out(2))
        real                            :: xi_out(nij_out(1)),xj_out(nij_out(2))
        real                            :: dxi_out,dxj_out
        integer                         :: i_out_,j_out_,k_out_

        ! for write

        real :: write_xijk(3)

        call MPI_Comm_rank(MPI_COMM_WORLD, MpiRank, ierr)        

        

        ! assign the indices for directions

        select case(direction)
        case('i')
            i_out_=2; j_out_=3; k_out_=1
        case('j')
            i_out_=1; j_out_=3; k_out_=2
        case('k')
            i_out_=1; j_out_=2; k_out_=3
        case default
            i_out_=1; j_out_=2; k_out_=3
        end select

        ! the range of save map

        if (present(xij_out_range)) then
            xij_out_range_here=xij_out_range
        else
            xij_out_range_here(1,:)=Tree%xijk_range(i_out_,:)
            xij_out_range_here(2,:)=Tree%xijk_range(j_out_,:)
        end if

        do i=1,nij_out(1)
            xi_out(i)=((i-0.5)*xij_out_range_here(1,2)+(nij_out(1)-i+0.5)*xij_out_range_here(1,1))/nij_out(1)
        end do
        do j=1,nij_out(2)
            xj_out(j)=((j-0.5)*xij_out_range_here(2,2)+(nij_out(2)-j+0.5)*xij_out_range_here(2,1))/nij_out(2)
        end do

        dxi_out=xi_out(2)-xi_out(1); dxj_out=xj_out(2)-xj_out(1)

        ! loop over all the blocks to fulfill the save map.
        ! initialize if_join to be MPI_UNDEFINED

        if_join=MPI_UNDEFINED

        save_primitive_local=0.
        save_primitive_global=0.

        do iLocalBlock=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iLocalBlock)

            call ModSave_2D_single(Block1,position,i_out_,j_out_,k_out_,&
                nij_out,xi_out,xj_out,dxi_out,dxj_out,save_primitive_local,if_inside_single)
            
            if (if_inside_single) if_join=1
        end do

        ! if if_inside is true then join the new mpi_common

        call MPI_Comm_split(MPI_COMM_WORLD, if_join, MpiRank, MPI_COMM_SAVE_SUBSET, ierr)

        if (if_join==1) Then
            call MPI_Comm_rank(MPI_COMM_SAVE_SUBSET,MpiSubRank,ierr)
            call MPI_Comm_size(MPI_COMM_SAVE_SUBSET,MpiSubSize,ierr)

            call MPI_Reduce(save_primitive_local,save_primitive_global,&
                Tree%LocalBlocks(1)%nvar*nij_out(1)*nij_out(2),&
                MPI_REAL, MPI_SUM, 0, MPI_COMM_SAVE_SUBSET, ierr)
        end if
        
        ! reduce

        ! Save it by the Rank0

        if (if_join==1 .and. MpiSubRank==0) then
            
            open(unit=logical_unit,file=filename, status='replace', action='write')
            write(logical_unit,*) nij_out

            do i=1,nij_out(1); do j=1,nij_out(2)
                write_xijk(i_out_)=xi_out(i)
                write_xijk(j_out_)=xj_out(j)
                write_xijk(k_out_)=position
                write(logical_unit,*)write_xijk,save_primitive_global(:,i,j)
            end do; end do

            close(logical_unit)
            print *,'File='//filename,', saved.'
        end if

    end subroutine ModSave_2D

    subroutine ModSave_2D_single(Block1,position,i_out_,j_out_,k_out_,&
        nij_out,xi_out,xj_out,dxi_out,dxj_out,save_primitive,if_inside_single)
        implicit none

        type(Block),intent(in)          :: Block1
        real,intent(in)                 :: position
        integer,intent(in)              :: i_out_,j_out_,k_out_
        integer,intent(in)              :: nij_out(2)
        real,intent(in)                 :: xi_out(nij_out(1)),xj_out(nij_out(2))
        real,intent(in)                 :: dxi_out,dxj_out

        real,intent(inout)              :: save_primitive(:,:,:)
        logical,intent(inout)           :: if_inside_single

        real                            :: xi_out1,xj_out1
        integer                         :: i,j
        real :: coord_ijk(3),posi_ijk(3),w_ijk(3)
        integer :: posi_ijk_int(3)

        if_inside_single=.false.

        if (position .ge. Block1%xijk_range(k_out_,1) .and. position .lt. Block1%xijk_range(k_out_,2)) then

            do i=floor((Block1%xijk_range(i_out_,1)-xi_out(1))/dxi_out+1),&
                floor((Block1%xijk_range(i_out_,2)-xi_out(1))/dxi_out)
    
                do j=floor((Block1%xijk_range(j_out_,1)-xj_out(1))/dxj_out+1),&
                    floor((Block1%xijk_range(j_out_,2)-xj_out(1))/dxj_out)

                    if_inside_single=.True.

                    if (i .gt. 0 .and. j .gt. 0 .and. i .le. nij_out(1) .and. j .le. nij_out(2)) then
                        xi_out1=xi_out(i)
                        xj_out1=xj_out(j)

                        coord_ijk=[xi_out1,xj_out1,position]
                        coord_ijk([i_out_,j_out_,k_out_])=coord_ijk

                        posi_ijk(1)=-Block1%ng+1.+(coord_ijk(1)-Block1%xi(-Block1%ng+1))/Block1%dxi
                        posi_ijk(2)=-Block1%ng+1.+(coord_ijk(2)-Block1%xj(-Block1%ng+1))/Block1%dxj
                        posi_ijk(3)=-Block1%ng+1.+(coord_ijk(3)-Block1%xk(-Block1%ng+1))/Block1%dxk

                        posi_ijk_int=floor(posi_ijk)
                        w_ijk=posi_ijk-posi_ijk_int

                        save_primitive(:,i,j)=&
                            Block1%primitive(:,posi_ijk_int(1)  ,posi_ijk_int(2)  ,posi_ijk_int(3)  )*&
                            (1.-w_ijk(1))*(1.-w_ijk(2))*(1.-w_ijk(3))+&
                            Block1%primitive(:,posi_ijk_int(1)+1,posi_ijk_int(2)  ,posi_ijk_int(3)  )*&
                            (   w_ijk(1))*(1.-w_ijk(2))*(1.-w_ijk(3))+&
                            Block1%primitive(:,posi_ijk_int(1)  ,posi_ijk_int(2)+1,posi_ijk_int(3)  )*&
                            (1.-w_ijk(1))*(   w_ijk(2))*(1.-w_ijk(3))+&
                            Block1%primitive(:,posi_ijk_int(1)  ,posi_ijk_int(2)  ,posi_ijk_int(3)+1)*&
                            (1.-w_ijk(1))*(1.-w_ijk(2))*(   w_ijk(3))+&
                            Block1%primitive(:,posi_ijk_int(1)+1,posi_ijk_int(2)+1,posi_ijk_int(3)  )*&
                            (   w_ijk(1))*(   w_ijk(2))*(1.-w_ijk(3))+&
                            Block1%primitive(:,posi_ijk_int(1)  ,posi_ijk_int(2)+1,posi_ijk_int(3)+1)*&
                            (1.-w_ijk(1))*(   w_ijk(2))*(   w_ijk(3))+&
                            Block1%primitive(:,posi_ijk_int(1)+1,posi_ijk_int(2)  ,posi_ijk_int(3)+1)*&
                            (   w_ijk(1))*(1.-w_ijk(2))*(   w_ijk(3))+&
                            Block1%primitive(:,posi_ijk_int(1)+1,posi_ijk_int(2)+1,posi_ijk_int(3)+1)*&
                            (   w_ijk(1))*(   w_ijk(2))*(   w_ijk(3))
                    end if
                    
                end do
            end do
        end if

    end subroutine ModSave_2D_single


end module ModSavePlot