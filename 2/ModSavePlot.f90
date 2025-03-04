module ModSavePlot

    use ModYinYangTree
    use ModMath
    use ModConst,only: dpi
    use MPI

    contains

    subroutine ModSave_Globe(Tree,r_save,ntp_out,filename,logical_unit)
        implicit none
        type(YYTree),intent(in),target ::  Tree
        real,intent(in)                     ::  r_save
        integer,intent(in)                  ::  ntp_out(2)
        character(len=*),intent(in)         ::  filename
        integer,intent(in)                  ::  logical_unit
        integer                             ::  MpiRank,MpiSubRank,MpiSubSize

        integer                             ::  iLocalBlock
        type(Block),pointer                 ::  Block1
        logical                             ::  if_inside_single
        integer                             ::  if_join,MPI_COMM_SAVE_SUBSET,ierr

        real                                ::  save_primitive_local(&
                                                Tree%LocalBlocks(1)%nvar,&
                                                ntp_out(1),ntp_out(2))
        real                                ::  save_primitive_global(&
                                                Tree%LocalBlocks(1)%nvar,&
                                                ntp_out(1),ntp_out(2))
        real                                ::  t_out(ntp_out(1)),p_out(ntp_out(2))
        integer                             ::  it,ip

        real                                ::  write_tp(2)
        
        ! MPI set up
        call MPi_Comm_rank(MPI_COMM_WORLD,MpiRank, ierr)

        ! Get the th and ph lists
        do it=1,ntp_out(1)
            t_out(it)=(it-0.50)/ntp_out(1)*dpi
        end do
        do ip=1,ntp_out(2)
            p_out(ip)=(ip-0.50)/ntp_out(2)*dpi*2-dpi
        end do

        if_join=MPI_UNDEFINED

        save_primitive_local=-1.e30
        save_primitive_global=0.

        do iLocalBlock=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iLocalBlock)
            call ModSave_Globe_single(Block1,r_save,ntp_out,t_out,p_out,&
                save_primitive_local,if_inside_single)
            
            if (if_inside_single) if_join=1
        end do

        ! Get the new communicator

        call MPI_Comm_split(MPI_COMM_WORLD, if_join, MpiRank, MPI_COMM_SAVE_SUBSET, ierr)

        ! Reduce

        if (if_join==1) Then
            call MPI_Comm_rank(MPI_COMM_SAVE_SUBSET,MpiSubRank,ierr)
            call MPI_Comm_size(MPI_COMM_SAVE_SUBSET,MpiSubSize,ierr)

            call MPI_Reduce(save_primitive_local,save_primitive_global,&
                Tree%LocalBlocks(1)%nvar*ntp_out(1)*ntp_out(2),&
                MPI_REAL, MPI_MAX, 0, MPI_COMM_SAVE_SUBSET, ierr)
        end if

        ! Save it by the SubRank=0

        if (if_join==1 .and. MpiSubRank==0) then
            
            open(unit=logical_unit,file=filename, status='replace', action='write')
            write(logical_unit,*) ntp_out

            do it=1,ntp_out(1); do ip=1,ntp_out(2)
                write_tp=[t_out(it),p_out(ip)]
                write(logical_unit,*)write_tp,save_primitive_global(:,it,ip)
            end do; end do

            close(logical_unit)
            print *,'File='//filename,', saved.'
        end if
    end subroutine ModSave_Globe

    subroutine ModSave_Globe_single(Block1,r_save,ntp_out,t_out,p_out,save_primitive,if_inside_single)
        implicit none

        type(Block),intent(in)              ::  Block1
        real,intent(in)                     ::  r_save
        integer,intent(in)                  ::  ntp_out(2)
        real,intent(in)                     ::  t_out(ntp_out(1)),p_out(ntp_out(2))
        real,intent(inout)                  ::  save_primitive(:,:,:)
        logical,intent(inout)               ::  if_inside_single

        integer                             ::  it,ip,rtp_posi_int(3)
        real                                ::  rtp(3),rtp_posi(3),rtp_weight(3)

        ! See if the Block1 has some points of radius r_save
        ! If yes, then do interpolation
        if ((Block1%xijk_range(1,1)-r_save)*(Block1%xijk_range(1,2)-r_save)<0.0) then
            if_inside_single=.true.

            ! Loop each point in the map
            ! if Yang, then first convert point coord to Yang and then
            ! see if it's inside the block.
            do it=1,ntp_out(1); do ip=1,ntp_out(2)
                rtp=[r_save,t_out(it),p_out(ip)]
                if (.not. Block1%if_yin) rtp=ModYinYang_CoordConv_0D(rtp)

                if ((Block1%xj(-Block1%ng+1)-rtp(2))*(Block1%xj(Block1%ng+Block1%nj)-rtp(2)).le.0.0 .and. &
                    (Block1%xk(-Block1%ng+1)-rtp(3))*(Block1%xk(Block1%ng+Block1%nk)-rtp(3)).le.0.0) then
                    
                    rtp_posi(1)=-Block1%ng+1.+(rtp(1)-Block1%xi(-Block1%ng+1))/Block1%dxi
                    rtp_posi(2)=-Block1%ng+1.+(rtp(2)-Block1%xj(-Block1%ng+1))/Block1%dxj
                    rtp_posi(3)=-Block1%ng+1.+(rtp(3)-Block1%xk(-Block1%ng+1))/Block1%dxk

                    rtp_posi_int=floor(rtp_posi)
                    rtp_weight=rtp_posi-rtp_posi_int

                    save_primitive(:,it,ip)=&
                        Block1%primitive(rtp_posi_int(1)  ,rtp_posi_int(2)  ,rtp_posi_int(3)  ,:)*&
                        (1.-rtp_weight(1))*(1.-rtp_weight(2))*(1.-rtp_weight(3))+&
                        Block1%primitive(rtp_posi_int(1)+1,rtp_posi_int(2)  ,rtp_posi_int(3)  ,:)*&
                        (   rtp_weight(1))*(1.-rtp_weight(2))*(1.-rtp_weight(3))+&
                        Block1%primitive(rtp_posi_int(1)  ,rtp_posi_int(2)+1,rtp_posi_int(3)  ,:)*&
                        (1.-rtp_weight(1))*(   rtp_weight(2))*(1.-rtp_weight(3))+&
                        Block1%primitive(rtp_posi_int(1)  ,rtp_posi_int(2)  ,rtp_posi_int(3)+1,:)*&
                        (1.-rtp_weight(1))*(1.-rtp_weight(2))*(   rtp_weight(3))+&
                        Block1%primitive(rtp_posi_int(1)+1,rtp_posi_int(2)+1,rtp_posi_int(3)  ,:)*&
                        (   rtp_weight(1))*(   rtp_weight(2))*(1.-rtp_weight(3))+&
                        Block1%primitive(rtp_posi_int(1)  ,rtp_posi_int(2)+1,rtp_posi_int(3)+1,:)*&
                        (1.-rtp_weight(1))*(   rtp_weight(2))*(   rtp_weight(3))+&
                        Block1%primitive(rtp_posi_int(1)+1,rtp_posi_int(2)  ,rtp_posi_int(3)+1,:)*&
                        (   rtp_weight(1))*(1.-rtp_weight(2))*(   rtp_weight(3))+&
                        Block1%primitive(rtp_posi_int(1)+1,rtp_posi_int(2)+1,rtp_posi_int(3)+1,:)*&
                        (   rtp_weight(1))*(   rtp_weight(2))*(   rtp_weight(3))
                    
                    if (.not. Block1%if_yin) save_primitive(2:4,it,ip)=&
                        ModYinYang_VecConv_0D(rtp,save_primitive(2:4,it,ip))
                    
                    !if (abs(save_primitive(3,it,ip)-1.0)>0.05) then
                    !    print *,'================'
                    !    print *,save_primitive(2:4,it,ip)
                    !    print *,rtp_posi_int
                    !    print *,Block1%primitive(2:4,rtp_posi_int(1),rtp_posi_int(2),rtp_posi_int(3))
                    !    print *,Block1%iBlock,it,ip
                    !end if
                end if
            end do; end do
        else
            if_inside_single=.false.
        end if
        
    end subroutine ModSave_Globe_single

end module ModSavePlot