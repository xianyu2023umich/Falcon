module ModSavePlot

    use ieee_arithmetic, only:  ieee_is_nan
    use ModBlock,       only:   BlockType
    use ModYinYang,     only:   ModYinYang_CoordConv_0D,&
                                ModYinYang_VecConv_0D
    
    use ModYinYangTree, only:   YYTree
    use ModConst,       only:   dpi
    use ModParameters,  only:   nvar,ni,nj,nk,ng,r_range,PlotType,Plots,nPlots
    use ModVariables,   only:   vr_,vt_,vp_,br_,bt_,bp_
    use ModMath,        only:   ModMath_IfLinesInterSect
    use MPI

    contains

    subroutine ModSave_DoAll(Tree,iStep)
        implicit none
        type(YYTree),intent(in),target      ::  Tree
        integer,intent(in)                  ::  iStep
        character(len=8)                    ::  iStep_char

        integer                             ::  iPlot
        type(PlotType),pointer              ::  Plot1
        integer                             ::  i

        ! Loop each plot
        do iPlot=1,nPlots
            Plot1=>Plots(iPlot)

            ! First see if it's time to save the plot
            if (mod(iStep,Plot1%nStepsSavePlot)==0) then

                ! Get the step char
                write(iStep_char,'(I8)')iStep
                do i=1,len_trim(iStep_char)
                    if (iStep_char(i:i)==' ') iStep_char(i:i)='0'
                end do

                ! See which type it is
                select case(Plot1%iType)
                case(0)
                    call ModSave_Globe_MPI_Accumulate(Tree,&
                        Plot1%rtp_SavePlot(1),[Plot1%nrtp_SavePlot(2),Plot1%nrtp_SavePlot(3)],&
                        'Sphere_'//iStep_char//'.dat',Plot1%logical_unit)
                case(1)
                    call ModSave_Fan(Tree,&
                        Plot1%rtp_SavePlot(3),[Plot1%nrtp_SavePlot(1),Plot1%nrtp_SavePlot(2)],&
                        'Fan_'//iStep_char//'.dat',Plot1%logical_unit)
                case(2)
                    call ModSave_Cube(Tree,&
                        Plot1%rtp_range_SavePlot(1,:),Plot1%nrtp_SavePlot,&
                        'Cube_'//iStep_char//'.dat',Plot1%logical_unit)
                end select
            end if
        end do
    end subroutine ModSave_DoAll

    subroutine ModSave_Cube(Tree,r_save_range,nrtp_out,filename,logical_unit)
        implicit none
        type(YYTree),intent(in),target      ::  Tree
        real,intent(in)                     ::  r_save_range(2)
        integer,intent(in)                  ::  nrtp_out(3)
        character(len=*),intent(in)         ::  filename
        integer,intent(in)                  ::  logical_unit
        type(BlockType),pointer             ::  Block1
        integer                             ::  iLocalBlock

        integer                             ::  MpiRank
        integer(kind=MPI_ADDRESS_KIND)      ::  bytes
        integer                             ::  disp_unit
        integer                             ::  MPI_window,ierr

        integer                             ::  ir,it,ip,rtp_posi_int(3)
        real                                ::  save_primitive_IV(nvar,nrtp_out(1),nrtp_out(2),nrtp_out(3))
        real                                ::  primitives_I(nvar)
        real                                ::  write_xyz(3)
        real                                ::  r_out(nrtp_out(1)),t_out(nrtp_out(2)),p_out(nrtp_out(3))
        real                                ::  rtp(3),rtp_posi(3),rtp_weight(3)

        ! MPI set up
        call MPi_Comm_rank(MPI_COMM_WORLD,MpiRank, ierr)

        ! Get the r, th and ph lists

        do ir=1,nrtp_out(1)
            r_out(ir)=((ir-1.0)*r_save_range(2)+&
                (nrtp_out(1)-ir)*r_save_range(1))/(nrtp_out(1)-1.0)
        end do
        do it=1,nrtp_out(2)
            t_out(it)=(it-0.50)/nrtp_out(2)*dpi
        end do
        do ip=1,nrtp_out(3)
            p_out(ip)=(ip-0.50)/nrtp_out(3)*dpi*2-dpi
        end do

        ! Get the bytes size of buffer
        disp_unit=sizeof(save_primitive_IV(1,1,1,1))
        bytes=nvar*nrtp_out(1)*nrtp_out(2)*nrtp_out(3)*disp_unit
        
        ! Create the window
        call MPI_Win_create(save_primitive_IV, bytes, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, MPI_window, ierr)
        call MPI_Win_lock_all(0, MPI_window, ierr)
        call MPI_Barrier(MPI_COMM_WORLD,ierr)

        ! Loop each block
        do iLocalBlock=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iLocalBlock)

            if (ModMath_IfLinesInterSect(Block1%xijk_range(1,:),r_save_range)) then
    
                ! Loop each point in the map
                ! if Yang, then first convert point coord to Yang and then
                ! see if it's inside the block.
                do ir=1,nrtp_out(1)
                    do it=1,nrtp_out(2)
                        do ip=1,nrtp_out(3)
                            rtp=[r_out(ir),t_out(it),p_out(ip)]
                            if (.not. Block1%if_yin) rtp=ModYinYang_CoordConv_0D(rtp)
        
                            if ((Block1%xi_I(-ng+1)-rtp(1))*(Block1%xi_I(ng+ni)-rtp(1)).le.0.0 .and. &
                                (Block1%xj_I(-ng+1)-rtp(2))*(Block1%xj_I(ng+nj)-rtp(2)).le.0.0 .and. &
                                (Block1%xk_I(-ng+1)-rtp(3))*(Block1%xk_I(ng+nk)-rtp(3)).le.0.0) then
                                
                                rtp_posi(1)=-ng+1.+(rtp(1)-Block1%xi_I(-ng+1))/Block1%dxi
                                rtp_posi(2)=-ng+1.+(rtp(2)-Block1%xj_I(-ng+1))/Block1%dxj
                                rtp_posi(3)=-ng+1.+(rtp(3)-Block1%xk_I(-ng+1))/Block1%dxk
        
                                rtp_posi_int=floor(rtp_posi)
                                rtp_weight=rtp_posi-rtp_posi_int
        
                                primitives_I(:)=&
                                    Block1%primitive_IV(rtp_posi_int(1)  ,rtp_posi_int(2)  ,rtp_posi_int(3)  ,:)*&
                                    (1.-rtp_weight(1))*(1.-rtp_weight(2))*(1.-rtp_weight(3))+&
                                    Block1%primitive_IV(rtp_posi_int(1)+1,rtp_posi_int(2)  ,rtp_posi_int(3)  ,:)*&
                                    (   rtp_weight(1))*(1.-rtp_weight(2))*(1.-rtp_weight(3))+&
                                    Block1%primitive_IV(rtp_posi_int(1)  ,rtp_posi_int(2)+1,rtp_posi_int(3)  ,:)*&
                                    (1.-rtp_weight(1))*(   rtp_weight(2))*(1.-rtp_weight(3))+&
                                    Block1%primitive_IV(rtp_posi_int(1)  ,rtp_posi_int(2)  ,rtp_posi_int(3)+1,:)*&
                                    (1.-rtp_weight(1))*(1.-rtp_weight(2))*(   rtp_weight(3))+&
                                    Block1%primitive_IV(rtp_posi_int(1)+1,rtp_posi_int(2)+1,rtp_posi_int(3)  ,:)*&
                                    (   rtp_weight(1))*(   rtp_weight(2))*(1.-rtp_weight(3))+&
                                    Block1%primitive_IV(rtp_posi_int(1)  ,rtp_posi_int(2)+1,rtp_posi_int(3)+1,:)*&
                                    (1.-rtp_weight(1))*(   rtp_weight(2))*(   rtp_weight(3))+&
                                    Block1%primitive_IV(rtp_posi_int(1)+1,rtp_posi_int(2)  ,rtp_posi_int(3)+1,:)*&
                                    (   rtp_weight(1))*(1.-rtp_weight(2))*(   rtp_weight(3))+&
                                    Block1%primitive_IV(rtp_posi_int(1)+1,rtp_posi_int(2)+1,rtp_posi_int(3)+1,:)*&
                                    (   rtp_weight(1))*(   rtp_weight(2))*(   rtp_weight(3))
                                
                                if (.not. Block1%if_yin) primitives_I(vr_:vp_)=&
                                    ModYinYang_VecConv_0D(rtp,primitives_I(vr_:vp_))

                                call MPI_Accumulate(primitives_I, nvar, MPI_REAL, 0,&
                                    int(nvar*(ir-1+nrtp_out(1)*(it-1+nrtp_out(2)*(ip-1))),MPI_ADDRESS_KIND),&
                                    nvar, MPI_REAL, MPI_REPLACE, MPI_window, ierr)
                            end if
                        end do
                    end do
                end do
            end if
        end do

        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_Win_flush(0, MPI_window, ierr)
        call MPI_Win_unlock_all(MPI_window, ierr)

        if (MpiRank==0) then
            
            open(unit=logical_unit,file=filename, status='replace', action='write')
            write(logical_unit,*) 'TITLE = "My Data on Sphere"'
            write(logical_unit,*) 'VARIABLES = "K", "J", "I", "X", "Y", "Z", "Rho1", "Vr", "Vt", "Vp", "S1"'
            write(logical_unit,*) 'ZONE T="STRUCTURE GRID", I=', nrtp_out(1),&
             ', J=', nrtp_out(2), ', K=', nrtp_out(3), ', F=POINT'

            do ip=1,nrtp_out(3)
                do it=1,nrtp_out(2)
                    do ir=1,nrtp_out(1)
                        write_xyz(1)=r_out(ir)*sin(t_out(it))*cos(p_out(ip))
                        write_xyz(2)=r_out(ir)*sin(t_out(it))*sin(p_out(ip))
                        write_xyz(3)=r_out(ir)*cos(t_out(it))

                        write(logical_unit,'(3I5,8ES16.8)') ir,it,ip,write_xyz,save_primitive_IV(:,ir,it,ip)
                    end do
                end do
            end do
            close(logical_unit)
            print *,'File='//filename,', saved.'
        end if
    end subroutine ModSave_Cube

    subroutine ModSave_Globe_MPI_Accumulate(Tree,r_save,ntp_out,filename,logical_unit)
        implicit none
        type(YYTree),intent(in),target      ::  Tree
        real,intent(in)                     ::  r_save
        integer,intent(in)                  ::  ntp_out(2)
        character(len=*),intent(in)         ::  filename
        integer,intent(in)                  ::  logical_unit
        type(BlockType),pointer             ::  Block1
        integer                             ::  iLocalBlock

        integer                             ::  MpiRank
        integer(kind=MPI_ADDRESS_KIND)      ::  bytes
        integer                             ::  disp_unit
        integer                             ::  MPI_window,ierr

        integer                             ::  it,ip,rtp_posi_int(3)
        real                                ::  save_primitive_III(nvar,ntp_out(1),ntp_out(2))
        real                                ::  primitives_I(nvar)
        real                                ::  write_tp(2),write_xyz(3)
        real                                ::  t_out(ntp_out(1)),p_out(ntp_out(2))
        real                                ::  rtp(3),rtp_posi(3),rtp_weight(3)

        ! MPI set up
        call MPi_Comm_rank(MPI_COMM_WORLD,MpiRank, ierr)

        ! Get the th and ph lists
        do it=1,ntp_out(1)
            t_out(it)=(it-0.50)/ntp_out(1)*dpi
        end do
        do ip=1,ntp_out(2)
            p_out(ip)=(ip-0.50)/ntp_out(2)*dpi*2-dpi
        end do

        ! Get the bytes size of buffer
        disp_unit=sizeof(save_primitive_III(1,1,1))
        bytes=nvar*ntp_out(1)*ntp_out(2)*disp_unit
        
        ! Create the window
        call MPI_Win_create(save_primitive_III, bytes, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, MPI_window, ierr)
        call MPI_Win_lock_all(0, MPI_window, ierr)
        call MPI_Barrier(MPI_COMM_WORLD,ierr)

        ! Loop each block
        do iLocalBlock=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iLocalBlock)

            if ((Block1%xijk_range(1,1)-r_save)*(Block1%xijk_range(1,2)-r_save)<0.0) then
    
                ! Loop each point in the map
                ! if Yang, then first convert point coord to Yang and then
                ! see if it's inside the block.
                do it=1,ntp_out(1)
                    do ip=1,ntp_out(2)
                        rtp=[r_save,t_out(it),p_out(ip)]
                        if (.not. Block1%if_yin) rtp=ModYinYang_CoordConv_0D(rtp)
    
                        if ((Block1%xj_I(-ng+1)-rtp(2))*(Block1%xj_I(ng+nj)-rtp(2)).le.0.0 .and. &
                            (Block1%xk_I(-ng+1)-rtp(3))*(Block1%xk_I(ng+nk)-rtp(3)).le.0.0) then
                            
                            rtp_posi(1)=-ng+1.+(rtp(1)-Block1%xi_I(-ng+1))/Block1%dxi
                            rtp_posi(2)=-ng+1.+(rtp(2)-Block1%xj_I(-ng+1))/Block1%dxj
                            rtp_posi(3)=-ng+1.+(rtp(3)-Block1%xk_I(-ng+1))/Block1%dxk
    
                            rtp_posi_int=floor(rtp_posi)
                            rtp_weight=rtp_posi-rtp_posi_int
    
                            primitives_I(:)=&
                                Block1%primitive_IV(rtp_posi_int(1)  ,rtp_posi_int(2)  ,rtp_posi_int(3)  ,:)*&
                                (1.-rtp_weight(1))*(1.-rtp_weight(2))*(1.-rtp_weight(3))+&
                                Block1%primitive_IV(rtp_posi_int(1)+1,rtp_posi_int(2)  ,rtp_posi_int(3)  ,:)*&
                                (   rtp_weight(1))*(1.-rtp_weight(2))*(1.-rtp_weight(3))+&
                                Block1%primitive_IV(rtp_posi_int(1)  ,rtp_posi_int(2)+1,rtp_posi_int(3)  ,:)*&
                                (1.-rtp_weight(1))*(   rtp_weight(2))*(1.-rtp_weight(3))+&
                                Block1%primitive_IV(rtp_posi_int(1)  ,rtp_posi_int(2)  ,rtp_posi_int(3)+1,:)*&
                                (1.-rtp_weight(1))*(1.-rtp_weight(2))*(   rtp_weight(3))+&
                                Block1%primitive_IV(rtp_posi_int(1)+1,rtp_posi_int(2)+1,rtp_posi_int(3)  ,:)*&
                                (   rtp_weight(1))*(   rtp_weight(2))*(1.-rtp_weight(3))+&
                                Block1%primitive_IV(rtp_posi_int(1)  ,rtp_posi_int(2)+1,rtp_posi_int(3)+1,:)*&
                                (1.-rtp_weight(1))*(   rtp_weight(2))*(   rtp_weight(3))+&
                                Block1%primitive_IV(rtp_posi_int(1)+1,rtp_posi_int(2)  ,rtp_posi_int(3)+1,:)*&
                                (   rtp_weight(1))*(1.-rtp_weight(2))*(   rtp_weight(3))+&
                                Block1%primitive_IV(rtp_posi_int(1)+1,rtp_posi_int(2)+1,rtp_posi_int(3)+1,:)*&
                                (   rtp_weight(1))*(   rtp_weight(2))*(   rtp_weight(3))
                            
                            if (.not. Block1%if_yin) primitives_I(vr_:vp_)=&
                                ModYinYang_VecConv_0D(rtp,primitives_I(vr_:vp_))

                            call MPI_Accumulate(primitives_I, nvar, MPI_REAL, 0,&
                                int(nvar*(it-1+ntp_out(1)*(ip-1)),MPI_ADDRESS_KIND),&
                                nvar, MPI_REAL, MPI_REPLACE, MPI_window, ierr)
                        end if
                    end do
                end do
            end if
        end do

        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_Win_flush(0, MPI_window, ierr)
        call MPI_Win_unlock_all(MPI_window, ierr)

        if (MpiRank==0) then
            
            open(unit=logical_unit,file=filename, status='replace', action='write')
            write(logical_unit,*) 'TITLE = "My Data on Sphere"'
            write(logical_unit,*) 'VARIABLES = "K", "J", "I", "X", "Y", "Z", "Rho1", "Vr", "Vt", "Vp", "S1"'
            write(logical_unit,*) 'ZONE T="STRUCTURE GRID", I=2, J=', ntp_out(1), ', K=', ntp_out(2), ', F=POINT'

            do ip=1,ntp_out(2)
                do it=1,ntp_out(1)
                    write_tp=[t_out(it),p_out(ip)]
                    write_xyz(1)=r_save*sin(t_out(it))*cos(p_out(ip))
                    write_xyz(2)=r_save*sin(t_out(it))*sin(p_out(ip))
                    write_xyz(3)=r_save*cos(t_out(it))

                    write(logical_unit,'(3I5,8ES16.8)') ip,it,1,write_xyz,save_primitive_III(:,it,ip)
                    write(logical_unit,'(3I5,8ES16.8)') ip,it,2,write_xyz,save_primitive_III(:,it,ip)
                end do
            end do
            close(logical_unit)
            print *,'File='//filename,', saved.'
        end if

    end subroutine ModSave_Globe_MPI_Accumulate

    subroutine ModSave_Globe(Tree,r_save,ntp_out,filename,logical_unit)
        implicit none
        type(YYTree),intent(in),target ::  Tree
        real,intent(in)                     ::  r_save
        integer,intent(in)                  ::  ntp_out(2)
        character(len=*),intent(in)         ::  filename
        integer,intent(in)                  ::  logical_unit
        integer                             ::  MpiRank,MpiSubRank,MpiSubSize

        integer                             ::  iLocalBlock
        type(BlockType),pointer             ::  Block1
        logical                             ::  if_inside_single
        integer                             ::  if_join,MPI_COMM_SAVE_SUBSET,ierr

        real                                ::  save_primitive_local(&
                                                nvar,ntp_out(1),ntp_out(2))
        real                                ::  save_primitive_global(&
                                                nvar,ntp_out(1),ntp_out(2))
        real                                ::  t_out(ntp_out(1)),p_out(ntp_out(2))
        integer                             ::  it,ip

        real                                ::  write_tp(2),write_xyz(3)
        
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
                nvar*ntp_out(1)*ntp_out(2),&
                MPI_REAL, MPI_MAX, 0, MPI_COMM_SAVE_SUBSET, ierr)
        end if

        ! Save it by the SubRank=0

        if (if_join==1 .and. MpiSubRank==0) then
            
            open(unit=logical_unit,file=filename, status='replace', action='write')
            write(10,*) 'TITLE = "My Data on Sphere"'
            write(10,*) 'VARIABLES = "K", "J", "I", "X", "Y", "Z", "Rho1", "Vr", "Vt", "Vp", "S1"'
            write(10,*) 'ZONE T="STRUCTURE GRID", I=2, J=', ntp_out(1), ', K=', ntp_out(2), ', F=POINT'

            do ip=1,ntp_out(2)
                do it=1,ntp_out(1)
                    write_tp=[t_out(it),p_out(ip)]
                    write_xyz(1)=r_save*sin(t_out(it))*cos(p_out(ip))
                    write_xyz(2)=r_save*sin(t_out(it))*sin(p_out(ip))
                    write_xyz(3)=r_save*cos(t_out(it))

                    write(logical_unit,'(3I5,8ES16.8)') ip,it,1,write_xyz,save_primitive_global(:,it,ip)
                    write(logical_unit,'(3I5,8ES16.8)') ip,it,2,write_xyz,save_primitive_global(:,it,ip)
                end do
            end do
            close(logical_unit)
            print *,'File='//filename,', saved.'
        end if
    end subroutine ModSave_Globe

    subroutine ModSave_Globe_single(Block1,r_save,ntp_out,t_out,p_out,save_primitive,if_inside_single)
        implicit none

        type(BlockType),intent(in)          ::  Block1
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
            do it=1,ntp_out(1)
                do ip=1,ntp_out(2)
                    rtp=[r_save,t_out(it),p_out(ip)]
                    if (.not. Block1%if_yin) rtp=ModYinYang_CoordConv_0D(rtp)

                    if ((Block1%xj_I(-ng+1)-rtp(2))*(Block1%xj_I(ng+nj)-rtp(2)).le.0.0 .and. &
                        (Block1%xk_I(-ng+1)-rtp(3))*(Block1%xk_I(ng+nk)-rtp(3)).le.0.0) then
                        
                        rtp_posi(1)=-ng+1.+(rtp(1)-Block1%xi_I(-ng+1))/Block1%dxi
                        rtp_posi(2)=-ng+1.+(rtp(2)-Block1%xj_I(-ng+1))/Block1%dxj
                        rtp_posi(3)=-ng+1.+(rtp(3)-Block1%xk_I(-ng+1))/Block1%dxk

                        rtp_posi_int=floor(rtp_posi)
                        rtp_weight=rtp_posi-rtp_posi_int

                        save_primitive(:,it,ip)=&
                            Block1%primitive_IV(rtp_posi_int(1)  ,rtp_posi_int(2)  ,rtp_posi_int(3)  ,:)*&
                            (1.-rtp_weight(1))*(1.-rtp_weight(2))*(1.-rtp_weight(3))+&
                            Block1%primitive_IV(rtp_posi_int(1)+1,rtp_posi_int(2)  ,rtp_posi_int(3)  ,:)*&
                            (   rtp_weight(1))*(1.-rtp_weight(2))*(1.-rtp_weight(3))+&
                            Block1%primitive_IV(rtp_posi_int(1)  ,rtp_posi_int(2)+1,rtp_posi_int(3)  ,:)*&
                            (1.-rtp_weight(1))*(   rtp_weight(2))*(1.-rtp_weight(3))+&
                            Block1%primitive_IV(rtp_posi_int(1)  ,rtp_posi_int(2)  ,rtp_posi_int(3)+1,:)*&
                            (1.-rtp_weight(1))*(1.-rtp_weight(2))*(   rtp_weight(3))+&
                            Block1%primitive_IV(rtp_posi_int(1)+1,rtp_posi_int(2)+1,rtp_posi_int(3)  ,:)*&
                            (   rtp_weight(1))*(   rtp_weight(2))*(1.-rtp_weight(3))+&
                            Block1%primitive_IV(rtp_posi_int(1)  ,rtp_posi_int(2)+1,rtp_posi_int(3)+1,:)*&
                            (1.-rtp_weight(1))*(   rtp_weight(2))*(   rtp_weight(3))+&
                            Block1%primitive_IV(rtp_posi_int(1)+1,rtp_posi_int(2)  ,rtp_posi_int(3)+1,:)*&
                            (   rtp_weight(1))*(1.-rtp_weight(2))*(   rtp_weight(3))+&
                            Block1%primitive_IV(rtp_posi_int(1)+1,rtp_posi_int(2)+1,rtp_posi_int(3)+1,:)*&
                            (   rtp_weight(1))*(   rtp_weight(2))*(   rtp_weight(3))
                        
                        if (.not. Block1%if_yin) save_primitive(vr_:vp_,it,ip)=&
                            ModYinYang_VecConv_0D(rtp,save_primitive(vr_:vp_,it,ip))
                    end if
                end do
            end do
        else
            if_inside_single=.false.
        end if
        
    end subroutine ModSave_Globe_single

    subroutine ModSave_Fan(Tree,p_save,nrt_out,filename,logical_unit)
        implicit none
        type(YYTree),intent(in),target      ::  Tree
        real,intent(in)                     ::  p_save
        integer,intent(in)                  ::  nrt_out(2)
        character(len=*),intent(in)         ::  filename
        integer,intent(in)                  ::  logical_unit
        integer                             ::  MpiRank,MpiSubRank,MpiSubSize

        integer                             ::  iLocalBlock
        type(BlockType),pointer             ::  Block1
        logical                             ::  if_inside_single
        integer                             ::  if_join,MPI_COMM_SAVE_SUBSET,ierr

        real                                ::  save_primitive_local(&
                                                nvar,nrt_out(1),nrt_out(2))
        real                                ::  save_primitive_global(&
                                                nvar,nrt_out(1),nrt_out(2))
        real                                ::  r_out(nrt_out(1)),t_out(nrt_out(2))
        integer                             ::  ir,it

        real                                ::  write_rt(2)
        
        ! MPI set up
        call MPi_Comm_rank(MPI_COMM_WORLD,MpiRank, ierr)

        ! Get the th and ph lists
        do ir=1,nrt_out(1)
            r_out(ir)=(r_range(1)*(nrt_out(1)-ir+0.5)+r_range(2)*(ir-0.5))/nrt_out(1)
        end do
        do it=1,nrt_out(2)
            t_out(it)=(it-0.50)/nrt_out(2)*dpi
        end do

        if_join=MPI_UNDEFINED

        save_primitive_local=-1.e30
        save_primitive_global=0.

        do iLocalBlock=1,size(Tree%LocalBlocks)
            Block1=>Tree%LocalBlocks(iLocalBlock)
            call ModSave_Fan_single(Block1,p_save,nrt_out,r_out,t_out,&
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
                nvar*nrt_out(1)*nrt_out(2),&
                MPI_REAL, MPI_MAX, 0, MPI_COMM_SAVE_SUBSET, ierr)
        end if

        ! Save it by the SubRank=0

        if (if_join==1 .and. MpiSubRank==0) then
            
            open(unit=logical_unit,file=filename, status='replace', action='write')
            write(logical_unit,*) nrt_out

            do ir=1,nrt_out(1)
                do it=1,nrt_out(2)
                    write_rt=[r_out(ir),t_out(it)]
                    write(logical_unit,*)write_rt,save_primitive_global(:,ir,it)
                end do
            end do

            close(logical_unit)
            print *,'File='//filename,', saved.'
        end if
    end subroutine ModSave_Fan

    subroutine ModSave_Fan_single(Block1,p_save,nrt_out,r_out,t_out,save_primitive,if_inside_single)
        implicit none

        type(BlockType),intent(in)          ::  Block1
        real,intent(in)                     ::  p_save
        integer,intent(in)                  ::  nrt_out(2)
        real,intent(in)                     ::  r_out(nrt_out(1)),t_out(nrt_out(2))
        real,intent(inout)                  ::  save_primitive(:,:,:)
        logical,intent(inout)               ::  if_inside_single

        integer                             ::  ir,it,rtp_posi_int(3)
        real                                ::  rtp(3),rtp_posi(3),rtp_weight(3)

        if_inside_single=.false.
        ! Loop each point in the map
        ! if Yang, then first convert point coord to Yang and then
        ! see if it's inside the block.
        do ir=1,nrt_out(1)
            do it=1,nrt_out(2)
                rtp=[r_out(ir),t_out(it),p_save]
                if (.not. Block1%if_yin) rtp=ModYinYang_CoordConv_0D(rtp)

                if ((Block1%xj_I(-ng+1)-rtp(2))*(Block1%xj_I(ng+nj)-rtp(2)).le.0.0 .and. &
                    (Block1%xk_I(-ng+1)-rtp(3))*(Block1%xk_I(ng+nk)-rtp(3)).le.0.0 .and. &
                    (Block1%xi_I(-ng+1)-rtp(1))*(Block1%xi_I(ng+ni)-rtp(1)).le.0.0) then

                    if_inside_single=.true.
                    
                    rtp_posi(1)=-ng+1.+(rtp(1)-Block1%xi_I(-ng+1))/Block1%dxi
                    rtp_posi(2)=-ng+1.+(rtp(2)-Block1%xj_I(-ng+1))/Block1%dxj
                    rtp_posi(3)=-ng+1.+(rtp(3)-Block1%xk_I(-ng+1))/Block1%dxk

                    rtp_posi_int=floor(rtp_posi)
                    rtp_weight=rtp_posi-rtp_posi_int

                    save_primitive(:,ir,it)=&
                        Block1%primitive_IV(:,rtp_posi_int(1)  ,rtp_posi_int(2)  ,rtp_posi_int(3)  )*&
                        (1.-rtp_weight(1))*(1.-rtp_weight(2))*(1.-rtp_weight(3))+&
                        Block1%primitive_IV(:,rtp_posi_int(1)+1,rtp_posi_int(2)  ,rtp_posi_int(3)  )*&
                        (   rtp_weight(1))*(1.-rtp_weight(2))*(1.-rtp_weight(3))+&
                        Block1%primitive_IV(:,rtp_posi_int(1)  ,rtp_posi_int(2)+1,rtp_posi_int(3)  )*&
                        (1.-rtp_weight(1))*(   rtp_weight(2))*(1.-rtp_weight(3))+&
                        Block1%primitive_IV(:,rtp_posi_int(1)  ,rtp_posi_int(2)  ,rtp_posi_int(3)+1)*&
                        (1.-rtp_weight(1))*(1.-rtp_weight(2))*(   rtp_weight(3))+&
                        Block1%primitive_IV(:,rtp_posi_int(1)+1,rtp_posi_int(2)+1,rtp_posi_int(3)  )*&
                        (   rtp_weight(1))*(   rtp_weight(2))*(1.-rtp_weight(3))+&
                        Block1%primitive_IV(:,rtp_posi_int(1)  ,rtp_posi_int(2)+1,rtp_posi_int(3)+1)*&
                        (1.-rtp_weight(1))*(   rtp_weight(2))*(   rtp_weight(3))+&
                        Block1%primitive_IV(:,rtp_posi_int(1)+1,rtp_posi_int(2)  ,rtp_posi_int(3)+1)*&
                        (   rtp_weight(1))*(1.-rtp_weight(2))*(   rtp_weight(3))+&
                        Block1%primitive_IV(:,rtp_posi_int(1)+1,rtp_posi_int(2)+1,rtp_posi_int(3)+1)*&
                        (   rtp_weight(1))*(   rtp_weight(2))*(   rtp_weight(3))
                    
                    if (.not. Block1%if_yin) save_primitive(vr_:vp_,ir,it)=&
                        ModYinYang_VecConv_0D(rtp,save_primitive(vr_:vp_,ir,it))
                
                end if
            end do
        end do
        
    end subroutine ModSave_Fan_single

end module ModSavePlot