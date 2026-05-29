module ModSavePlot

    use ieee_arithmetic, only:  ieee_is_nan
    use ModBlock,       only:   BlockType
    use ModYinYang,     only:   ModYinYang_CoordConv_0D,&
                                ModYinYang_VecConv_0D
    
    use ModYinYangTree, only:   YYTree
    use ModConst,       only:   dpi,R_sun__CGS
    use ModParameters,  only:   nvar,ni,nj,nk,ng,r_range,r_range_Rsun,PlotType,Plots,nPlots,iEquation
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
                    call ModSave_Globe(Tree,&
                        Plot1%rtp_SavePlot(1),[Plot1%nrtp_SavePlot(2),Plot1%nrtp_SavePlot(3)],&
                        'Sphere_'//iStep_char//'.dat',Plot1%logical_unit)
                case(1)
                    call ModSave_Fan(Tree,&
                        Plot1%rtp_SavePlot(3),[Plot1%nrtp_SavePlot(1),Plot1%nrtp_SavePlot(2)],&
                        'Fan_'//iStep_char//'.dat',Plot1%logical_unit)
                case(2)
                    call ModSave_Cube(Tree,&
                        Plot1%rtp_range_SavePlot,Plot1%nrtp_SavePlot,&
                        Plot1%if_cube_grid_points,&
                        'Cube_'//iStep_char//'.dat',Plot1%logical_unit)
                case(3)
                    call ModSave_AllCells(Tree,iStep,Plot1%logical_unit)
                end select
            end if
        end do
    end subroutine ModSave_DoAll

    !subroutine ModSave_Energy_Flux(Tree,nr_out,filename,logical_unit)
    !    implicit none
    !    type(YYTree),intent(in),target      ::  Tree
    !    integer,intent(in)                  ::  nr_out
    !    character(len=*),intent(in)         ::  filename
    !    integer,intent(in)                  ::  logical_unit
    !    integer                             ::  MpiRank,MpiSubRank,MpiSubSize
    !    integer                             ::  iLocalBlock
    !    type(BlockType),pointer             ::  Block1
    !    logical                             ::  if_inside_single





    !end subroutine ModSave_Energy_Flux

    ! It's better to use MPI_Reduce to save the data.

    subroutine ModSave_Cube(Tree,rtp_range_Rsun,nrtp_out,if_grid_points,filename,logical_unit)
        implicit none
        type(YYTree),intent(in),target      ::  Tree
        real(8),intent(in)                  ::  rtp_range_Rsun(3,2)
        integer,intent(in)                  ::  nrtp_out(3)
        logical,intent(in)                  ::  if_grid_points
        character(len=*),intent(in)         ::  filename
        integer,intent(in)                  ::  logical_unit
        integer                             ::  nrtp_save(3)
        integer                             ::  iLocalBlock
        type(BlockType),pointer             ::  Block1

        integer                             ::  MpiRank,MpiSubRank,MpiSubSize
        logical                             ::  if_inside_single
        integer                             ::  if_join,MPI_COMM_SAVE_SUBSET,ierr

        integer                             ::  ir,it,ip
        real(8),allocatable                 ::  save_primitive_layer_III_local(:,:,:),&
                                                save_primitive_layer_III_global(:,:,:)
        real(8),allocatable                 ::  save_primitive_IV(:,:,:,:)
        real(8)                             ::  write_xyz_Rsun(3)
        real(8)                             ::  r_out(nrtp_out(1))
        real(8),allocatable                 ::  t_out(:),p_out(:)
        character(len=100)                  ::  line_format

        nrtp_save = nrtp_out
        if (if_grid_points) then
            nrtp_save(2) = nrtp_out(2) + 1
            nrtp_save(3) = nrtp_out(3) + 1
        end if

        allocate(save_primitive_layer_III_local(nvar,nrtp_save(2),nrtp_save(3)))
        allocate(save_primitive_layer_III_global(nvar,nrtp_save(2),nrtp_save(3)))
        allocate(t_out(nrtp_save(2)))
        allocate(p_out(nrtp_save(3)))

        ! MPI set up
        call MPi_Comm_rank(MPI_COMM_WORLD,MpiRank, ierr)

        ! Get the r, th and ph lists
        do ir=1,nrtp_out(1)
            r_out(ir)=((ir-1.0)*rtp_range_Rsun(1,2)*R_sun__CGS+&
                (nrtp_out(1)-ir)*rtp_range_Rsun(1,1)*R_sun__CGS)/(nrtp_out(1)-1.0)
        end do
        do it=1,nrtp_save(2)
            if (if_grid_points) then
                t_out(it)=rtp_range_Rsun(2,1)+(it-1.0d0)/(nrtp_save(2)-1.0d0)*(rtp_range_Rsun(2,2)-rtp_range_Rsun(2,1))
            else
                t_out(it)=rtp_range_Rsun(2,1)+((it-0.50)/nrtp_save(2))*(rtp_range_Rsun(2,2)-rtp_range_Rsun(2,1))
            end if
        end do
        do ip=1,nrtp_save(3)
            if (if_grid_points) then
                p_out(ip)=rtp_range_Rsun(3,1)+(ip-1.0d0)/(nrtp_save(3)-1.0d0)*(rtp_range_Rsun(3,2)-rtp_range_Rsun(3,1))
            else
                p_out(ip)=rtp_range_Rsun(3,1)+((ip-0.50)/nrtp_save(3))*(rtp_range_Rsun(3,2)-rtp_range_Rsun(3,1))
            end if
        end do

        ! Allocate the IV at the first layer by rank 0.
        if (MpiRank==0) then
            allocate(save_primitive_IV(nvar,nrtp_out(1),nrtp_save(2),nrtp_save(3)))
            save_primitive_IV=-1.e30
        end if

        ! Mpi_reduce each layer so that each time it's not that large.
        do ir=1,nrtp_out(1)

            ! Initialize the variables
            ! Rank 0 ALWAYS join, because it's the rank to perform saving.
            if_join=MPI_UNDEFINED
            if (MpiRank==0) if_join=1

            save_primitive_layer_III_local=-1.e30
            save_primitive_layer_III_global=0.

            ! Loop each block to 1. fulfill the layer III and 2. see if it's inside the block.
            do iLocalBlock=1,size(Tree%LocalBlocks)
                Block1=>Tree%LocalBlocks(iLocalBlock)
                call ModSave_Globe_single(Block1,r_out(ir),nrtp_save(2:3),t_out,p_out,&
                    save_primitive_layer_III_local,if_inside_single)
                
                if (if_inside_single) if_join=1
            end do

            ! Get the new communicator

            call MPI_Comm_split(MPI_COMM_WORLD, if_join, MpiRank, MPI_COMM_SAVE_SUBSET, ierr)

            ! Reduce this layer to the root

            if (if_join==1) then
                call MPI_Comm_rank(MPI_COMM_SAVE_SUBSET,MpiSubRank,ierr)
                call MPI_Comm_size(MPI_COMM_SAVE_SUBSET,MpiSubSize,ierr)

                call MPI_Reduce(save_primitive_layer_III_local,save_primitive_layer_III_global,&
                    nvar*nrtp_save(2)*nrtp_save(3),&
                    MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_SAVE_SUBSET, ierr)

                ! Save the data to the IV
                if (MpiRank==0 .and. allocated(save_primitive_IV)) &
                    save_primitive_IV(:,ir,:,:)=save_primitive_layer_III_global

                ! Free the communicator
                call MPI_Comm_free(MPI_COMM_SAVE_SUBSET,ierr)
            end if            
        end do

        if (MpiRank==0 .and. allocated(save_primitive_IV)) then
            if (if_grid_points) then
                nrtp_save(2) = nrtp_out(2) + 1
                nrtp_save(3) = nrtp_out(3) + 1
            end if
            open(unit=logical_unit,file=filename, status='replace', action='write')
            write(logical_unit,*) 'TITLE = "Cube"'
            write(logical_unit,*) 'VARIABLES = "K", "J", "I", "X [Rs]", "Y [Rs]", "Z [Rs]", "Rho1", "Vr", "Vt", "Vp", "S1"'
            write(logical_unit,*) 'ZONE T="STRUCTURE GRID", I=', nrtp_out(1),&
             ', J=', nrtp_save(2), ', K=', nrtp_save(3), ', F=POINT'

            select case(iEquation)
            case(0)
                line_format='(3I5,8ES16.8)'
            case(1)
                line_format='(3I5,11ES16.8)'
            end select

            !print *,save_primitive_IV(:,nrtp_out(1),nrtp_save(2),nrtp_save(3))
            
            do ip=1,nrtp_save(3)
                do it=1,nrtp_save(2)
                    do ir=1,nrtp_out(1)
                        write_xyz_Rsun(1)=r_out(ir)*sin(t_out(it))*cos(p_out(ip))/R_sun__CGS
                        write_xyz_Rsun(2)=r_out(ir)*sin(t_out(it))*sin(p_out(ip))/R_sun__CGS
                        write_xyz_Rsun(3)=r_out(ir)*cos(t_out(it))/R_sun__CGS

                        write(logical_unit,line_format) ir,it,ip,write_xyz_Rsun,save_primitive_IV(:,ir,it,ip)
                    end do
                end do
            end do
            close(logical_unit)
            print *,'File='//filename,', saved.'
        end if
    end subroutine ModSave_Cube

    subroutine ModSave_Globe(Tree,r_save_Rsun,ntp_out,filename,logical_unit)
        implicit none
        type(YYTree),intent(in),target      ::  Tree
        real(8),intent(in)                  ::  r_save_Rsun
        real(8)                             ::  r_save_cgs
        integer,intent(in)                  ::  ntp_out(2)
        character(len=*),intent(in)         ::  filename
        integer,intent(in)                  ::  logical_unit
        integer                             ::  MpiRank,MpiSubRank,MpiSubSize

        integer                             ::  iLocalBlock
        type(BlockType),pointer             ::  Block1
        logical                             ::  if_inside_single
        integer                             ::  if_join,MPI_COMM_SAVE_SUBSET,ierr

        real(8)                             ::  save_primitive_local(&
                                                nvar,ntp_out(1),ntp_out(2))
        real(8)                             ::  save_primitive_global(&
                                                nvar,ntp_out(1),ntp_out(2))
        real(8)                             ::  t_out(ntp_out(1)),p_out(ntp_out(2))
        integer                             ::  it,ip

        real(8)                             ::  write_tp(2),write_xyz_Rsun(3)

        ! unit change
        r_save_cgs=r_save_Rsun*R_sun__CGS
        
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
            call ModSave_Globe_single(Block1,r_save_cgs,ntp_out,t_out,p_out,&
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
                MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_SAVE_SUBSET, ierr)
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
                    write_xyz_Rsun(1)=r_save_Rsun*sin(t_out(it))*cos(p_out(ip))
                    write_xyz_Rsun(2)=r_save_Rsun*sin(t_out(it))*sin(p_out(ip))
                    write_xyz_Rsun(3)=r_save_Rsun*cos(t_out(it))

                    write(logical_unit,'(3I5,8ES16.8)') ip,it,1,write_xyz_Rsun,save_primitive_global(:,it,ip)
                    write(logical_unit,'(3I5,8ES16.8)') ip,it,2,write_xyz_Rsun,save_primitive_global(:,it,ip)
                end do
            end do
            close(logical_unit)
            print *,'File='//filename,', saved.'
        end if
    end subroutine ModSave_Globe

    subroutine ModSave_Globe_single(Block1,r_save,ntp_out,t_out,p_out,save_primitive,if_inside_single)
        implicit none

        type(BlockType),intent(in)          ::  Block1
        real(8),intent(in)                  ::  r_save
        integer,intent(in)                  ::  ntp_out(2)
        real(8),intent(in)                  ::  t_out(ntp_out(1)),p_out(ntp_out(2))
        real(8),intent(inout)               ::  save_primitive(:,:,:)
        logical,intent(inout)               ::  if_inside_single

        integer                             ::  it,ip,rtp_posi_int(3)
        real(8)                             ::  rtp(3),rtp_posi(3),rtp_weight(3)
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
                        if (.not. Block1%if_yin) save_primitive(Block1%vr_:Block1%vp_,it,ip)=&
                            ModYinYang_VecConv_0D(rtp,save_primitive(Block1%vr_:Block1%vp_,it,ip))
                        
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
        real(8),intent(in)                  ::  p_save
        integer,intent(in)                  ::  nrt_out(2)
        character(len=*),intent(in)         ::  filename
        integer,intent(in)                  ::  logical_unit
        integer                             ::  MpiRank,MpiSubRank,MpiSubSize

        integer                             ::  iLocalBlock
        type(BlockType),pointer             ::  Block1
        logical                             ::  if_inside_single
        integer                             ::  if_join,MPI_COMM_SAVE_SUBSET,ierr

        real(8)                             ::  save_primitive_local(&
                                                nvar,nrt_out(1),nrt_out(2))
        real(8)                             ::  save_primitive_global(&
                                                nvar,nrt_out(1),nrt_out(2))
        real(8)                             ::  r_out(nrt_out(1)),t_out(nrt_out(2))
        integer                             ::  ir,it

        real(8)                             ::  write_rt(2)
        
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
                MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_SAVE_SUBSET, ierr)
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
        real(8),intent(in)                  ::  p_save
        integer,intent(in)                  ::  nrt_out(2)
        real(8),intent(in)                  ::  r_out(nrt_out(1)),t_out(nrt_out(2))
        real(8),intent(inout)               ::  save_primitive(:,:,:)
        logical,intent(inout)               ::  if_inside_single

        integer                             ::  ir,it,rtp_posi_int(3)
        real(8)                             ::  rtp(3),rtp_posi(3),rtp_weight(3)
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
                    
                    if (.not. Block1%if_yin) save_primitive(Block1%vr_:Block1%vp_,ir,it)=&
                        ModYinYang_VecConv_0D(rtp,save_primitive(Block1%vr_:Block1%vp_,ir,it))
                
                end if
            end do
        end do
        
    end subroutine ModSave_Fan_single

    subroutine ModSave_AllCells(Tree, iStep, logical_unit)
        implicit none
        type(YYTree),intent(in),target      ::  Tree
        integer,intent(in)                  ::  iStep
        integer,intent(in)                  ::  logical_unit

        character(len=8)                    ::  iStep_char, iBlock_char
        character(len=30)                   ::  dirname
        character(len=60)                   ::  filename
        integer                             ::  iLocalBlock, i, j, k, MpiRank, ierr
        type(BlockType),pointer             ::  Block1
        real(8)                             ::  r, theta, phi, sth, cth, sph, cph
        real(8)                             ::  xyz(3), vx0, vy0, vz0
        real(8)                             ::  vars_out(nvar)

        ! Build zero-padded step string
        write(iStep_char,'(I8)') iStep
        do i = 1, 8
            if (iStep_char(i:i)==' ') iStep_char(i:i)='0'
        end do
        dirname = 'AllCells_' // iStep_char

        ! Rank 0 creates the directory, then all ranks wait
        call MPI_Comm_rank(MPI_COMM_WORLD, MpiRank, ierr)
        if (MpiRank==0) call execute_command_line('mkdir -p ' // trim(dirname))
        call MPI_Barrier(MPI_COMM_WORLD, ierr)

        ! Rank 0 writes a metadata file with global block count and block dimensions
        if (MpiRank==0) then
            open(unit=logical_unit, file=trim(dirname)//'/info.dat', status='replace', action='write')
            write(logical_unit,'(A,I0)') 'nBlocks = ', Tree%NumLeafNodes
            write(logical_unit,'(A,I0)') 'ni = ', ni
            write(logical_unit,'(A,I0)') 'nj = ', nj
            write(logical_unit,'(A,I0)') 'nk = ', nk
            write(logical_unit,'(A,I0)') 'nVarsTotal = ', nvar + 3
            select case(iEquation)
            case(0)
                write(logical_unit,'(A)') 'VARIABLES = "X [Rs]", "Y [Rs]", "Z [Rs]", '//&
                    '"Rho1", "Vx", "Vy", "Vz", "S1"'
            case(1)
                write(logical_unit,'(A)') 'VARIABLES = "X [Rs]", "Y [Rs]", "Z [Rs]", '//&
                    '"Rho1", "Vx", "Vy", "Vz", "Bx", "By", "Bz", "S1", "PHI"'
            case(2)
                write(logical_unit,'(A)') 'VARIABLES = "X [Rs]", "Y [Rs]", "Z [Rs]", '//&
                    '"Rho", "Vx", "Vy", "Vz", "Bx", "By", "Bz", "Te", "W_plus", "W_minus"'
            end select
            close(logical_unit)
        end if

        ! Loop over blocks owned by this rank
        do iLocalBlock = 1, size(Tree%LocalBlocks)
            Block1 => Tree%LocalBlocks(iLocalBlock)

            ! Allocate and compute corner values by averaging 8 surrounding cell centers
            call ModSavePlot_AllocateCorner(Block1)
            do k = 1, nk+1
                do j = 1, nj+1
                    do i = 1, ni+1
                        Block1%primitive_corner_IV(i,j,k,:) = 0.125d0 * ( &
                            Block1%primitive_IV(i-1,j-1,k-1,:) + Block1%primitive_IV(i  ,j-1,k-1,:) + &
                            Block1%primitive_IV(i-1,j  ,k-1,:) + Block1%primitive_IV(i  ,j  ,k-1,:) + &
                            Block1%primitive_IV(i-1,j-1,k  ,:) + Block1%primitive_IV(i  ,j-1,k  ,:) + &
                            Block1%primitive_IV(i-1,j  ,k  ,:) + Block1%primitive_IV(i  ,j  ,k  ,:) )
                    end do
                end do
            end do

            ! Open the file for this block
            write(iBlock_char,'(I8)') Block1%iBlock
            do i = 1, 8
                if (iBlock_char(i:i)==' ') iBlock_char(i:i)='0'
            end do
            filename = trim(dirname) // '/Block_' // iBlock_char // '.dat'
            open(unit=logical_unit, file=filename, status='replace', action='write')

            do k = 1, nk+1
                do j = 1, nj+1
                    do i = 1, ni+1
                        r     = Block1%xi_F(i)
                        theta = Block1%xj_F(j)
                        phi   = Block1%xk_F(k)
                        sth = sin(theta);  cth = cos(theta)
                        sph = sin(phi);    cph = cos(phi)

                        ! Position: rtp → Cartesian xyz in units of Rs
                        ! Yin: standard;  Yang: (x,y,z) = (-x0,z0,y0)
                        if (Block1%if_yin) then
                            xyz(1) = r*sth*cph
                            xyz(2) = r*sth*sph
                            xyz(3) = r*cth
                        else
                            xyz(1) = -r*sth*cph
                            xyz(2) = r*cth
                            xyz(3) = r*sth*sph
                        end if
                        xyz = xyz / R_sun__CGS

                        vars_out = Block1%primitive_corner_IV(i,j,k,:)

                        ! Velocity: (vr,vt,vp) → (vx,vy,vz)
                        vx0 = vars_out(Block1%vr_)*sth*cph + vars_out(Block1%vt_)*cth*cph - vars_out(Block1%vp_)*sph
                        vy0 = vars_out(Block1%vr_)*sth*sph + vars_out(Block1%vt_)*cth*sph + vars_out(Block1%vp_)*cph
                        vz0 = vars_out(Block1%vr_)*cth     - vars_out(Block1%vt_)*sth
                        if (Block1%if_yin) then
                            vars_out(Block1%vr_) = vx0;  vars_out(Block1%vt_) = vy0;  vars_out(Block1%vp_) = vz0
                        else
                            vars_out(Block1%vr_) = -vx0;  vars_out(Block1%vt_) = vz0;  vars_out(Block1%vp_) = vy0
                        end if

                        ! B field: (br,bt,bp) → (bx,by,bz), only if present
                        if (Block1%br_ > 0) then
                            vx0 = vars_out(Block1%br_)*sth*cph + vars_out(Block1%bt_)*cth*cph - vars_out(Block1%bp_)*sph
                            vy0 = vars_out(Block1%br_)*sth*sph + vars_out(Block1%bt_)*cth*sph + vars_out(Block1%bp_)*cph
                            vz0 = vars_out(Block1%br_)*cth     - vars_out(Block1%bt_)*sth
                            if (Block1%if_yin) then
                                vars_out(Block1%br_) = vx0;  vars_out(Block1%bt_) = vy0;  vars_out(Block1%bp_) = vz0
                            else
                                vars_out(Block1%br_) = -vx0;  vars_out(Block1%bt_) = vz0;  vars_out(Block1%bp_) = vy0
                            end if
                        end if

                        write(logical_unit,'(*(ES16.8))') xyz, vars_out
                    end do
                end do
            end do

            close(logical_unit)
            call ModSavePlot_DeallocateCorner(Block1)
        end do

    end subroutine ModSave_AllCells

    subroutine ModSavePlot_AllocateCorner(Block1)
        implicit none
        type(BlockType),intent(inout)   ::  Block1
        allocate(Block1%primitive_corner_IV(1:ni+1,1:nj+1,1:nk+1,1:nvar))
        Block1%primitive_corner_IV = 0.d0
    end subroutine ModSavePlot_AllocateCorner

    subroutine ModSavePlot_DeallocateCorner(Block1)
        implicit none
        type(BlockType),intent(inout)   ::  Block1
        if (allocated(Block1%primitive_corner_IV)) deallocate(Block1%primitive_corner_IV)
    end subroutine ModSavePlot_DeallocateCorner

end module ModSavePlot
