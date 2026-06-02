module ModSaveLog

    use ModControl,     only:   iStep,t,dt
    use ModBlock,       only:   BlockType
    use ModYinYang,     only:   ModYinYang_CoordConv_0D,&
                                ModYinYang_VecConv_0D
    use ModYinYangTree, only:   YYTree
    use ModConst,       only:   dpi,R_sun__CGS
    use ModParameters,  only:   nvar,ni,nj,nk,ng,r_range,r_range_Rsun,&
                                iEquation,nLogs,Logs,LogType,MpiRank
    use ModMath,        only:   ModMath_IfLinesInterSect
    use ModEquation,    only:   ModEquation_Dynamo_Get_p1
    use MPI

    implicit none

    ! The flag determining whether initialized.
    ! Turn on after DoAll called first time. 
    logical             ::  if_initialized = .false.

    contains

    subroutine ModSaveLog_DoAll(Tree)
        implicit none
        type(YYTree),target     ::  Tree
        type(LogType),pointer   ::  Log1
        integer                 ::  iLog,ir
        character(len=100)      ::  filename
        logical                 ::  if_intersect
        integer                 ::  ierr

        ! The first time to call this subroutine
        ! we need to initialize.

        if (.not. if_initialized) then
            do iLog=1,nLogs
                Log1=>Logs(iLog)

                ! Initial set for each itype.
                select case (Log1%iType)
                case (0)
                case (1)
                    ! If the maximum of r_range is negative
                    ! then it means it reads r_range
                    if (maxval(Log1%r_range_Savelog)<=0.0d0) then
                        Log1%r_range_Savelog = r_range
                    end if

                    ! Get r_list
                    if (Log1%nr_SaveLog>1) then
                        allocate(Log1%r_list_SaveLog(Log1%nr_SaveLog))
                        
                        do ir=1,Log1%nr_SaveLog
                            Log1%r_list_SaveLog(ir) = Log1%r_range_Savelog(1)       +   &
                                (Log1%r_range_Savelog(2)-Log1%r_range_Savelog(1))   *   &
                                real(ir-1,8) / real(Log1%nr_SaveLog-1,8)
                        end do

                        allocate(Log1%data_SaveLog(Log1%nr_SaveLog))
                        Log1%data_SaveLog = 0.0d0
                    else
                        if (MpiRank==0) print *, 'Warning: &
                            nr_SaveLog should be larger than 1 for layer log.'
                    end if

                    ! Open the logical unit for each log at first time.
                    ! filename should be something like Log/Log_VarName.dat
                    if (iLog==1) call execute_command_line('mkdir -p Log', wait=.true.)
                    filename = 'Log/Log_'//trim(adjustl(Log1%VarName))//'.dat'
                    open(unit=Log1%logical_unit,file=filename, &
                        status='replace', action='write')
                    
                    ! Write nr at first line and r_list at the second line for later use.
                    ! Third name is 'iStep dt t {Varname}'
                    write(Log1%logical_unit,*) Log1%nr_SaveLog
                    write(Log1%logical_unit,*) Log1%r_list_SaveLog
                    write(Log1%logical_unit,*) 'iStep dt t '//trim(adjustl(Log1%VarName))
                case (2)
                end select
            end do

            ! Set the flag to true after initialization.
            if_initialized = .true.
        end if

        ! Now call the actual saving.
        call ModSaveLog_AllLayers(Tree)
    end subroutine ModSaveLog_DoAll

    ! Close all logical units for logs at the end of the program.
    subroutine ModSaveLog_CloseAll()
        implicit none
        type(LogType),pointer   ::  Log1
        integer                 ::  iLog

        do iLog=1,nLogs
            Log1=>Logs(iLog)
            close(Log1%logical_unit)
        end do
    end subroutine ModSaveLog_CloseAll

    subroutine ModSaveLog_AllLayers(Tree)
        implicit none
        type(YYTree),target     ::  Tree
        integer                 ::  iLog
        type(LogType),pointer   ::  Log1
        integer                 ::  iLocalBlock
        type(BlockType),pointer ::  Block1
        real(8),allocatable     ::  data_SaveLog(:)
        logical                 ::  if_intersect
        integer                 ::  ierr

        ! Outer loop: logs — so reset/reduce/write happen once per log per timestep.
        do iLog=1,nLogs
            Log1=>Logs(iLog)

            if (Log1%iType==1 .and. mod(iStep,Log1%nStepsSaveLog)==0) then
                Log1%data_SaveLog(:)=0.0d0

                ! Inner loop: accumulate contributions from all local blocks.
                do iLocalBlock=1,Tree%nLocalBlocks
                    Block1=>Tree%LocalBlocks(iLocalBlock)
                    if_intersect = ModMath_IfLinesInterSect(&
                        Block1%xijk_range(1,:),Log1%r_range_Savelog)
                    if (if_intersect) call ModSaveLog_OneLayer_SingleBlock(Block1, Log1)
                end do

                ! Reduce across all ranks, then rank 0 writes.
                allocate(data_SaveLog(Log1%nr_SaveLog))
                call MPI_Reduce(Log1%data_SaveLog, data_SaveLog, &
                    Log1%nr_SaveLog, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
                if (MpiRank==0) write(Log1%logical_unit,*) iStep,dt,t,data_SaveLog
                deallocate(data_SaveLog)
            end if
        end do
    end subroutine ModSaveLog_AllLayers

    subroutine ModSaveLog_OneLayer_SingleBlock(Block1, Log1)
        implicit none
        type(BlockType),pointer ::  Block1
        type(LogType),pointer   ::  Log1
        real(8)                 ::  r,weight
        integer                 ::  ir_pos_int,ir,j,k
        real(8)                 ::  T0,p0,rho0,e1,p1,enthapy,dS
        real(8)                 ::  kinetic_energy,radiative_flux,cooling_flux
        real(8)                 ::  primitive_I(nvar)

        ! Loop each ir to see whether
        ! it is within the block's r range.
        do ir=1,Log1%nr_SaveLog
            if (Log1%r_list_SaveLog(ir)>=Block1%xijk_range(1,1) .and. &
                Log1%r_list_SaveLog(ir)<=Block1%xijk_range(1,2)) then
                
                ! get r, and the position index of this r
                r=Log1%r_list_SaveLog(ir)
                ir_pos_int=int((r-Block1%xi_I(1))/Block1%dxi)+1

                ! weight is the fraction part
                weight=(r-Block1%xi_I(1))/Block1%dxi+1-ir_pos_int

                ! Select which to calculate
                select case (Log1%VarName)
                case ('Le')
                    ! Perturbed enthalpy flux
                    ! = \int (rho0 * e1 + p1 - p0 * rho1/rho0) v_r dS
                    ! Which means we need e1.
                    ! Since de=T ds - p d(1/rho)
                    ! = T ds + p/rho^2 drho, we have e1=T s1 + p/rho^2 rho1
                    ! The coefficients T and p/rho^2 can be estimated by
                    ! background values at this r.

                    T0 = Block1%te0_I(ir_pos_int)*(1.0d0-weight)+&
                        Block1%te0_I(ir_pos_int+1)*weight
                    p0 = Block1%p0_I(ir_pos_int)*(1.0d0-weight)+&
                        Block1%p0_I(ir_pos_int+1)*weight
                    rho0 = Block1%rho0_I(ir_pos_int)*(1.0d0-weight)+&
                        Block1%rho0_I(ir_pos_int+1)*weight
                    
                    ! Call calculation of p1
                    Block1%primitive=>Block1%primitive_IV
                    call ModEquation_Dynamo_Get_p1(Block1)

                    ! Loop all nj*nk points at this ir.
                    do j=1,nj
                        do k=1,nk
                            ! get primitive variables at this point by interpolation in r.
                            primitive_I(:) = &
                                Block1%primitive_IV(ir_pos_int,j,k,:)*(1.0d0-weight)+&
                                Block1%primitive_IV(ir_pos_int+1,j,k,:)*weight

                            ! e1=T s1 + p/rho^2 rho1
                            e1 =  T0*primitive_I(Block1%s1_)   &
                                + p0/rho0**2              &
                                * primitive_I(Block1%rho1_)

                            p1 = Block1%p1_III(ir_pos_int,j,k)*(1.0d0-weight)+&
                                Block1%p1_III(ir_pos_int+1,j,k)*weight
                            
                            ! enthalpy = rho0 * e1 + p1 - p0 * rho1/rho0
                            enthapy = rho0 * e1 + p1  &
                                - p0 * primitive_I(Block1%rho1_)/rho0
                            
                            ! Get the area of this piece at r.
                            dS = r ** 2 * (cos(Block1%xj_F(j))-cos(Block1%xj_F(j+1)))*&
                                (Block1%xk_F(k+1)-Block1%xk_F(k))
                            
                            ! Add this piece to data.
                            log1%data_SaveLog(ir) = log1%data_SaveLog(ir) &
                                + enthapy * primitive_I(Block1%vr_)*dS
                        end do
                    end do
                case ('Lk')
                    ! Kinetic energy flux
                    ! = \int 0.5 rho0 v^2 v_r dS
                    ! First get background rho0 at this r.
                    rho0 = Block1%rho0_I(ir_pos_int)*(1.0d0-weight)+&
                        Block1%rho0_I(ir_pos_int+1)*weight
                    
                    do j=1,nj
                        do k=1,nk
                            ! First get primitive
                            primitive_I(:) = &
                                Block1%primitive_IV(ir_pos_int,j,k,:)*(1.0d0-weight)+&
                                Block1%primitive_IV(ir_pos_int+1,j,k,:)*weight
                            
                            ! Get kinetic energy
                            kinetic_energy = 0.5d0 * rho0 * (primitive_I(Block1%vr_)**2 + &
                                primitive_I(Block1%vt_)**2 + &
                                primitive_I(Block1%vp_)**2)
                            
                            ! Get the area of this piece at r.
                            dS = r ** 2 * (cos(Block1%xj_F(j))-cos(Block1%xj_F(j+1)))*&
                                (Block1%xk_F(k+1)-Block1%xk_F(k))
                            
                            ! Add this piece to data.
                            log1%data_SaveLog(ir) = log1%data_SaveLog(ir) &
                                + kinetic_energy * primitive_I(Block1%vr_)*dS
                        end do
                    end do
                case ('Lr')
                    ! Since the radiative flux is uniform, 
                    ! just do one interpolation and get the
                    ! total area and then get the flux.

                    ! The radiative flux
                    radiative_flux = Block1%diffusion_flux_I(ir_pos_int)*(1.0d0-weight)+&
                        Block1%diffusion_flux_I(ir_pos_int+1)*(weight)
                    
                    ! Get the total area of this block at r.
                    dS = r ** 2 * (cos(Block1%xj_F(1))-cos(Block1%xj_F(nj+1)))*&
                        (Block1%xk_F(nk+1)-Block1%xk_F(1))
                    
                    ! Add this piece to data.
                    log1%data_SaveLog(ir) = log1%data_SaveLog(ir) &
                        + radiative_flux * dS
                case ('Ls')
                    ! Ls is similar to Lr.
                    ! The cooling flux
                    cooling_flux = Block1%cooling_flux_I(ir_pos_int)*(1.0d0-weight)+&
                        Block1%cooling_flux_I(ir_pos_int+1)*(weight)
                    
                    ! Get the total area of this block at r.
                    dS = r ** 2 * (cos(Block1%xj_F(1))-cos(Block1%xj_F(nj+1)))*&
                        (Block1%xk_F(nk+1)-Block1%xk_F(1))
                    
                    ! Add this piece to data.
                    log1%data_SaveLog(ir) = log1%data_SaveLog(ir) &
                        + cooling_flux * dS
                end select
            end if
        end do
    end subroutine ModSaveLog_OneLayer_SingleBlock


end module ModSaveLog