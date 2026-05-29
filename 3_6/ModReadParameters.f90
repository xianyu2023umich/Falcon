module ModReadParameters

    use ModConst,           only:   R_sun__CGS
    use ModParameters,      only:   mpirank,&
                                    r_range,r_range_Rsun,&
                                    ni,nj,nk,ng,nvar,&
                                    Artificial_heating_ratio,&
                                    nSteps,&
                                    CFL,&
                                    NameEquation,iEquation,&
                                    Initiation_type,Initiation_type_index,randVelocity_rms,&
                                    InitiationB_type,InitiationB_type_index,Bphi_uniform,&
                                    iGeometry,&
                                    DoCheck,&
                                    Plots,nPlots,PlotType,&
                                    rLevelInitial,&
                                    nAMRs,AMRs,AMRType,&
                                    if_do_echo,nStepsEcho,&
                                    if_involve_B,DivB_method,DivB_option
    use ModControl,         only:   if_param_file_opened
    use ModLookUpTable,     only:   ModLookUpTable_Read
    use ModEOS,             only:   ModEOS_init
    use ModOpacity,         only:   ModOpacity_init
    use ModStratification_new, only: ModStratification_new_Init

    contains

    subroutine ModReadParameters_read(filename,logical_unit)
        implicit none
        character(len=*),intent(in)     ::  filename            ! Parameter input filename
        character(len=256)              ::  line,var            ! Which line it is
        character(len=22)               ::  name_sub='ModReadParameters_read'
        integer,intent(in)              ::  logical_unit
        integer                         ::  ios                 ! For reading
        logical                         ::  if_close

        ! Default if_close to true.
        if_close=.true.

        ! Open the input file if not opened yet.
        if (.not. if_param_file_opened) then
            open(logical_unit, file=filename, status='old', action='read', iostat=ios)
            if (ios /= 0) then
                write(*,*) "Error from ",name_sub,": Error opening file: ", filename
                stop
            end if
            if_param_file_opened = .true.
        end if

        ! The main loop to read the file.
        do  
            read(logical_unit, *, iostat=ios) line
            if (ios /= 0) exit  ! End of file
            
            if (line(1:1) == "#") then              ! If starts with # then enter this command block
                select case (trim(adjustl(line)))   ! See which command it is

                ! Read the domain range. 
                case("#DOMAIN")
                    read(logical_unit, *, iostat=ios) r_range_Rsun(1)
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading r_range(1)"
                        stop 1
                    end if
                    read(logical_unit, *, iostat=ios) r_range_Rsun(2)
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading r_range(2)"
                        stop 1
                    end if

                r_range=r_range_Rsun*R_sun__CGS

                ! The grid within the Block
                ! Read ni,nj,nk,ng
                case("#GRID")
                    read(logical_unit, *, iostat=ios) ni
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading ni"
                        stop 1
                    end if
                    read(logical_unit, *, iostat=ios) nj
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading nj"
                        stop 1
                    end if
                    read(logical_unit, *, iostat=ios) nk
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading nk"
                        stop 1
                    end if
                    read(logical_unit, *, iostat=ios) ng
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading ng"
                        stop 1
                    end if

                case("#RADIALGRID")
                    read(logical_unit, *, iostat=ios) rLevelInitial
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading rLevelInitial"
                        stop 1
                    end if
                
                ! AMR grid
                case("#AMR")
                    call ModReadParameters_read_AMR(logical_unit)

                ! When to stop
                case("#STOP")
                    read(logical_unit, *, iostat=ios) nSteps
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading nSteps"
                        stop 1
                    end if

                case("#SAVEPLOT")
                    call ModReadParameters_read_SavePlot(logical_unit)
                    

                case("#TIMESTEPPING")
                    read(logical_unit, *, iostat=ios) CFL
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading CFL"
                        stop 1
                    end if

                case("#HEATING")
                    read(logical_unit, *, iostat=ios) Artificial_heating_ratio
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading Artificial_heating_ratio"
                        stop 1
                    end if

                ! Read three parameters about ModelS:
                ! 1. ModelS_delta 2. ModelS_c_sound__CGS 3. ModelS_filename
                case("#MODELS")
                    call ModReadParameters_read_Models(logical_unit)

                case("#EQUATION")
                    call ModReadParameters_read_Equation(logical_unit)

                case("#INITIATION")
                    read(logical_unit, *, iostat=ios) Initiation_type
                    select case(Initiation_type)
                    case("none","None","NONE")
                        Initiation_type_index=0
                    case("Harmonics","harmonics")
                        Initiation_type_index=1
                    case("RandomV","randomv","randomvelocity","RandomVelocity")
                        Initiation_type_index=3
                        read(logical_unit, *, iostat=ios) randVelocity_rms
                        if (ios/=0) then
                            write(*,*) "Error from ",name_sub,": Error reading randVelocity_rms"
                            stop 1
                        end if
                    end select

                case("#INITIATIONB")
                    read(logical_unit, *, iostat=ios) InitiationB_type
                    select case(InitiationB_type)
                    case("none","None","NONE")
                        InitiationB_type_index=0
                    case("Bp","bp","BP")
                        InitiationB_type_index=1
                        read(logical_unit, *, iostat=ios) Bphi_uniform
                        if (ios/=0) then
                            write(*,*) "Error from ",name_sub,": Error reading Bphi_uniform"
                            stop 1
                        end if
                    case default
                        write(*,*) "Error from ",name_sub,": Unknown InitiationB_type: ",trim(adjustl(InitiationB_type))
                        stop 1
                    end select

                ! Read the geometry type
                case("#GEOMETRY")
                    read(logical_unit, *, iostat=ios) iGeometry
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading iGeometry"
                        stop 1
                    end if

                case("#DOCHECK")
                    read(logical_unit, *, iostat=ios) DoCheck
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading DoCheck"
                        stop 1
                    end if

                case("#ECHO")
                    read(logical_unit, *, iostat=ios) if_do_echo
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading if_do_echo"
                        stop 1
                    end if

                    read(logical_unit, *, iostat=ios) nStepsEcho
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading nStepsEcho"
                        stop 1
                    end if
                
                case("#UPDATEB")
                    call ModReadParameters_read_UpdateB(logical_unit)
                    

                case ("#CHECKPOINT")
                    if_close = .false.
                    exit

                case default
                    write(*,*) "Error from ",name_sub,": Unknown command: ",trim(adjustl(line))
                    stop 1
                end select
            end if
        end do

        ! Several things to check at the end:
        ! 1. If iequation=0, e.g., hd, but if_involve_B is true, then turn it off
        if (iEquation==0 .and. if_involve_B) then
            write(*,*) "Warning from ",name_sub,": if_involve_B is true but iEquation=0 (HD). Turn off if_involve_B."
            if_involve_B=.false.
        end if

        ! If there's no checkpoint then it means no more to read.
        if(if_close) then
            close(logical_unit)
            if_param_file_opened = .False.
        end if
    end subroutine ModReadParameters_read

    subroutine ModReadParameters_read_UpdateB(logical_unit)
        implicit none
        character(len=31)               ::  name_sub='ModReadParameters_read_DivB'
        integer,intent(in)              ::  logical_unit
        integer                         ::  ios                 ! For reading

        read(logical_unit, *, iostat=ios) if_involve_B
        if (ios/=0) then
            write(*,*) "Error from ",name_sub,": Error reading if_involve_B"
            stop 1
        end if

        if (if_involve_B) then
            read(logical_unit, *, iostat=ios) DivB_method
            if (ios/=0) then
                write(*,*) "Error from ",name_sub,": Error reading DivB_method"
                stop 1
            end if

            select case(trim(adjustl(DivB_method)))
            case('None','none','NONE')
                DivB_option=0
            case('GLM','glm')
                DivB_option=1
            case('DivBSource','divbsource','DIVBSOURCE')
                DivB_option=2
            end select
        end if
    end subroutine ModReadParameters_read_UpdateB

    subroutine ModReadParameters_read_SavePlot(logical_unit)
        implicit none
        character(len=31)               ::  name_sub='ModReadParameters_read_SavePlot'
        integer,intent(in)              ::  logical_unit
        integer                         ::  ios                 ! For reading
        integer                         ::  iPlot               ! For reading Plots
        type(PlotType),pointer          ::  Plot1
        integer                         ::  idirection          ! For reading rtp_range_SavePlot

        ! Read nPlots
        read(logical_unit, *, iostat=ios) nPlots
        if (ios/=0) then
            write(*,*) "Error from ",name_sub,": Error reading nPlots"
            stop 1
        end if

        ! Allocate Plots
        allocate(Plots(nPlots))

        ! Read each plot
        do iPlot=1,nPlots
            Plot1=>Plots(iPlot)
            Plot1%logical_unit=10+iPlot

            ! First see which type it is
            read(logical_unit, *, iostat=ios) Plot1%charType
            if (ios/=0) then
                write(*,*) "Error from ",name_sub,": Error reading charType=",trim(adjustl(Plot1%charType))
                stop 1
            end if

            ! Then see if this type is supported. If so then set the iType
            ! and continue reading.
            select case(Plot1%charType)


            ! For the sphere, we need the nStepsSavePlot, rSave, nthSavePlot, nphSavePlot
            case("sphere","Sphere","SPHERE")
                Plot1%iType=0
                
                read(logical_unit, *, iostat=ios) Plot1%nStepsSavePlot
                if (ios/=0) then
                    write(*,*) "Error from ",name_sub,": Error reading nStepsSavePlot"
                    stop 1
                end if
                read(logical_unit, *, iostat=ios) Plot1%rtp_SavePlot(1)
                if (ios/=0) then
                    write(*,*) "Error from ",name_sub,": Error reading rSave"
                    stop 1
                end if
                read(logical_unit, *, iostat=ios) Plot1%nrtp_SavePlot(2)
                if (ios/=0) then
                    write(*,*) "Error from ",name_sub,": Error reading nthSavePlot"
                    stop 1
                end if
                read(logical_unit, *, iostat=ios) Plot1%nrtp_SavePlot(3)
                if (ios/=0) then
                    write(*,*) "Error from ",name_sub,": Error reading nphSavePlot"
                    stop 1
                end if

            ! For the fan, we need the nStepsSavePlot, r_range, nrSavePlot, nthSavePlot, phSave
            case("fan","Fan","FAN")
                Plot1%iType=1
                read(logical_unit, *, iostat=ios) Plot1%nStepsSavePlot
                if (ios/=0) then
                    write(*,*) "Error from ",name_sub,": Error reading nStepsSavePlot"
                    stop 1
                end if
                read(logical_unit, *, iostat=ios) Plot1%rtp_range_SavePlot(1,1)
                if (ios/=0) then
                    write(*,*) "Error from ",name_sub,": Error reading r_range_SavePlot(1)"
                    stop 1
                end if
                read(logical_unit, *, iostat=ios) Plot1%rtp_range_SavePlot(1,2)
                if (ios/=0) then
                    write(*,*) "Error from ",name_sub,": Error reading r_range_SavePlot(2)"
                    stop 1
                end if
                read(logical_unit, *, iostat=ios) Plot1%nrtp_SavePlot(1)
                if (ios/=0) then
                    write(*,*) "Error from ",name_sub,": Error reading nrSavePlot"
                    stop 1
                end if
                read(logical_unit, *, iostat=ios) Plot1%nrtp_SavePlot(2)
                if (ios/=0) then
                    write(*,*) "Error from ",name_sub,": Error reading nthSavePlot"
                    stop 1
                end if
                read(logical_unit, *, iostat=ios) Plot1%rtp_SavePlot(3)
                if (ios/=0) then
                    write(*,*) "Error from ",name_sub,": Error reading phSave"
                    stop 1
                end if

            ! For the cube, we need:
            ! nStepsSavePlot, 
            ! r_range, nrSavePlot
            ! th_range, nthSavePlot
            ! ph_range, nphSavePlot
            case("cube","Cube","CUBE")
                Plot1%iType=2
                read(logical_unit, *, iostat=ios) Plot1%nStepsSavePlot
                if (ios/=0) then
                    write(*,*) "Error from ",name_sub,": Error reading nStepsSavePlot"
                    stop 1
                end if
                do idirection=1,3
                    read(logical_unit, *, iostat=ios) Plot1%rtp_range_SavePlot(idirection,1)
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading rtp_range_SavePlot(",idirection,",1)"
                        stop 1
                    end if
                    read(logical_unit, *, iostat=ios) Plot1%rtp_range_SavePlot(idirection,2)
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading rtp_range_SavePlot(",idirection,",2)"
                        stop 1
                    end if
                    read(logical_unit, *, iostat=ios) Plot1%nrtp_SavePlot(idirection)
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading nrtp_SavePlot(",idirection,")"
                        stop 1
                    end if
                end do
                read(logical_unit, *, iostat=ios) Plot1%if_cube_grid_points
                if (ios/=0) then
                    write(*,*) "Error from ",name_sub,": Error reading if_cube_grid_points"
                    stop 1
                end if


            case("3D","3d",'AllCells','allcells')
                Plot1%iType=3
                read(logical_unit, *, iostat=ios) Plot1%nStepsSavePlot
                if (ios/=0) then
                    write(*,*) "Error from ",name_sub,": Error reading nStepsSavePlot"
                    stop 1
                end if
            case default
                write(*,*) "Error from ",name_sub,": Unknown plot type: ",trim(adjustl(Plots(iPlot)%charType))
                stop 1
            end select
        end do

    end subroutine ModReadParameters_read_SavePlot

    subroutine ModReadParameters_read_AMR(logical_unit)
        implicit none
        character(len=31)               ::  name_sub='ModReadParameters_read_AMR'
        integer,intent(in)              ::  logical_unit
        integer                         ::  ios                 ! For reading
        integer                         ::  iAMR                ! For reading AMR
        integer                         ::  idirection          ! For reading rtp_range
        type(AMRType),pointer           ::  AMR1
        
        ! Read nAMRs
        read(logical_unit, *, iostat=ios) nAMRs
        if (ios/=0) then
            write(*,*) "Error from ",name_sub,": Error reading nAMRs"
            stop 1
        end if
        if (nAMRs<0) then
            write(*,*) "Error from ",name_sub,": Incorrect nAMRs=",nAMRs
            stop 1
        end if
        if (nAMRs>16) then
            write(*,*) "Error from ",name_sub,": Too many nAMRs=",nAMRs
            stop 1
        end if

        ! For a proper nAMRs loop to read the AMR range.
        if (nAMRs>=1) then
            ! First allocate
            allocate(AMRs(nAMRs))

            ! Loop each AMR
            do iAMR=1,nAMRs
                AMR1=>AMRs(iAMR)

                ! First read the rtp_if_divide
                read(logical_unit, *, iostat=ios) AMR1%rtp_if_divide(1:3)
                if (ios/=0) then
                    write(*,*) "Error from ",name_sub,": Error reading AMR1%rtp_if_divide(1:3)"
                    stop 1
                end if

                ! Then read if_global
                read(logical_unit, *, iostat=ios) AMR1%if_global
                if (ios/=0) then
                    write(*,*) "Error from ",name_sub,": Error reading AMR1%if_global"
                    stop 1
                end if

                ! If it's global, then we only need the r_range. Otherwise, we need the rtp_range.
                if (AMR1%if_global) then
                    read(logical_unit, *, iostat=ios) AMR1%rtp_range_Rsun(1,1)
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading AMR1%rtp_range(1,1)"
                        stop 1
                    end if
                    read(logical_unit, *, iostat=ios) AMR1%rtp_range_Rsun(1,2)
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading AMR1%rtp_range(1,2)"
                        stop 1
                    end if
                else
                    do idirection=1,3
                        read(logical_unit, *, iostat=ios) AMR1%rtp_range_Rsun(idirection,1)
                        if (ios/=0) then
                            write(*,*) "Error from ",name_sub,": Error reading AMR1%rtp_range(",idirection,",1)"
                            stop 1
                        end if
                        read(logical_unit, *, iostat=ios) AMR1%rtp_range_Rsun(idirection,2)
                        if (ios/=0) then
                            write(*,*) "Error from ",name_sub,": Error reading AMR1%rtp_range(",idirection,",2)"
                            stop 1
                        end if
                    end do
                end if
            end do
        end if
    end subroutine ModReadParameters_read_AMR

    subroutine ModReadParameters_read_Equation(logical_unit)
        implicit none
        ! character(len=31)               ::  name_sub='ModReadParameters_read_Equation'
        integer,intent(in)              ::  logical_unit
        integer                         ::  ios                 ! For reading

        ! Read the equation type
        read(logical_unit, *, iostat=ios) NameEquation
        select case(NameEquation)
        case("MHD","mhd")
            iEquation=1
            nvar=9
        case("HD","hd")
            iEquation=0
            nvar=5
        end select
    end subroutine ModReadParameters_read_Equation

    subroutine ModReadParameters_read_Models(logical_unit)
        implicit none
        integer,intent(in)              ::  logical_unit
        integer                         ::  ios                 ! For reading
        character(len=256)              ::  eos_file
        character(len=256)              ::  entropy_file
        character(len=256)              ::  opacity_file

        read(logical_unit, *, iostat=ios) eos_file
        if (ios/=0) then
            write(*,*) "Error reading eos_file."
            stop 1
        end if
        read(logical_unit, *, iostat=ios) entropy_file
        if (ios/=0) then
            write(*,*) "Error reading entropy_file."
            stop 1
        end if
        read(logical_unit, *, iostat=ios) opacity_file
        if (ios/=0) then
            write(*,*) "Error reading opacity_file."
            stop 1
        end if

        call ModLookUpTable_Read(eos_file, logical_unit=2)
        call ModLookUpTable_Read(entropy_file, logical_unit=2)
        call ModLookUpTable_Read(opacity_file, logical_unit=2)

        call ModEOS_init
        call ModOpacity_init
        call ModStratification_new_Init
    end subroutine ModReadParameters_read_Models
end module ModReadParameters
