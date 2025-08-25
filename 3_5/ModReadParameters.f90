module ModReadParameters

    use ModParameters,      only:   r_range,ni,nj,nk,ng,nvar,&
                                    ModelS_delta,ModelS_c_sound__CGS,&
                                    ModelS_rmax,ModelS_dc_type,ModelS_dc_rmax,ModelS_filename,&
                                    nSteps,CFL,ModelS_heating_ratio,NameEquation,iEquation,Initiation_type,&
                                    Initiation_type_index,rLevelInitial,iGeometry,DoCheck,&
                                    Plots,nPlots,PlotType,nAMRs,AMRs,AMRType
    use ModStratification,  only:   ModStratification_read_lookuptable,ModStratification_DoAll
    use ModVariables,       only:   rho1_,vr_,vt_,vp_,br_,bt_,bp_,s1_

    contains

    subroutine ModReadParameters_read(filename,logical_unit)
        implicit none
        character(len=*),intent(in)     ::  filename            ! Parameter input filename
        character(len=256)              ::  line,var            ! Which line it is
        character(len=22)               ::  name_sub='ModReadParameters_read'
        integer,intent(in)              ::  logical_unit
        integer                         ::  ios                 ! For reading

        ! Open the input file
        open(logical_unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) "Error from ",name_sub,": Error opening file: ", filename
            stop
        end if

        ! The main loop to read the file.
        do  
            read(logical_unit, '(A)', iostat=ios) line
            if (ios /= 0) exit  ! End of file
            
            if (line(1:1) == "#") then              ! If starts with # then enter this command block
                select case (trim(adjustl(line)))   ! See which command it is

                ! Read the domain range. 
                case("#DOMAIN")
                    read(logical_unit, *, iostat=ios) r_range(1)
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading r_range(1)"
                        stop 1
                    end if
                    read(logical_unit, *, iostat=ios) r_range(2)
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading r_range(2)"
                        stop 1
                    end if

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
                    
                ! About artificial cooling
                ! Read ModelS_rmax. Then read the type of dc
                ! If user, then read an extra line
                case("#ARTIFICIALCOOLING")
                    read(logical_unit, *, iostat=ios) ModelS_rmax
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading ModelS_rmax"
                        stop 1
                    end if

                    read(logical_unit, *, iostat=ios) ModelS_dc_type
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading ModelS_rmax"
                        stop 1
                    end if

                    ! After reading ModelS_dc_type,
                    ! if it's read then read one more time
                    ! If it's scale then calculate ModelS_dc_rmax later.
                    ! Otherwise just stop.
                    select case(ModelS_dc_type)
                    case("read","Read","READ")
                        read(logical_unit, *, iostat=ios) ModelS_dc_rmax
                        if (ios/=0) then
                            write(*,*) "Error from ",name_sub,": Error reading RandomScale"
                            stop 1
                        end if
                    case("scale")
                    case default
                        write(*,*) "Error from ",name_sub,": Unknown ModelS_dc_type=",trim(adjustl(var))
                        STOP
                    end select
                
                ! Read three parameters about ModelS:
                ! 1. ModelS_delta 2. ModelS_c_sound__CGS 3. ModelS_filename
                case("#MODELS")
                    read(logical_unit, *, iostat=ios) ModelS_delta
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading ModelS_delta"
                        stop 1
                    end if

                    read(logical_unit, *, iostat=ios) ModelS_c_sound__CGS
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading ModelS_c_sound__CGS"
                        stop 1
                    end if
                    
                    read(logical_unit, *, iostat=ios) ModelS_filename
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading ModelS_filename"
                        stop 1
                    end if

                    call ModStratification_read_lookuptable(ModelS_filename,logical_unit=2,head_size=1)
                    call ModStratification_DoAll
                
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
                    read(logical_unit, *, iostat=ios) ModelS_heating_ratio
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading ModelS_heating_ratio"
                        stop 1
                    end if

                case("#EQUATION")
                    call ModReadParameters_read_Equation(logical_unit)

                case("#INITIATION")
                    read(logical_unit, *, iostat=ios) Initiation_type
                    select case(Initiation_type)
                    case("none","None","NONE")
                        Initiation_type_index=0
                    case("Harmonics","harmonics")
                        Initiation_type_index=1
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

                case default
                    write(*,*) "Error from ",name_sub,": Unknown command: ",trim(adjustl(line))
                    stop 1
                end select
            end if
        end do

        close(logical_unit)
    end subroutine ModReadParameters_read

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
                    read(logical_unit, *, iostat=ios) AMR1%rtp_range(1,1)
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading AMR1%rtp_range(1,1)"
                        stop 1
                    end if
                    read(logical_unit, *, iostat=ios) AMR1%rtp_range(1,2)
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading AMR1%rtp_range(1,2)"
                        stop 1
                    end if
                else
                    do idirection=1,3
                        read(logical_unit, *, iostat=ios) AMR1%rtp_range(idirection,1)
                        if (ios/=0) then
                            write(*,*) "Error from ",name_sub,": Error reading AMR1%rtp_range(",idirection,",1)"
                            stop 1
                        end if
                        read(logical_unit, *, iostat=ios) AMR1%rtp_range(idirection,2)
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
        character(len=31)               ::  name_sub='ModReadParameters_read_Equation'
        integer,intent(in)              ::  logical_unit
        integer                         ::  ios                 ! For reading

        ! Read the equation type
        read(logical_unit, *, iostat=ios) NameEquation
        select case(NameEquation)
        case("MHD","mhd")
            iEquation=1
            nvar=8
            rho1_=1
            vr_=2
            vt_=3
            vp_=4
            br_=5
            bt_=6
            bp_=7
            s1_=8
        case("HD","hd")
            iEquation=0
            nvar=5
            rho1_=1
            vr_=2
            vt_=3
            vp_=4
            s1_=5
        end select
    end subroutine ModReadParameters_read_Equation
end module ModReadParameters
