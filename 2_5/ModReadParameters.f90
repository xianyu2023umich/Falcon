module ModReadParameters

    use ModParameters,      only:   r_range,ni,nj,nk,ng,nvar,&
                                    ModelS_delta,ModelS_c_sound__CGS,&
                                    ModelS_rmax,ModelS_dc_type,ModelS_dc_rmax,ModelS_filename,&
                                    nStepsSavePlot,nthSavePlot,nphSavePlot,nSteps,rSave,CFL,&
                                    ModelS_heating_ratio,NameEquation,Initiation_type,&
                                    Initiation_type_index
    use ModStratification,  only:   ModStratification_DoAll
    use ModAMR,             only:   AMR_nLevels,AMR_r_ranges,AMR_rtp_if_divide
    use ModVariables,       only:   rho1_,vr_,vt_,vp_,br_,bt_,bp_,s1_

    contains

    subroutine ModReadParameters_read(filename,logical_unit)
        implicit none
        character(len=*),intent(in)     ::  filename            ! Parameter input filename
        character(len=256)              ::  line,var            ! Which line it is
        character(len=22)               ::  name_sub='ModReadParameters_read'
        integer,intent(in)              ::  logical_unit
        integer                         ::  ios                 ! For reading
        integer                         ::  AMR_iLevel          ! For reading AMR

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
                
                ! AMR grid
                case("#AMR")
                    read(logical_unit, *, iostat=ios) AMR_nLevels
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading AMR_nLevels"
                        stop 1
                    end if
                    if (AMR_nLevels<0) then
                        write(*,*) "Error from ",name_sub,": Incorrect AMR_nLevels=",AMR_nLevels
                        stop 1
                    end if
                    if (AMR_nLevels>16) then
                        write(*,*) "Error from ",name_sub,": Too large AMR_nLevels=",AMR_nLevels
                        stop 1
                    end if

                    ! For a proper AMR_nLevels loop to read the AMR range.
                    if (AMR_nLevels>=1) then
                        ! First allocate
                        allocate(AMR_r_ranges(2,AMR_nLevels))
                        allocate(AMR_rtp_if_divide(AMR_nLevels,3))

                        do AMR_iLevel=1,AMR_nLevels
                            read(logical_unit, *, iostat=ios) AMR_rtp_if_divide(AMR_iLevel,1:3)
                            if (ios/=0) then
                                write(*,*) "Error from ",name_sub,": Error reading AMR_if_divide_r(AMR_iLevel)"
                                stop 1
                            end if
                            read(logical_unit, *, iostat=ios) AMR_r_ranges(1,AMR_iLevel)
                            if (ios/=0) then
                                write(*,*) "Error from ",name_sub,": Error reading AMR_r_ranges(1,AMR_iLevel)"
                                stop 1
                            end if
                            read(logical_unit, *, iostat=ios) AMR_r_ranges(2,AMR_iLevel)
                            if (ios/=0) then
                                write(*,*) "Error from ",name_sub,": Error reading AMR_r_ranges(2,AMR_iLevel)"
                                stop 1
                            end if
                        end do
                    end if

                ! When to stop
                case("#STOP")
                    read(logical_unit, *, iostat=ios) nSteps
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading nSteps"
                        stop 1
                    end if

                case("#SAVEPLOT")
                    read(logical_unit, *, iostat=ios) nStepsSavePlot
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading nStepsSavePlot"
                        stop 1
                    end if
                    read(logical_unit, *, iostat=ios) rSave
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading rSave"
                        stop 1
                    end if
                    read(logical_unit, *, iostat=ios) nthSavePlot
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading nthSavePlot"
                        stop 1
                    end if
                    read(logical_unit, *, iostat=ios) nphSavePlot
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading nphSavePlot"
                        stop 1
                    end if

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
                    read(logical_unit, *, iostat=ios) NameEquation
                    select case(NameEquation)
                    case("MHD","mhd")
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
                        nvar=5
                        rho1_=1
                        vr_=2
                        vt_=3
                        vp_=4
                        s1_=5
                    end select

                case("#INITIATION")
                    read(logical_unit, *, iostat=ios) Initiation_type
                    select case(Initiation_type)
                    case("none","None","NONE")
                        Initiation_type_index=0
                    case("Harmonics","harmonics")
                        Initiation_type_index=1
                    end select
                end select
            end if
        end do

        close(logical_unit)
    end subroutine ModReadParameters_read
end module ModReadParameters
