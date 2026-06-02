module ModParameters

    type LogType
        character(len=100)  ::  VarName                     ! Volume or layers
        character(len=100)  ::  charType                    ! Volume or layers
        integer         ::      iType                       ! 0: volume, 1: layers
        integer         ::      nStepsSaveLog
        integer         ::      logical_unit
        real(8)         ::      r_range_SaveLog(2)          ! in cgs
        integer         ::      nr_SaveLog         
        real(8),allocatable ::  r_list_SaveLog(:)           ! in cgs         
        real(8),allocatable ::  data_SaveLog(:)             ! the data to be saved, which will be allocated in ModSaveLog_DoAll 
    end type LogType

    type PlotType
        character(len=100)  ::  charType
        integer         ::      iType
        integer         ::      nStepsSavePlot
        integer         ::      logical_unit
        logical         ::      if_cube_grid_points=.false.
        real(8)         ::      rtp_SavePlot(3)        
        real(8)         ::      rtp_range_SavePlot(3,2)    
        integer         ::      nrtp_SavePlot(3)          
    end type PlotType

    type AMRType
        integer         ::      iLevel
        logical         ::      rtp_if_divide(3)
        logical         ::      if_global
        real(8)         ::      rtp_range_Rsun(3,2)
    end type AMRType

    integer             ::      MpiSize,MpiRank             !   MPI

    integer             ::      iGeometry                   !   0 is cartesian, 1 is spherical
    real(8)             ::      r_range(2)                  !   Domain range
    real(8)             ::      r_range_Rsun(2)             !   Same thing but in Rsun

    integer             ::      ni                          !   Grid in a block
    integer             ::      nj
    integer             ::      nk
    integer             ::      ng                          !   Ghost cell grid
    integer             ::      nvar=5                      !   number of variables

    real(8)             ::      Artificial_heating_ratio=1.0

    integer             ::      nPlots=0
    type(PlotType),&
        allocatable,&
        target          ::      Plots(:)
    
    integer             ::      nLogs=0
    type(LogType),&
        allocatable,&
        target          ::      Logs(:)

    integer             ::      nSteps                      !   Total n of steps
    real(8)             ::      CFL

    character(len=100)  ::      NameEquation="HD"
    integer             ::      iEquation=0                 !   0: HD, 1: MHD, 2: Corona    

    character(len=100)  ::      Initiation_type
    integer             ::      Initiation_type_index=0
    real(8)             ::      randVelocity_rms=1.0d3

    character(len=100)  ::      InitiationB_type
    integer             ::      InitiationB_type_index=0
    real(8)             ::      Bphi_uniform=1.0d2

    integer             ::      rLevelInitial
    integer             ::      rLevelOption

    integer             ::      Multigrid_nLevels

    logical             ::      DoCheck=.false.

    integer             ::      nAMRs
    type(AMRType),&
        allocatable,&
        target          ::      AMRs(:)

    ! Coronal heating parameters
    real(8)             ::      LperpSqrtB=1.5e9

    logical             ::      if_involve_B
    character(len=100)  ::      DivB_method                 ! "GLM" or "CT"
    integer             ::      DivB_option=0               ! 0: nothing. 1: GLM. 2: Powell's source term.
    
    logical             ::      if_do_echo=.false.
    integer             ::      nStepsEcho=10

    
end module ModParameters
