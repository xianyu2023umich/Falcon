module ModParameters

    type PlotType
        character(len=100)  ::  charType
        integer         ::      iType
        integer         ::      nStepsSavePlot
        integer         ::      logical_unit
        logical             ::      if_cube_grid_points=.false.
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

    real(8)             ::      ModelS_delta=1.0e-6         !   Used for scales
    real(8)             ::      ModelS_c_sound__CGS=1.8e5   !   constant c speed in rsst

    character(len=100)  ::      ModelS_dc_type              !   Type of artificial cooling
    real(8)             ::      ModelS_rmax                 !   Top of cooling
    real(8)             ::      ModelS_dc_rmax=-1.0
    character(len=100)  ::      ModelS_filename
    real(8)             ::      ModelS_heating_ratio=1.0

    integer             ::      nPlots=0
    type(PlotType),&
        allocatable,&
        target          ::      Plots(:)

    integer             ::      nSteps                      !   Total n of steps
    real(8)             ::      CFL

    character(len=100)  ::      NameEquation="HD"
    integer             ::      iEquation=0                 !   0: HD, 1: MHD, 2: Corona         
    character(len=100)  ::      Initiation_type
    integer             ::      Initiation_type_index=0

    integer             ::      rLevelInitial

    integer             ::      Multigrid_nLevels

    logical             ::      DoCheck=.false.

    integer             ::      nAMRs
    type(AMRType),&
        allocatable,&
        target          ::      AMRs(:)

    ! Coronal heating parameters
    real(8)             ::      LperpSqrtB=1.5e9

    character(len=100)  ::      DivB_method                 ! "GLM" or "CT"
    integer             ::      DivB_option=0               ! 0: nothing. 1: GLM. 2: Powell's source term.
    
end module ModParameters
