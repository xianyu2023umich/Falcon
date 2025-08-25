module ModParameters

    type PlotType
        character(len=100)  ::  charType
        integer         ::      iType
        integer         ::      nStepsSavePlot
        integer         ::      logical_unit
        real            ::      rtp_SavePlot(3)        
        real            ::      rtp_range_SavePlot(3,2)    
        integer         ::      nrtp_SavePlot(3)          
    end type PlotType

    type AMRType
        integer         ::      iLevel
        logical         ::      rtp_if_divide(3)
        logical         ::      if_global
        real            ::      rtp_range(3,2)
    end type AMRType

    integer             ::      MpiSize,MpiRank             !   MPI

    integer             ::      iGeometry                   !   0 is cartesian, 1 is spherical
    real                ::      r_range(2)                  !   Domain range

    integer             ::      ni                          !   Grid in a block
    integer             ::      nj
    integer             ::      nk
    integer             ::      ng                          !   Ghost cell grid
    integer             ::      nvar=5                      !   number of variables

    real                ::      ModelS_delta=1.0e-6         !   Used for scales
    real                ::      ModelS_c_sound__CGS=1.8e5   !   constant c speed in rsst

    character(len=100)  ::      ModelS_dc_type              !   Type of 
    real                ::      ModelS_rmax                 !   Top of cooling
    real                ::      ModelS_dc_rmax=-1.0
    character(len=100)  ::      ModelS_filename
    real                ::      ModelS_heating_ratio=1.0

    integer             ::      nPlots=0
    type(PlotType),&
        allocatable,&
        target          ::  Plots(:)

    integer             ::      nSteps                      !   Total n of steps
    real                ::      CFL

    character(len=100)  ::      NameEquation="HD"
    integer             ::      iEquation=0                 !   0: HD, 1: MHD             
    character(len=100)  ::      Initiation_type
    integer             ::      Initiation_type_index=0

    integer             ::      rLevelInitial

    integer             ::      Multigrid_nLevels

    logical             ::      DoCheck=.false.

    integer             ::      nAMRs
    type(AMRType),&
        allocatable,&
        target          ::  AMRs(:)
    
end module ModParameters