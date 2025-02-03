module ModParameters

    integer             ::      MpiSize,MpiRank             !   MPI

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

    integer             ::      nStepsSavePlot
    integer             ::      nthSavePlot
    integer             ::      nphSavePlot
    real                ::      rSave

    integer             ::      nSteps                      !   Total n of steps
    real                ::      CFL

    character(len=100)  ::      NameEquation="HD"
end module ModParameters