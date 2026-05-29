Module ModControl

    integer         :: nMaxBlocksPerRank = 150

    ! Time stepping part: time and steps
    integer         :: iStep = 0        ! current step number
    real(8)         :: t     = 0.0d0    ! simulation time
    real(8)         :: dt    = 0.0d0    ! current timestep

    ! Whether the param.in is opened currently.
    logical         :: if_param_file_opened = .false.

    ! Status of run
    ! 0: initiation
    ! 1: running
    ! 2: hold
    integer         :: iStatus = 0
end Module ModControl