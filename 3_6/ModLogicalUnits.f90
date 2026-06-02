Module ModLogicalUnits
    ! Central registry of Fortran logical (I/O) unit numbers.
    ! Claim a new constant here before using a new unit elsewhere.
    !
    ! Reserved ranges:
    !   1        — test programs (passed explicitly, not a named constant)
    !   2        — lookup tables during initiation (EOS, entropy, opacity)
    !   11–20    — plot output files  (iUnit_plot_base + iPlot, up to 10 plots)
    !   21–40    — log  output files  (iUnit_log_base  + iLog,  up to 10 logs)
    !   42       — parameter file (PARAM.in)

    ! Reading parameter input files (PARAM.in), opened/closed around each read.
    integer, parameter :: iUnit_param_file   = 42

    ! Reading lookup table files (EOS, entropy, opacity) during initiation.
    ! Opened and closed before any plot files are active.
    integer, parameter :: iUnit_lookup_table = 2

    ! Base unit for plot output files.  Plot #i opens unit (iUnit_plot_base + i).
    integer, parameter :: iUnit_plot_base    = 10

    ! Base unit for log  output files.  Log  #i opens unit (iUnit_log_base  + i).
    integer, parameter :: iUnit_log_base     = 20

end Module ModLogicalUnits
