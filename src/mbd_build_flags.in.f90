module mbd_build_flags

logical, parameter :: WITH_MPI = @WITH_MPI@
logical, parameter :: WITH_SCALAPACK = .false.
logical, parameter :: WITH_ELPA = .false.
integer, parameter :: LEGENDRE_PREC = 8

end module
