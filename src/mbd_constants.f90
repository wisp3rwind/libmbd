! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module mbd_constants
!! Constants used throughout.

implicit none

integer, parameter :: dp = kind(0.d0)
real(dp), parameter :: pi = acos(-1.d0)
real(dp), parameter :: ang = 1.8897259886d0
    !! Value of angstrom in atomic units

integer, parameter :: MBD_EXC_NEG_EIGVALS = 1
    !! Negative eigenvalue exception
integer, parameter :: MBD_EXC_NEG_POL = 2
    !! Negative polarizability exception
integer, parameter :: MBD_EXC_LINALG = 3
    !! Exception in LAPACK or ScaLAPACK
integer, parameter :: MBD_EXC_UNIMPL = 4
    !! Functionality is not implemented
integer, parameter :: MBD_EXC_DAMPING = 5
    !! Damping-function exception
integer, parameter :: MBD_EXC_INPUT = 6
    !! Invalid input

integer, parameter :: MBD_LOG_LVL_DEBUG = -1
integer, parameter :: MBD_LOG_LVL_INFO = 0
integer, parameter :: MBD_LOG_LVL_WARN = 1
integer, parameter :: MBD_LOG_LVL_ERROR = 2

real(dp), parameter :: ZERO_REAL = 0d0
complex(dp), parameter :: ZERO_COMPLEX = (0d0, 0d0)
complex(dp), parameter :: IMI = (0d0, 1d0)

integer, parameter :: DAMPING_UNKNOWN = 0
integer, parameter :: DAMPING_BARE = 1
integer, parameter :: DAMPING_DIP_1MEXP = 2
integer, parameter :: DAMPING_FERMI_DIP = 3
integer, parameter :: DAMPING_SQRTFERMI_DIP = 4
integer, parameter :: DAMPING_CUSTOM_DIP = 5
integer, parameter :: DAMPING_DIP_CUSTOM = 6
integer, parameter :: DAMPING_DIP_GG = 7
integer, parameter :: DAMPING_FERMI_DIP_GG = 8
integer, parameter :: DAMPING_SQRTFERMI_DIP_GG = 9
integer, parameter :: DAMPING_CUSTOM_DIP_GG = 10

end module
