! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
program mbd_grad_bench

use mbd_constants
use mbd_geom, only: geom_t
use mbd_benches
#ifdef WITH_MPI
use mbd_mpi
#endif

implicit none

character(len=:), allocatable :: bench_name, fname
integer :: repetitions = 0, N = 0
logical :: periodic = .false.

#ifdef WITH_MPI
integer :: rank, err
#endif

call parse_args()

#ifdef WITH_MPI
call MPI_INIT(err)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, err)
#endif

call exec_bench(bench_name, repetitions, N, fname)

#ifdef WITH_MPI
call MPI_FINALIZE(err)
#endif

contains

subroutine read_arg(pos, buf)
    integer, intent(inout) :: pos
    character(len=:), allocatable, intent(inout) :: buf

    integer :: length, status

    if (.not. allocated(buf)) then
        allocate(character(len=1) :: buf)
    end if

    call get_command_argument(pos, buf, length, status)
    if (status <= 0) then
        deallocate(buf)
        allocate(character(len=length) :: buf)
        call get_command_argument(pos, buf, length, status)
    end if
    if (status /= 0) then
        stop 2
    end if

    pos = pos + 1
end subroutine

subroutine parse_args()
    integer :: count_, pos
    character(len=:), allocatable :: buf

    count_ = command_argument_count()
    pos = 1
    do
        call read_arg(pos, buf)

        select case (buf)
        case ("--name")
            call read_arg(pos, bench_name)
        case ("-N")
            call read_arg(pos, buf)
            read(buf, *) N
        case ("--repetitions")
            call read_arg(pos, buf)
            read(buf, *) repetitions
        case ("--periodic")
            periodic = .true.
        case ("--fname")
            call read_arg(pos, fname)
        case default 
            print *, 'Unknown argument:'
            print *, buf
            stop 1
        end select

        if (pos > count_) then
            exit
        end if
    end do

    if (N <= 0) then
        print *, "Invalid or missing value for -N"
        stop 1
    end if

    if (repetitions <= 0) then
        print *, "Invalid or missing value for --repetitions"
        stop 1
    end if

    if (.not. allocated(bench_name)) then
        print *, "Missing value for --name"
        stop 1
    end if
end subroutine

subroutine exec_bench(bench_name, repetitions, N, fname)
    implicit none

    character(len=*), intent(in) :: bench_name
    integer, intent(in) :: repetitions, N
    character(len=:), allocatable, intent(in) :: fname
    type(geom_t) :: geom

    real(dp) :: res
    logical :: is_initialized
    integer :: i

    is_initialized = .false.

    !allocate (geom)
    geom%log%level = MBD_LOG_LVL_DEBUG
    geom%parallel_mode = 'none'

    do i = 1, repetitions
        select case (bench_name)
        case ('plain')
            res = res + bench_mbd_plain(geom, N, periodic, is_initialized)
        case ('rsscs')
            res = res + bench_mbd_rsscs(geom, N, periodic, is_initialized)
        case ('T-bare')
            res = res + bench_mbd_T(geom, N, periodic, is_initialized, 'fermi,dip')
        case ('T-erf')
            res = res + bench_mbd_T(geom, N, periodic, is_initialized, 'fermi,dip,gg')
        case default
            print *, 'Unknown benchmark!'
            exit
        end select

        if (geom%has_exc()) then
            print *, 'Exception!'
            print *, geom%exc%origin
            print *, geom%exc%msg
            print *, geom%exc%code
            exit
        end if

        is_initialized = .true.
    end do

    if (is_initialized) then
        call geom%timer%print()
        if (allocated(fname)) then
            call geom%timer%print_csv(fname)
        end if

        ! Be sure that the result is used and the compiler cannot just
        ! remove the calculation (do the same for forces!).
    end if
    print '(F16.6)', res

    call geom%destroy()

    if (geom%has_exc()) then
        !deallocate (geom)
        stop 1
    end if
    !deallocate (geom)
end subroutine

end program

