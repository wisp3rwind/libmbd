! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_benches

use mbd_constants
use mbd_damping, only: damping_t, damping_fermi
use mbd_dipole, only: dipole_matrix, T_bare, T_erf_coulomb, damping_grad, T_erfc
use mbd_geom, only: geom_t
use mbd_gradients, only: grad_t, grad_matrix_re_t, grad_request_t, grad_scalar_t
use mbd_hamiltonian, only: get_mbd_hamiltonian_energy
use mbd_matrix, only: matrix_re_t, matrix_cplx_t
use mbd_methods, only: get_mbd_scs_energy
use mbd_ts, only: get_ts_energy
use mbd_scs, only: run_scs
use mbd_utils, only: diff7, findval, tostr, result_t

implicit none

contains

subroutine build_coords(geom, N, periodic)
    type(geom_t), intent(inout) :: geom
    logical, intent(in) :: periodic
    integer, intent(in) :: N

    integer :: n_atoms, i, j
    real(dp) :: r0
    real(dp), allocatable :: coords(:, :)

    ! quadratic sheet in x-y plane
    n_atoms = N * N
    r0 = 4d0
    allocate (coords(3, n_atoms), source=0d0)
    do i = 1, N
        do j = 1, N
            coords(1, (i - 1) * N + j) = (i - 1) * r0
            coords(2, (i - 1) * N + j) = (j - 1) * r0
            coords(3, (i - 1) * N + j) = 1.23
        end do
    end do
    geom%coords = coords

    ! stack squares along z axis
    if (periodic) then
        if (.not. allocated(geom%lattice)) then
            allocate(geom%lattice(3, 3), source=0d0)
        end if
        geom%lattice(:, 1) = [real(dp) :: N * r0, 0_dp, 0_dp]
        geom%lattice(:, 2) = [real(dp) :: 0_dp, N * r0, 0_dp]
        geom%lattice(:, 3) = [real(dp) :: 0_dp, 0_dp, r0]

        ! FIXME: Determine in a sensible way for this system
        geom%k_grid = [5, 5, 5]
    end if
end subroutine


real(dp) function bench_mbd_plain(geom, N, periodic, is_initialized) result(E)
    type(geom_t), intent(inout) :: geom
    integer, intent(in) :: N
    logical, intent(in) :: periodic
    logical, intent(in) :: is_initialized

    type(damping_t) :: damp
    real(dp), allocatable :: alpha_0(:), omega(:)
    integer :: n_atoms
    type(result_t) :: res

    call build_coords(geom, N, periodic)
    if (.not. is_initialized) then
        call geom%init()
    end if

    damp%version = 'fermi,dip'

    n_atoms = geom%siz()

    allocate(damp%r_vdw(n_atoms), source=3.55d0)
    damp%beta = 0.83d0
    allocate(alpha_0(n_atoms), source=12d0)
    allocate(omega(n_atoms), source=.7d0)

    res = get_mbd_hamiltonian_energy(geom, alpha_0, omega, damp, &
        grad_request_t(dcoords=.true.))
    E = res%energy
end function


real(dp) function bench_mbd_rsscs(geom, N, periodic, is_initialized) result(E)
    type(geom_t), intent(inout) :: geom
    integer, intent(in) :: N
    logical, intent(in) :: periodic
    logical, intent(in) :: is_initialized

    type(damping_t) :: damp
    real(dp), allocatable :: alpha_0(:), omega(:), C6(:)
    integer :: n_atoms
    type(result_t) :: res

    call build_coords(geom, N, periodic)
    if (.not. is_initialized) then
        call geom%init()
    end if

    damp%version = 'fermi,dip'

    n_atoms = geom%siz()

    allocate(damp%r_vdw(n_atoms), source=3.55d0)
    damp%beta = 0.83d0
    allocate(alpha_0(n_atoms), source=12d0)
    allocate(omega(n_atoms), source=.7d0)
    allocate(C6(n_atoms), source=65d0)

    res = get_mbd_scs_energy(geom, 'rsscs', alpha_0, C6, damp, &
        grad_request_t(dcoords=.true.))
    E = res%energy
end function

real(dp) function bench_mbd_T(geom, N, periodic, is_initialized, version) result(res)
    type(geom_t), intent(inout) :: geom
    integer, intent(in) :: N
    logical, intent(in) :: periodic
    logical, intent(in) :: is_initialized
    character(len=*), intent(in) :: version

    type(damping_t) :: damp
    real(dp), allocatable :: alpha_0(:), omega(:)
    integer :: n_atoms
    type(matrix_re_t) :: T
    type(grad_matrix_re_t) :: dT

    call build_coords(geom, N, periodic)
    if (.not. is_initialized) then
        call geom%init()
    end if

    damp%version = version

    n_atoms = geom%siz()

    allocate(damp%r_vdw(n_atoms), source=3.55d0)
    damp%beta = 0.83d0
    allocate(alpha_0(n_atoms), source=12d0)
    allocate(omega(n_atoms), source=.7d0)

    T = dipole_matrix(geom, damp, dT, grad_request_t(dcoords=.true.))
    res = 0.0
end function

end
