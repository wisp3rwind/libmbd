! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DO_COMPLEX_TYPE

module mbd_dipole
!! Construction of dipole tensors and dipole matrices.

use mbd_constants
use mbd_matrix, only: matrix_re_t, matrix_cplx_t
use mbd_geom, only: geom_t, supercell_circum
use mbd_damping, only: damping_t, damping_fermi, damping_sqrtfermi, &
    op1minus_grad
use mbd_gradients, only: grad_t, grad_matrix_re_t, grad_matrix_cplx_t, &
    grad_scalar_t, grad_request_t
use mbd_lapack, only: eigvals, inverse
use mbd_linalg, only: outer
use mbd_tensor, only: tensor_3x3_t, tensor_3x3_grad_scalar_t, tensor_3x3_grad_vector_t
use mbd_utils, only: tostr, shift_idx

implicit none

private
public :: dipole_matrix, T_bare, T_erf_coulomb, damping_grad, T_erfc

interface dipole_matrix
    !! Form either a real or a complex dipole matrix.
    !!
    !! The real-typed version is equivalent to \(\mathbf q=0\).
    !!
    !! $$
    !! \boldsymbol{\mathcal A}:=[\mathbf a_1\mathbf a_2\mathbf
    !! a_3],\qquad\boldsymbol{\mathcal B}:=[\mathbf b_1\mathbf b_2\mathbf b_3]
    !! \\ \boldsymbol{\mathcal B}=2\pi(\boldsymbol{\mathcal A}^{-1})^\mathrm
    !! T,\qquad \partial\boldsymbol{\mathcal B}=-\big((\partial\boldsymbol{\mathcal
    !! A})\boldsymbol{\mathcal A}^{-1}\big)^\mathrm T\boldsymbol{\mathcal B}
    !! \\ \mathbf R_\mathbf n=\boldsymbol{\mathcal A}\mathbf
    !! n,\qquad\partial\mathbf R_\mathbf n=(\partial\boldsymbol{\mathcal
    !! A})\mathbf n,
    !! \\ \mathbf G_\mathbf m=\boldsymbol{\mathcal B}\mathbf m,\qquad
    !! \partial\mathbf G_\mathbf m=-\big((\partial\boldsymbol{\mathcal
    !! A})\boldsymbol{\mathcal A}^{-1}\big)^\mathrm T\mathbf G_\mathbf m,
    !! \\ \frac{\partial G_{\mathbf ma}}{\partial A_{bc}}=-\mathcal A^{-1}_{ca}G_{\mathbf mb}
    !! $$
    !!
    !! $$
    !! \begin{gathered}
    !! \mathbf T_{ij}(\mathbf q)=\sum_{\mathbf n}\mathbf T(\mathbf R_{\mathbf
    !! nij})\mathrm e^{-\mathrm i\mathbf q\cdot\mathbf R_{\mathbf nij}},\quad\mathbf
    !! R_{\mathbf nij}=\mathbf R_j+\mathbf R_\mathbf n-\mathbf R_i
    !! \\ \frac{\mathrm d\mathbf R_{\mathbf nij}}{\mathrm d\mathbf
    !! R_k}=(\delta_{jk}-\delta_{ik})\mathbf I
    !! \\ \mathbf{T}_{ij}(\mathbf{q})\approx\mathbf{T}^\text{Ew}_{ij}(\mathbf{q})
    !! =\sum_\mathbf n^{|\mathbf R_{\mathbf nij}|<R_\text c}\mathbf
    !! T^\text{erfc}(\mathbf R_{\mathbf nij};\gamma)\mathrm e^{-\mathrm i\mathbf
    !! q\cdot\mathbf R_{\mathbf nij}} +\frac{4\pi}{V_\text{uc}}\sum_{\mathbf
    !! m}^{0<|\mathbf k_\mathbf m|<k_\text c}\mathbf{\hat k}_\mathbf
    !! m\otimes\mathbf{\hat k}_\mathbf m\,\mathrm e^{-\frac{k_\mathbf
    !! m^2}{4\gamma^2}-\mathrm i\mathbf G_\mathbf m\cdot\mathbf R_{ij}}
    !! \\ -\frac{4\gamma^3}{3\sqrt\pi}\delta_{ij}\mathbf I +\delta(\mathbf q)\frac{4
    !! \pi}{3 V_\text{uc}}\mathbf I,\qquad \mathbf k_\mathbf m=\mathbf G_\mathbf
    !! m+\mathbf q
    !! \end{gathered}
    !! $$
    !!
    !!
    !! $$
    !! \partial\left(\frac{4\pi}{V_\text{uc}}\right)=-\frac{4\pi}{V_\text{uc}}\frac{\partial
    !! V_\text{uc}}{V_\text{uc}},\qquad\frac{\partial
    !! V_\text{uc}}{\partial\boldsymbol{\mathcal
    !! A}}=V_\text{uc}(\boldsymbol{\mathcal A}^{-1})^\mathrm T
    !! \\ \partial(k^2)=2\mathbf k\cdot\partial\mathbf k
    !! \\ \mathbf{\hat k}\otimes\partial\mathbf{\hat k}=\frac{\mathbf
    !! k\otimes\partial\mathbf k}{k^2}-\frac{\mathbf k\otimes\mathbf
    !! k}{2k^4}\partial(k^2)
    !! $$
    !!
    !! $$
    !! \gamma:=\frac{2.5}{\sqrt[3]{V_\text{uc}}},\quad R_\text
    !! c:=\frac6\gamma,\quad k_\text c:=10\gamma
    !! $$
    module procedure dipole_matrix_real
    module procedure dipole_matrix_complex
end interface

contains

function parse_damping_version(damp) result(version)
    type(damping_t), intent(in) :: damp
    integer :: version

    select case (damp%version)
        case ("bare")
            version = DAMPING_BARE
        case ("dip,1mexp")
            version = DAMPING_DIP_1MEXP
        case ("fermi,dip")
            version = DAMPING_FERMI_DIP
        case ("sqrtfermi,dip")
            version = DAMPING_SQRTFERMI_DIP
        case ("custom,dip")
            version = DAMPING_CUSTOM_DIP
        case ("dip,custom")
            version = DAMPING_DIP_CUSTOM
        case ("dip,gg")
            version = DAMPING_DIP_GG
        case ("fermi,dip,gg")
            version = DAMPING_FERMI_DIP_GG
        case ("sqrtfermi,dip,gg")
            version = DAMPING_SQRTFERMI_DIP_GG
        case ("custom,dip,gg")
            version = DAMPING_CUSTOM_DIP_GG
        case default
            version = DAMPING_UNKNOWN
    end select
end function

#endif

#ifndef DO_COMPLEX_TYPE
type(matrix_re_t) function dipole_matrix_real( &
        geom, damp, ddipmat, grad) result(dipmat)
    use mbd_constants, only: ZERO => ZERO_REAL
#else
type(matrix_cplx_t) function dipole_matrix_complex( &
        geom, damp, ddipmat, grad, q) result(dipmat)
    use mbd_constants, only: ZERO => ZERO_COMPLEX
#endif

    type(geom_t), intent(inout) :: geom
    type(damping_t), intent(in) :: damp
    type(grad_request_t), intent(in), optional :: grad
#ifndef DO_COMPLEX_TYPE
    type(grad_matrix_re_t), intent(out), optional :: ddipmat
#else
    type(grad_matrix_cplx_t), intent(out), optional :: ddipmat
    real(dp), intent(in) :: q(3)
#endif

    real(dp) :: Rn(3), Rnij(3), Rnij_norm, T(3, 3), f_damp, &
        sigma_ij, T0(3, 3), beta_R_vdw
    integer :: i_atom, j_atom, i_cell, n(3), range_n(3), i, j, &
        my_i_atom, my_j_atom, i_latt, my_nr, my_nc, &
        damping_version
    logical :: do_ewald, is_periodic
    type(grad_matrix_re_t) :: dT, dT0, dTew
    type(grad_scalar_t) :: df
    type(grad_request_t) :: grad_ij
#ifndef DO_COMPLEX_TYPE
    real(dp) :: Tij(3, 3)
    type(grad_matrix_re_t) :: dTij
#else
    real(dp) :: qR
    complex(dp) :: Tij(3, 3), exp_qR
    type(grad_matrix_cplx_t) :: dTij
#endif

    ! Allocate dipole matrix
#ifdef WITH_SCALAPACK
    call dipmat%init(geom%idx, geom%blacs)
#else
    call dipmat%init(geom%idx)
#endif

    ! Check & prepare periodic calculation
    do_ewald = .false.
    is_periodic = allocated(geom%lattice)
    if (is_periodic) then
        do_ewald = geom%gamm > 0d0
        range_n = supercell_circum(geom%lattice, geom%real_space_cutoff)
    else
        range_n(:) = 0
    end if

    ! Prepare gradient arrays
    if (present(grad)) then
        grad_ij = grad
        grad_ij%dcoords = grad%dcoords .or. grad%dlattice
    end if
    if (grad_ij%dcoords) then
        allocate (dTij%dr(3, 3, 3))
        allocate (dT%dr(3, 3, 3))
        allocate (dT0%dr(3, 3, 3))
        allocate (dTew%dr(3, 3, 3))
        allocate (df%dr(3))
    end if
    if (grad_ij%dr_vdw) allocate (dT%dvdw(3, 3), source=0d0)
    my_nr = size(dipmat%idx%i_atom)
    my_nc = size(dipmat%idx%j_atom)
    allocate (dipmat%val(3 * my_nr, 3 * my_nc), source=ZERO)
    if (present(grad)) then
        if (grad%dcoords) then
            allocate (ddipmat%dr(3 * my_nr, 3 * my_nc, 3), source=ZERO)
        end if
        if (grad%dlattice) then
            allocate (ddipmat%dlattice(3 * my_nr, 3 * my_nc, 3, 3), source=ZERO)
        end if
        if (grad%dr_vdw) then
            allocate (ddipmat%dvdw(3 * my_nr, 3 * my_nc), source=ZERO)
            allocate (dTij%dvdw(3, 3))
        end if
        if (grad%dsigma) then
            allocate (ddipmat%dsigma(3 * my_nr, 3 * my_nc), source=ZERO)
            allocate (dTij%dsigma(3, 3))
        end if
#ifdef DO_COMPLEX_TYPE
        if (grad%dq) then
            allocate (ddipmat%dq(3 * my_nr, 3 * my_nc, 3), source=ZERO)
            allocate (dTij%dq(3, 3, 3))
        end if
#endif
    end if

    ! Build dipole matrix, including real-space summation if periodic
    damping_version = parse_damping_version(damp)
    call geom%clock(11)
    n = [0, 0, -1]
    each_cell: do i_cell = 1, product(1 + 2 * range_n)
        call shift_idx(n, -range_n, range_n)

        if (is_periodic) then
            Rn = matmul(geom%lattice, n)
        else
            Rn(:) = 0d0
        end if

        each_atom: do my_i_atom = 1, size(dipmat%idx%i_atom)
            i_atom = dipmat%idx%i_atom(my_i_atom)
            each_atom_pair: do my_j_atom = 1, size(dipmat%idx%j_atom)
                j_atom = dipmat%idx%j_atom(my_j_atom)

                ! i_cell == 1 => n = [0, 0, 0] => (i == j both atoms identical)
                if (i_cell == 1) then
                    if (i_atom == j_atom) cycle
                end if

                Rnij = geom%coords(:, i_atom) - geom%coords(:, j_atom) - Rn
                Rnij_norm = sqrt(sum(Rnij**2))
                if (is_periodic .and. Rnij_norm > geom%real_space_cutoff) cycle
                if (allocated(damp%R_vdw)) then
                    beta_R_vdw = damp%beta * sum(damp%R_vdw([i_atom, j_atom]))
                end if
                if (allocated(damp%sigma)) then
                    sigma_ij = damp%mayer_scaling &
                        * sqrt(sum(damp%sigma([i_atom, j_atom])**2))
                end if

                ! Evaluate i,j dipole interaction
                select case (damping_version)
                    case (DAMPING_BARE)
                        T = T_bare(Rnij, dT, grad_ij%dcoords)
                    case (DAMPING_DIP_1MEXP)
                        T = T_1mexp_coulomb(Rnij, beta_R_vdw, damp%a)
                    case (DAMPING_FERMI_DIP)
                        f_damp = damping_fermi(Rnij, beta_R_vdw, damp%a, df, grad_ij)
                        T0 = T_bare(Rnij, dT0, grad_ij%dcoords)
                        T = damping_grad(f_damp, df, T0, dT0, dT, grad_ij)
                    case (DAMPING_SQRTFERMI_DIP)
                        f_damp = damping_sqrtfermi(Rnij, beta_R_vdw, damp%a, df, grad_ij)
                        T0 = T_bare(Rnij, dT0, grad_ij%dcoords)
                        T = damping_grad(f_damp, df, T0, dT0, dT, grad_ij)
                    case (DAMPING_CUSTOM_DIP)
                        T = damp%damping_custom(i_atom, j_atom) * T_bare(Rnij)
                    case (DAMPING_DIP_CUSTOM)
                        T = damp%potential_custom(:, :, i_atom, j_atom)
                    case (DAMPING_DIP_GG)
                        T = T_erf_coulomb(Rnij, sigma_ij, dT, grad_ij)
                    case (DAMPING_FERMI_DIP_GG)
                        f_damp = damping_fermi(Rnij, beta_R_vdw, damp%a, df, grad_ij)
                        call op1minus_grad(f_damp, df)
                        T0 = T_erf_coulomb(Rnij, sigma_ij, dT0, grad_ij)
                        T = damping_grad(f_damp, df, T0, dT0, dT, grad_ij)
                        ! No Ewald summation for short-range interaction
                        do_ewald = .false.
                    case (DAMPING_SQRTFERMI_DIP_GG)
                        f_damp = damping_sqrtfermi(Rnij, beta_R_vdw, damp%a, df, grad_ij)
                        call op1minus_grad(f_damp, df)
                        T0 = T_erf_coulomb(Rnij, sigma_ij, dT0, grad_ij)
                        T = damping_grad(f_damp, df, T0, dT0, dT, grad_ij)
                        ! No Ewald summation for short-range interaction
                        do_ewald = .false.
                    case (DAMPING_CUSTOM_DIP_GG)
                        T = (1d0 - damp%damping_custom(i_atom, j_atom)) * &
                            T_erf_coulomb(Rnij, sigma_ij)
                        ! No Ewald summation for short-range interaction
                        do_ewald = .false.
                end select
                if (grad_ij%dr_vdw) dT%dvdw = damp%beta * dT%dvdw

                ! Cut off bare long-range interaction if using Ewald
                if (do_ewald) then
                    T = T + T_erfc_minus_bare(Rnij, geom%gamm, dTew, grad_ij)
                    if (grad_ij%dcoords) dT%dr = dT%dr + dTew%dr
                end if

                ! Copy over to potentially complex-valued tensor
                Tij = T
                if (grad_ij%dcoords) dTij%dr = dT%dr
                if (grad_ij%dr_vdw) dTij%dvdw = dT%dvdw
                if (grad_ij%dsigma) dTij%dsigma = dT%dsigma

                ! Multiply exp(-i q * Rnij), update gradients accordingly
#ifdef DO_COMPLEX_TYPE
                qR = dot_product(q, Rnij)
                exp_qR = complex(cos(qR), -sin(qR))
                Tij = T * exp_qR
                if (grad_ij%dcoords) then
#ifndef WITHOUT_DO_CONCURRENT
                    do concurrent(i=1:3)
#else
                    do i = 1, 3
#endif
                        dTij%dr(:, :, i) = dT%dr(:, :, i) * exp_qR - IMI * q(i) * Tij
                    end do
                end if
                if (grad_ij%dsigma) dTij%dsigma = dT%dsigma * exp_qR
                if (grad_ij%dr_vdw) dTij%dvdw = dT%dvdw * exp_qR
                if (grad_ij%dq) then
                    do concurrent(i=1:3)
                        dTij%dq(:, :, i) = -IMI * Rnij(i) * Tij
                    end do
                end if
#endif

                ! Store T
                i = 3 * (my_i_atom - 1)
                j = 3 * (my_j_atom - 1)
                associate (T_sub => dipmat%val(i + 1:i + 3, j + 1:j + 3))
                    T_sub = T_sub + Tij
                end associate

                ! Store gradients
                if (.not. present(grad)) cycle
                if (grad%dcoords .and. i_atom /= j_atom) then
                    associate (dTdR_sub => ddipmat%dr(i + 1:i + 3, j + 1:j + 3, :))
                        dTdR_sub = dTdR_sub + dTij%dr
                    end associate
                end if
                if (grad%dlattice) then
                     do i_latt = 1, 3
                        associate ( &
                            dTda_sub => ddipmat%dlattice(i + 1:i + 3, j + 1:j + 3, i_latt, :) &
                        )
                            dTda_sub = dTda_sub - dTij%dr * (n(i_latt))
                        end associate
                    end do
                end if
                if (grad%dr_vdw) then
                    associate (dTdRvdw_sub => ddipmat%dvdw(i + 1:i + 3, j + 1:j + 3))
                        dTdRvdw_sub = dTdRvdw_sub + dTij%dvdw
                    end associate
                end if
                if (grad%dsigma) then
                    associate (dTdsigma_sub => ddipmat%dsigma(i + 1:i + 3, j + 1:j + 3))
                        dTdsigma_sub = dTdsigma_sub + dTij%dsigma
                    end associate
                end if
#ifdef DO_COMPLEX_TYPE
                if (grad%dq) then
                    associate (dTdq_sub => ddipmat%dq(i + 1:i + 3, j + 1:j + 3, :))
                        dTdq_sub = dTdq_sub + dTij%dq
                    end associate
                end if
#endif
            end do each_atom_pair
        end do each_atom
    end do each_cell

    call geom%clock(-11)

    ! Add long-range Ewald contributions
    if (do_ewald) then
#ifndef DO_COMPLEX_TYPE
        call add_ewald_dipole_parts_real(geom, dipmat, ddipmat, grad)
#else
        call add_ewald_dipole_parts_complex(geom, dipmat, ddipmat, grad, q)
#endif
    end if
end function

#ifndef DO_COMPLEX_TYPE
subroutine add_ewald_dipole_parts_real(geom, dipmat, ddipmat, grad)
    type(matrix_re_t), intent(inout) :: dipmat
    type(grad_matrix_re_t), intent(inout), optional :: ddipmat
#else
subroutine add_ewald_dipole_parts_complex(geom, dipmat, ddipmat, grad, q)
    type(matrix_cplx_t), intent(inout) :: dipmat
    type(grad_matrix_cplx_t), intent(inout), optional :: ddipmat
#endif
    type(geom_t), intent(inout) :: geom
    type(grad_request_t), intent(in), optional :: grad
#ifdef DO_COMPLEX_TYPE
    real(dp), intent(in) :: q(3)
#endif

    logical :: do_surface
    real(dp) :: rec_latt(3, 3), volume, G(3), Rij(3), k(3), &
        k_sq, k_sq_inv, G_Rij, latt_inv(3, 3), &
        dGdA(3), dk_sqdA, dkk_dA(3, 3), vol_prefactor, sin_GR, cos_GR, &
        k_otimes_k(3, 3), exp_k_sq_gamma
    integer :: &
        i_atom, j_atom, i, j, i_xyz, m(3), i_m, &
        range_m(3), my_i_atom, my_j_atom, i_latt, a, b
#ifndef DO_COMPLEX_TYPE
    real(dp) :: exp_GR, iexp_GR
#else
    complex(dp) :: exp_GR, iexp_GR
    integer :: c
    real(dp) :: dkk_dq(3, 3, 3)
#endif

    latt_inv = inverse(geom%lattice)
    rec_latt = 2 * pi * transpose(latt_inv)
    volume = abs(dble(product(eigvals(geom%lattice))))
    vol_prefactor = 4 * pi / volume
    range_m = supercell_circum(rec_latt, geom%rec_space_cutoff)
    call geom%clock(12)
    m = [0, 0, -1]
    each_recip_vec: do i_m = 1, product(1 + 2 * range_m)
        call shift_idx(m, -range_m, range_m)
        G = matmul(rec_latt, m)
#ifdef DO_COMPLEX_TYPE
        k = G + q
#else
        k = G
#endif
        k_sq = sum(k**2)
        k_sq_inv = 1 / k_sq
        if (sqrt(k_sq) > geom%rec_space_cutoff .or. sqrt(k_sq) < 1d-15) cycle

        exp_k_sq_gamma = exp(-k_sq / (4 * geom%gamm**2))
        do concurrent(a=1:3, b=1:3)
            k_otimes_k(a, b) = k(a) * k(b) * k_sq_inv
        end do
        each_atom: do my_i_atom = 1, size(dipmat%idx%i_atom)
            i_atom = dipmat%idx%i_atom(my_i_atom)
            each_atom_pair: do my_j_atom = 1, size(dipmat%idx%j_atom)
                j_atom = dipmat%idx%j_atom(my_j_atom)
                Rij = geom%coords(:, i_atom) - geom%coords(:, j_atom)
                G_Rij = dot_product(G, Rij)
                sin_GR = sin(G_Rij)
                cos_GR = cos(G_Rij)
#ifdef DO_COMPLEX_TYPE
                exp_GR = complex(cos_GR, sin_GR)
                iexp_GR = complex(-sin_GR, cos_GR)
#else
                exp_GR = cos_GR
                iexp_GR = -sin_GR
#endif
                i = 3 * (my_i_atom - 1)
                j = 3 * (my_j_atom - 1)
                associate (T_sub => dipmat%val(i + 1:i + 3, j + 1:j + 3))
                    T_sub = T_sub + (vol_prefactor * exp_k_sq_gamma * exp_GR) * k_otimes_k
                end associate
                if (.not. present(grad)) cycle
                if (grad%dcoords .and. i_atom /= j_atom) then
                    associate (dTdR_sub => ddipmat%dr(i + 1:i + 3, j + 1:j + 3, :))
                        ! TODO should be do-concurrent, but this crashes IBM XL
                        ! 16.1.1, see issue #16
                        do i_xyz = 1, 3
                            dTdR_sub(:, :, i_xyz) = dTdR_sub(:, :, i_xyz) &
                                + (vol_prefactor * exp_k_sq_gamma * G(i_xyz) * iexp_GR) * k_otimes_k
                        end do
                    end associate
                end if
                if (grad%dlattice) then
                    do i_latt = 1, 3
                        do i_xyz = 1, 3
                            dGdA = -latt_inv(i_latt, :) * G(i_xyz)
                            dk_sqdA = dot_product(k, dGdA)
                            do concurrent(a=1:3, b=1:3)
                                dkk_dA(a, b) = k(a) * k_sq_inv * dGdA(b)
                            end do
                            dkk_dA = dkk_dA + transpose(dkk_dA)
                            ! Using associate here was causing weird seg faults
                            ! with some Intel compilers, reporting i_xyz being
                            ! zero in the index, even though it printed as 1
                            ! associate ( &
                            !     dTda_sub => ddipmat%dlattice(i+1:i+3, j+1:j+3, i_latt, i_xyz) &
                            ! )
                            ddipmat%dlattice(i + 1:i + 3, j + 1:j + 3, i_latt, i_xyz) = &
                                ddipmat%dlattice(i + 1:i + 3, j + 1:j + 3, i_latt, i_xyz) &
                                + vol_prefactor * exp_k_sq_gamma * ( &
                                    - exp_GR * latt_inv(i_latt, i_xyz) &
                                    - exp_GR * dk_sqdA * (1 + 1 / (2 * geom%gamm**2)) &
                                    + iexp_GR * dot_product(dGdA, Rij) &
                                ) * k_otimes_k &
                                + vol_prefactor * exp_k_sq_gamma * exp_GR * dkk_dA
                            ! end associate
                        end do
                    end do
                end if
#ifdef DO_COMPLEX_TYPE
                if (grad%dq) then
                    do concurrent(a=1:3, b=1:3, c=1:3)
                        dkk_dq(a, b, c) = -2 * k(a) * k(b) * k(c) * k_sq_inv**2
                    end do
                    do concurrent(a=1:3, b=1:3)
                        dkk_dq(b, a, a) = dkk_dq(b, a, a) + k(b) * k_sq_inv
                    end do
                    do concurrent(a=1:3, b=1:3)
                        dkk_dq(a, b, a) = dkk_dq(a, b, a) + k(b) * k_sq_inv
                    end do
                    associate (dTdq_sub => ddipmat%dq(i + 1:i + 3, j + 1:j + 3, :))
                        dTdq_sub = dTdq_sub + vol_prefactor * exp_k_sq_gamma * exp_GR * dkk_dq
                        ! TODO should be do-concurrent, but this crashes IBM XL
                        ! 16.1.1, see issue #16
                        do a = 1, 3
                            dTdq_sub(:, :, a) = dTdq_sub(:, :, a) &
                                - (1 / (2 * geom%gamm**2) * vol_prefactor * exp_k_sq_gamma * exp_GR) * k(a) * k_otimes_k
                        end do
                    end associate
                end if
#endif
            end do each_atom_pair
        end do each_atom
    end do each_recip_vec

    ! self energy
    call dipmat%add_diag_scalar(-4 * geom%gamm**3 / (3 * sqrt(pi)))

    ! surface term
#ifdef DO_COMPLEX_TYPE
    do_surface = sqrt(sum(q**2)) < 1d-15
#else
    do_surface = .true.
#endif
    if (do_surface) then
        do my_i_atom = 1, size(dipmat%idx%i_atom)
            do my_j_atom = 1, size(dipmat%idx%j_atom)
                do i_xyz = 1, 3
                    i = 3 * (my_i_atom - 1) + i_xyz
                    j = 3 * (my_j_atom - 1) + i_xyz
                    dipmat%val(i, j) = dipmat%val(i, j) + vol_prefactor / 3
                    if (.not. present(grad)) cycle
                    if (grad%dlattice) then
                        ddipmat%dlattice(i, j, :, :) = ddipmat%dlattice(i, j, :, :) &
                            - vol_prefactor / 3 * latt_inv
                    end if
                end do
            end do
        end do
    end if
    call geom%clock(-12)
end subroutine

#ifndef DO_COMPLEX_TYPE
#   define DO_COMPLEX_TYPE
#   include "mbd_dipole.F90"

function build_T(r, d_inv, B, C, D, request_dr, request_dparam, dr, dB, dC, dparam) result(T)
    !! $$
    !! T_{ab}(\mathbf r)=\frac{\partial^2}{\partial r_a\partial r_b}\frac1r=
    !! \frac{-3r_ar_b+r^2\delta_{ab}}{r^5},\qquad
    !! \frac{\partial T_{ab}(\mathbf r)}{\partial r_c}=-3\left(
    !! \frac{r_a\delta_{bc}+r_b\delta_{ca}+r_c\delta_{ab}}{r^5}-
    !! \frac{5r_ar_br_c}{r^7}
    !! \right)
    !! $$
    real(dp), intent(in) :: r(3), d_inv, B, C, D
    logical, intent(in) :: request_dr, request_dparam
    real(dp), intent(in), optional :: dB, dC
    real(dp), intent(out), optional :: dr(3, 3, 3), dparam(3, 3)
    real(dp) :: T(3, 3)

    real(dp) :: C_over_d
    real(dp) :: ux, uy, uz
    real(dp) :: xx, yy, zz, xy, xz, yz
    real(dp) :: xxx, yyy, zzz, yxx, yyx, zxx, zyy, zzx, zzy, xyz

    ux = d_inv * r(1)
    uy = d_inv * r(2)
    uz = d_inv * r(3)

    xx = C * ux**2 + B
    yy = C * uy**2 + B
    zz = C * uz**2 + B
    xy = C * ux * uy
    xz = C * ux * uz
    yz = C * uy * uz

    T(1, 1) = xx
    T(1, 2) = xy
    T(1, 3) = xz
    T(2, 1) = xy
    T(2, 2) = yy
    T(2, 3) = yz
    T(3, 1) = xz
    T(3, 2) = yz
    T(3, 3) = zz

    if (request_dr .and. present(dr)) then
        C_over_d = C * d_inv

        ! Diagonal
        xxx = ux * (3 * C_over_d + ux**2 * D)
        yyy = uy * (3 * C_over_d + uy**2 * D)
        zzz = uz * (3 * C_over_d + uz**2 * D)

        ! Terms with two identical indices
        yxx = uy * (C_over_d + ux**2 * D)
        yyx = ux * (C_over_d + uy**2 * D)
        zxx = uz * (C_over_d + ux**2 * D)
        zyy = uz * (C_over_d + uy**2 * D)
        zzx = ux * (C_over_d + uz**2 * D)
        zzy = uy * (C_over_d + uz**2 * D)

        ! Terms with all indices distinct
        xyz = ux * uy * uz * D

        dr(1, 1, 1) = xxx
        dr(2, 1, 1) = yxx
        dr(3, 1, 1) = zxx
        dr(1, 2, 1) = yxx
        dr(2, 2, 1) = yyx
        dr(3, 2, 1) = xyz
        dr(1, 3, 1) = zxx
        dr(2, 3, 1) = xyz
        dr(3, 3, 1) = zzx

        dr(1, 1, 2) = yxx
        dr(2, 1, 2) = yyx
        dr(3, 1, 2) = xyz
        dr(1, 2, 2) = yyx
        dr(2, 2, 2) = yyy
        dr(3, 2, 2) = zyy
        dr(1, 3, 2) = xyz
        dr(2, 3, 2) = zyy
        dr(3, 3, 2) = zzy

        dr(1, 1, 3) = zxx
        dr(2, 1, 3) = xyz
        dr(3, 1, 3) = zzx
        dr(1, 2, 3) = xyz
        dr(2, 2, 3) = zyy
        dr(3, 2, 3) = zzy
        dr(1, 3, 3) = zzx
        dr(2, 3, 3) = zzy
        dr(3, 3, 3) = zzz
    end if

    if (request_dparam .and. present(dparam) .and. present(dB) .and. present(dC)) then
        xx = dC * ux**2 + dB
        yy = dC * uy**2 + dB
        zz = dC * uz**2 + dB
        xy = dC * ux * uy
        xz = dC * ux * uz
        yz = dC * uy * uz

        dparam(1, 1) = xx
        dparam(2, 1) = xy
        dparam(3, 1) = xz
        dparam(1, 2) = xy
        dparam(2, 2) = yy
        dparam(3, 2) = yz
        dparam(1, 3) = xz
        dparam(2, 3) = yz
        dparam(3, 3) = zz
    end if
end function

function T_bare(r, dT, grad) result(T)
    !! $$
    !! T_{ab}(\mathbf r)=\frac{\partial^2}{\partial r_a\partial r_b}\frac1r=
    !! \frac{-3r_ar_b+r^2\delta_{ab}}{r^5},\qquad
    !! \frac{\partial T_{ab}(\mathbf r)}{\partial r_c}=-3\left(
    !! \frac{r_a\delta_{bc}+r_b\delta_{ca}+r_c\delta_{ab}}{r^5}-
    !! \frac{5r_ar_br_c}{r^7}
    !! \right)
    !! $$
    real(dp), intent(in) :: r(3)
    type(grad_matrix_re_t), intent(inout), optional :: dT
    logical, intent(in), optional :: grad
    real(dp) :: T(3, 3)

    real(dp) :: d_inv, B, C, D

    d_inv = 1_dp / sqrt(sum(r**2))

    B = d_inv**3
    C = -3 * B
    D = 15 * d_inv * B

    T = build_T(r, d_inv, B, C, D, grad, .false., dT%dr)
end function

function T_erfc(r, gamm, dT, grad) result(T)
    !! $$
    !! T_{ab}^\text{erfc}(\mathbf r,\gamma)
    !! =-3\frac{r_ar_b}{r^5}C(r,\gamma)+\frac{\delta_{ab}}{r^3}B(r,\gamma)
    !! $$
    !!
    !! $$
    !! \begin{aligned}
    !! \frac{\partial T_{ab}^\text{erfc}(\mathbf r,\gamma)}{\partial r_c}
    !! &=-\left(
    !! \frac{r_a\delta_{bc}+r_b\delta_{ca}}{r^5}-
    !! 5\frac{r_ar_br_c}{r^7}
    !! \right)C(r,\gamma)-3\frac{r_c\delta_{ab}}{r^5}B(r,\gamma)
    !! \\ &-\frac{r_ar_br_c}{r^6}\frac{\partial C(r,\gamma)}{\partial
    !! r}+\frac{r_c\delta_{ab}}{r^4}\frac{\partial B(r,\gamma)}{\partial r}
    !! \end{aligned}
    !! $$
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: gamm
    type(grad_matrix_re_t), intent(inout), optional :: dT
    type(grad_request_t), intent(in), optional :: grad
    real(dp) :: T(3, 3)

    real(dp) :: r_1, r_2, d_inv
    real(dp) :: x, theta, rho, B, C, D, dB, dC

    r_2 = sum(r**2)
    r_1 = sqrt(r_2)
    d_inv = 1_dp / r_1

    x = gamm * r_1
    theta = (2 / sqrt(pi)) * x * exp(-x**2)
    B = d_inv**3 * (erfc(x) + theta)
    rho = 2 * gamm**2 * theta
    C = -3 * B - d_inv * rho
    D = (2 * gamm**2 * rho - 5 * d_inv * C)

    dB = 0_dp
    dC = 0_dp
    if (grad%dgamma) then
        dB = -2 * d_inv * gamm * theta
        dC = -2 * x**2 * dB
    end if

    T = build_T(r, d_inv, B, C, D, grad%dcoords, grad%dgamma, dT%dr, dB, dC, dT%dgamma)
end function

! FIXME: Add unit test to compare against T_erfc - T_bare
function T_erfc_minus_bare(r, gamm, dT, grad) result(T)
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: gamm
    type(grad_matrix_re_t), intent(inout), optional :: dT
    type(grad_request_t), intent(in), optional :: grad
    real(dp) :: T(3, 3)

    real(dp) :: r_1, r_2, d_inv
    real(dp) :: x, theta, rho, B, C, D, dB, dC

    r_2 = sum(r**2)
    r_1 = sqrt(r_2)
    d_inv = 1_dp / r_1

    x = gamm * r_1
    theta = (2 / sqrt(pi)) * x * exp(-x**2)
    B = d_inv**3 * (erfc(x) + theta - 1)
    rho = 2 * gamm**2 * theta
    C = -3 * B - d_inv * rho
    D = (2 * gamm**2 * rho - 5 * d_inv * C)

    dB = 0_dp
    dC = 0_dp
    if (grad%dgamma) then
        dB = -2 * d_inv * gamm * theta
        dC = -2 * x**2 * dB
    end if

    T = build_T(r, d_inv, B, C, D, grad%dcoords, grad%dgamma, dT%dr, dB, dC, dT%dgamma)
end function

function T_erf_coulomb(r, sigma, dT, grad) result(T)
    !! $$
    !! \begin{aligned}
    !! T^\text{GG}_{ab}(\mathbf r,\sigma)&=
    !! \frac{\partial^2}{\partial r_a\partial r_b}\frac{\operatorname{erf}(\zeta)}r=
    !! \big(\operatorname{erf}(\zeta)-\Theta(\zeta)\big)T_{ab}(\mathbf r)+
    !! 2\zeta^2\Theta(\zeta)\frac{r_ar_b}{r^5}
    !! \\ \Theta(\zeta)&=\frac{2\zeta}{\sqrt\pi}\exp(-\zeta^2),\qquad
    !! \zeta=\frac r\sigma
    !! \\ \frac{\mathrm d T_{ab}^\text{GG}(\mathbf r,\sigma)}{\mathrm dr_c}&=
    !! 2\zeta\Theta(\zeta)\left(T_{ab}(\mathbf r)+(3-2\zeta^2)\frac{r_ar_b}{r^5}\right)
    !! \frac{\mathrm d\zeta}{\mathrm dr_c}
    !! \\ &+\big(\operatorname{erf}(\zeta)-\Theta(\zeta)\big)
    !! \frac{\partial T_{ab}(\mathbf r)}{\partial r_c}-
    !! 2\zeta^2\Theta(\zeta)\left(
    !! \frac13\frac{\partial T_{ab}(\mathbf r)}{\partial r_c}+
    !! \frac{r_c\delta_{ab}}{r^5}
    !! \right)
    !! \\ \qquad\frac{\mathrm d\zeta}{\mathrm dr_c}&=
    !! \frac{r_c}{r\sigma}-\frac r{\sigma^2}\frac{\mathrm d\sigma}{\mathrm dr_c}
    !! \end{aligned}
    !! $$
    real(dp), intent(in) :: r(3)
    real(dp), intent(in) :: sigma
    type(grad_matrix_re_t), intent(inout), optional :: dT
    type(grad_request_t), intent(in), optional :: grad
    real(dp) :: T(3, 3)

    real(dp) :: r_1, d_inv, sigma_inv, x, theta, rho, B, C, D, dB, dC

    r_1 = sqrt(sum(r**2))
    d_inv = 1_dp / r_1
    sigma_inv = 1 / sigma

    x = r_1 * sigma_inv
    theta = (2 / sqrt(pi)) * x * exp(-x**2)
    B = d_inv**3 * (erf(x) - theta)
    rho = 2 * sigma_inv**2 * theta
    C = -3 * B + d_inv * rho
    D = -2 * sigma_inv**2 * rho - 5 * d_inv * C

    dB = 0_dp
    dC = 0_dp
    if (grad%dsigma) then
        dB = -2 * d_inv * sigma_inv**3 * theta
        dC = -2 * x**2 * dB
    end if

    T = build_T(r, d_inv, B, C, D, grad%dcoords, grad%dsigma, dT%dr, dB, dC, dT%dsigma)
end function

function T_1mexp_coulomb(r, sigma, a, dT, grad) result(T)
    real(dp), intent(in) :: r(3), sigma, a
    type(grad_matrix_re_t), intent(inout), optional :: dT
    type(grad_request_t), intent(in), optional :: grad
    real(dp) :: T(3, 3)

    real(dp) :: r_1, d_inv, x, B, C, D, tmp, dB, dC

    r_1 = sqrt(sum(r**2))
    d_inv = 1 / r_1

    x = (r_1 / sigma)**a
    B = -x * a * exp(-x) * d_inv**2
    C = -B * (2 + a * (x - 1_dp))
    D = B * d_inv * (8 + a * (6 * (x - 1) + a * (x**2 - 3 * x + 1)))

    dB = 0_dp
    dC = 0_dp
    if (grad%dsigma) then
        tmp = B * a * d_inv * x**(1 / a)
        dB = tmp * (x - 1)
        dC = -tmp * (2 * (x - 1) + a * (x**2 - 3 * x + 1))
    end if

    T = build_T(r, d_inv, B, C, D, grad%dcoords, grad%dsigma, dT%dr, dB, dC, dT%dsigma)
end function

function damping_grad(f, df, T, dT, dfT, grad) result(fT)
    real(dp), intent(in) :: f
    type(grad_scalar_t), intent(in) :: df
    real(dp), intent(in) :: T(3, 3)
    type(grad_matrix_re_t), intent(inout) :: dT
    type(grad_matrix_re_t), intent(inout) :: dfT
    type(grad_request_t), intent(in) :: grad
    real(dp) :: fT(3, 3)

    integer :: c, b, a

    do concurrent(c = 1:3, b = 1:3)
        fT(b, c) = f * T(b, c)
    end do
    if (grad%dcoords) then
        dfT%dr = 0_dp
        if (allocated(df%dr)) then
            do concurrent(c = 1:3, b = 1:3, a = 1:3)
                dfT%dr(a, b, c) = df%dr(c) * T(a, b)
            end do
        end if
        if (allocated(dT%dr)) then
            do concurrent(c = 1:3, b = 1:3, a = 1:3)
                dfT%dr(a, b, c) = dfT%dr(a, b, c) + f * dT%dr(a, b, c)
            end do
        end if
    end if
    if (grad%dr_vdw) then
        dfT%dvdw = 0_dp
        if (allocated(df%dvdw)) dfT%dvdw = df%dvdw * T
        if (allocated(dT%dvdw)) dfT%dvdw = dfT%dvdw + f * dT%dvdw
    end if
    if (grad%dsigma) dfT%dsigma = f * dT%dsigma
end function

end module

#endif
