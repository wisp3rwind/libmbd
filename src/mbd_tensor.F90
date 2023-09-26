#ifndef DO_COMPLEX_TYPE

module mbd_tensor

use mbd_constants

implicit none

private

type, public :: tensor_3x3_re_grad_scalar_t
    real(dp) :: xx, yy, zz, xy, xz, yz
contains
    procedure :: copy_to => tensor_3x3_re_grad_scalar_copy_to
    procedure :: add_to => tensor_3x3_re_grad_scalar_add_to
end type

type, public :: tensor_3x3_cplx_grad_scalar_t
    complex(dp) :: xx, yy, zz, xy, xz, yz
contains
    procedure :: copy_to => tensor_3x3_cplx_grad_scalar_copy_to
    procedure :: add_to => tensor_3x3_cplx_grad_scalar_add_to
end type

type, public :: tensor_3x3_re_grad_vector_t
    real(dp) :: xx_dx, xx_dy, xx_dz
    real(dp) :: yy_dx, yy_dy, yy_dz
    real(dp) :: zz_dx, zz_dy, zz_dz
    real(dp) :: xy_dx, xy_dy, xy_dz
    real(dp) :: xz_dx, xz_dy, xz_dz
    real(dp) :: yz_dx, yz_dy, yz_dz
contains
    procedure :: copy_to => tensor_3x3_re_grad_vector_copy_to
    procedure :: add_to => tensor_3x3_re_grad_vector_add_to
end type

type, public :: tensor_3x3_cplx_grad_vector_t
    complex(dp) :: xx_dx, xx_dy, xx_dz
    complex(dp) :: yy_dx, yy_dy, yy_dz
    complex(dp) :: zz_dx, zz_dy, zz_dz
    complex(dp) :: xy_dx, xy_dy, xy_dz
    complex(dp) :: xz_dx, xz_dy, xz_dz
    complex(dp) :: yz_dx, yz_dy, yz_dz
contains
    procedure :: copy_to => tensor_3x3_cplx_grad_vector_copy_to
    procedure :: add_to => tensor_3x3_cplx_grad_vector_add_to
end type

type, public :: tensor_3x3_re_t
    real(dp) :: xx, yy, zz, xy, xz, yz
contains
    procedure :: copy_to => tensor_3x3_re_copy_to
    procedure :: add_to => tensor_3x3_re_add_to
    procedure, pass(this) :: mul_re => tensor_3x3_re_mul_re
    procedure, pass(this) :: mul_cplx => tensor_3x3_re_mul_cplx
    procedure, pass(this) :: add_re => tensor_3x3_re_add_re
    procedure, pass(this) :: sub_re => tensor_3x3_re_sub_re
    procedure, pass(this) :: assign_scalar_re => tensor_3x3_re_assign_re
    procedure, pass(this) :: clear => tensor_3x3_re_clear

    generic :: operator(*) => mul_re, mul_cplx
    generic :: operator(+) => add_re
    generic :: operator(-) => sub_re
    generic :: assignment(=) => assign_scalar_re
end type

type, public :: tensor_3x3_cplx_t
    complex(dp) :: xx, yy, zz, xy, xz, yz
contains
    procedure :: copy_to => tensor_3x3_cplx_copy_to
    procedure :: add_to => tensor_3x3_cplx_add_to
    procedure, pass(this) :: mul_re => tensor_3x3_cplx_mul_re
    procedure, pass(this) :: mul_cplx => tensor_3x3_cplx_mul_cplx
    procedure, pass(this) :: add_cplx => tensor_3x3_cplx_add_cplx
    procedure, pass(this) :: sub_cplx => tensor_3x3_cplx_sub_cplx
    procedure, pass(this) :: assign_scalar_cplx => tensor_3x3_cplx_assign_cplx
    procedure, pass(this) :: clear => tensor_3x3_cplx_clear

    generic :: operator(*) => mul_re, mul_cplx
    generic :: operator(+) => add_cplx
    generic :: operator(-) => sub_cplx
    generic :: assignment(=) => assign_scalar_cplx
end type

contains

#endif

#ifndef DO_COMPLEX_TYPE
pure subroutine tensor_3x3_re_copy_to(this, dest)
    class(tensor_3x3_re_t), intent(in) :: this
    real(dp), intent(inout) :: dest(:, :)
#else
pure subroutine tensor_3x3_cplx_copy_to(this, dest)
    class(tensor_3x3_cplx_t), intent(in) :: this
    complex(dp), intent(inout) :: dest(:, :)
#endif

    dest(1, 1) = this%xx
    dest(1, 2) = this%xy
    dest(1, 3) = this%xz
    dest(2, 1) = this%xy
    dest(2, 2) = this%yy
    dest(2, 3) = this%yz
    dest(3, 1) = this%xz
    dest(3, 2) = this%yz
    dest(3, 3) = this%zz
end subroutine

#ifndef DO_COMPLEX_TYPE
pure subroutine tensor_3x3_re_add_to(this, dest)
    class(tensor_3x3_re_t), intent(in) :: this
    real(dp), intent(inout) :: dest(:, :)
#else
pure subroutine tensor_3x3_cplx_add_to(this, dest)
    class(tensor_3x3_cplx_t), intent(in) :: this
    complex(dp), intent(inout) :: dest(:, :)
#endif

    dest(1, 1) = dest(1, 1) + this%xx
    dest(1, 2) = dest(1, 2) + this%xy
    dest(1, 3) = dest(1, 3) + this%xz
    dest(2, 1) = dest(2, 1) + this%xy
    dest(2, 2) = dest(2, 2) + this%yy
    dest(2, 3) = dest(2, 3) + this%yz
    dest(3, 1) = dest(3, 1) + this%xz
    dest(3, 2) = dest(3, 2) + this%yz
    dest(3, 3) = dest(3, 3) + this%zz
end subroutine

#ifndef DO_COMPLEX_TYPE
pure function tensor_3x3_re_add_re(this, other) result(res)
    class(tensor_3x3_re_t), intent(in) :: this, other
    type(tensor_3x3_re_t) :: res
#else
pure function tensor_3x3_cplx_add_cplx(this, other) result(res)
    class(tensor_3x3_cplx_t), intent(in) :: this, other
    type(tensor_3x3_cplx_t) :: res
#endif

    res%xx = this%xx + other%xx
    res%yy = this%yy + other%yy
    res%zz = this%zz + other%zz
    res%xy = this%xy + other%xy
    res%xz = this%xz + other%xz
    res%yz = this%yz + other%yz
end function

#ifndef DO_COMPLEX_TYPE
pure function tensor_3x3_re_sub_re(this, other) result(res)
    class(tensor_3x3_re_t), intent(in) :: this, other
    type(tensor_3x3_re_t) :: res
#else
pure function tensor_3x3_cplx_sub_cplx(this, other) result(res)
    class(tensor_3x3_cplx_t), intent(in) :: this, other
    type(tensor_3x3_cplx_t) :: res
#endif

    res%xx = this%xx - other%xx
    res%yy = this%yy - other%yy
    res%zz = this%zz - other%zz
    res%xy = this%xy - other%xy
    res%xz = this%xz - other%xz
    res%yz = this%yz - other%yz
end function

#ifndef DO_COMPLEX_TYPE
pure subroutine tensor_3x3_re_assign_re(this, c)
    class(tensor_3x3_re_t), intent(inout) :: this
    real(dp), intent(in) :: c
#else
pure subroutine tensor_3x3_cplx_assign_cplx(this, c)
    class(tensor_3x3_cplx_t), intent(inout) :: this
    complex(dp), intent(in) :: c
#endif

    this%xx = c
    this%yy = c
    this%zz = c
    this%xy = c
    this%xz = c
    this%yz = c
end subroutine

#ifndef DO_COMPLEX_TYPE
pure subroutine tensor_3x3_re_clear(this)
    class(tensor_3x3_re_t), intent(inout) :: this
#else
pure subroutine tensor_3x3_cplx_clear(this)
    class(tensor_3x3_cplx_t), intent(inout) :: this
#endif

    this%xx = 0.0
    this%yy = 0.0
    this%zz = 0.0
    this%xy = 0.0
    this%xz = 0.0
    this%yz = 0.0
end subroutine

#ifndef DO_COMPLEX_TYPE
pure function tensor_3x3_re_mul_re(c, this) result(res)
    real(dp), intent(in) :: c
    class(tensor_3x3_re_t), intent(in) :: this
    type(tensor_3x3_re_t) :: res
#else
pure function tensor_3x3_cplx_mul_re(c, this) result(res)
    real(dp), intent(in) :: c
    class(tensor_3x3_cplx_t), intent(in) :: this
    type(tensor_3x3_cplx_t) :: res
#endif

    res%xx = c * this%xx
    res%yy = c * this%yy
    res%zz = c * this%zz
    res%xy = c * this%xy
    res%xz = c * this%xz
    res%yz = c * this%yz
end function

#ifndef DO_COMPLEX_TYPE
pure function tensor_3x3_re_mul_cplx(c, this) result(res)
    complex(dp), intent(in) :: c
    class(tensor_3x3_re_t), intent(in) :: this
    type(tensor_3x3_cplx_t) :: res
#else
pure function tensor_3x3_cplx_mul_cplx(c, this) result(res)
    complex(dp), intent(in) :: c
    class(tensor_3x3_cplx_t), intent(in) :: this
    type(tensor_3x3_cplx_t) :: res
#endif

    res%xx = c * this%xx
    res%yy = c * this%yy
    res%zz = c * this%zz
    res%xy = c * this%xy
    res%xz = c * this%xz
    res%yz = c * this%yz
end function

#ifndef DO_COMPLEX_TYPE
pure subroutine tensor_3x3_re_grad_scalar_copy_to(this, dest)
    class(tensor_3x3_re_grad_scalar_t), intent(inout) :: this
    real(dp), intent(inout) :: dest(:, :)
#else
pure subroutine tensor_3x3_cplx_grad_scalar_copy_to(this, dest)
    class(tensor_3x3_cplx_grad_scalar_t), intent(inout) :: this
    complex(dp), intent(inout) :: dest(:, :)
#endif

    dest(1, 1) = this%xx
    dest(1, 2) = this%xy
    dest(1, 3) = this%xz
    dest(2, 1) = this%xy
    dest(2, 2) = this%yy
    dest(2, 3) = this%yz
    dest(3, 1) = this%xz
    dest(3, 2) = this%yz
    dest(3, 3) = this%zz
end subroutine

#ifndef DO_COMPLEX_TYPE
pure subroutine tensor_3x3_re_grad_scalar_add_to(this, dest)
    class(tensor_3x3_re_grad_scalar_t), intent(in) :: this
    real(dp), intent(inout) :: dest(:, :)
#else
pure subroutine tensor_3x3_cplx_grad_scalar_add_to(this, dest)
    class(tensor_3x3_cplx_grad_scalar_t), intent(in) :: this
    complex(dp), intent(inout) :: dest(:, :)
#endif

    dest(1, 1) = dest(1, 1) + this%xx
    dest(1, 2) = dest(1, 2) + this%xy
    dest(1, 3) = dest(1, 3) + this%xz
    dest(2, 1) = dest(2, 1) + this%xy
    dest(2, 2) = dest(2, 2) + this%yy
    dest(2, 3) = dest(2, 3) + this%yz
    dest(3, 1) = dest(3, 1) + this%xz
    dest(3, 2) = dest(3, 2) + this%yz
    dest(3, 3) = dest(3, 3) + this%zz
end subroutine

#ifndef DO_COMPLEX_TYPE
pure subroutine tensor_3x3_re_grad_vector_copy_to(this, dest)
    class(tensor_3x3_re_grad_vector_t), intent(inout) :: this
    real(dp), intent(inout) :: dest(:, :, :)
#else
pure subroutine tensor_3x3_cplx_grad_vector_copy_to(this, dest)
    class(tensor_3x3_cplx_grad_vector_t), intent(inout) :: this
    complex(dp), intent(inout) :: dest(:, :, :)
#endif

    dest(1, 1, 1) = this%xx_dx
    dest(2, 1, 1) = this%xy_dx
    dest(3, 1, 1) = this%xz_dx
    dest(1, 2, 1) = this%xy_dx
    dest(2, 2, 1) = this%yy_dx
    dest(3, 2, 1) = this%yz_dx
    dest(1, 3, 1) = this%xz_dx
    dest(2, 3, 1) = this%yz_dx
    dest(3, 3, 1) = this%zz_dx

    dest(1, 1, 2) = this%xx_dy
    dest(2, 1, 2) = this%xy_dy
    dest(3, 1, 2) = this%xz_dy
    dest(1, 2, 2) = this%xy_dy
    dest(2, 2, 2) = this%yy_dy
    dest(3, 2, 2) = this%yz_dy
    dest(1, 3, 2) = this%xz_dy
    dest(2, 3, 2) = this%yz_dy
    dest(3, 3, 2) = this%zz_dy

    dest(1, 1, 3) = this%xx_dz
    dest(2, 1, 3) = this%xy_dz
    dest(3, 1, 3) = this%xz_dz
    dest(1, 2, 3) = this%xy_dz
    dest(2, 2, 3) = this%yy_dz
    dest(3, 2, 3) = this%yz_dz
    dest(1, 3, 3) = this%xz_dz
    dest(2, 3, 3) = this%yz_dz
    dest(3, 3, 3) = this%zz_dz
end subroutine

#ifndef DO_COMPLEX_TYPE
pure subroutine tensor_3x3_re_grad_vector_add_to(this, dest)
    class(tensor_3x3_re_grad_vector_t), intent(in) :: this
    real(dp), intent(inout) :: dest(:, :, :)
#else
pure subroutine tensor_3x3_cplx_grad_vector_add_to(this, dest)
    class(tensor_3x3_cplx_grad_vector_t), intent(in) :: this
    complex(dp), intent(inout) :: dest(:, :, :)
#endif

    dest(1, 1, 1) = dest(1, 1, 1) + this%xx_dx
    dest(2, 1, 1) = dest(2, 1, 1) + this%xy_dx
    dest(3, 1, 1) = dest(3, 1, 1) + this%xz_dx
    dest(1, 2, 1) = dest(1, 2, 1) + this%xy_dx
    dest(2, 2, 1) = dest(2, 2, 1) + this%yy_dx
    dest(3, 2, 1) = dest(3, 2, 1) + this%yz_dx
    dest(1, 3, 1) = dest(1, 3, 1) + this%xz_dx
    dest(2, 3, 1) = dest(2, 3, 1) + this%yz_dx
    dest(3, 3, 1) = dest(3, 3, 1) + this%zz_dx

    dest(1, 1, 2) = dest(1, 1, 2) + this%xx_dy
    dest(2, 1, 2) = dest(2, 1, 2) + this%xy_dy
    dest(3, 1, 2) = dest(3, 1, 2) + this%xz_dy
    dest(1, 2, 2) = dest(1, 2, 2) + this%xy_dy
    dest(2, 2, 2) = dest(2, 2, 2) + this%yy_dy
    dest(3, 2, 2) = dest(3, 2, 2) + this%yz_dy
    dest(1, 3, 2) = dest(1, 3, 2) + this%xz_dy
    dest(2, 3, 2) = dest(2, 3, 2) + this%yz_dy
    dest(3, 3, 2) = dest(3, 3, 2) + this%zz_dy

    dest(1, 1, 3) = dest(1, 1, 3) + this%xx_dz
    dest(2, 1, 3) = dest(2, 1, 3) + this%xy_dz
    dest(3, 1, 3) = dest(3, 1, 3) + this%xz_dz
    dest(1, 2, 3) = dest(1, 2, 3) + this%xy_dz
    dest(2, 2, 3) = dest(2, 2, 3) + this%yy_dz
    dest(3, 2, 3) = dest(3, 2, 3) + this%yz_dz
    dest(1, 3, 3) = dest(1, 3, 3) + this%xz_dz
    dest(2, 3, 3) = dest(2, 3, 3) + this%yz_dz
    dest(3, 3, 3) = dest(3, 3, 3) + this%zz_dz
end subroutine

#ifndef DO_COMPLEX_TYPE
#   define DO_COMPLEX_TYPE
#   include "mbd_tensor.F90"

end module

#endif
