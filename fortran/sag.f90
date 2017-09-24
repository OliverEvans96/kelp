module sag
use utils
use fastgl

implicit none

! Spatial grids do not include upper endpoints.
! Angular grids do include upper endpoints.
! Both include lower endpoints.

! To use:
! call grid%set_bounds(...)
! call grid%set_num(...) (or set_spacing)
! call grid%init()
! ...
! call grid%deinit()

!integer, parameter :: pi = 3.141592653589793D+00

type index_list
   integer i, j, k, l, m
 contains
   procedure :: init => index_list_init
   procedure :: print => index_list_print
end type index_list

type angle_dim
   integer num
   double precision minval, maxval, prefactor
   double precision, dimension(:), allocatable :: vals, weights, sin, cos
 contains
   procedure :: set_bounds => angle_set_bounds
   procedure :: set_num => angle_set_num
   procedure :: deinit => angle_deinit
   procedure :: integrate_points => angle_integrate_points
   procedure :: integrate_func => angle_integrate_func
   procedure :: assign_legendre
end type angle_dim

type space_dim
   integer num
   double precision minval, maxval, spacing
   double precision, dimension(:), allocatable :: vals
 contains
   procedure :: integrate_points => space_integrate_points
   procedure :: trapezoid_rule
   procedure :: set_bounds => space_set_bounds
   procedure :: set_num => space_set_num
   procedure :: set_spacing => space_set_spacing
   procedure :: set_num_from_spacing
   procedure :: set_spacing_from_num
   procedure :: deinit => space_deinit
   procedure :: assign_linspace
end type space_dim

type space_angle_grid !(sag)
  type(space_dim) :: x, y, z
  type(angle_dim) :: theta, phi
contains
  procedure :: set_bounds => sag_set_bounds
  procedure :: set_num => sag_set_num
  procedure :: integrate_angle_2d => sag_integrate_angle_2d
  procedure :: init => sag_init
  procedure :: deinit => sag_deinit
  procedure :: set_num_from_spacing => sag_set_num_from_spacing
  procedure :: set_spacing_from_num => sag_set_spacing_from_num
end type space_angle_grid

contains

  subroutine index_list_init(indices)
    class(index_list) indices
    indices%i = 1
    indices%j = 1
    indices%k = 1
    indices%l = 1
    indices%m = 1
  end subroutine

  subroutine index_list_print(indices)
    class(index_list) indices

    write(*,*) 'i, j, k, l, m =', indices%i, indices%j, indices%k, indices%l, indices%m
  end subroutine index_list_print

  subroutine angle_set_bounds(angle, minval, maxval)
    class(angle_dim) :: angle
    double precision minval, maxval
    angle%minval = minval
    angle%maxval = maxval
  end subroutine angle_set_bounds

  subroutine angle_set_num(angle, num)
    class(angle_dim) :: angle
    integer num
    angle%num = num
  end subroutine angle_set_num

  ! To calculate \int_{xmin}^{xmax} f(x) dx :
  ! int = prefactor * sum(weights * f(roots))
  subroutine assign_legendre(angle)
    class(angle_dim) :: angle
    double precision root, weight, theta
    integer i
    ! glpair produces both x and theta, where x=cos(theta). We'll throw out theta.

    allocate(angle%vals(angle%num))
    allocate(angle%weights(angle%num))
    allocate(angle%sin(angle%num))
    allocate(angle%cos(angle%num))

    ! Prefactor for integration
    ! From change of variables
    angle%prefactor = (angle%maxval - angle%minval) / 2.d0

    do i = 1, angle%num
       call glpair(angle%num, i, theta, weight, root)
       call affine_transform(root, -1.d0, 1.d0, angle%minval, angle%maxval)
       angle%vals(i) = root
       angle%weights(i) = weight
       angle%sin(i) = sin(root)
       angle%cos(i) = cos(root)
    end do

  end subroutine assign_legendre

  ! Integrate callable function over angle via Gauss-Legendre quadrature

  function angle_integrate_func(angle, func_callable) result(integral)
    class(angle_dim) :: angle
    double precision, external :: func_callable
    double precision, dimension(:), allocatable :: func_vals
    double precision integral
    integer i

    allocate(func_vals(angle%num))

    do i=1, angle%num
       func_vals(i) = func_callable(angle%vals(i))
    end do

    integral = angle%integrate_points(func_vals)

    deallocate(func_vals)
  end function angle_integrate_func

  ! Integrate function given function values sampled at legendre theta values
  function angle_integrate_points(angle, func_vals) result(integral)
    class(angle_dim) :: angle
    double precision, dimension(angle%num) :: func_vals
    double precision integral

    integral = angle%prefactor * sum(angle%weights * func_vals)
  end function angle_integrate_points

  subroutine angle_deinit(angle)
    class(angle_dim) :: angle
    deallocate(angle%vals)
    deallocate(angle%weights)
    deallocate(angle%sin)
    deallocate(angle%cos)
  end subroutine angle_deinit


  !! SPACE !!

  ! Integrate function given function values sampled at even grid points
  function space_integrate_points(space, func_vals) result(integral)
    class(space_dim) :: space
    double precision, dimension(space%num) :: func_vals
    double precision integral

    ! Encapsulate actual method for easy switching
    integral = space%trapezoid_rule(func_vals)

  end function space_integrate_points

  function trapezoid_rule(space, func_vals) result(integral)
    class(space_dim) :: space
    double precision, dimension(space%num) :: func_vals
    double precision integral

    integral = 0.5d0 * space%spacing * sum(func_vals)
  end function

  subroutine space_set_bounds(space, minval, maxval)
    class(space_dim) :: space
    double precision minval, maxval
    space%minval = minval
    space%maxval = maxval
  end subroutine space_set_bounds

  subroutine space_set_num(space, num)
    class(space_dim) :: space
    integer num
    space%num = num
  end subroutine space_set_num

  subroutine space_set_spacing(space, spacing)
    class(space_dim) :: space
    integer spacing
    space%spacing = spacing
  end subroutine space_set_spacing

  subroutine assign_linspace(space)
    class(space_dim) :: space
    integer i

    allocate(space%vals(space%num))

    space%spacing = spacing_from_num(space%minval, space%maxval, space%num)

    do i=1, space%num
       space%vals(i) = space%minval + dble(i-1) * space%spacing
    end do
  end subroutine assign_linspace

  subroutine set_spacing_from_num(space)
    class(space_dim) :: space
    space%spacing = spacing_from_num(space%minval, space%maxval, space%num)
  end subroutine set_spacing_from_num

  subroutine set_num_from_spacing(space)
    class(space_dim) :: space
    space%num = num_from_spacing(space%minval, space%maxval, space%spacing)
  end subroutine set_num_from_spacing

  subroutine space_deinit(space)
    class(space_dim) :: space
    deallocate(space%vals)
  end subroutine space_deinit

  !! SAG !!

  subroutine sag_set_bounds(grid, xmin, xmax, ymin, ymax, zmin, zmax)
    class(space_angle_grid) :: grid
    double precision xmin, xmax, ymin, ymax, zmin, zmax

    double precision, parameter :: thetamin = 0, thetamax = 2*pi
    double precision, parameter :: phimin = 0, phimax = pi

    call grid%x%set_bounds(xmin, xmax)
    call grid%y%set_bounds(ymin, ymax)
    call grid%z%set_bounds(zmin, zmax)
    call grid%theta%set_bounds(thetamin, thetamax)
    call grid%phi%set_bounds(phimin, phimax)
  end subroutine sag_set_bounds

  subroutine sag_set_spacing(grid, dx, dy, dz)
    class(space_angle_grid) :: grid
    double precision dx, dy, dz
    grid%x%spacing = dx
    grid%y%spacing = dy
    grid%z%spacing = dz
  end subroutine sag_set_spacing

  subroutine sag_set_num(grid, nx, ny, nz, ntheta, nphi)
    class(space_angle_grid) :: grid
    integer nx, ny, nz, ntheta, nphi
    call grid%x%set_num(nx)
    call grid%y%set_num(ny)
    call grid%z%set_num(nz)
    call grid%theta%set_num(ntheta)
    call grid%phi%set_num(nphi)
  end subroutine sag_set_num

  subroutine sag_init(grid)
    class(space_angle_grid) :: grid

    call grid%x%assign_linspace()
    call grid%y%assign_linspace()
    call grid%z%assign_linspace()

    call grid%theta%set_bounds(0.d0, 2.d0*pi)
    call grid%phi%set_bounds(0.d0, pi)
    call grid%theta%assign_legendre()
    call grid%phi%assign_legendre()

  end subroutine sag_init

  function sag_integrate_angle_2d(grid, func_vals) result(integral)
    class(space_angle_grid) grid
    double precision, dimension(:,:) :: func_vals
    double precision integral
    double precision prefactor
    integer lp, mp

    prefactor = grid%theta%prefactor * grid%phi%prefactor
    integral = 0

    do lp=1, grid%theta%num
       do mp=1, grid%phi%num
          integral = integral + prefactor * func_vals(lp, mp)
       end do
    end do

  end function sag_integrate_angle_2d

  subroutine sag_set_spacing_from_num(grid)
    class(space_angle_grid) :: grid
    call grid%x%set_spacing_from_num()
    call grid%y%set_spacing_from_num()
    call grid%z%set_spacing_from_num()
  end subroutine sag_set_spacing_from_num

  subroutine sag_set_num_from_spacing(grid)
    class(space_angle_grid) :: grid
    call grid%x%set_num_from_spacing()
    call grid%y%set_num_from_spacing()
    call grid%z%set_num_from_spacing()

  end subroutine sag_set_num_from_spacing
  
  subroutine sag_deinit(grid)
    class(space_angle_grid) :: grid
    call grid%x%deinit()
    call grid%y%deinit()
    call grid%z%deinit()
    call grid%theta%deinit()
    call grid%phi%deinit()
  end subroutine sag_deinit

  ! Affine shift on x from [xmin, xmax] to [ymin, ymax]
  subroutine affine_transform(x, xmin, xmax, ymin, ymax)
    double precision x, xmin, xmax, ymin, ymax
    x = ymin + (ymax-ymin)/(xmax-xmin) * (x-xmin)
  end subroutine affine_transform

  function num_from_spacing(xmin, xmax, dx) result(n)
    double precision xmin, xmax, dx
    integer n
    n = floor( (xmax - xmin) / dx )
  end function num_from_spacing

  function spacing_from_num(xmin, xmax, nx) result(dx)
    double precision xmin, xmax, dx
    integer nx
    dx = (xmax - xmin) / dble(nx)
  end function spacing_from_num
end module sag
