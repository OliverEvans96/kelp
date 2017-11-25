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

type, extends(space_dim) :: angle_dim
   double precision, dimension(:), allocatable :: sin, cos, half_cos, half_cos_sum
 contains
   procedure :: calculate_trig => angle_calculate_trig
   procedure :: deinit_trig => angle_deinit_trig
end type angle_dim

type space_angle_grid !(sag)
  type(space_dim) :: x, y, z
  type(angle_dim) :: theta, phi
  double precision, dimension(:,:), allocatable :: x_factor, y_factor
contains
  procedure :: set_bounds => sag_set_bounds
  procedure :: set_num => sag_set_num
  procedure :: integrate_angle_2d => sag_integrate_angle_2d
  procedure :: init => sag_init
  procedure :: deinit => sag_deinit
  procedure :: set_num_from_spacing => sag_set_num_from_spacing
  procedure :: set_spacing_from_num => sag_set_spacing_from_num
  procedure :: calculate_factors => sag_calculate_factors
  procedure :: integrate_pole_ray
  procedure :: integrate_interior_ray
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


  !! ANGLE !!

  subroutine angle_calculate_trig(angle)
    class(angle_dim) :: angle
    integer i

    allocate(angle%sin(angle%num))
    allocate(angle%cos(angle%num))
    allocate(angle%half_cos(angle%num-1))
    allocate(angle%half_cos_sum(2:angle%num-1))

    do i=1, angle%num
       angle%sin(i) = sin(angle%vals(i))
       angle%cos(i) = cos(angle%vals(i))
    end do

    do i=1, angle%num-1
       angle%half_cos(i) = cos(angle%vals(i)+angle%spacing/2)
    end do

    do i=2, angle%num-1
       angle%half_cos_sum(i) = angle%half_cos(i-1) + angle%half_cos(i)
    end do
  end subroutine angle_calculate_trig

  subroutine angle_deinit_trig(angle)
    class(angle_dim) angle
    deallocate(angle%sin)
    deallocate(angle%cos)
    deallocate(angle%half_cos)
    deallocate(angle%half_cos_sum)
  end subroutine angle_deinit_trig


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

    write(*,*) 'Assigning x linspace'
    call grid%x%assign_linspace()
    write(*,*) 'Assigning y linspace'
    call grid%y%assign_linspace()
    write(*,*) 'Assigning z linspace'
    call grid%z%assign_linspace()

    call grid%theta%set_bounds(0.d0, 2.d0*pi)
    call grid%phi%set_bounds(0.d0, pi)
    call grid%theta%assign_linspace()
    call grid%phi%assign_linspace()

    call grid%theta%calculate_trig()
    call grid%phi%calculate_trig()

    call grid%calculate_factors()

  end subroutine sag_init

  subroutine sag_calculate_factors(grid)
    ! Factors by which depth difference is multiplied
    ! in order to calculate distance traveled in the
    ! (x, y) direction along a ray in the (theta, phi)
    ! direction
    class(space_angle_grid) :: grid
    integer l, m
    double precision theta, phi
    integer ntheta, nphi

    ntheta = grid%theta%num
    nphi = grid%phi%num

    allocate(grid%x_factor(ntheta, nphi))
    allocate(grid%y_factor(ntheta, nphi))

    do l=1, ntheta
       theta = grid%theta%vals(l)
       do m=1, nphi
          phi = grid%phi%vals(m)
          grid%x_factor(l, m) = tan(phi) * cos(theta)
          grid%y_factor(l, m) = tan(phi) * sin(theta)
       end do
    end do
  end subroutine sag_calculate_factors

  function sag_integrate_angle_2d(grid, func_vals) result(integral)
    ! Integrate a function over the unit sphere,
    ! assuming that function values are constant over cells
    ! in the uniform rectangular 2d angular grid
    class(space_angle_grid) grid
    double precision, dimension(:,:) :: func_vals
    double precision integral
    integer l, m, nphi, ntheta
    double precision dtheta, phi_slice

    nphi = grid%phi%num
    ntheta = grid%theta%num
    dtheta = grid%theta%spacing

    ! No theta dependence at poles
    integral = (1-grid%phi%half_cos(1)) &
         * (func_vals(1,1) + func_vals(1,nphi))

    ! Loop over interior
    do m=2, nphi-1
       phi_slice = 0
       do l=1, ntheta
          phi_slice = phi_slice + func_vals(l, m)
       end do
       integral = integral + phi_slice * grid%phi%half_cos_sum(m)
    end do
  end function sag_integrate_angle_2d

  function integrate_pole_ray(grid, m, func_val) result(integral)
    ! Integrate a function over a singular interior spherical rectangle
    ! m is the phi index of the pole
    ! (must be either 1 or nphi)
    class(space_angle_grid) grid
    double precision func_val
    double precision integral
    integer m
    double precision dtheta

    dtheta = grid%theta%spacing

    ! No theta dependence at poles
    integral = (1-grid%phi%half_cos(1)) * func_val
  end function integrate_pole_ray

  function integrate_interior_ray(grid, l, m, func_val) result(integral)
    ! Integrate a function over a singular interior spherical rectangle
    ! l, m are the theta and phi indices respectively
    ! of the rectangle.
    class(space_angle_grid) grid
    double precision func_val
    double precision integral
    integer l, m
    double precision dtheta

    integral = func_val * grid%phi%half_cos_sum(m)
  end function integrate_interior_ray

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
    call grid%theta%deinit_trig()
    call grid%phi%deinit_trig()

    deallocate(grid%x_factor)
    deallocate(grid%y_factor)
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
