module sag
use utils

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

type space_angle_grid !(sag)
  integer nx, ny, nz, ntheta, nphi
  double precision :: xmin, xmax, ymin, ymax, zmin, zmax
  double precision :: dx, dy, dz, dtheta, dphi
  ! Azimuthal angle
  double precision :: thetamin = 0, thetamax = 2*pi
  ! Polar angle
  double precision :: phimin = 0, phimax = pi

  double precision, dimension(:), allocatable :: x, y, z, theta, phi
  ! Gauss-Legendre weights for numerical integration
  double precision, dimension(:), allocatable :: theta_weights, phi_weights
  double precision theta_prefactor, phi_prefactor
contains
  procedure :: set_bounds => sag_set_bounds
  procedure :: set_num => sag_set_num
  procedure :: init => sag_init
  procedure :: deinit => sag_deinit
  !procedure :: set_num_from_spacing
  procedure :: set_spacing_from_num
  procedure :: theta_integrate_points
  procedure :: theta_integrate_func
  procedure :: phi_integrate_points
  procedure :: phi_integrate_func
end type space_angle_grid


contains
  subroutine sag_set_bounds(grid, xmin, xmax, ymin, ymax, zmin, zmax)
    class(space_angle_grid) :: grid
    double precision xmin, xmax, ymin, ymax, zmin, zmax
    grid%xmin = xmin
    grid%ymin = ymin
    grid%zmin = zmin
    grid%xmax = xmax
    grid%ymax = ymax
    grid%zmax = zmax
  end subroutine sag_set_bounds

  ! subroutine sag_set_spacing(grid, dx, dy, dz)
  !   class(space_angle_grid) :: grid
  !   grid%dx = dx
  !   grid%dy = dy
  !   grid%dz = dz
  ! 
  !   call grid%set_num_from_spacing()
  ! end subroutine sag_set_spacing

  subroutine sag_set_num(grid, nx, ny, nz, ntheta, nphi)
    class(space_angle_grid) :: grid
    integer nx, ny, nz, ntheta, nphi
    grid%nx = nx
    grid%ny = ny
    grid%nz = nz
    grid%ntheta = ntheta
    grid%nphi = nphi

    call grid%set_spacing_from_num()
  end subroutine sag_set_num

  subroutine sag_init(grid)
    class(space_angle_grid) :: grid

    call set_spacing_from_num(grid)

    call assign_linspace(grid%x, grid%xmin, grid%xmax, grid%nx)
    call assign_linspace(grid%y, grid%ymin, grid%ymax, grid%ny)
    call assign_linspace(grid%z, grid%zmin, grid%zmax, grid%nz)

    call assign_legendre(grid%theta, grid%theta_weights, grid%theta_prefactor, grid%thetamin, grid%thetamax, grid%ntheta)
    call assign_legendre(grid%phi, grid%phi_weights, grid%phi_prefactor, grid%phimin, grid%phimax, grid%nphi)

  end subroutine sag_init

  subroutine sag_deinit(grid)
    class(space_angle_grid) :: grid
    deallocate(grid%x)
    deallocate(grid%y)
    deallocate(grid%z)
    deallocate(grid%theta)
    deallocate(grid%phi)
  end subroutine sag_deinit

  ! subroutine set_num_from_spacing(grid)
  !   class(space_angle_grid) :: grid
  !   grid%nx = num_from_spacing(grid%xmin, grid%xmax, grid%dx)
  !   grid%ny = num_from_spacing(grid%ymin, grid%ymax, grid%dy)
  !   grid%nz = num_from_spacing(grid%zmin, grid%zmax, grid%dz)
  ! end subroutine set_num_from_spacing

  subroutine set_spacing_from_num(grid)
    class(space_angle_grid) :: grid
    grid%dx = spacing_from_num(grid%xmin, grid%xmax, grid%nx)
    grid%dy = spacing_from_num(grid%ymin, grid%ymax, grid%ny)
    grid%dz = spacing_from_num(grid%zmin, grid%zmax, grid%nz)
  end subroutine set_spacing_from_num

  ! subroutine assign_linspace_spacing(arr, xmin, xmax, dx)
  !   double precision, dimension(:), allocatable :: arr
  !   double precision xmin, xmax, dx
  !   integer n, i

  !   n = num_from_spacing(xmin, xmax, dx)

  !   allocate(arr(n))

  !   do i=1, n
  !      arr = xmin + dble(i-1) * dx
  !   end do
  ! end subroutine assign_linspace_spacing

  subroutine assign_linspace(arr, xmin, xmax, n)
    double precision, dimension(:), allocatable :: arr
    double precision xmin, xmax, dx
    integer n, i

    allocate(arr(n))

    dx = spacing_from_num(xmin, xmax, n)
    write(*,*) 'dx = ', dx

    do i=1, n
       arr(i) = xmin + dble(i-1) * dx
    end do
  end subroutine assign_linspace

  ! To calculate \int_{xmin}^{xmax} f(x) dx :
  ! int = prefactor * sum(weights * f(roots))
  subroutine assign_legendre(roots, weights, prefactor, xmin, xmax, nx)
    double precision, dimension(:), allocatable :: roots, weights
    double precision xmin, xmax, prefactor, root, weight, theta
    integer nx, i
    ! glpair produces both x and theta, where x=cos(theta). We'll throw out theta.

    allocate(roots(nx))
    allocate(weights(nx))

    ! Prefactor for integration
    ! From change of variables
    prefactor = (xmax - xmin) / 2

    do i = 1, nx
       call glpair(nx, i, theta, weight, root)
       call affine_transform(root, -1.d0, 1.d0, xmin, xmax)
       roots(i) = root
       weights(i) = weight
    end do
  end subroutine assign_legendre
  
  ! Integrate callable function over theta
  function theta_integrate_func(grid, func_callable) result(integral)
    class(space_angle_grid) :: grid
    double precision, external :: func_callable
    double precision, dimension(:), allocatable :: func_vals
    double precision integral
    integer i

    allocate(func_vals(grid%ntheta))

    do i=1, grid%ntheta
       func_vals(i) = func_callable(grid%theta(i))
    end do

    integral = grid%theta_integrate_points(func_vals)

    deallocate(func_vals)
  end function theta_integrate_func

  ! Integrate function given function values sampled at legendre theta values
  function theta_integrate_points(grid, func_vals) result(integral)
    class(space_angle_grid) :: grid
    double precision, dimension(grid%ntheta) :: func_vals
    double precision integral

    integral = grid%theta_prefactor * sum(grid%theta_weights * func_vals)
  end function theta_integrate_points
  
  ! Integrate callable function over phi
  function phi_integrate_func(grid, func_callable) result(integral)
    class(space_angle_grid) :: grid
    double precision, external :: func_callable
    double precision, dimension(:), allocatable :: func_vals
    double precision integral
    integer i

    allocate(func_vals(grid%nphi))

    do i=1, grid%nphi
       func_vals(i) = func_callable(grid%phi(i))
    end do

    integral = grid%phi_integrate_points(func_vals)

    deallocate(func_vals)
  end function phi_integrate_func

  ! Integrate function given function values sampled at legendre phi values
  function phi_integrate_points(grid, func_vals) result(integral)
    class(space_angle_grid) :: grid
    double precision, dimension(grid%nphi) :: func_vals
    double precision integral

    integral = grid%phi_prefactor * sum(grid%phi_weights * func_vals)
  end function phi_integrate_points

  ! Affine shift on x from [xmin, xmax] to [ymin, ymax]
  subroutine affine_transform(x, xmin, xmax, ymin, ymax)
    double precision x, xmin, xmax, ymin, ymax
    x = ymin + (ymax-ymin)/(xmax-xmin) * (x-xmin)
  end subroutine affine_transform

  ! function num_from_spacing(xmin, xmax, dx) result(n)
  !   double precision xmin, xmax, dx
  !   n = floor( (xmax - xmin) / dx )
  ! end function num_from_spacing

  function spacing_from_num(xmin, xmax, nx) result(dx)
    double precision xmin, xmax, dx
    integer nx
    dx = (xmax - xmin) / dble(nx)
  end function spacing_from_num
end module sag
