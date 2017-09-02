module sag
! NONE OF THESE GRIDS WILL INCLUDE THE UPPER ENDPOINT.
! To use:
! call grid%set_bounds(...)
! call grid%set_num(...) (or set_spacing)
! call grid%init()

integer, parameter :: pi = 3.141592653589793D+00

type space_angle_grid !(sag)
  integer nx, ny, nz, ntheta, nphi
  double precision :: xmin, xmax, ymin, ymax, zmin, zmax
  double precision :: dx, dy, dz, dtheta, dphi
  ! Azimuthal angle
  double precision, private :: thetamin = 0, thetamax = 2*pi
  ! Polar angle
  double precision, private :: phimin = 0, phimax = pi

  double precision, dimension(:), allocatable :: x, y, z, theta, phi
contains
  procedure :: set_bounds => sag_set_bounds
  procedure :: set_num => sag_set_num
  procedure :: init => sag_init
  procedure :: set_num_from_spacing
  procedure :: set_spacing_from_num
end type space_angle_grid


contains
  subroutine sag_set_bounds(grid, xmin, xmax, ymin, ymax, zmin, zmax)
    class(space_angle_grid) :: grid
    grid%xmin = xmin
    grid%zmin = ymin
    grid%ymin = zmin
    grid%xmax = xmax
    grid%zmax = ymax
    grid%ymax = zmax
  end subroutine sag_set_bounds

  subroutine sag_set_spacing(grid, dx, dy, dz, dtheta, dphi)
    class(space_angle_grid) :: grid
    grid%dx = dx
    grid%dy = dy
    grid%dz = dz
    grid%dtheta = dtheta
    grid%dphi = dphi

    call grid%set_num_from_spacing()
  end subroutine sag_set_spacing

  subroutine sag_set_num(grid, nx, ny, nz, ntheta, nphi)
    class(space_angle_grid) :: grid
    grid%nx = nx
    grid%ny = ny
    grid%nz = nz
    grid%ntheta = ntheta
    grid%nphi = nphi

    call grid%set_spacing_from_num()
  end subroutine sag_set_num

  subroutine sag_init(grid)
    class(space_angle_grid) :: grid

    call assign_linspace(grid%x, grid%xmin, grid%xmax, grid%dx)
    call assign_linspace(grid%y, grid%ymin, grid%ymax, grid%dy)
    call assign_linspace(grid%z, grid%zmin, grid%zmax, grid%dz)
    call assign_linspace(grid%theta, grid%thetamin, grid%thetamax, grid%dtheta)
    call assign_linspace(grid%phi, grid%phimin, grid%phimax, grid%dphi)
  end subroutine sag_init
  
  subroutine set_num_from_spacing(grid)
    class(space_angle_grid) :: grid
    grid%nx = num_from_spacing(grid%xmin, grid%xmax, grid%dx)
    grid%ny = num_from_spacing(grid%ymin, grid%ymax, grid%dy)
    grid%nz = num_from_spacing(grid%zmin, grid%zmax, grid%dz)
    grid%ntheta = num_from_spacing(grid%thetamin, grid%thetamax, grid%dtheta)
    grid%nphi = num_from_spacing(grid%phimin, grid%phimax, grid%dphi)
  end subroutine set_num_from_spacing

  subroutine set_spacing_from_num(grid)
    class(space_angle_grid) :: grid
    grid%dx= spacing_from_num(grid%xmin, grid%xmax, grid%nx)
    grid%dy= spacing_from_num(grid%ymin, grid%ymax, grid%ny)
    grid%dz= spacing_from_num(grid%zmin, grid%zmax, grid%nz)
    grid%dtheta= spacing_from_num(grid%thetamin, grid%thetamax, grid%ntheta)
    grid%dphi= spacing_from_num(grid%phimin, grid%phimax, grid%nphi)
  end subroutine set_spacing_from_num

  subroutine assign_linspace(arr, xmin, xmax, dx)
    double precision, dimension(:), allocatable :: arr
    double precision xmin, xmax, dx
    integer n

    n = num_from_spacing(xmin, xmax, dx)

    allocate(arr(n))

    do i=1, n
       arr(i) = xmin + dble(i-1) * dx
    end do
  end subroutine assign_linspace

  function num_from_spacing(xmin, xmax, dx) result(n)
    double precision xmin, xmax, dx
    n = floor( (xmax - xmin) / dx )
  end function num_from_spacing

  function spacing_from_num(xmin, xmax, nx) result(dx)
    double precision xmin, xmax
    integer n
    dx = floor( (xmax - xmin) / dble(n) )
  end function spacing_from_num
end module sag
