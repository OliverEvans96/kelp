module sag
! NONE OF THESE GRIDS WILL INCLUDE THE UPPER ENDPOINT.
! To use:
! call grid%set_bounds(...)
! call grid%set_num(...) (or set_spacing)
! call grid%init()
type space_angle_grid !(sag)
  integer nx, ny, nz, ntheta, nphi
  double precision :: xmin, xmax, ymin, ymax, zmin, zmax
  double precision :: dx, dy, dz, dtheta, dphi
  ! Azimuthal angle
  double precision, parameter :: thetamin = 0, thetamax = 2*pi
  ! Polar angle
  double precision, parameter :: phimin = 0, phimax = pi

  double precision, dimension(:), allocatable :: x, y, z, theta, phi
contains
  procedure :: set_bounds => sag_set_bounds
  procedure :: set_num => sag_set_num
  procedure :: set_angle => sag_set_angle
  procedure :: init => sag_init
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
  end subroutine init_grid
  
  subroutine set_num_from_spacing(grid)
    class(space_angle_grid) :: grid
    grid%nx = num_from_spacing(grid%dx)
    grid%ny = num_from_spacing(grid%dy)
    grid%nz = num_from_spacing(grid%dz)
    grid%ntheta = num_from_spacing(grid%dtheta)
    grid%nphi = num_from_spacing(grid%dphi)
  end subroutine set_num_from_spacing

  subroutine set_spacing_from_num(grid)
    class(space_angle_grid) :: grid
    grid%dx= num_from_spacing(grid%nx)
    grid%dy= num_from_spacing(grid%ny)
    grid%dz= num_from_spacing(grid%nz)
    grid%dtheta= num_from_spacing(grid%ntheta)
    grid%dphi= num_from_spacing(grid%nphi)
  end subroutine set_spacing_from_num

  subroutine assign_linspace(arr, xmin, xmax, dx)
    double precision, dimension(:), allocatable :: arr
    double precision xmin, xmax
    integer n

    n = n_from_dx(xmin, xmax, dx)

    allocate(arr(n))

    do i=1, n
       arr(i) = xmin + (i-1) * dx
    end do
  end subroutine assign_linspace

  function num_from_spacing(xmin, xmax, dx) result(n)
    double precision xmin, xmax, dx
    n = floor( (xmax - xmin) / dx )
  end function n_from_dx

  function spacing_from_num(xmin, xmax, nx) results dx
    double precision xmin, xmax
    integer n
    dx = floor( (xmax - xmin) / float(n) )
  end function spacing_from_num
end module sag
