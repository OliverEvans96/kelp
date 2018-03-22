module sag
use utils
use fastgl

implicit none

! Spatial grids do not include upper endpoints.
! Angular grids do include upper endpoints.
! Both include lower endpoints.

! To use:
! call grid%set_bounds(...)
! call grid%set_num(...) (or set_uniform_spacing)
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

type angle2d
   integer ntheta, nphi, nomega
   double precision dtheta, dphi
   double precision, dimension(:), allocatable :: theta, phi, theta_edge, phi_edge
   double precision, dimension(:), allocatable ::  theta_p, phi_p, theta_edge_p, phi_edge_p
   double precision, dimension(:), allocatable :: cos_theta, sin_theta, cos_phi, sin_phi
   double precision, dimension(:), allocatable :: cos_theta_edge, sin_theta_edge, cos_phi_edge, sin_phi_edge
 contains
   procedure :: set_num => angle_set_num
   procedure :: phat, lhat, mhat
   procedure :: init => angle_init ! Call after set_num
   procedure :: integrate_points => angle_integrate_points
   procedure :: integrate_func => angle_integrate_func
   procedure :: deinit => angle_deinit
end type angle2d

type space_dim
   integer num
   double precision minval, maxval
   double precision, dimension(:), allocatable :: vals, edges, spacing
 contains
   procedure :: integrate_points => space_integrate_points
   procedure :: trapezoid_rule
   procedure :: set_bounds => space_set_bounds
   procedure :: set_num => space_set_num
   procedure :: set_uniform_spacing => space_set_uniform_spacing
   !procedure :: set_num_from_spacing
   procedure :: set_uniform_spacing_from_num
   procedure :: set_spacing_array => space_set_spacing_array
   procedure :: deinit => space_deinit
   procedure :: assign_linspace
end type space_dim

type space_angle_grid !(sag)
  type(space_dim) :: x, y, z
  type(angle2d) :: angles
  double precision, dimension(:), allocatable :: x_factor, y_factor
contains
  procedure :: set_bounds => sag_set_bounds
  procedure :: set_num => sag_set_num
  procedure :: init => sag_init
  procedure :: deinit => sag_deinit
  !procedure :: set_num_from_spacing => sag_set_num_from_spacing
  procedure :: set_uniform_spacing_from_num => sag_set_uniform_spacing_from_num
  procedure :: calculate_factors => sag_calculate_factors
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

  subroutine angle_set_num(angles, ntheta, nphi)
    class(angle2d) :: angles
    integer ntheta, nphi
    angles%ntheta = ntheta
    angles%nphi = nphi
    angles%nomega = ntheta*(nphi-2) + 2
  end subroutine angle_set_num

  function lhat(angles, p) result(l)
    class(angle2d) :: angles
    integer l, p
    l = mod1(p, angles%ntheta)
  end function lhat

  function mhat(angles, p) result(m)
    class(angle2d) :: angles
    integer m, p
    m = ceiling(dble(p)/dble(angles%ntheta)) + 1
  end function mhat

  function phat(angles, l, m) result(p)
    class(angle2d) :: angles
    integer l, m, p

    p = (m-2)*angles%ntheta + l
  end function phat

  subroutine angle_init(angles)
    class(angle2d) :: angles
    integer l, m, p

    allocate(angles%theta(angles%ntheta))
    allocate(angles%phi(angles%nphi))
    allocate(angles%theta_edge(angles%ntheta))
    allocate(angles%phi_edge(angles%nphi-1))
    allocate(angles%theta_p(angles%nomega))
    allocate(angles%phi_p(angles%nomega))
    allocate(angles%theta_edge_p(angles%nomega))
    allocate(angles%phi_edge_p(angles%nomega))
    allocate(angles%cos_theta(angles%nomega))
    allocate(angles%sin_theta(angles%nomega))
    allocate(angles%cos_phi(angles%nomega))
    allocate(angles%sin_phi(angles%nomega))
    allocate(angles%cos_theta_edge(angles%nomega))
    allocate(angles%sin_theta_edge(angles%nomega))
    allocate(angles%cos_phi_edge(angles%nomega-1))
    allocate(angles%sin_phi_edge(angles%nomega-1))

    ! Calculate spacing
    angles%dtheta = 2.d0*pi/dble(angles%ntheta)
    angles%dphi = pi/dble(angles%nphi-1)

    ! Create grids
    do l=1, angles%ntheta
       angles%theta(l) = dble(l-1)*angles%dtheta
       angles%theta_edge(l) = dble(l-1)*angles%dtheta
       angles%cos_theta_edge(l) = cos(angles%theta_edge(l))
       angles%sin_theta_edge(l) = sin(angles%theta_edge(l))
    end do

    do m=1, angles%nphi
       angles%phi(m) = dble(m-1.d0)*angles%dphi
       if(m<angles%nphi) then
          angles%phi_edge(m) = dble(m-0.5d0)*angles%dphi
          angles%cos_phi_edge(m) = cos(angles%phi_edge(m))
          angles%sin_phi_edge(m) = sin(angles%phi_edge(m))
       end if
    end do

    ! Create p arrays
    do m=2, angles%nphi-1
       do l=1, angles%ntheta
          p = angles%phat(l, m)

          angles%theta_p(p) = angles%theta(l)
          angles%theta_edge_p(p) = angles%theta_edge(l)
          angles%phi_p(p) = angles%phi(m)
          angles%phi_edge_p(p) = angles%phi_edge(m)

       end do
    end do

    ! Poles
    l=1
    ! North Pole
    p = angles%nomega-1
    m=1
    angles%theta_p(p) = angles%theta(l)
    angles%theta_edge_p(p) = angles%theta_edge(l)
    angles%phi_p(p) = angles%phi(m)
    ! phi_edge_p only defined up to nphi-1.
    angles%phi_edge_p(p) = angles%phi_edge(m)

    ! South Pole
    p = angles%nomega
    m = angles%nphi
    angles%theta_p(p) = angles%theta(l)
    angles%theta_edge_p(p) = angles%theta_edge(l)
    angles%phi_p(p) = angles%phi(m)

  end subroutine angle_init

  ! Integrate function given function values at grid cells
  function angle_integrate_points(angles, func_vals) result(integral)
    class(angle2d) :: angles
    double precision, dimension(angles%nomega) :: func_vals
    double precision integral
    integer p, m

    integral = 0

    ! Interior
    do p=1, angles%nomega-2
       m = angles%mhat(p)
       integral = integral + func_vals(p) * ( &
       angles%cos_phi_edge(m-1) &
       - angles%cos_phi_edge(m))
    end do

    ! Scale
    integral = angles%dtheta * integral

    ! Include poles
    integral = integral + 2.d0*pi*(1.d0-cos(angles%dphi/2.d0)) * (&
         func_vals(angles%nomega-1) &
         + func_vals(angles%nomega))

  end function angle_integrate_points

  function angle_integrate_func(angles, func_callable) result(integral)
    class(angle2d) :: angles
    double precision, external :: func_callable
    double precision, dimension(:), allocatable :: func_vals
    double precision integral
    integer p
    double precision theta, phi

    allocate(func_vals(angles%nomega))

    do p=1, angles%nomega
       theta = angles%theta_p(p)
       phi = angles%phi_p(p)
       func_vals(p) = func_callable(theta, phi)
    end do

    integral = angles%integrate_points(func_vals)

    deallocate(func_vals)
  end function angle_integrate_func


  subroutine angle_deinit(angles)
    class(angle2d) :: angles
    deallocate(angles%theta)
    deallocate(angles%phi)
    deallocate(angles%theta_edge)
    deallocate(angles%phi_edge)
    deallocate(angles%theta_p)
    deallocate(angles%phi_p)
    deallocate(angles%theta_edge_p)
    deallocate(angles%phi_edge_p)
    deallocate(angles%cos_theta)
    deallocate(angles%sin_theta)
    deallocate(angles%cos_phi)
    deallocate(angles%sin_phi)
    deallocate(angles%cos_theta_edge)
    deallocate(angles%sin_theta_edge)
    deallocate(angles%cos_phi_edge)
    deallocate(angles%sin_phi_edge)
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

    integral = 0.5d0 * sum(func_vals * space%spacing)
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

  subroutine space_set_uniform_spacing(space, spacing)
    class(space_dim) :: space
    double precision spacing
    integer k
    do k=1, space%num
      space%spacing(k) = spacing
   end do
  end subroutine space_set_uniform_spacing

  subroutine space_set_spacing_array(space, spacing)
    class(space_dim) :: space
    double precision, dimension(space%num) :: spacing
    space%spacing = spacing
  end subroutine space_set_spacing_array

  subroutine assign_linspace(space)
    class(space_dim) :: space
    double precision spacing
    integer i

    allocate(space%vals(space%num))
    allocate(space%edges(space%num))
    allocate(space%spacing(space%num))

    spacing = spacing_from_num(space%minval, space%maxval, space%num)
    call space%set_uniform_spacing(spacing)

    do i=1, space%num
       space%edges(i) = space%minval + dble(i-1) * space%spacing(i)
       space%vals(i) = space%minval + dble(i-0.5d0) * space%spacing(i)
    end do

  end subroutine assign_linspace

  subroutine set_uniform_spacing_from_num(space)
    ! Create evenly spaced grid (linspace)
    class(space_dim) :: space
    integer k
    double precision spacing

    spacing = spacing_from_num(space%minval, space%maxval, space%num)
    call space%set_uniform_spacing(spacing)

  end subroutine set_uniform_spacing_from_num

 !  subroutine set_num_from_spacing(space)
 !    class(space_dim) :: space
 !    !space%num = num_from_spacing(space%minval, space%maxval, space%spacing)

 !  end subroutine set_num_from_spacing

  subroutine space_deinit(space)
    class(space_dim) :: space
    deallocate(space%vals)
    deallocate(space%edges)
    deallocate(space%spacing)
  end subroutine space_deinit

  !! SAG !!

  subroutine sag_set_bounds(grid, xmin, xmax, ymin, ymax, zmin, zmax)
    class(space_angle_grid) :: grid
    double precision xmin, xmax, ymin, ymax, zmin, zmax

    call grid%x%set_bounds(xmin, xmax)
    call grid%y%set_bounds(ymin, ymax)
    call grid%z%set_bounds(zmin, zmax)
  end subroutine sag_set_bounds

  subroutine sag_set_uniform_spacing(grid, dx, dy, dz)
    class(space_angle_grid) :: grid
    double precision dx, dy, dz
    call grid%x%set_uniform_spacing(dx)
    call grid%y%set_uniform_spacing(dy)
    call grid%z%set_uniform_spacing(dz)
  end subroutine sag_set_uniform_spacing

  subroutine sag_set_num(grid, nx, ny, nz, ntheta, nphi)
    class(space_angle_grid) :: grid
    integer nx, ny, nz, ntheta, nphi
    call grid%x%set_num(nx)
    call grid%y%set_num(ny)
    call grid%z%set_num(nz)
    call grid%angles%set_num(ntheta, nphi)
  end subroutine sag_set_num

  subroutine sag_init(grid)
    class(space_angle_grid) :: grid

    call grid%x%assign_linspace()
    call grid%y%assign_linspace()
    call grid%z%assign_linspace()

    call grid%angles%init()
    call grid%calculate_factors()

  end subroutine sag_init

  subroutine sag_calculate_factors(grid)
    ! Factors by which depth difference is multiplied
    ! in order to calculate distance traveled in the
    ! (x, y) direction along a ray in the (theta, phi)
    ! direction
    class(space_angle_grid) :: grid
    integer p, nomega
    double precision theta, phi

    nomega = grid%angles%nomega

    allocate(grid%x_factor(nomega))
    allocate(grid%y_factor(nomega))

    do p=1, nomega
       theta = grid%angles%theta_p(p)
       phi = grid%angles%phi_p(p)
       grid%x_factor(p) = tan(phi) * cos(theta)
       grid%y_factor(p) = tan(phi) * sin(theta)
    end do

  end subroutine sag_calculate_factors

  subroutine sag_set_uniform_spacing_from_num(grid)
    class(space_angle_grid) :: grid
    call grid%x%set_uniform_spacing_from_num()
    call grid%y%set_uniform_spacing_from_num()
    call grid%z%set_uniform_spacing_from_num()
  end subroutine sag_set_uniform_spacing_from_num

  ! subroutine sag_set_num_from_spacing(grid)
  !   class(space_angle_grid) :: grid
  !   call grid%x%set_num_from_spacing()
  !   call grid%y%set_num_from_spacing()
  !   call grid%z%set_num_from_spacing()

  ! end subroutine sag_set_num_from_spacing
  
  subroutine sag_deinit(grid)
    class(space_angle_grid) :: grid
    call grid%x%deinit()
    call grid%y%deinit()
    call grid%z%deinit()
    call grid%angles%deinit()

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
