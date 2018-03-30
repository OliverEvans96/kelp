module test_grid
  use sag
  implicit none

contains

  function test_2d_angular_integration(func, ntheta, nphi) result(integral)
    type(space_angle_grid) grid
    double precision, external :: func
    integer ntheta, nphi
    double precision integral

    call grid%set_bounds(-1.d0,1.d0,-2.d0,2.d0,-3.d0,3.d0)
    call grid%set_num(4,4,4, ntheta, nphi)
    call grid%init()

    integral = grid%angles%integrate_func(func)

    !if(passed) then
    !   write(*,'(A,I2,A,I2)') '* PASS - test_2d_angular_integration: ', ntheta, 'x', nphi
    !else
    !   write(*,'(A,I2,A,I2)') '* FAIL - test_2d_angular_integration: ', ntheta, 'x', nphi
    !end if

    call grid%deinit()
  end function test_2d_angular_integration

  function test_angle_p_conversions(ntheta, nphi) result(passed)
  integer ntheta, nphi
  logical passed
  integer l, m, p
  type(space_angle_grid) grid

  call grid%set_bounds(-1.d0,1.d0,-2.d0,2.d0,-3.d0,3.d0)
  call grid%set_num(4,4,4, ntheta, nphi)
  call grid%init()

  passed = .true.

  ! Make sure that conversions between l,m and p
  ! are invertible on interior of angular grid
  do l=1, ntheta
     do m=2, nphi-1
        p = grid%angles%phat(l, m)
        if(grid%angles%lhat(p) .ne. l) then
           passed = .false.
           exit
        end if
        if(grid%angles%mhat(p) .ne. m) then
           passed = .false.
           exit
        end if
     end do
     if(.not. passed) then
        exit
     end if
  end do

  !if(passed) then
  !   write(*,'(A,I2,A,I2)') '* PASS - test_p_conversions: ', ntheta, 'x', nphi
  !else
  !   write(*,'(A,I2,A,I2)') '* FAIL - test_p_conversions: ', ntheta, 'x', nphi
  !end if

  call grid%deinit()

  end function test_angle_p_conversions

  function make_fortran_grid(xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz, ntheta, nphi) result(grid)
    double precision, intent(in) ::  xmin, xmax, ymin, ymax, zmin, zmax
    integer, intent(in) :: nx, ny, nz, ntheta, nphi
    type(space_angle_grid) grid

    call grid%set_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
    call grid%set_num(nx, ny, nz, ntheta, nphi)
    call grid%init()
  end function make_fortran_grid

  subroutine destruct_fortran_grid(grid)
    type(space_angle_grid) grid
    call grid%deinit()
  end subroutine destruct_fortran_grid

  subroutine make_grid(&
       xmin, xmax, ymin, ymax, zmin, zmax,&
       nx, ny, nz, ntheta, nphi,&
       dx, dy, dz, dtheta, dphi,&
       x, y, z, theta, phi,&
       x_edge, y_edge, z_edge, theta_edge, phi_edge)
    double precision, intent(in) ::  xmin, xmax, ymin, ymax, zmin, zmax
    integer, intent(in) :: nx, ny, nz, ntheta, nphi
    double precision, intent(out) :: dtheta, dphi
    double precision, intent(out), dimension(nx) :: x, x_edge, dx
    double precision, intent(out), dimension(ny) :: y, y_edge, dy
    double precision, intent(out), dimension(nz) :: z, z_edge, dz
    double precision, intent(out), dimension(ntheta) :: theta
    double precision, intent(out), dimension(ntheta) :: theta_edge
    double precision, intent(out), dimension(nphi) :: phi
    double precision, intent(out), dimension(nphi-1) :: phi_edge

    type(space_angle_grid) grid

    call grid%set_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
    call grid%set_num(nx, ny, nz, ntheta, nphi)
    call grid%init()

    dx = grid%x%spacing
    dy = grid%y%spacing
    dz = grid%z%spacing
    dtheta = grid%angles%dtheta
    dphi = grid%angles%dphi

    x = grid%x%vals
    y = grid%y%vals
    z = grid%z%vals
    theta = grid%angles%theta
    phi = grid%angles%phi

    x_edge = grid%x%edges
    y_edge = grid%y%edges
    z_edge = grid%z%edges
    theta_edge = grid%angles%theta_edge
    phi_edge = grid%angles%phi_edge
    call grid%deinit()
  end subroutine make_grid

end module test_grid
