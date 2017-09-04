program test_grid
  use kelp_context

  double precision thetaint, phiint
  double precision, external :: f1, f2
  type(space_angle_grid) grid

  call grid%set_bounds(-1.d0,1.d0,-2.d0,2.d0,-3.d0,3.d0)
  call grid%set_num(20,40,60, 5, 5)
  call grid%init()

  write(*,*) 'x'
  call print_array(grid%x, grid%nx, 1)
  write(*,*) 'y'
  call print_array(grid%y, grid%ny, 1)
  write(*,*) 'z'
  call print_array(grid%z, grid%nz, 1)
  write(*,*) 'theta'
  call print_array(grid%theta, grid%ntheta, 1)
  write(*,*) 'phi'
  call print_array(grid%phi, grid%nphi, 1)

  ! Test integration
  ! \int_0^2pi sin(theta/2) dtheta = 4
  ! \int_0^pi sin(phi) dphi = 2

  ! Hard
  thetaint = pi * sum(grid%theta_weights * sin(grid%theta/2))
  phiint = pi/2.d0 * sum(grid%phi_weights * sin(grid%phi))
  write(*,*) 'thetaint = ', thetaint
  write(*,*) 'phiint = ', phiint

  ! Easier
  write(*,*) 'thetaint = ', grid%theta_integrate_points(sin(grid%theta/2))
  write(*,*) 'phiint = ', grid%phi_integrate_points(sin(grid%phi))

  ! Easiest
  write(*,*) 'thetaint = ', grid%theta_integrate_func(f1)
  write(*,*) 'phiint = ', grid%phi_integrate_func(f2)


  call grid%deinit()
end program test_grid

function f1(x)
  double precision x, f1
  f1 = sin(x/2)
end function f1

function f2(x)
  double precision x, f2
  f2 = sin(x)
end function
