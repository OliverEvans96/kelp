program test_grid
  use kelp_context

  double precision thetaint, phiint
  double precision, external :: f
  type(space_angle_grid) grid

  call grid%set_bounds(-1.d0,1.d0,-2.d0,2.d0,-3.d0,3.d0)
  call grid%set_num(2,4,6, 10, 10)
  call grid%init()

  call print_array(grid%x%vals, grid%x%num, 1)
  call print_array(grid%y%vals, grid%y%num, 1)
  call print_array(grid%z%vals, grid%z%num, 1)
  call print_array(grid%angles%theta, grid%angles%ntheta, 1)
  call print_array(grid%angles%phi, grid%angles%nphi, 1)

  ! Test integration
  ! \int_0^2pi sin(theta/2) dtheta = 4
  ! \int_0^pi sin(phi) dphi = 2

  ! Hard
  ! thetaint = pi * sum(grid%theta%weights * sin(grid%theta%vals/2))
  ! phiint = pi/2.d0 * sum(grid%phi%weights * sin(grid%phi%vals))
  ! write(*,*) 'thetaint = ', thetaint
  ! write(*,*) 'phiint = ', phiint

  ! ! Easier
  ! write(*,*) 'thetaint = ', grid%theta%integrate_points(sin(grid%theta%vals/2))
  ! write(*,*) 'phiint = ', grid%phi%integrate_points(sin(grid%phi%vals))

  ! ! Easiest
  ! write(*,*) 'thetaint = ', grid%theta%integrate_func(f1)
  ! write(*,*) 'phiint = ', grid%phi%integrate_func(f2)


  write(*,*) 'angle_int = ', grid%angles%integrate_func(f)

  call grid%deinit()
end program test_grid

function f(theta, phi)
  double precision theta, phi, f
  f = 1.d0 + 0.1d0*(sin(10.d0*phi)+cos(10.d0*theta))
end function f


! function f1(x)
!   double precision x, f1
!   f1 = sin(x/2)
! end function f1
! 
! function f2(x)
!   double precision x, f2
!   f2 = sin(x)
! end function
