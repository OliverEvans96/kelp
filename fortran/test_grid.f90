program test_grid
  use test_grid_mod
  implicit none
  double precision, external :: f_ang

  call test_2d_angular_integration(f_ang, 10, 10, 13.823007d0)
  call test_2d_angular_integration(f_ang, 20, 20, 12.570662d0)

  call test_angle_p_conversions(10, 10)
  call test_angle_p_conversions(15, 25)
  call test_angle_p_conversions(30, 15)
end program test_grid

function f_ang(theta, phi)
  double precision theta, phi, f_ang
  f_ang = 1.d0 + 0.1d0*(sin(10.d0*phi)+cos(10.d0*theta))
end function f_ang
