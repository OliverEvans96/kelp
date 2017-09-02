! Kelp 3D
! Oliver Evans
! 8/31/2017

! Given superindividual/water current data at each depth, generate kelp distribution at each point in 3D space

use prob
use cube_grid
use intlib

subroutine generate_grid(xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz)
end subroutine generate_grid

subroutine calculate_kelp_on_grid()
end subroutine calculate_kelp_on_grid

function length_distribution_cdf(L, L_mean, L_std)
! C_L(L)
  double precision, intent(in) :: L, L_mean, L_std
  double precision, intent(out) :: length_distribution_cdf

  ! Not 100% sure about the order of arguments. Probably good.
  call normal_cdf(L, L_mean, L_std)

end function length_distribution_cdf

function angle_distribution_pdf(theta_f, theta_w, v_w)
! P_{\theta_f}(\theta_f)
  double precision, intent(in) :: theta_f, v_w, theta_w
  double precision, intent(out) :: angle_distribution_pdf
  !!! NORMAL DISTRIBUTION PDF !!!
  call von_mises_pdf(theta_f, theta_w, v_w, angle_distribution_pdf)
end function angle_distribution_pdf

function alpha(fs, fr)
  double precision, intent(in) :: fs, fr
  double precision, intent(out) :: alpha

  alpha = atan( 2*fs*fr / (1 + fs))

end function alpha

subroutine polar_to_cartesian(r, theta, x, y)
  double precision, intent(in) :: r, theta
  double precision, intent(out) :: x, y
  x = r*cos(theta)
  y = r*sin(theta)
end subroutine polar_to_cartesian

subroutine cartesian_to_polar(x, y, r, theta)
  double precision, intent(in) :: x, y
  double precision, intent(out) :: r, theta
  r = sqrt(x**2 + y**2)
  theta = atan2(x, y)
end subroutine cartesian_to_polar

function integrate_length_angle_dist()
end function integrate_length_angle_dist

function shading_region_limits()
end function shading_region_limits

function prob_kelp(theta_p, r_p)
! P_s(theta_p, r_p)
  double precision, intent(in) :: theta_p, r_p,
  double precision, intent(out) :: prob_kelp

  !!! NUMERICAL INTEGRATION !!!

  call shading_region_limits()
  call integrate_length_angle_dist()

  double precision, intent(in), external :: 

end function prob_kelp

function min_shading_length(r_p, theta, L)
! L_min(\theta)
  double precision, intent(in) :: r_p, theta, L
  double precision, intent(out) :: min_shading_length

  min_shading_length = r_p * L / frond_edge(theta)
end function min_shading_length

function frond_edge(theta, theta_f L, fs, fr)
! r_f(\theta)
  double precision, intent(in) :: theta, theta_f, L, fs, fr
  double precision, intent(out) :: frond_edge

  frond_edge = relative_frond_edge(theta - theta_f + pi/2.d0)

end function frond_edge

function relative_frond_edge(theta_prime, L, fs, fr)
! r_f'(\theta')
  double precision, intent(in) :: theta_prime, L, fs, fr
  double precision, intent(out) :: relative_frond_edge

  relative_frond_edge = L / (sin(theta_prime) + angular_sign(theta_prime * alpha(fs, fr) * cos(theta_prime)))
end function relative_frond_edge

function angular_sign(theta_prime)
! S(\theta')
  double precision, intent(in) :: theta_prime
  double precision, intent(out) :: angular_sign

  angular_sign = sgn(theta_prime - pi/2.d0)
end function angular_sign

function sgn(x)
! Standard signum function
  sgn = sign(1,x) 
  if(x .eq. 0.) sgn = 0 
end function signum


