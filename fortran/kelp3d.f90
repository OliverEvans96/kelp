! Kelp 3D
! Oliver Evans
! 8/31/2017

! Given superindividual/water current data at each depth, generate kelp distribution at each point in 3D space

module kelp3d

use sag
use prob

subroutine generate_grid(xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz, ntheta, nphi, grid, p_kelp)
  double precision, intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax 
  integer, intent(in) :: nx, ny, nz, ntheta, nphi
  type(space_angle_grid), intent(out) :: grid
  double precision, dimension(:,:,:), allocatable, intent(out) :: p_kelp

  grid%set_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
  grid%set_num(nx, ny, nz, ntheta, nphi)

  allocate(p_kelp(nx,ny,nz))

end subroutine generate_grid

subroutine deinit(grid, p_kelp)
  grid%deinit()
  deallocate(p_kelp)
end subroutine

subroutine calculate_kelp_on_grid(grid, p_kelp, fs, fr)
  type(space_angle_grid) :: grid
  double precision, dimension(grid%nx, grid%ny, grid%nz) :: p_kelp
  double precision theta_p, r_p
  integer i, j, k

  do i=1, grid%nx
     do j=1, grid%ny
        do k=1, grid%nz
           cartesian_to_polar(grid%x(i), grid%y(j), theta_p, r_p)
           p_kelp = prob_kelp(theta_p, r_p, fs, fr)
        end do
     end do
  end do
end subroutine calculate_kelp_on_grid

function length_distribution_cdf(L, L_mean, L_std) result(output)
! C_L(L)
  double precision, intent(in) :: L, L_mean, L_std
  double precision, intent(out) :: output

  call normal_cdf(L, L_mean, L_std, output)
end function length_distribution_cdf

function angle_distribution_pdf(theta_f, theta_w, v_w) result output
! P_{\theta_f}(\theta_f)
  double precision, intent(in) :: theta_f, v_w, theta_w
  double precision, intent(out) :: output
  call von_mises_pdf(theta_f, theta_w, v_w, output
end function angle_distribution_pdf

function alpha(fs, fr)
  double precision, intent(in) :: fs, fr
  double precision, intent(out) :: alpha

  alpha = atan(tan_alpha(fs, fr))

end function alpha

function tan_alpha(fs, fr)
  double precision, intent(in) :: fs, fr
  double precision, intent(out) :: tan_alpha

  tan_alpha = 2*fs*fr / (1 + fs)

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

function shading_region_limits(theta_low_lim, theta_high_lim, theta_p, fs, fr)
  double precision, intent(in) :: theta_p, fs, fr
  double precision, intent(out) :: theta_low_lim, theta_high_lim

  theta_low_lim = theta_p - alpha(fs, fr)
  theta_high_lim = theta_p + alpha(fs, fr)
end function shading_region_limits

function prob_kelp(theta_p, r_p, fs, fr)
! P_s(theta_p, r_p)
  double precision, intent(in) :: theta_p, r_p, fs, fr
  integer, intent(in) :: quadrature_degree
  double precision, intent(out) :: prob_kelp

  call shading_region_limits(theta_low_lim, theta_high_lim, theta_p, fs, fr)
  prob_kelp = integrate_ps(theta_low_lim, theta_high_lim, quadrature_degree, theta_p, r_p, fs, fr)
end function prob_kelp

function integrate_ps(theta_low_lim, theta_high_lim, quadrature_degree, theta_p, r_p, fs, fr) result(integral)
  double precision, intent(in) :: theta_low_lim, theta_high_lim
  integer, intent(in) :: quadrature_degree
  double precision, intent(in) :: theta_p, r_p, fs, fr
  double precision, intent(out) :: integral
  double precision, dimension(:), allocatable :: integrand_vals
  integer i

  type(angle_dim) :: theta_f

  allocate(integrand_vals(theta_f%num))

  theta_f%set_bounds(theta_low_lim, theta_high_lim)
  theta_f%set_num(quadrature_degree)

  do i=1, theta_f%num
    integrand_vals(i) = ps_integrand(theta_f(i), theta_p, r_p, fs, fr)
  end do

  integral = theta_f%integrate_points(integrand_vals)

  deallocate(integrand_vals)
  theta_f%deinit()

end function integrate_ps


function ps_integrand(theta_f, theta_p, r_p, fs, fr)
  double precision theta_f, theta_p, l_min, r_p
  double precision angular_part, length_part
  double precision ps_integrand

  l_min = min_shading_length(theta_f, theta_p, r_p, fs, fr)

  angular_part = angle_distribution_pdf(theta_f) 
  length_part =  1 - length_distribution_cdf(l_min)

  ps_integrand = angular_part * length_part
end function ps_integrand


function min_shading_length(theta_f, theta_p, r_p, fs, fr) result(l_min)
! L_min(\theta)
  double precision, intent(in) :: theta_f, theta_p, r_p, fs, fr
  double precision, intent(out) :: l_min
  double precision tpp

  ! tpp = theta_p_prime
  tpp = theta_p - theta_f + pi / 2.d0
  l_min = r_p * (sin(tpp) + angular_sign(tpp) * tan_alpha(fs, fr) * cos(tpp))
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


end module
