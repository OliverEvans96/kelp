! Kelp 3D
! Oliver Evans
! 8/31/2017

! Given superindividual/water current data at each depth, generate kelp distribution at each point in 3D space

module kelp3d

use kelp_context

implicit none

contains

subroutine generate_grid(xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz, ntheta, nphi, grid, p_kelp)
  double precision xmin, xmax, ymin, ymax, zmin, zmax
  integer nx, ny, nz, ntheta, nphi
  type(space_angle_grid) grid
  double precision, dimension(:,:,:), allocatable :: p_kelp

  call grid%set_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
  call grid%set_num(nx, ny, nz, ntheta, nphi)

  allocate(p_kelp(nx,ny,nz))

end subroutine generate_grid

subroutine kelp3d_deinit(grid, rope, p_kelp)
  type(space_angle_grid) grid
  type(rope_state) rope
  double precision, dimension(:,:,:), allocatable :: p_kelp
  call rope%deinit()
  call grid%deinit()
  deallocate(p_kelp)
end subroutine kelp3d_deinit

subroutine calculate_kelp_on_grid(grid, p_kelp, frond, rope, quadrature_degree)
  type(space_angle_grid), intent(in) :: grid
  type(frond_shape), intent(in) :: frond
  type(rope_state), intent(in) :: rope
  type(point3d) point
  integer, intent(in) :: quadrature_degree
  double precision, dimension(grid%x%num, grid%y%num, grid%z%num) :: p_kelp
  type(depth_state) depth

  integer i, j, k, nx, ny, nz
  double precision x, y, z
  ! Number of periodic images
  ! to consider in each horizontal direction
  ! for kelp distribution
  ! n_images=1 => 3x3 meta-grid
  ! n_images=2 => 5x5 meta-grid (probably not necessary)
  integer n_images
  integer im_i, im_j
  double precision x_width, y_width

  x_width = grid%x%maxval - grid%x%minval
  y_width = grid%y%maxval - grid%y%minval

  n_images = 1

  nx = grid%x%num
  ny = grid%y%num
  nz = grid%z%num

  p_kelp(:,:,:) = 0

  ! !$omp parallel do default(none) private(point,depth) &
  ! !$omp private(i,j,k,im_i,im_j) shared(nx,ny,nz,n_images) &
  ! !$omp shared(frond,rope,grid,quadrature_degree) &
  ! !$omp num_threads(num_threads) !collapse(2)
  do k=1, nz
    z = grid%z%vals(k)
    call depth%set_depth(rope, grid, k)
    do im_i=-n_images, n_images
      do im_j=-n_images, n_images
        do i=1, nx
          x = im_i*x_width + grid%x%vals(i)
          do j=1, ny
            y = im_j*y_width + grid%y%vals(j)
            call point%set_cart(x, y, z)
            p_kelp(i, j, k) = p_kelp(i,j,k) &
                  + kelp_proportion(point, frond, grid, depth, quadrature_degree)
         end do
        end do
      end do
    end do
  end do
end subroutine calculate_kelp_on_grid

subroutine shading_region_limits(theta_low_lim, theta_high_lim, point, frond)
  type(point3d), intent(in) :: point
  type(frond_shape), intent(in) :: frond
  double precision, intent(out) :: theta_low_lim, theta_high_lim

  theta_low_lim = point%theta - frond%alpha
  theta_high_lim = point%theta + frond%alpha
end subroutine shading_region_limits

function prob_kelp(point, frond, depth, quadrature_degree)
! P_s(theta_p, r_p) - This is the proportion of the population of this depth layer which can be found in this Cartesian grid cell.
  type(point3d), intent(in) :: point
  type(frond_shape), intent(in) :: frond
  type(depth_state), intent(in) :: depth
  integer, intent(in) :: quadrature_degree
  double precision prob_kelp
  double precision theta_low_lim, theta_high_lim

  call shading_region_limits(theta_low_lim, theta_high_lim, point, frond)
  prob_kelp = integrate_ps(theta_low_lim, theta_high_lim, quadrature_degree, point, frond, depth)
end function prob_kelp

function kelp_proportion(point, frond, grid, depth, quadrature_degree)
  ! This is the proportion of the volume of the Cartesian grid cell occupied by kelp
  type(point3d), intent(in) :: point
  type(frond_shape), intent(in) :: frond
  type(depth_state), intent(in) :: depth
  type(space_angle_grid), intent(in) :: grid
  integer, intent(in) :: quadrature_degree
  double precision p_k, n, t, dz
  double precision kelp_proportion

  n = depth%num_fronds
  dz = grid%z%spacing(depth%depth_layer)
  t = frond%ft
  !write(*,*) 'KELP PROPORTION'
  !write(*,*) 'n=', n
  !write(*,*) 'dz=', dz
  !write(*,*) 't=', t
  !write(*,*) 'coef=', n*t/dz
  p_k = prob_kelp(point, frond, depth, quadrature_degree)
  kelp_proportion = n*t/dz * p_k
end function kelp_proportion

function integrate_ps(theta_low_lim, theta_high_lim, quadrature_degree, point, frond, depth) result(integral)
  type(point3d), intent(in) :: point
  type(frond_shape), intent(in) :: frond
  double precision, intent(in) :: theta_low_lim, theta_high_lim
  integer, intent(in) :: quadrature_degree
  type(depth_state), intent(in) :: depth
  double precision integral
  double precision, dimension(:), allocatable :: integrand_vals
  integer i

  type(angle_dim) :: theta_f
  call theta_f%set_bounds(theta_low_lim, theta_high_lim)
  call theta_f%set_num(quadrature_degree)
  call theta_f%assign_legendre()

  allocate(integrand_vals(theta_f%num))

  do i=1, theta_f%num
    integrand_vals(i) = ps_integrand(theta_f%vals(i), point, frond, depth)
  end do

  integral = theta_f%integrate_points(integrand_vals)

  deallocate(integrand_vals)
  call theta_f%deinit()

end function integrate_ps

function ps_integrand(theta_f, point, frond, depth)
  type(point3d), intent(in) :: point
  type(frond_shape), intent(in) :: frond
  type(depth_state), intent(in) :: depth
  double precision theta_f, l_min
  double precision angular_part, length_part
  double precision ps_integrand

  l_min = min_shading_length(theta_f, point, frond)

  angular_part = depth%angle_distribution_pdf(theta_f)
  length_part =  1 - depth%length_distribution_cdf(l_min)

  ps_integrand = angular_part * length_part
end function ps_integrand


function min_shading_length(theta_f, point, frond) result(l_min)
! L_min(\theta)
  type(point3d), intent(in) :: point
  type(frond_shape), intent(in) :: frond
  double precision, intent(in) :: theta_f
  double precision l_min
  double precision tpp
  double precision frond_frac

  ! tpp === theta_p_prime
  tpp = point%theta - theta_f + pi / 2.d0
  frond_frac = 2.d0 * frond%fr / (1.d0 + frond%fs)
  l_min = point%r * (sin(tpp) + angular_sign(tpp) * frond_frac * cos(tpp))
end function min_shading_length

! function frond_edge(theta, theta_f, L, fs, fr)
! ! r_f(\theta)
!   double precision, intent(in) :: theta, theta_f, L, fs, fr
!   double precision, intent(out) :: frond_edge
!
!   frond_edge = relative_frond_edge(theta - theta_f + pi/2.d0)
!
! end function frond_edge
!
! function relative_frond_edge(theta_prime, L, fs, fr)
! ! r_f'(\theta')
!   double precision, intent(in) :: theta_prime, L, fs, fr
!   double precision, intent(out) :: relative_frond_edge
!
!   relative_frond_edge = L / (sin(theta_prime) + angular_sign(theta_prime * alpha(fs, fr) * cos(theta_prime)))
! end function relative_frond_edge

function angular_sign(theta_prime)
! S(\theta')
  double precision, intent(in) :: theta_prime
  double precision angular_sign

  ! This seems to be incorrect in summary.pdf as of 9/9/18
  ! In the report, it's written as sgn(theta_print - pi/2.d0)
  ! This results in L_min < 0 - not good!
  angular_sign = sgn(pi/2.d0 - theta_prime)
end function angular_sign

end module kelp3d

