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

subroutine calculate_kelp_on_grid(grid, p_kelp, frond, rope, quadrature_degree, n_images, num_threads)
  type(space_angle_grid), intent(in) :: grid
  type(frond_shape), intent(in) :: frond
  type(rope_state), intent(in) :: rope
  type(point3d) point
  integer, intent(in) :: quadrature_degree
  integer, optional :: n_images
  double precision, dimension(grid%x%num, grid%y%num, grid%z%num) :: p_kelp
  type(depth_state) depth
  integer num_threads

  integer i, j, k, nx, ny, nz
  double precision x, y, z
  ! Number of periodic images
  ! to consider in each horizontal direction
  ! for kelp distribution
  ! n_images=1 => 3x3 meta-grid
  ! n_images=2 => 5x5 meta-grid (only necessary for very dense kelp ropes)
  integer im_i, im_j
  double precision x_width, y_width

  x_width = grid%x%maxval - grid%x%minval
  y_width = grid%y%maxval - grid%y%minval

  if(.not. present(n_images)) then
    n_images = 1
  end if

  nx = grid%x%num
  ny = grid%y%num
  nz = grid%z%num

  p_kelp(:,:,:) = 0

  !$omp parallel do default(shared) private(x,y,z) &
  !$omp firstprivate(point,depth) &
  !$omp private(i,j,k,im_i,im_j) shared(nx,ny,nz,n_images) &
  !$omp shared(frond,rope,grid,quadrature_degree) &
  !$omp shared(p_kelp,x_width,y_width) &
  !$omp num_threads(num_threads) collapse(3) &
  !$omp schedule(dynamic, 10) ! 10 grid points per thread
  do k=1, nz
    do i=1, nx
      do j=1, ny
        z = grid%z%vals(k)
        call depth%set_depth(rope, grid, k)
        do im_i=-n_images, n_images
          x = im_i*x_width + grid%x%vals(i)
          do im_j=-n_images, n_images
            y = im_j*y_width + grid%y%vals(j)
            call point%set_cart(x, y, z)
            p_kelp(i, j, k) = p_kelp(i,j,k) &
                  + kelp_proportion(point, frond, grid, depth, quadrature_degree)
         end do
        end do
      end do
    end do
  end do
  !omp end do
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

subroutine gaussian_blur_2d(A, sigma, dx, dy, nk, num_threads)
  ! 2D Gaussian blur (periodic BC) with std sigma
  ! with kernel radius of nk (full size (2*nk+1)x(2*nk+1))
  ! applied to matrix A with element spacings dx and dy.
  double precision, intent(inout), dimension(:, :) :: A
  double precision, intent(in) :: sigma, dx, dy
  ! kernel half width
  integer, intent(in) :: nk
  ! kernel full width
  integer kw
  integer num_threads

  ! A matrix size
  integer nx, ny

  ! indices
  integer i1, j1
  integer i2, j2
  integer i, j
  ! kernel
  double precision, dimension(:,:), allocatable :: k
  ! output matrix
  double precision, dimension(:,:), allocatable :: B
  ! kernel independent variables
  double precision x, y

  if(sigma > 0) then
    nx = size(A, 1)
    ny = size(A, 2)

    kw = 2*nk + 1

    allocate(B(nx, ny))
    allocate(k(kw, kw))
    write(*,*) 'creating kernel', sigma, nk
    ! Create kernel
    do i1=-nk, nk
      x = i1*dx
      i = i1+nk+1
      do j1=-nk, nk
          y = j1*dy
          j = j1+nk+1
          k(i,j) = exp(-(x**2+y**2)/(2*sigma**2))
      end do

    end do
    ! normalize kernel
    k = k / sum(k)

    write(*,*) 'convolving'
    ! convolve
    !$omp parallel do default(private) private(x,y) &
    !$omp private(i,j,i1,j1,i2,j2) shared(nx,ny,nk,kw) &
    !$omp shared(A,B,k) &
    !$omp num_threads(num_threads) collapse(2) &
    !$omp schedule(dynamic, 10) ! 10 grid points per thread
    do i1=1, nx
      do j1=1, ny
          B(i1, j1) = 0
          do i2=1, kw
            do j2=1, kw
                i = mod1(i1 - nk + i2 - 1, nx)
                j = mod1(j1 - nk + j2 - 1, ny)
                B(i1, j1) = B(i1, j1) + k(i2, j2) * A(i, j)
            end do
          end do
      end do
    end do
    !omp end parallel do
    write(*,*) 'done convolving'

    ! Update original matrix
    A(:,:) = B(:,:)
    deallocate(k)
    deallocate(B)
    write(*,*) 'gb2d done.'
 end if
end subroutine gaussian_blur_2d

end module kelp3d

