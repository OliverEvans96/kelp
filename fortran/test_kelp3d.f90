module test_kelp3d_mod
  use kelp3d
  use hdf5_utils
  implicit none
contains

subroutine init_grid(grid)
  type(space_angle_grid) grid

  grid%x%minval = -1.d0
  grid%x%maxval = 1.d0
  grid%x%spacing = .5d-1

  grid%y%minval = -1.d0
  grid%y%maxval = 1.d0
  grid%y%spacing = .5d-1

  grid%z%minval = 0.d0
  grid%z%maxval = 1.d1
  grid%z%spacing = 1.d0

  grid%theta%num = 20
  grid%phi%num = 10

  call grid%set_num_from_spacing()
  call grid%init()

end subroutine init_grid

subroutine init_rope(rope, grid)
  type(rope_state) rope
  type(space_angle_grid) grid
  double precision, dimension(:), allocatable :: z, ones
  integer nz
  double precision zmax
  double precision a, b, c

  integer i

  call rope%init(grid)

  nz = grid%z%num
  zmax = grid%z%maxval

  allocate(z(nz))
  allocate(ones(nz))

  z = grid%z%vals
  ones = z / z

  ! Make up some goofy shape
  a = 2.d-1
  b = 1.d-1
  c = 5.d-1
  rope%frond_lengths = exp(-a*z) * sin(z) ** 2
  rope%frond_stds = b * ones
  rope%water_speeds = c * ones
  rope%water_angles = 2*pi / z(nz) * z

  deallocate(z)
  deallocate(ones)
end subroutine init_rope

subroutine init_kelp(p_kelp, grid)
  double precision, dimension(:,:,:), allocatable :: p_kelp
  type(space_angle_grid) grid
  integer nx, ny, nz

  nx = grid%x%num
  ny = grid%y%num
  nz = grid%z%num

  allocate(p_kelp(nx,ny,nz))

end subroutine init_kelp

subroutine init_frond(frond)
  type(frond_shape) frond
  double precision fs, fr
  fs = .5d0
  fr = 2.d0
  call frond%set_shape(fs, fr)
end subroutine

subroutine init_params(quadrature_degree)
  integer quadrature_degree
  quadrature_degree = 5
end subroutine init_params

subroutine calculate_memory_usage(grid)
  type(space_angle_grid) grid
  integer nx, ny, nz, ntheta, nphi, space_bytes, total_bytes
  double precision space_megabytes, total_megabytes

  nx = grid%x%num
  ny = grid%y%num
  nz = grid%z%num
  ntheta = grid%theta%num
  nphi = grid%phi%num

  space_bytes = 8 * nx * ny * nz
  total_bytes = space_bytes * ntheta * nphi

  space_megabytes = dble(space_bytes) / dble(1024 ** 2)
  total_megabytes = dble(total_bytes) / dble(1024 ** 2)

  write(*,'(A,F10.3,A)') 'Space: ', space_megabytes, ' MB'
  write(*,'(A,F10.3,A)') 'Total: ', total_megabytes, ' MB'
end subroutine


subroutine run_test_kelp_3d(grid, rope, p_kelp, theta)
  type(space_angle_grid) grid
  type(rope_state) rope
  type(frond_shape) frond
  integer quadrature_degree
  double precision, dimension(:,:,:), allocatable :: p_kelp

  type(depth_state) depth
  type(angle_dim) theta
  double precision thetamin, thetamax
  integer ntheta
  double precision pth
  integer i

  call init_grid(grid)
  call init_rope(rope, grid)
  call init_kelp(p_kelp, grid)
  call init_frond(frond)
  call init_params(quadrature_degree)

  call calculate_memory_usage(grid)
  call calculate_kelp_on_grid(grid, p_kelp, frond, rope, quadrature_degree)

  thetamin = -3 * pi
  thetamax = 3 * pi
  ntheta = 61

  call theta%set_bounds(thetamin, thetamax)
  call theta%set_num(ntheta)
  call theta%assign_legendre()

  call depth%set_depth(rope, grid, 1)

  do i=1, theta%num
    pth = depth%angle_distribution_pdf(theta%vals(i))
    write(*,*) 'Theta = ', theta%vals(i)
    write(*,*) 'P_th = ', pth
    write(*,*)
  end do

end subroutine run_test_kelp_3d
end module test_kelp3d_mod

!! Main Program !!

program test_kelp3d
  use test_kelp3d_mod
  implicit none
  double precision, dimension(:,:,:), allocatable :: p_kelp
  type(space_angle_grid) grid
  type(rope_state) rope
  type(angle_dim) theta

  call run_test_kelp_3d(grid, rope, p_kelp, theta)

  call theta%deinit()
  call deinit(grid, rope, p_kelp)

end program test_kelp3d
