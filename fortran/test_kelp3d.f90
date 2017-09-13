program test_kelp3d
  use test_kelp3d_mod
  implicit none
  double precision, dimension(:,:,:), allocatable :: p_kelp
  type(space_angle_grid) grid
  type(rope_state) rope
  type(angle_dim) theta

  call test_vm_dist(grid, rope, p_kelp, theta)

  call theta%deinit()
  call deinit(grid, rope, p_kelp)

end program test_kelp3d
