program test_rte3d
  use test_kelp3d_mod
  use rte3d
  implicit none
  double precision, dimension(:,:,:), allocatable :: p_kelp
  type(space_angle_grid) grid
  type(rope_state) rope

  call gen_kelp(grid, rope, p_kelp)

  call deinit(grid, rope, p_kelp)

end program test_rte3d
