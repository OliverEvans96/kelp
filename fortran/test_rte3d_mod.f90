module test_rte3d_mod
  use test_kelp3d_mod
  use rte3d
  implicit none
contains

  subroutine run_test_rte3d
    double precision, dimension(:,:,:), allocatable :: p_kelp
    type(space_angle_grid) grid
    type(rope_state) rope
    call gen_kelp(grid, rope, p_kelp)
    call kelp3d_deinit(grid, rope, p_kelp)
  end subroutine run_test_rte3d

end module test_rte3d_mod
