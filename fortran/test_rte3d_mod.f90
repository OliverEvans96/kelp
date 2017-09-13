module test_rte3d_mod
  use test_kelp3d_mod
  use rte3d
  implicit none
contains

  subroutine run_test_rte3d
    double precision, dimension(:,:,:), allocatable :: p_kelp
    type(space_angle_grid) grid
    type(rope_state) rope
    type(rte_mat) mat
    type(optical_properties) iops

    write(*,*) 'Start'
    call gen_kelp(grid, rope, p_kelp)
    write(*,*) 'Gen matrix'
    call gen_matrix(grid, mat, iops)

    write(*,*) 'Solve'
    call mat%solve()

    write(*,*) 'Kelp deinit'
    call kelp3d_deinit(grid, rope, p_kelp)
    write(*,*) 'RTE3D deinit'
    call rte3d_deinit(mat)

  end subroutine run_test_rte3d

end module test_rte3d_mod
