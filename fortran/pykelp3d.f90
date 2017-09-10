program pykelp3d
  use kelp3d
  use hdf5_utils
  implicit none

  type(space_angle_grid) grid
  type(rope_state) rope
  type(frond_shape) frond
  integer quadrature_degree
  double precision, dimension(:,:,:), allocatable :: p_kelp

  character(len=23), parameter :: gridfile = 'hdf5/kelp3d/grid.hdf5'
  character(len=23), parameter :: ropefile = 'hdf5/kelp3d/rope.hdf5'
  character(len=23), parameter :: kelpfile = 'hdf5/kelp3d/kelp.hdf5'
  character(len=25), parameter :: frondfile = 'hdf5/kelp3d/frond.hdf5'
  character(len=25), parameter :: paramfile = 'hdf5/kelp3d/param.hdf5'

  call hdf_read_grid(gridfile, grid)
  call hdf_read_rope(ropefile, rope, grid)
  call hdf_read_frond(frondfile, frond)
  call hdf_read_params(paramfile, quadrature_degree)

  allocate(p_kelp(grid%x%num, grid%y%num, grid%z%num))

  call calculate_kelp_on_grid(grid, p_kelp, frond, rope, quadrature_degree)
  call hdf_write_kelp(kelpfile, p_kelp, grid)

  call deinit(grid, rope, p_kelp)

end program pykelp3d
