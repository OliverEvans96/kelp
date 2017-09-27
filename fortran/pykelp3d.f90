program pykelp3d
  use kelp3d
  implicit none

  type(space_angle_grid) grid
  type(rope_state) rope
  type(frond_shape) frond
  type(light_state) light
  integer quadrature_degree
  double precision, dimension(:,:,:), allocatable :: p_kelp

  character(len=23), parameter :: gridfile = 'hdf5/kelp3d/grid.hdf5'
  character(len=23), parameter :: ropefile = 'hdf5/kelp3d/rope.hdf5'
  character(len=23), parameter :: kelpfile = 'hdf5/kelp3d/kelp.hdf5'
  character(len=25), parameter :: frondfile = 'hdf5/kelp3d/frond.hdf5'
  character(len=25), parameter :: paramfile = 'hdf5/kelp3d/param.hdf5'
  character(len=23), parameter :: radfile = 'hdf5/kelp3d/rad.hdf5'
  character(len=24), parameter :: irradfile = 'hdf5/kelp3d/irrad.hdf5'

  call hdf_read_grid(gridfile, grid)
  call hdf_read_rope(ropefile, rope, grid)
  call hdf_read_frond(frondfile, frond)
  call hdf_read_params(paramfile, quadrature_degree)

  allocate(p_kelp(grid%x%num, grid%y%num, grid%z%num))

  call calculate_kelp_on_grid(grid, p_kelp, frond, rope, quadrature_degree)
  call hdf_write_kelp(kelpfile, p_kelp, grid)

  call light%calculate_radiance()
  call light%calculate_irradance()

  call hdf_write_rad(radfile, light%radiance, grid)
  call hdf_write_irrad(irradfile, light%irradiance, grid)

  call kelp3d_deinit(grid, rope, p_kelp)

end program pykelp3d
