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

  write(*,*) 'A'
  call hdf_read_grid(gridfile, grid )
  write(*,*) 'C'
  call hdf_read_rope(ropefile, rope, grid)
  write(*,*) 'D'
  !call hdf_read_kelp(kelpfile, p_kelp)
  call hdf_read_frond(frondfile, frond)
  write(*,*) 'E'
  call hdf_read_params(paramfile, quadrature_degree)
  write(*,*) 'F'

  allocate(p_kelp(grid%x%num, grid%y%num, grid%z%num))

  call calculate_kelp_on_grid(grid, p_kelp, frond, rope, quadrature_degree)
  write(*,*) 'G'
  call hdf_write_kelp(kelpfile, p_kelp, grid)
  write(*,*) 'G'

  call deinit(grid, rope, p_kelp)
  write(*,*) 'I'

end program pykelp3d
