module pykelp3d_wrap
  use kelp3d
  implicit none

contains
  subroutine py_gen_kelp(xmin, xmax, nx, ymin, ymax, ny, zmax, nz, &
      frond_lengths, frond_stds, num_fronds, water_speeds, water_angles, &
      fs, fr, ft, p_kelp)

    integer nx, ny, nz
    double precision xmin, xmax, ymin, ymax, zmax
    double precision, dimension(nz) :: frond_lengths, frond_stds, &
        water_speeds, water_angles, num_fronds
    double precision fs, fr, ft
    double precision, dimension(nx, ny, nz) :: p_kelp
    integer quadrature_degree

    type(space_angle_grid) grid
    type(rope_state) rope
    type(frond_shape) frond

    quadrature_degree = 5

    ! if(.not. present(quadrature_degree)) then
    !    quadrature_degree = 5
    ! endif

    ! INIT GRID
    grid%x%minval = xmin
    grid%x%maxval = xmax
    grid%x%num = nx

    grid%y%minval = ymin
    grid%y%maxval = ymax
    grid%y%num = ny

    grid%z%minval = 0.d0
    grid%z%maxval = zmax
    grid%z%num = nz

    ! Doesn't actually matter, since we'll have to regenerate
    ! the grid for the RTE part since we can't pass derived types here.
    grid%theta%num = 10
    grid%phi%num = 10

    call grid%set_spacing_from_num()
    call grid%init()

    ! INIT ROPE
    call rope%init(grid)

    rope%frond_lengths = frond_lengths
    rope%frond_stds = frond_stds
    rope%num_fronds = num_fronds
    rope%water_speeds = water_speeds
    rope%water_angles = water_angles

    ! INIT FROND
    call frond%set_shape(fs, fr, ft)

    ! CALCULATE KELP
    call calculate_kelp_on_grid(grid, p_kelp, frond, rope, quadrature_degree)

    call rope%deinit()
    call grid%deinit()

  end subroutine py_gen_kelp

end module pykelp3d_wrap
