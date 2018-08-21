module light_utils_wrap
  use light_context

contains

  subroutine calculate_irradiance(radiance, irradiance, nx, ny, nz, ntheta, nphi)
    integer, intent(in) :: nx, ny, nz, ntheta, nphi
    double precision, intent(in), dimension(nx, ny, nz, (nphi-2)*ntheta+2) :: radiance
    double precision, intent(inout), dimension(:,:,:) :: irradiance

    type(space_angle_grid) grid
    type(optical_properties) iops
    type(rte_mat) mat
    type(light_state) light

    ! Not important for calculating irradiance
    double precision xmin, xmax, ymin, ymax, zmin, zmax
    xmin = 0.d0
    xmax = 1.d0
    ymin = 0.d0
    ymax = 1.d0
    zmin = 0.d0
    zmax = 1.d0

    ! INIT GRID
    call grid%set_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
    call grid%set_num(nx, ny, nz, ntheta, nphi)
    call grid%init()

    call light%init_grid(grid)
    light%radiance = radiance
    call light%calculate_irradiance()
    irradiance = light%irradiance

    call grid%deinit()
    call light%deinit()
  end subroutine calculate_irradiance

end module light_utils_wrap
