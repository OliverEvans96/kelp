module pyrte3d_wrap
  use rte3d
  implicit none

contains
  subroutine py_calculate_light_field( &
       xmin, xmax, nx, ymin, ymax, ny, zmax, nz, ntheta, nphi, &
       a_w, a_k, b_w, b_k, num_vsf, vsf_angles, vsf_vals, &
       theta_s, phi_s, max_rad, decay, &
       tol_abs, tol_rel, maxiter_inner, maxiter_outer, &
       p_kelp, radiance, irradiance)

    integer nx, ny, nz, ntheta, nphi
    double precision xmin, xmax, ymin, ymax, zmax
    double precision a_w, a_k, b_w, b_k
    integer num_vsf
    double precision, dimension(num_vsf) :: vsf_angles, vsf_vals
    double precision theta_s, phi_s, max_rad, decay
    double precision, dimension(nx, ny, nz) :: p_kelp
    double precision, dimension(nx, ny, nz, ntheta, nphi) :: radiance
    double precision, dimension(nx, ny, nz) :: irradiance
    double precision tol_abs, tol_rel
    integer maxiter_inner, maxiter_outer

    type(space_angle_grid) grid
    type(rte_mat) mat
    type(optical_properties) iops
    type(light_state) light

    ! INIT GRID
    grid%x%minval = xmin
    grid%x%maxval = xmax
    grid%x%num = nx

    grid%y%minval = ymin
    grid%y%maxval = ymax
    grid%y%num = ny

    grid%z%maxval = zmax
    grid%z%num = nz

    grid%theta%num = ntheta
    grid%phi%num = nphi

    call grid%set_spacing_from_num()
    call grid%init()

    ! INIT IOPS
    call iops%init(grid)
    iops%abs_kelp = a_k
    iops%scat_kelp = b_k
    iops%abs_water = a_w
    iops%scat_water = b_w
    iops%num_vsf = num_vsf
    iops%vsf_angles = vsf_angles
    iops%vsf_vals = vsf_vals
    call iops%calc_vsf_on_grid()
    call iops%calculate_coef_grids(p_kelp)

    ! INIT MAT
    call gen_matrix(grid, mat, iops)
    mat%params%maxiter_outer = maxiter_outer
    mat%params%maxiter_inner = maxiter_inner
    mat%params%tol_abs = tol_abs
    mat%params%tol_rel = tol_rel

    ! Solve system
    call light%init(mat)
    call light%calculate_radiance()
    call light%calculate_irradiance()

    radiance = light%radiance
    irradiance = light%irradiance

  end subroutine py_calculate_light_field

end module pyrte3d_wrap
