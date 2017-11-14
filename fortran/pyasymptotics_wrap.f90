module pyasymptotics_wrap
  use asymptotics
  implicit none

contains
  subroutine py_calculate_asymptotic_light_field( &
       xmin, xmax, nx, ymin, ymax, ny, zmax, nz, ntheta, nphi, &
       a_w, a_k, b_w, b_k, num_vsf, vsf_angles, vsf_vals, &
       theta_s, phi_s, max_rad, decay, &
       tol_abs, tol_rel, maxiter_inner, maxiter_outer, &
       p_kelp, radiance, irradiance, num_scatters, gmres_flag)

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

    integer num_scatters
    logical gmres_flag

    type(space_angle_grid) grid
    type(rte_mat) mat
    type(optical_properties) iops
    type(light_state) light
    type(boundary_condition) bc

    ! INIT GRID
    write(*,*) 'Grid'
    grid%x%minval = xmin
    grid%x%maxval = xmax
    grid%x%num = nx

    grid%y%minval = ymin
    grid%y%maxval = ymax
    grid%y%num = ny

    grid%z%minval = 0.d0
    grid%z%maxval = zmax
    grid%z%num = nz

    grid%theta%num = ntheta
    grid%phi%num = nphi

    call grid%set_spacing_from_num()
    call grid%init()

    ! INIT IOPS
    write(*,*) 'IOPs'
    iops%abs_kelp = a_k
    iops%scat_kelp = b_k
    iops%abs_water = a_w
    iops%scat_water = b_w
    iops%num_vsf = num_vsf

    call iops%init(grid)
    iops%vsf_angles = vsf_angles
    iops%vsf_vals = vsf_vals

    call iops%calc_vsf_on_grid()
    call iops%calculate_coef_grids(p_kelp)

    write(*,*) 'BC'
    call bc%init(grid, theta_s, phi_s, decay, max_rad)

    write(*,*) 'Scatter'
    call calculate_light_with_scattering(grid, bc, iops, radiance, num_scatters)

    !! I THINK WE CAN DO WITHOUT COPYING RADIANCE !!

    if(gmres_flag) then

      ! INIT MAT
       write(*,*) 'GMRES Mat'
      ! Set boundary condition
      call mat%init(grid, iops)
      call mat%set_bc(bc, radiance)
      call gen_matrix(mat)
      mat%params%maxiter_outer = maxiter_outer
      mat%params%maxiter_inner = maxiter_inner
      mat%params%tol_abs = tol_abs
      mat%params%tol_rel = tol_rel

      ! Initialize & set initial guess
      write(*,*) 'Light'
      light%radiance = radiance

      ! Solve system
      write(*,*) 'Calculate Radiance'
      call light%calculate_radiance()

      call mat%deinit()
    else
       call light%init_grid(grid)
       light%radiance = radiance
    endif

    write(*,*) 'Irrad'
    call light%calculate_irradiance()

    radiance = light%radiance
    irradiance = light%irradiance

    write(*,*) 'deinit'
    call bc%deinit()
    call iops%deinit()
    call light%deinit()
    call grid%deinit()

    write(*,*) 'done'
  end subroutine py_calculate_asymptotic_light_field

end module pyasymptotics_wrap
