module pykelp3d_wrap
  !use kelp_context
  !use light_context
  !use rte_sparse_matrices
  use rte3d
  use kelp3d
  use asymptotics
  use light_interface
  implicit none

contains
  subroutine gen_kelp(xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz, &
      frond_lengths, frond_stds, num_fronds, water_speeds, water_angles, &
      fs, fr, ft, p_kelp)

    integer nx, ny, nz
    double precision xmin, xmax, ymin, ymax, zmin, zmax
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

    call grid%set_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
    ! ntheta, nphi don't actually matter, since we'll have to regenerate
    ! the grid for the RTE part since we can't pass derived types here.
    call grid%set_num(nx, ny, nz, 4, 4)
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

  end subroutine gen_kelp

  ! TODO: max_rad -> surface_irrad
  ! TODO: Calculate avg_irrad, perceied_irrad
  ! using functions from light_interface.f90
  subroutine calculate_light_field( &
       xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz, ntheta, nphi, &
       a_w, a_k, b, num_vsf, vsf_angles, vsf_vals, &
       theta_s, phi_s, max_rad, decay, &
       p_kelp, radiance, irradiance, avg_irrad, perceived_irrad, &
       num_scatters, fd_flag, lis_opts, lis_iter, lis_time, lis_resid)

    integer nx, ny, nz, ntheta, nphi
    double precision xmin, xmax, ymin, ymax, zmin, zmax
    double precision a_w, a_k, b
    integer num_vsf
    double precision, dimension(num_vsf) :: vsf_angles, vsf_vals
    double precision theta_s, phi_s, max_rad, decay
    double precision, dimension(nx, ny, nz) :: p_kelp
    double precision, dimension(nx, ny, nz, ntheta*(nphi-2)+2) :: radiance
    double precision, dimension(nx, ny, nz) :: irradiance
    double precision, dimension(nz) :: avg_irrad, perceived_irrad
    character*(*) :: lis_opts
    integer lis_iter
    double precision :: lis_time, lis_resid

    integer num_scatters
    ! Rely on external function to solve sparse matrix
    ! (e.g. from Julia or Python)
    logical fd_flag

    type(space_angle_grid) grid
    type(rte_mat) mat
    type(optical_properties) iops
    type(light_state) light
    type(boundary_condition) bc
    integer k

    ! INIT GRID
    write(*,*) 'Grid'
    call grid%set_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
    call grid%set_num(nx, ny, nz, ntheta, nphi)
    call grid%init()

    ! INIT IOPS
    write(*,*) 'IOPs'
    iops%num_vsf = num_vsf

    call iops%init(grid)
    iops%vsf_angles = vsf_angles
    iops%vsf_vals = vsf_vals
    call iops%calc_vsf_on_grid()

    iops%abs_kelp = a_k
    do k=1, grid%z%num
       iops%abs_water(k) = a_w
    end do
    iops%scat = b
    call iops%calculate_coef_grids(p_kelp)

    write(*,*) 'max abs =', maxval(iops%abs_grid)
    write(*,*) 'max loc =', maxloc(iops%abs_grid)

    write(*,*) 'BC'
    call bc%init(grid, theta_s, phi_s, decay, max_rad)

    write(*,*) 'Scatter'
    call calculate_light_with_scattering(grid, bc, iops, radiance, num_scatters)

    if(fd_flag) then

      ! INIT MAT
      write(*,*) 'Sparse Matrix'
      ! Set boundary condition
      call mat%init(grid, iops)
      call mat%set_bc(bc)
      call gen_matrix(mat)

      ! Set sparse solver and params
      !if(present(solver_callback)) then
      !   mat%solver => solver_callback
      !end if
      ! call mat%set_solver_params( &
      !      maxiter_outer, maxiter_inner, &
      !      tol_abs, tol_rel)

      call mat%set_solver_opts(lis_opts)

      ! Initialize & set initial guess
      write(*,*) 'Light'
      call light%init(mat)
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

    ! TODO: Make sure this works
    call lis_solver_get_iters(mat%solver, lis_iter)
    call lis_solver_get_time(mat%solver, lis_time)
    call lis_solver_get_residualnorm(mat%solver, lis_resid)

    call calculate_average_irradiacnce(grid, light, avg_irrad)
    call calculate_perceived_irradiance(grid, p_kelp, &
         perceived_irrad, irradiance)

    write(*,*) 'deinit'
    call bc%deinit()
    call iops%deinit()
    call light%deinit()
    call grid%deinit()

    write(*,*) 'done'
  end subroutine calculate_light_field

end module pykelp3d_wrap
