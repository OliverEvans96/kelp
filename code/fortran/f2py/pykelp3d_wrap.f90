module pykelp3d_wrap
  !use kelp_context
  !use light_context
  !use rte_sparse_matrices
  use rte3d
  use kelp3d
  use asymptotics
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

  subroutine calculate_light_field( &
       xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz, ntheta, nphi, &
       a_w, a_k, b, num_vsf, vsf_angles, vsf_vals, &
       theta_s, phi_s, max_rad, decay, &
       tol_abs, tol_rel, maxiter_inner, maxiter_outer, &
       p_kelp, radiance, irradiance, num_scatters, &
       sparse_flag, solver_callback)

    integer nx, ny, nz, ntheta, nphi
    double precision xmin, xmax, ymin, ymax, zmin, zmax
    double precision a_w, a_k, b
    integer num_vsf
    double precision, dimension(num_vsf) :: vsf_angles, vsf_vals
    double precision theta_s, phi_s, max_rad, decay
    double precision, dimension(nx, ny, nz) :: p_kelp
    double precision, dimension(nx, ny, nz, ntheta*(nphi-2)+2) :: radiance
    double precision, dimension(nx, ny, nz) :: irradiance
    double precision tol_abs, tol_rel
    integer maxiter_inner, maxiter_outer

    integer num_scatters
    ! Rely on external function to solve sparse matrix
    ! (e.g. from Julia or Python)
    logical sparse_flag
    procedure(solver_interface), optional :: solver_callback

    type(space_angle_grid) grid
    type(rte_mat) mat
    type(optical_properties) iops
    type(light_state) light
    type(boundary_condition) bc
    integer k


    !! DUMMY VARIABLES
    ! Explicit constants won't do.
    ! f2py needs typed variables.
    integer, parameter ::  tmp_n_total=1
    integer, parameter ::  tmp_nonzero=1
    integer, dimension(tmp_nonzero) :: tmp_row
    integer, dimension(tmp_nonzero) :: tmp_col
    double precision, dimension(tmp_nonzero) :: tmp_data
    double precision, dimension(tmp_nonzero) :: tmp_sol
    double precision, dimension(tmp_n_total) :: tmp_rhs
    integer :: tmp_maxiter_outer
    integer :: tmp_maxiter_inner
    double precision :: tmp_tol_abs
    double precision :: tmp_tol_rel
    ! TODO: FINISH THIS, ADD UNDERSCORES
    ! Have to explicitly have a call to solver_callback
    ! in this function so that f2py treats it as a callback.
    tmp_row = (/ 1 /)
    tmp_col = (/ 1 /)
    tmp_data = (/ 0.d0 /)
    tmp_sol = (/ 0.d0 /)
    tmp_rhs = (/ 0.d0 /)
    tmp_maxiter_inner = 0
    tmp_maxiter_outer = 0
    tmp_tol_abs = 0.d0
    tmp_tol_rel = 0.d0
    if(.false.) then
       call solver_callback(tmp_n_total, tmp_nonzero, &
           tmp_row, tmp_col, tmp_data, tmp_sol, tmp_rhs, &
           tmp_maxiter_outer, tmp_maxiter_inner, &
           tmp_tol_abs, tmp_tol_rel)
   end if


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

    write(*,*) 'BC'
    call bc%init(grid, theta_s, phi_s, decay, max_rad)

    write(*,*) 'Scatter'
    call calculate_light_with_scattering(grid, bc, iops, radiance, num_scatters)

    if(sparse_flag) then

      ! INIT MAT
      write(*,*) 'Sparse Matrix'
      ! Set boundary condition
      call mat%init(grid, iops)
      call mat%set_bc(bc)
      call gen_matrix(mat)

      ! Set sparse solver and params
      if(present(solver_callback)) then
         mat%solver => solver_callback
      end if
      call mat%set_solver_params( &
           maxiter_outer, maxiter_inner, &
           tol_abs, tol_rel)

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

    write(*,*) 'deinit'
    call bc%deinit()
    call iops%deinit()
    call light%deinit()
    call grid%deinit()

    write(*,*) 'done'
  end subroutine calculate_light_field

end module pykelp3d_wrap
