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
    type(boundary_condition) bc

    open(unit=7, file='frad_in.txt')
    open(unit=8, file='frad_out.txt')

    ! Report arguments
    write(*,*) 'xmin =', xmin
    write(*,*) 'xmax =', xmax
    write(*,*) 'nx =', nx
    write(*,*) 'ymin =', ymin
    write(*,*) 'ymax =', ymax
    write(*,*) 'ny =', ny
    write(*,*) 'zmax =', zmax
    write(*,*) 'nz =', nz
    write(*,*) 'ntheta =', ntheta
    write(*,*) 'nphi =', nphi
    write(*,*) 'a_w =', a_w
    write(*,*) 'a_k =', a_k
    write(*,*) 'b_w =', b_w
    write(*,*) 'b_k =', b_k
    write(*,*) 'num_vsf =', num_vsf
    write(*,*) 'vsf_angles =', vsf_angles
    write(*,*) 'vsf_vals =', vsf_vals
    write(*,*) 'theta_s =', theta_s
    write(*,*) 'phi_s =', phi_s
    write(*,*) 'max_rad =', max_rad
    write(*,*) 'decay =', decay
    write(*,*) 'tol_abs =', tol_abs
    write(*,*) 'tol_rel =', tol_rel
    write(*,*) 'maxiter_inner =', maxiter_inner
    write(*,*) 'maxiter_outer =', maxiter_outer
    !write(*,*) 'p_kelp =', p_kelp
    write(7,*) radiance
    !write(*,*) 'irradiance =', irradiance


    ! INIT GRID
    write(*,*) 'Grid'
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
    write(*,*) 'IOPs'
    write(*,*) 'coefs'
    iops%abs_kelp = a_k
    iops%scat_kelp = b_k
    iops%abs_water = a_w
    iops%scat_water = b_w
    iops%num_vsf = num_vsf
    write(*,*) 'init'
    call iops%init(grid)
    iops%vsf_angles = vsf_angles
    iops%vsf_vals = vsf_vals
    write(*,*) 'VSF'
    call iops%calc_vsf_on_grid()
    write(*,*) 'Coef grid'
    call iops%calculate_coef_grids(p_kelp)

    ! INIT MAT
    write(*,*) 'MAT'

    ! Set boundary condition
    bc%max_rad = max_rad
    bc%decay = decay
    bc%theta_s = theta_s
    bc%phi_s = phi_s

    write(*,*) 'mat init'
    call mat%init(grid, iops)
    call mat%set_bc(bc, radiance)
    call gen_matrix(mat)
    write(*,*) 'RHS_SUM =', sum(mat%rhs)
    write(*,*) 'Gen matrix'
    write(*,*) 'Set params'
    mat%params%maxiter_outer = maxiter_outer
    mat%params%maxiter_inner = maxiter_inner
    mat%params%tol_abs = tol_abs
    mat%params%tol_rel = tol_rel

    ! Initialize & set initial guess
    write(*,*) 'Light init'
    call light%init(mat)

    write(*,*) 'set radiance'
    light%radiance = radiance

    ! Solve system
    write(*,*) 'calculate radiance'
    call light%calculate_radiance()
    write(*,*) 'calculate irradiance'
    call light%calculate_irradiance()

    write(*,*) 'save radiance'
    radiance = light%radiance
    write(*,*) 'save irradiance'
    irradiance = light%irradiance


    write(*,*) 'deinit'
    write(*,*) 'ready'
    call iops%deinit()
    write(*,*) 'a'
    call light%deinit()
    write(*,*) 'b'
    call mat%deinit()
    write(*,*) 'c'
    call grid%deinit()

    write(*,*) 'writing new radiance to file 8'
    write(8,*) radiance

    close(7)
    close(8)

    write(*,*) 'done'
  end subroutine py_calculate_light_field

end module pyrte3d_wrap
