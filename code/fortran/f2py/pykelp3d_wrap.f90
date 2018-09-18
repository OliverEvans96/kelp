module pykelp3d_wrap
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

  ! TODO: Error on odd nphi
  subroutine calculate_light_field( &
       xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz, ntheta, nphi, &
       a_w, a_k, b, num_vsf, vsf_angles, vsf_vals, &
       theta_s, phi_s, I0, decay, &
       p_kelp, radiance, irradiance, &
       num_scatters, num_threads, fd_flag, &
       lis_opts, lis_iter, lis_time, lis_resid)

    integer, intent(in) :: nx, ny, nz, ntheta, nphi
    double precision, intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax
    double precision, intent(in) :: a_w, a_k, b
    integer, intent(in) ::num_vsf
    double precision, dimension(num_vsf), intent(in) :: vsf_angles, vsf_vals
    double precision, intent(in) :: theta_s, phi_s, I0, decay
    double precision, dimension(nx, ny, nz), intent(in) :: p_kelp
    double precision, dimension(nx, ny, nz, ntheta*(nphi-2)+2), intent(inout) :: radiance
    double precision, dimension(nx, ny, nz), intent(inout) :: irradiance

    integer, intent(in) :: num_scatters
    integer, intent(in) :: num_threads
    logical, intent(in) :: fd_flag
    character*(*), intent(in) :: lis_opts

    integer, intent(inout) :: lis_iter
    double precision, intent(inout) :: lis_time, lis_resid

    type(space_angle_grid) grid
    type(rte_mat) mat
    type(optical_properties) iops
    type(light_state) light
    type(boundary_condition) bc
    integer k

    double precision, dimension(:,:,:,:), allocatable :: source

    ! INIT GRID
    write(*,*) 'Grid'
    call grid%set_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
    call grid%set_num(nx, ny, nz, ntheta, nphi)
    call grid%init()

    allocate(source( &
         grid%x%num, &
         grid%y%num, &
         grid%z%num, &
         grid%angles%nomega))

    ! Initialize source to zero
    source(:,:,:,:) = 0.d0

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
    write(*,*) 'min abs =', maxval(iops%abs_grid)
    write(*,*) 'mean abs =', sum(iops%abs_grid)/size(iops%abs_grid)

    write(*,*) 'BC'
    call bc%init(grid, theta_s, phi_s, decay, I0)

    !write(*,*) 'bc_grid = ', bc%bc_grid

    write(*,*) 'Calculate asymptotic light field'
    call calculate_asymptotic_light_field(&
         grid, bc, iops, source, &
         radiance, num_scatters, num_threads)

    if(fd_flag) then

      ! INIT MAT
      write(*,*) 'Sparse Matrix'
      ! Set boundary condition
      call mat%init(grid, iops)
      call mat%set_bc(bc)
      call gen_matrix(mat, num_threads)

      ! Set solver options
      call mat%set_solver_opts(lis_opts)

      ! Initialize & set initial guess
      write(*,*) 'Light'
      call light%init(mat)
      light%radiance = radiance

      ! Solve system
      write(*,*) 'Calculate Radiance'
      call light%calculate_radiance()

      call mat%get_solver_stats(lis_iter, lis_time, lis_resid)
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

    deallocate(source)

    write(*,*) 'done'
  end subroutine calculate_light_field

  ! This subroutine allows manual specification
  ! of the following from callback arguments.
  ! - absorption coefficient (x, y, z)
  ! - source function (x, y, z, theta, phi)
  ! - surface boundary condition (theta, phi)
  ! - volume scattering function (delta)
  !
  ! This is useful for Code Verification
  ! of the Radiative Transfer solution algorithms
  ! via the Method of Manufactured Solutions.
  !
  ! NOTE: The manufactured solution must still meet the
  ! "no upwelling light from below" boundary condition.
  subroutine solve_rte_with_callbacks( &
       xmin, xmax, nx, &
       ymin, ymax, ny, &
       zmin, zmax, nz, &
       ntheta, nphi, &
       b, abs_func, source_func, source_expansion_func, bc_func, vsf_func, &
       radiance, irradiance, &
       num_scatters, num_threads, fd_flag, &
       lis_opts, lis_iter, lis_time, lis_resid)
    integer, intent(in) :: nx, ny, nz, ntheta, nphi
    double precision, intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax

    double precision, intent(in) :: b
    external :: abs_func, source_func, bc_func, vsf_func

    double precision, dimension(nx, ny, nz, ntheta*(nphi-2)+2), intent(inout) :: radiance
    double precision, dimension(nx, ny, nz), intent(inout) :: irradiance

    integer, intent(in) :: num_scatters
    integer, intent(in) :: num_threads
    logical, intent(in) :: fd_flag
    character*(*), intent(in) :: lis_opts

    integer, intent(inout) :: lis_iter
    double precision, intent(inout) :: lis_time, lis_resid

    ! End of procedure arguments

    type(space_angle_grid) grid
    type(rte_mat) mat
    type(optical_properties) iops
    type(light_state) light
    type(boundary_condition) bc
    integer i, j, k, p

    double precision dvsf

    double precision, dimension(:,:,:,:,:), allocatable :: source_expansion
    double precision, dimension(:,:,:,:), allocatable :: source

    ! Arrays for evaluating callbacks
    ! via vectorized numpy functions to improve speed
    double precision, dimension(:,:,:,:), allocatable :: x
    double precision, dimension(:,:,:,:), allocatable :: y
    double precision, dimension(:,:,:,:), allocatable :: z
    double precision, dimension(:,:,:,:), allocatable :: theta
    double precision, dimension(:,:,:,:), allocatable :: phi
    ! These temporary arrays are necessary because
    ! f2py can't infer the datatypes of members of derived types.
    ! and therefore expects the wrong kind of callbacks if they are used.
    double precision, dimension(:,:,:), allocatable :: x1
    double precision, dimension(:,:,:), allocatable :: y1
    double precision, dimension(:,:,:), allocatable :: z1
    double precision, dimension(:), allocatable :: theta1
    double precision, dimension(:), allocatable :: phi1
    double precision, dimension(:), allocatable :: tmp_vsf_cos, tmp_vsf_vals
    double precision, dimension(:,:,:), allocatable :: tmp_spatial
    double precision, dimension(:), allocatable :: tmp_angular
    integer num_vsf, nomega
    integer n

    ! The following line is an important f2py directive,
    ! not a comment.
    !f2py intent(out) tmp_vsf_vals, tmp_spatial, tmp_angular, source

    ! INIT GRID
    write(*,*) 'Grid'
    write(*,*) 'nx =', nx
    write(*,*) 'ny =', ny
    write(*,*) 'nz =', nz
    write(*,*) 'ntheta =', ntheta
    write(*,*) 'nphi =', nphi
    call grid%set_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
    call grid%set_num(nx, ny, nz, ntheta, nphi)
    call grid%init()

    allocate(source( &
         grid%x%num, &
         grid%y%num, &
         grid%z%num, &
         grid%angles%nomega))
    ! Last dimension iterates over expansion terms
    allocate(source_expansion( &
         grid%x%num, &
         grid%y%num, &
         grid%z%num, &
         grid%angles%nomega, &
         num_scatters+1))

    allocate(x(grid%x%num, 1, 1, 1))
    allocate(y(1, grid%y%num, 1, 1))
    allocate(z(1, 1, grid%z%num, 1))
    allocate(theta(1, 1, 1, grid%angles%nomega))
    allocate(phi(1, 1, 1, grid%angles%nomega))

    allocate(x1(grid%x%num, 1, 1))
    allocate(y1(1, grid%y%num, 1))
    allocate(z1(1, 1, grid%z%num))
    allocate(theta1(grid%angles%nomega))
    allocate(phi1(grid%angles%nomega))

    allocate(tmp_spatial(grid%x%num, grid%y%num, grid%z%num))
    allocate(tmp_angular(grid%angles%nomega))

    ! Create spatial/angular grid variables
    do i=1, grid%x%num
       x(i,:,:,:) = grid%x%vals(i)
       x1(i,:,:) = grid%x%vals(i)
    end do
    do j=1, grid%y%num
       y(:,j,:,:) = grid%y%vals(j)
       y1(:,j,:) = grid%y%vals(j)
    end do
    do k=1, grid%z%num
       z(:,:,k,:) = grid%z%vals(k)
       z1(:,:,k) = grid%z%vals(k)
    end do
    do p=1, grid%angles%nomega
       theta(:,:,:,p) = grid%angles%theta_p(p)
       phi(:,:,:,p) = grid%angles%phi_p(p)
       theta1(p) = grid%angles%theta_p(p)
       phi1(p) = grid%angles%phi_p(p)
    end do

    write(*,*) 'BC'
    ! Allocate bc_grid & initialize
    call bc%init(grid, 0.d0, 0.d0, 0.d0, 1.d0)

    ! INIT IOPS
    write(*,*) 'IOPs'

    iops%num_vsf = 101
    allocate(tmp_vsf_cos(iops%num_vsf))
    allocate(tmp_vsf_vals(iops%num_vsf))
    call iops%init(grid)
    ! Evaluate VSF on discrete grid
    ! Angles evenly spaced, endpoints included

    dvsf = pi/(dble(iops%num_vsf-1))
    do i=1, iops%num_vsf
       ! Defined on [0, pi]
       iops%vsf_angles(i) = dble(i-1)*dvsf
    end do
    ! Defined on [-1, 1]
    tmp_vsf_cos = cos(iops%vsf_angles)

    iops%scat = b
    num_vsf = iops%num_vsf

    write(*,*) 'calling vsf'
    write(*,*) 'num_vsf = ', num_vsf
    write(*,*) 'shape(tmp_vsf_cos) = ', shape(tmp_vsf_cos)
    write(*,*) 'shape(tmp_vsf_vals) = ', shape(tmp_vsf_vals)
    call vsf_func(tmp_vsf_cos, tmp_vsf_vals, num_vsf)
    write(*,*) 'vsf called'
    iops%vsf_vals = tmp_vsf_vals
    call iops%calc_vsf_on_grid()

    ! Calculate absorption coefficient
    ! and source term on discrete grids
    write(*,*) 'calling abs_func on x1, y1, z1, tmp_spatial'
    write(*,*) 'x1 = ', x1
    write(*,*) 'y1 = ', y1
    write(*,*) 'z1 = ', z1
    write(*,*) 'shape(tmp_spatial) = ', shape(tmp_spatial)
    call abs_func(x1, y1, z1, tmp_spatial, nx, ny, nz)
    write(*,*) 'abs called'
    iops%abs_grid = tmp_spatial
    write(*,*) 'assigned'

    nomega = grid%angles%nomega
    write(*,*) 'Evaluating source expansion'
    write(*,*) 'shape(x) =', shape(x)
    write(*,*) 'shape(y) =', shape(y)
    write(*,*) 'shape(z) =', shape(z)
    write(*,*) 'shape(theta) =', shape(theta)
    write(*,*) 'shape(phi) =', shape(phi)
    write(*,*) 'shape(source) =', shape(source)
    do n=0, num_scatters
       write(*,*) 'n =', n, '/', num_scatters
       call source_expansion_func(x, y, z, theta, phi, n, source, nx, ny, nz, nomega, num_scatters)
       write(*,*) 'source expansion called.'
       source_expansion(:,:,:,:,n+1) = source
       write(*,*) 'source expansion copied'
    end do
    write(*,*) 'Evaluating exact source'
    call source_func(x, y, z, theta, phi, n, source, nx, ny, nz, nomega)

    write(*,*) 'source called'
    call iops%set_source(source)
    write(*,*) 'source assigned'

    write(*,*) 'calling bc'
    call bc_func(theta1, phi1, tmp_angular, nomega)
    write(*,*) 'bc called'
    bc%bc_grid(1:grid%angles%nomega/2) = tmp_angular(1:grid%angles%nomega/2)
    write(*,*) 'bc assigned'

    write(*,*) 'min abs =', minval(iops%abs_grid)
    write(*,*) 'max abs =', maxval(iops%abs_grid)
    write(*,*) 'mean abs =', sum(iops%abs_grid)/size(iops%abs_grid)
    write(*,*) ''

    write(*,*) 'min source =', minval(source)
    write(*,*) 'max source =', maxval(source)
    write(*,*) 'mean source =', sum(source)/size(source)

    write(*,*) 'min bc =', minval(bc%bc_grid)
    write(*,*) 'max bc =', maxval(bc%bc_grid)
    write(*,*) 'mean bc =', sum(bc%bc_grid)/size(bc%bc_grid)

    write(*,*) 'Calculate asymptotic light field'
    call calculate_asymptotic_light_field_expanded_source(&
         grid, bc, iops, source, source_expansion, &
         radiance, num_scatters, num_threads)

    if(fd_flag) then

      ! INIT MAT
      write(*,*) 'Sparse Matrix'
      ! Set boundary condition
      call mat%init(grid, iops)
      call mat%set_bc(bc)
      call gen_matrix(mat, num_threads)

      ! Set solver options
      call mat%set_solver_opts(lis_opts)

      ! Initialize & set initial guess
      write(*,*) 'Light'
      call light%init(mat)
      light%radiance = radiance

      ! Solve system
      write(*,*) 'Calculate Radiance'
      call light%calculate_radiance()

      call mat%get_solver_stats(lis_iter, lis_time, lis_resid)
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

    deallocate(source_expansion)
    deallocate(source)
    deallocate(x)
    deallocate(y)
    deallocate(z)
    deallocate(theta)
    deallocate(phi)
    deallocate(tmp_vsf_cos)
    deallocate(tmp_vsf_vals)
    deallocate(tmp_spatial)
    deallocate(tmp_angular)

    write(*,*) 'done'
  end subroutine solve_rte_with_callbacks

end module pykelp3d_wrap
