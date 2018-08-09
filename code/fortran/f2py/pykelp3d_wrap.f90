module pykelp3d_wrap
  use rte3d
  use kelp3d
  use asymptotics
  use light_interface
  implicit none

  interface
     subroutine func1d(delta, output)
       double precision, intent(in) :: delta
       double precision, intent(out) :: output
     end subroutine func1d

     subroutine func2d(theta, phi, output)
       double precision, intent(in) :: theta, phi
       double precision, intent(out) :: output
     end subroutine func2d

     subroutine func3d(x, y, z, output)
       double precision, intent(in) :: x, y, z
       double precision, intent(out) :: output
     end subroutine func3d

     subroutine func5d(x, y, z, theta, phi, output)
       double precision, intent(in) :: x, y, z, theta, phi
       double precision, intent(out) :: output
     end subroutine func5d
  end interface

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
       theta_s, phi_s, I0, decay, &
       p_kelp, radiance, irradiance, avg_irrad, perceived_irrad, &
       num_scatters, fd_flag, lis_opts, &
       lis_iter, lis_time, lis_resid)

    integer, intent(in) :: nx, ny, nz, ntheta, nphi
    double precision, intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax
    double precision, intent(in) :: a_w, a_k, b
    integer, intent(in) ::num_vsf
    double precision, dimension(num_vsf), intent(in) :: vsf_angles, vsf_vals
    double precision, intent(in) :: theta_s, phi_s, I0, decay
    double precision, dimension(nx, ny, nz), intent(in) :: p_kelp
    double precision, dimension(nx, ny, nz, ntheta*(nphi-2)+2), intent(inout) :: radiance
    double precision, dimension(nx, ny, nz), intent(inout) :: irradiance
    ! This is sloppy, but these functions need reals.
    real, dimension(nz), intent(inout) :: avg_irrad, perceived_irrad

    integer, intent(in) :: num_scatters
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
    !write(*,*) 'full abs:'
    !write(*,*) iops%abs_grid

    write(*,*) 'BC'
    call bc%init(grid, theta_s, phi_s, decay, I0)

    !write(*,*) 'bc_grid = ', bc%bc_grid

    write(*,*) 'Calculate asymptotic light field'
    call calculate_asymptotic_light_field(grid, bc, iops, source, radiance, num_scatters)

    if(fd_flag) then

      ! INIT MAT
      write(*,*) 'Sparse Matrix'
      ! Set boundary condition
      call mat%init(grid, iops)
      call mat%set_bc(bc)
      call gen_matrix(mat)

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
       b, &
       !abs_func, source_func, bc_func, vsf_func, &
       radiance, irradiance, avg_irrad, perceived_irrad, &
       num_scatters, fd_flag, lis_opts, &
       lis_iter, lis_time, lis_resid)
    integer, intent(in) :: nx, ny, nz, ntheta, nphi
    double precision, intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax

    double precision, intent(in) :: b
    ! procedure(func3d) :: abs_func
    ! procedure(func5d) :: source_func
    ! procedure(func2d) :: bc_func
    ! procedure(func1d) :: vsf_func

    double precision, dimension(nx, ny, nz, ntheta*(nphi-2)+2), intent(inout) :: radiance
    double precision, dimension(nx, ny, nz), intent(inout) :: irradiance
    ! This is sloppy, but these functions need reals.
    real, dimension(nz), intent(inout) :: avg_irrad, perceived_irrad

    integer, intent(in) :: num_scatters
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
    double precision x, y, z, theta, phi
    double precision delta, tmp

    double precision dvsf

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

    ! INIT IOPS
    write(*,*) 'IOPs'

    iops%num_vsf = 101
    call iops%init(grid)
    ! Evaluate VSF on discrete grid
    ! Angles evenly spaced, endpoints included
    dvsf = pi/(dble(iops%num_vsf-1))
    do i=1, iops%num_vsf
       iops%vsf_angles(i) = dble(i-1)*dvsf
       delta = iops%vsf_angles(i)
       ! call vsf_func(delta, tmp)
       ! iops%vsf_vals(i) = tmp
    end do
    call iops%calc_vsf_on_grid()

    iops%scat = b

    ! Calculate absorption coefficient
    ! and source term on discrete grids
    do i=1, grid%x%num
       x = grid%x%vals(i)
       do j=1, grid%y%num
          y = grid%x%vals(j)
          do k=1, grid%z%num
             z = grid%x%vals(k)
             ! call abs_func(x, y, z, tmp)
             ! iops%abs_grid(i,j,k) = tmp

             do p=1, grid%angles%nomega
                theta = grid%angles%theta_p(p)
                phi = grid%angles%phi_p(p)
                ! call source_func(x, y, z, theta, phi, tmp)
                ! source(i,j,k,p) = tmp
             end do
          end do
       end do
    end do

    write(*,*) 'BC'
    ! Allocate bc_grid & initialize
    call bc%init(grid, 0.d0, 0.d0, 0.d0, 1.d0)
    ! Apply correct BC
    do p=1, grid%angles%nomega/2
       theta = grid%angles%theta_p(p)
       phi = grid%angles%phi_p(p)
       ! call bc_func(theta, phi, tmp)
       ! bc%bc_grid(p) = tmp
    end do

    write(*,*) 'Calculate asymptotic light field'
    call calculate_asymptotic_light_field(grid, bc, iops, source, radiance, num_scatters)

    if(fd_flag) then

      ! INIT MAT
      write(*,*) 'Sparse Matrix'
      ! Set boundary condition
      call mat%init(grid, iops)
      call mat%set_bc(bc)
      call gen_matrix(mat)

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
  end subroutine solve_rte_with_callbacks

end module pykelp3d_wrap
