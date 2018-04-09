module test_asymptotics
  use asymptotics
  use light_context
  implicit none

contains

  function test_max_cells(&
       xmin, xmax, ymin, ymax, zmin, zmax,&
       nx, ny, nz, ntheta, nphi) result(max_cells)
    type(space_angle_grid) grid
    integer max_cells
    double precision, intent(in) ::  xmin, xmax, ymin, ymax, zmin, zmax
    integer, intent(in) :: nx, ny, nz, ntheta, nphi

    call grid%set_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
    call grid%set_num(nx, ny, nz, ntheta, nphi)
    call grid%init()

    max_cells = calculate_max_cells(grid)
    call grid%deinit()
  end function test_max_cells

  subroutine test_traverse(&
       xmin, xmax, ymin, ymax, zmin, zmax,&
       nx, ny, nz, ntheta, nphi,&
       i, j, k, l, m,&
       s_array, ds, a_tilde, gn, rad_scatter,&
       num_cells, pkelp_fun)
    type(space_angle_grid) grid
    type(optical_properties) iops
    integer i, j, k, l, m, p
    integer ip, jp, kp
    integer num_cells
    double precision xmin, xmax, ymin, ymax, zmin, zmax
    double precision, dimension(:,:,:), allocatable :: p_kelp
    ! Fortran doesn't seem to enforce upper limits
    ! Couldn't find any reason not to do this
    double precision, dimension(1) :: s_array, ds, a_tilde, gn
    double precision, dimension(1,1,1,1) :: rad_scatter
    integer nx, ny, nz, ntheta, nphi, nomega
    double precision, allocatable, dimension(:,:,:,:) :: scatter_integral, source
    double precision, external :: pkelp_fun

    call grid%set_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
    call grid%set_num(nx, ny, nz, ntheta, nphi)
    call grid%init()
    nomega = grid%angles%nomega

    allocate(p_kelp(nx,ny,nz))
    allocate(scatter_integral(nx, ny, nz, nomega))
    allocate(source(nx,ny,nz,nomega))

    ! Just set vsf to be even, same angles as theta for convenience
    iops%num_vsf = ntheta
    call iops%init(grid)
    ! Copy arrays to avoid confusing allocation issues
    iops%vsf_angles(:) = grid%angles%theta(:)
    iops%vsf_vals(:) = 0*grid%angles%theta(:) + 1

    iops%abs_kelp = 1.0
    iops%scat = 0.0

    do kp=1, grid%z%num
      iops%abs_water(kp) = 0.0
    end do

    do ip=1, nx
       do jp=1, ny
          do kp=1, nz
             p_kelp(ip,jp,kp) = pkelp_fun(&
                  grid%x%vals(ip),&
                  grid%y%vals(jp),&
                  grid%z%vals(kp))
          end do
       end do
    end do

    call iops%calc_vsf_on_grid()
    call iops%calculate_coef_grids(p_kelp)

    p = grid%angles%phat(l,m)

    call calculate_source(grid, iops, rad_scatter, source, scatter_integral)

    call traverse_ray(grid, iops, source, i, j, k, p, s_array, ds, a_tilde, gn, num_cells)

    call iops%deinit()
    call grid%deinit()

    deallocate(p_kelp)
    deallocate(scatter_integral)
    deallocate(source)
  end subroutine test_traverse

  subroutine test_asymptotics_1d(I0, a, b, vsf_func, zmin, zmax, nz, z, Lp, Lm, num_scatters)
    !character(len=5), parameter :: fmtstr = 'E13.4'
    !character(len=56) :: vsf_file
    double precision, external :: vsf_func
    double precision, intent(in) :: I0, a, b, zmin, zmax
    integer, intent(in) :: nz
    integer, intent(in) :: num_scatters
    ! Downwelling (Lp) and upwelling (Lm) radiance
    double precision, intent(out), dimension(nz) :: Lp, Lm, z
    type(space_angle_grid) grid
    type(optical_properties) iops
    type(boundary_condition) bc
    type(light_state) light
    integer k
    double precision, dimension(1,1,nz) :: p_kelp

    call grid%set_bounds(0.d0, 1.d0, 0.d0, 1.d0, zmin, zmax)
    call grid%set_num(1, 1, nz, 1, 2)
    call grid%init()

    ! INIT IOPS
    iops%num_vsf = 55
    call iops%init(grid)
    iops%abs_kelp = a
    do k=1, nz
      iops%abs_water(k) = a
      p_kelp(1,1,k) = 0.d0
    end do

    iops%scat = b

    call iops%calculate_coef_grids(p_kelp)

    ! Will be called on cos vals [-1, 1]
    call iops%vsf_from_function(vsf_func)
    !vsf_file = '/home/oliver/academic/research/kelp/data/vsf/nuc_vsf.txt'
    !call iops%load_vsf(vsf_file, fmtstr)

    ! Straight downward light
    call bc%init(grid, 0.d0, 0.d0, 0.d0, I0)

    ! Rescale surface radiance to match surface irradiance
    bc%bc_grid = bc%bc_grid * I0 / grid%angles%integrate_points(bc%bc_grid)

    call light%init_grid(grid)

    call calculate_light_with_scattering(grid, bc, iops, light%radiance, num_scatters)

    ! Extract radiance
    Lp(:) = light%radiance(1,1,:,1)
    Lm(:) = light%radiance(1,1,:,2)

    ! Extract positions
    z(:) = grid%z%vals(:)

    call bc%deinit()
    call iops%deinit()
    call light%deinit()
    call grid%deinit()
  end subroutine test_asymptotics_1d

  subroutine test_asymptotics_3d(&
       I0, theta_s, phi_s, decay, &
       a_func, b, vsf_func, &
       nx, ny, nz, ntheta, nphi, &
       num_scatters, rad, irrad)

    double precision, external :: a_func, vsf_func
    double precision, intent(in) :: I0, b, theta_s, phi_s, decay
    integer, intent(in) :: nx, ny, nz, ntheta, nphi
    integer, intent(in) :: num_scatters
    double precision, dimension(nx, ny, nz, ntheta*(nphi-2)+2), intent(out) :: rad
    double precision, dimension(nx, ny, nz), intent(out) :: irrad

    type(space_angle_grid) grid
    type(optical_properties) iops
    type(boundary_condition) bc
    type(light_state) light

    double precision x, y, z
    integer i, j, k

    call grid%set_bounds(0.d0, 1.d0, 0.d0, 1.d0, 0.d0, 1.d0)
    call grid%set_num(nx, ny, nz, ntheta, nphi)
    call grid%init()

    ! INIT IOPS
    iops%num_vsf = 55
    call iops%init(grid)

    do i=1, nx
       x = grid%x%vals(i)
       do j=1, ny
          y = grid%y%vals(j)
          do k=1, nz
             z = grid%z%vals(k)
             iops%abs_grid(i,j,k) = a_func(x, y, z)
          end do
       end do
    end do

    iops%scat = b

    ! Will be called on cos vals [-1, 1]
    call iops%vsf_from_function(vsf_func)

    ! Straight downward light
    call bc%init(grid, theta_s, phi_s, decay, I0)

    ! Rescale surface radiance to match surface irradiance
    bc%bc_grid = bc%bc_grid * I0 / grid%angles%integrate_points(bc%bc_grid)

    call light%init_grid(grid)
    call calculate_light_with_scattering(grid, bc, iops, light%radiance, num_scatters)
    call light%calculate_irradiance()

    rad = light%radiance
    irrad = light%irradiance

    write(*,*) 'irrad'
    write(*,*) irrad

    call bc%deinit()
    call iops%deinit()
    call light%deinit()
    call grid%deinit()
  end subroutine test_asymptotics_3d

end module test_asymptotics
