module test_asymptotics
  use asymptotics
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
    double precision, allocatable, dimension(:) :: scatter_integrand
    double precision, external :: pkelp_fun

    call grid%set_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
    call grid%set_num(nx, ny, nz, ntheta, nphi)
    call grid%init()
    nomega = grid%angles%nomega

    ! write(*,*) 'nx = ', nx
    ! write(*,*) 'ny = ', ny
    ! write(*,*) 'nz = ', nz
    ! write(*,*) 'nomega = ', nomega

    ! write(*,*) 'x =', grid%x%vals
    ! write(*,*) 'y =', grid%y%vals
    ! write(*,*) 'z =', grid%z%vals
    ! write(*,*) 'theta =', grid%angles%theta
    ! write(*,*) 'phi =', grid%angles%phi

    allocate(p_kelp(nx,ny,nz))
    allocate(scatter_integral(nx, ny, nz, nomega))
    allocate(scatter_integrand(nomega))
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
    !write(*,*) 'PHAT l, m, p =', l, m, p

    call calculate_source(grid, iops, rad_scatter, source, scatter_integral, scatter_integrand)

    call traverse_ray(grid, iops, source, i, j, k, p, s_array, ds, a_tilde, gn, num_cells, .false.)

    call iops%deinit()
    call grid%deinit()

    deallocate(p_kelp)
    deallocate(scatter_integral)
    deallocate(scatter_integrand)
    deallocate(source)
  end subroutine test_traverse

  subroutine test_asymptotics_1d(I0, a, b, vsf_func, zmin, zmax, nz, z, Lp, Lm, num_scatters)
    !character(len=5), parameter :: fmtstr = 'E13.4'
    !character(len=56) :: vsf_file
    double precision, external :: vsf_func
    double precision, intent(in) :: I0, a, b, zmin, zmax
    integer, intent(in) :: nz, num_scatters
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
    write(*,*) 'IOPs'
    iops%abs_kelp = a
    do k=1, nz
      iops%abs_water(k) = a
      p_kelp(1,1,k) = 0.d0
    end do

    iops%scat = b

    call iops%calculate_coef_grids(p_kelp)

    !write(*,*) 'iop init'
    !iops%vsf_angles = vsf_angles
    !iops%vsf_vals = vsf_vals

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
    !write(*,*) 'a'
    call iops%deinit()
    !write(*,*) 'b'
    call light%deinit()
    !write(*,*) 'c'
    call grid%deinit()
    !write(*,*) 'e'
  end subroutine test_asymptotics_1d

end module test_asymptotics
