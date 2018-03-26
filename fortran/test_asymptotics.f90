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
       s_array, ds, a_tilde, gn, rad_scatter, num_cells)
    type(space_angle_grid) grid
    type(optical_properties) iops
    integer i, j, k, l, m, p
    integer ip, jp, kp
    integer num_cells
    double precision xmin, xmax, ymin, ymax, zmin, zmax
    double precision, dimension(:,:,:), allocatable :: p_kelp
    ! Fortran doesn't seem to enforce upper limits
    ! Couldn't find any reason not to do this
    double precision, dimension(:) :: s_array, ds, a_tilde, gn
    double precision, dimension(:,:,:,:) :: rad_scatter
    integer nx, ny, nz, ntheta, nphi, nomega
    double precision, allocatable, dimension(:,:,:,:) :: scatter_integral, source
    double precision, allocatable, dimension(:) :: scatter_integrand

    call grid%set_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
    call grid%set_num(nx, ny, nz, ntheta, nphi)
    call grid%init()
    nomega = grid%angles%nomega

    write(*,*) 'nx = ', nx
    write(*,*) 'ny = ', ny
    write(*,*) 'nz = ', nz
    write(*,*) 'nomega = ', nomega

    write(*,*) 'x =', grid%x%vals
    write(*,*) 'y =', grid%y%vals
    write(*,*) 'z =', grid%z%vals
    write(*,*) 'theta =', grid%angles%theta
    write(*,*) 'phi =', grid%angles%phi

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

    iops%scat_kelp = 0.0
    iops%abs_kelp = 1.0

    do kp=1, grid%z%num
      iops%abs_water(kp) = 1.0
      iops%scat_water(kp) = 0.0
    end do

    do ip=1, nx
       do jp=1, ny
          do kp=1, nz
             p_kelp(ip,jp,kp) = 0.5d0
          end do
       end do
    end do

    call iops%calc_vsf_on_grid()
    call iops%calculate_coef_grids(p_kelp)

    p = grid%angles%phat(l,m)
    write(*,*) 'PHAT l, m, p =', l, m, p

    call calculate_source(grid, iops, rad_scatter, source, scatter_integral, scatter_integrand)

    call traverse_ray(grid, iops, source, i, j, k, p, s_array, ds, a_tilde, gn, num_cells)

    call iops%deinit()
    call grid%deinit()

    deallocate(p_kelp)
    deallocate(scatter_integral)
    deallocate(scatter_integrand)
    deallocate(source)
  end subroutine test_traverse

end module test_asymptotics
