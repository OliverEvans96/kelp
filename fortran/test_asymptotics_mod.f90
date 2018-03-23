module test_asymptotics_mod
  use asymptotics
  implicit none

contains

  subroutine test_traverse
    type(space_angle_grid) grid
    type(optical_properties) iops
    integer i, j, k, l, m, p
    integer max_cells
    double precision, dimension(:,:,:), allocatable :: p_kelp
    double precision, dimension(:), allocatable :: ds, f
    integer nx, ny, nz, ntheta, nphi
    logical passed

    nx = 10
    ny = 10
    nz = 10
    ntheta = 10
    nphi = 10

    allocate(p_kelp(nx,ny,nz))

    call grid%set_bounds(-1.d0,1.d0,-2.d0,2.d0,-3.d0,3.d0)
    call grid%set_num(nx, ny, nz, ntheta, nphi)
    call grid%init()

    call iops%init(grid)

    iops%abs_kelp = 1.0
    iops%scat_kelp = 0.0
    iops%abs_water = 1.0
    iops%scat_water = 0.0

    ! Just set vsf to be even, same angles as theta for convenience
    iops%num_vsf = ntheta
    iops%vsf_angles = grid%angles%theta
    iops%vsf_vals = 0*grid%angles%theta + 1

    do i=1, nx
       do j=1, ny
          do k=1, nz
             p_kelp(i,j,k) = 0.5d0
          end do
       end do
    end do

    call iops%calc_vsf_on_grid()
    call iops%calculate_coef_grids(p_kelp)

    i = 3
    j = 2
    k = 1
    l = 3
    m = 2

    p = grid%angles%phat(l,m)
    ! TODO: Formalize this into function in asymptotics
    ! This is definitely not right as is.
    max_cells = ceiling(grid%z%spacing(1)/grid%angles%cos_phi_p(p))

    allocate(ds(max_cells))
    allocate(f(max_cells))

    call traverse(grid, iops, i, j, k, p, ds, f)

    deallocate(ds)
    deallocate(f)
    deallocate(p_kelp)

    call iops%deinit()
    call grid%deinit()

    passed = .true.
    if(passed) then
       write(*,*) '* PASS - test_traverse'
    else
       write(*,*) '* FAIL - test_traverse'
    end if

  end subroutine test_traverse

end module test_asymptotics_mod
