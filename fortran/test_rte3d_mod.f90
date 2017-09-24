module test_rte3d_mod
  use test_kelp3d_mod
  use rte3d
  implicit none
contains

  subroutine run_test_rte3d()
    double precision, dimension(:,:,:), allocatable :: p_kelp
    type(space_angle_grid) grid
    type(rope_state) rope
    type(rte_mat) mat
    type(light_state) light
    type(optical_properties) iops

    character(len=21), parameter :: matfile = 'hdf5/kelp3d/mat.hdf5'

    write(*,*) 'Start'
    call gen_kelp(grid, rope, p_kelp)
    write(*,*) 'Test Domain:'
    write(*,*) 'nx =', grid%x%num
    write(*,*) 'ny =', grid%y%num
    write(*,*) 'nz =', grid%z%num
    write(*,*) 'ntheta =', grid%theta%num
    write(*,*) 'nphi =', grid%phi%num

    write(*,*) 'IOPs'
    call set_iops(iops, grid, p_kelp)
    write(*,*) 'Gen matrix'
    call gen_matrix(grid, mat, iops)

    call light%init(mat)

    write(*,*) 'Solver Params'
    call set_solver_params(mat)

    write(*,*) 'ent =', mat%ent, '/', mat%nonzero

    write(*,*) 'Write'
    call mat%to_hdf(matfile)

    write(*,*) 'Radiance'
    call light%calculate_radiance()

    write(*,*) 'Irradiance'
    call light%calculate_irradiance()

    write(*,*) ' deinit'
    call kelp3d_deinit(grid, rope, p_kelp)
    write(*,*) 'RTE3D deinit'
    call rte3d_deinit(mat, iops, light)

  end subroutine run_test_rte3d

  subroutine set_iops(iops, grid, p_kelp)
    type(optical_properties) iops
    type(space_angle_grid) grid
    double precision, dimension(:,:,:) :: p_kelp
    character(len=23) vsf_file
    character(len=6) vsf_fmt

    call iops%init(grid)

    ! Number of data points from Petzold
    iops%num_vsf = 55
    vsf_file = 'data/vsf/bahama_vsf.txt'
    vsf_fmt = 'ES13.5'

    iops%abs_kelp = 5
    iops%scat_kelp = 5

    iops%abs_water = 1
    iops%scat_water = 1

    call iops%load_vsf(vsf_file, vsf_fmt)
    call iops%calculate_coef_grids(p_kelp)

  end subroutine set_iops

  subroutine set_solver_params(mat)
    type(rte_mat) mat

    mat%params%maxiter_outer = 10
    mat%params%maxiter_inner = 10
    mat%params%tol_abs = 1.d-3
    mat%params%tol_rel = 1.d-3

  end subroutine set_solver_params

end module test_rte3d_mod
