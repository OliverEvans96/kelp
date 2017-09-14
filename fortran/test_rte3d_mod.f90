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
    type(optical_properties) iops

    write(*,*) 'Start'
    call gen_kelp(grid, rope, p_kelp)
    write(*,*) 'Test Domain:'
    write(*,*) 'nx =', grid%x%num
    write(*,*) 'ny =', grid%y%num
    write(*,*) 'nz =', grid%z%num
    write(*,*) 'ntheta =', grid%theta%num
    write(*,*) 'nphi =', grid%phi%num

    write(*,*) 'IOPs'
    call set_iops(iops, grid)
    write(*,*) 'Gen matrix'
    call gen_matrix(grid, mat, iops)

    write(*,*) 'Solver Params'
    call set_solver_params(mat)

    write(*,*) 'Solve'
    !call mat%solve()

    write(*,*) 'Kelp deinit'
    call kelp3d_deinit(grid, rope, p_kelp)
    write(*,*) 'RTE3D deinit'
    call rte3d_deinit(mat, iops)

  end subroutine run_test_rte3d

  subroutine set_iops(iops, grid)
    type(optical_properties) iops
    type(space_angle_grid) grid
    character(len=23) vsf_file
    character(len=6) vsf_fmt

    call iops%init(grid)

    ! Number of data points from Petzold
    iops%num_vsf = 55
    vsf_file = 'data/vsf/bahama_vsf.txt'
    vsf_fmt = 'ES13.5'
    call iops%load_vsf(vsf_file, vsf_fmt)

    iops%abs_kelp = 5
    iops%scat_kelp = 5

    iops%abs_water = 1
    iops%scat_water = 1

  end subroutine set_iops

  subroutine set_solver_params(mat)
    type(rte_mat) mat

    mat%params%maxiter_outer = 10
    mat%params%maxiter_inner = 10
    mat%params%tol_abs = 1.d-3
    mat%params%tol_rel = 1.d-3

  end subroutine set_solver_params

end module test_rte3d_mod
