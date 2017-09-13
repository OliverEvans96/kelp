module rte_sparse_matrices
use sag
use kelp_context
use mgmres

type solver_params
   integer maxiter_inner, maxiter_outer
   double precision tol_abs, tol_rel
end type solver_params

type rte_mat
   type(space_angle_grid) grid
   type(optical_properties) iops
   type(solver_params) params
   integer nx, ny, nz, ntheta, nphi
   integer ent, i, j, k, l, m
   integer nonzero, n_total

   ! A stored in coordinate form in row, col, data
   integer, dimension(:), allocatable :: row, col
   double precision, dimension(:), allocatable :: data
   ! b and x stored in rhs in full form
   double precision, dimension(:), allocatable :: rhs, sol
 contains
   procedure :: init => mat_init
   procedure :: deinit => mat_deinit
   procedure calculate_size
   procedure set_solver_params
   procedure :: assign => mat_assign
   procedure :: solve => mat_solve
   procedure ind
   procedure attenuate
   procedure angular_integral

   ! Derivative subroutines
   procedure x_cd2
   procedure x_cd2_first
   procedure x_cd2_last
   procedure y_cd2
   procedure y_cd2_first
   procedure y_cd2_last
   procedure z_cd2
   procedure z_fd2
   procedure z_bd2
   procedure z_surface_bc
   procedure z_bottom_bc
end type rte_mat

contains

  subroutine mat_init(mat, grid, iops)
    class(rte_mat) mat
    type(space_angle_grid) grid
    type(optical_properties) iops
    integer nnz, n_total

    mat%grid = grid
    mat%iops = iops

    call mat%calculate_size()

    n_total = mat%n_total
    nnz = mat%nonzero
    allocate(mat%row(nnz))
    allocate(mat%col(nnz))
    allocate(mat%data(nnz))
    allocate(mat%rhs(n_total))
    allocate(mat%sol(n_total))

  end subroutine mat_init

  subroutine mat_deinit(mat)
    class(rte_mat) mat
    deallocate(mat%row)
    deallocate(mat%col)
    deallocate(mat%data)
    deallocate(mat%rhs)
    deallocate(mat%sol)
  end subroutine mat_deinit

  subroutine calculate_size(mat)
    class(rte_mat) mat
    type(space_angle_grid) grid
    integer nx, ny, nz, ntheta, nphi
    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    ntheta = grid%theta%num
    nphi = grid%phi%num

    mat%nonzero = nx * ny * (nz-2) * (6 + ntheta * nphi)
    mat%n_total = nx * ny * nz * ntheta * nphi
  
  end subroutine calculate_size

  subroutine mat_solve(mat)
    class(rte_mat) mat
    type(solver_params) params

    params = mat%params

    call mgmres_st(mat%n_total, mat%nonzero, mat%row, mat%col, mat%data, &
         mat%sol, mat%rhs, params%maxiter_outer, params%maxiter_inner, &
         params%tol_abs, params%tol_rel)

  end subroutine mat_solve

  subroutine set_solver_params(mat, maxiter_outer, &
       maxiter_inner, tol_abs, tol_rel)
    class(rte_mat) mat

    mat%params%maxiter_outer = maxiter_outer
    mat%params%maxiter_inner = maxiter_inner
    mat%params%tol_abs = tol_abs
    mat%params%tol_rel = tol_rel
  end subroutine set_solver_params

  function ind(mat, i, j, k, l, m)
    ! Assuming var ordering: z, y, x, phi, theta
    class(rte_mat) mat
    type(space_angle_grid) grid
    integer i, j, k, l, m, ind
    integer z_block_size, y_block_size, x_block_size
    integer phi_block_size, theta_block_size
    integer nx, ny, nz, ntheta, nphi
    grid = mat%grid

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    ntheta = grid%theta%num
    nphi = grid%phi%num

    theta_block_size = 1
    phi_block_size = theta_block_size * ntheta
    x_block_size = phi_block_size * nphi
    y_block_size = x_block_size * nx
    z_block_size = y_block_size * ny

    ind = (i-1) * x_block_size + (j-1) * y_block_size + (k-1) * z_block_size + (l) * theta_block_size + (m-1) * phi_block_size
  end function ind

  subroutine mat_assign(mat, data, i, j, k, l, m)
    class(rte_mat) mat
    double precision data
    integer i, j, k, l, m
    integer row_num, col_num

    row_num = mat%ind(mat%i, mat%j, mat%k, mat%l, mat%m)
    col_num = mat%ind(i, j, k, l, m)

    mat%row(mat%ent) = row_num
    mat%col(mat%ent) = col_num
    mat%data(mat%ent) = data

    mat%ent = mat%ent + 1
  end subroutine mat_assign

  subroutine attenuate(mat, indices)
    class(rte_mat) mat
    type(optical_properties) iops
    double precision attenuation
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m

    iops = mat%iops

    aa = iops%abs_grid(i, j, k)
    bb = iops%scat_grid(i, j, k)

    attenuation = aa + bb
    call mat%assign(attenuation, i, j, k, l, m)
  end subroutine attenuate

  subroutine x_cd2(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision theta, phi, dx
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    mat%grid = grid

    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dx = grid%x%spacing

    val = sin(phi) * cos(theta) / (2.d0 * dx)

    call mat%assign(-val,i-1,j,k,l,m)
    call mat%assign(val,i+1,j,k,l,m)
  end subroutine x_cd2
  
  subroutine x_cd2_first(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision theta, phi, dx
    integer nx
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    mat%grid = grid

    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dx = grid%x%spacing
    nx = grid%x%num

    val = sin(phi) * cos(theta) / (2.d0 * dx)

    call mat%assign(-val,nx,j,k,l,m)
    call mat%assign(val,i+1,j,k,l,m)
  end subroutine x_cd2_first
  
  subroutine x_cd2_last(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision theta, phi, dx
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    mat%grid = grid

    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dx = grid%x%spacing
    nx = grid%x%num

    val = sin(phi) * cos(theta) / (2.d0 * dx)

    call mat%assign(-val,i-1,j,k,l,m)
    call mat%assign(val,1,j,k,l,m)
  end subroutine x_cd2_last
  
  subroutine y_cd2(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision theta, phi, dy
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    mat%grid = grid


    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dy = grid%y%spacing

    val = sin(phi) * sin(theta) / (2.d0 * dy)

    call mat%assign(-val,i,j-1,k,l,m)
    call mat%assign(val,i,j+1,k,l,m)
  end subroutine y_cd2
  
  subroutine y_cd2_first(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision theta, phi, dy
    integer ny
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    mat%grid = grid


    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dy = grid%y%spacing
    ny = grid%y%num

    val = sin(phi) * sin(theta) / (2.d0 * dy)

    call mat%assign(-val,i,ny,k,l,m)
    call mat%assign(val,i,j+1,k,l,m)
  end subroutine y_cd2_first

  subroutine y_cd2_last(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision theta, phi, dy
    integer ny
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    mat%grid = grid


    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dy = grid%y%spacing

    val = sin(phi) * sin(theta) / (2.d0 * dy)

    call mat%assign(-val,i,j-1,k,l,m)
    call mat%assign(val,i,1,k,l,m)
  end subroutine y_cd2_last

  subroutine z_cd2(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision theta, phi, dz
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    mat%grid = grid


    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dz = grid%z%spacing

    val = cos(phi) / (2.d0 * dz)

    call mat%assign(-val,i,j,k-1,l,m)
    call mat%assign(val,i,j,k+1,l,m)
  end subroutine z_cd2
  
  subroutine z_fd2(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, val1, val2, val3
    double precision theta, phi, dz
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    mat%grid = grid


    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dz = grid%z%spacing

    val = cos(phi) / (2.d0 * dz)

    val1 = -3.d0 * val
    val2 = 2.d0 * val
    val3 = -val

    call mat%assign(val1,i,j,k,l,m)
    call mat%assign(val2,i,j,k+1,l,m)
    call mat%assign(val3,i,j,k+2,l,m)
  end subroutine z_fd2

  subroutine z_bd2(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, val1, val2, val3
    double precision theta, phi, dz
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    mat%grid = grid


    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dz = grid%z%spacing

    val = cos(phi) / (2.d0 * dz)

    val1 = 3.d0 * val
    val2 = -2.d0 * val
    val3 = val

    call mat%assign(val1,i,j,k,l,m)
    call mat%assign(val2,i,j,k-1,l,m)
    call mat%assign(val3,i,j,k-2,l,m)
  end subroutine z_bd2

  subroutine angular_integral(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    type(optical_properties) iops
    ! Primed integration variables
    integer lp, mp
    double precision val
    double precision prefactor
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid
    iops = mat%iops

    prefactor = grid%theta%prefactor * grid%phi%prefactor

    do lp=1, grid%theta%num
       do mp=1, grid%phi%num
          val = prefactor * iops%vsf(l,m,lp,mp)
       end do
    end do

  end subroutine angular_integral

  subroutine z_surface_bc(mat, indices)
    class(rte_mat) mat
    type(index_list) indices
  end subroutine z_surface_bc

  subroutine z_bottom_bc(mat, indices)
    class(rte_mat) mat
    type(index_list) indices
  end subroutine z_bottom_bc

  ! Finite difference wrappers

  subroutine wrap_x_cd2(mat, indices)
    type(rte_mat) mat
    type(index_list) indices
    call mat%x_cd2(indices)
  end subroutine wrap_x_cd2
  
  subroutine wrap_x_cd2_last(mat, indices)
    type(rte_mat) mat
    type(index_list) indices
    call mat%x_cd2_last(indices)
  end subroutine wrap_x_cd2_last

  subroutine wrap_x_cd2_first(mat, indices)
    type(rte_mat) mat
    type(index_list) indices
    call mat%x_cd2_first(indices)
  end subroutine wrap_x_cd2_first

  subroutine wrap_y_cd2(mat, indices)
    type(rte_mat) mat
    type(index_list) indices
    call mat%y_cd2(indices)
  end subroutine wrap_y_cd2

  subroutine wrap_y_cd2_last(mat, indices)
    type(rte_mat) mat
    type(index_list) indices
    call mat%y_cd2_last(indices)
  end subroutine wrap_y_cd2_last

  subroutine wrap_y_cd2_first(mat, indices)
    type(rte_mat) mat
    type(index_list) indices
    call mat%y_cd2_first(indices)
  end subroutine wrap_y_cd2_first

  subroutine wrap_z_cd2(mat, indices)
    type(rte_mat) mat
    type(index_list) indices
    call mat%z_cd2(indices)
  end subroutine wrap_z_cd2

end module rte_sparse_matrices

