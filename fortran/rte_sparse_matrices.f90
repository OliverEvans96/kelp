module rte_sparse_matrices
use sag
use kelp_context
use mgmres
implicit none

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
   procedure :: set_solver_params => mat_set_solver_params
   procedure :: set_ind => mat_set_ind
   procedure :: assign => mat_assign
   procedure :: add => mat_add
   procedure :: assign_rhs => mat_assign_rhs
   procedure :: find_index => mat_find_index
   procedure :: solve => mat_solve
   procedure :: ind => mat_ind
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

    call zeros(mat%rhs, n_total)
    call zeros(mat%sol, n_total)

  end subroutine mat_init

  subroutine mat_deinit(mat)
    class(rte_mat) mat
    write(*,*) 'md.1'
    write(*,*) allocated(mat%row)
    deallocate(mat%row)
    write(*,*) 'md.2'
    deallocate(mat%col)
    write(*,*) 'md.3'
    deallocate(mat%data)
    write(*,*) 'md.4'
    deallocate(mat%rhs)
    write(*,*) 'md.5'
    deallocate(mat%sol)
    write(*,*) 'md.6'
  end subroutine mat_deinit

  subroutine calculate_size(mat)
    class(rte_mat) mat
    integer nx, ny, nz, ntheta, nphi

    nx = mat%grid%x%num
    ny = mat%grid%y%num
    nz = mat%grid%z%num
    ntheta = mat%grid%theta%num
    nphi = mat%grid%phi%num

    mat%nonzero = 2 * nx * ny * ntheta * nphi * (nz-1) * (6 + ntheta * nphi)
    mat%n_total = nx * ny * nz * ntheta * nphi

  end subroutine calculate_size

  subroutine mat_solve(mat)
    class(rte_mat) mat
    type(solver_params) params

    params = mat%params

    write(*,*) 'mat%n_total =', mat%n_total
    write(*,*) 'mat%nonzero =', mat%nonzero
    write(*,*) 'size(mat%row) =', size(mat%row)
    write(*,*) 'size(mat%col) =', size(mat%col)
    write(*,*) 'size(mat%data) =', size(mat%data)
    write(*,*) 'size(mat%sol) =', size(mat%sol)
    write(*,*) 'size(mat%rhs) =', size(mat%rhs)
    write(*,*) 'params%maxiter_outer =', params%maxiter_outer
    write(*,*) 'params%maxiter_inner =', params%maxiter_inner
    write(*,*) 'params%tol_rel =', params%tol_rel
    write(*,*) 'params%tol_abs =', params%tol_abs

    call mgmres_st(mat%n_total, mat%nonzero, mat%row, mat%col, mat%data, &
         mat%sol, mat%rhs, params%maxiter_outer, params%maxiter_inner, &
         params%tol_abs, params%tol_rel)

  end subroutine mat_solve

  subroutine mat_set_solver_params(mat, maxiter_outer, &
       maxiter_inner, tol_abs, tol_rel)
    class(rte_mat) mat
    integer maxiter_outer, maxiter_inner
    double precision tol_abs, tol_rel

    mat%params%maxiter_outer = maxiter_outer
    mat%params%maxiter_inner = maxiter_inner
    mat%params%tol_abs = tol_abs
    mat%params%tol_rel = tol_rel
  end subroutine mat_set_solver_params

  function mat_ind(mat, i, j, k, l, m) result(ind)
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
  end function mat_ind

  subroutine mat_set_ind(mat, indices)
    ! These indices act as a row counter
    class(rte_mat) mat
    type(index_list) indices

    !write(*,*) 'Setting Mat Ind'
    !call indices%print()

    mat%i = indices%i
    mat%j = indices%j
    mat%k = indices%k
    mat%l = indices%l
    mat%m = indices%m
  end subroutine mat_set_ind

  subroutine mat_assign(mat, data, i, j, k, l, m)
    ! It's assumed that this is the only time this entry is defined
    class(rte_mat) mat
    double precision data
    integer i, j, k, l, m
    integer row_num, col_num

    row_num = mat%ind(mat%i, mat%j, mat%k, mat%l, mat%m)
    col_num = mat%ind(i, j, k, l, m)

    !write(*,*) 'inds =', i, j, k, l, m
    !write(*,*) 'rcd =', row_num, col_num, data
    !write(*,*) 'ent =', mat%ent

    if( ((row_num .gt. mat%n_total) .or. (col_num .gt. mat%nonzero)) &
         .or. ((row_num .lt. 1) .or. (col_num .lt. 1)) ) then
       write(*,*) 'OUT OF BOUNDS: ', row_num, col_num
       write(*,*) 'Using ROW i,j,k,l,m =', mat%i, mat%j, mat%k, mat%l, mat%m
       write(*,*) 'Using COL i,j,k,l,m =', i, j, k, l, m
       write(*,*)
       !!! END PROGRAM !!!
       write(*,*) 'Exiting'
       stop
       !write(*,*) 'whereas N =', mat%n_total
    end if

    mat%row(mat%ent) = row_num
    !write(*,*) 'a.1'
    mat%col(mat%ent) = col_num
    !write(*,*) 'a.2'
    mat%data(mat%ent) = data
    !write(*,*) 'a.3'

    mat%ent = mat%ent + 1
    !write(*,*) 'a.4'
  end subroutine mat_assign

  subroutine mat_add(mat, data, i, j, k, l, m)
    ! Use this when you know that this entry has already been assigned
    ! and you'd like to add this value to the existing value.
    class(rte_mat) mat
    double precision data
    integer i, j, k, l, m
    integer row_num, col_num
    integer index

    row_num = mat%ind(mat%i, mat%j, mat%k, mat%l, mat%m)
    col_num = mat%ind(i, j, k, l, m)

    index = mat%find_index(row_num, col_num)

    mat%data(index) = mat%data(index) + data
  end subroutine mat_add

  subroutine mat_assign_rhs(mat, data, i, j, k, l, m)
    class(rte_mat) mat
    double precision data
    integer i, j, k, l, m
    integer row_num

    row_num = mat%ind(mat%i, mat%j, mat%k, mat%l, mat%m)

    mat%rhs(row_num) = data
  end subroutine mat_assign_rhs

  function mat_find_index(mat, row_num, col_num) result(index)
    ! Find the position in row, col, data where this entry
    ! is defined.
    class(rte_mat) mat
    integer row_num, col_num, index

    ! Only search up to most recently assigned index
    do index=1, mat%ent-1
       if( (mat%row(index) .eq. row_num) .and. (mat%col(index) .eq. col_num)) then
          exit
       end if
    end do
  end function mat_find_index

  subroutine attenuate(mat, indices)
    ! Has to be called after angular_integral
    ! Because they both write to the same matrix entry
    ! And adding here is more efficient than a conditional
    ! in the angular loop.
    class(rte_mat) mat
    type(optical_properties) iops
    double precision attenuation
    type(index_list) indices
    integer i, j, k, l, m
    double precision aa, bb
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m

    iops = mat%iops

    aa = iops%abs_grid(i, j, k)
    bb = iops%scat_grid(i, j, k)

    attenuation = aa + bb
    call mat%add(attenuation, i, j, k, l, m)
  end subroutine attenuate

  subroutine x_cd2(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision sintheta, costheta, sinphi, cosphi, dx
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid
    write(*,*) 'x cd2'

    sintheta = grid%theta%sin(k)
    costheta = grid%theta%cos(k)
    sinphi = grid%phi%sin(l)
    cosphi = grid%phi%cos(l)
    dx = grid%x%spacing

    val = sinphi * costheta / (2.d0 * dx)

    write(*,*) 'VAL =', val

    write(*,*) '1'
    call mat%set_ind(indices)
    !call indices%print()
    call mat%assign(-val,i-1,j,k,l,m)
    write(*,*) '2'
    call mat%assign(val,i+1,j,k,l,m)
  end subroutine x_cd2
  
  subroutine x_cd2_first(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision sintheta, costheta, sinphi, cosphi, dx
    integer nx
    type(index_list) indices
    integer i, j, k, l, m
    write(*,*) 'x cd2 first'

    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid

    sintheta = grid%theta%sin(k)
    costheta = grid%theta%cos(k)
    sinphi = grid%phi%sin(l)
    cosphi = grid%phi%cos(l)
    dx = grid%x%spacing
    nx = grid%x%num

    val = sinphi * costheta / (2.d0 * dx)

    call mat%set_ind(indices)
    !call indices%print()
    call mat%assign(-val,nx,j,k,l,m)
    call mat%assign(val,i+1,j,k,l,m)
  end subroutine x_cd2_first
  
  subroutine x_cd2_last(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision sintheta, costheta, sinphi, cosphi, dx
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid
    write(*,*) 'x cd2 last'

    sintheta = grid%theta%sin(k)
    costheta = grid%theta%cos(k)
    sinphi = grid%phi%sin(l)
    cosphi = grid%phi%cos(l)
    dx = grid%x%spacing

    val = sinphi * costheta / (2.d0 * dx)

    call mat%set_ind(indices)
    !call indices%print()
    call mat%assign(-val,i-1,j,k,l,m)
    call mat%assign(val,1,j,k,l,m)
  end subroutine x_cd2_last
  
  subroutine y_cd2(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision sintheta, costheta, sinphi, cosphi, dy
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid

    write(*,*) 'ycd2'


    sintheta = grid%theta%sin(k)
    costheta = grid%theta%cos(k)
    sinphi = grid%phi%sin(l)
    cosphi = grid%phi%cos(l)
    dy = grid%y%spacing

    val = sinphi * sintheta / (2.d0 * dy)

    call mat%set_ind(indices)
    !call indices%print()
    call mat%assign(-val,i,j-1,k,l,m)
    call mat%assign(val,i,j+1,k,l,m)
  end subroutine y_cd2
  
  subroutine y_cd2_first(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision sintheta, costheta, sinphi, cosphi, dy
    integer ny
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid

    write(*,*) 'y cd2 first'

    sintheta = grid%theta%sin(k)
    costheta = grid%theta%cos(k)
    sinphi = grid%phi%sin(l)
    cosphi = grid%phi%cos(l)
    dy = grid%y%spacing
    ny = grid%y%num

    val = sinphi * sintheta / (2.d0 * dy)

    call mat%set_ind(indices)
    !call indices%print()
    call mat%assign(-val,i,ny,k,l,m)
    call mat%assign(val,i,j+1,k,l,m)
  end subroutine y_cd2_first

  subroutine y_cd2_last(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision sintheta, costheta, sinphi, cosphi, dy
    integer ny
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid
    write(*,*) 'y cd2 last'


    sintheta = grid%theta%sin(k)
    costheta = grid%theta%cos(k)
    sinphi = grid%phi%sin(l)
    cosphi = grid%phi%cos(l)
    dy = grid%y%spacing

    val = sinphi * sintheta / (2.d0 * dy)

    call mat%set_ind(indices)
    !call indices%print()
    call mat%assign(-val,i,j-1,k,l,m)
    call mat%assign(val,i,1,k,l,m)
  end subroutine y_cd2_last

  subroutine z_cd2(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision sintheta, costheta, sinphi, cosphi, dz
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid
    write(*,*) 'z cd2'

    sintheta = grid%theta%sin(k)
    costheta = grid%theta%cos(k)
    sinphi = grid%phi%sin(l)
    cosphi = grid%phi%cos(l)
    dz = grid%z%spacing

    val = cosphi / (2.d0 * dz)

    call mat%set_ind(indices)
    write(*,*) 'zcd2 ind:'
    !call indices%print()
    call mat%assign(-val,i,j,k-1,l,m)
    call mat%assign(val,i,j,k+1,l,m)
  end subroutine z_cd2
  
  subroutine z_fd2(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, val1, val2, val3
    double precision sintheta, costheta, sinphi, cosphi, dz
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid
    write(*,*) 'z fd2'


    sintheta = grid%theta%sin(k)
    costheta = grid%theta%cos(k)
    sinphi = grid%phi%sin(l)
    cosphi = grid%phi%cos(l)
    dz = grid%z%spacing

    val = cosphi / (2.d0 * dz)

    val1 = -3.d0 * val
    val2 = 2.d0 * val
    val3 = -val

    call mat%set_ind(indices)
    !call indices%print()
    call mat%assign(val1,i,j,k,l,m)
    call mat%assign(val2,i,j,k+1,l,m)
    call mat%assign(val3,i,j,k+2,l,m)
  end subroutine z_fd2

  subroutine z_bd2(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, val1, val2, val3
    double precision sintheta, costheta, sinphi, cosphi, dz
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid
    write(*,*) 'z bd2'

    sintheta = grid%theta%sin(k)
    costheta = grid%theta%cos(k)
    sinphi = grid%phi%sin(l)
    cosphi = grid%phi%cos(l)
    dz = grid%z%spacing

    val = cosphi / (2.d0 * dz)

    val1 = 3.d0 * val
    val2 = -2.d0 * val
    val3 = val

    call mat%set_ind(indices)
    write(*,*) 'Z BD2'
    call indices%print()
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
    call mat%set_ind(indices)
    !call indices%print()

    do lp=1, grid%theta%num
       do mp=1, grid%phi%num
          val = prefactor * iops%vsf(l,m,lp,mp)
          call mat%assign(val, i, j, k, lp, mp)
       end do
    end do

  end subroutine angular_integral

  subroutine z_surface_bc(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision bc_val
    double precision sintheta, costheta, sinphi, cosphi, dx
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m

    grid = mat%grid

    sintheta = grid%theta%sin(k)
    costheta = grid%theta%cos(k)
    sinphi = grid%phi%sin(l)
    cosphi = grid%phi%cos(l)
    dx = grid%x%spacing

    ! Constant light from above in all directions
    bc_val = 1.d0

    call mat%set_ind(indices)
    !call indices%print()
    call mat%assign(1.d0,i,j,k,l,m)
    call mat%assign_rhs(bc_val, i, j, k, l, m)
  end subroutine z_surface_bc

  subroutine z_bottom_bc(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision bc_val
    double precision sintheta, costheta, sinphi, cosphi, dx
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m

    grid = mat%grid

    sintheta = grid%theta%sin(k)
    costheta = grid%theta%cos(k)
    sinphi = grid%phi%sin(l)
    cosphi = grid%phi%cos(l)
    dx = grid%x%spacing

    ! No light from below
    bc_val = 0.d0

    call mat%set_ind(indices)
    !call indices%print()
    call mat%assign(1.d0,i,j,k,l,m)
    call mat%assign_rhs(bc_val, i, j, k, l, m)
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

