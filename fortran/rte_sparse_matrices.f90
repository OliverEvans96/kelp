module rte_sparse_matrices
use sag
use kelp_context
use mgmres
!use hdf5_utils
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
   integer repeat_index
   integer nonzero, n_total
   integer x_block_size, y_block_size, z_block_size, theta_block_size, phi_block_size

   double precision, dimension(:,:,:,:), allocatable :: surface_vals

   ! A stored in coordinate form in row, col, data
   integer, dimension(:), allocatable :: row, col
   double precision, dimension(:), allocatable :: data
   ! b and x stored in rhs in full form
   double precision, dimension(:), allocatable :: rhs, sol

   ! Given row and column number, determine corresponding position in row, col, data
   !integer, dimension(:,:), allocatable :: index_map
 contains
   procedure :: init => mat_init
   procedure :: deinit => mat_deinit
   procedure calculate_size
   procedure :: set_solver_params => mat_set_solver_params
   procedure :: set_ind => mat_set_ind
   procedure :: assign => mat_assign
   procedure :: add => mat_add
   procedure :: assign_rhs => mat_assign_rhs
   !procedure :: store_index => mat_store_index
   !procedure :: find_index => mat_find_index
   procedure :: initial_guess => mat_initial_guess
   procedure :: solve => mat_solve
   procedure :: ind => mat_ind
   procedure :: calculate_repeat_index => mat_calculate_repeat_index
   !procedure :: to_hdf => mat_to_hdf
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
    allocate(mat%surface_vals(grid%x%num, grid%y%num, grid%theta%num, grid%phi%num))
    allocate(mat%row(nnz))
    allocate(mat%col(nnz))
    allocate(mat%data(nnz))
    allocate(mat%rhs(n_total))
    allocate(mat%sol(n_total))

    call zeros(mat%rhs, n_total)
    call zeros(mat%sol, n_total)

    ! Start at first entry in row, col, data vectors
    mat%ent = 1

  end subroutine mat_init

  subroutine mat_deinit(mat)
    class(rte_mat) mat
    deallocate(mat%row)
    deallocate(mat%col)
    deallocate(mat%data)
    deallocate(mat%rhs)
    deallocate(mat%sol)
    deallocate(mat%surface_vals)
  end subroutine mat_deinit

  subroutine calculate_size(mat)
    class(rte_mat) mat
    integer nx, ny, nz, ntheta, nphi

    nx = mat%grid%x%num
    ny = mat%grid%y%num
    nz = mat%grid%z%num
    ntheta = mat%grid%theta%num
    nphi = mat%grid%phi%num

    mat%nonzero = nx * ny * ntheta * nphi * ( (nz-1) * (6 + ntheta * nphi) + 1)
    mat%n_total = nx * ny * nz * ntheta * nphi

    mat%theta_block_size = 1
    mat%phi_block_size = mat%theta_block_size * ntheta
    mat%x_block_size = mat%phi_block_size * nphi
    mat%y_block_size = mat%x_block_size * nx
    mat%z_block_size = mat%y_block_size * ny

  end subroutine calculate_size

!  subroutine mat_to_hdf(mat,filename)
!    class(rte_mat) mat
!    character(len=*) filename
!    call write_coo(filename, mat%row, mat%col, mat%data, mat%nonzero)
!  end subroutine mat_to_hdf

  subroutine mat_initial_guess(mat)
    class(rte_mat) mat
    integer i, j, k, l, m
    double precision dz
    double precision aa, bb, atten
    integer index

    ! Initial guess: Exponential decay downward using absorption coefficient

    index = 1
    dz = mat%grid%z%spacing

    do k=1, mat%grid%z%num
       do j=1, mat%grid%y%num
          do i=1, mat%grid%x%num
             ! Absorption coefficient
             aa = mat%iops%abs_grid(i,j,k)
             ! Scattering coefficient
             bb = mat%iops%scat_grid(i,j,k)
             ! Attenuation factor
             atten = 1.d0 - aa * dz
             ! Downwelling light
             do m=1, mat%grid%phi%num / 2
                do l=1, mat%grid%theta%num
                   mat%sol(index) = mat%surface_vals(i,j,l,m)
                   write(1,*) 'sol(', index, ') =', mat%sol(index)
                   index = index + 1
                   mat%surface_vals(i,j,l,m) = atten * mat%surface_vals(i,j,l,m)
                end do
             end do
             ! Upwelling light
             ! Still counting m from 1 since surface_vals is only defined up to nphi/2
             ! However index has incremented to the correct position for downwelling light.
             do m=1, mat%grid%phi%num / 2
                do l=1, mat%grid%theta%num
                   mat%sol(index) = mat%surface_vals(i,j,l,m)
                   write(1,*) 'sol(', index, ') =', mat%sol(index)
                   index = index + 1
                   mat%surface_vals(i,j,l,m) = bb * atten * mat%surface_vals(i,j,l,m)
                end do
             end do
          end do
       end do
    end do
  end subroutine mat_initial_guess


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
    write(1,*) mat%row
    write(2,*) mat%col
    write(3,*) mat%data
    write(4,*) mat%rhs

    call mat%initial_guess()
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

  subroutine mat_calculate_repeat_index(mat, indices)
    ! Must be called from angular loop
    ! Assuming outer loop is theta, inner is phi
    class(rte_mat) mat
    type(index_list) indices

    mat%repeat_index = mat%ent + mat%grid%phi%num * (indices%l-1) &
         + indices%m - 1
  end subroutine mat_calculate_repeat_index

  function mat_ind(mat, i, j, k, l, m) result(ind)
    ! Assuming var ordering: z, y, x, phi, theta
    class(rte_mat) mat
    type(space_angle_grid) grid
    integer i, j, k, l, m
    integer ind
    grid = mat%grid

    ind = (i-1) * mat%x_block_size + (j-1) * mat%y_block_size + &
         (k-1) * mat%z_block_size + (l) * mat%theta_block_size + &
         (m-1) * mat%phi_block_size
  end function mat_ind

  subroutine mat_set_ind(mat, indices)
    ! These indices act as a row counter
    class(rte_mat) mat
    type(index_list) indices

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

    mat%row(mat%ent) = row_num
    mat%col(mat%ent) = col_num
    mat%data(mat%ent) = data

    ! Remember where we stored information for this matrix element
    !call mat%store_index(row_num, col_num)

    mat%ent = mat%ent + 1
  end subroutine mat_assign

  subroutine mat_add(mat, data)
    ! Use this when you know that this entry has already been assigned
    ! and you'd like to add this value to the existing value.
    class(rte_mat) mat
    double precision data
    integer index

    index = mat%repeat_index

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

  ! subroutine mat_store_index(mat, row_num, col_num)
  !   ! Remember where we stored information for this matrix element
  !   class(rte_mat) mat
  !   integer row_num, col_num
  !   !mat%index_map(row_num, col_num) = mat%ent
  ! end subroutine

  ! function mat_find_index(mat, row_num, col_num) result(index)
  !   ! Find the position in row, col, data where this entry
  !   ! is defined.
  !   class(rte_mat) mat
  !   integer row_num, col_num, index

  !   index = mat%index_map(row_num, col_num)

  !   ! This took up 95% of execution time.
  !   ! Only search up to most recently assigned index
  !   ! do index=1, mat%ent-1
  !   !    if( (mat%row(index) .eq. row_num) .and. (mat%col(index) .eq. col_num)) then
  !   !       exit
  !   !    end if
  !   ! end do
  ! end function mat_find_index

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
    call mat%add(attenuation)
  end subroutine attenuate

  subroutine x_cd2(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, dx
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid

    dx = grid%x%spacing

    val = grid%phi%sin(m) * grid%theta%cos(l) / (2.d0 * dx)

    call mat%set_ind(indices)
    call mat%assign(-val,i-1,j,k,l,m)
    call mat%assign(val,i+1,j,k,l,m)
  end subroutine x_cd2
  
  subroutine x_cd2_first(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, dx
    integer nx
    type(index_list) indices
    integer i, j, k, l, m

    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid

    dx = grid%x%spacing
    nx = grid%x%num

    val = grid%phi%sin(m) * grid%theta%cos(l) / (2.d0 * dx)

    call mat%set_ind(indices)
    call mat%assign(-val,nx,j,k,l,m)
    call mat%assign(val,i+1,j,k,l,m)
  end subroutine x_cd2_first
  
  subroutine x_cd2_last(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, dx
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid

    dx = grid%x%spacing

    val = grid%phi%sin(m) * grid%theta%cos(l) / (2.d0 * dx)

    call mat%set_ind(indices)
    call mat%assign(-val,i-1,j,k,l,m)
    call mat%assign(val,1,j,k,l,m)
  end subroutine x_cd2_last
  
  subroutine y_cd2(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, dy
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid

    dy = grid%y%spacing

    val = grid%phi%sin(m) * grid%theta%sin(l) / (2.d0 * dy)

    call mat%set_ind(indices)
    call mat%assign(-val,i,j-1,k,l,m)
    call mat%assign(val,i,j+1,k,l,m)
  end subroutine y_cd2
  
  subroutine y_cd2_first(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, dy
    integer ny
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid

    dy = grid%y%spacing
    ny = grid%y%num

    val = grid%phi%sin(m) * grid%theta%sin(l) / (2.d0 * dy)

    call mat%set_ind(indices)
    call mat%assign(-val,i,ny,k,l,m)
    call mat%assign(val,i,j+1,k,l,m)
  end subroutine y_cd2_first

  subroutine y_cd2_last(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, dy
    integer ny
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid

    dy = grid%y%spacing

    val = grid%phi%sin(m) * grid%theta%sin(l) / (2.d0 * dy)

    call mat%set_ind(indices)
    call mat%assign(-val,i,j-1,k,l,m)
    call mat%assign(val,i,1,k,l,m)
  end subroutine y_cd2_last

  subroutine z_cd2(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, dz
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid

    dz = grid%z%spacing

    val = grid%phi%cos(m) / (2.d0 * dz)

    call mat%set_ind(indices)
    call mat%assign(-val,i,j,k-1,l,m)
    call mat%assign(val,i,j,k+1,l,m)
  end subroutine z_cd2
  
  subroutine z_fd2(mat, indices)
    ! Has to be called after angular_integral
    ! Because they both write to the same matrix entry
    ! And adding here is more efficient than a conditional
    ! in the angular loop.
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, val1, val2, val3, dz
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid

    dz = grid%z%spacing

    val = grid%phi%cos(m) / (2.d0 * dz)

    val1 = -3.d0 * val
    val2 = 2.d0 * val
    val3 = -val

    call mat%set_ind(indices)
    call mat%add(val1)
    call mat%assign(val2,i,j,k+1,l,m)
    call mat%assign(val3,i,j,k+2,l,m)
  end subroutine z_fd2

  subroutine z_bd2(mat, indices)
    ! Has to be called after angular_integral
    ! Because they both write to the same matrix entry
    ! And adding here is more efficient than a conditional
    ! in the angular loop.
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, val1, val2, val3, dz
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m
    grid = mat%grid

    dz = grid%z%spacing

    val = grid%phi%cos(m) / (2.d0 * dz)

    val1 = 3.d0 * val
    val2 = -2.d0 * val
    val3 = val

    call mat%set_ind(indices)
    call mat%add(val1)
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

    ! Store entry number of element at i,j,k,l,m
    ! to allow other functions to add to it
    call mat%calculate_repeat_index(indices)

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
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m

    grid = mat%grid

    ! Constant light from above in all directions
    bc_val = 1.d1

    mat%surface_vals(i,j,l,m) = bc_val

    call mat%set_ind(indices)
    call mat%assign(1.d0,i,j,k,l,m)
    call mat%assign_rhs(bc_val, i, j, k, l, m)
  end subroutine z_surface_bc

  subroutine z_bottom_bc(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision bc_val
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m

    grid = mat%grid

    ! No light from below
    bc_val = 0.d0

    call mat%set_ind(indices)
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

