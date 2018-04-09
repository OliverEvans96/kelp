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
   integer nx, ny, nz, nomega
   integer ent, i, j, k, p
   integer repeat_ent
   integer nonzero, n_total
   integer x_block_size, y_block_size, z_block_size, omega_block_size

   double precision, dimension(:), allocatable :: surface_vals

   ! A stored in coordinate form in row, col, data
   integer, dimension(:), allocatable :: row, col
   double precision, dimension(:), allocatable :: data
   ! b and x stored in rhs in full form
   double precision, dimension(:), allocatable :: rhs, sol

   ! Pointer to solver subroutine
   ! Set to mgmres by default
   procedure(solver_interface), pointer, nopass :: solver => mgmres_st

 contains
   procedure :: init => mat_init
   procedure :: deinit => mat_deinit
   procedure :: calculate_size
   procedure :: set_solver_params => mat_set_solver_params
   procedure :: set_row => mat_set_row
   procedure :: assign => mat_assign
   procedure :: add => mat_add
   procedure :: assign_rhs => mat_assign_rhs
   !procedure :: store_index => mat_store_index
   !procedure :: find_index => mat_find_index
   procedure :: set_bc => mat_set_bc
   procedure :: solve => mat_solve
   procedure :: ind => mat_ind
   procedure :: calculate_repeat_ent => mat_calculate_repeat_ent
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

abstract interface
   ! Define interface for external procedure
   ! https://stackoverflow.com/questions/8549415/how-to-declare-the-interface-section-for-a-procedure-argument-which-in-turn-ref
   subroutine solver_interface(n_total, nonzero, row, col, data, &
        sol, rhs, maxiter_outer, maxiter_inner, &
        tol_abs, tol_rel)
     integer ::  n_total, nonzero
     integer, dimension(nonzero) :: row, col
     double precision, dimension(nonzero) :: data
     double precision, dimension(nonzero) :: sol
     double precision, dimension(n_total) :: rhs
     integer :: maxiter_outer, maxiter_inner
     double precision :: tol_abs, tol_rel
   end subroutine solver_interface
end interface

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
    allocate(mat%surface_vals(grid%angles%nomega))
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
    integer nx, ny, nz, nomega

    nx = mat%grid%x%num
    ny = mat%grid%y%num
    nz = mat%grid%z%num
    nomega = mat%grid%angles%nomega

    !mat%nonzero = nx * ny * ntheta * nphi * ( (nz-1) * (6 + ntheta * nphi) + 1)
    mat%nonzero = nx * ny * nomega * (nz * (nomega + 6) - 1)
    mat%n_total = nx * ny * nz * nomega

    !mat%theta_block_size = 1
    !mat%phi_block_size = mat%theta_block_size * ntheta
    mat%omega_block_size = 1
    mat%y_block_size = mat%omega_block_size * nomega
    mat%x_block_size = mat%y_block_size * ny
    mat%z_block_size = mat%x_block_size * nx

  end subroutine calculate_size

!  subroutine mat_to_hdf(mat,filename)
!    class(rte_mat) mat
!    character(len=*) filename
!    call write_coo(filename, mat%row, mat%col, mat%data, mat%nonzero)
!  end subroutine mat_to_hdf

  subroutine mat_set_bc(mat, bc)
    class(rte_mat) mat
    class(boundary_condition) bc
    integer p

    do p=1, mat%grid%angles%nomega/2
       mat%surface_vals(p) = bc%bc_grid(p)
    end do
  end subroutine mat_set_bc

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
    open(unit=1, file='row.txt')
    open(unit=2, file='col.txt')
    open(unit=3, file='data.txt')
    open(unit=4, file='rhs.txt')
    write(1,*) mat%row
    write(2,*) mat%col
    write(3,*) mat%data
    write(4,*) mat%rhs
    close(1)
    close(2)
    close(3)
    close(4)

    open(unit=5, file='sol.txt')
    call mat%solver(mat%n_total, mat%nonzero, &
         mat%row, mat%col, mat%data, mat%sol, mat%rhs, &
         params%maxiter_outer, params%maxiter_inner, &
         params%tol_abs, params%tol_rel)

    write(5,*) mat%sol
    close(5)

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

  subroutine mat_calculate_repeat_ent(mat)
    ! Should be called when incrementing row
    class(rte_mat) mat

    ! Index of L_{ijklp}
    ! whose coefficient will be modified
    ! several times per row
    mat%repeat_ent = mat%ent + mat%p - 1

  end subroutine mat_calculate_repeat_ent

  function mat_ind(mat, i, j, k, p) result(ind)
    ! Assuming var ordering: z, x, y, omega
    class(rte_mat) mat
    type(space_angle_grid) grid
    integer i, j, k, p
    integer ind
    grid = mat%grid

    ind = (i-1) * mat%x_block_size + (j-1) * mat%y_block_size + &
         (k-1) * mat%z_block_size + p * mat%omega_block_size
  end function mat_ind

  subroutine mat_set_row(mat, indices)
    ! These indices act as a row counter
    ! Row should always be incremented in
    ! angular_integral, which should be called
    ! before derivatives, bcs, and attenuation
    class(rte_mat) mat
    type(index_list) indices

    mat%i = indices%i
    mat%j = indices%j
    mat%k = indices%k
    mat%p = indices%p

    call mat%calculate_repeat_ent()

  end subroutine mat_set_row

  subroutine mat_assign(mat, val, i, j, k, p)
    ! It's assumed that this is the only time this entry is defined
    class(rte_mat) mat
    double precision val
    integer i, j, k, p
    integer row_num, col_num

    row_num = mat%ind(mat%i, mat%j, mat%k,  mat%p)
    col_num = mat%ind(i, j, k, p)

    mat%row(mat%ent) = row_num
    mat%col(mat%ent) = col_num
    if(isnan(val)) then
       write(*,*) 'ISNAN'
       write(*,*) 'row = ', row_num
       write(*,*) 'col = ', col_num
       write(*,*) 'mat_index =', mat%i, mat%j, mat%k, mat%p
       write(*,*) 'index =', i, j, k, p
       write(*,*) 'entry =', mat%ent
    endif

    ! if(i.eq.mat%i .and. j.eq.mat%j .and. k.eq.mat%j .and. l.eq.mat%l .and. p.eq.mat%p) then
    !    write(*,*) 'diag: ', val
    ! endif

    mat%data(mat%ent) = val

    ! Remember where we stored information for this matrix element
    !call mat%store_index(row_num, col_num)

    mat%ent = mat%ent + 1
  end subroutine mat_assign

  subroutine mat_add(mat, val)
    ! Use this when you know that this entry has already been assigned
    ! and you'd like to add this value to the existing value.

    class(rte_mat) mat
    double precision val
    integer index

    ! Entry number where value is already stored
    index = mat%repeat_ent

    mat%data(index) = mat%data(index) + val
  end subroutine mat_add

  subroutine mat_assign_rhs(mat, data)
    class(rte_mat) mat
    double precision data
    integer row_num

    row_num = mat%ind(mat%i, mat%j, mat%k, mat%p)
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
    double precision aa, bb
    iops = mat%iops

    aa = iops%abs_grid(indices%i, indices%j, indices%k)
    bb = iops%scat

    attenuation = aa + bb*(1-iops%vsf_integral(indices%p, indices%p))
    call mat%add(attenuation)
  end subroutine attenuate

  subroutine x_cd2(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, dx
    type(index_list) indices
    integer i, j, k, p
    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p
    grid = mat%grid

    dx = grid%x%spacing(indices%k)

    val = grid%angles%sin_phi_p(p) * grid%angles%cos_theta_p(p) / (2.d0 * dx)

    call mat%assign(-val,i-1,j,k,p)
    call mat%assign(val,i+1,j,k,p)
  end subroutine x_cd2

  subroutine x_cd2_first(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, dx
    integer nx
    type(index_list) indices
    integer i, j, k, p

    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p
    grid = mat%grid

    dx = grid%x%spacing(indices%k)
    nx = grid%x%num

    val = grid%angles%sin_phi_p(p) * grid%angles%cos_theta_p(p) / (2.d0 * dx)

    call mat%assign(-val,nx,j,k,p)
    call mat%assign(val,i+1,j,k,p)
  end subroutine x_cd2_first

  subroutine x_cd2_last(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, dx
    type(index_list) indices
    integer i, j, k, p
    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p
    grid = mat%grid

    dx = grid%x%spacing(indices%k)

    val = grid%angles%sin_phi_p(p) * grid%angles%cos_theta_p(p) / (2.d0 * dx)

    call mat%assign(-val,i-1,j,k,p)
    call mat%assign(val,1,j,k,p)
  end subroutine x_cd2_last

  subroutine y_cd2(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, dy
    type(index_list) indices
    integer i, j, k, p
    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p
    grid = mat%grid

    dy = grid%y%spacing(indices%k)

    val = grid%angles%sin_phi_p(p) * grid%angles%sin_theta_p(p) / (2.d0 * dy)

    call mat%assign(-val,i,j-1,k,p)
    call mat%assign(val,i,j+1,k,p)
  end subroutine y_cd2

  subroutine y_cd2_first(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, dy
    integer ny
    type(index_list) indices
    integer i, j, k, p
    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p
    grid = mat%grid

    dy = grid%y%spacing(indices%k)
    ny = grid%y%num

    val = grid%angles%sin_phi_p(p) * grid%angles%sin_theta_p(p) / (2.d0 * dy)

    call mat%assign(-val,i,ny,k,p)
    call mat%assign(val,i,j+1,k,p)
  end subroutine y_cd2_first

  subroutine y_cd2_last(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, dy
    type(index_list) indices
    integer i, j, k, p
    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p
    grid = mat%grid

    dy = grid%y%spacing(indices%k)

    val = grid%angles%sin_phi_p(p) * grid%angles%sin_theta_p(p) / (2.d0 * dy)

    call mat%assign(-val,i,j-1,k,p)
    call mat%assign(val,i,1,k,p)
  end subroutine y_cd2_last

  subroutine z_cd2(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, dz
    type(index_list) indices
    integer i, j, k, p
    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p
    grid = mat%grid

    dz = grid%z%spacing(indices%k)

    val = grid%angles%cos_phi_p(p) / (2.d0 * dz)

    call mat%assign(-val,i,j,k-1,p)
    call mat%assign(val,i,j,k+1,p)
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
    integer i, j, k, p
    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p
    grid = mat%grid

    dz = grid%z%spacing(indices%k)

    val = grid%angles%cos_phi_p(p) / (2.d0 * dz)

    val1 = -3.d0 * val
    val2 = 4.d0 * val
    val3 = -val

    call mat%add(val1)
    call mat%assign(val2,i,j,k+1,p)
    call mat%assign(val3,i,j,k+2,p)
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
    integer i, j, k, p
    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p
    grid = mat%grid


    dz = grid%z%spacing(indices%k)

    val = grid%angles%cos_phi_p(p) / (2.d0 * dz)

    val1 = 3.d0 * val
    val2 = -4.d0 * val
    val3 = val

    call mat%add(val1)
    call mat%assign(val2,i,j,k-1,p)
    call mat%assign(val3,i,j,k-2,p)
  end subroutine z_bd2

  subroutine angular_integral(mat, indices)
    class(rte_mat) mat
    ! Primed angular integration variables
    integer pp
    double precision val
    type(index_list) indices

    call mat%set_row(indices)

    ! Interior
    do pp=1, mat%grid%angles%nomega
       ! TODO: Make sure I don't have p and pp backwards
       val = mat%iops%scat * mat%iops%vsf_integral(indices%p, pp)
       call mat%assign(val, indices%i, indices%j, indices%k, pp)
    end do

  end subroutine angular_integral

  subroutine z_surface_bc(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision bc_val
    type(index_list) indices
    double precision val1, val2

    grid = mat%grid

    val1 = grid%angles%cos_phi_p(indices%p) / (5.d0 * grid%z%spacing(1))
    val2 = 7.d0 * val1

    call mat%assign(val1,indices%i,indices%j,2,indices%p)
    call mat%add(val2)

    bc_val = 8.d0 * mat%surface_vals(indices%p) / (5.d0 * grid%z%spacing(1))
    call mat%assign_rhs(bc_val)

  end subroutine z_surface_bc

  subroutine z_bottom_bc(mat, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    type(index_list) indices
    double precision val1, val2

    grid = mat%grid

    val1 = -grid%angles%cos_phi_p(indices%p) / (5.d0 * grid%z%spacing(1))
    val2 = 7.d0 * val1

    call mat%assign(val1,indices%i,indices%j,grid%z%num-1,indices%p)
    call mat%add(val2)

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
