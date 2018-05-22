module rte_sparse_matrices
use sag
use kelp_context
use mgmres
!use hdf5_utils
#include "lisf.h"
implicit none

type solver_opts
   integer maxiter_inner, maxiter_outer
   double precision tol_abs, tol_rel
end type solver_opts

type rte_mat
   type(space_angle_grid) grid
   type(optical_properties) iops
   type(solver_opts) params
   integer nx, ny, nz, nomega
   integer i, j, k, p
   integer nonzero, n_total
   integer x_block_size, y_block_size, z_block_size, omega_block_size

   double precision, dimension(:), allocatable :: surface_vals

   ! A stored in coordinate form in row, col, data
   integer, dimension(:), allocatable :: row, col
   double precision, dimension(:), allocatable :: data

   ! Lis Matrix and vectors
   LIS_MATRIX A
   LIS_VECTOR b, x
   LIS_SOLVER solver
   LIS_INTEGER ierr
   character(len=256) solver_opts

   ! Pointer to solver subroutine
   ! Set to mgmres by default
   !procedure(solver_interface), pointer, nopass :: solver => mgmres_st

 contains
   procedure :: init => mat_init
   procedure :: deinit => mat_deinit
   procedure :: calculate_size
   procedure :: set_solver_opts => mat_set_solver_opts
   procedure :: assign => mat_assign
   procedure :: add => mat_add
   procedure :: assign_rhs => mat_assign_rhs
   !procedure :: store_index => mat_store_index
   !procedure :: find_index => mat_find_index
   procedure :: set_bc => mat_set_bc
   procedure :: solve => mat_solve
   procedure :: ind => mat_ind
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

interface
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

    LIS_INTEGER comm_world

    comm_world = LIS_COMM_WORLD

    mat%grid = grid
    mat%iops = iops

    call mat%calculate_size()

    mat%solver_opts = ''
    mat%ierr = 0

    n_total = mat%n_total
    nnz = mat%nonzero

    call lis_initialize(mat%ierr)

    call lis_solver_create(mat%solver, mat%ierr)

    call lis_matrix_create(comm_world, mat%A, mat%ierr)
    call lis_vector_create(comm_world, mat%b, mat%ierr)
    call lis_vector_create(comm_world, mat%x, mat%ierr)

    call lis_matrix_set_size(mat%A, 0, n_total, mat%ierr)
    call lis_vector_set_size(mat%b, 0, n_total, mat%ierr)
    call lis_vector_set_size(mat%x, 0, n_total, mat%ierr)

    call lis_vector_set_all(0.0d0, mat%x, mat%ierr)
    call lis_vector_set_all(0.0d0, mat%b, mat%ierr)

    if(mat%ierr .ne. 0) then
       write(*,*) 'INIT ERR: ', mat%ierr
       call exit(1)
    end if

    allocate(mat%row(nnz))
    allocate(mat%col(nnz))
    allocate(mat%data(nnz))
    allocate(mat%surface_vals(grid%angles%nomega))
  end subroutine mat_init

  subroutine mat_deinit(mat)
    class(rte_mat) mat

    call lis_matrix_destroy(mat%A, mat%ierr)
    call lis_vector_destroy(mat%b, mat%ierr)
    call lis_vector_destroy(mat%x, mat%ierr)
    call lis_solver_destroy(mat%solver, mat%ierr)
    call lis_finalize(mat%ierr)

    if(mat%ierr .ne. 0) then
       write(*,*) 'DEINIT ERR: ', mat%ierr
       call exit(1)
    end if

    deallocate(mat%row)
    deallocate(mat%col)
    deallocate(mat%data)
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
    character(len=24) init_opt
    LIS_VECTOR res

    ! type(solver_opts) params

    !params = mat%params

    write(*,*) 'mat%n_total =', mat%n_total
    write(*,*) 'mat%nonzero =', mat%nonzero
    ! write(*,*) 'size(mat%row) =', size(mat%row)
    ! write(*,*) 'size(mat%col) =', size(mat%col)
    ! write(*,*) 'size(mat%data) =', size(mat%data)
    ! write(*,*) 'size(mat%sol) =', size(mat%sol)
    ! write(*,*) 'size(mat%rhs) =', size(mat%rhs)
    ! write(*,*) 'params%maxiter_outer =', params%maxiter_outer
    ! write(*,*) 'params%maxiter_inner =', params%maxiter_inner
    ! write(*,*) 'params%tol_rel =', params%tol_rel
    ! write(*,*) 'params%tol_abs =', params%tol_abs
    open(unit=1, file='row.txt')
    open(unit=2, file='col.txt')
    open(unit=3, file='data.txt')
    !open(unit=4, file='rhs.txt')
    open(unit=5, file='sol.txt')
    write(1,*) mat%row
    write(2,*) mat%col
    write(3,*) mat%data
    !write(4,*) mat%rhs

    close(1)
    close(2)
    close(3)
    !close(4)

    ! call mat%solver(mat%n_total, mat%nonzero, &
    !      mat%row, mat%col, mat%data, mat%sol, mat%rhs, &
    !      params%maxiter_outer, params%maxiter_inner, &
    !      params%tol_abs, params%tol_rel)

    call lis_matrix_set_coo(mat%nonzero, mat%row, mat%col, mat%data, &
         mat%A, mat%ierr)
    call lis_matrix_assemble(mat%A, mat%ierr)

    write(*,*) 'MM'
    call lis_matvec(mat%A, mat%x, res, mat%ierr)
    write(*,*) 'Q'
    call lis_vector_print(res, mat%ierr)

    write(*,*) 'row =', mat%row(1:10)
    write(*,*) 'col =', mat%col(1:10)
    write(*,*) 'data =', mat%data(1:10)
    write(*,*) 'b'
    call lis_vector_print(mat%b, mat%ierr)
    write(*,*) 'x'
    call lis_vector_print(mat%x, mat%ierr)

    write(*,*) ''

    if(len(trim(mat%solver_opts)) .gt. 0) then
       call lis_solver_set_option(mat%solver_opts, mat%solver, mat%ierr)
    end if
    init_opt = "-initx_zeros false"
    call lis_solver_set_option(init_opt, mat%solver, mat%ierr)
    call lis_solve(mat%A, mat%b, mat%x, mat%solver, mat%ierr)

    ! write(5,*) mat%sol
    ! close(5)

  end subroutine mat_solve

  subroutine mat_set_solver_opts(mat, solver_opts)
    class(rte_mat) mat
    character(len=*) solver_opts
    mat%solver_opts = solver_opts
  end subroutine mat_set_solver_opts

  function mat_ind(mat, i, j, k, p) result(ind)
    ! Assuming var ordering: z, x, y, omega
    class(rte_mat) mat
    integer i, j, k, p
    integer ind

    ind = (i-1) * mat%x_block_size + (j-1) * mat%y_block_size + &
         (k-1) * mat%z_block_size + p * mat%omega_block_size
  end function mat_ind

  subroutine mat_assign(mat, row_num, ent, val, i, j, k, p)
    ! It's assumed that this is the only time this entry is defined
    class(rte_mat) mat
    double precision val
    integer i, j, k, p
    integer row_num, ent

    mat%row(ent) = row_num
    mat%col(ent) = mat%ind(i, j, k, p)
    mat%data(ent) = val

    if(ent .eq. 1) then
       write(*,*) 'rcd =', mat%row(ent), mat%col(ent), mat%data(ent)
       write(*,*) ':)a'
    end if

    ent = ent + 1
  end subroutine mat_assign

  subroutine mat_add(mat, repeat_ent, val)
    ! Use this when you know that this entry has already been assigned
    ! and you'd like to add this value to the existing value.

    class(rte_mat) mat
    double precision val
    integer repeat_ent

    ! Entry number where value is already stored
    mat%data(repeat_ent) = mat%data(repeat_ent) + val

    if(repeat_ent .eq. 1) then
       write(*,*) 'rcvd =', mat%row(1), mat%col(1), val, mat%data(1)
       write(*,*) ':)b'
    end if

  end subroutine mat_add

  subroutine mat_assign_rhs(mat, row_num, data)
    class(rte_mat) mat
    double precision data
    integer row_num

    call lis_vector_set_value(LIS_INS_VALUE, row_num, data, mat%b, mat%ierr)
    if(mat%ierr .ne. 0) then
       write(*,*) 'RHS ERR: ', mat%ierr
       call exit(1)
    end if

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

  subroutine attenuate(mat, indices, repeat_ent)
    ! Has to be called after angular_integral
    ! Because they both write to the same matrix entry
    ! And adding here is more efficient than a conditional
    ! in the angular loop.
    class(rte_mat) mat
    double precision attenuation
    type(index_list) indices
    double precision aa, bb
    integer repeat_ent

    aa = mat%iops%abs_grid(indices%i, indices%j, indices%k)
    bb = mat%iops%scat
    attenuation = aa + bb

    call mat%add(repeat_ent, attenuation)
  end subroutine attenuate

  subroutine x_cd2(mat, indices, row_num, ent)
    class(rte_mat) mat
    double precision val, dx
    type(index_list) indices
    integer i, j, k, p
    integer row_num, ent

    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p

    dx = mat%grid%x%spacing(1)

    val = mat%grid%angles%sin_phi_p(p) &
         * mat%grid%angles%cos_theta_p(p) / (2.d0 * dx)

    call mat%assign(row_num,ent,-val,i-1,j,k,p)
    call mat%assign(row_num,ent,val,i+1,j,k,p)
  end subroutine x_cd2

  subroutine x_cd2_first(mat, indices, row_num, ent)
    class(rte_mat) mat
    double precision val, dx
    integer nx
    type(index_list) indices
    integer i, j, k, p
    integer row_num, ent

    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p

    dx = mat%grid%x%spacing(1)
    nx = mat%grid%x%num

    val = mat%grid%angles%sin_phi_p(p) &
         * mat%grid%angles%cos_theta_p(p) / (2.d0 * dx)

    call mat%assign(row_num,ent,-val,nx,j,k,p)
    call mat%assign(row_num,ent,val,i+1,j,k,p)
  end subroutine x_cd2_first

  subroutine x_cd2_last(mat, indices, row_num, ent)
    class(rte_mat) mat
    double precision val, dx
    type(index_list) indices
    integer i, j, k, p
    integer row_num, ent

    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p

    dx = mat%grid%x%spacing(1)

    val = mat%grid%angles%sin_phi_p(p) &
         * mat%grid%angles%cos_theta_p(p) / (2.d0 * dx)

    call mat%assign(row_num,ent,-val,i-1,j,k,p)
    call mat%assign(row_num,ent,val,1,j,k,p)
  end subroutine x_cd2_last

  subroutine y_cd2(mat, indices, row_num, ent)
    class(rte_mat) mat
    double precision val, dy
    type(index_list) indices
    integer i, j, k, p
    integer row_num, ent

    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p

    dy = mat%grid%y%spacing(1)

    val = mat%grid%angles%sin_phi_p(p) &
         * mat%grid%angles%sin_theta_p(p) / (2.d0 * dy)

    call mat%assign(row_num,ent,-val,i,j-1,k,p)
    call mat%assign(row_num,ent,val,i,j+1,k,p)
  end subroutine y_cd2

  subroutine y_cd2_first(mat, indices, row_num, ent)
    class(rte_mat) mat
    double precision val, dy
    integer ny
    type(index_list) indices
    integer i, j, k, p
    integer row_num, ent

    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p

    dy = mat%grid%y%spacing(1)
    ny = mat%grid%y%num

    val = mat%grid%angles%sin_phi_p(p) &
         * mat%grid%angles%sin_theta_p(p) / (2.d0 * dy)

    call mat%assign(row_num,ent,-val,i,ny,k,p)
    call mat%assign(row_num,ent,val,i,j+1,k,p)
  end subroutine y_cd2_first

  subroutine y_cd2_last(mat, indices, row_num, ent)
    class(rte_mat) mat
    double precision val, dy
    type(index_list) indices
    integer i, j, k, p
    integer row_num, ent

    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p

    dy = mat%grid%y%spacing(1)

    val = mat%grid%angles%sin_phi_p(p) &
         * mat%grid%angles%sin_theta_p(p) / (2.d0 * dy)

    call mat%assign(row_num,ent,-val,i,j-1,k,p)
    call mat%assign(row_num,ent,val,i,1,k,p)
  end subroutine y_cd2_last

  subroutine z_cd2(mat, indices, row_num, ent)
    class(rte_mat) mat
    double precision val, dz
    type(index_list) indices
    integer i, j, k, p
    integer row_num, ent

    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p

    dz = mat%grid%z%spacing(indices%k)

    val = mat%grid%angles%cos_phi_p(p) / (2.d0 * dz)

    call mat%assign(row_num,ent,-val,i,j,k-1,p)
    call mat%assign(row_num,ent,val,i,j,k+1,p)
  end subroutine z_cd2

  subroutine z_fd2(mat, indices, row_num, ent, repeat_ent)
    ! Has to be called after angular_integral
    ! Because they both write to the same matrix entry
    ! And adding here is more efficient than a conditional
    ! in the angular loop.
    class(rte_mat) mat
    double precision val, val1, val2, val3, dz
    type(index_list) indices
    integer i, j, k, p
    integer row_num, ent, repeat_ent

    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p

    dz = mat%grid%z%spacing(indices%k)

    val = mat%grid%angles%cos_phi_p(p) / (2.d0 * dz)

    val1 = -3.d0 * val
    val2 = 4.d0 * val
    val3 = -val

    call mat%add(repeat_ent, val1)
    call mat%assign(row_num,ent,val2,i,j,k+1,p)
    call mat%assign(row_num,ent,val3,i,j,k+2,p)
  end subroutine z_fd2

  subroutine z_bd2(mat, indices, row_num, ent, repeat_ent)
    ! Has to be called after angular_integral
    ! Because they both write to the same matrix entry
    ! And adding here is more efficient than a conditional
    ! in the angular loop.
    class(rte_mat) mat
    double precision val, val1, val2, val3, dz
    type(index_list) indices
    integer i, j, k, p
    integer row_num, ent, repeat_ent

    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p

    dz = mat%grid%z%spacing(indices%k)

    val = mat%grid%angles%cos_phi_p(p) / (2.d0 * dz)

    val1 = 3.d0 * val
    val2 = -4.d0 * val
    val3 = val

    call mat%add(repeat_ent, val1)
    call mat%assign(row_num,ent,val2,i,j,k-1,p)
    call mat%assign(row_num,ent,val3,i,j,k-2,p)
  end subroutine z_bd2

  subroutine angular_integral(mat, indices, row_num, ent)
    class(rte_mat) mat
    ! Primed angular integration variables
    integer pp
    double precision val
    type(index_list) indices
    integer row_num, ent

    ! Interior
    do pp=1, mat%grid%angles%nomega
       ! TODO: Make sure I don't have p and pp backwards
       val = -mat%iops%scat * mat%iops%vsf_integral(indices%p, pp)
       call mat%assign(row_num, ent, val, indices%i, indices%j, indices%k, pp)
    end do
  end subroutine angular_integral

  subroutine z_surface_bc(mat, indices, row_num, ent, repeat_ent)
    class(rte_mat) mat
    double precision bc_val
    type(index_list) indices
    double precision val1, val2, dz
    integer row_num, ent, repeat_ent

    dz = mat%grid%z%spacing(1)

    val1 = mat%grid%angles%cos_phi_p(indices%p) / (5.d0 * dz)
    val2 = 7.d0 * val1
    bc_val = 8.d0 * val1 * mat%surface_vals(indices%p)

    call mat%assign(row_num,ent,val1,indices%i,indices%j,2,indices%p)
    call mat%add(repeat_ent, val2)
    call mat%assign_rhs(row_num, bc_val)

  end subroutine z_surface_bc

    subroutine z_bottom_bc(mat, indices, row_num, ent, repeat_ent)
    class(rte_mat) mat
    type(index_list) indices
    double precision val1, val2, dz
    integer nz
    integer row_num, ent, repeat_ent

    dz = mat%grid%z%spacing(1)
    nz = mat%grid%z%num

    val1 = -mat%grid%angles%cos_phi_p(indices%p) / (5.d0 * dz)
    val2 = 7.d0 * val1

    call mat%assign(row_num,ent,val1,indices%i,indices%j,nz-1,indices%p)
    call mat%add(repeat_ent, val2)

  end subroutine z_bottom_bc

  ! Finite difference wrappers

  ! subroutine wrap_x_cd2(mat, indices)
  !   type(rte_mat) mat
  !   type(index_list) indices
  !   call mat%x_cd2(indices)
  ! end subroutine wrap_x_cd2

  ! subroutine wrap_x_cd2_last(mat, indices)
  !   type(rte_mat) mat
  !   type(index_list) indices
  !   call mat%x_cd2_last(indices)
  ! end subroutine wrap_x_cd2_last

  ! subroutine wrap_x_cd2_first(mat, indices)
  !   type(rte_mat) mat
  !   type(index_list) indices
  !   call mat%x_cd2_first(indices)
  ! end subroutine wrap_x_cd2_first

  ! subroutine wrap_y_cd2(mat, indices)
  !   type(rte_mat) mat
  !   type(index_list) indices
  !   call mat%y_cd2(indices)
  ! end subroutine wrap_y_cd2

  ! subroutine wrap_y_cd2_last(mat, indices)
  !   type(rte_mat) mat
  !   type(index_list) indices
  !   call mat%y_cd2_last(indices)
  ! end subroutine wrap_y_cd2_last

  ! subroutine wrap_y_cd2_first(mat, indices)
  !   type(rte_mat) mat
  !   type(index_list) indices
  !   call mat%y_cd2_first(indices)
  ! end subroutine wrap_y_cd2_first

  ! subroutine wrap_z_cd2(mat, indices)
  !   type(rte_mat) mat
  !   type(index_list) indices
  !   call mat%z_cd2(indices)
  ! end subroutine wrap_z_cd2

end module rte_sparse_matrices
