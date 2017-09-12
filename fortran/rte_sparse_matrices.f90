module rte_sparse_matrices
use sag
use kelp_context

type rte_mat
   type(space_angle_grid) grid
   integer nx, ny, nz, ntheta, nphi
   integer ent, i, j, k, l, m

   ! A stored in coordinate form in row, col, data
   integer, dimension(:), allocatable :: row, col
   double precision, dimension(:), allocatable :: data
   ! b stored in rhs in full form
   double precision, dimension(:), allocatable :: rhs
 contains
   procedure init, deinit
   procedure :: assign => mat_assign
   procedure nonzero
   procedure ind

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

  subroutine init(mat)
    class(rte_mat) mat
    type(space_angle_grid) grid
    integer nnz

    mat%grid = grid
    nnz = mat%nonzero()
    allocate(mat%row(nnz))
    allocate(mat%col(nnz))
    allocate(mat%data(nnz))
  end subroutine init

  subroutine deinit(mat)
    class(rte_mat) mat
    deallocate(mat%row)
    deallocate(mat%col)
    deallocate(mat%data)
  end subroutine deinit

  function nonzero(mat)
    class(rte_mat) mat
    double precision nonzero
    !!! *** !!!
    nonzero = 0
  end function nonzero

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

  subroutine attenuate(mat, iops, indices)
    class(rte_mat) mat
    type(space_angle_grid) grid
    type(optical_properties) iops
    double precision attenuation
    type(index_list) indices
    integer i, j, k, l, m
    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m

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

  subroutine angular_integral(mat, iops, indices)
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
    mat%grid = grid

    prefactor = grid%theta%prefactor * grid%phi%prefactor

    do lp=1, grid%theta%num
       do mp=1, grid%phi%num
          val = prefactor * iops%vsf(l,m,lp,mp)
       end do
    end do

  end subroutine angular_integral

  subroutine z_surface_bc(mat)
    class(rte_mat) mat
  end subroutine z_surface_bc

  subroutine z_bottom_bc(mat)
    class(rte_mat) mat
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

