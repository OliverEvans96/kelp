module rte_sparse_matrices
use sag

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
   procedure assign => mat_assign
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

  subroutine init(mat, grid)
    class(rte_matrix) mat
    type(space_angle_grid) grid
    integer nnz

    nnz = mat%nonzero()
    allocate(mat%row(nnz))
    allocate(mat%col(nnz))
    allocate(mat%data(nnz))
  end subroutine init

  subroutine deinit(mat)
    class(rte_matrix) mat
    deallocate(mat%row)
    deallocate(mat%col)
    deallocate(mat%data)
  end subroutine deinit

  function nonzero(mat)
    !!! *** !!!
  end function nonzero

  function ind(mat, grid, i, j, k, l, m)
    ! Assuming var ordering: z, y, x, phi, theta
    class(rte_mat) mat
    type(space_angle_grid) grid
    integer i, j, k, l, m, ind
    integer z_block_size, y_block_size, x_block_size,
    integer phi_block_size, theta_block_size

    integer nx, ny, nz, ntheta, nphi
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
  end subroutine assign

  subroutine attenuate(mat, iops)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision attenuation
    attenuation = aa(mat%i, mat%j, mat%k) + bb(mat%i, mat%j, mat%k)
    mat%assign(attenuation, i, j, k, l, m)
  end subroutine attenuate

  subroutine x_cd2(mat, grid)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision theta, phi, dx

    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dx = grid%x%spacing

    val = sin(phi) * cos(theta) / (2.d0 * dx)

    mat%assign(-val,i-1,j,k,l,m)
    mat%assign(val,i+1,j,k,l,m)
  end subroutine x_cd2
  
  subroutine x_cd2_first(mat, grid)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision theta, phi, dx
    integer nx

    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dx = grid%x%spacing
    nx = grid%x%num

    val = sin(phi) * cos(theta) / (2.d0 * dx)

    mat%assign(-val,nx,j,k,l,m)
    mat%assign(val,i+1,j,k,l,m)
  end subroutine x_cd2
  
  subroutine x_cd2_last(mat, grid)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision theta, phi, dx

    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dx = grid%x%spacing
    nx = grid%x%num

    val = sin(phi) * cos(theta) / (2.d0 * dx)

    mat%assign(-val,i-1,j,k,l,m)
    mat%assign(val,1,j,k,l,m)
  end subroutine x_cd2_first
  
  subroutine y_cd2(mat, grid)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision theta, phi, dy

    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dy = grid%y%spacing

    val = sin(phi) * sin(theta) / (2.d0 * dy)

    mat%assign(-val,i,j-1,k,l,m)
    mat%assign(val,i,j+1,k,l,m)
  end subroutine y_cd2
  
  subroutine y_cd2_first(mat, grid)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision theta, phi, dy
    integer ny

    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dy = grid%y%spacing
    ny = grid%y%num

    val = sin(phi) * sin(theta) / (2.d0 * dy)

    mat%assign(-val,i,ny,k,l,m)
    mat%assign(val,i,j+1,k,l,m)
  end subroutine y_cd2_first

  subroutine y_cd2_last(mat, grid)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision theta, phi, dy
    integer ny

    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dy = grid%y%spacing

    val = sin(phi) * sin(theta) / (2.d0 * dy)

    mat%assign(-val,i,j-1,k,l,m)
    mat%assign(val,i,1,k,l,m)
  end subroutine y_cd2_last

  subroutine z_cd2(mat, grid)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val
    double precision theta, phi, dz

    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dz = grid%z%spacing

    val = cos(phi) / (2.d0 * dz)

    mat%assign(-val,i,j,k-1,l,m)
    mat%assign(val,i,j,k+1,l,m)
  end subroutine z_cd2
  
  subroutine z_fd2(mat, grid)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, val1, val2, val3
    double precision theta, phi, dz

    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dz = grid%z%spacing

    val = cos(phi) / (2.d0 * dz)

    val1 = -3.d0 * val
    val2 = 2.d0 * val
    val3 = -val

    mat%assign(val1,i,j,k,l,m)
    mat%assign(val2,i,j,k+1,l,m)
    mat%assign(val3,i,j,k+2,l,m)
  end subroutine z_fd2

  subroutine z_bd2(mat, grid)
    class(rte_mat) mat
    type(space_angle_grid) grid
    double precision val, val1, val2, val3
    double precision theta, phi, dz

    theta = grid%theta%vals(mat%k)
    phi = grid%theta%vals(mat%l)
    dz = grid%z%spacing

    val = cos(phi) / (2.d0 * dz)

    val1 = 3.d0 * val
    val2 = -2.d0 * val
    val3 = val

    mat%assign(val1,i,j,k,l,m)
    mat%assign(val2,i,j,k-1,l,m)
    mat%assign(val3,i,j,k-2,l,m)
  end subroutine z_bd2

  subroutine angular_integral(mat, grid, iops)
    class(rte_mat) mat
    type(space_angle_grid) grid
    type(optical_properties) iops
    ! Primed integration variables
    integer lp, mp

    double precision prefactor
    double precision angdiff

    do lp=1, grid%theta%num
       do mp=1, grid%phi%num
          prefactor * iops%vsf(angdiff)
       end do
    end do

  end subroutine angular_integral

end module
