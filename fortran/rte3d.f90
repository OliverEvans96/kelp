module rte3d
use kelp_context
implicit none
contains

  subroutine gen_matrix(grid, mat)
    type(space_angle_grid) grid
    type(rte_matrix) mat
    integer nonzero
    double precision phi, sinphi, theta, sintheta, costheta 
    integer i, j, k, l, m

    nonzero = mat%nonzero()

    mat%init(grid)

    do k=2 , grid%z%num - 1
       do i=1, grid%x%num
          do j=1, grid%y%num
             do m=1, grid%phi%num
                phi = grid%phi%vals(m)
                sinphi = sin(phi)
                do l=1, grid%theta%num
                   call mat%attenuate(grid)
                   call mat%x_cd2(grid)
                   call mat%y_cd2(grid)
                   call mat%z_cd2(grid)
                   call mat%angular_integral(grid)
                end do
             end do
          end do
       end do
    end do

    mat%deinit()

end subroutine gen_matrix

subroutine matrix_deinit()

end module rte3d
