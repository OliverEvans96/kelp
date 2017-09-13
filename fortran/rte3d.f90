module rte3d
use kelp_context
use rte_sparse_matrices
implicit none
contains

subroutine interior_space_loop(mat, indices)
  type(rte_mat) mat
  type(space_angle_grid) grid
  type(index_list) indices
  integer i, j, k

  ! z interior
  do k=2 , grid%z%num - 1
     ! x first
     indices%i=1
       ! y first
       indices%j=1
       call interior_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2_first)
       ! y interior
       do j=2, grid%y%num - 1
       indices%j = j
           call interior_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2)
       end do
       ! y last
       indices%j=grid%y%num
       call interior_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2_last)

     ! x interior
     do i=2, grid%x%num - 1
     indices%i = i
        ! y first
        indices%j=1
        call interior_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2_first)
        ! y interior
        do j=2, grid%y%num - 1
        indices%j = j
           call interior_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2)
        end do
        ! y last
        indices%j=grid%y%num
        call interior_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2_last)
     end do

     ! x last
     indices%i=grid%x%num
       ! y first
       indices%j=1
       call interior_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2_first)
       ! y interior
       do j=2, grid%y%num - 1
       indices%j = j
          call interior_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2)
       end do
       ! y last
       indices%j=grid%y%num
       call interior_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2_last)

  end do

end subroutine


subroutine surface_space_loop(mat, indices)
  type(rte_mat) mat
  type(space_angle_grid) grid
  type(index_list) indices
  integer i, j, k

  ! z surface
  indices%k=1
     ! x first
     indices%i=1
       ! y first
       indices%j=1
       call surface_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2_first)
       ! y surface
       do j=2, grid%y%num - 1
       indices%j = j
           call surface_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2)
       end do
       ! y last
       indices%j=grid%y%num
       call surface_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2_last)

     ! x surface
     do i=2, grid%x%num - 1
     indices%i = i
        ! y first
        indices%j=1
        call surface_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2_first)
        ! y surface
        do j=2, grid%y%num - 1
        indices%j = j
           call surface_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2)
        end do
        ! y last
        indices%j=grid%y%num
        call surface_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2_last)
     end do

     ! x last
     indices%i=grid%x%num
       ! y first
       indices%j=1
       call surface_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2_first)
       ! y surface
       do j=2, grid%y%num - 1
       indices%j = j
          call surface_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2)
       end do
       ! y last
       indices%j=grid%y%num
       call surface_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2_last)

end subroutine surface_space_loop

subroutine bottom_space_loop(mat, indices)
  type(rte_mat) mat
  type(space_angle_grid) grid
  type(index_list) indices
  integer i, j, k

  ! z bottom
  indices%k=grid%z%num
     ! x first
     indices%i=1
       ! y first
       indices%j=1
       call bottom_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2_first)
       ! y bottom
       do j=2, grid%y%num - 1
       indices%j = j
           call bottom_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2)
       end do
       ! y last
       indices%j=grid%y%num
       call bottom_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2_last)

     ! x bottom
     do i=2, grid%x%num - 1
     indices%i = i
        ! y first
        indices%j=1
        call bottom_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2_first)
        ! y bottom
        do j=2, grid%y%num - 1
        indices%j = j
           call bottom_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2)
        end do
        ! y last
        indices%j=grid%y%num
        call bottom_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2_last)
     end do

     ! x last
     indices%i=grid%x%num
       ! y first
       indices%j=1
       call bottom_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2_first)
       ! y bottom
       do j=2, grid%y%num - 1
       indices%j = j
          call bottom_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2)
       end do
       ! y last
       indices%j=grid%y%num
       call bottom_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2_last)

end subroutine

subroutine interior_angle_loop(mat, indices, ddx, ddy)
  type(space_angle_grid) grid
  type(rte_mat) mat
  type(index_list) indices
  integer l, m

  ! Allow derivative subroutines to be passed as arguments
  interface
     subroutine ddx(mat, indices)
       use sag
       use rte_sparse_matrices
       type(rte_mat) mat
       type(index_list) indices
     end subroutine ddx
     subroutine ddy(mat, indices)
       use sag
       use rte_sparse_matrices
       type(rte_mat) mat
       type(index_list) indices
     end subroutine ddy
  end interface

  do m=1, grid%phi%num
  indices%m = m
     do l=1, grid%theta%num
     indices%l = l
        call ddx(mat, indices)
        call ddy(mat, indices)
        call mat%z_cd2(indices)
        call mat%angular_integral(indices)
        call mat%attenuate(indices)
     end do
  end do
end subroutine


subroutine surface_angle_loop(mat, indices, ddx, ddy)
  type(space_angle_grid) grid
  type(rte_mat) mat
  type(index_list) indices
  integer l, m
  
  ! Allow derivative subroutines to be passed as arguments
  interface
     subroutine ddx(mat, indices)
       use sag
       use rte_sparse_matrices
       type(rte_mat) mat
       type(index_list) indices
     end subroutine ddx
     subroutine ddy(mat, indices)
       use sag
       use rte_sparse_matrices
       type(rte_mat) mat
       type(index_list) indices
     end subroutine ddy
  end interface

  ! Downwelling
  do m=1, grid%phi%num/2
  indices%m = m
     do l=1, grid%theta%num
     indices%l = l
        call mat%z_surface_bc(indices)
     end do
  end do

  ! Upwelling
  do m=grid%phi%num/2+1, grid%phi%num
  indices%m = m
     do l=1, grid%theta%num
     indices%l = l
        call ddx(mat, indices)
        call ddy(mat, indices)
        call mat%z_fd2(indices)
        call mat%angular_integral(indices)
        call mat%attenuate(indices)
     end do
  end do
end subroutine surface_angle_loop

subroutine bottom_angle_loop(mat, indices, ddx, ddy)
  type(space_angle_grid) grid
  type(rte_mat) mat
  type(index_list) indices
  integer l, m
  
  ! Allow derivative subroutines to be passed as arguments
  interface
     subroutine ddx(mat, indices)
       use sag
       use rte_sparse_matrices
       type(rte_mat) mat
       type(index_list) indices
     end subroutine ddx
     subroutine ddy(mat, indices)
       use sag
       use rte_sparse_matrices
       type(rte_mat) mat
       type(index_list) indices
     end subroutine ddy
  end interface

  ! Downwelling
  do m=1, grid%phi%num/2
  indices%m = m
     do l=1, grid%theta%num
     indices%l = l
        call ddx(mat, indices)
        call ddy(mat, indices)
        call mat%z_bd2(indices)
        call mat%angular_integral(indices)
        call mat%attenuate(indices)
     end do
  end do

  ! Upwelling
  do m=grid%phi%num/2+1, grid%phi%num
  indices%m = m
     do l=1, grid%theta%num
     indices%l = l
        call mat%z_bottom_bc(indices)
     end do
  end do
end subroutine bottom_angle_loop

subroutine gen_matrix(grid, mat, iops)
  type(rte_mat) mat
  type(space_angle_grid) grid
  type(optical_properties) iops
  type(index_list) indices

  call mat%init(grid, iops)

  call surface_space_loop(mat, indices)
  call interior_space_loop(mat, indices)
  call bottom_space_loop(mat, indices)

  call mat%deinit()

end subroutine gen_matrix

end module rte3d
