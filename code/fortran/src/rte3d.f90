module rte3d
use kelp_context
use rte_sparse_matrices
use light_context
implicit none

contains

subroutine interior_space_loop(mat, indices)
  type(rte_mat) mat
  type(space_angle_grid) grid
  type(index_list) indices
  integer i, j, k

  grid = mat%grid

  ! z interior
  !$OMP PARALLEL DO
  do k=2 , grid%z%num - 1
     indices%k = k
     write(*,*) 'k =', indices%k, '/', grid%z%num
     ! x first
     indices%i=1
       ! y first
       indices%j=1
       call interior_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2_first)

       ! y interior
       !$OMP PARALLEL DO
       do j=2, grid%y%num - 1
       indices%j = j
           call interior_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2)
       end do
       !$OMP END PARALLEL DO
       ! y last
       indices%j=grid%y%num
       call interior_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2_last)

     ! x interior
     !$OMP PARALLEL DO
     do i=2, grid%x%num - 1
     indices%i = i
        ! y first
        indices%j=1
        call interior_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2_first)
        ! y interior
        !$OMP PARALLEL DO
        do j=2, grid%y%num - 1
        indices%j = j
           call interior_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2)
        end do
        !$OMP END PARALLEL DO
        ! y last
        indices%j=grid%y%num
           call interior_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2_last)
     end do
     !$OMP END PARALLEL DO

     ! x last
     indices%i=grid%x%num
       ! y first
       indices%j=1
          call interior_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2_first)
       ! y interior
       !$OMP PARALLEL DO
       do j=2, grid%y%num - 1
          indices%j = j
          call interior_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2)
       end do
       !$OMP END PARALLEL DO
       ! y last
       indices%j=grid%y%num
          call interior_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2_last)

  end do
  !$OMP END PARALLEL DO

end subroutine


subroutine surface_space_loop(mat, indices)
  type(rte_mat) mat
  type(space_angle_grid) grid
  type(index_list) indices
  integer i, j

  grid = mat%grid

  ! z surface
  indices%k=1
     write(*,*) 'k =', indices%k, '/', grid%z%num
     ! x first
     indices%i=1
        ! y first
        indices%j=1
           call surface_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2_first)
        ! y interior
        !$OMP PARALLEL DO
        do j=2, grid%y%num - 1
           indices%j = j
           call surface_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2)
        end do
        !$OMP END PARALLEL DO
        ! y last
        indices%j=grid%y%num
           call surface_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2_last)

     ! x interior
     !$OMP PARALLEL DO
     do i=2, grid%x%num - 1
     indices%i = i
        ! y first
        indices%j=1
           call surface_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2_first)
        ! y interior
        !$OMP PARALLEL DO
        do j=2, grid%y%num - 1
        indices%j = j
           call surface_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2)
        end do
        !$OMP END PARALLEL DO
        ! y last
        indices%j=grid%y%num
           call surface_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2_last)
     end do
     !$OMP END PARALLEL DO

     ! x last
     indices%i=grid%x%num
       ! y first
       indices%j=1
       call surface_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2_first)
       ! y surface
       !$OMP PARALLEL DO
       do j=2, grid%y%num - 1
       indices%j = j
          call surface_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2)
       end do
       !$OMP END PARALLEL DO
       ! y last
       indices%j=grid%y%num
       call surface_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2_last)

end subroutine surface_space_loop

subroutine bottom_space_loop(mat, indices)
  type(rte_mat) mat
  type(space_angle_grid) grid
  type(index_list) indices
  integer i, j

  grid = mat%grid

  ! z bottom
  indices%k=grid%z%num
     write(*,*) 'k =', indices%k, '/', grid%z%num
     ! x first
     indices%i=1
       ! y first
       indices%j=1
          call bottom_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2_first)
       ! y interior
       !$OMP PARALLEL DO
       do j=2, grid%y%num - 1
          indices%j = j
          call bottom_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2)
       end do
       !$OMP END PARALLEL DO
       ! y last
       indices%j=grid%y%num
          call bottom_angle_loop(mat, indices, wrap_x_cd2_first, wrap_y_cd2_last)

     ! x interior
     !$OMP PARALLEL DO
     do i=2, grid%x%num - 1
     indices%i = i
        ! y first
        indices%j=1
           call bottom_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2_first)
        ! y bottom
        !$OMP PARALLEL DO
        do j=2, grid%y%num - 1
           indices%j = j
           call bottom_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2)
        end do
        !$OMP END PARALLEL DO
        ! y last
        indices%j=grid%y%num
           call bottom_angle_loop(mat, indices, wrap_x_cd2, wrap_y_cd2_last)
     end do
     !$OMP END PARALLEL DO

     ! x last
     indices%i=grid%x%num
       ! y first
       indices%j=1
       call bottom_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2_first)
       ! y interior
       !$OMP PARALLEL DO
       do j=2, grid%y%num - 1
          indices%j = j
          call bottom_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2)
       end do
       !$OMP END PARALLEL DO
       ! y last
       indices%j=grid%y%num
          call bottom_angle_loop(mat, indices, wrap_x_cd2_last, wrap_y_cd2_last)

end subroutine

subroutine interior_angle_loop(mat, indices, ddx, ddy)
  type(space_angle_grid) grid
  type(rte_mat) mat
  type(index_list) indices
  integer p

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

  grid = mat%grid

  do p=1, grid%angles%nomega
      indices%p = p
      call mat%angular_integral(indices)
      call ddx(mat, indices)
      call ddy(mat, indices)
      call mat%z_cd2(indices)
      call mat%attenuate(indices)
  end do
end subroutine


subroutine surface_angle_loop(mat, indices, ddx, ddy)
  type(space_angle_grid) grid
  type(rte_mat) mat
  type(index_list) indices
  integer p

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

  grid = mat%grid

  ! Downwelling
  ! TODO: Downwelling vs upwelling
  do p=1, grid%angles%nomega / 2
     call mat%z_surface_bc(indices)
  end do

  ! Upwelling
  do p=grid%angles%nomega/2+1, grid%angles%nomega
     indices%p = p
     call mat%angular_integral(indices)
     call ddx(mat, indices)
     call ddy(mat, indices)
     call mat%z_fd2(indices)
     call mat%attenuate(indices)
  end do

end subroutine surface_angle_loop

subroutine bottom_angle_loop(mat, indices, ddx, ddy)
  type(space_angle_grid) grid
  type(rte_mat) mat
  type(index_list) indices
  integer p

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

  grid = mat%grid

  ! Downwelling
  do p=1, grid%angles%nomega
     indices%p = p
     call mat%angular_integral(indices)
     call ddx(mat, indices)
     call ddy(mat, indices)
     call mat%z_bd2(indices)
     call mat%attenuate(indices)
  end do

  ! Upwelling
  do p=grid%angles%nomega/2+1, grid%angles%nomega
     indices%p = p
     call mat%z_bottom_bc(indices)
  end do

end subroutine bottom_angle_loop

subroutine gen_matrix(mat)
  type(rte_mat) mat
  type(index_list) indices

  call indices%init()

  call surface_space_loop(mat, indices)
  call interior_space_loop(mat, indices)
  call bottom_space_loop(mat, indices)
end subroutine gen_matrix

subroutine rte3d_deinit(mat, iops, light)
  type(rte_mat) mat
  type(optical_properties) iops
  type(light_state) light

  call mat%deinit()
  call iops%deinit()
  call light%deinit()
end subroutine

end module rte3d
