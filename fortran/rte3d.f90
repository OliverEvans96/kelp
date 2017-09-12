module rte3d
use kelp_context
use rte_matrix
implicit none
contains

subroutine interior(grid, mat, kelp, indices)
  class(rte_mat) rte_matrix
  type(space_angle_grid) grid
  type(kelp_state) kelp
  type(index_list) indices

  ! z interior
  do indices%k=2 , grid%z%num - 1
     ! x first
     indices%i=1
       ! y first
       indices%j=1
       call interior_angle_loop(grid, mat, indices, mat%x_cd2_first, mat%y_cd2_first)
       ! y interior
       do indices%j=2, grid%y%num - 1
           call interior_angle_loop(grid, mat, indices, mat%x_cd2_first, mat%y_cd2)
       end do
       ! y last
       indices%j=grid%y%num
       call interior_angle_loop(grid, mat, indices, mat%x_cd2_first, mat%y_cd2_last)

     ! x interior
     do indices%i=2, grid%x%num - 1
        ! y first
        indices%j=1
        call interior_angle_loop(grid, mat, indices, mat%x_cd2, mat%y_cd2_first)
        ! y interior
        do indices%j=2, grid%y%num - 1
           call interior_angle_loop(grid, mat, indices, mat%x_cd2, mat%y_cd2)
        end do
        ! y last
        indices%j=grid%y%num
        call interior_angle_loop(grid, mat, indices, mat%x_cd2, mat%y_cd2_last)
     end do

     ! x last
     indices%i=grid%x%num
       ! y first
       indices%j=1
       call interior_angle_loop(grid, mat, indices, mat%x_cd2_last, mat%y_cd2_first)
       ! y interior
       do indices%j=2, grid%y%num - 1
          call interior_angle_loop(grid, mat, indices, mat%x_cd2_last, mat%y_cd2)
       end do
       ! y last
       indices%j=grid%y%num
       call interior_angle_loop(grid, mat, indices, mat%x_cd2_last, mat%y_cd2_last)

  end do

end subroutine


subroutine surface(grid, mat, kelp, indices)
  class(rte_mat) rte_matrix
  type(space_angle_grid) grid
  type(kelp_state) kelp
  type(index_list) indices

  ! z surface
  indices%k=1
     ! x first
     indices%i=1
       ! y first
       indices%j=1
       call surface_angle_loop(grid, mat, indices, mat%x_cd2_first, mat%y_cd2_first)
       ! y surface
       do indices%j=2, grid%y%num - 1
           call surface_angle_loop(grid, mat, indices, mat%x_cd2_first, mat%y_cd2)
       end do
       ! y last
       indices%j=grid%y%num
       call surface_angle_loop(grid, mat, indices, mat%x_cd2_first, mat%y_cd2_last)

     ! x surface
     do indices%i=2, grid%x%num - 1
        ! y first
        indices%j=1
        call surface_angle_loop(grid, mat, indices, mat%x_cd2, mat%y_cd2_first)
        ! y surface
        do indices%j=2, grid%y%num - 1
           call surface_angle_loop(grid, mat, indices, mat%x_cd2, mat%y_cd2)
        end do
        ! y last
        indices%j=grid%y%num
        call surface_angle_loop(grid, mat, indices, mat%x_cd2, mat%y_cd2_last)
     end do

     ! x last
     indices%i=grid%x%num
       ! y first
       indices%j=1
       call surface_angle_loop(grid, mat, indices, mat%x_cd2_last, mat%y_cd2_first)
       ! y surface
       do indices%j=2, grid%y%num - 1
          call surface_angle_loop(grid, mat, indices, mat%x_cd2_last, mat%y_cd2)
       end do
       ! y last
       indices%j=grid%y%num
       call surface_angle_loop(grid, mat, indices, mat%x_cd2_last, mat%y_cd2_last)

end subroutine

subroutine bottom(grid, mat, kelp, indices)
  class(rte_mat) rte_matrix
  type(space_angle_grid) grid
  type(kelp_state) kelp
  type(index_list) indices

  ! z bottom
  indices%k=grid%z%num
     ! x first
     indices%i=1
       ! y first
       indices%j=1
       call bottom_angle_loop(grid, mat, indices, mat%x_cd2_first, mat%y_cd2_first)
       ! y bottom
       do indices%j=2, grid%y%num - 1
           call bottom_angle_loop(grid, mat, indices, mat%x_cd2_first, mat%y_cd2)
       end do
       ! y last
       indices%j=grid%y%num
       call bottom_angle_loop(grid, mat, indices, mat%x_cd2_first, mat%y_cd2_last)

     ! x bottom
     do indices%i=2, grid%x%num - 1
        ! y first
        indices%j=1
        call bottom_angle_loop(grid, mat, indices, mat%x_cd2, mat%y_cd2_first)
        ! y bottom
        do indices%j=2, grid%y%num - 1
           call bottom_angle_loop(grid, mat, indices, mat%x_cd2, mat%y_cd2)
        end do
        ! y last
        indices%j=grid%y%num
        call bottom_angle_loop(grid, mat, indices, mat%x_cd2, mat%y_cd2_last)
     end do

     ! x last
     indices%i=grid%x%num
       ! y first
       indices%j=1
       call bottom_angle_loop(grid, mat, indices, mat%x_cd2_last, mat%y_cd2_first)
       ! y bottom
       do indices%j=2, grid%y%num - 1
          call bottom_angle_loop(grid, mat, indices, mat%x_cd2_last, mat%y_cd2)
       end do
       ! y last
       indices%j=grid%y%num
       call bottom_angle_loop(grid, mat, indices, mat%x_cd2_last, mat%y_cd2_last)

end subroutine

subroutine interior_angle_loop(grid, mat, indices, ddx, ddy)
  type(space_angle_grid) grid
  type(rte_mat) mat
  type(index_list) indices

  ! Allow derivative subroutines to be passed as arguments
  interface
     subroutine ddx
       type(space_angle_grid) grid
       type(index_list) indices
     end subroutine ddx
     subroutine ddy
       type(space_angle_grid) grid
       type(indey_list) indices
     end subroutine ddy
  end interface

  do indices%m=1, grid%phi%num
     do indices%l=1, grid%theta%num
        call mat%attenuate(grid, kelp, indices)
        call ddx(grid, indices)
        call ddy(grid, indices)
        call mat%z_cd2(grid, indices)
        call mat%angular_integral(grid, kelp, indices)
     end do
  end do
end subroutine


subroutine surface_angle_loop(grid, mat, indices, ddx, ddy)
  type(space_angle_grid) grid
  type(rte_mat) mat
  type(index_list) indices
  
  ! Allow derivative subroutines to be passed as arguments
  interface
     subroutine ddx
       type(space_angle_grid) grid
       type(index_list) indices
     end subroutine ddx
     subroutine ddy
       type(space_angle_grid) grid
       type(indey_list) indices
     end subroutine ddy
  end interface

  ! Downwelling
  do indices%m=1, grid%phi%num/2
     do indices%l=1, grid%theta%num
        call mat%z_surface_bc(grid, indices)
     end do
  end do

  ! Upwelling
  do indices%m=grid%phi%num/2+1, grid%phi%num
     do indices%l=1, grid%theta%num
        call mat%attenuate(grid, kelp, indices)
        call ddx(grid, indices)
        call ddy(grid, indices)
        call mat%z_fd2(grid, indices)
        call mat%angular_integral(grid, kelp, indices)
     end do
  end do
end subroutine surface_angle_loop

subroutine bottom_angle_loop(grid, mat, indices, ddx, ddy)
  type(space_angle_grid) grid
  type(rte_mat) mat
  type(index_list) indices
  
  ! Allow derivative subroutines to be passed as arguments
  interface
     subroutine ddx
       type(space_angle_grid) grid
       type(index_list) indices
     end subroutine ddx
     subroutine ddy
       type(space_angle_grid) grid
       type(indey_list) indices
     end subroutine ddy
  end interface

  ! Downwelling
  do indices%m=1, grid%phi%num/2
     do indices%l=1, grid%theta%num
        call mat%attenuate(grid, kelp, indices)
        call ddx(grid, indices)
        call ddy(grid, indices)
        call mat%z_bd2(grid, indices)
        call mat%angular_integral(grid, kelp, indices)
     end do
  end do

  ! Upwelling
  do indices%m=grid%phi%num/2+1, grid%phi%num
     do indices%l=1, grid%theta%num
        call mat%z_bottom_bc(grid, indices)
     end do
  end do
end subroutine bottom_angle_loop

subroutine gen_matrix(grid, kelp, mat)
  type(rte_mat) mat
  type(space_angle_grid) grid
  type(kelp_state) kelp
  type(index_list) indices

  mat%init(grid)

  call surface(mat, grid, kelp, indices)
  call interior(mat, grid, kelp, indices)
  call bottom(mat, grid, kelp, indices)

  mat%deinit()

end subroutine gen_matrix

end module rte3d
