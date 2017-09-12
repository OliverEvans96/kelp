module rte3d
use kelp_context
use rte_matrix
implicit none
contains

subroutine interior
  do indices%k=2 , grid%z%num - 1
     do indices%i=1, grid%x%num
        do indices%j=1, grid%y%num
           call interior_angle_loop(grid, mat, indices)
        end do
     end do
  end do
end subroutine

subroutine surface(grid, mat, kelp, indices)
  type(space_angle_grid) grid
  type(rte_mat) mat
  type(index_list) indices

  indices%k = 1
  do indices%i=1, grid%x%num
     do indices%j=1, grid%y%num
        call surface_angle_loop(grid, mat, indices)
     end do
  end do
end subroutine surface

subroutine bottom(grid, mat, kelp, indices)
  type(space_angle_grid) grid
  type(rte_mat) mat
  type(index_list) indices

  indices%k = 1
  do indices%i=1, grid%x%num
     do indices%j=1, grid%y%num
        call bottom_angle_loop(grid, mat, indices)
     end do
  end do
end subroutine bottom

subroutine interior_angle_loop(grid, mat, indices)
  type(space_angle_grid) grid
  type(rte_mat) mat
  type(index_list) indices
  do indices%m=1, grid%phi%num
     do indices%l=1, grid%theta%num
        call mat%attenuate(grid, kelp, indices)
        call mat%x_cd2(grid, indices)
        call mat%y_cd2(grid, indices)
        call mat%z_cd2(grid, indices)
        call mat%angular_integral(grid, kelp, indices)
     end do
  end do
end subroutine


subroutine surface_angle_loop(grid, mat, indices)
  type(space_angle_grid) grid
  type(rte_mat) mat
  type(index_list) indices
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
        call mat%x_cd2(grid, indices)
        call mat%y_cd2(grid, indices)
        call mat%z_fd2(grid, indices)
        call mat%angular_integral(grid, kelp, indices)
     end do
  end do
end subroutine surface_angle_loop

subroutine bottom_angle_loop(grid, mat, indices)
  type(space_angle_grid) grid
  type(rte_mat) mat
  type(index_list) indices
  ! Downwelling
  do indices%m=1, grid%phi%num/2
     do indices%l=1, grid%theta%num
        call mat%attenuate(grid, kelp, indices)
        call mat%x_cd2(grid, indices)
        call mat%y_cd2(grid, indices)
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

subroutine gen_matrix(grid, mat)
  type(space_angle_grid) grid
  type(rte_mat) mat
  type(index_list) indices

  mat%init(grid)

  do indices%i=1, grid%x%num
     do indices%j=1, grid%y%num
        do indices%m=1, grid%phi%num/2
           do indices%l=1, grid%theta%num
              call mat%attenuate(grid, kelp, indices)
              call mat%x_cd2(grid, indices)
              call mat%y_cd2(grid, indices)
              call mat%z_cd2(grid, indices)
              call mat%angular_integral(grid, kelp, indices)
           end do
        end do
     end do
  end do

  do indices%k=2 , grid%z%num - 1
      do indices%i=1, grid%x%num
        do indices%j=1, grid%y%num
            do indices%m=1, grid%phi%num/2
              do indices%l=1, grid%theta%num
                  call mat%attenuate(grid, kelp, indices)
                  call mat%x_cd2(grid, indices)
                  call mat%y_cd2(grid, indices)
                  call mat%z_cd2(grid, indices)
                  call mat%angular_integral(grid, kelp, indices)
              end do
            end do
        end do
      end do
  end do

  do indices%i=1, grid%x%num
     do indices%j=1, grid%y%num
        do indices%m=1, grid%phi%num/2
           do indices%l=1, grid%theta%num
              call mat%attenuate(grid, kelp, indices)
              call mat%x_cd2(grid, indices)
              call mat%y_cd2(grid, indices)
              call mat%z_cd2(grid, indices)
              call mat%angular_integral(grid, kelp, indices)
           end do
        end do
     end do
  end do

  mat%deinit()

end subroutine gen_matrix

end module rte3d
