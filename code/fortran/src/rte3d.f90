module rte3d
use kelp_context
use rte_sparse_matrices
use light_context
implicit none

interface
   subroutine deriv_interface(mat, indices, ent)
     use rte_sparse_matrices
     class(rte_mat) mat
     type(index_list) indices
     integer(8) ent
   end subroutine deriv_interface
   subroutine angle_loop_interface(mat, indices, ddx, ddy)
     use rte_sparse_matrices
     import deriv_interface
     type(space_angle_grid) grid
     type(rte_mat) mat
     type(index_list) indices
     procedure(deriv_interface) :: ddx, ddy
   end subroutine angle_loop_interface
end interface

contains

subroutine whole_space_loop(mat, indices, num_threads)
  type(rte_mat) mat
  type(index_list) indices
  integer i, j, k
  integer num_threads

  procedure(deriv_interface), pointer :: ddx, ddy
  procedure(angle_loop_interface), pointer :: angle_loop

  !$omp parallel do default(none) shared(mat) &
  !$omp private(ddx,ddy,angle_loop, k, i, j) private(indices) &
  !$omp num_threads(num_threads) collapse(3)
  do k=1, mat%grid%z%num
     do i=1, mat%grid%x%num
        do j=1, mat%grid%y%num
           indices%k = k
           if(k .eq. 1) then
              angle_loop => surface_angle_loop
           else if(k .eq. mat%grid%z%num) then
              angle_loop => bottom_angle_loop
           else
              angle_loop => interior_angle_loop
           end if

           indices%i = i
           if(indices%i .eq. 1) then
              ddx => x_cd2_first
           else if(indices%i .eq. mat%grid%x%num) then
              ddx => x_cd2_last
           else
              ddx => x_cd2
           end if

           indices%j = j
           if(indices%j .eq. 1) then
              ddy => y_cd2_first
           else if(indices%j .eq. mat%grid%y%num) then
              ddy => y_cd2_last
           else
              ddy => y_cd2
           end if

           call angle_loop(mat, indices, ddx, ddy)
        end do
     end do
  end do
  !$omp end parallel do
end subroutine whole_space_loop

function calculate_start_ent(grid, indices) result(ent)
  type(space_angle_grid) grid
  type(index_list) indices
  integer(8) ent
  integer(8) boundary_nnz, interior_nnz
  integer(8) num_boundary, num_interior
  integer(8) num_this_x, num_this_z

  ! Nonzero matrix entries for an surface or bottom spatial grid cell
  ! Definitely an integer since nomega is even
  boundary_nnz = grid%angles%nomega * (2 * grid%angles%nomega + 11) / 2
  ! Nonzero matrix entries for an interior spatial grid cell
  interior_nnz = grid%angles%nomega * (grid%angles%nomega + 6)

  ! Order: z, x, y, omega
  ! Total number traversed so far in each spatial category
  ! row
  num_this_x = indices%j - 1
  ! depth layer
  num_this_z = (indices%i - 1) * grid%y%num + num_this_x

  ! Calculate number of spatial grid cells of each type which have
  ! already been traversed up to this point
  if(indices%k .eq. 1) then
     num_boundary = num_this_z
     num_interior = 0
  else if(indices%k .eq. grid%z%num) then
     num_boundary = (grid%x%num * grid%y%num) + num_this_z
     num_interior = (grid%z%num-2) * grid%x%num * grid%y%num
  else
     num_boundary = grid%x%num * grid%y%num
     num_interior = num_this_z + (indices%k-2) * grid%x%num * grid%y%num
  end if

  ent = num_boundary * boundary_nnz + num_interior * interior_nnz + 1
end function calculate_start_ent

function calculate_repeat_ent(ent, p) result(repeat_ent)
  integer p
  integer(8) ent, repeat_ent
  ! Entry number for row=mat%ind(i,j,k,p), col=mat%ind(i,j,k,p),
  ! which will be modified multiple times in this matrix row
  repeat_ent = ent + p - 1
end function calculate_repeat_ent

subroutine interior_angle_loop(mat, indices, ddx, ddy)
  type(rte_mat) mat
  type(index_list) indices
  procedure(deriv_interface) :: ddx, ddy
  integer p
  integer(8) ent, repeat_ent
  integer(8) row_num

  ! Determine which matrix row to start at
  ent = calculate_start_ent(mat%grid, indices)
  indices%p = 1
  row_num = mat%ind(indices%i, indices%j, indices%k, indices%p)

  do p=1, mat%grid%angles%nomega
     indices%p = p
     repeat_ent = calculate_repeat_ent(ent, p)
     call mat%set_row(ent, row_num)
     call mat%angular_integral(indices, ent)
     call ddx(mat, indices, ent)
     call ddy(mat, indices, ent)
     call mat%z_cd2(indices, ent)
     call mat%attenuate(indices, repeat_ent)
     call mat%add_source(indices, row_num)
     row_num = row_num + 1
  end do
end subroutine

subroutine surface_angle_loop(mat, indices, ddx, ddy)
  type(rte_mat) mat
  type(index_list) indices
  integer p
  procedure(deriv_interface) :: ddx, ddy
  integer(8) ent, repeat_ent
  integer(8) row_num

  ! Determine which matrix row to start at
  ent = calculate_start_ent(mat%grid, indices)
  indices%p = 1
  row_num = mat%ind(indices%i, indices%j, indices%k, indices%p)

  ! Downwelling
  do p=1, mat%grid%angles%nomega / 2
     indices%p = p
     repeat_ent = calculate_repeat_ent(ent, p)
     call mat%set_row(ent, row_num)
     call mat%angular_integral(indices, ent)
     call ddx(mat, indices, ent)
     call ddy(mat, indices, ent)
     call mat%z_surface_bc(indices, row_num, ent, repeat_ent)
     call mat%attenuate(indices, repeat_ent)
     call mat%add_source(indices, row_num)
     row_num = row_num + 1
  end do
  ! Upwelling
  do p=mat%grid%angles%nomega/2+1, mat%grid%angles%nomega
     indices%p = p
     repeat_ent = calculate_repeat_ent(ent, p)
     call mat%set_row(ent, row_num)
     call mat%angular_integral(indices, ent)
     call ddx(mat, indices, ent)
     call ddy(mat, indices, ent)
     call mat%z_fd2(indices, ent, repeat_ent)
     call mat%attenuate(indices, repeat_ent)
     call mat%add_source(indices, row_num)
     row_num = row_num + 1
  end do
end subroutine surface_angle_loop

subroutine bottom_angle_loop(mat, indices, ddx, ddy)
  type(rte_mat) mat
  type(index_list) indices
  integer p
  integer(8) row_num, ent, repeat_ent
  procedure(deriv_interface) :: ddx, ddy

  ! Determine which matrix row to start at
  ent = calculate_start_ent(mat%grid, indices)
  indices%p = 1
  row_num = mat%ind(indices%i, indices%j, indices%k, indices%p)

  ! Downwelling
  do p=1, mat%grid%angles%nomega/2
     indices%p = p
     repeat_ent = calculate_repeat_ent(ent, p)
     call mat%set_row(ent, row_num)
     call mat%angular_integral(indices, ent)
     call ddx(mat, indices, ent)
     call ddy(mat, indices, ent)
     call mat%z_bd2(indices, ent, repeat_ent)
     call mat%attenuate(indices, repeat_ent)
     call mat%add_source(indices, row_num)
     row_num = row_num + 1
  end do
  ! Upwelling
  do p=mat%grid%angles%nomega/2+1, mat%grid%angles%nomega
     indices%p = p
     repeat_ent = calculate_repeat_ent(ent, p)
     call mat%set_row(ent, row_num)
     call mat%angular_integral(indices, ent)
     call ddx(mat, indices, ent)
     call ddy(mat, indices, ent)
     call mat%z_bottom_bc(indices, ent, repeat_ent)
     call mat%attenuate(indices, repeat_ent)
     call mat%add_source(indices, row_num)
     row_num = row_num + 1
  end do
end subroutine bottom_angle_loop

subroutine gen_matrix(mat, num_threads)
  type(rte_mat) mat
  type(index_list) indices
  integer num_threads

  call indices%init()

  call whole_space_loop(mat, indices, num_threads)
  ! call surface_space_loop(mat, indices)
  ! call interior_space_loop(mat, indices)
  ! call bottom_space_loop(mat, indices)
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
