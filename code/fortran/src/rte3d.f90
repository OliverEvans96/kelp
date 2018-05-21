module rte3d
use kelp_context
use rte_sparse_matrices
use light_context
implicit none

interface
   subroutine deriv_interface(mat, indices, row_num)
     use rte_sparse_matrices
     class(rte_mat) mat
     type(index_list) indices
     integer row_num
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

subroutine whole_space_loop(mat, indices)
  type(rte_mat) mat
  type(index_list) indices
  integer i, j, k

  procedure(deriv_interface), pointer :: ddx, ddy
  procedure(angle_loop_interface), pointer :: angle_loop

  !$ integer omp_get_max_threads
  !$ integer num_threads_z, num_threads_x, num_threads_y

  ! Enable nested parallelism
  !$ call omp_set_nested(.true.)

  ! Use nz procs for outer loop,
  ! or num_procs if num_procs < nz
  ! Divide the rest of the tasks as appropriate

  !$ num_threads_z = min(omp_get_max_threads(), mat%grid%z%num)
  !$ num_threads_x = min( &
  !$    omp_get_max_threads()/num_threads_z, &
  !$    mat%grid%x%num)
  !$ num_threads_y = min( &
  !$    omp_get_max_threads()/(num_threads_z*num_threads_x), &
  !$    mat%grid%y%num)

  !$ write(*,*) 'num_procs =', omp_get_max_threads()
  !$ write(*,*) 'ntz =', num_threads_z
  !$ write(*,*) 'ntx =', num_threads_x
  !$ write(*,*) 'nty =', num_threads_y

  !$omp parallel do default(none) shared(mat) &
  !$omp private(ddx,ddy,angle_loop, k, i, j) private(indices) &
  !$omp shared(num_threads_x,num_threads_y,num_threads_z) &
  !$omp num_threads(num_threads_z) if(num_threads_z .gt. 1)
  do k=1, mat%grid%z%num
     write(*,*) 'k =', k
     indices%k = k
     if(k .eq. 1) then
        angle_loop => surface_angle_loop
     else if(k .eq. mat%grid%z%num) then
        angle_loop => bottom_angle_loop
     else
        angle_loop => interior_angle_loop
     end if

     !$omp parallel do default(none) shared(mat) private(i,j) &
     !$omp firstprivate(indices,angle_loop, k) private(ddx,ddy) &
     !$omp shared(num_threads_x,num_threads_y,num_threads_z) &
     !$omp num_threads(num_threads_x) if(num_threads_x .gt. 1)
     do i=1, mat%grid%x%num
        indices%i = i
        if(indices%i .eq. 1) then
           ddx => x_cd2_first
        else if(indices%i .eq. mat%grid%x%num) then
           ddx => x_cd2_last
        else
           ddx => x_cd2
        end if
        !$omp parallel do default(none) shared(mat) private(j) &
        !$omp firstprivate(indices,ddx,ddy,angle_loop, i, k) &
        !$omp shared(num_threads_x,num_threads_y,num_threads_z) &
        !$omp num_threads(num_threads_y) if(num_threads_y .gt. 1)
        do j=1, mat%grid%y%num
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
        !$omp end parallel do
     end do
     !$omp end parallel do
  end do
  !$omp end parallel do
end subroutine whole_space_loop
subroutine interior_angle_loop(mat, indices, ddx, ddy)
  type(rte_mat) mat
  type(index_list) indices
  procedure(deriv_interface) :: ddx, ddy
  integer p
  integer row_num

  ! Determine which matrix row to start at
  indices%p = 1
  row_num = mat%ind(indices%i, indices%j, indices%k, indices%p)

  do p=1, mat%grid%angles%nomega
     indices%p = p
     call mat%angular_integral(indices, row_num)
     call ddx(mat, indices, row_num)
     call ddy(mat, indices, row_num)
     call mat%z_cd2(indices, row_num)
     call mat%attenuate(indices, row_num)
     row_num = row_num + 1
  end do
end subroutine

subroutine surface_angle_loop(mat, indices, ddx, ddy)
  type(rte_mat) mat
  type(index_list) indices
  integer p
  procedure(deriv_interface) :: ddx, ddy
  integer row_num

  ! Determine which matrix row to start at
  indices%p = 1
  row_num = mat%ind(indices%i, indices%j, indices%k, indices%p)

  ! Downwelling
  do p=1, mat%grid%angles%nomega / 2
     indices%p = p
     call mat%angular_integral(indices, row_num)
     call ddx(mat, indices, row_num)
     call ddy(mat, indices, row_num)
     call mat%z_surface_bc(indices, row_num)
     call mat%attenuate(indices, row_num)
     row_num = row_num + 1
  end do
  ! Upwelling
  do p=mat%grid%angles%nomega/2+1, mat%grid%angles%nomega
     indices%p = p
     call mat%angular_integral(indices, row_num)
     call ddx(mat, indices, row_num)
     call ddy(mat, indices, row_num)
     call mat%z_fd2(indices, row_num)
     call mat%attenuate(indices, row_num)
     row_num = row_num + 1
  end do
end subroutine surface_angle_loop

subroutine bottom_angle_loop(mat, indices, ddx, ddy)
  type(rte_mat) mat
  type(index_list) indices
  integer p
  integer row_num
  procedure(deriv_interface) :: ddx, ddy

  ! Determine which matrix row to start at
  indices%p = 1
  row_num = mat%ind(indices%i, indices%j, indices%k, indices%p)

  ! Downwelling
  do p=1, mat%grid%angles%nomega/2
     indices%p = p
     call mat%angular_integral(indices, row_num)
     call ddx(mat, indices, row_num)
     call ddy(mat, indices, row_num)
     call mat%z_bd2(indices, row_num)
     call mat%attenuate(indices, row_num)
     row_num = row_num + 1
  end do
  ! Upwelling
  do p=mat%grid%angles%nomega/2+1, mat%grid%angles%nomega
     indices%p = p
     call mat%angular_integral(indices, row_num)
     call ddx(mat, indices, row_num)
     call ddy(mat, indices, row_num)
     call mat%z_bottom_bc(indices, row_num)
     call mat%attenuate(indices, row_num)
     row_num = row_num + 1
  end do
end subroutine bottom_angle_loop

subroutine gen_matrix(mat)
  type(rte_mat) mat
  type(index_list) indices

  call indices%init()

  call whole_space_loop(mat, indices)
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
