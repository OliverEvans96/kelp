module hdf5_utils
use hdf5
use utils
use kelp_context
implicit none

contains

  subroutine hdf_read_grid(filename, grid)
    integer error
    type(space_angle_grid) grid
    double precision, dimension(6) :: bounds
    double precision, dimension(5) :: spacings

    integer(hid_t) :: file_id
    character(len=*) :: filename

    call h5open_f(error)
    call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)

    call grid_read_bounds(bounds, file_id, error)
    call grid_read_spacings(spacings, file_id, error)

    call h5fclose_f(file_id, error)
    call h5close_f(error)

    call grid%set_num_from_spacing()
    call grid%init()
  end subroutine hdf_read_grid

  subroutine hdf_read_rope(filename, rope)
    integer error
    type(rope_state) rope
    double precision, dimension(:) :: kelp_lengths, kelp_stds, water_speeds, water_angles

    integer(hid_t) :: file_id
    character(len=*) :: filename

    call h5open_f(error)
    call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)
    integer(hsize_t), dimension(1) :: data_dims

    data_dims(1) = rope%nz

    call dset_read_1d_double_array(file_id, 'kelp_lengths', data_dims, kelp_lengths, error)
    call dset_read_1d_double_array(file_id, 'kelp_stds', data_dims, kelp_stds, error)
    call dset_read_1d_double_array(file_id, 'water_speeds', data_dims, water_speeds, error)
    call dset_read_1d_double_array(file_id, 'water_angles', data_dims, water_angles, error)

    call h5fclose_f(file_id, error)
    call h5close_f(error)

    call rope%init()
  end subroutine hdf_read_rope

  subroutine hdf_read_kelp(filename, p_kelp)
  end subroutine hdf_read_kelp

  !subroutine hdf_read_light()
  !end subroutine hdf_read_light

  !- hdf WRITE kelp

  subroutine hdf_write_kelp(filename)
  end subroutine hdf_write_kelp

  ! subroutine hdf_read_rad(filename, rad)
  ! end subroutine hdf_read_rad()

  ! subroutine hdf_read_irrad()
  ! end subroutine hdf_read_irrad

  ! - grid utils
  subroutine grid_read_bounds(bounds, file_id, error)
    integer, parameter :: num_bounds = 6
    integer(hid_t) :: file_id
    integer error
    double precision, dimension(num_bounds) :: bounds 
    integer(hsize_t), dimension(1) :: data_dims

    data_dims(1) = num_bounds

    call dset_read_1d_double_array(file_id, 'bounds', data_dims, bounds, error)
  end subroutine grid_read_bounds

  subroutine grid_read_nums(nums, file_id, error)
    integer, parameter :: num_nums = 5
    integer(hid_t) :: file_id
    integer error
    integer, dimension(num_nums) :: nums 
    integer(hsize_t), dimension(1) :: data_dims

    data_dims(1) = num_nums

    call dset_read_1d_double_array(file_id, 'nums', data_dims, nums, error)
  end subroutine grid_read_nums

  subroutine grid_read_spacings(spacings, file_id, error)
    integer, parameter :: num_spacings = 5
    integer(hid_t) :: file_id
    integer error
    double precision, dimension(num_spacings) :: spacings 
    integer(hsize_t), dimension(1) :: data_dims

    array_dim(1) = num_spacings

    call dset_read_1d_double_array(file_id, 'spacings', data_dims, spacings, error)
  end subroutine grid_read_spacings


  !- dset READ int

  subroutine dset_read_int(file_id, dsetname, output, error)

    integer(hid_t) :: file_id, dset_id
    integer(hsize_t), dimension(1) :: int_data_dim
    integer, dimension(1) :: tmp_int_arr
    character(len=*) :: dsetname
    integer error
    integer output

    int_data_dim(1) = 1

    call h5dopen_f(file_id, dsetname, dset_id, error)
    call h5dread_f(dset_id, h5t_native_integer, tmp_int_arr, int_data_dim, error)
    call h5dclose_f(dset_id, error)

    output = tmp_int_arr(1)

  end subroutine dset_read_int

  subroutine dset_read_1d_int_array(file_id, dsetname, data_dims, output, error)
    integer(hid_t) :: file_id, dset_id
    integer(hsize_t), dimension(1) :: data_dims
    character(len=*) :: dsetname
    integer, dimension(n) :: array
    integer error

    call h5dopen_f(file_id, dsetname, dset_id, error)
    call h5dread_f(dset_id, h5t_native_integer, array, data_dims, error)
    call h5dclose_f(dset_id, error)
  end subroutine dset_read_1d_int_array


  !- dset READ double
  subroutine dset_read_1d_double_array(file_id, dsetname, data_dims, output, error)
    integer(hid_t) :: file_id, dset_id
    integer(hsize_t), dimension(1) :: data_dims
    character(len=*) :: dsetname
    double precision, dimension(n) :: array
    integer error

    call h5dopen_f(file_id, dsetname, dset_id, error)
    call h5dread_f(dset_id, h5t_native_double, array, data_dims, error)
    call h5dclose_f(dset_id, error)
  end subroutine dset_read_1d_double_array

  ! subroutine dset_read_3d_double_array()
  ! end subroutine dset_read_3d_double_array

  ! subroutine dset_read_5d_double_array()
  ! end subroutine dset_read_5d_double_array


  !- dset WRITE int


  !- dset WRITE double
  subroutine dset_write_1d_double_array(file_id, dsetname, data_dims, input, error)
    integer(hid_t) :: file_id, dset_id
    integer(hsize_t), dimension(1) :: data_dims
    character(len=*) :: dsetname
    integer n
    integer, dimension(n) :: array
    integer error

    data_dims(1) = n

    call h5dopen_f(file_id, dsetname, dset_id, error)
    call h5dread_f(dset_id, h5t_native_integer, array, data_dims, error)
    call h5dclose_f(dset_id, error)
  end subroutine dset_write_1d_double_array

  subroutine dset_write_3d_double_array(file_id, dsetname, data_dims, input, error)
  end subroutine dset_write_3d_double_array

  !subroutine dset_write_5d_double_array(file_id, dsetname, data_dims, input, error)
  !end subroutine dset_write_5d_double_array




  !- hdf READ double
subroutine hdf_read_1d_double_array(filename, dsetname, arr, n)
  integer error
  integer n
  character(len=*) filename, dsetname
  integer(hid_t) :: file_id
  double precision, dimension(n) :: arr

  call h5open_f(error)
  call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)

  call read_double_1d_dset(arr, dsetname, file_id, n, error)

  call h5fclose_f(file_id, error)
  call h5close_f(error)
end subroutine hdf_read_1d_double_array


  !- hdf COO

subroutine read_1d_array_from_coo(filename, arr, n, nnz)
  integer error
  integer n, m, nnz

  integer, dimension(nnz) :: row, col
  double precision, dimension(nnz) :: data
  double precision, dimension(n) :: arr

  integer(hid_t) :: file_id
  character(len=*) :: filename

  m = 1

  call h5open_f(error)
  call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)

  call read_int_array(row, 'row', file_id, nnz, error)
  call read_int_array(col, 'col', file_id, nnz, error)
  call read_double_1d_dset(data, 'data', file_id, nnz, error)

  call h5fclose_f(file_id, error)
  call h5close_f(error)

  ! write(*,*) 'row'
  ! call print_int_array(row, nnz, 1)
  ! write(*,*) 'col'
  ! call print_int_array(col, nnz, 1)
  ! write(*,*) 'data'
  ! call print_array(data, nnz, 1)

  call dense_from_coo(row, col, data, arr, n, m, nnz)

  !write(*,*) 'arr'
  !call print_array(arr, n, m)

end subroutine read_1d_array_from_coo

subroutine read_coo(filename, row, col, data, nnz)
  integer error
  integer n, m, nnz

  integer, dimension(nnz) :: row, col
  double precision, dimension(nnz) :: data

  integer(hid_t) :: file_id
  character(len=*) :: filename

  call h5open_f(error)
  call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)

  call read_int_array(row, 'row', file_id, nnz, error)
  call read_int_array(col, 'col', file_id, nnz, error)
  call read_double_1d_dset(data, 'data', file_id, nnz, error)

  call h5fclose_f(file_id, error)
  call h5close_f(error)
end subroutine

subroutine read_info(filename, n, m, nnz)
  integer error
  integer n, m, nnz

  integer(hid_t) :: file_id

  character(len=*) :: filename

  call h5open_f(error)
  call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)

  nnz = read_integer(file_id, "nnz", error)
  m = read_integer(file_id, "m", error)
  n = read_integer(file_id, "n", error)

  call h5fclose_f(file_id, error)
  call h5close_f(error)
end subroutine read_info

!- COO utils

subroutine coo_from_dense(a, n, m, row, col, data)
  integer n, m
  double precision, dimension(n, m) :: a
  integer, dimension(n*m) :: row, col
  double precision, dimension(n*m) :: data
  integer i, j, k

  k = 1
  do i=1, m
     do j=1, n
        row(k) = i
        col(k) = j
        data(k) = a(i,j)

        k = k+1
     end do
  end do
end subroutine coo_from_dense

subroutine dense_from_coo(row, col, data, a, n, m, nnz)
  integer n, m, nnz
  double precision, dimension(n, m) :: a
  integer, dimension(nnz) :: row, col
  double precision, dimension(nnz) :: data
  integer i, j, k

  ! This could definitely be improved.
  do i = 1, n
     do j = 1, m
        a(i,j) = 0
     end do
  end do

   do k=1, nnz
      i = row(k)
      j = col(k)
      a(i,j) = data(k)
   end do

end subroutine dense_from_coo


end module hdf5_utils

