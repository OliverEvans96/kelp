module hdf5_utils
use hdf5
use utils
implicit none

contains

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

subroutine read_1d_double_array(filename, dsetname, arr, n)
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
end subroutine read_1d_double_array


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


subroutine read_int_array(array, dsetname, file_id, n, error)
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
end subroutine read_int_array

subroutine read_double_1d_dset(array, dsetname, file_id, n, error)
  integer(hid_t) :: file_id, dset_id
  integer(hsize_t), dimension(1) :: data_dims
  character(len=*) :: dsetname
  integer n
  double precision, dimension(n) :: array
  integer error

  data_dims(1) = n

  call h5dopen_f(file_id, dsetname, dset_id, error)
  call h5dread_f(dset_id, h5t_native_double, array, data_dims, error)
  call h5dclose_f(dset_id, error)
end subroutine read_double_1d_dset

function read_integer(file_id, dsetname, error) result(output)

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

end function read_integer

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

