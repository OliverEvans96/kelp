module hdf5_utils
use hdf5

contains

subroutine read_coo(filename, allocatable_array)
  implicit none
  integer error
  integer(hsize_t), dimension(2) :: data_dims
  integer m, n, nnz
  double precision, dimension(:,:), allocatable :: a
  double precision, dimension(:,:), allocatable :: b, x

  integer, dimension(:), allocatable :: row, col, data
  double precision, dimension(:,:), allocatable :: allocatable_array

  integer(hid_t) :: file_id

  character(len=*) :: filename
  !character(len=3) :: dsetname = "nnz"

  call h5open_f(error)
  call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)

  nnz = read_integer(file_id, "nnz", error)
  m = read_integer(file_id, "m", error)
  n = read_integer(file_id, "n", error)

  allocate(row(nnz))
  allocate(col(nnz))
  allocate(data(nnz))

  !!! READ 1D ARRAY HERE !!!

  write(*,*) "nnz = ", nnz
  write(*,*) "m = ", m
  write(*,*) "n = ", n

  call h5fclose_f(file_id, error)
  call h5close_f(error)
end subroutine

function read_integer(file_id, dsetname, error) result(output)
  implicit none

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


end module hdf5_utils

