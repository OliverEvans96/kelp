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
    !double precision, dimension(5) :: spacings
    integer, dimension(5) :: nums

    integer(hid_t), dimension(1) :: data_dims
    integer(hid_t) :: file_id
    character(len=*) :: filename

    call h5open_f(error)
    call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)

    data_dims(1) = 6
    call dset_read_1d_double(file_id, 'bounds', data_dims, bounds, error)
    data_dims(1) = 5
    call dset_read_1d_int(file_id, 'nums', data_dims, nums, error)

    call h5fclose_f(file_id, error)
    call h5close_f(error)

    grid%x%minval = bounds(1)
    grid%x%maxval = bounds(2)
    grid%y%minval = bounds(3)
    grid%y%maxval = bounds(4)
    grid%z%minval = bounds(5)
    grid%z%maxval = bounds(6)

    grid%x%num = nums(1)
    grid%y%num = nums(2)
    grid%z%num = nums(3)
    grid%theta%num = nums(4)
    grid%phi%num = nums(5)

    call grid%set_spacing_from_num()
    call grid%init()
  end subroutine hdf_read_grid

  subroutine hdf_read_rope(filename, rope, grid)
    integer error
    type(rope_state) rope
    type(space_angle_grid) grid

    integer(hid_t) :: file_id
    character(len=*) :: filename
    integer(hsize_t), dimension(1) :: data_dims

    call h5open_f(error)
    call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)

    call rope%init(grid)
    data_dims(1) = rope%nz

    call dset_read_1d_double(file_id, 'frond_lengths', data_dims, rope%frond_lengths, error)
    call dset_read_1d_double(file_id, 'frond_stds', data_dims, rope%frond_stds, error)
    call dset_read_1d_double(file_id, 'water_speeds', data_dims, rope%water_speeds, error)
    call dset_read_1d_double(file_id, 'water_angles', data_dims, rope%water_angles, error)

    call h5fclose_f(file_id, error)
    call h5close_f(error)
  end subroutine hdf_read_rope

  ! subroutine hdf_read_kelp(filename, p_kelp, grid)
  !   integer error
  !   double precision, dimension(:,:,:) :: p_kelp
  !   type(space_angle_grid) grid

  !   integer(hid_t) :: file_id
  !   character(len=*) :: filename
  !   integer(hsize_t), dimension(3) :: data_dims

  !   call h5open_f(error)
  !   call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)

  !   data_dims(1) = grid%x%num
  !   data_dims(2) = grid%y%num
  !   data_dims(3) = grid%z%num

  !   call dset_read_3d_double(file_id, 'p_kelp', data_dims, p_kelp, error)

  !   call h5fclose_f(file_id, error)
  !   call h5close_f(error)
  ! end subroutine hdf_read_kelp

  subroutine hdf_read_frond(filename, frond)
    integer error
    type(frond_shape) frond

    integer(hid_t) :: file_id
    character(len=*) :: filename
    integer(hsize_t), dimension(1) :: data_dims
    double precision, dimension(3) :: frond_arr
    double precision fs, fr, ft

    call h5open_f(error)
    call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)

    data_dims(1) = 2

    call dset_read_1d_double(file_id, 'frond_arr', data_dims, frond_arr, error)
    fs = frond_arr(1)
    fr = frond_arr(2)
    ft = frond_arr(3)

    call frond%set_shape(fs, fr, ft)

    call h5fclose_f(file_id, error)
    call h5close_f(error)
  end subroutine hdf_read_frond

  subroutine hdf_read_params(filename, quadrature_degree)
    integer error
    integer quadrature_degree

    integer(hid_t) :: file_id
    character(len=*) :: filename
    integer(hsize_t), dimension(1) :: data_dims
    double precision, dimension(1) :: param_arr

    call h5open_f(error)
    call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)

    data_dims(1) = 1

    call dset_read_1d_double(file_id, 'param_arr', data_dims, param_arr, error)
    quadrature_degree = param_arr(1)

    call h5fclose_f(file_id, error)
    call h5close_f(error)
  end subroutine hdf_read_params
  !subroutine hdf_read_light()
  !end subroutine hdf_read_light

  !- hdf WRITE kelp

  subroutine hdf_write_kelp(filename, p_kelp, grid)
    integer error
    type(space_angle_grid) grid
    double precision, dimension(:,:,:) :: p_kelp

    integer(hid_t) :: file_id
    character(len=*) :: filename
    integer(hsize_t), dimension(3) :: data_dims

    call h5open_f(error)
    ! Create file
    call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, error)

    data_dims(1) = grid%x%num
    data_dims(2) = grid%y%num
    data_dims(3) = grid%z%num
    call dset_write_3d_double(file_id, 'p_kelp', data_dims, p_kelp, error)

    call h5fclose_f(file_id, error)
    call h5close_f(error)
  end subroutine hdf_write_kelp

  subroutine hdf_write_radiance(filename, radiance, grid)
    integer error
    type(space_angle_grid) grid

    integer(hid_t) :: file_id
    character(len=*) :: filename
    integer(hsize_t), dimension(5) :: data_dims
    double precision, dimension(:,:,:,:,:) :: radiance

    call h5open_f(error)
    ! Create file
    call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, error)

    data_dims(1) = grid%x%num
    data_dims(2) = grid%y%num
    data_dims(3) = grid%z%num
    data_dims(4) = grid%theta%num
    data_dims(5) = grid%phi%num
    call dset_write_5d_double(file_id, 'radiance', data_dims, radiance, error)

    call h5fclose_f(file_id, error)
    call h5close_f(error)
  end subroutine hdf_write_radiance

  subroutine hdf_write_irradiance(filename, irradiance, grid)
    integer error
    type(space_angle_grid) grid

    integer(hid_t) :: file_id
    character(len=*) :: filename
    integer(hsize_t), dimension(3) :: data_dims
    double precision, dimension(:,:,:) :: irradiance

    call h5open_f(error)
    ! Create file
    call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, error)

    data_dims(1) = grid%x%num
    data_dims(2) = grid%y%num
    data_dims(3) = grid%z%num
    call dset_write_3d_double(file_id, 'irradiance', data_dims, irradiance, error)

    call h5fclose_f(file_id, error)
    call h5close_f(error)
  end subroutine hdf_write_irradiance


  ! subroutine hdf_read_irrad()
  ! end subroutine hdf_read_irrad

  ! - grid utils
  ! subroutine grid_read_bounds(bounds, file_id, error)
  !   integer, parameter :: num_bounds = 6
  !   integer(hid_t) :: file_id
  !   integer error
  !   double precision, dimension(num_bounds) :: bounds 
  !   integer(hsize_t), dimension(1) :: data_dims

  !   data_dims(1) = num_bounds

  !   call dset_read_1d_double(file_id, 'bounds', data_dims, bounds, error)
  ! end subroutine grid_read_bounds

  ! subroutine grid_read_nums(nums, file_id, error)
  !   integer, parameter :: num_nums = 5
  !   integer(hid_t) :: file_id
  !   integer error
  !   integer, dimension(num_nums) :: nums 
  !   integer(hsize_t), dimension(1) :: data_dims

  !   data_dims(1) = num_nums

  !   call dset_read_1d_double(file_id, 'nums', data_dims, nums, error)
  ! end subroutine grid_read_nums

  ! subroutine grid_read_spacings(spacings, file_id, error)
  !   integer, parameter :: num_spacings = 5
  !   integer(hid_t) :: file_id
  !   integer error
  !   double precision, dimension(num_spacings) :: spacings 
  !   integer(hsize_t), dimension(1) :: data_dims

  !   array_dim(1) = num_spacings

  !   call dset_read_1d_double(file_id, 'spacings', data_dims, spacings, error)
  ! end subroutine grid_read_spacings


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

  subroutine dset_read_1d_int(file_id, dsetname, data_dims, output, error)
    integer(hid_t) :: file_id, dset_id
    integer(hsize_t), dimension(1) :: data_dims
    character(len=*) :: dsetname
    integer, dimension(:) :: output
    integer error

    call h5dopen_f(file_id, dsetname, dset_id, error)
    call h5dread_f(dset_id, h5t_native_integer, output, data_dims, error)
    call h5dclose_f(dset_id, error)
  end subroutine dset_read_1d_int


  !- dset READ double
  subroutine dset_read_1d_double(file_id, dsetname, data_dims, output, error)
    integer(hid_t) :: file_id, dset_id
    integer(hsize_t), dimension(1) :: data_dims
    character(len=*) :: dsetname
    double precision, dimension(:) :: output
    integer error

    call h5dopen_f(file_id, dsetname, dset_id, error)
    call h5dread_f(dset_id, h5t_native_double, output, data_dims, error)
    call h5dclose_f(dset_id, error)
  end subroutine dset_read_1d_double

  subroutine dset_read_3d_double(file_id, dsetname, data_dims, output, error)
    integer(hid_t) :: file_id, dset_id
    integer(hsize_t), dimension(3) :: data_dims
    character(len=*) :: dsetname
    double precision, dimension(:,:,:) :: output
    integer error

    call h5dopen_f(file_id, dsetname, dset_id, error)
    call h5dread_f(dset_id, h5t_native_double, output, data_dims, error)
    call h5dclose_f(dset_id, error)
  end subroutine dset_read_3d_double

  ! subroutine dset_read_5d_double()
  ! end subroutine dset_read_5d_double


  !- dset WRITE int
  subroutine dset_write_1d_int(file_id, dsetname, data_dims, output, error)
    integer(hid_t) :: file_id, dset_id, dspace_id
    integer, parameter :: rank = 1
    integer(hsize_t), dimension(rank) :: data_dims
    character(len=*) :: dsetname
    integer, dimension(:) :: output
    integer error

    call h5screate_simple_f(rank, data_dims, dspace_id, error)
    ! Create dataset
    call h5dcreate_f(file_id, dsetname, h5t_native_integer, dspace_id, dset_id, error)
    call h5dwrite_f(dset_id, h5t_native_integer, output, data_dims, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)
  end subroutine dset_write_1d_int


  !- dset WRITE double
  subroutine dset_write_1d_double(file_id, dsetname, data_dims, input, error)
    integer(hid_t) :: file_id, dset_id, dspace_id
    integer, parameter :: rank = 1
    integer(hsize_t), dimension(rank) :: data_dims
    character(len=*) :: dsetname
    double precision, dimension(:) :: input
    integer error

    ! Create dataspace
    call h5screate_simple_f(rank, data_dims, dspace_id, error)
    ! Create dataset
    call h5dcreate_f(file_id, dsetname, h5t_native_double, dspace_id, dset_id, error)
    call h5dwrite_f(dset_id, h5t_native_double, input, data_dims, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)
  end subroutine dset_write_1d_double

  subroutine dset_write_3d_double(file_id, dsetname, data_dims, input, error)
    integer(hid_t) :: file_id, dspace_id, dset_id
    integer, parameter :: rank = 3
    integer(hsize_t), dimension(rank) :: data_dims
    character(len=*) :: dsetname
    double precision, dimension(:,:,:) :: input
    integer error

    ! Create dataspace
    call h5screate_simple_f(rank, data_dims, dspace_id, error)
    ! Create dataset
    call h5dcreate_f(file_id, dsetname, h5t_native_double, dspace_id, dset_id, error)

    call h5dwrite_f(dset_id, h5t_native_double, input, data_dims, error)

    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)
  end subroutine dset_write_3d_double

  subroutine dset_write_5d_double(file_id, dsetname, data_dims, input, error)
    integer(hid_t) :: file_id, dspace_id, dset_id
    integer, parameter :: rank = 5
    integer(hsize_t), dimension(rank) :: data_dims
    character(len=*) :: dsetname
    double precision, dimension(:,:,:,:,:) :: input
    integer error
    integer i

    ! Create dataspace
    call h5screate_simple_f(rank, data_dims, dspace_id, error)
    ! Create dataset
    call h5dcreate_f(file_id, dsetname, h5t_native_double, dspace_id, dset_id, error)

    call h5dwrite_f(dset_id, h5t_native_double, input, data_dims, error)

    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)
  end subroutine dset_write_5d_double

  !- hdf READ int

  subroutine hdf_read_int(file_id, dsetname, output, error)

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

  end subroutine hdf_read_int

  !- hdf READ double
subroutine hdf_read_1d_double(filename, dsetname, output, n)
  integer error
  integer n
  character(len=*) filename, dsetname
  integer(hid_t) :: file_id
  integer(hid_t), dimension(1) :: data_dims
  double precision, dimension(n) :: output

  data_dims(1) = n

  call h5open_f(error)
  call h5fopen_f(filename, h5f_acc_rdonly_f, file_id, error)

  call dset_read_1d_double(file_id, dsetname, data_dims, output, error)

  call h5fclose_f(file_id, error)
  call h5close_f(error)
end subroutine hdf_read_1d_double


  !- hdf COO

subroutine read_1d_array_from_coo(filename, arr, n, nnz)
  integer error
  integer n, m, nnz

  integer, dimension(nnz) :: row, col
  double precision, dimension(nnz) :: data
  double precision, dimension(n) :: arr

  integer(hid_t) :: file_id
  integer(hid_t), dimension(1) :: data_dims

  character(len=*) :: filename

  m = 1

  call h5open_f(error)
  call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)

  call dset_read_1d_int(file_id, 'row', data_dims, row, error)
  call dset_read_1d_int(file_id, 'col', data_dims, col, error)
  call dset_read_1d_double(file_id, 'data', data_dims, data, error)

  call h5fclose_f(file_id, error)
  call h5close_f(error)

  call dense_from_coo(row, col, data, arr, n, m, nnz)
end subroutine read_1d_array_from_coo

subroutine read_coo(filename, row, col, data, nnz)
  integer error
  integer n, m, nnz

  integer, dimension(nnz) :: row, col
  double precision, dimension(nnz) :: data

  integer(hid_t), dimension(1) :: data_dims
  integer(hid_t) :: file_id
  character(len=*) :: filename

  data_dims(1) = nnz

  call h5open_f(error)
  call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)

  call dset_read_1d_int(file_id, 'row', data_dims, row, error)
  call dset_read_1d_int(file_id, 'col', data_dims, col, error)
  call dset_read_1d_double(file_id, 'data', data_dims, data, error)

  call h5fclose_f(file_id, error)
  call h5close_f(error)
end subroutine

subroutine write_coo(filename, row, col, data, nnz)
  integer error
  integer n, m, nnz

  integer, dimension(nnz) :: row, col
  double precision, dimension(nnz) :: data

  integer(hid_t), dimension(1) :: data_dims
  integer(hid_t) :: file_id
  character(len=*) :: filename

  data_dims(1) = nnz

  call h5open_f(error)
  call h5fcreate_f(filename, h5f_acc_trunc_f, file_id, error)

  call dset_write_1d_int(file_id, 'row', data_dims, row, error)
  call dset_write_1d_int(file_id, 'col', data_dims, col, error)
  call dset_write_1d_double(file_id, 'data', data_dims, data, error)

  call h5fclose_f(file_id, error)
  call h5close_f(error)
end subroutine write_coo

subroutine read_info(filename, n, m, nnz)
  integer error
  integer n, m, nnz

  integer(hid_t) :: file_id

  character(len=*) :: filename

  call h5open_f(error)
  call h5fopen_f(filename, h5f_acc_rdwr_f, file_id, error)

  call hdf_read_int(file_id, "nnz", nnz, error)
  call hdf_read_int(file_id, "m", m, error)
  call hdf_read_int(file_id, "n", n, error)

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

