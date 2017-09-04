program test_gmres
  use hdf5_utils

  character(len=24) :: filename = "fortran/download/h5/A.h5"
  double precision, dimension(:,:), allocatable :: arr

  call read_coo(filename, arr)
  
end program test_gmres

