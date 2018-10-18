module type_consts
  use iso_fortran_env
  ! Type of integer to use for matrix index
  ! int32 max size is 2^31-1 ~= 2e9
  integer, parameter :: index_kind = int64
  ! LIS_INTEGER length defined by #define LONG__LONG
  ! in light_context.f90 and rte_sparse_matrices.f90
end module type_consts

