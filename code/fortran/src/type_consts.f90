module type_consts
  use iso_fortran_env
  ! Type of integer to use for matrix index
  ! int32 max size is 2^31-1 ~= 2e9
  integer, parameter :: index_kind = int64
end module type_consts

