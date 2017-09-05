program test_gmres

  implicit none

  call test_easy()
  call test_real()

end program test_gmres

subroutine test_real()
  use hdf5_utils
  use utils
  use mgmres

  implicit none

  integer error
  character(len=24) :: a_fn = "fortran/hdf5test/A.h5"
  character(len=24) :: rhs_fn = "fortran/hdf5test/b.h5"
  character(len=24) :: x_fn = "fortran/hdf5test/x.h5"
  character(len=27) :: test_fn = "fortran/hdf5test/test.h5"
  integer, dimension(:), allocatable :: row, col
  double precision, dimension(:), allocatable :: data, x, x_true, rhs, test
  integer n, m, nnz_a, nnz_rhs, nnz_x_true
  integer n_test, m_test, nnz_test

  ! Max # of outer iterations
  integer itr_max
  ! Max # of inner iterations
  integer mr
  double precision tol_abs, tol_rel

  itr_max = 10
  mr = 10
  tol_abs = 1.d-5
  tol_rel = 1.d-5

  call read_info(rhs_fn, n, m, nnz_rhs)
  call read_info(x_fn, n, m, nnz_x_true)
  call read_info(a_fn, n, m, nnz_a)
  call read_info(test_fn, n_test, m_test, nnz_test)
  write(*,*) "nnz = ", nnz_test
  write(*,*) "n = ", n_test
  write(*,*) "m = ", m_test

  allocate(row(nnz_a))
  allocate(col(nnz_a))
  allocate(data(nnz_a))

  allocate(x(n))
  allocate(x_true(n))
  allocate(rhs(n))
  allocate(test(n_test))

  call zeros(x, n)

  call read_coo(a_fn, row, col, data, nnz_a)
  call read_1d_array_from_coo(rhs_fn, rhs, n, nnz_rhs)
  call read_1d_array_from_coo(x_fn, x_true, n, nnz_x_true)

  call mgmres_st(n, nnz_a, row, col, data, x, rhs, itr_max, mr, tol_abs, tol_rel)

  write(*,*) 'ERROR = ', sum((x-x_true)**2)
  write(*,*) 'If this is small, then it worked!'

  deallocate(row)
  deallocate(col)
  deallocate(data)
  deallocate(x)
  deallocate(x_true)
  deallocate(rhs)
  deallocate(test)

end subroutine test_real


subroutine test_easy()
  use hdf5_utils
  use utils
  integer error
  integer, dimension(:), allocatable :: row, col
  double precision, dimension(:), allocatable :: data, x, x_true, rhs
  integer, parameter :: n = 10, nnz = 100

  ! Max # of outer iterations
  integer itr_max
  ! Max # of inner iterations
  integer mr
  double precision tol_abs, tol_rel

  itr_max = 100
  mr = 10
  tol_abs = 1.d-3
  tol_rel = 1.d-3

  allocate(row(nnz))
  allocate(col(nnz))
  allocate(data(nnz))

  allocate(x(n))
  allocate(x_true(n))
  allocate(rhs(n))

  call zeros(x, n)
  
  call create_test_equation(row, col, data, rhs, x_true)
  
  call mgmres_st(n, nnz, row, col, data, x, rhs, itr_max, mr, tol_abs, tol_rel)

  call print_array(x, n, 1)

  write(*,*) 'True solution:'
  call print_array(x_true, n, 1)
  ! Correct solution is [1 2 3 4 5 6 7 8 9 10]

  deallocate(row)
  deallocate(col)
  deallocate(data)
  deallocate(x)
  deallocate(x_true)
  deallocate(rhs)
end subroutine test_easy


subroutine zeros(x, n)
  integer n, i
  double precision, dimension(n) :: x

  do i=1, n
     x(i) = 0
  end do
end subroutine zeros

subroutine create_test_equation(row, col, data, b, x_true)
  use utils
  use hdf5_utils
  implicit none
  integer, parameter :: n = 10
  double precision, dimension(n,n) :: a
  double precision, dimension(n) :: b, x_true
  integer, dimension(n*n) :: row, col
  double precision, dimension(n*n) :: data
  integer i, j

  b(1) = -23.d0
  b(2) = -6.d0
  b(3) = -9.d0
  b(4) = -12.d0
  b(5) = -15.d0
  b(6) = -18.d0
  b(7) = -21.d0
  b(8) = -24.d0
  b(9) = -27.d0
  b(10) = -10.d0

  do i=1, n
     x_true(i) = i

     do j=1, n
        if(i .eq. j) then
           a(i,j) = 1.d0
        else if((j .eq. mod1(i-1, n)) .or. (j .eq. mod1(i+1, n))) then
           a(i,j) = -2.d0
        else
           a(i,j) = 0.d0
        end if
     end do
  end do

  write(*,*) 'A'
  call print_array(a, n, n)
  write(*,*) 'rhs'
  call print_array(b, n, 1)

  call coo_from_dense(a, n, n, row, col, data)
end subroutine create_test_equation


