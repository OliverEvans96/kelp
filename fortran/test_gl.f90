program test_gl
  use fastgl
  integer n, k
  double precision theta, weight, x

  n = 10

  do k = 1, n
    call glpair(n, k, theta, weight, x)
    write(*,*) theta, weight, x
  end do

  ! Looks good. x = cos(theta).

end program test_gl
