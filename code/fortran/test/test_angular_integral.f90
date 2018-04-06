program test_angular_integral
  use sag

  type(space_angle_grid) grid
  integer i, j, k, l, m
  double precision theta, phi
  double precision, dimension(:,:), allocatable :: func
  double precision integral

  grid%x%minval = 0
  grid%x%maxval = 1
  grid%x%num = 10

  grid%y%minval = 0
  grid%y%maxval = 1
  grid%y%num = 10

  grid%z%minval = 0
  grid%z%maxval = 1
  grid%z%num = 10

  grid%theta%num = 10
  grid%phi%num = 10

  call grid%init()

  allocate(func(grid%theta%num, grid%phi%num))

  do l=1, grid%theta%num
     theta = grid%theta%vals(l)
     do m=1, grid%phi%num
        phi = grid%phi%vals(m)
        func(l,m) = sin(theta/2) * sin(phi)
     end do
  end do

  integral = grid%integrate_angle_2d(func)

  write(*,*) 'integral = ', integral

  deallocate(func)
end program test_angular_integral

