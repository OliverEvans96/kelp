module test_grid
  use sag
  implicit none

contains

  function test_2d_angular_integration(func, ntheta, nphi) result(integral)
    type(space_angle_grid) grid
    double precision, external :: func
    integer ntheta, nphi
    double precision integral

    call grid%set_bounds(-1.d0,1.d0,-2.d0,2.d0,-3.d0,3.d0)
    call grid%set_num(4,4,4, ntheta, nphi)
    call grid%init()

    integral = grid%angles%integrate_func(func)

    !if(passed) then
    !   write(*,'(A,I2,A,I2)') '* PASS - test_2d_angular_integration: ', ntheta, 'x', nphi
    !else
    !   write(*,'(A,I2,A,I2)') '* FAIL - test_2d_angular_integration: ', ntheta, 'x', nphi
    !end if

    call grid%deinit()
  end function test_2d_angular_integration

  function test_angle_p_conversions(ntheta, nphi) result(passed)
  integer ntheta, nphi
  logical passed
  integer l, m, p
  type(space_angle_grid) grid

  call grid%set_bounds(-1.d0,1.d0,-2.d0,2.d0,-3.d0,3.d0)
  call grid%set_num(4,4,4, ntheta, nphi)
  call grid%init()

  passed = .true.

  ! Make sure that conversions between l,m and p
  ! are invertible on interior of angular grid
  do l=1, ntheta
     do m=2, nphi-1
        p = grid%angles%phat(l, m)
        if(grid%angles%lhat(p) .ne. l) then
           passed = .false.
           exit
        end if
        if(grid%angles%mhat(p) .ne. m) then
           passed = .false.
           exit
        end if
     end do
     if(.not. passed) then
        exit
     end if
  end do

  !if(passed) then
  !   write(*,'(A,I2,A,I2)') '* PASS - test_p_conversions: ', ntheta, 'x', nphi
  !else
  !   write(*,'(A,I2,A,I2)') '* FAIL - test_p_conversions: ', ntheta, 'x', nphi
  !end if

  call grid%deinit()

  end function test_angle_p_conversions

end module test_grid
