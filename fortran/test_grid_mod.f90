module test_grid_mod
  use sag
  implicit none

contains

  subroutine test_2d_angular_integration(func, ntheta, nphi, ans)
    logical passed
    type(space_angle_grid) grid
    double precision, external :: func
    integer ntheta, nphi
    double precision ans

    call grid%set_bounds(-1.d0,1.d0,-2.d0,2.d0,-3.d0,3.d0)
    call grid%set_num(4,4,4, ntheta, nphi)
    call grid%init()

    if(abs(grid%angles%integrate_func(func) - ans) .lt. 1.d-5) then
      passed = .true.
   else
      passed = .false.
   end if

    if(passed) then
       write(*,'(A,I2,A,I2)') '* PASS - test_2d_angular_integration: ', ntheta, 'x', nphi
    else
       write(*,'(A,I2,A,I2)') '* FAIL - test_2d_angular_integration: ', ntheta, 'x', nphi
    end if

    call grid%deinit()
  end subroutine test_2d_angular_integration

  subroutine test_angle_p_conversions(ntheta, nphi)
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

  if(passed) then
     write(*,'(A,I2,A,I2)') '* PASS - test_p_conversions: ', ntheta, 'x', nphi
  else
     write(*,'(A,I2,A,I2)') '* FAIL - test_p_conversions: ', ntheta, 'x', nphi
  end if

  call grid%deinit()

  end subroutine test_angle_p_conversions

end module
