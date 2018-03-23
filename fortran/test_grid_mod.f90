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

end module
