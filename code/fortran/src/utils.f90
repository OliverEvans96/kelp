! General utilities which might be useful in other settings
module utils
implicit none

! Constants
double precision, parameter :: pi = 4.D0 * datan(1.D0)

contains

! Determine base directory relative to current directory
! by looking for Makefile, which is in the base dir
! Assuming that this is executed from within the git repo.
function getbasedir()
    implicit none

    ! INPUTS:
    ! Number of paths to check
    integer, parameter :: numpaths = 3
    ! Maximum length of path names
    integer, parameter :: maxlength = numpaths * 2 - 1
    ! Paths to check for Makefile
    character(len=maxlength), parameter, dimension(numpaths) :: check_paths &
            = (/ '.    ', '..   ', '../..' /)
    ! Temporary path string
    character(len=maxlength) tmp_path
    ! Whether Makefile has been found yet
    logical found
    ! Path counter
    integer ii
    ! Lengths of paths
    integer, dimension(numpaths) ::  pathlengths

    ! OUTPUT:
    ! getbasedir - relative path to base directory
    ! Will either return '.', '..', or '../..'
    character(len=maxlength) getbasedir


    ! Determine length of each path
    pathlengths(1) = 1
    do ii = 2, numpaths
        pathlengths(ii) = 2 + 3 * (ii - 2)
    end do

    ! Loop through paths
    do ii = 1, numpaths
        ! Determine this path
        tmp_path = check_paths(ii)

        ! Check whether Makefile is in this directory
        !write(*,*) 'Checking "', tmp_path(1:pathlengths(ii)), '"'
        inquire(file=tmp_path(1:pathlengths(ii)) // '/Makefile', exist=found)
        ! If so, stop. Otherwise, keep looking.
        if(found) then
            getbasedir = tmp_path(1:pathlengths(ii))
            exit
        end if
    end do

    ! If it hasn't been found, then this script was probably called
    ! from outside of the repository.
    if(.not. found) then
        write(*,*) 'BASE DIR NOT FOUND.'
    end if

end function

! Determine array size from min, max and step
! If alignment is off, array will overstep the maximum
function bnd2max(xmin,xmax,dx)
    implicit none

    ! INPUTS:
    ! xmin - minimum x value in array
    ! xmax - maximum x value in array (inclusive)
    ! dx - step size
    double precision, intent(in) :: xmin, xmax, dx

    ! OUTPUT:
    ! step2max - maximum index of array
    integer bnd2max

    ! Calculate array size
    bnd2max = int(ceiling((xmax-xmin)/dx))
end function

! Create array from bounds and number of elements
! xmax is not included in array
function bnd2arr(xmin,xmax,imax)
    implicit none

    ! INPUTS:
    ! xmin - minimum x value in array
    ! xmax - maximum x value in array (exclusive)
    double precision, intent(in) :: xmin, xmax
    ! imax - number of elements in array
    integer imax

    ! OUTPUT:
    ! bnd2arr - array to generate
    double precision, dimension(imax) :: bnd2arr

    ! BODY:

    ! Counter
    integer ii
    ! Step size
    double precision dx

    ! Calculate step size
    dx = (xmax - xmin) / imax

    ! Generate array
    do ii = 1, imax
        bnd2arr(ii) = xmin + (ii-1) * dx
    end do

end function

function mod1(i, n)
  implicit none
  integer i, n, m
  integer mod1

  m = modulo(i, n)

  if(m .eq. 0) then
     mod1 = n
  else
     mod1 = m
  end if

end function mod1

function sgn_int(x)
  integer x, sgn_int
  ! Standard signum function
  sgn_int = sign(1,x) 
  if(x .eq. 0.) sgn_int = 0 
end function sgn_int

function sgn(x)
  double precision x, sgn
  ! Standard signum function
  sgn = sign(1.d0,x) 
  if(x .eq. 0.) sgn = 0 
end function sgn

! Interpolate single point from 1D data
function interp(x0,xx,yy,nn)
    implicit none

    ! INPUTS:
    ! x0 - x value at which to interpolate
    double precision, intent(in) :: x0
    ! xx - ordered x values at which y data is sampled
    ! yy - corresponding y values to interpolate
    double precision, dimension (nn), intent(in) :: xx,yy
    ! nn - length of data
    integer, intent(in) :: nn

    ! OUTPUT:
    ! interp - interpolated y value
    double precision interp

    ! BODY:

    ! Index of lower-adjacent data (xx(i) < x0 < xx(i+1))
    integer ii
    ! Slope of liine between (xx(ii),yy(ii)) and (xx(ii+1),yy(ii+1))
    double precision mm

    ! If out of bounds, then return endpoint value
    if (x0 < xx(1)) then
       interp = yy(1)
    else if (x0 > xx(nn)) then
       interp = yy(nn)
    else

      ! Determine ii
      do ii = 1, nn
          if (xx(ii) > x0) then
              ! We've now gone one index too far.
              exit
          end if
      end do

      ! Determine whether we're on the right endpoint
      if(ii-1 < nn) then
          ! If this is a legitimate interpolation, then
          ! subtract since we went one index too far
          ii = ii - 1

          ! Calculate slope
          mm = (yy(ii+1) - yy(ii)) / (xx(ii+1) - xx(ii))

          ! Return interpolated value
          interp = yy(ii) + mm * (x0 - xx(ii))
      else
          ! If we're actually interpolating the right endpoint,
          ! then just return it.
          interp = yy(nn)
      end if

   end if

end function

! Calculate unshifted position of periodic image
! Assuming xmin, xmax are extreme attainable values of x
function shift_mod(x, xmin, xmax)
  double precision x, xmin, xmax
  double precision mod_part, shift_mod
  mod_part = mod(x-xmin, xmax-xmin)
  if(mod_part .ge. 0) then
     ! In this case, mod_part is distance between image & lower bound
     shift_mod = xmin + mod_part
  else
     ! In this case, mod_part is distance between image & upper bound
     shift_mod = xmax + mod_part
  endif
end function shift_mod

! Bilinear interpolation on evenly spaced 2D grid
! Assume upper endpoint is not included and is identical
! to the lower endpoint, which is included.
function bilinear_array_periodic(x, y, nx, ny, x_vals, y_vals, fun_vals)
  implicit none
  double precision x, y
  integer nx, ny
  double precision, dimension(:) :: x_vals, y_vals
  double precision, dimension(:,:) :: fun_vals

  double precision dx, dy, xmin, ymin
  integer i0, j0, i1, j1
  double precision x0, x1, y0, y1
  double precision z00, z10, z01, z11

  double precision bilinear_array_periodic

  xmin = x_vals(1)
  ymin = y_vals(1)
  dx = x_vals(2) - x_vals(1)
  dy = y_vals(2) - y_vals(1)

  ! Add 1 for one-indexing
  i0 = int(floor((x-xmin)/dx))+1
  j0 = int(floor((y-ymin)/dy))+1

  x0 = x_vals(i0)
  y0 = y_vals(j0)

  ! Periodic wrap
  if(i0 .lt. nx) then
     i1 = i0 + 1
     x1 = x_vals(i1)
  else
     i1 = 1
     x1 = x_vals(nx) + dx
  endif

  if(j0 .lt. ny) then
     j1 = j0 + 1
     y1 = y_vals(j1)
  else
     j1 = 1
     y1 = y_vals(ny) + dy
  endif

  z00 = fun_vals(i0,j0)
  z10 = fun_vals(i1,j0)
  z01 = fun_vals(i0,j1)
  z11 = fun_vals(i1,j1)

  bilinear_array_periodic = bilinear(x, y, x0, y0, x1, y1, z00, z01, z10, z11)
end function bilinear_array_periodic

! Bilinear interpolation on evenly spaced 2D grid
! Assume upper and lower endpoints are included
function bilinear_array(x, y, x_vals, y_vals, fun_vals)
  implicit none
  double precision x, y
  double precision, dimension(:) :: x_vals, y_vals
  double precision, dimension(:,:) :: fun_vals

  double precision dx, dy, xmin, ymin
  integer i0, j0, i1, j1
  double precision x0, x1, y0, y1
  double precision z00, z10, z01, z11

  double precision bilinear_array

  xmin = x_vals(1)
  ymin = y_vals(1)
  dx = x_vals(2) - x_vals(1)
  dy = y_vals(2) - y_vals(1)

  ! Add 1 for one-indexing
  i0 = int(floor((x-xmin)/dx))+1
  j0 = int(floor((y-ymin)/dy))+1
  i1 = i0 + 1
  j1 = j0 + 1

  ! Bounds checking
  ! if(i0 .lt. 1) then
  !    i0 = 1
  !    i1 = 1
  ! else if(i1 .gt. nx) then
  !    i0 = nx
  !    i1 = nx
  ! endif
  ! if(j0 .lt. 1) then
  !    j0 = 1
  !    j1 = 1
  ! else if(j1 .gt. ny) then
  !    j0 = ny
  !    j1 = ny
  ! endif

  x0 = x_vals(i0)
  x1 = x_vals(i1)
  y0 = y_vals(j0)
  y1 = y_vals(j1)

  z00 = fun_vals(i0,j0)
  z10 = fun_vals(i1,j0)
  z01 = fun_vals(i0,j1)
  z11 = fun_vals(i1,j1)

  bilinear_array = bilinear(x, y, x0, y0, x1, y1, z00, z01, z10, z11)
end function bilinear_array

! ilinear interpolation of a function of two variables
! over a rectangle of points.
! Weight each point by the area of the sub-rectangle involving
! the point (x,y) and the point diagonally across the rectangle
function bilinear(x, y, x0, y0, x1, y1, z00, z01, z10, z11)
  implicit none
  double precision x, y
  double precision x0, y0, x1, y1, z00, z01, z10, z11
  double precision a, b, c, d
  double precision bilinear

  a = (x-x0)*(y-y0)
  b = (x1-x)*(y-y0)
  c = (x-x0)*(y1-y)
  d = (x1-x)*(y1-y)

  bilinear = (a*z11 + b*z01 + c*z10 + d*z00) / (a + b + c + d)
end function bilinear

! Integrate using left endpoint rule
! Assuming the right endpoint is not included in arr
function lep_rule(arr,dx,nn)
    implicit none

    ! INPUTS:
    ! arr - array to integrate
    double precision, dimension(nn) :: arr
    ! dx - array spacing (mesh size)
    double precision dx
    ! nn - length of arr
    integer, intent(in) :: nn

    ! OUTPUT:
    ! lep_rule - integral w/ left endpoint rule
    double precision lep_rule

    ! BODY:

    ! Counter
    integer ii

    ! Set output to zero
    lep_rule = 0.0d0

    ! Accumulate integral
    do ii = 1, nn
        lep_rule = lep_rule + arr(ii) * dx
    end do

end function

! Integrate using trapezoid rule
! Assuming both endpoints are included in arr
function trap_rule_dx(arr, dx, nn)
  implicit none
  double precision, dimension(nn) :: arr
  double precision dx
  integer ii, nn
  double precision trap_rule_dx

  trap_rule_dx = 0.0d0

  do ii=1, nn-1
     trap_rule_dx = trap_rule_dx + 0.5d0 * dx * (arr(ii) + arr(ii+1))
  end do

end function trap_rule_dx

! Integrate using trapezoid rule
! Assuming both endpoints are included in arr
function trap_rule_uneven(xx, yy, nn)
  implicit none
  double precision, dimension(nn) :: xx
  double precision, dimension(nn) :: yy
  integer ii, nn
  double precision trap_rule_uneven

  trap_rule_uneven = 0.0d0

  do ii=1, nn-1
     trap_rule_uneven = trap_rule_uneven + 0.5d0 * (xx(ii+1)-xx(ii)) * (yy(ii) + yy(ii+1))
  end do
end function trap_rule_uneven

function trap_rule_dx_uneven(dx, yy, nn)
  implicit none
  double precision, dimension(nn-1) :: dx
  double precision, dimension(nn) :: yy
  integer ii, nn
  double precision trap_rule_dx_uneven

  trap_rule_dx_uneven = 0.0d0

  do ii=1, nn-1
     trap_rule_dx_uneven = trap_rule_dx_uneven + 0.5d0 * dx(ii) * (yy(ii) + yy(ii+1))
  end do
end function trap_rule_dx_uneven

! Integrate using midpoint rule
! First and last bins, only use inner half
function midpoint_rule_halfends(dx, yy, nn) result(integral)
  implicit none
  integer ii, nn
  double precision, dimension(nn) :: dx, yy
  double precision integral

  if(nn > 1) then
    integral = .5d0 * (dx(1)*yy(1) + dx(nn)*yy(nn))

    do ii=2, nn-1
      integral = integral + dx(ii)*yy(ii)
    end do
  else
    integral = 0.d0
  end if
end function midpoint_rule_halfends

! Normalize 1D array and return integral w/ left endpoint rule
function normalize_dx(arr,dx,nn)
    implicit none

    ! INPUTS:
    ! arr - array to normalize
    double precision, dimension(nn) :: arr
    ! dx - array spacing (mesh size)
    double precision dx
    ! nn - length of arr
    integer, intent(in) :: nn

    ! OUTPUT:
    ! normalize - integral before normalization (left endpoint rule)
    double precision normalize_dx

    ! BODY:

    ! Calculate integral
    normalize_dx = lep_rule(arr,dx,nn)

    ! Normalize array
    arr = arr / normalize_dx

  end function normalize_dx

! Normalize 1D unevenly-spaced array and
! return integral w/ trapezoid rule
! Will not be quite accurate if rightmost endpoint is not included
! (Very small for VSF, so not a big deal there)
! Modifies yy in place
function normalize_uneven(xx, yy, nn) result(norm)
  implicit none

  ! INPUTS:
  ! xx, yy - array values of data to normalize
  double precision, dimension(nn) :: xx, yy
  ! nn - length of arr
  integer, intent(in) :: nn

  ! OUTPUT:
  ! normalize - integral before normalization (left endpoint rule)
  double precision norm

  ! BODY:

  ! Calculate integral
  ! PERHAPS WE SHOULD USE TRAPEZOID RULE
  norm = trap_rule_uneven(xx, yy, nn)

  ! Normalize array
  yy(:) = yy(:) / norm

end function normalize_uneven

! Read 2D array from file
function read_array(filename,fmtstr,nn,mm,skiplines_in)
    implicit none

    ! INPUTS:
    ! filename - path to file to be read
    ! fmtstr - input format (no parentheses, don't specify columns)
    ! e.g. 'E10.2', not '(2E10.2)'
    character(len=*), intent(in) :: filename, fmtstr
    ! nn - Number of data rows in file
    ! mm - number of data columns in file
    integer, intent(in) :: nn, mm
    ! skiplines - optional - number of lines to skip from header
    integer, optional :: skiplines_in
    integer skiplines

    ! OUTPUT:
    double precision, dimension(nn,mm) :: read_array

    ! BODY:

    ! Row counter
    integer ii
    ! File unit number
    integer, parameter :: un = 10
    ! Final format to use
    character(len=256) finfmt

    ! Generate final format string
    write(finfmt,'(A,I1,A,A)') '(', mm, fmtstr, ')'

    ! Print message
    !write(*,*) 'Reading data from "', trim(filename), '"'
    !write(*,*) 'using format "', trim(finfmt), '"'

    ! Open file
    open(unit=un, file=trim(filename), status='old', form='formatted')

    ! Skip lines if desired
    if(present(skiplines_in)) then
        skiplines = skiplines_in
        do ii = 1, skiplines
            ! Read without variable ignores the line
            read(un,*)
        end do
    else
        skiplines = 0
    end if

    ! Loop through lines
    do ii = 1, nn
        ! Read one row at a time
        read(unit=un, fmt=trim(finfmt)) read_array(ii,:)
    end do

    ! Close file
    close(unit=un)

end function

! Print 2D array to stdout
subroutine print_int_array(arr,nn,mm,fmtstr_in)
  implicit none

  ! INPUTS:
  ! arr - array to print
  integer, dimension (nn,mm), intent(in) :: arr
  ! nn - number of data rows in file
  ! nn - number of data columns in file
  integer, intent(in) :: nn, mm
  ! fmtstr - output format (no parentheses, don't specify columns)
  ! e.g. 'E10.2', not '(2E10.2)'
  character(len=*), optional :: fmtstr_in
  character(len=256) fmtstr

  ! NO OUTPUTS

  ! BODY

  ! Row counter
  integer ii
  ! Final format to use
  character(len=256) finfmt

  ! Determine string format
  if(present(fmtstr_in)) then
     fmtstr = fmtstr_in
  else
     fmtstr = 'I10'
  end if

  ! Generate final format string
  write(finfmt,'(A,I4,A,A)') '(', mm, trim(fmtstr), ')'

  ! Loop through rows
  do ii = 1, nn
     ! Print one row at a time
     write(*,finfmt) arr(ii,:)
  end do

  ! Print blank line after
  write(*,*) ' '

end subroutine print_int_array

subroutine print_array(arr,nn,mm,fmtstr_in)
    implicit none

    ! INPUTS:
    ! arr - array to print
    double precision, dimension (nn,mm), intent(in) :: arr
    ! nn - number of data rows in file
    ! nn - number of data columns in file
    integer, intent(in) :: nn, mm
    ! fmtstr - output format (no parentheses, don't specify columns)
    ! e.g. 'E10.2', not '(2E10.2)'
    character(len=*), optional :: fmtstr_in
    character(len=256) fmtstr

    ! NO OUTPUTS

    ! BODY

    ! Row counter
    integer ii
    ! Final format to use
    character(len=256) finfmt

    ! Determine string format
    if(present(fmtstr_in)) then
        fmtstr = fmtstr_in
    else
        fmtstr = 'ES10.2'
    end if

    ! Generate final format string
    write(finfmt,'(A,I4,A,A)') '(', mm, trim(fmtstr), ')'

    ! Loop through rows
    do ii = 1, nn
        ! Include row number
        !write(*,'(I10)', advance='no') ii
        ! Print one row at a time
        write(*,finfmt) arr(ii,:)
    end do

    ! Print blank line after
    write(*,*) ' '

end subroutine

! Write 1D array to file
subroutine write_vec(arr,nn,filename,fmtstr_in)
    implicit none

    ! INPUTS:
    ! arr - array to print
    double precision, dimension (nn), intent(in) :: arr
    ! nn - number of data rows in file
    ! nn - number of data columns in file
    integer, intent(in) :: nn
    ! filename - file to write to
    character(len=*) filename
    ! fmtstr - output format (no parentheses, don't specify columns)
    ! e.g. 'E10.2', not '(2E10.2)'
    character(len=*), optional :: fmtstr_in
    character(len=256) fmtstr

    ! NO OUTPUTS

    ! BODY

    ! Row counter
    integer ii
    ! Final format to use
    character(len=256) finfmt
    ! Dummy file unit to use
    integer, parameter :: un = 20

    ! Open file for writing
    open(unit=un, file=trim(filename), status='replace', form='formatted')

    ! Determine string format
    if(present(fmtstr_in)) then
        fmtstr = fmtstr_in
    else
        fmtstr = 'E10.2'
    end if

    ! Generate final format string
    write(finfmt,'(A,A,A)') '(', trim(fmtstr), ')'

    ! Loop through rows
    do ii = 1, nn
        ! Print entry per row
        write(un,finfmt) arr(ii)
    end do

    ! Close file
    close(unit=un)

end subroutine

! Write 2D array to file
subroutine write_array(arr,nn,mm,filename,fmtstr_in)
    implicit none

    ! INPUTS:
    ! arr - array to print
    double precision, dimension (nn,mm), intent(in) :: arr
    ! nn - number of data rows in file
    ! nn - number of data columns in file
    integer, intent(in) :: nn, mm
    ! filename - file to write to
    character(len=*) filename
    ! fmtstr - output format (no parentheses, don't specify columns)
    ! e.g. 'E10.2', not '(2E10.2)'
    character(len=*), optional :: fmtstr_in
    character(len=256) fmtstr

    ! NO OUTPUTS

    ! BODY

    ! Row counter
    integer ii
    ! Final format to use
    character(len=256) finfmt
    ! Dummy file unit to use
    integer, parameter :: un = 20

    ! Open file for writing
    open(unit=un, file=trim(filename), status='replace', form='formatted')

    ! Determine string format
    if(present(fmtstr_in)) then
        fmtstr = fmtstr_in
    else
        fmtstr = 'E10.2'
    end if

    ! Generate final format string
    write(finfmt,'(A,I4,A,A)') '(', mm, trim(fmtstr), ')'

    ! Loop through rows
    do ii = 1, nn
        ! Print one row at a time
        write(un,finfmt) arr(ii,:)
    end do

    ! Close file
    close(unit=un)

end subroutine

subroutine zeros(x, n)
  implicit none
  integer n, i
  double precision, dimension(n) :: x

  do i=1, n
     x(i) = 0
  end do
end subroutine zeros

end module
