module prob
contains
subroutine von_mises_pdf ( x, a, b, pdf )

!*****************************************************************************80
!
!! VON_MISES_PDF evaluates the von Mises PDF.
!
!  Discussion:
!
!    PDF(A,B;X) = EXP ( B * COS ( X - A ) ) / ( 2 * PI * I0(B) )
!
!    where:
!
!      I0(*) is the modified Bessel function of the first
!      kind of order 0.
!
!    The von Mises distribution for points on the unit circle is
!    analogous to the normal distribution of points on a line.
!    The variable X is interpreted as a deviation from the angle A,
!    with B controlling the amount of dispersion.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jerry Banks, editor,
!    Handbook of Simulation,
!    Engineering and Management Press Books, 1998, page 160.
!
!    Donald Best, Nicholas Fisher,
!    Efficient Simulation of the von Mises Distribution,
!    Applied Statistics,
!    Volume 28, Number 2, pages 152-157.
!
!    Merran Evans, Nicholas Hastings, Brian Peacock,
!    Statistical Distributions,
!    Wiley, 2000,
!    LC: QA273.6.E92, pages 189-191.
!
!    Kanti Mardia, Peter Jupp,
!    Directional Statistics,
!    Wiley, 2000,
!    LC: QA276.M335
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!    A - PI <= X <= A + PI.
!
!    Input, real ( kind = 8 ) A, a parameter of the PDF.
!    A is the preferred direction, in radians.
!    -PI <= A <= PI.
!
!    Input, real ( kind = 8 ) B, a parameter of the PDF.
!    B measures the "concentration" of the distribution around the
!    angle A.  B = 0 corresponds to a uniform distribution
!    (no concentration).  Higher values of B cause greater concentration
!    of probability near A.
!    0.0 <= B.
!
!    Output, real ( kind = 8 ) PDF, the value of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
!123  real ( kind = 8 ) bessel_i0
  real ( kind = 8 ) pdf
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  if ( x < a - r8_pi ) then
    pdf = 0.0D+00
  else if ( x <= a + r8_pi ) then
    pdf = exp ( b * cos ( x - a ) ) / ( 2.0D+00 * r8_pi * bessel_i0 ( b ) )
  else
    pdf = 0.0D+00
  end if

  return
end subroutine von_mises_pdf


subroutine normal_01_cdf ( x, cdf )

!*****************************************************************************80
!
!! NORMAL_01_CDF evaluates the Normal 01 CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    AG Adams,
!    Algorithm 39,
!    Areas Under the Normal Curve,
!    Computer Journal,
!    Volume 12, pages 197-198, 1969.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ), parameter :: a1 = 0.398942280444D+00
  real ( kind = 8 ), parameter :: a2 = 0.399903438504D+00
  real ( kind = 8 ), parameter :: a3 = 5.75885480458D+00
  real ( kind = 8 ), parameter :: a4 = 29.8213557808D+00
  real ( kind = 8 ), parameter :: a5 = 2.62433121679D+00
  real ( kind = 8 ), parameter :: a6 = 48.6959930692D+00
  real ( kind = 8 ), parameter :: a7 = 5.92885724438D+00
  real ( kind = 8 ), parameter :: b0 = 0.398942280385D+00
  real ( kind = 8 ), parameter :: b1 = 3.8052D-08
  real ( kind = 8 ), parameter :: b2 = 1.00000615302D+00
  real ( kind = 8 ), parameter :: b3 = 3.98064794D-04
  real ( kind = 8 ), parameter :: b4 = 1.98615381364D+00
  real ( kind = 8 ), parameter :: b5 = 0.151679116635D+00
  real ( kind = 8 ), parameter :: b6 = 5.29330324926D+00
  real ( kind = 8 ), parameter :: b7 = 4.8385912808D+00
  real ( kind = 8 ), parameter :: b8 = 15.1508972451D+00
  real ( kind = 8 ), parameter :: b9 = 0.742380924027D+00
  real ( kind = 8 ), parameter :: b10 = 30.789933034D+00
  real ( kind = 8 ), parameter :: b11 = 3.99019417011D+00
  real ( kind = 8 ) cdf
  real ( kind = 8 ) q
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  |X| <= 1.28.
!
  if ( abs ( x ) <= 1.28D+00 ) then

    y = 0.5D+00 * x * x

    q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
      + a6 / ( y + a7 ) ) ) )
!
!  1.28 < |X| <= 12.7
!
  else if ( abs ( x ) <= 12.7D+00 ) then

    y = 0.5D+00 * x * x

    q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
      + b2 / ( abs ( x ) + b3 &
      + b4 / ( abs ( x ) - b5 &
      + b6 / ( abs ( x ) + b7 &
      - b8 / ( abs ( x ) + b9 &
      + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
!
!  12.7 < |X|
!
  else

    q = 0.0D+00

  end if
!
!  Take account of negative X.
!
  if ( x < 0.0D+00 ) then
    cdf = q
  else
    cdf = 1.0D+00 - q
  end if

  return
end subroutine normal_01_cdf


subroutine normal_cdf ( x, a, b, cdf )

!*****************************************************************************80
!
!! NORMAL_CDF evaluates the Normal CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Input, real ( kind = 8 ) A, B, the parameters of the PDF.
!    0.0 < B.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  y = ( x - a ) / b

  call normal_01_cdf ( y, cdf )

  return
end subroutine normal_cdf


function bessel_i0 ( arg )

!*****************************************************************************80
!
!! BESSEL_I0 evaluates the modified Bessel function I0(X).
!
!  Discussion:
!
!    The main computation evaluates slightly modified forms of
!    minimax approximations generated by Blair and Edwards, Chalk
!    River (Atomic Energy of Canada Limited) Report AECL-4928,
!    October, 1974.  This transportable program is patterned after
!    the machine dependent FUNPACK packet NATSI0, but cannot match
!    that version for efficiency or accuracy.  This version uses
!    rational functions that theoretically approximate I-SUB-0(X)
!    to at least 18 significant decimal digits.
!
!  Machine dependent constants:
!
!    beta   = Radix for the floating-point system
!    maxexp = Smallest power of beta that overflows
!    XMAX =   Largest argument acceptable to BESI0;  Solution to
!             equation:
!               W(X) * (1+1/(8*X)+9/(128*X^2) = beta^maxexp
!             where  W(X) = EXP(X)/sqrt(2*PI*X)
!
!    Approximate values for some important machines are:
!
!                             beta       maxexp       XMAX
!
!    CRAY-1        (S.P.)       2         8191       5682.810
!    Cyber 180/855
!      under NOS   (S.P.)       2         1070        745.893
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)       2          128         91.900
!    IEEE (IBM/XT,
!      SUN, etc.)  (D.P.)       2         1024        713.986
!    IBM 3033      (D.P.)      16           63        178.182
!    VAX           (S.P.)       2          127         91.203
!    VAX D-Format  (D.P.)       2          127         91.203
!    VAX G-Format  (D.P.)       2         1023        713.293
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 October 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.
!
!    Output, real ( kind = 8 ) BESSEL_I0, the value of the modified
!    Bessel function of the first kind.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ) bessel_i0
  real ( kind = 8 ), parameter :: exp40 = 2.353852668370199854D+17
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter, dimension ( 15 ) :: p = (/ &
    -5.2487866627945699800D-18, &
    -1.5982226675653184646D-14, &
    -2.6843448573468483278D-11, &
    -3.0517226450451067446D-08, &
    -2.5172644670688975051D-05, &
    -1.5453977791786851041D-02, &
    -7.0935347449210549190D+00, &
    -2.4125195876041896775D+03, &
    -5.9545626019847898221D+05, &
    -1.0313066708737980747D+08, &
    -1.1912746104985237192D+10, &
    -8.4925101247114157499D+11, &
    -3.2940087627407749166D+13, &
    -5.5050369673018427753D+14, &
    -2.2335582639474375249D+15 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: pp = (/ &
    -3.9843750000000000000D-01, &
     2.9205384596336793945D+00, &
    -2.4708469169133954315D+00, &
     4.7914889422856814203D-01, &
    -3.7384991926068969150D-03, &
    -2.6801520353328635310D-03, &
     9.9168777670983678974D-05, &
    -2.1877128189032726730D-06 /)
  real ( kind = 8 ), parameter, dimension ( 5 ) :: q = (/ &
    -3.7277560179962773046D+03, &
     6.5158506418655165707D+06, &
    -6.5626560740833869295D+09, &
     3.7604188704092954661D+12, &
    -9.7087946179594019126D+14 /)
  real ( kind = 8 ), parameter, dimension ( 7 ) :: qq = (/ &
    -3.1446690275135491500D+01, &
     8.5539563258012929600D+01, &
    -6.0228002066743340583D+01, &
     1.3982595353892851542D+01, &
    -1.1151759188741312645D+00, &
     3.2547697594819615062D-02, &
    -5.5194330231005480228D-04 /)
  real ( kind = 8 ), parameter :: rec15 = 6.6666666666666666666D-02
  real ( kind = 8 ) sump
  real ( kind = 8 ) sumq
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xmax = 91.9D+00
  real ( kind = 8 ) xx

  x = abs ( arg )

  if ( x < epsilon ( arg ) ) then
    value = 1.0D+00
  else if ( x < 15.0D+00 ) then
!
!  EPSILON ( ARG ) <= ABS(ARG) < 15.0D+00
!
    xx = x * x
    sump = p(1)
    do i = 2, 15
      sump = sump * xx + p(i)
    end do

    xx = xx - 225.0D+00
    sumq = (((( &
           xx + q(1) ) &
         * xx + q(2) ) &
         * xx + q(3) ) &
         * xx + q(4) ) &
         * xx + q(5)

    value = sump / sumq

  else if ( 15.0D+00 <= x ) then

    if ( xmax < x ) then
      value = huge ( value )
    else
!
!  15.0D+00 <= ABS(ARG)
!
      xx = 1.0D+00 / x - rec15

      sump = ((((((  &
                  pp(1) &
           * xx + pp(2) ) &
           * xx + pp(3) ) &
           * xx + pp(4) ) &
           * xx + pp(5) ) &
           * xx + pp(6) ) &
           * xx + pp(7) ) &
           * xx + pp(8)

      sumq = (((((( &
             xx + qq(1) ) &
           * xx + qq(2) ) &
           * xx + qq(3) ) &
           * xx + qq(4) ) &
           * xx + qq(5) ) &
           * xx + qq(6) ) &
           * xx + qq(7)

      value = sump / sumq
!
!  Calculation reformulated to avoid premature overflow.
!
      if ( x <= xmax - 15.0D+00 ) then
        a = exp ( x )
        b = 1.0D+00
      else
        a = exp ( x - 40.0D+00 )
        b = exp40
      end if

      value = ( ( value * a - pp(1) * a ) / sqrt ( x ) ) * b

    end if

  end if

  bessel_i0 = value

  return
end function bessel_i0

end module