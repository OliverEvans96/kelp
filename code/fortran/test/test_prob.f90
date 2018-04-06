program test_prob
  ! Need VM PDF and Normal CDF
  use prob

  integer n, k
  double precision, parameter :: thetamin=0.d0, thetamax=3.1415926d0, xmin=-1.d0, xmax=1.d0
  double precision theta, x, vm_mu, normal_mu, vm_kappa, normal_sigma
  double precision normal, vm
  double precision dx, dtheta

  ! Number of points to sample, including both endpoints
  n = 10

  normal_mu = 5.d-1
  normal_sigma = 1.d0

  vm_mu = 1.d0
  vm_kappa = 3.d0

  dx = (xmax - xmin) / (n-1)
  dtheta = (thetamax - thetamin) / (n-1)

  write(*,'(A1, A20,A20,A20,A20)') "#", "X","NORMAL","THETA","VON MISES"

  do k = 1, n
     x = xmin + (k-1) * dx
     theta = thetamin + (k-1) * dtheta

     call von_mises_pdf(theta, vm_mu, vm_kappa, vm)
     call normal_cdf(x, normal_mu, normal_sigma, normal)

     write(*,*) x, normal, theta, vm
  end do

  ! Looks good.

end program
