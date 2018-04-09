function vsf_func(costheta)
  use utils
  double precision a, vsf_func, costheta
  ! vsf_func decay parameter
  a = 10.d0
  vsf_func = exp(a*(costheta)) * a / (2*pi*(exp(a)-exp(-a)))
end function vsf_func

function a_func(x, y, z)
  double precision x, y, z, a_func
  a_func = x + y + z ** 2
end function a_func

program test_test_asymptotics
  !use kelp_context
  !use asymptotics
  !use light_context
  use test_asymptotics
  implicit none

  double precision, external :: vsf_func, a_func
  double precision :: I0, b, theta_s, phi_s, decay
  integer :: nx, ny, nz, ntheta, nphi, nomega
  integer :: num_scatters

  double precision, dimension(:,:,:,:), allocatable :: rad
  double precision, dimension(:,:,:), allocatable :: irrad

  nx = 10
  ny = 10
  nz = 10
  ntheta = 10
  nphi = 10

  b = 0.5

  theta_s = 0.0
  phi_s = 0.0
  I0 = 1.0
  decay = 1.0

  num_scatters = 1

  nomega = ntheta * (nphi-2) + 2

  allocate(rad(nx,ny,nz,nomega))
  allocate(irrad(nx,ny,nz))

  call test_asymptotics_3d(&
       I0, theta_s, phi_s, decay, &
       a_func, b, vsf_func, &
       nx, ny, nz, ntheta, nphi, &
       num_scatters, rad, irrad)

  deallocate(irrad)
  deallocate(rad)

end program test_test_asymptotics
