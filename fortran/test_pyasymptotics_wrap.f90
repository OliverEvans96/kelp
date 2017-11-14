program test_pyasymptotics_wrap
use pyasymptotics_wrap

double precision xmin, xmax, ymin, ymax, zmax
integer nx, ny, nz, ntheta, nphi
double precision a_w, a_k, b_w, b_k
integer num_vsf
double precision, dimension(:), allocatable :: vsf_angles, vsf_vals
double precision theta_s, phi_s, max_rad, decay
double precision tol_abs, tol_rel
integer maxiter_inner, maxiter_outer
double precision, dimension(:,:,:), allocatable :: p_kelp
double precision, dimension(:,:,:,:,:), allocatable :: radiance
double precision, dimension(:,:,:), allocatable :: irradiance
integer num_scatters
logical gmres_flag

integer i, j, k, l, m

xmin = -1
xmax = 1
ymin = -1
ymax = 1
zmax = 1

nx = 6
ny = 6
nz = 6
ntheta = 6
nphi = 6

num_vsf = 10

allocate(vsf_angles(num_vsf))
allocate(vsf_vals(num_vsf))
allocate(p_kelp(nx, ny, nz))
allocate(radiance(nx, ny, nz, ntheta, nphi))
allocate(irradiance(nx, ny, nz))

a_w = 1
a_k = 0
b_w = 0
b_k = 0

do i=1, num_vsf
   vsf_angles(i) = pi/num_vsf * (i-1)
   vsf_vals(i) = pi - vsf_angles(i)
end do

theta_s = 0
phi_s = 0
max_rad = 1
decay = 1

tol_abs = 1.d-3
tol_rel = 1.d-3
maxiter_inner = 10
maxiter_outer = 10

num_scatters = 3
gmres_flag = .true.

do i=1, nx
   do j=1, ny
      do k=1, nz
         p_kelp(i,j,k) = 0
         irradiance(i,j,k) = 0
         do l=1, ntheta
            do m=1, nphi
               radiance(i,j,k,l,m) = 0
            end do
         end do
      end do
   end do
end do

call py_calculate_asymptotic_light_field( &
     xmin, xmax, nx, ymin, ymax, ny, zmax, nz, ntheta, nphi, &
     a_w, a_k, b_w, b_k, num_vsf, vsf_angles, vsf_vals, &
     theta_s, phi_s, max_rad, decay, &
     tol_abs, tol_rel, maxiter_inner, maxiter_outer, &
     p_kelp, radiance, irradiance, num_scatters, gmres_flag)

!write(*,*) 'radiance =', radiance

deallocate(vsf_angles)
deallocate(vsf_vals)
deallocate(p_kelp)
deallocate(radiance)
deallocate(irradiance)

end program test_pyasymptotics_wrap

