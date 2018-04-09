program test_pyasymptotics_wrap
use pykelp3d_wrap

double precision xmin, xmax, ymin, ymax, zmin, zmax
integer nx, ny, nz, ntheta, nphi, nomega
double precision a_w, a_k, b
integer num_vsf
double precision, dimension(:), allocatable :: vsf_angles, vsf_vals
double precision theta_s, phi_s, max_rad, decay
double precision tol_abs, tol_rel
integer maxiter_inner, maxiter_outer
double precision, dimension(:,:,:), allocatable :: p_kelp
double precision, dimension(:,:,:,:), allocatable :: radiance
double precision, dimension(:,:,:), allocatable :: irradiance
integer num_scatters
logical gmres_flag

integer i, j, k, p

xmin = -1
xmax = 1
ymin = -1
ymax = 1
zmin = 0
zmax = 2

nx = 10
ny = 10
nz = 10
ntheta = 10
nphi = 10

nomega = ntheta*(nphi-2)+2

num_vsf = 30

allocate(vsf_angles(num_vsf))
allocate(vsf_vals(num_vsf))
allocate(p_kelp(nx, ny, nz))
allocate(radiance(nx, ny, nz, nomega))
allocate(irradiance(nx, ny, nz))

a_w = 1.d0
a_k = 0.d0
b = 0.2d0

do i=1, num_vsf
   vsf_angles(i) = pi/num_vsf * (i-1)
   vsf_vals(i) = pi - vsf_angles(i)
end do

theta_s = 0.d0
phi_s = 0.d0
max_rad = 1.d0
decay = 1.d0

tol_abs = 1.d-3
tol_rel = 1.d-3
maxiter_inner = 100
maxiter_outer = 100

num_scatters = 1
gmres_flag = .true.

do i=1, nx
   do j=1, ny
      do k=1, nz
         p_kelp(i,j,k) = 0
         irradiance(i,j,k) = 0
         do p=1, nomega
               radiance(i,j,k,p) = 0
         end do
      end do
   end do
end do

call calculate_light_field( &
     xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz, ntheta, nphi, &
     a_w, a_k, b, num_vsf, vsf_angles, vsf_vals, &
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

