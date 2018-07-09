program nc_pykelp3d_wrap
use pykelp3d_wrap
use netcdf
implicit none

double precision xmin, xmax, ymin, ymax, zmin, zmax
double precision rope_spacing
integer nx, ny, nz, ntheta, nphi, nomega
double precision a_w, a_k, b
double precision fr, fs, ft
double precision num_dens
integer num_vsf
double precision, dimension(:), allocatable :: vsf_angles, vsf_vals
double precision theta_s, phi_s, max_rad, I0, decay
double precision, dimension(:,:,:), allocatable :: p_kelp
double precision, dimension(:,:,:,:), allocatable :: radiance
double precision, dimension(:,:,:), allocatable :: irradiance
real, dimension(:), allocatable :: avg_irrad, perceived_irrad
integer num_scatters
integer fd_flag_int
logical fd_flag
character(len=256) :: lis_opts
integer lis_iter
double precision :: lis_time, lis_resid

double precision, dimension(:), allocatable :: frond_lengths
double precision, dimension(:), allocatable :: frond_stds
double precision, dimension(:), allocatable :: water_speeds
double precision, dimension(:), allocatable :: water_angles
double precision, dimension(:), allocatable :: num_fronds
character(len=256) kelp_dist
double precision max_length

character(len=256) :: nc_path
integer ncid, varid
integer nc_err
double precision dz

integer i, j, k, p


! Path to .nc file as command line argument
if(command_argument_count() .ne. 1) then
   write(*,*) "Please supplt .nc file as the only command line argument."
   call exit(1)
end if

call get_command_argument(1, nc_path)

write(*,*) 'Reading netCDF: ', trim(nc_path)
nc_err = nf90_open(nc_path, nf90_nowrite, ncid)

nc_err = nf90_inq_varid(ncid, "rope_spacing", varid)
nc_err = nf90_get_var(ncid, varid, rope_spacing)
write(*,*) "READ: rope_spacing = ", rope_spacing

nc_err = nf90_inq_varid(ncid, "nx", varid)
nc_err = nf90_get_var(ncid, varid, nx)
write(*,*) "READ: nx = ", nx

nc_err = nf90_inq_varid(ncid, "ny", varid)
nc_err = nf90_get_var(ncid, varid, ny)
write(*,*) "READ: ny = ", ny

nc_err = nf90_inq_varid(ncid, "nz", varid)
nc_err = nf90_get_var(ncid, varid, nz)
write(*,*) "READ: nz = ", nz

nc_err = nf90_inq_varid(ncid, "zmax", varid)
nc_err = nf90_get_var(ncid, varid, zmax)
write(*,*) "READ: zmax = ", zmax

nc_err = nf90_inq_varid(ncid, "ntheta", varid)
nc_err = nf90_get_var(ncid, varid, ntheta)
write(*,*) "READ: ntheta = ", ntheta

nc_err = nf90_inq_varid(ncid, "nphi", varid)
nc_err = nf90_get_var(ncid, varid, nphi)
write(*,*) "READ: nphi = ", nphi

nc_err = nf90_inq_varid(ncid, "a_w", varid)
nc_err = nf90_get_var(ncid, varid, a_w)
write(*,*) "READ: a_w = ", a_w

nc_err = nf90_inq_varid(ncid, "a_k", varid)
nc_err = nf90_get_var(ncid, varid, a_k)
write(*,*) "READ: a_k = ", a_k

nc_err = nf90_inq_varid(ncid, "b", varid)
nc_err = nf90_get_var(ncid, varid, b)
write(*,*) "READ: b = ", b

nc_err = nf90_inq_varid(ncid, "theta_s", varid)
nc_err = nf90_get_var(ncid, varid, theta_s)
write(*,*) "READ: theta_s = ", theta_s

nc_err = nf90_inq_varid(ncid, "phi_s", varid)
nc_err = nf90_get_var(ncid, varid, phi_s)
write(*,*) "READ: phi_s = ", phi_s

nc_err = nf90_inq_varid(ncid, "I0", varid)
nc_err = nf90_get_var(ncid, varid, I0)
write(*,*) "READ: I0 = ", I0

nc_err = nf90_inq_varid(ncid, "decay", varid)
nc_err = nf90_get_var(ncid, varid, decay)
write(*,*) "READ: decay = ", decay

nc_err = nf90_inq_varid(ncid, "fs", varid)
nc_err = nf90_get_var(ncid, varid, fs)
write(*,*) "READ: fs = ", fs

nc_err = nf90_inq_varid(ncid, "fr", varid)
nc_err = nf90_get_var(ncid, varid, fr)
write(*,*) "READ: fr = ", fr

nc_err = nf90_inq_varid(ncid, "ft", varid)
nc_err = nf90_get_var(ncid, varid, ft)
write(*,*) "READ: ft = ", ft

nc_err = nf90_inq_varid(ncid, "kelp_dist", varid)
nc_err = nf90_get_var(ncid, varid, kelp_dist, count=(/ 9 /))
write(*,*) "READ: kelp_dist = ", kelp_dist

nc_err = nf90_inq_varid(ncid, "max_length", varid)
nc_err = nf90_get_var(ncid, varid, max_length)
write(*,*) "READ: max_length = ", max_length

nc_err = nf90_inq_varid(ncid, "num_dens", varid)
nc_err = nf90_get_var(ncid, varid, num_dens)
write(*,*) "READ: num_dens = ", num_dens

nc_err = nf90_inq_varid(ncid, "num_scatters", varid)
nc_err = nf90_get_var(ncid, varid, num_scatters)
write(*,*) "READ: num_scatters = ", num_scatters

nc_err = nf90_inq_varid(ncid, "fd_flag", varid)
nc_err = nf90_get_var(ncid, varid, fd_flag_int)
write(*,*) "READ: fd_flag = ", fd_flag

nc_err = nf90_inq_varid(ncid, "lis_opts", varid)
nc_err = nf90_get_var(ncid, varid, lis_opts)
write(*,*) "READ: lis_opts = ", lis_opts

! Close dataset
nc_err = nf90_close(ncid)

! Construct options for `calculate_light_field` from .nc values

! fd_flag (convert integer -> logical)
if(fd_flag_int .eq. 0) then
   fd_flag = .false.
else
   fd_flag = .true.
end if

! x-y limits
xmin = -rope_spacing / 2
xmax = rope_spacing / 2
ymin = -rope_spacing / 2
ymax = rope_spacing / 2
zmin = 0

! ntheta, nphi specified.
nomega = ntheta*(nphi-2)+2

allocate(frond_lengths(nz))
allocate(frond_stds(nz))
allocate(water_speeds(nz))
allocate(water_angles(nz))
allocate(num_fronds(nz))

allocate(vsf_angles(num_vsf))
allocate(vsf_vals(num_vsf))
allocate(p_kelp(nx, ny, nz))
allocate(radiance(nx, ny, nz, nomega))
allocate(irradiance(nx, ny, nz))
allocate(avg_irrad(nz))
allocate(perceived_irrad(nz))

! Mimicking what kelp_compute.py does
num_vsf = ntheta
do i=1, num_vsf
   vsf_angles(i) = dble(i-1) * pi / dble(num_vsf - 1)
   vsf_vals(i) = 1.d0/(4.d0*pi)
end do

! Initialize everything to zero
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

! TODO: Do I need to convert here?
max_rad = I0

dz = (zmax-zmin)/nz
do k=1, nz
  num_fronds(k) = num_dens * dz
end do

call get_kelp_dist(kelp_dist, max_length, zmin, zmax, nz, frond_lengths, frond_stds, water_speeds, water_angles)

call gen_kelp(xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz, &
     frond_lengths, frond_stds, num_fronds, water_speeds, water_angles, &
     fs, fr, ft, p_kelp)

call calculate_light_field( &
     xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz, ntheta, nphi, &
     a_w, a_k, b, num_vsf, vsf_angles, vsf_vals, &
     theta_s, phi_s, max_rad, decay, &
     p_kelp, radiance, irradiance, avg_irrad, perceived_irrad, &
     num_scatters, fd_flag, lis_opts, lis_iter, lis_time, lis_resid)

!write(*,*) 'radiance =', radiance

! write(*,*) 'rad min: ', minval(radiance)
! write(*,*) 'rad max: ', maxval(radiance)
! write(*,*) 'rad mean: ', sum(radiance)/size(radiance)
! write(*,*)
! write(*,*) 'irrad min: ', minval(irradiance)
! write(*,*) 'irrad max: ', maxval(irradiance)
! write(*,*) 'irrad mean: ', sum(irradiance)/size(irradiance)

deallocate(vsf_angles)
deallocate(vsf_vals)
deallocate(p_kelp)
deallocate(radiance)
deallocate(irradiance)
deallocate(avg_irrad)
deallocate(perceived_irrad)
deallocate(frond_lengths)
deallocate(frond_stds)
deallocate(water_angles)
deallocate(water_speeds)
deallocate(num_fronds)
! TODO: Double-check deallocation

end program nc_pykelp3d_wrap

subroutine get_kelp_dist(kelp_dist, max_length, zmin, zmax, nz, frond_lengths, frond_stds, water_speeds, water_angles)
  use utils

  character(len=*), intent(in) :: kelp_dist
  double precision, intent(in) :: max_length
  double precision, intent(in) :: zmin, zmax
  integer, intent(in) :: nz
  double precision, dimension(nz), intent(out) :: frond_lengths, frond_stds
  double precision, dimension(nz), intent(out) :: water_speeds, water_angles

  integer k
  double precision maxval
  double precision dz
  ! TODO: Have to allocate?
  double precision :: z

  maxval = 12 * exp(-2.d0) + 0.5

  dz = (zmax - zmin) / dble(nz)
  do k=1, nz
  end do

  ! TODO: Not sure about string comparison
  do k=1, nz
     z = (k-0.5d0) * dz
     if(kelp_dist .eq. "top-heavy") then
        frond_lengths(k) = max_length * (3.d0 * z**2 * exp(-z) + 0.5d0) / maxval
     else if(kelp_dist .eq. "bottom-heavy") then
        frond_lengths(k) = max_length * (3.d0 * (zmax-z)**2 * exp(z-zmax) + 0.5) / maxval
     else if(kelp_dist .eq. "uniform") then
        frond_lengths(k) = max_length
     else if(kelp_dist .eq. "none") then
        frond_lengths(k) = 0.d0
     end if

     frond_stds(k) = 0.d0
     water_speeds(k) = 0.5d0
     water_angles(k) = 2*pi / zmax * (z-zmin)
  end do
end subroutine get_kelp_dist
