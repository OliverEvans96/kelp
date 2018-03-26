program test_traverse_prog
  use test_asymptotics
  implicit none

  integer i, j, k, l, m, p
  integer ip, jp, kp, pp
  integer num_cells
  double precision xmin, xmax, ymin, ymax, zmin, zmax
  integer nx, ny, nz, ntheta, nphi, nomega
  double precision, dimension(:), allocatable :: path_lengths, s_array, ds, a_tilde, gn
  double precision, dimension(:,:,:,:), allocatable :: rad_scatter

  integer max_cells

  xmin = 0
  ymin = 0
  zmin = 0
  xmax = 10
  ymax = 10
  zmax = 10

  nx = 10
  ny = 10
  nz = 10
  ntheta = 10
  nphi = 10

  nomega = ntheta*(nphi-2)+2

  i = 4
  j = 3
  k = 6
  l = 5
  m = 2

  max_cells = 100

  allocate(path_lengths(max_cells))
  allocate(ds(max_cells))
  allocate(a_tilde(max_cells))
  allocate(gn(max_cells))
  allocate(rad_scatter(nx,ny,nz,nomega))

  do ip=1, max_cells
     path_lengths(ip) = 0
     ds(ip) = 0
     a_tilde(ip) = 0
     gn(ip) = 0
  end do

  do ip=1, nx
     do jp=1, ny
        do kp=1, nz
           do pp=1, nomega
              rad_scatter(ip,jp,kp,pp) = 0
           end do
        end do
     end do
  end do

  call test_traverse(&
       xmin, xmax, ymin, ymax, zmin, zmax,&
       nx, ny, nz, ntheta, nphi,&
       i, j, k, l, m,&
       path_lengths, ds, a_tilde, gn, rad_scatter, num_cells)

  deallocate(path_lengths)
  deallocate(ds)
  deallocate(a_tilde)
  deallocate(gn)
  deallocate(rad_scatter)
end program test_traverse_prog
