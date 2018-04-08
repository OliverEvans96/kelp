function vsf(costheta)
  double precision a, vsf, costheta
  a = 5.d0
  vsf = exp(-a*(costheta+1)) * a / (exp(a)-exp(-a))
end function vsf

program test_test_asymptotics
  use kelp_context
  use asymptotics
  use light_context
  double precision, external :: vsf
  double precision :: I0, b, theta_s, phi_s, decay
  integer :: nx, ny, nz, ntheta, nphi, nomega
  integer :: num_scatters

  type(space_angle_grid) grid
  type(optical_properties) iops
  type(boundary_condition) bc
  type(light_state) light

  double precision x, y, z
  integer i, j, k

  nx = 10
  ny = 10
  nz = 10
  ntheta = 10
  nphi = 10

  a = 1.0
  b = 0.1

  theta_s = 0.0
  phi_s = 0.0
  I0 = 1.0
  decay = 1.0

  num_scatters = 1

  nomega = ntheta * (nphi-2) + 2

  call grid%set_bounds(0.d0, 1.d0, 0.d0, 1.d0, 0.d0, 1.d0)
  call grid%set_num(nx, ny, nz, ntheta, nphi)
  call grid%init()

  ! INIT IOPS
  iops%num_vsf = 55
  call iops%init(grid)

  do i=1, nx
      x = grid%x%vals(i)
      do j=1, ny
        y = grid%y%vals(j)
        do k=1, nz
            z = grid%z%vals(k)
            iops%abs_grid(i,j,k) = a
        end do
      end do
  end do

  iops%scat = b

  ! Will be called on cos vals [-1, 1]
  call iops%vsf_from_function(vsf)

  ! Straight downward light
  call bc%init(grid, theta_s, phi_s, decay, I0)

  ! Rescale surface radiance to match surface irradiance
  bc%bc_grid = bc%bc_grid * I0 / grid%angles%integrate_points(bc%bc_grid)

  call light%init_grid(grid)
  call calculate_light_with_scattering(grid, bc, iops, light%radiance, num_scatters)
  call light%calculate_irradiance()

  ! write(*,*) 'irrad'
  ! write(*,*) light%irradiance

  call bc%deinit()
  call iops%deinit()
  call light%deinit()
  call grid%deinit()

end program test_test_asymptotics
