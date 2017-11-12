module asymptotics
  use rte3d

  subroutine calculate_light_before_scattering(grid, bc, iops, radiance)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:,:) :: radiance
    integer i, j, k, l, m
    type(index_list) indices
    double precision surface_val


    do l=1, grid%theta%num
       indices%l = l
       do m=1, grid%phi%num
          indices%m = m
          surface_val = bc%bc_grid(indices%l,indices%m)
          do i=1, grid%x%num
             indices%i = i
             do j=1, grid%y%num
                indices%j = j
                do k=1, grid%z%num
                   indices%k = k
                   radiance(i,j,k,l,m) = absorb_along_path(grid, bc, iops, indices, 1)
                end do
             end do
          end do
       end do
    end do
  end subroutine calculate_light_before_scattering

  ! Perform one scattering event
  subroutine scatter(grid, bc, iops, rad_prescatter, rad_postscatter)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:,:) :: rad_prescatter, rad_postscatter
    integer i, j, k, l, m, kp
    type(index_list) indices
    double precision surface_val
    double precision dpath

    double precision, dimension(:), allocatable :: integrand

    ! Allocate largest that will be needed.
    ! Most cases will only use a subset of the array.
    allocate(integrand(grid%z%num))

    calculate_source(grid, iops, rad_prescatter, source)

    do indices%i=1, grid%x%num
       do indices%j=1, grid%y%num
          do indices%k=1, grid%z%num
             do indices%l=1, grid%theta%num
                do indices%m=1, grid%phi%num

                   do kp = 1, k
                      integrand(kp) = source(i,j,k,l,m) * absorb_along_path(grid, bc, iops, indices, kp)
                   end do

                   dpath = grid%z%spacing * sec(phi)
                   rad_postscatter(i,j,k,l,m) = trap_rule(integrand, dpath, k)

                end do
             end do
          end do
       end do
    end do
  end subroutine scatter

  ! Calculate source from no-scatter or previous scattering layer
  subroutine calculate_source(grid, iops, rad_prescatter, source)
    type(space_angle_grid) grid
    type(optical_properties) iops
    double precision, dimension(:,:,:,:,:) :: rad_prescatter, source
    double precision, dimension(:,:,:,:,:), allocatable :: scatter_integral
    integer nx, ny, nz, ntheta, nphi

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    ntheta = grid%theta%num
    nphi = grid%phi%num

    allocate(scatter_integral(nx, ny, nz, ntheta, nphi))

    call calculate_scatter_integral(grid, iops, rad_prescatter, scatter_integral)
    source = -rad_prescatter + scatter_integral

  end subroutine calculate_source

  subroutine calculate_scatter_integral(grid, iops, rad_prescatter, scatter_integral
    type(space_angle_grid) grid
    type(optical_properties) iops
    double precision, dimension(:,:,:,:,:) :: rad_prescatter, scatter_integral

    double precision, dimension(:,:), allocatable :: integrand
    double precision angle_diff

    integer i, j, k, l, m
    integer lp, mp

    allocate(integrand(grid%theta%num, grid%phi%num))

    do i=1, grid%x%num
       do j=1, grid%y%num
          do l=1, grid%z%num
             do l=1, grid%theta%num
                do m=1, grid%phi%num
                   do lp=1, grid%theta%num
                      do mp=1, grid%phi%num
                         integrand(lp, mp) = iops%vsf(l,m,lp,mp) * rad_prescatter(i,j,k,lp,mp)
                      end do
                   end do

                   scatter_integral(i,j,k,l,m) = grid.integrate_angle_2d(integrand)
                end do
             end do
          end do
       end do
    end do

    deallocate(integrand)

  end subroutine calculate_scatter_integral

  function absorb_along_path(grid, bc, iops, indices, kmin)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    type(index_list) indices
    double precision x, y, z, theta, phi
    integer nx, ny

    double precision, dimension(:), allocatable :: abs_along_path
    double precision dpath, total_abs, absorb_along_path

    ! Primed variables
    ! x, y are interpolated
    double precision xp, yp, zp
    ! z is evaluated on grid
    integer kp

    ! Interpolation variables
    integer i0, j0,
    double precision x0, y0, x1, y1

    ! Depth at which to begin integration
    integer kmin

    nx = grid%x%num
    ny = grid%y%num

    x = grid%x%vals(indices%i)
    y = grid%y%vals(indices%j)
    z = grid%z%vals(indices%k)
    theta = grid%theta%vals(indices%l)
    phi = grid%phi%vals(indices%m)

    allocate(abs_coef_along_path(indices%k))

    do kp=kmin, indices%k
       zp = grid%z%vals(kp)
       xp = x + zp * tan(phi) * cos(theta)
       yp = y + zp * tan(phi) * sin(theta)

       abs_coef_along_path(kp) = bilinear_array(xp, yp, nx, ny grid%x%vals, grid%y%vals, iops%abs_coef(:,:,kp))
    end do

    dpath = grid%z%spacing * sec(phi)
    total_abs = trap_rule(abs_coef_along_path, dpath, indices%k)

    deallocate(abs_coef_along_path)

    absorb_along_path = bc%bc_grid(indices%l, indices%m) * exp(-total_abs)

  end function absorb_along_path
end module asymptotics
