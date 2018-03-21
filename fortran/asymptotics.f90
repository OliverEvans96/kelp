module asymptotics
  use rte3d
  implicit none
  contains

  subroutine calculate_light_with_scattering(grid, bc, iops, radiance, num_scatters)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:,:) :: radiance
    double precision, dimension(:,:,:,:,:), allocatable :: rad_scatter
    integer num_scatters
    integer i, j, k, l, m
    integer nx, ny, nz, ntheta, nphi
    type(index_list) indices
    double precision surface_val
    double precision bb
    integer n

    double precision, dimension(:), allocatable :: path_spacing
    double precision, dimension(:), allocatable :: outer_path_integrand

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    ntheta = grid%theta%num
    nphi = grid%phi%num

    allocate(path_spacing(nz))
    allocate(outer_path_integrand(nz))

    ! Scattering coefficient
    ! Allowed to vary over depth
    ! Reset radiance
    !write(*,*) 'Before Scattering'
    call calculate_light_before_scattering(grid, bc, iops, radiance, path_spacing, outer_path_integrand)

    if (num_scatters .gt. 0) then
       call calculate_light_after_scattering(grid, bc, iops, radiance, num_scatters, path_spacing, outer_path_integrand)
    end if

    deallocate(path_spacing)
    deallocate(outer_path_integrand)
  end subroutine calculate_light_with_scattering

  subroutine calculate_light_before_scattering(grid, bc, iops, radiance, path_spacing, outer_path_integrand)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:,:) :: radiance
    double precision, dimension(:) :: path_spacing
    double precision, dimension(:) :: outer_path_integrand
    integer i, j, k, l, m
    type(index_list) indices
    double precision surface_val

    ! Downwelling light
   do m=1, grid%phi%num/2
      indices%m = m
      do l=1, grid%theta%num
         indices%l = l
         surface_val = bc%bc_grid(indices%l,indices%m)
         do i=1, grid%x%num
            indices%i = i
            do j=1, grid%y%num
               indices%j = j
               do k=1, grid%z%num
                  indices%k = k
                  radiance(i,j,k,l,m) = absorb_along_path(grid, bc, iops, indices, 1,&
                       surface_val, path_spacing, outer_path_integrand)
                end do
             end do
          end do
       end do
    end do

    ! No upwelling light before scattering
    do m=grid%phi%num/2+1, grid%phi%num
       do l=1, grid%theta%num
          do i=1, grid%x%num
             do j=1, grid%y%num
                do k=1, grid%z%num
                   radiance(i,j,k,l,m) = 0.d0
                end do
             end do
          end do
       end do
    end do
  end subroutine calculate_light_before_scattering

  subroutine calculate_light_after_scattering(grid, bc, iops, radiance, num_scatters, path_spacing, outer_path_integrand)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:,:) :: radiance
    integer num_scatters
    double precision, dimension(:) :: path_spacing
    double precision, dimension(:) :: outer_path_integrand
    double precision, dimension(:,:,:,:,:), allocatable :: rad_scatter
    integer n, k
    double precision bb

    allocate(rad_scatter(grid%x%num, grid%y%num, grid%z%num, grid%theta%num, grid%phi%num))
    rad_scatter = radiance

    ! For now, just use first entry from scattering array everywhere
    bb = iops%scat_water(1)
    write(*,*) 'prescatter = ', rad_scatter(:,:,:,2,2)

    do n=1, num_scatters
       write(*,*) 'scatter #', n
       call scatter(grid, bc, iops, rad_scatter, path_spacing, outer_path_integrand)
       radiance = radiance + bb**n * rad_scatter
       write(*,*) 'rad_scatter = ', rad_scatter(:,:,:,2,2)
    end do

    deallocate(rad_scatter)
  end subroutine calculate_light_after_scattering

  ! Perform one scattering event
  subroutine scatter(grid, bc, iops, rad_scatter, path_spacing, outer_path_integrand)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:,:) :: rad_scatter
    integer i, j, k, l, m, kp
    double precision, dimension(:,:,:,:,:), allocatable :: source, scatter_integral
    double precision, dimension(:,:), allocatable :: scatter_integrand
    double precision, dimension(:) :: path_spacing
    double precision, dimension(:), allocatable :: inner_path_integrand
    double precision, dimension(:) :: outer_path_integrand
    double precision path_source
    integer nx, ny, nz, ntheta, nphi

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    ntheta = grid%theta%num
    nphi = grid%phi%num

    allocate(source(nx, ny, nz, ntheta, nphi))
    allocate(scatter_integral(nx, ny, nz, ntheta, nphi))
    allocate(scatter_integrand(ntheta, nphi))
    ! Allocate largest that will be needed.
    ! Most cases will only use a subset of these arrays.
    allocate(inner_path_integrand(nz))

    write(*,*) 'abs = ', iops%abs_grid

    call calculate_source(grid, iops, rad_scatter, source, scatter_integral, scatter_integrand)
    call advect_light(grid, iops, source, rad_scatter, path_spacing,&
         inner_path_integrand, outer_path_integrand)

    !write(*,*) 'source'
    ! do i=1, nx
    !    do j=1, ny
    !       do k=1, nz
    !         !write(*,'(10E15.3)') sum(rad_scatter(i,j,k,:,:))
    !          !write(*,'(10E15.3)') sum(source(i,j,k,:,:))
    !      end do
    !    end do
    !    !write(*,*) ''
    ! end do

    deallocate(source)
    deallocate(scatter_integral)
    deallocate(scatter_integrand)
    deallocate(inner_path_integrand)
  end subroutine scatter

  ! Calculate source from no-scatter or previous scattering layer
  subroutine calculate_source(grid, iops, rad_scatter, source, scatter_integral, scatter_integrand)
    type(space_angle_grid) grid
    type(optical_properties) iops
    double precision, dimension(:,:,:,:,:) :: rad_scatter, source
    double precision, dimension(:,:,:,:,:) :: scatter_integral
    double precision, dimension(:,:) :: scatter_integrand
    type(index_list) indices
    integer nx, ny, nz, ntheta, nphi
    integer i, j, k, l, m

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    ntheta = grid%theta%num
    nphi = grid%phi%num

    do i=1, nx
       indices%i = i
       do j=1, ny
          indices%j = j
          do k=1, nz
             indices%k = k
             do l=1, ntheta
                indices%l = l
                do m=1, nphi
                   indices%m = m
                   call calculate_scatter_integral(&
                        grid, iops, rad_scatter,&
                        scatter_integral,&
                        scatter_integrand, indices)
                   !write(*,*) 'integrand ', i, j, k, l, m
                   !write(*,*) scatter_integrand(:,:)
                   !write(*,*) 'integral', scatter_integral(i,j,k,l,m)
                   !write(*,*) ''
                end do
             end do
          end do
       end do
    end do

    ! THINK ABOUT WHAT VALUE VSF(0) SHOULD HAVE
    source = -rad_scatter + scatter_integral

  end subroutine calculate_source

  subroutine calculate_scatter_integral(grid, iops, rad_scatter, scatter_integral, scatter_integrand, indices)
    type(space_angle_grid) grid
    type(optical_properties) iops
    double precision, dimension(:,:,:,:,:) :: rad_scatter, scatter_integral
    double precision, dimension(:,:) :: scatter_integrand
    type(index_list) indices
    integer lp, mp

    do lp=1, grid%theta%num
      do mp=1, grid%phi%num
          scatter_integrand(lp, mp) = iops%vsf(&
                indices%l,&
                indices%m,&
                lp,&
                mp) &
            * rad_scatter(&
                indices%i,&
                indices%j,&
                indices%k,&
                lp,&
                mp)
      end do
    end do

    scatter_integral(indices%i,indices%j,indices%k,indices%l,indices%m) = grid%integrate_angle_2d(scatter_integrand)

   if(indices%i .eq. 1 .and. indices%j .eq. 5 .and. indices%k .eq. 3 .and. indices%l .eq. 3 .and. indices%m .eq. 3) then
       write(*,*) ''
       write(*,*) 'theta'
       write(*,*) grid%theta%vals
       write(*,*) 'phi'
       write(*,*) grid%phi%vals
       write(*,*) 'integrand ', indices%i, indices%j, indices%k, indices%l, indices%m
       write(*,*) scatter_integrand(:,:)
       write(*,*) 'integral', scatter_integral(indices%i,indices%j,indices%k,indices%l,indices%m)
       write(*,*) ''
       write(*,*) ''
       write(*,*) ''
    end if

  end subroutine calculate_scatter_integral

  subroutine advect_light(grid, iops, source, rad_scatter, path_spacing,&
       inner_path_integrand, outer_path_integrand)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:,:) :: rad_scatter
    integer i, j, k, l, m, kp
    type(index_list) indices
    double precision, dimension(:,:,:,:,:) :: source

    double precision, dimension(:) :: path_spacing
    double precision, dimension(:) :: inner_path_integrand, outer_path_integrand
    double precision path_source

    integer nx, ny, nz, ntheta, nphi
    double precision x, y, z, theta, phi

    integer kp_start, kp_stop, direction

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    ntheta = grid%theta%num
    nphi = grid%phi%num

    ! Downwelling (project ray to surface)
    kp_start = 1
    direction = 1
    do m=1, nphi / 2
       indices%m = m
       phi = grid%phi%vals(m)
       ! k=1 is already determined by surface BC
       ! Although k=1 would just be a no-op
       do k=2, nz
          indices%k = k
          kp_stop = k
          ! CHECK THIS
          path_spacing(:k) = grid%z%spacing(:k) / grid%phi%cos(m)
          do l=1, ntheta
              indices%l = l
              do i=1, nx
                indices%i = i
                x = grid%x%vals(i)
                do j=1, ny
                   indices%j = j
                   y = grid%y%vals(j)
                   call integrate_ray(&
                        x, y, k, source, grid, iops,&
                        indices, rad_scatter, path_spacing,&
                        inner_path_integrand, outer_path_integrand,&
                        kp_start, kp_stop, direction)
                end do
             end do
          end do
       end do
    end do

    ! Upwelling (ray project to bottom)
    kp_start = nz
    direction = -1
    do m=grid%phi%num/2 + 1, grid%phi%num
       indices%m = m
       phi = grid%phi%vals(m)
       ! k=nz is already determined by bottom BC
       ! Although k=nz would just be a no-op
       do k=1, grid%z%num-1
          indices%k = k
          kp_stop = k
          ! Minus sign since cos(phi) < 0 for upwelling
          ! CHECK THIS
          path_spacing(k:) = - grid%z%spacing(k:) / grid%phi%cos(m)
          do l=1, grid%theta%num
              indices%l = l
              do i=1, grid%x%num
                indices%i = i
                x = grid%x%vals(i)
                do j=1, grid%y%num
                    indices%j = j
                    y = grid%y%vals(j)
                    call integrate_ray(&
                         x, y, k, source, grid, iops,&
                         indices, rad_scatter, path_spacing,&
                         inner_path_integrand, outer_path_integrand,&
                         kp_start, kp_stop, direction)
                end do
             end do
          end do
       end do
    end do
  end subroutine advect_light

  subroutine integrate_ray(x, y, k, source,&
       grid, iops, indices, rad_scatter, path_spacing, inner_path_integrand,&
       outer_path_integrand, kp_start, kp_stop, direction)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:,:) :: rad_scatter
    integer i, j, k, l, m, kp
    type(index_list) indices
    double precision, dimension(:,:,:,:,:) :: source

    double precision, dimension(:) :: path_spacing
    double precision, dimension(:) :: inner_path_integrand, outer_path_integrand
    double precision path_source

    double precision xmin, xmax, ymin, ymax
    integer nx, ny, nz, ntheta, nphi
    double precision x, y, z, theta, phi

    integer kp_start, kp_stop, direction
    integer path_length

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    ntheta = grid%theta%num
    nphi = grid%phi%num

    xmin = grid%x%minval
    xmax = grid%x%maxval
    ymin = grid%y%minval
    ymax = grid%y%maxval

    i = indices%i
    j = indices%j
    k = indices%k
    l = indices%l
    m = indices%m

    path_length = abs(kp_stop - kp_start) + 1

    z = grid%z%vals(k)

    ! Downwelling: kp_start = 1, kp_stop = k, direction = 1
    ! Upwelling: kp_start = nz, kp_stop = k, direction = -1
    do kp = kp_start, kp_stop, direction
       path_source = interpolate_ray_at_depth(&
            x, y, z, kp, nx, ny, nz, xmin, xmax,&
            ymin, ymax, grid%x%vals, grid%y%vals,&
            grid%z%vals, grid%x_factor(l,m),&
            grid%y_factor(l,m), source(:,:,:,l,m))
       !path_integrand(kp) = path_source * absorb_along_path(&
       !     grid, bc, iops, indices, kp)
       ! path_spacing defined in advect_light
       outer_path_integrand(kp) = absorb_along_path(grid, bc, iops, indices,&
            kp, path_source, path_spacing, inner_path_integrand)
    end do
    rad_scatter(i,j,k,l,m) = trap_rule_uneven(path_spacing, inner_path_integrand, path_length)
  end subroutine integrate_ray

  ! Calculate the percent of radiance remaining after passing
  ! through depth layers beginning with kmin, and ending before k.
  ! Path may be downwelling or upwelling.
  ! Absorption is applied at the end of the path, but not the beginning
  ! This is so as not to double absorb upon scattering.
  ! Also, scattering does not occur at the surface for downwelling light
  ! (BC is specified once and for all)
  ! i.e., radiance from nth scattering also considers nth absorption.
  function absorb_along_path(grid, bc, iops, indices,&
       kmin, path_source, path_spacing, inner_path_integrand)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    type(index_list) indices
    double precision x, y, z, theta, phi
    integer nx, ny, nz
    double precision xmin, xmax, ymin, ymax
    double precision x_mod, y_mod

    double precision path_source
    double precision, dimension(:) :: path_spacing, inner_path_integrand
    double precision total_abs

    double precision absorb_along_path

    ! Primed variables
    ! x, y are interpolated
    double precision xp, yp, zp
    ! z is evaluated on grid
    integer kp

    ! Depth at which to begin integration
    integer kmin
    integer path_length, direction, index

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num

    xmin = grid%x%minval
    xmax = grid%x%maxval
    ymin = grid%y%minval
    ymax = grid%y%maxval

    x = grid%x%vals(indices%i)
    y = grid%y%vals(indices%j)
    z = grid%z%vals(indices%k)
    theta = grid%theta%vals(indices%l)
    phi = grid%phi%vals(indices%m)

    ! Number of coefficients to integrate over
    path_length = abs(indices%k - kmin) + 1

    ! Whether the ray is moving up or down
    ! 1 for down, -1 for up, 0 for neither
    direction = sgn_int(indices%k-kmin)

    total_abs = 0

    index = 0
    if(direction .ne. 0) then
       do kp=kmin, indices%k, direction
         index = index + 1
         inner_path_integrand(index) = interpolate_ray_at_depth(&
              x, y, z, kp, nx, ny, nz, xmin, xmax, ymin, ymax,&
              grid%x%vals, grid%y%vals, grid%z%vals,&
              grid%x_factor(indices%l, indices%m),&
              grid%y_factor(indices%l, indices%m),&
              iops%abs_grid)
       end do
    end if

    ! Only pass relevant portion of path_spacing to integrator
    path_spacing(:path_length) = grid%z%spacing(kmin:indices%k:direction) / abs(grid%phi%cos(indices%m))
    total_abs = midpoint_rule_halfends(path_spacing(:path_length), inner_path_integrand, path_length)

    absorb_along_path = path_source * exp(-total_abs)

    ! if(indices%k<3) then
    !    write(*,*) 'ind = ', indices%i, indices%j, indices%k, indices%l, indices%m
    !    write(*,*) 'path_spacing      = ', path_spacing(:path_length)
    !    write(*,*) 'path_abs          = ', inner_path_integrand
    !    write(*,*) 'path_source       = ', path_source
    !    write(*,*) 'absorb_along_path = ', absorb_along_path
    ! end if

  end function absorb_along_path

  ! Given a point (x, y, z) and a direction (theta, phi) such that
  ! x_factor = abs(tan(phi) * cos(theta)),
  ! y_factor = abs(tan(phi) * sin(theta)),
  ! Project a ray to the depth layer kp and interpolate fun_vals
  ! at the unshifted periodic image of the projection.
  function interpolate_ray_at_depth(x, y, z, kp, nx, ny, nz,&
       xmin, xmax, ymin, ymax, x_vals, y_vals, z_vals,&
       x_factor, y_factor, fun_vals)
    double precision x, y, z
    double precision xp, yp, zp
    integer nx, ny, nz, kp
    double precision xmin, xmax, ymin, ymax
    double precision, dimension(nx) :: x_vals
    double precision, dimension(ny) :: y_vals
    double precision, dimension(nz) :: z_vals
    double precision, dimension(:,:,:) :: fun_vals
    double precision x_factor, y_factor
    double precision z_diff, x_mod, y_mod
    double precision interpolate_ray_at_depth
    logical debug_flag

    zp = z_vals(kp)
    z_diff = zp - z
    xp = x + z_diff * x_factor
    yp = y + z_diff * y_factor

    x_mod = shift_mod(xp, xmin, xmax)
    y_mod = shift_mod(yp, ymin, ymax)

    interpolate_ray_at_depth = bilinear_array_periodic(x_mod, y_mod, nx, ny, x_vals, y_vals, fun_vals(:,:,kp))

  end function interpolate_ray_at_depth
end module asymptotics
