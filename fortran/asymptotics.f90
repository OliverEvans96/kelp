module asymptotics
  use rte3d
  implicit none
  contains

  subroutine calculate_light_with_scattering(grid, bc, iops, radiance, num_scatters)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: radiance
    double precision, dimension(:,:,:,:), allocatable :: rad_scatter
    integer num_scatters
    integer nx, ny, nz, nomega
    type(index_list) indices
    double precision surface_val
    double precision bb
    integer n
    integer max_cells

    double precision, dimension(:), allocatable :: path_length, path_spacing, a_tilde, gn

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    nomega = grid%angles%nomega

    max_cells = calculate_max_cells(grid)

    allocate(path_length(max_cells))
    allocate(path_spacing(max_cells))
    allocate(a_tilde(max_cells))
    allocate(gn(max_cells))


    ! Scattering coefficient
    ! Allowed to vary over depth
    ! Reset radiance
    !write(*,*) 'Before Scattering'
    call calculate_light_before_scattering(grid, bc, iops, radiance, path_spacing, a_tilde)

    if (num_scatters .gt. 0) then
       call calculate_light_after_scattering(grid, bc, iops, radiance, num_scatters, path_length, path_spacing, a_tilde, gn)
    end if

    allocate(path_length(max_cells))
    allocate(path_spacing(max_cells))
    deallocate(path_spacing)
    deallocate(a_tilde)
  end subroutine calculate_light_with_scattering

  subroutine calculate_light_before_scattering(grid, bc, iops, radiance, path_spacing, a_tilde)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: radiance
    double precision, dimension(:) :: path_spacing, a_tilde
    integer i, j, k, p
    type(index_list) indices
    double precision surface_val

    ! Downwelling light
    ! TODO: Check indices
    do p=1, grid%angles%nomega/2
       indices%p = p
       surface_val = bc%bc_grid(indices%p)
       do i=1, grid%x%num
          indices%i = i
          do j=1, grid%y%num
            indices%j = j
            do k=1, grid%z%num
                indices%k = k
                radiance(i,j,k,p) = absorb_along_path(grid, bc, iops, indices, 1,&
                    surface_val, path_spacing, a_tilde)
            end do
          end do
       end do
     end do

     ! No upwelling light before scattering
     do p = grid%angles%nomega/2+1, grid%angles%nomega
        do i=1, grid%x%num
          do j=1, grid%y%num
              do k=1, grid%z%num
                radiance(i,j,k,p) = 0.d0
              end do
          end do
        end do
     end do
  end subroutine calculate_light_before_scattering

  subroutine calculate_light_after_scattering(grid, bc, iops, radiance, num_scatters, path_length, path_spacing, a_tilde, gn)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: radiance
    integer num_scatters
    double precision, dimension(:) :: path_length, path_spacing, a_tilde, gn
    double precision, dimension(:,:,:,:), allocatable :: rad_scatter
    integer n, k
    double precision bb

    allocate(rad_scatter(grid%x%num, grid%y%num, grid%z%num, grid%angles%nomega))
    rad_scatter = radiance

    ! For now, just use first entry from scattering array everywhere
    bb = iops%scat_water(1)
    write(*,*) 'prescatter = ', rad_scatter(:,:,:,2)

    do n=1, num_scatters
       write(*,*) 'scatter #', n
       call scatter(grid, bc, iops, rad_scatter, path_length, path_spacing, a_tilde, gn)
       radiance = radiance + bb**n * rad_scatter
       write(*,*) 'rad_scatter = ', rad_scatter(:,:,:,2)
    end do

    deallocate(rad_scatter)
  end subroutine calculate_light_after_scattering

  ! Perform one scattering event
  subroutine scatter(grid, bc, iops, rad_scatter, path_length, path_spacing, a_tilde, gn)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: rad_scatter
    double precision, dimension(:,:,:,:), allocatable :: source, scatter_integral
    double precision, dimension(:), allocatable :: scatter_integrand
    double precision, dimension(:) :: path_length, path_spacing, a_tilde, gn
    double precision path_source
    integer max_cells
    integer nx, ny, nz, nomega

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    nomega = grid%angles%nomega


    allocate(source(nx, ny, nz, nomega))
    allocate(scatter_integral(nx, ny, nz, nomega))
    allocate(scatter_integrand(nomega))
    ! Allocate largest that will be needed.
    ! Most cases will only use a subset of these arrays.

    !write(*,*) 'abs = ', iops%abs_grid

    call calculate_source(grid, iops, rad_scatter, source, scatter_integral, scatter_integrand)
    call advect_light(grid, iops, source, rad_scatter, path_length, path_spacing, a_tilde, gn)

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
  end subroutine scatter

  ! Calculate source from no-scatter or previous scattering layer
  subroutine calculate_source(grid, iops, rad_scatter, source, scatter_integral, scatter_integrand)
    type(space_angle_grid) grid
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: rad_scatter, source
    double precision, dimension(:,:,:,:) :: scatter_integral
    double precision, dimension(:) :: scatter_integrand
    type(index_list) indices
    integer nx, ny, nz, nomega
    integer i, j, k, p

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    nomega = grid%angles%nomega

    do i=1, nx
       indices%i = i
       do j=1, ny
          indices%j = j
          do k=1, nz
             indices%k = k
             ! TODO: DOUBLE CHECK INDICES
             do p=1, nomega
                indices%p = p
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

    ! TODO: THINK ABOUT WHAT VALUE VSF(0) SHOULD HAVE
    source = -rad_scatter + scatter_integral

  end subroutine calculate_source

  subroutine calculate_scatter_integral(grid, iops, rad_scatter, scatter_integral, scatter_integrand, indices)
    type(space_angle_grid) grid
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: rad_scatter, scatter_integral
    double precision, dimension(:) :: scatter_integrand
    type(index_list) indices
    integer pp

    ! TODO: Exclude forward direction from scattering integral?
    do pp=1, grid%angles%nomega
       scatter_integrand(pp) = iops%vsf(&
             indices%p,&
             pp) &
         * rad_scatter(&
             indices%i,&
             indices%j,&
             indices%k,&
             pp)
    end do

    scatter_integral(indices%i,indices%j,indices%k,indices%p) = grid%angles%integrate_points(scatter_integrand)

   ! if(indices%i .eq. 1 .and. indices%j .eq. 5 .and. indices%k .eq. 3 .and. indices%p .eq. 13) then
   !     write(*,*) ''
   !     write(*,*) 'theta'
   !     write(*,*) grid%angles%theta
   !     write(*,*) 'phi'
   !     write(*,*) grid%angles%phi
   !     write(*,*) 'integrand ', indices%i, indices%j, indices%k, indices%p
   !     write(*,*) scatter_integrand(:)
   !     write(*,*) 'integral', scatter_integral(indices%i,indices%j,indices%k,indices%p)
   !     write(*,*) ''
   !     write(*,*) ''
   !     write(*,*) ''
   !  end if

  end subroutine calculate_scatter_integral

  subroutine advect_light(grid, iops, source, rad_scatter, path_length, path_spacing, a_tilde, gn)
    type(space_angle_grid) grid
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: rad_scatter
    double precision, dimension(:) :: path_length, path_spacing, a_tilde, gn
    double precision path_source

    integer i, j, k, p, kp
    type(index_list) indices
    double precision, dimension(:,:,:,:) :: source


    integer nx, ny, nz, nomega
    double precision x, y, z, theta, phi

    integer kp_start, kp_stop, direction

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    nomega = grid%angles%nomega

    ! Downwelling (project ray to surface)
    kp_start = 1
    direction = 1
    do p=1, nomega/2
       indices%p = p
       phi = grid%angles%phi_p(p)
       ! k=1 is already determined by surface BC
       ! Although k=1 would just be a no-op
       do k=2, nz
          indices%k = k
          kp_stop = k
          ! TODO: CHECK THIS
          path_spacing(:k) = grid%z%spacing(:k) / grid%angles%cos_phi_p(p)
          ! TODO: CHECK INDICES
          do i=1, nx
            indices%i = i
            x = grid%x%vals(i)
            do j=1, ny
                indices%j = j
                y = grid%y%vals(j)
                call integrate_ray(grid, iops, source,&
                     rad_scatter, path_length, path_spacing,&
                     a_tilde, gn, i, j, k, p)
            end do
          end do
       end do
    end do

    ! Upwelling (ray project to bottom)
    kp_start = nz
    direction = -1
    do p=nomega+1, nomega
       indices%p = p
       ! k=nz is already determined by bottom BC
       ! Although k=nz would just be a no-op
       do k=1, grid%z%num-1
          indices%k = k
          kp_stop = k
          ! Minus sign since cos(phi) < 0 for upwelling
          ! TODO: CHECK THIS
          path_spacing(k:) = - grid%z%spacing(k:) / grid%angles%cos_phi_p(p)
          do i=1, grid%x%num
            indices%i = i
            x = grid%x%vals(i)
            do j=1, grid%y%num
                indices%j = j
                y = grid%y%vals(j)
                call integrate_ray(grid, iops, source,&
                     rad_scatter, path_length, path_spacing,&
                     a_tilde, gn, i, j, k, p)
             end do
          end do
       end do
    end do
  end subroutine advect_light

  ! New algorithm, double integral over piecewise-constant 1d funcs
  subroutine integrate_ray(grid, iops, source, rad_scatter, path_length, path_spacing, a_tilde, gn, i, j, k, p) 
    type(space_angle_grid) grid
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: source, rad_scatter
    double precision, dimension(:) :: path_length, path_spacing, a_tilde, gn
    integer i, j, k, p
    integer num_cells
    double precision integral

    call traverse_ray(grid, iops, source, i, j, k, p, path_length, path_spacing, a_tilde, gn, num_cells)
    rad_scatter(i,j,k,p) = calculate_ray_integral(num_cells, path_length, path_spacing, a_tilde, gn)
  end subroutine integrate_ray

  function calculate_ray_integral(num_cells, s, ds, a_tilde, gn) result(integral)
    integer num_cells
    double precision, dimension(num_cells) :: s, ds, a_tilde, gn
    double precision integral, bi, di
    integer i, j

    integral = 0
    do i=1, num_cells
       bi = -a_tilde(i)*s(i+1)
       do j=1, num_cells-1
          bi = bi - a_tilde(j)*ds(j)
       end do

       if(a_tilde(i) .eq. 0) then
          di = ds(i)*exp(bi)
       else
          di = (exp(a_tilde(i)*s(i+1))-exp(a_tilde(i)*s(i)))/a_tilde(i)
       end if

       integral = integral + gn(i)*di
    end do

  end function calculate_ray_integral


  subroutine old_integrate_ray(x, y, k, source,&
       grid, iops, indices, rad_scatter, path_spacing, inner_path_integrand,&
       outer_path_integrand, kp_start, kp_stop, direction)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: rad_scatter
    integer i, j, k, p, kp
    type(index_list) indices
    double precision, dimension(:,:,:,:) :: source

    double precision, dimension(:) :: path_spacing
    double precision, dimension(:) :: inner_path_integrand, outer_path_integrand
    double precision path_source

    double precision xmin, xmax, ymin, ymax
    integer nx, ny, nz, nomega
    double precision x, y, z, theta, phi

    integer kp_start, kp_stop, direction
    integer path_length

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    nomega = grid%angles%nomega

    xmin = grid%x%minval
    xmax = grid%x%maxval
    ymin = grid%y%minval
    ymax = grid%y%maxval

    i = indices%i
    j = indices%j
    k = indices%k
    p = indices%p

    path_length = abs(kp_stop - kp_start) + 1

    z = grid%z%vals(k)

    ! Downwelling: kp_start = 1, kp_stop = k, direction = 1
    ! Upwelling: kp_start = nz, kp_stop = k, direction = -1
    do kp = kp_start, kp_stop, direction
       path_source = interpolate_ray_at_depth(&
            x, y, z, kp, nx, ny, nz, xmin, xmax,&
            ymin, ymax, grid%x%vals, grid%y%vals,&
            grid%z%vals, grid%x_factor(p),&
            grid%y_factor(p), source(:,:,:,p))
       !path_integrand(kp) = path_source * absorb_along_path(&
       !     grid, bc, iops, indices, kp)
       ! path_spacing defined in advect_light
       outer_path_integrand(kp) = absorb_along_path(grid, bc, iops, indices,&
            kp, path_source, path_spacing, inner_path_integrand)
    end do
    rad_scatter(i,j,k,p) = trap_rule_uneven(path_spacing, inner_path_integrand, path_length)
  end subroutine old_integrate_ray


  ! Calculate maximum number of cells a path through the grid could take
  ! This is a loose upper bound
  function calculate_max_cells(grid) result(max_cells)
    type(space_angle_grid) grid
    integer max_cells
    double precision dx, dy, zrange, phi_middle

    ! Angle that will have the longest ray
    phi_middle = grid%angles%phi(grid%angles%nphi/2)
    dx = grid%x%spacing(1)
    dy = grid%y%spacing(1)
    zrange = grid%z%maxval - grid%z%minval

    max_cells = grid%z%num + ceiling((1/dx+1/dy)*zrange*tan(phi_middle))
  end function calculate_max_cells

  ! Traverse from surface or bottom to point (xi, yj, zk)
  ! in the direction omega_p, extracting path lengths (ds) and
  ! function values (f) along the way,
  ! as well as number of cells traversed (n).
  subroutine traverse_ray(grid, iops, source, i, j, k, p, s_array, ds, a_tilde, gn, num_cells)
    type(space_angle_grid) grid
    type(optical_properties) iops
    integer i, j, k, p
    double precision, dimension(:,:,:,:) :: source
    double precision, dimension(:), intent(out) :: s_array, ds, a_tilde, gn
    integer t
    integer, intent(out) :: num_cells

    double precision p0x, p0y, p0z
    double precision p1x, p1y, p1z
    double precision z0
    double precision s_tilde
    integer dir_x, dir_y, dir_z
    integer shift_x, shift_y, shift_z
    integer cell_x, cell_y, cell_z
    integer edge_x, edge_y, edge_z
    double precision s_next_x, s_next_y, s_next_z, s_next
    double precision s, integral

    ! Destination point
    p1x = grid%x%vals(i)
    p1y = grid%y%vals(j)
    p1z = grid%z%vals(k)

    ! Direction
    if(p < (grid%angles%nomega-2)/2) then
       ! Upwelling light originates from bottom
       z0 = grid%z%maxval
    else
       ! Downwelling light originates from surface
       z0 = 0.d0
    end if

    ! Total path length from origin to destination
    ! (sign is correct for upwelling and downwelling)
    s_tilde = (p1z - z0)/grid%angles%cos_phi_p(p)

    ! Origin point
    p0x = p1x - s_tilde * grid%angles%sin_phi_p(p) * grid%angles%cos_theta_p(p)
    p0y = p1y - s_tilde * grid%angles%sin_phi_p(p) * grid%angles%sin_theta_p(p)
    p0z = p1z - s_tilde * grid%angles%cos_phi_p(p)

    ! Direction of ray in each dimension. 1 => increasing. -1 => decreasing.
    dir_x = sgn(p1x-p0x)
    dir_y = sgn(p1y-p0y)
    dir_z = sgn(p1z-p0z)

    ! Shifts
    ! Conversion from cell_inds to edge_inds
    ! merge is fortran's ternary operator
    shift_x = merge(1,0,dir_x>0)
    shift_y = merge(1,0,dir_y>0)
    shift_z = merge(1,0,dir_z>0)

    ! Indices for cell containing origin point
    ! TODO: Is this a good idea? Assuming that x spacing is constant
    cell_x = ceiling((p0x-grid%x%minval)/grid%x%spacing(1))
    cell_y = ceiling((p0y-grid%y%minval)/grid%y%spacing(1))
    cell_z = merge(1, grid%z%num, dir_z>0)

    ! Edge indices preceeding starting cells
    edge_x = cell_x + shift_x
    edge_y = cell_y + shift_y
    edge_z = cell_z + shift_z

    ! Path length to next edge plane in each dimension
    s_next_x = (grid%x%edges(edge_x) - p0x)/(p1x-p0x)
    s_next_y = (grid%y%edges(edge_y) - p0y)/(p1y-p0y)
    s_next_z = (grid%z%edges(edge_z) - p0z)/(p1z-p0z)

    ! Initialize loop variables
    s = 0.d0
    integral = 0.d0

    ! Start with t=0 so that we can increment before storing,
    ! so that t will be the number of grid cells at the end of the loop.
    t=0

    do while (s .lt. s_tilde)
       ! Move cell counter
       t = t + 1

       ! Store path length
       s_array(t) = s

       ! Extract function values
       a_tilde(t) = iops%abs_grid(cell_x, cell_y, cell_z)
       gn(t) = source(cell_x, cell_y, cell_z, p)

       ! Move to next cell in path
       if (s_next_x .le. min(s_next_y, s_next_z)) then
          ! x edge is closest
          s_next = s_next_x
          cell_x = cell_x + dir_x
          edge_x = edge_x + dir_x
          s_next_x = (grid%x%edges(edge_x) - p0x)/(p1x-p0x)
       else if (s_next_y .le. min(s_next_x, s_next_z)) then
          ! y edge is closest
          s_next = s_next_y
          cell_y = cell_y + dir_y
          edge_y = edge_y + dir_y
          s_next_y = (grid%y%edges(edge_y) - p0y)/(p1y-p0y)
       else
          ! z edge is closest
          s_next = s_next_z
          cell_z = cell_z + dir_z
          edge_z = edge_z + dir_z
          s_next_z = (grid%z%edges(edge_z) - p0z)/(p1z-p0z)
       end if

       ! Extract path length from same cell as function vals
       ds(t) = s_next - s

       ! Update path length
       s = s_next
    end do

    ! Return number of cells traversed
    num_cells = t
  end subroutine traverse_ray

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
    theta = grid%angles%theta_p(indices%p)
    phi = grid%angles%phi_p(indices%p)

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
              grid%x_factor(indices%p),&
              grid%y_factor(indices%p),&
              iops%abs_grid)
       end do
    end if

    ! Only pass relevant portion of path_spacing to integrator
    path_spacing(:path_length) = grid%z%spacing(kmin:indices%k:direction) / abs(grid%angles%cos_phi_p(indices%p))
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
