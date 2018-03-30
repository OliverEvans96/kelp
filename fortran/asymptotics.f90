module asymptotics
  use rte3d
  implicit none
  contains

  subroutine calculate_light_with_scattering(grid, bc, iops, radiance, num_scatters)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: radiance
    double precision, dimension(:,:,:,:), allocatable :: rad_scatter, source
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
    allocate(source(nx, ny, nz, nomega))


    ! Scattering coefficient
    ! Allowed to vary over depth
    ! Reset radiance
    !write(*,*) 'Before Scattering'
    ! write(*,*) '-1'
    call calculate_light_before_scattering(grid, bc, iops, source, radiance, path_length, path_spacing, a_tilde, gn)
    ! write(*,*) '0'

    if (num_scatters .gt. 0) then
       call calculate_light_after_scattering(grid, bc, iops, source, radiance, num_scatters, path_length, path_spacing, a_tilde, gn)
    end if

    ! write(*,*) '1'
    deallocate(path_length)
    ! write(*,*) '2'
    deallocate(path_spacing)
    ! write(*,*) '3'
    deallocate(a_tilde)
    ! write(*,*) '4'
    deallocate(gn)
    ! write(*,*) '5'
    deallocate(source)
    ! write(*,*) '6'
  end subroutine calculate_light_with_scattering

  subroutine calculate_light_before_scattering(grid, bc, iops, source, radiance, path_length, path_spacing, a_tilde, gn)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: radiance, source 
    double precision, dimension(:) :: path_length, path_spacing, a_tilde, gn
    integer i, j, k, p
    double precision surface_val
    integer num_cells

    ! Downwelling light
    do p=1, grid%angles%nomega/2
       surface_val = bc%bc_grid(p)
       do i=1, grid%x%num
          do j=1, grid%y%num
            do k=1, grid%z%num
               call attenuate_light_from_surface(&
                    grid, iops, source, i, j, k, p,&
                    radiance, path_length, path_spacing,&
                    a_tilde, gn, bc)
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

  subroutine attenuate_light_from_surface(grid, iops, source, i, j, k, p,&
       radiance, path_length, path_spacing, a_tilde, gn, bc)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: radiance, source
    double precision, dimension(:) :: path_length, path_spacing, a_tilde, gn
    integer i, j, k, p
    integer num_cells
    double precision atten
    logical print_flag
    integer ip

    if(i .eq. 1 .and. j .eq. 1 .and. k .eq. 1 .and. p .eq. 7) then
       print_flag = .false.
    else
       print_flag = .false.
    end if


    ! Don't need gn here, so just ignore it
    call traverse_ray(grid, iops, source, i, j, k, p, path_length, path_spacing, a_tilde, gn, num_cells, print_flag)
    ! Start with surface bc and attenuate along path

    ! write(*,*)
    ! write(*,*) 'ds =', path_spacing(1:num_cells)
    ! write(*,*) 'a_tilde = ', a_tilde(1:num_cells)

    ! do ip=1, num_cells
    !    write(*,*) 'ip, prod =', ip, path_spacing(ip) * a_tilde(ip)
    ! end do

    ! write(*,*) 'atten_arg = ', dot_product(path_spacing(1:num_cells), a_tilde(1:num_cells))
    atten = sum(path_spacing(1:num_cells) * a_tilde(1:num_cells))
    radiance(i,j,k,p) = bc%bc_grid(p) * exp(-atten)
    !if(atten > 50.d0) then
    !   ! Prevent underflow
    !   radiance(i,j,k,p) = 0.d0
    !else
    !end if

    ! if(i .eq. 1 .and. j .eq. 1 .and. k .eq. 1 .and. p .eq. 7) then
    !   write(*,*) 'i, j, k, p =', i, j, k, p
    !   write(*,*) 'phi = ', grid%angles%phi_p(p) * 180/pi
    !   write(*,*) 'num_cells = ', num_cells
    !   write(*,*) 'path_spacing =', path_spacing(1:num_cells)
    !   write(*,*) 'bc, exp_arg, rad = ', bc%bc_grid(p), -sum(path_spacing(1:num_cells) * a_tilde(1:num_cells)), radiance(i,j,k,p)
    !   write(*,*)
    ! end if

  end subroutine attenuate_light_from_surface

  subroutine calculate_light_after_scattering(grid, bc, iops, source, radiance,&
       num_scatters, path_length, path_spacing, a_tilde, gn)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: radiance, source
    integer num_scatters
    double precision, dimension(:) :: path_length, path_spacing, a_tilde, gn
    double precision, dimension(:,:,:,:), allocatable :: rad_scatter
    integer n, k
    double precision bb

    allocate(rad_scatter(grid%x%num, grid%y%num, grid%z%num, grid%angles%nomega))
    rad_scatter = radiance

    ! For now, just use first entry from scattering array everywhere
    bb = iops%scat_water(1)
    write(*,*) 'prescatter = ', rad_scatter(2,2,:,2)

    do n=1, num_scatters
       write(*,*) 'scatter #', n
       call scatter(grid, bc, iops, source, rad_scatter, path_length, path_spacing, a_tilde, gn)
       radiance = radiance + bb**n * rad_scatter
       write(*,*) 'rad_scatter = ', rad_scatter(2,2,:,2)
    end do

    write(*,*) 'deallocate rad_scatter'
    deallocate(rad_scatter)
    write(*,*) 'deallocated'
  end subroutine calculate_light_after_scattering

  ! Perform one scattering event
  subroutine scatter(grid, bc, iops, source, rad_scatter, path_length, path_spacing, a_tilde, gn)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: rad_scatter, source
    double precision, dimension(:,:,:,:), allocatable :: scatter_integral
    double precision, dimension(:), allocatable :: scatter_integrand
    double precision, dimension(:) :: path_length, path_spacing, a_tilde, gn
    double precision path_source
    integer max_cells
    integer nx, ny, nz, nomega

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    nomega = grid%angles%nomega

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

    write(*,*) 'A'
    deallocate(scatter_integral)
    write(*,*) 'B'
    deallocate(scatter_integrand)
    write(*,*) 'C'
  end subroutine scatter

  ! Calculate source from no-scatter or previous scattering layer
  subroutine calculate_source(grid, iops, rad_scatter, source, scatter_integral, scatter_integrand)
    type(space_angle_grid) grid
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: rad_scatter, source, scatter_integral
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
             do p=1, nomega
                indices%p = p
                call calculate_scatter_integral(&
                    grid, iops, rad_scatter,&
                    scatter_integral,&
                    scatter_integrand, indices)
                ! write(*,*) 'integrand ', i, j, k, p
                ! write(*,*) scatter_integrand(:)
                ! write(*,*) 'integral', scatter_integral(i,j,k,p)
                ! write(*,*) ''
             end do
          end do
       end do
    end do

    source = -rad_scatter + scatter_integral

  end subroutine calculate_source

  subroutine calculate_scatter_integral(grid, iops, rad_scatter, scatter_integral, scatter_integrand, indices)
    type(space_angle_grid) grid
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: rad_scatter, scatter_integral
    double precision, dimension(:) :: scatter_integrand
    type(index_list) indices
    integer pp

    ! Current direction is already excluded by VSF
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

   if(indices%i .eq. 2 .and. indices%j .eq. 5 .and. indices%k .eq. 3 .and. indices%p .eq. 5) then
       write(*,*) ''
       write(*,*) 'theta'
       write(*,*) grid%angles%theta
       write(*,*) 'phi'
       write(*,*) grid%angles%phi
       write(*,*) 'integrand ', indices%i, indices%j, indices%k, indices%p
       write(*,*) scatter_integrand(:)
       write(*,*) 'integral', scatter_integral(indices%i,indices%j,indices%k,indices%p)
       write(*,*) ''
       write(*,*) ''
       write(*,*) ''
    end if

  end subroutine calculate_scatter_integral

  subroutine advect_light(grid, iops, source, rad_scatter, path_length, path_spacing, a_tilde, gn)
    type(space_angle_grid) grid
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: rad_scatter, source
    double precision, dimension(:) :: path_length, path_spacing, a_tilde, gn
    double precision path_source

    integer i, j, k, p, kp
    type(index_list) indices


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

    call traverse_ray(grid, iops, source, i, j, k, p, path_length, path_spacing, a_tilde, gn, num_cells, .false.)
    rad_scatter(i,j,k,p) = calculate_ray_integral(num_cells, path_length, path_spacing, a_tilde, gn)
  end subroutine integrate_ray

  function calculate_ray_integral(num_cells, s, ds, a_tilde, gn) result(integral)
    ! Double integral which accumulates all scattered light along the path
    ! via an angular integral and attenuates it by integrating along the path
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
  subroutine traverse_ray(grid, iops, source, i, j, k, p, s_array, ds, a_tilde, gn, num_cells, print_flag)
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
    double precision s_tilde, s
    integer dir_x, dir_y, dir_z
    integer shift_x, shift_y, shift_z
    integer cell_x, cell_y, cell_z
    integer edge_x, edge_y, edge_z
    integer first_x, last_x, first_y, last_y, last_z
    double precision s_next_x, s_next_y, s_next_z, s_next
    double precision x_factor, y_factor, z_factor
    double precision ds_x, ds_y
    double precision, dimension(grid%z%num) :: ds_z
    logical print_flag
    double precision test, zero
    double precision smx, smy
    zero = 0.d0

    ! Divide by these numbers to get path separation
    ! from separation in individual dimensions
    x_factor = grid%angles%sin_phi_p(p) * grid%angles%cos_theta_p(p)
    y_factor = grid%angles%sin_phi_p(p) * grid%angles%sin_theta_p(p)
    z_factor = grid%angles%cos_phi_p(p)

    if(print_flag) then
       write(*,*) 'x_factor = ', x_factor
       write(*,*) 'y_factor = ', y_factor
       write(*,*) 'z_factor = ', z_factor
       write(*,*) 'sin(phi) = ', grid%angles%sin_phi_p(p)
       write(*,*) 'sin(theta) = ', grid%angles%sin_theta_p(p)
    end if

    ! This one is an array because z spacing can vary
    ! z_factor should never be 0, because the ray will never
    ! reach the surface or bottom.
    ds_z(1:grid%z%num) = grid%z%spacing(1:grid%z%num)/z_factor

    ! Destination point
    p1x = grid%x%vals(i)
    p1y = grid%y%vals(j)
    p1z = grid%z%vals(k)

    !write(*,*) 'T1'
    ! Direction
    if(p .le. grid%angles%nomega/2) then
       ! Downwelling light originates from surface
       z0 = grid%z%minval
       dir_z = 1
    else
       ! Upwelling light originates from bottom
       z0 = grid%z%maxval
       dir_z = -1
    end if

    !write(*,*) 'T2'
    ! Total path length from origin to destination
    ! (sign is correct for upwelling and downwelling)
    s_tilde = (p1z - z0)/grid%angles%cos_phi_p(p)

    ! Path spacings between edge intersections in each dimension
    ! Set to 2*s_tilde if infinite in this dimension so that it's unreachable
    ! Assume x & y spacings are even, so it's okay to just use the first value.
    if(x_factor .eq. 0) then
       ds_x = 2*s_tilde
    else
       ds_x = abs(grid%x%spacing(1)/x_factor)
    end if
    if(y_factor .eq. 0) then
       ds_y = 2*s_tilde
    else
       ds_y = abs(grid%y%spacing(1)/y_factor)
    end if
    !write(*,*) 'T3'

    ! Origin point
    p0x = p1x - s_tilde * x_factor
    p0y = p1y - s_tilde * y_factor
    p0z = p1z - s_tilde * z_factor

    ! Direction of ray in each dimension. 1 => increasing. -1 => decreasing.
    dir_x = sgn(p1x-p0x)
    dir_y = sgn(p1y-p0y)

    ! Shifts
    ! Conversion from cell_inds to edge_inds
    ! merge is fortran's ternary operator
    shift_x = merge(1,0,dir_x>0)
    shift_y = merge(1,0,dir_y>0)
    shift_z = merge(1,0,dir_z>0)

    ! Indices for cell containing origin point
    cell_x = floor((p0x-grid%x%minval)/grid%x%spacing(1)) + 1
    cell_y = floor((p0y-grid%y%minval)/grid%y%spacing(1)) + 1
    ! x and y may be in periodic image, so shift back.
    cell_x = mod1(cell_x, grid%x%num)
    cell_y = mod1(cell_y, grid%y%num)

    !write(*,*) 'T4'
    ! z starts at top or bottom depending on direction.
    if(dir_z > 0) then
       cell_z = 1
    else
       cell_z = grid%z%num
    end if

    ! Edge indices preceeding starting cells
    edge_x = mod1(cell_x + shift_x, grid%x%num)
    edge_y = mod1(cell_y + shift_y, grid%y%num)
    edge_z = mod1(cell_z + shift_z, grid%z%num)

    ! First and last cells given direction
    if(dir_x .gt. 0) then
       first_x = 1
       last_x = grid%x%num
    else
       first_x = grid%x%num
       last_x = 1
    end if
    if(dir_y .gt. 0) then
       first_y = 1
       last_y = grid%y%num
    else
       first_y = grid%y%num
       last_y = 1
    end if
    if(dir_z .gt. 0) then
       last_z = grid%z%num
    else
       last_z = 1
    end if

    ! Calculate periodic images
    smx = shift_mod(p0x, grid%x%minval, grid%x%maxval)
    smy = shift_mod(p0y, grid%y%minval, grid%y%maxval)

    !write(*,*) 'T5'
    ! Path length to next edge plane in each dimension
    if(x_factor .eq. 0) then
       ! Will never cross, so set above total path length
       s_next_x = 2*s_tilde
    else if(cell_x .eq. last_x) then
       ! If starts out at last cell, then compare to periodic image
       s_next_x = (grid%x%edges(first_x) + dir_x * (grid%x%maxval - grid%x%minval)&
            - smx) / x_factor
    else
       ! Otherwise, just compare to next cell
       s_next_x = (grid%x%edges(edge_x) - smx) / x_factor
    end if

    ! Path length to next edge plane in each dimension
    if(y_factor .eq. 0) then
       ! Will never cross, so set above total path length
       s_next_y = 2*s_tilde
    else if(cell_y .eq. last_y) then
       ! If starts out at last cell, then compare to periodic image
       s_next_y = (grid%y%edges(first_y) + dir_y * (grid%y%maxval - grid%y%minval)&
            - smy) / y_factor
    else
       ! Otherwise, just compare to next cell
       s_next_y = (grid%y%edges(edge_y) - smy) / y_factor
    end if

    s_next_z = ds_z(cell_z)

    ! This should not happen
    if(s_next_x .lt. 0 .or. s_next_y .lt. 0 .or. s_next_z .lt. 0 .or. s_tilde .lt. 0) then
       write(*,*) 'NEGATIVE S NEXT'
    end if


    !write(*,*) 'T6'
    ! Initialize path
    s = 0.d0

    ! Start with t=0 so that we can increment before storing,
    ! so that t will be the number of grid cells at the end of the loop.
    t=0

    if(print_flag) then
       write(*,*) 'i, j, k, p =', i, j, k, p
       write(*,*) 'l, m =', grid%angles%lhat(p), grid%angles%mhat(p)
       write(*,*) 'theta, phi = ', 180/pi*grid%angles%theta_p(p), 180/pi*grid%angles%phi_p(p)

       write(*,*) 'z0 =', z0

      write(*,*) 'p0 = ', p0x, p0y, p0z
      write(*,*) 'p1 = ', p1x, p1y, p1z

      write(*,*) 'dir_x =', dir_x
      write(*,*) 'dir_y =', dir_y
      write(*,*) 'dir_z =', dir_z

      write(*,*) 'shift_x = ', shift_x
      write(*,*) 'shift_y = ', shift_y
      write(*,*) 'shift_z = ', shift_z

      write(*,*) 'xmin = ', grid%x%minval
      write(*,*) 'ymin = ', grid%y%minval
      write(*,*) 'zmin = ', grid%z%minval

      write(*,*) 'dx = ', grid%x%spacing(1)
      write(*,*) 'dy = ', grid%y%spacing(1)
      write(*,*) 'dz = ', grid%z%spacing(1)

      write(*,*) 'init cell_x = ', cell_x
      write(*,*) 'init cell_y = ', cell_y
      write(*,*) 'init cell_z = ', cell_z

      write(*,*) 's_next_x =', s_next_x
      write(*,*) 's_next_y =', s_next_y
      write(*,*) 's_next_z =', s_next_z

      write(*,*) 's, s_tilde = ', s, s_tilde

      write(*,*) ''
      write(*,*) '---'
   end if

   !write(*,*) 'T7'
    ! s is the beginning of the current cell,
    ! s_next is the end of the current cell.
    do while (s .lt. s_tilde)
       ! Move cell counter
       t = t + 1

       if(print_flag .and. t .lt. 100) then
          write(*,*)
          write(*,*) 't = ', t
          !write(*,*) 'cell_x = ', cell_x
          !write(*,*) 'cell_y = ', cell_y
          !write(*,*) 'cell_z = ', cell_z
          write(*,*) 'traverse cell ', cell_x, cell_y, cell_z

          write(*,*) 'x plane: ', grid%x%edges(edge_x)
          write(*,*) 'y plane: ', grid%y%edges(edge_y)
          write(*,*) 'z plane: ', grid%z%edges(edge_z)

          write(*,*) 'x center: ', grid%x%vals(cell_x)
          write(*,*) 'y center: ', grid%y%vals(cell_y)
          write(*,*) 'z center: ', grid%z%vals(cell_z)

          write(*,*) 's_next_x =', s_next_x
          write(*,*) 's_next_y =', s_next_y
          write(*,*) 's_next_z =', s_next_z
       end if

       ! Extract function values
       if(t > size(a_tilde)) then
          write(*,*) 'TOO BIG: t, max_cells = ', t, size(a_tilde)
          !write(*,*) 'i, j, k, p =', i, j, k, p
       end if

       a_tilde(t) = iops%abs_grid(cell_x, cell_y, cell_z)
       gn(t) = source(cell_x, cell_y, cell_z, p)

       !write(*,*) 'a_tilde = ', a_tilde(t)

       ! Move to next cell in path
       if(s_next_x .le. min(s_next_y, s_next_z)) then
          ! x edge is closest
          s_next = s_next_x

          ! Increment indices (periodic)
          cell_x = mod1(cell_x + dir_x, grid%x%num)
          edge_x = mod1(edge_x + dir_x, grid%x%num)

          ! x intersection after the one at s=s_next
          s_next_x = s_next + ds_x

       else if (s_next_y .le. min(s_next_x, s_next_z)) then
          ! y edge is closest
          s_next = s_next_y

          ! Increment indices (periodic)
          cell_y = mod1(cell_y + dir_y, grid%y%num)
          edge_y = mod1(edge_y + dir_y, grid%y%num)

          ! y intersection after the one at s=s_next
          s_next_y = s_next + ds_y

       else
          ! z edge is closest
          s_next = s_next_z

          ! Increment indices
          cell_z = cell_z + dir_z
          edge_z = edge_z + dir_z

          ! z intersection after the one at s=s_next
          if(cell_z .ne. last_z) then
             ! Only look ahead if we aren't at the end
             s_next_z = s_next + ds_z(cell_z)
          else
             ! Otherwise, no need to continue.
             ! this is our final destination.
             s_next_z = 2*s_tilde
          end if

       end if

       ! Cut off early if this is the end
       ! This will be the last cell traversed if s_next >= s_tilde
       s_next = min(s_tilde, s_next)


       ! Store path length
       s_array(t) = s_next
       ! Extract path length from same cell as function vals
       ds(t) = s_next - s

       if(print_flag .and. t .lt. 100) then
        write(*,*) 's, s_next = ', s, s_next
        write(*,*) 's_next/s_tilde = ', s_next / s_tilde
        write(*,*) 'ds = ', ds(t)
       end if

       ! Update path length
       s = s_next
    end do
    !write(*,*) 'T8'

    ! Return number of cells traversed
    num_cells = t

    !write(*,*) 'final position', p1x, p1y, p1z
    !write(*,*) 'T9'

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
