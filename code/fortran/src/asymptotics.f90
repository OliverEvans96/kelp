module asymptotics
  use kelp_context
  !use rte_sparse_matrices
  !use light_context
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

    call calculate_light_before_scattering(grid, bc, iops, source, radiance, path_length, path_spacing, a_tilde, gn)

    if (num_scatters .gt. 0) then
       call calculate_light_after_scattering(&
            grid, bc, iops, source, radiance, &
            num_scatters, path_length, path_spacing, &
            a_tilde, gn)
    end if

    deallocate(path_length)
    deallocate(path_spacing)
    deallocate(a_tilde)
    deallocate(gn)
    deallocate(source)
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
    integer ip

    ! Don't need gn here, so just ignore it
    call traverse_ray(grid, iops, source, i, j, k, p, path_length, path_spacing, a_tilde, gn, num_cells)

    ! Start with surface bc and attenuate along path
    atten = sum(path_spacing(1:num_cells) * a_tilde(1:num_cells))
    radiance(i,j,k,p) = bc%bc_grid(p) * exp(-atten)
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
    bb = iops%scat

    do n=1, num_scatters
       write(*,*) 'scatter #', n
       call scatter(grid, bc, iops, source, rad_scatter, path_length, path_spacing, a_tilde, gn)
       radiance = radiance + bb**n * rad_scatter
    end do

    deallocate(rad_scatter)
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

    call calculate_source(grid, iops, rad_scatter, source, scatter_integral, scatter_integrand)
    call advect_light(grid, iops, source, rad_scatter, path_length, path_spacing, a_tilde, gn)

    deallocate(scatter_integral)
    deallocate(scatter_integrand)
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
             end do
          end do
       end do
    end do

    source(:,:,:,:) = -rad_scatter(:,:,:,:) + scatter_integral(:,:,:,:)

  end subroutine calculate_source

  subroutine calculate_scatter_integral(grid, iops, rad_scatter, scatter_integral, scatter_integrand, indices)
    type(space_angle_grid) grid
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: rad_scatter, scatter_integral
    double precision, dimension(:) :: scatter_integrand
    type(index_list) indices
    integer pp
    double precision s

    scatter_integral(indices%i,indices%j,indices%k,indices%p) &
         = sum(iops%vsf_integral(indices%p, :) &
         * rad_scatter(&
           indices%i,&
           indices%j,&
           indices%k,:))

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

    do p=1, nomega
       do k=1, nz
          do i=1, nx
            do j=1, ny
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
    ! Double integral which accumulates all scattered light along the path
    ! via an angular integral and attenuates it by integrating along the path
    integer num_cells
    double precision, dimension(num_cells) :: ds, a_tilde, gn
    double precision, dimension(num_cells+1) :: s
    double precision integral, bi, di
    integer i, j

    integral = 0
    do i=1, num_cells
       bi = -a_tilde(i)*s(i+1)
       do j=i+1, num_cells
          bi = bi - a_tilde(j)*ds(j)
       end do

       if(a_tilde(i) .eq. 0) then
          di = ds(i)
       else
          di = (exp(a_tilde(i)*s(i+1))-exp(a_tilde(i)*s(i)))/a_tilde(i)
       end if

       integral = integral + gn(i)*di * exp(bi)
    end do

  end function calculate_ray_integral

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
    double precision smx, smy

    ! Divide by these numbers to get path separation
    ! from separation in individual dimensions
    x_factor = grid%angles%sin_phi_p(p) * grid%angles%cos_theta_p(p)
    y_factor = grid%angles%sin_phi_p(p) * grid%angles%sin_theta_p(p)
    z_factor = grid%angles%cos_phi_p(p)

    ! Destination point
    p1x = grid%x%vals(i)
    p1y = grid%y%vals(j)
    p1z = grid%z%vals(k)

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

    ! Total path length from origin to destination
    ! (sign is correct for upwelling and downwelling)
    s_tilde = (p1z - z0)/grid%angles%cos_phi_p(p)

    ! Path spacings between edge intersections in each dimension
    ! Set to 2*s_tilde if infinite in this dimension so that it's unreachable
    ! Assume x & y spacings are uniform,
    ! so it's okay to just use the first value.
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

    ! This one is an array because z spacing can vary
    ! z_factor should never be 0, because the ray will never
    ! reach the surface or bottom.
    ds_z(1:grid%z%num) = dir_z * grid%z%spacing(1:grid%z%num)/z_factor

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

    ! Path length to next edge plane in each dimension
    if(abs(x_factor) .lt. 1.d-10) then
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
    if(abs(y_factor) .lt. 1.d-10) then
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

    ! Initialize path
    s = 0.d0
    s_array(1) = 0.d0

    ! Start with t=0 so that we can increment before storing,
    ! so that t will be the number of grid cells at the end of the loop.
    t=0

    ! s is the beginning of the current cell,
    ! s_next is the end of the current cell.
    do while (s .lt. s_tilde)
       ! Move cell counter
       t = t + 1

       ! Extract function values
       a_tilde(t) = iops%abs_grid(cell_x, cell_y, cell_z)
       gn(t) = source(cell_x, cell_y, cell_z, p)

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
       s_array(t+1) = s_next
       ! Extract path length from same cell as function vals
       ds(t) = s_next - s

       ! Update path length
       s = s_next
    end do

    ! Return number of cells traversed
    num_cells = t

  end subroutine traverse_ray
end module asymptotics
