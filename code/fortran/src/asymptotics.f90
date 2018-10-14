module asymptotics
  use kelp_context
  !use rte_sparse_matrices
  !use light_context
  implicit none
  contains

  subroutine calculate_asymptotic_light_field(grid, bc, iops, source, radiance, num_scatters, num_threads)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: radiance
    double precision, dimension(:,:,:,:), allocatable :: rad_scatter
    double precision, dimension(:,:,:,:) :: source
    integer num_scatters
    integer nx, ny, nz, nomega
    integer max_cells
    integer n
    logical bc_flag
    integer num_threads

    double precision bb

    double precision, dimension(:), allocatable :: path_length, path_spacing, a_tilde, gn

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    nomega = grid%angles%nomega

    max_cells = calculate_max_cells(grid)

    allocate(path_length(max_cells+1))
    allocate(path_spacing(max_cells))
    allocate(a_tilde(max_cells))
    allocate(gn(max_cells))
    allocate(rad_scatter(grid%x%num, grid%y%num, grid%z%num, grid%angles%nomega))

    write(*,*) 'before'
    write(*,*) 'min radiance =', minval(radiance)
    write(*,*) 'max radiance =', maxval(radiance)
    write(*,*) 'mean radiance =', sum(radiance)/size(radiance)

    write(*,*) 'advect source + bc'
    bc_flag = .true.
    call advect_light( &
         grid, iops, source, radiance, &
         path_length, path_spacing, &
         a_tilde, gn, bc_flag, num_threads, bc)

    write(*,*) 'after'
    write(*,*) 'min radiance =', minval(radiance)
    write(*,*) 'max radiance =', maxval(radiance)
    write(*,*) 'mean radiance =', sum(radiance)/size(radiance)

    rad_scatter = radiance
    bb = iops%scat

    do n=1, num_scatters
       write(*,*) 'scatter #', n
       call scatter(grid, iops, source, rad_scatter, path_length, path_spacing, a_tilde, gn, num_threads)
       radiance = radiance + bb**n * rad_scatter
    end do

    write(*,*) 'asymptotics complete'

    deallocate(path_length)
    deallocate(path_spacing)
    deallocate(a_tilde)
    deallocate(gn)
    deallocate(rad_scatter)
  end subroutine calculate_asymptotic_light_field

  subroutine calculate_asymptotic_light_field_expanded_source(&
       grid, bc, iops, source, &
       source_expansion, radiance, &
       num_scatters, num_threads)
    type(space_angle_grid) grid
    type(boundary_condition) bc
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: radiance
    double precision, dimension(:,:,:,:,:) :: source_expansion
    integer num_scatters
    integer nx, ny, nz, nomega
    integer max_cells
    integer n
    logical bc_flag
    integer num_threads

    double precision bb

    double precision, dimension(:), allocatable :: path_length, path_spacing, a_tilde, gn
    double precision, dimension(:,:,:,:), allocatable :: source
    double precision, dimension(:,:,:,:), allocatable :: rad_scatter
    double precision, dimension(:,:,:,:), allocatable :: scatter_integral

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    nomega = grid%angles%nomega

    max_cells = calculate_max_cells(grid)

    allocate(path_length(max_cells+1))
    allocate(path_spacing(max_cells))
    allocate(a_tilde(max_cells))
    allocate(gn(max_cells))
    allocate(rad_scatter(grid%x%num, grid%y%num, grid%z%num, grid%angles%nomega))
    allocate(scatter_integral(nx, ny, nz, nomega))

    write(*,*) 'advect source + bc'
    bc_flag = .true.
    call advect_light( &
         grid, iops, source_expansion(:,:,:,:,1), radiance, &
         path_length, path_spacing, &
         a_tilde, gn, bc_flag, num_threads, bc)
    ! Disable BC for scattering advection
    bc_flag = .false.

    rad_scatter = radiance
    bb = iops%scat

    do n=1, num_scatters
       write(*,*) 'scatter #', n
       call calculate_scatter_source(grid, iops, rad_scatter, source, scatter_integral, num_threads)
       source = source + source_expansion(:,:,:,:,n+1)
       call advect_light(grid, iops, source, rad_scatter, path_length, path_spacing, a_tilde, gn, bc_flag, num_threads)

       radiance = radiance + bb**n * rad_scatter

    end do

    write(*,*) 'asymptotics complete'

    deallocate(path_length)
    deallocate(path_spacing)
    deallocate(a_tilde)
    deallocate(gn)
    deallocate(rad_scatter)
    deallocate(scatter_integral)
  end subroutine calculate_asymptotic_light_field_expanded_source

  ! Add attenuated surface light to existing radiance
  subroutine advect_surface_bc(&
       i, j, k, p, radiance, &
       path_spacing, num_cells, a_tilde, bc)
    type(boundary_condition) bc
    double precision, dimension(:,:,:,:) :: radiance
    double precision, dimension(:) :: path_spacing, a_tilde
    integer i, j, k, p
    integer num_cells
    double precision atten

    atten = sum(path_spacing(1:num_cells) * a_tilde(1:num_cells))
    ! Avoid underflow
    if(atten .lt. 100.d0) then
       radiance(i,j,k,p) = radiance(i,j,k,p) + bc%bc_grid(p) * exp(-atten)
    end if
  end subroutine advect_surface_bc

  ! Perform one scattering event
  subroutine scatter(grid, iops, source, rad_scatter, path_length, path_spacing, a_tilde, gn, num_threads)
    type(space_angle_grid) grid
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: rad_scatter, source
    double precision, dimension(:,:,:,:), allocatable :: scatter_integral
    double precision, dimension(:) :: path_length, path_spacing, a_tilde, gn
    integer nx, ny, nz, nomega
    logical bc_flag
    integer num_threads

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    nomega = grid%angles%nomega
    bc_flag = .false.

    allocate(scatter_integral(nx, ny, nz, nomega))

    call calculate_scatter_source(grid, iops, rad_scatter, source, scatter_integral, num_threads)
    call advect_light(grid, iops, source, rad_scatter, path_length, path_spacing, a_tilde, gn, bc_flag, num_threads)

    deallocate(scatter_integral)
  end subroutine scatter

  ! Calculate source from no-scatter or previous scattering layer
  subroutine calculate_scatter_source(grid, iops, rad_scatter, source, scatter_integral, num_threads)
    type(space_angle_grid) grid
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: rad_scatter
    double precision, dimension(:,:,:,:) :: source
    double precision, dimension(:,:,:,:) :: scatter_integral
    type(index_list) indices
    integer nx, ny, nz, nomega
    integer i, j, k, p
    integer num_threads

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num
    nomega = grid%angles%nomega

    ! Use collapse to combine outer two loops
    ! to avoid manually deciding nesting details.

    !$omp parallel do default(none) private(indices) &
    !$omp private(i,j,k,p) shared(nx,ny,nz,nomega) &
    !$omp shared(iops, rad_scatter, scatter_integral) &
    !$omp num_threads(num_threads) collapse(2)

    do k=1, nz
       do i=1, nx
          indices%k = k
          indices%i = i
          do j=1, ny
             indices%j = j
             do p=1, nomega
                indices%p = p
                call calculate_scatter_integral(&
                    iops, rad_scatter,&
                    scatter_integral,&
                    indices)
             end do
          end do
       end do
    end do
    !$omp end parallel do

    source(:,:,:,:) = -rad_scatter(:,:,:,:) + scatter_integral(:,:,:,:)

    write(*,*) 'source min: ', minval(source)
    write(*,*) 'source max: ', maxval(source)
    write(*,*) 'source mean: ', sum(source)/size(source)

  end subroutine calculate_scatter_source

  subroutine calculate_scatter_integral(iops, rad_scatter, scatter_integral, indices)
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: rad_scatter, scatter_integral
    type(index_list) indices

    scatter_integral(indices%i,indices%j,indices%k,indices%p) &
         = sum(iops%vsf_integral(indices%p, :) &
         * rad_scatter(&
           indices%i,&
           indices%j,&
           indices%k,:))

  end subroutine calculate_scatter_integral

  subroutine advect_light(grid, iops, source, rad_scatter, path_length, path_spacing, a_tilde, gn, bc_flag, num_threads, bc)
    type(space_angle_grid) grid
    type(optical_properties) iops
    double precision, dimension(:,:,:,:) :: rad_scatter, source
    double precision, dimension(:) :: path_length, path_spacing, a_tilde, gn
    logical bc_flag
    type(boundary_condition), intent(in), optional :: bc
    integer i, j, k, p
    integer num_threads

    !$omp parallel do default(none) &
    !$omp private(i,j,k,p) &
    !$omp shared(rad_scatter,source,grid,iops,bc_flag,bc) &
    !$omp private(path_length,path_spacing,a_tilde,gn) &
    !$omp num_threads(num_threads) collapse(2)
    do k=1, grid%z%num
       do i=1, grid%x%num
          do j=1, grid%y%num
             do p=1, grid%angles%nomega
                call integrate_ray(grid, iops, source,&
                     rad_scatter, path_length, path_spacing,&
                     a_tilde, gn, i, j, k, p, bc_flag, bc)
            end do
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine advect_light

  ! New algorithm, double integral over piecewise-constant 1d funcs
  subroutine integrate_ray(grid, iops, source, rad_scatter, path_length, path_spacing, a_tilde, gn, i, j, k, p, bc_flag, bc)
    type(space_angle_grid) :: grid
    type(optical_properties) :: iops
    double precision, dimension(:,:,:,:) :: source
    double precision, dimension(:,:,:,:) :: rad_scatter
    integer :: i, j, k, p
    ! The following are only passed to avoid unnecessary allocation
    double precision, dimension(:) :: path_length, path_spacing, a_tilde, gn
    logical bc_flag
    type(boundary_condition), intent(in), optional :: bc

    integer num_cells

    call traverse_ray(grid, iops, source, i, j, k, p, path_length, path_spacing, a_tilde, gn, num_cells)
    rad_scatter(i,j,k,p) = calculate_ray_integral(num_cells, path_length, path_spacing, a_tilde, gn)

    if(bc_flag .and. p .le. grid%angles%nomega/2) then
       call advect_surface_bc(&
            i, j, k, p, rad_scatter, &
            path_spacing, num_cells, &
            a_tilde, bc)
    end if

    ! if((i .eq. 1) &
    !      .and. (j .eq. 1) &
    !      !.and. (k .eq. grid%z%num/2) &
    !      .and. ( &
    !      (p .eq. 1) .or. (p .eq. grid%angles%nomega) &
    !      )) then

    !    write(*,*) 'ray (', i, ', ', j, ', ', k, ', ', p, ')'
    !    write(*,*) 'num_cells = ', num_cells
    !    write(*,*) 'path_spacing:'
    !    write(*,*) path_spacing(1:num_cells)
    !    write(*,*) 'path_length:'
    !    write(*,*) path_length(1:num_cells+1)
    !    write(*,*) 'a_tilde:'
    !    write(*,*) a_tilde(1:num_cells)
    !    write(*,*) 'gn:'
    !    write(*,*) gn(1:num_cells)
    !    write(*,*)
    ! end if

  end subroutine integrate_ray

  function calculate_ray_integral(num_cells, s, ds, a_tilde, gn) result(integral)
    ! Double integral which accumulates all scattered light along the path
    ! via an angular integral and attenuates it by integrating along the path
    integer :: num_cells
    double precision, dimension(num_cells) :: ds, a_tilde, gn
    double precision, dimension(num_cells+1) :: s
    double precision :: integral
    double precision bi, di_exp_bi
    double precision cutoff
    integer i, j

    ! Maximum absorption coefficient suitable for numerical computation
    cutoff = 10.d0

    integral = 0
    do i=1, num_cells
       bi = -a_tilde(i)*s(i+1)
       do j=i+1, num_cells
          bi = bi - a_tilde(j)*ds(j)
       end do

       ! In this case, so much absorption has occurred
       ! previously on the path that we don't need
       ! to continue, and we might get underflow if we do.
       if(bi .lt. -100.d0) then
          di_exp_bi = 0.d0
       else

         ! Without this conditional, overflow occurs.
         ! Which is unnecessary, because large absorption
         ! means very small light added to the ray
         ! at this grid cell.
         if(a_tilde(i) .lt. cutoff) then
             if(a_tilde(i) .eq. 0) then
                 di_exp_bi = ds(i) * exp(bi)
             else
                ! In an attempt to avoid overflow
                ! and reduce compute time,
                ! I'm combining exponentials.
                ! di*exp(bi) -> di_exp_bi
                di_exp_bi = (exp(a_tilde(i)*s(i+1) + bi) - exp(a_tilde(i)*s(i) + bi))/a_tilde(i)
             end if
             integral = integral + gn(i)*di_exp_bi
         end if
       end if
    end do

  end function calculate_ray_integral

  ! Calculate maximum number of cells a path through the grid could take
  ! This is a loose upper bound
  function calculate_max_cells(grid) result(max_cells)
    type(space_angle_grid) :: grid
    integer :: max_cells
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
    type(space_angle_grid) :: grid
    type(optical_properties) :: iops
    double precision, dimension(:,:,:,:) :: source
    integer :: i, j, k, p
    double precision, dimension(:) :: s_array, ds, a_tilde, gn
    integer :: num_cells

    integer t
    double precision p0x, p0y, p0z
    double precision p1x, p1y, p1z
    double precision z0
    double precision s_tilde, s
    integer dir_x, dir_y, dir_z
    integer shift_x, shift_y
    integer cell_x, cell_y, cell_z
    integer edge_x, edge_y
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

    !write(*,*) 'START PATH.'
    !write(*,*) 'ijk = ', i, j, k

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
    ! (e.g., if ray is parallel to x axis, then no x intersection will occur.)
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
    ! z_factor should never be 0,
    ! because the ray is then horizontal
    ! and infinite in length.
    ! z_factor != 0 is ensured when nphi is even.
    ds_z(1:grid%z%num) = dir_z * grid%z%spacing(1:grid%z%num)/z_factor

    ! Origin point
    p0x = p1x - s_tilde * x_factor
    p0y = p1y - s_tilde * y_factor
    p0z = p1z - s_tilde * z_factor

    ! Direction of ray in each dimension. 1 => increasing. -1 => decreasing.
    dir_x = int(sgn(p1x-p0x))
    dir_y = int(sgn(p1y-p0y))

    ! Shifts
    ! Conversion from cell_inds to edge_inds
    ! merge is fortran's ternary operator
    shift_x = merge(1,0,dir_x>0)
    shift_y = merge(1,0,dir_y>0)

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

    ! First and last cells in each
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
    t = 0

    ! s is the beginning of the current cell,
    ! s_next is the end of the current cell.
    do while (s .lt. s_tilde)
       ! Move cell counter
       t = t + 1

       ! Extract function values
       a_tilde(t) = iops%abs_grid(cell_x, cell_y, cell_z)
       gn(t) = source(cell_x, cell_y, cell_z, p)

       !write(*,*) ''
       !write(*,*) 's_next_x = ', s_next_x
       !write(*,*) 's_next_y = ', s_next_y
       !write(*,*) 's_next_z = ', s_next_z
       !write(*,*) 'theta, phi =', grid%angles%theta_p(p)*180.d0/pi, grid%angles%phi_p(p)*180.d0/pi
       !write(*,*) 's = ', s, '/', s_tilde
       !write(*,*) 'cell_z =', cell_z, '/', grid%z%num
       !write(*,*) 's_next_z =', s_next_z
       !write(*,*) 'last_z =', last_z
       !write(*,*) 'new'

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

       else if(s_next_z .le. min(s_next_x, s_next_y)) then
          ! z edge is closest
          s_next = s_next_z

          ! Increment indices
          cell_z = cell_z + dir_z

          !write(*,*) 'z edge, s_next =', s_next

          ! z intersection after the one at s=s_next
          if(dir_z * (last_z - cell_z) .gt. 0) then
             ! Only look ahead if we aren't at the end
             s_next_z = s_next + ds_z(cell_z)
          else
             ! Otherwise, no need to continue.
             ! this is our final destination.
          !   exit
             s_next_z = 2*s_tilde
             !write(*,*) 'end. s_next_z =', s_next_z
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
