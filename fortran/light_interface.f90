module light_interface_module
  use rte3d
  use kelp3d
  use asymptotics
  implicit none

contains
  subroutine full_light_calculations( &
    ! OPTICAL PROPERTIES
    abs_kelp, &
    scat_kelp, &
    abs_water, &
    scat_water, &
    num_vsf, &
    vsf_file, &
    ! SUNLIGHT
    solar_zenith, &
    solar_azimuthal, &
    surface_irrad, &
    ! KELP &
    num_si, &
    si_area, &
    si_ind, &
    frond_thickness, &
    frond_aspect_ratio, &
    frond_shape_ratio, &
    ! WATER CURRENT
    current_speeds, &
    current_angles, &
    ! SPACING
    rope_spacing, &
    depth_spacing, &
    ! SOLVER PARAMETERS
    nx, &
    ny, &
    nz, &
    ntheta, &
    nphi, &
    num_scatters, &
    ! FINAL RESULT
    irrad)

    implicit none

    ! OPTICAL PROPERTIES
    integer, intent(in) :: nx, ny, nz, ntheta, nphi
    ! Absorption and scattering coefficients
    double precision, intent(in) :: abs_kelp
    double precision, intent(in) :: scat_kelp
    double precision, dimension(nz), intent(in) :: abs_water
    double precision, dimension(nz), intent(in) :: scat_water
    ! Volume scattering function
    integer, intent(in) :: num_vsf
    character(len=*) :: vsf_file
    !double precision, dimension(num_vsf), intent(int) :: vsf_angles
    !double precision, dimension(num_vsf), intent(int) :: vsf_vals

    ! SUNLIGHT
    double precision, intent(in) :: solar_zenith
    double precision, intent(in) :: solar_azimuthal
    double precision, intent(in) :: surface_irrad

    ! KELP
    ! Number of Superindividuals in each depth level
    integer, intent(in) :: num_si
    ! si_area(i,j) = area of superindividual j at depth i
    double precision, dimension(nz, num_si), intent(in) :: si_area
    ! si_ind(i,j) = number of inidividuals represented by superindividual j at depth i
    integer, dimension(nz, num_si), intent(in) :: si_ind
    ! Thickness of each frond
    double precision, intent(in) :: frond_thickness
    ! Ratio of length to width (0,infty)
    double precision, intent(in) :: frond_aspect_ratio
    ! Rescaled position of greatest width (0=base, 1=tip)
    double precision, intent(in) :: frond_shape_ratio

    ! WATER CURRENT
    double precision, dimension(nz), intent(in) :: current_speeds
    double precision, dimension(nz), intent(in) :: current_angles

    ! SPACING
    double precision, intent(in) :: rope_spacing
    double precision, dimension(nz), intent(in) :: depth_spacing
    ! SOLVER PARAMETERS
    integer, intent(in) :: num_scatters

    ! FINAL RESULT
    real, dimension(nz), intent(out) :: irrad

    !-------------!

    double precision xmin, xmax, ymin, ymax, zmin, zmax
    character(len=5), parameter :: fmtstr = 'E13.4'
    !double precision, dimension(num_vsf) :: vsf_angles, vsf_vals
    double precision max_rad, decay
    integer quadrature_degree

    type(space_angle_grid) grid
    type(rte_mat) mat
    type(optical_properties) iops
    type(light_state) light
    type(rope_state) rope
    type(frond_shape) frond
    type(boundary_condition) bc

    double precision, dimension(:), allocatable :: pop_length_means, pop_length_stds
    ! Number of fronds in each depth layer
    double precision, dimension(:), allocatable :: num_fronds

    integer k

    double precision, dimension(:,:,:), allocatable :: p_kelp

    allocate(pop_length_means(nz))
    allocate(pop_length_stds(nz))
    allocate(num_fronds(nz))
    allocate(p_kelp(nx, ny, nz))

    xmin = -rope_spacing/2
    xmax = rope_spacing/2

    ymin = -rope_spacing/2
    ymax = rope_spacing/2

    zmin = 0.d0
    zmax = sum(depth_spacing)

    ! INIT GRID
    !write(*,*) 'Grid'
    grid%x%minval = xmin
    grid%x%maxval = xmax
    grid%x%num = nx

    grid%y%minval = ymin
    grid%y%maxval = ymax
    grid%y%num = ny

    grid%z%minval = zmin
    grid%z%maxval = zmax
    grid%z%num = nz

    grid%theta%num = ntheta
    grid%phi%num = nphi

    call grid%init()
    !call grid%set_uniform_spacing_from_num()
    call grid%z%set_spacing_array(depth_spacing)

    call rope%init(grid)

    ! Calculate kelp distribution
    call calculate_length_dist_from_superinds( &
    nz, &
    num_si, &
    si_area, &
    si_ind, &
    frond_aspect_ratio, &
    num_fronds, &
    pop_length_means, &
    pop_length_stds)

    rope%frond_lengths = pop_length_means
    rope%frond_stds = pop_length_stds
    rope%num_fronds = num_fronds
    rope%water_speeds = current_speeds
    rope%water_angles = current_angles

    ! INIT FROND
    call frond%set_shape(frond_shape_ratio, frond_aspect_ratio, frond_thickness)

    ! CALCULATE KELP
    quadrature_degree = 5
    call calculate_kelp_on_grid(grid, p_kelp, frond, rope, quadrature_degree)
    ! INIT IOPS
    !write(*,*) 'IOPs'
    iops%abs_kelp = abs_kelp
    iops%scat_kelp = scat_kelp
    iops%abs_water = abs_water
    iops%scat_water = scat_water
    iops%num_vsf = num_vsf


    !write(*,*) 'iop init'
    call iops%init(grid)
    !iops%vsf_angles = vsf_angles
    !iops%vsf_vals = vsf_vals
    call iops%load_vsf(vsf_file, fmtstr)

    ! load_vsf already calls calc_vsf_on_grid
    !call iops%calc_vsf_on_grid()
    call iops%calculate_coef_grids(p_kelp)

    !write(*,*) 'BC'
    max_rad = 1 ! Doesn't matter because we'll rescale
    decay = 1 ! Does matter, but maybe not much. Determines drop-off from angle
    call bc%init(grid, solar_zenith, solar_azimuthal, decay, max_rad)
    ! Rescale surface radiance to match surface irradiance
    bc%bc_grid = bc%bc_grid * surface_irrad / grid%integrate_angle_2d(bc%bc_grid)

    call light%init_grid(grid)

    !write(*,*) 'Scatter'
    call calculate_light_with_scattering(grid, bc, iops, light%radiance, num_scatters)

    !write(*,*) 'Irrad'
    call light%calculate_irradiance()

    ! Calculate average irradiance
    do k=1, nz
       irrad(k) = real(sum(light%irradiance(:,:,k)) / nx / ny)
    end do

    !write(*,*) 'deinit'
    call bc%deinit()
    call iops%deinit()
    call light%deinit()
    call rope%deinit()
    call grid%deinit()

    deallocate(pop_length_means)
    deallocate(pop_length_stds)
    deallocate(num_fronds)
    deallocate(p_kelp)

    write(*,*) 'abs_water = ', abs_water
    write(*,*) 'scat_water = ', scat_water
    write(*,*) 'irrad = ', irrad

    !write(*,*) 'done'
  end subroutine full_light_calculations

  subroutine calculate_length_dist_from_superinds( &
    nz, &
    num_si, &
    si_area, &
    si_ind, &
    frond_aspect_ratio, &
    num_fronds, &
    pop_length_means, &
    pop_length_stds)

    implicit none

    ! Number of depth levels
    integer, intent(in) :: nz
    ! Number of Superindividuals in each depth level
    integer, intent(in) :: num_si
    ! si_area(i,j) = area of superindividual j at depth i
    double precision, dimension(nz, num_si), intent(in) :: si_area
    ! si_area(i,j) = number of inidividuals represented by superindividual j at depth i
    integer, dimension(nz, num_si), intent(in) :: si_ind
    double precision, intent(in) :: frond_aspect_ratio

    double precision, dimension(nz), intent(out) :: num_fronds
    ! Population mean area at each depth level
    double precision, dimension(nz), intent(out) :: pop_length_means
    ! Population area standard deviation at each depth level
    double precision, dimension(nz), intent(out) :: pop_length_stds

    !---------------!

    integer i, k
    ! Numerators for mean and std
    double precision mean_num, std_num
    ! Convert area to length
    double precision, dimension(num_si) :: si_length

    do k=1, nz
       mean_num = 0.d0
       std_num = 0.d0
       num_fronds(k) = 0

       do i=1, num_si
          si_length(i) = sqrt(2*frond_aspect_ratio*si_area(k,i))
          mean_num = mean_num + si_length(i)
          num_fronds(k) = num_fronds(k) + si_ind(k,i)
       end do

       pop_length_means(k) = mean_num / num_fronds(k)

       do i=1, num_si
          std_num = std_num + (si_length(i) - pop_length_means(k))**2
       end do

       pop_length_stds(k) = std_num / (num_fronds(k) - 1)

    end do

  end subroutine calculate_length_dist_from_superinds

end module light_interface_module
