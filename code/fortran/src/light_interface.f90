module light_interface
  use rte3d
  use kelp3d
  use asymptotics
  implicit none

contains
  subroutine full_light_calculations( &
    ! OPTICAL PROPERTIES
    absorptance_kelp, & ! NOT THE SAME AS ABSORPTION COEFFICIENT
    abs_water, &
    scat, &
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
    ! FINAL RESULTS
    perceived_irrad, &
    avg_irrad)

    implicit none

    ! OPTICAL PROPERTIES
    integer, intent(in) :: nx, ny, nz, ntheta, nphi
    ! Absorption and scattering coefficients
    double precision, intent(in) :: absorptance_kelp, scat
    double precision, dimension(nz), intent(in) :: abs_water
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
    double precision, dimension(nz, num_si), intent(in) :: si_ind
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
    real, dimension(nz), intent(out) :: avg_irrad, perceived_irrad

    !-------------!

    double precision xmin, xmax, ymin, ymax, zmin, zmax
    character(len=5), parameter :: fmtstr = 'E13.4'
    !double precision, dimension(num_vsf) :: vsf_angles, vsf_vals
    double precision max_rad, decay
    integer quadrature_degree

    type(space_angle_grid) grid
    type(optical_properties) iops
    type(light_state) light
    type(rope_state) rope
    type(frond_shape) frond
    type(boundary_condition) bc

    double precision, dimension(:), allocatable :: pop_length_means, pop_length_stds
    ! Number of fronds in each depth layer
    double precision, dimension(:), allocatable :: num_fronds
    double precision, dimension(:,:,:), allocatable :: p_kelp

    write(*,*) 'Light calculation'

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

    write(*,*) 'Grid'
    call grid%set_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
    call grid%set_num(nx, ny, nz, ntheta, nphi)
    call grid%init()
    !call grid%set_uniform_spacing_from_num()
    call grid%z%set_spacing_array(depth_spacing)

    call rope%init(grid)

    write(*,*) 'Rope'
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

    write(*,*) 'frond_lengths  =', rope%frond_lengths
    write(*,*) 'frond_stds  =', rope%frond_stds
    write(*,*) 'num_fronds  =', rope%num_fronds
    write(*,*) 'water_speeds  =', rope%water_speeds
    write(*,*) 'water_angles  =', rope%water_angles

    write(*,*) 'Frond'
    ! INIT FROND
    call frond%set_shape(frond_shape_ratio, frond_aspect_ratio, frond_thickness)
    ! CALCULATE KELP
    quadrature_degree = 5
    call calculate_kelp_on_grid(grid, p_kelp, frond, rope, quadrature_degree)
    ! INIT IOPS
    iops%num_vsf = num_vsf
    call iops%init(grid)
    write(*,*) 'IOPs'
    iops%abs_kelp = absorptance_kelp / frond_thickness
    iops%abs_water = abs_water
    iops%scat = scat

    !write(*,*) 'iop init'
    !iops%vsf_angles = vsf_angles
    !iops%vsf_vals = vsf_vals
    call iops%load_vsf(vsf_file, fmtstr)

    ! load_vsf already calls calc_vsf_on_grid
    !call iops%calc_vsf_on_grid()
    call iops%calculate_coef_grids(p_kelp)

    !write(*,*) 'BC'
    max_rad = 1.d0 ! Doesn't matter because we'll rescale
    decay = 1.d0 ! Does matter, but maybe not much. Determines drop-off from angle
    call bc%init(grid, solar_zenith, solar_azimuthal, decay, max_rad)
    ! Rescale surface radiance to match surface irradiance
    bc%bc_grid = bc%bc_grid * surface_irrad / grid%angles%integrate_points(bc%bc_grid)

    write(*,*) 'bc'
    write(*,*) bc%bc_grid

    ! write(*,*) 'bc'
    ! do i=1, grid%y%num
    !     write(*,'(10F15.3)') bc%bc_grid(i,:)
    ! end do

    call light%init_grid(grid)

    write(*,*) 'Scatter'
    call calculate_light_with_scattering(grid, bc, iops, light%radiance, num_scatters)

    write(*,*) 'Irrad'
    call light%calculate_irradiance()

    ! Calculate output variables
    call calculate_average_irradiance(grid, light, avg_irrad)
    call calculate_perceived_irradiance(grid, p_kelp, &
         perceived_irrad, light%irradiance)

    !write(*,*) 'vsf_angles = ', iops%vsf_angles
    !write(*,*) 'vsf_vals = ', iops%vsf_vals
    !write(*,*) 'vsf norm  = ', grid%integrate_angle_2d(iops%vsf(1,1,:,:))

    ! write(*,*) 'abs_water = ', abs_water
    ! write(*,*) 'scat_water = ', scat_water
    write(*,*) 'kelp '
    write(*,*) p_kelp(:,:,:)
    write(*,*) 'ft =', frond%ft

    write(*,*) 'irrad'
    write(*,*) light%irradiance

    write(*,*) 'avg_irrad = ', avg_irrad
    write(*,*) 'perceived_irrad = ', perceived_irrad

    write(*,*) 'deinit'
    call bc%deinit()
    !write(*,*) 'a'
    call iops%deinit()
    !write(*,*) 'b'
    call light%deinit()
    !write(*,*) 'c'
    call rope%deinit()
    !write(*,*) 'd'
    call grid%deinit()
    !write(*,*) 'e'

    deallocate(pop_length_means)
    deallocate(pop_length_stds)
    deallocate(num_fronds)
    deallocate(p_kelp)

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
    double precision, dimension(nz, num_si), intent(in) :: si_ind
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
          si_length(i) = sqrt(2.d0*frond_aspect_ratio*si_area(k,i))
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

  subroutine calculate_average_irradiance(grid, light, avg_irrad)
    type(space_angle_grid) grid
    type(light_state) light
    real, dimension(:) :: avg_irrad
    integer k, nx, ny, nz

    nx = grid%x%num
    ny = grid%y%num
    nz = grid%z%num

    do k=1, nz
       avg_irrad(k) = real(sum(light%irradiance(:,:,k)) / nx / ny)
    end do
  end subroutine calculate_average_irradiance

  subroutine calculate_perceived_irradiance(grid, p_kelp, &
       perceived_irrad, irradiance)
    type(space_angle_grid) grid
    double precision, dimension(:,:,:) :: p_kelp
    real, dimension(:) :: perceived_irrad
    double precision, dimension(:,:,:) :: irradiance
    double precision total_kelp
    integer center_i1, center_i2, center_j1, center_j2

    integer k

    ! Calculate the average irradiance experienced over the frond.
    ! Has same units as irradiance.
    ! If no kelp, then just take the irradiance at the center
    ! of the grid.
    do k=1, grid%z%num
       total_kelp = sum(p_kelp(:,:,k))
       if(total_kelp .eq. 0) then
          center_i1 = int(ceiling(grid%x%num / 2.d0))
          center_j1 = int(ceiling(grid%y%num / 2.d0))
          ! For even grid, use average of center two cells
          ! For odd grid, just use center cell
          if(mod(grid%x%num, 2) .eq. 0) then
             center_i2 = center_i1 + 1
          else
             center_i2 = center_i1
          end if
          if(mod(grid%y%num, 2) .eq. 0) then
             center_j2 = center_j1 + 1
          else
             center_j2 = center_j1
          end if


          ! Irradiance at the center of the grid (at the rope)
          perceived_irrad(k) = real(sum(irradiance( &
               center_i1:center_i2, &
               center_j1:center_j2, k)) &
               / ((center_i2-center_i1+1) * (center_j2-center_j1+1)))
       else
          ! Average irradiance weighted by kelp distribution
          perceived_irrad(k) = real( &
              sum(p_kelp(:,:,k)*irradiance(:,:,k)) &
              / total_kelp)
       end if
    end do

  end subroutine calculate_perceived_irradiance

end module light_interface
