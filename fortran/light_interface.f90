module light_interface
contains
  subroutine full_light_calculations(
    ! OPTICAL PROPERTIES
    abs_kelp,
    abs_water,
    scat_kelp,
    scat_water,
    num_vsf,
    vsf_angles,
    vsf_vals,
    ! SUNLIGHT
    solar_zenith,
    solar_azimuthal,
    surface_irrad,
    ! KELP
    num_si,
    si_area,
    si_ind,
    frond_thickness,
    frond_aspect,
    frond_shape,
    ! WATER CURRENT
    current_speed,
    current_angle
    ! SPACING
    rope_spacing,
    depth_spacing,
    num_depth_levels,
    ! SOLVER PARAMETERS
    spatial_res,
    angular_res,
    num_scatters,
    ! FINAL RESULT
    absorbed_light)

    ! OPTICAL PROPERTIES
    ! Absorption and scattering coefficients
    double precision, intent(in) :: abs_kelp
    double precision, intent(in) :: abs_water
    double precision, intent(in) :: scat_kelp
    double precision, intent(in) :: scat_water
    ! Volume scattering function
    integer, intent(in) :: num_vsf
    double precision, dimension(num_vsf), intent(int) :: vsf_angles
    double precision, dimension(num_vsf), intent(int) :: vsf_vals

    ! SUNLIGHT
    double precision, intent(in) :: solar_zenith, solar_azimuthal
    double precision, intent(in) :: surface_irrad

    ! KELP
    ! Number of Superindividuals in each depth level
    integer, intent(in) :: num_si
    ! si_area(i,j) = area of superindividual j at depth i
    double precision, dimension(num_depth_levels, num_si), intent(in) :: si_area
    ! si_area(i,j) = number of inidividuals represented by superindividual j at depth i
    integer, dimension(num_depth_levels, num_si), intent(in) :: si_area
    ! Thickness of each frond
    double precision, intent(in) :: frond_thickness
    ! Ratio of length to width (0,infty)
    double precision, intent(in) :: frond_aspect
    ! Rescaled position of greatest width (0=base, 1=tip)
    double precision, intent(in) :: frond_shape


    ! WATER CURRENT
    double precision, dimension(num_depth_levels) :: current_speed
    double precision, dimension(num_depth_levels) :: current_angle

    ! SPACING
    double precision, intent(in) :: rope_spacing
    double precision, intent(in) :: depth_spacing
    integer, intent(in) :: num_depth_levels

    ! SOLVER PARAMETERS
    double precision, intent(in) :: spatial_res
    double precision, intent(in) :: angular_res
    integer, intent(in) :: num_scatters

    ! FINAL RESULT
    double precision, dimension(num_depth_levels) :: absorbed_light
  end subroutine full_light_calculations

  subroutine calculate_length_dist_from_superinds(
    num_depth_levels,
    num_si,
    si_area,
    si_ind,
    frond_aspect,
    num_fronds,
    pop_length_means,
    pop_length_stds)

    ! Number of depth levels
    integer, intent(in) :: num_depth_levels
    ! Number of Superindividuals in each depth level
    integer, intent(in) :: num_si
    ! si_area(i,j) = area of superindividual j at depth i
    double precision, dimension(num_depth_levels, num_si), intent(in) :: si_area
    ! si_area(i,j) = number of inidividuals represented by superindividual j at depth i
    integer, dimension(num_depth_levels, num_si), intent(in) :: si_area
    double precision, intent(in) :: frond_aspect

    integer, dimension(num_depth_levels), intent(out) :: num_fronds
    ! Population mean area at each depth level
    double precision, dimension(num_depth_levels), intent(out) :: pop_length_means
    ! Population area standard deviation at each depth level
    double precision, dimension(num_depth_levels), intent(out) :: pop_length_stds

    !---------------!

    integer i, k
    ! Numerators for mean and std
    double precision mean_num, std_num
    ! Convert area to length
    double precision, dimension(num_si) :: si_length

    ! PROBABLY NEED EXPLICIT ALLOCATE STATEMENTS HERE

    do k=1, num_depth_levels
       mean_num = 0.d0
       std_num = 0.d0
       num_fronds(k) = 0

       do i=1, num_si
          si_length(i) = sqrt(2*frond_aspect*si_area(k,i))
          mean_num = mean_num + si_length(i)
          num_fronds(k) = num_fronds + si_ind(k,i)
       end do

       pop_mean(k) = mean_num / dble(num_fronds(k))

       do i=1, num_si
          std_num = std_num + (si_length(i) = pop_mean(k))**2
       end do

       pop_std(k) = std_num / dble(num_fronds(k) - 1)

    end do

  end subroutine calculate_area_dist_from_superinds

  subroutine calculate_

end module light_interface
