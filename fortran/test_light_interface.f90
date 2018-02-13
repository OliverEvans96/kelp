program test_light_interface
  use light_interface_module
    double precision abs_kelp
    double precision abs_water
    double precision scat_kelp
    double precision scat_water
    integer num_vsf
    character(len=56) vsf_file
    double precision solar_zenith
    double precision solar_azimuthal
    integer num_si
    double precision, dimension(:), allocatable ::  si_area
    integer, dimension(:), allocatable ::  si_ind
    double precision frond_thickness
    double precision frond_aspect_ratio
    double precision frond_shape_ratio
    double precision, dimension(:), allocatable ::  current_speeds
    double precision, dimension(:), allocatable ::  current_angles
    double precision rope_spacing
    double precision depth_spacing
    integer nx
    integer ny
    integer nz
    integer ntheta
    integer nphi
    integer num_scatters
    double precision, dimension(:), allocatable ::  pre_kelp_irrad
    double precision, dimension(:), allocatable ::  post_kelp_irrad

    integer k

    abs_kelp = 0.d0
    abs_water = 0.d0
    scat_kelp = 0.d0
    scat_water = 0.d0
    num_vsf = 0
    vsf_file = '/home/oliver/academic/research/kelp/data/vsf/nuc_vsf.txt'
    solar_zenith = 0.d0
    solar_azimuthal = 0.d0
    num_si = 10

    nx = 10
    ny = 10
    nz = 10
    ntheta = 10
    nphi = 10
    num_scatters = 1

    allocate(si_area(nz))
    allocate(si_ind(nz))
    allocate(current_speeds(nz))
    allocate(current_angles(nz))
    allocate(pre_kelp_irrad(nz))
    allocate(post_kelp_irrad(nz))

    do k=1, nz
       si_area(k) = 0.d0
       si_ind(k) = 0.d0
       current_speeds(k) = 0.d0
       current_angles(k) = 0.d0
       pre_kelp_irrad(k) = 0.d0
       post_kelp_irrad(k) = 0.d0
    end do

  call full_light_calculations( &
       abs_kelp,        &
       abs_water,       &
       scat_kelp,      &
       scat_water,      &
       num_vsf,         &
       vsf_file,        &
                                ! SUNLIGHT       &
       solar_zenith,    &
       solar_azimuthal, &
                                ! KELP           &
       num_si,          &
       si_area,         &
       si_ind,          &
       frond_thickness, &
       frond_aspect_ratio,    &
       frond_shape_ratio,     &
                                ! WATER CURRENT  &
       current_speeds,  &
       current_angles,   &
                                ! SPACING        &
       rope_spacing,    &
       depth_spacing,   &
                                ! SOLVER PARAMETERS
       nx, &
       ny, &
       nz, &
       ntheta, &
       nphi, &
       num_scatters,    &
                                ! LIGHT WITHOUT KELP
       pre_kelp_irrad,  &
                                ! FINAL RESULT   &
       post_kelp_irrad)

end program test_light_interface
