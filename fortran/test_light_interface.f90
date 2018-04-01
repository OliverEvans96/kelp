program test_light_interface
  use light_interface_module

  implicit none
    double precision abs_kelp, scat
    double precision, dimension(:), allocatable :: abs_water
    integer num_vsf
    character(len=56) vsf_file
    double precision solar_zenith
    double precision solar_azimuthal
    double precision surface_irrad
    integer num_si
    double precision, dimension(:,:), allocatable ::  si_area
    double precision, dimension(:,:), allocatable ::  si_ind
    double precision frond_thickness
    double precision frond_aspect_ratio
    double precision frond_shape_ratio
    double precision, dimension(:), allocatable ::  current_speeds
    double precision, dimension(:), allocatable ::  current_angles
    double precision rope_spacing
    double precision, dimension(:), allocatable ::  depth_spacing
    integer nx
    integer ny
    integer nz
    integer ntheta
    integer nphi
    integer num_scatters
    real, dimension(:,:), allocatable :: available_light
    real, dimension(:), allocatable ::  avg_irrad

    integer i, k

    abs_kelp = 1.d0
    scat = 0.01d0

    num_vsf = 55
    vsf_file = '/home/oliver/academic/research/kelp/data/vsf/nuc_vsf.txt'

    solar_zenith = 0.d0
    solar_azimuthal = 0.d0
    surface_irrad = 1.d0

    rope_spacing = 10.d0

    num_si = 10
    nx = 1
    ny = 1
    nz = 5
    ntheta = 2
    nphi = 2
    num_scatters = 1

    allocate(abs_water(nz))
    allocate(si_area(nz, num_si))
    allocate(si_ind(nz, num_si))
    allocate(current_speeds(nz))
    allocate(current_angles(nz))
    allocate(available_light(nz, num_si))
    allocate(avg_irrad(nz))
    allocate(depth_spacing(nz))

    ! Use kelp iops for water so that medium is homogeneous
    do k=1, nz
       abs_water(k) = abs_kelp

       do i=1, num_si
          si_area(k, i) = 3.d0
          si_ind(k, i) = 10
       end do

       current_speeds(k) = 0.d0
       current_angles(k) = 0.d0

       avg_irrad(k) = 0.d0

       depth_spacing(k) = 1.d0
    end do

    ! Frond properties
    frond_aspect_ratio = 2.d0
    frond_shape_ratio = 0.5d0
    frond_thickness = 0.05d0

    call full_light_calculations( &
         ! OPTICAL PROPERTIES
         abs_kelp, &
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
         available_light, &
         avg_irrad)

  deallocate(abs_water)
  deallocate(si_area)
  deallocate(si_ind)
  deallocate(current_speeds)
  deallocate(current_angles)
  deallocate(available_light)
  deallocate(avg_irrad)
  deallocate(depth_spacing)


end program test_light_interface
