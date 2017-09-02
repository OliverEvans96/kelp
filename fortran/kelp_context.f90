module context_module
use sag

type optical_properties
  integer num_vsf
  double precision, dimension(:), allocatable :: vsf_array
  double precision abs_water, abs_kelp, scat_water, scat_kelp

  contains
    procedure :: load_vsf
end type optical_properties


type rope_state
   type(space_angle_grid) :: grid

   double precision, dimension(:), allocatable :: frond_lengths, frond_stds, water_speed, water_angle

contains
    procedure :: init => init_rope_state
end type rope_state


contains
  subroutine load_vsf(iops)
    class(optical_properties) :: iops

  end subroutine load_vsf

  subroutine init_rope_state(rope)
    class(rope_state) :: rope
    allocate(frond_lengths(grid%nz))
    allocate(frond_stds(grid%nz))
    allocate(water_speed(grid%nz))
    allocate(water_angle(grid%nz))
  end subroutine init_kelp_state

  ! Read array from file
  subroutine read_array(arr, n, filename)

  end subroutine read_array



end module context_module

