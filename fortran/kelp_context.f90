module kelp_context
use sag
use utils

type optical_properties
  integer num_vsf
  double precision, dimension(:), allocatable :: vsf
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
  subroutine load_vsf(iops, filename, fmtstr)
    class(optical_properties) :: iops
    character(len=*) :: filename, fmtstr
    double precision, dimension(:,:), allocatable :: tmp_2d_arr

    allocate(tmp_2d_arr(iops%num_vsf, 1))
    allocate(iops%vsf(iops%num_vsf))

    tmp_2d_arr = read_array(filename, fmtstr, iops%num_vsf, 1, 0)
    iops%vsf = tmp_2d_arr(:,1)
  end subroutine load_vsf

  subroutine init_rope_state(rope)
    class(rope_state) :: rope
    allocate(rope%frond_lengths(rope%grid%z%num))
    allocate(rope%frond_stds(rope%grid%z%num))
    allocate(rope%water_speed(rope%grid%z%num))
    allocate(rope%water_angle(rope%grid%z%num))
  end subroutine init_rope_state

end module kelp_context
