module kelp_context
use sag
use prob

! Point in cylindrical coordinates
type point3d
   double precision x, y, z, theta, r
 contains
   procedure :: set_cart => point_set_cart
   procedure :: set_cyl => point_set_cyl
   procedure :: cartesian_to_polar
   procedure :: polar_to_cartesian
end type point3d

type frond_shape
  double precision fs, fr, tan_alpha, alpha
contains
  procedure :: set_shape => frond_set_shape
  procedure :: calculate_angles => frond_calculate_angles
end type frond_shape

type rope_state
   integer nz
   double precision, dimension(:), allocatable :: frond_lengths, frond_stds, water_speeds, water_angles
contains
    procedure :: init => rope_init
    procedure :: deinit => rope_deinit
end type rope_state

type depth_state
   double precision frond_length, frond_std, water_speed, water_angle
contains
   procedure :: set_depth
   procedure :: length_distribution_cdf
   procedure :: angle_distribution_pdf
end type depth_state

type optical_properties
   integer num_vsf
   double precision, dimension(:), allocatable :: vsf
   double precision abs_water, abs_kelp, scat_water, scat_kelp
 contains
   procedure :: load_vsf
   procedure :: deinit => iop_deinit
end type optical_properties

contains
  subroutine point_set_cart(point, x, y, z)
    class(point3d) :: point
    double precision x, y, z
    point%x = x
    point%y = y
    point%z = z
    call point%cartesian_to_polar()
  end subroutine point_set_cart

  subroutine point_set_cyl(point, theta, r, z)
    class(point3d) :: point
    double precision theta, r, z
    point%theta = theta
    point%r = r
    point%z = z
    call point%polar_to_cartesian()
  end subroutine point_set_cyl

  subroutine polar_to_cartesian(point)
    class(point3d) :: point
    point%x = point%r*cos(point%theta)
    point%y = point%r*sin(point%theta)
  end subroutine polar_to_cartesian

  subroutine cartesian_to_polar(point)
    class(point3d) :: point
    point%r = sqrt(point%x**2 + point%y**2)
    point%theta = atan2(point%x, point%y)
  end subroutine cartesian_to_polar

  subroutine frond_set_shape(frond, fs, fr)
    class(frond_shape) frond
    double precision fs, fr
    frond%fs = fs
    frond%fr = fr
    call frond%calculate_angles()
  end subroutine frond_set_shape

  subroutine frond_calculate_angles(frond)
    class(frond_shape) frond
    frond%tan_alpha = 2*frond%fs*frond%fr / (1 + frond%fs)
    frond%alpha = tan(frond%tan_alpha)
  end subroutine

  subroutine load_vsf(iops, filename, fmtstr)
    class(optical_properties) :: iops
    character(len=*) :: filename, fmtstr
    double precision, dimension(:,:), allocatable :: tmp_2d_arr

    allocate(tmp_2d_arr(iops%num_vsf, 1))
    allocate(iops%vsf(iops%num_vsf))

    tmp_2d_arr = read_array(filename, fmtstr, iops%num_vsf, 1, 0)
    iops%vsf = tmp_2d_arr(:,1)
  end subroutine load_vsf

  subroutine rope_init(rope, grid)
    class(rope_state) :: rope
    type(space_angle_grid) :: grid

    rope%nz = grid%z%num
    allocate(rope%frond_lengths(rope%nz))
    allocate(rope%frond_stds(rope%nz))
    allocate(rope%water_speeds(rope%nz))
    allocate(rope%water_angles(rope%nz))
  end subroutine rope_init
  
  subroutine rope_deinit(rope)
    class(rope_state) rope
    deallocate(rope%frond_lengths)
    deallocate(rope%frond_stds)
    deallocate(rope%water_speeds)
    deallocate(rope%water_angles)
  end subroutine rope_deinit

  subroutine set_depth(depth, rope, depth_layer)
    class(depth_state) depth
    type(rope_state) rope
    integer depth_layer

    depth%frond_length = rope%frond_lengths(depth_layer)
    depth%frond_std = rope%frond_stds(depth_layer)
    depth%water_speed = rope%water_speeds(depth_layer)
    depth%water_angle = rope%water_angles(depth_layer)
  end subroutine set_depth

  function length_distribution_cdf(depth, L) result(output)
    ! C_L(L)
    class(depth_state) depth
    double precision L, L_mean, L_std
    double precision output

    L_mean = depth%frond_length
    L_std = depth%frond_std

    call normal_cdf(L, L_mean, L_std, output)
  end function length_distribution_cdf

  function angle_distribution_pdf(depth, theta_f) result(output)
    ! P_{\theta_f}(\theta_f)
    class(depth_state) depth
    double precision theta_f, v_w, theta_w
    double precision output

    v_w = depth%water_speed
    theta_w = depth%water_angle

    call von_mises_pdf(theta_f, theta_w, v_w, output)
  end function angle_distribution_pdf

  subroutine iop_deinit(iops)
    class(optical_properties) iops
    deallocate(iops%vsf)
  end subroutine iop_deinit

end module kelp_context
