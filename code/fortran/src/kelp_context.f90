module kelp_context
use sag
use prob
implicit none

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
  double precision fs, fr, tan_alpha, alpha, ft
contains
  procedure :: set_shape => frond_set_shape
  procedure :: calculate_angles => frond_calculate_angles
end type frond_shape

type rope_state
   integer nz
   double precision, dimension(:), allocatable :: frond_lengths, frond_stds, num_fronds, water_speeds, water_angles
contains
    procedure :: init => rope_init
    procedure :: deinit => rope_deinit
end type rope_state

type depth_state
   double precision frond_length, frond_std, num_fronds, water_speeds, water_angles, depth
   integer depth_layer
contains
   procedure :: set_depth
   procedure :: length_distribution_cdf
   procedure :: angle_distribution_pdf
   procedure :: angle_mod
   procedure :: angle_diff_2d
end type depth_state

type optical_properties
   integer num_vsf
   type(space_angle_grid) grid
   double precision, dimension(:), allocatable :: vsf_angles, vsf_vals, vsf_cos
   double precision, dimension(:), allocatable :: abs_water, scat_water
   double precision abs_kelp, vsf_scat_coef, scat
   ! On x, y, z grid - including water & kelp.
   double precision, dimension(:,:,:), allocatable :: abs_grid
   ! On theta, phi, theta_prime, phi_prime grid
   double precision, dimension(:,:), allocatable :: vsf, vsf_integral
 contains
   procedure :: init => iop_init
   procedure :: calculate_coef_grids
   procedure :: load_vsf
   procedure :: eval_vsf
   procedure :: calc_vsf_on_grid
   procedure :: deinit => iop_deinit
   procedure :: vsf_from_function
end type optical_properties

type boundary_condition
   double precision max_rad, decay, theta_s, phi_s
   type(space_angle_grid) grid
   double precision, dimension(:), allocatable :: bc_grid
 contains
   procedure :: bc_gaussian
   procedure :: init => bc_init
   procedure :: deinit => bc_deinit
end type boundary_condition

contains

  function bc_gaussian(bc, theta, phi)
    class(boundary_condition) bc
    double precision theta, phi, diff
    double precision bc_gaussian
    diff = angle_diff_3d(theta, phi, bc%theta_s, bc%phi_s)
    bc_gaussian = bc%max_rad * exp(-bc%decay * diff)
  end function bc_gaussian

  subroutine bc_init(bc, grid, theta_s, phi_s, decay, max_rad)
    class(boundary_condition) bc
    type(space_angle_grid) grid
    double precision theta_s, phi_s, decay, max_rad
    integer p
    double precision theta, phi

    allocate(bc%bc_grid(grid%angles%nomega))

    bc%theta_s = theta_s
    bc%phi_s = phi_s
    bc%decay = decay
    bc%max_rad = max_rad

    ! Only set BC for downwelling light
    do p=1, grid%angles%nomega/2
       theta = grid%angles%theta_p(p)
       phi = grid%angles%phi_p(p)
       bc%bc_grid(p) = bc%bc_gaussian(theta, phi)
    end do
    ! Zero upwelling light specified at surface
    do p=grid%angles%nomega/2+1, grid%angles%nomega
       bc%bc_grid(p) = 0.d0
    end do

  end subroutine bc_init

  subroutine bc_deinit(bc)
    class(boundary_condition) bc
    deallocate(bc%bc_grid)
    end subroutine

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
    point%theta = atan2(point%y, point%x)
  end subroutine cartesian_to_polar

  subroutine frond_set_shape(frond, fs, fr, ft)
    class(frond_shape) frond
    double precision fs, fr, ft
    frond%fs = fs
    frond%fr = fr
    frond%ft = ft
    call frond%calculate_angles()
  end subroutine frond_set_shape

  subroutine frond_calculate_angles(frond)
    class(frond_shape) frond
    frond%tan_alpha = 2.d0*frond%fs*frond%fr / (1.d0 + frond%fs)
    frond%alpha = atan(frond%tan_alpha)
  end subroutine

  subroutine iop_init(iops, grid)
    class(optical_properties) iops
    type(space_angle_grid) grid

    iops%grid = grid

    ! Assume that these are preallocated and passed to function
    ! Nevermind, don't assume this.
    allocate(iops%abs_water(grid%z%num))
    allocate(iops%scat_water(grid%z%num))

    ! Assume that these must be allocated here
    allocate(iops%vsf_angles(iops%num_vsf))
    allocate(iops%vsf_vals(iops%num_vsf))
    allocate(iops%vsf_cos(iops%num_vsf))
    allocate(iops%vsf(grid%angles%nomega,grid%angles%nomega))
    allocate(iops%vsf_integral(grid%angles%nomega,grid%angles%nomega))
    allocate(iops%abs_grid(grid%x%num, grid%y%num, grid%z%num))
  end subroutine iop_init

  subroutine calculate_coef_grids(iops, p_kelp)
    class(optical_properties) iops
    double precision, dimension(:,:,:) :: p_kelp

    integer k

    ! Allow water IOPs to vary over depth
    do k=1, iops%grid%z%num
      iops%abs_grid(:,:,k) = (iops%abs_kelp - iops%abs_water(k)) * p_kelp(:,:,k) + iops%abs_water(k)
   end do

  end subroutine calculate_coef_grids


  subroutine load_vsf(iops, filename, fmtstr)
    class(optical_properties) :: iops
    character(len=*) :: filename, fmtstr
    double precision, dimension(:,:), allocatable :: tmp_2d_arr
    integer num_rows, num_cols, skiplines_in

    ! First column is the angle at which the measurement is taken
    ! Second column is the value of the VSF at that angle
    num_rows = iops%num_vsf
    num_cols = 2
    skiplines_in = 1 ! Ignore comment on first line

    allocate(tmp_2d_arr(num_rows, num_cols))

    tmp_2d_arr = read_array(filename, fmtstr, num_rows, num_cols, skiplines_in)
    iops%vsf_angles = tmp_2d_arr(:,1)
    iops%vsf_vals = tmp_2d_arr(:,2)

    ! write(*,*) 'vsf_angles = ', iops%vsf_angles
    ! write(*,*) 'vsf_vals = ', iops%vsf_vals

    ! Pre-evaluate for all pair of angles
    call iops%calc_vsf_on_grid()
  end subroutine load_vsf

  function eval_vsf(iops, theta)
    class(optical_properties) iops
    double precision theta
    double precision eval_vsf
    ! No need to set vsf(0) = 0.
    ! It's the area under the curve that matters, not the value.
    eval_vsf = interp(theta, iops%vsf_angles, iops%vsf_vals, iops%num_vsf)

  end function eval_vsf

  subroutine rope_init(rope, grid)
    class(rope_state) :: rope
    type(space_angle_grid) :: grid

    rope%nz = grid%z%num
    allocate(rope%frond_lengths(rope%nz))
    allocate(rope%frond_stds(rope%nz))
    allocate(rope%water_speeds(rope%nz))
    allocate(rope%water_angles(rope%nz))
    allocate(rope%num_fronds(rope%nz))
  end subroutine rope_init

  subroutine rope_deinit(rope)
    class(rope_state) rope
    deallocate(rope%frond_lengths)
    deallocate(rope%frond_stds)
    deallocate(rope%water_speeds)
    deallocate(rope%water_angles)
    deallocate(rope%num_fronds)
  end subroutine rope_deinit

  subroutine set_depth(depth, rope, grid, depth_layer)
    class(depth_state) depth
    type(rope_state) rope
    type(space_angle_grid) grid
    integer depth_layer

    depth%frond_length = rope%frond_lengths(depth_layer)
    depth%frond_std = rope%frond_stds(depth_layer)
    depth%water_speeds = rope%water_speeds(depth_layer)
    depth%water_angles = rope%water_angles(depth_layer)
    depth%num_fronds = rope%num_fronds(depth_layer)
    depth%depth_layer = depth_layer
    depth%depth = grid%z%vals(depth_layer)
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
    double precision diff

    v_w = depth%water_speeds
    theta_w = depth%water_angles

    ! von_mises_pdf is only defined on [-pi, pi]
    ! So take difference of angles and input into
    ! von_mises dist. centered & x=0.

    diff = depth%angle_diff_2d(theta_f, theta_w)

    call von_mises_pdf(diff, 0.d0, v_w, output)
  end function angle_distribution_pdf

  function angle_mod(depth, theta) result(mod_theta)
    ! Shift theta to the interval [-pi, pi]
    ! which is where von_mises_pdf is defined.

    class(depth_state) depth
    double precision theta, mod_theta

    mod_theta = mod(theta + pi, 2.d0*pi) - pi
  end function angle_mod

  function angle_diff_2d(depth, theta1, theta2) result(diff)
    ! Shortest difference between two angles which may be
    ! in different periods.
    class(depth_state) depth
    double precision theta1, theta2, diff
    double precision modt1, modt2

    ! Shift to [0, 2*pi]
    modt1 = mod(theta1, 2*pi)
    modt2 = mod(theta2, 2*pi)

    ! https://gamedev.stackexchange.com/questions/4467/comparing-angles-and-working-out-the-difference

    diff = pi - abs(abs(modt1-modt2) - pi)
  end function angle_diff_2d

  function angle_diff_3d(theta, phi, theta_prime, phi_prime) result(diff)
    ! Angle between two vectors in spherical coordinates
    double precision theta, phi, theta_prime, phi_prime
    double precision alpha, diff

    ! Faster, but produces lots of NaNs (at least in Python)
    !alpha = sin(theta)*sin(theta_prime)*cos(theta-theta_prime) + cos(phi)*cos(phi_prime)


    ! Slower, but more accurate
    alpha = (sin(phi)*sin(phi_prime) &
      * (cos(theta)*cos(theta_prime) + sin(theta)*sin(theta_prime)) &
      + cos(phi)*cos(phi_prime))

    ! Avoid out-of-bounds errors due to rounding
    alpha = min(1.d0, alpha)
    alpha = max(-1.d0, alpha)

    diff = acos(alpha)
  end function angle_diff_3d

  subroutine vsf_from_function(iops, func)
    class(optical_properties) iops
    double precision, external :: func
    integer i
    type(angle_dim) :: angle1d

    call angle1d%set_bounds(-1.d0, 1.d0)
    call angle1d%set_num(iops%num_vsf)
    call angle1d%assign_legendre()

    iops%vsf_angles = acos(angle1d%vals)
    do i=1, iops%num_vsf
       iops%vsf_vals(i) = func(iops%vsf_angles(i))
    end do

    call iops%calc_vsf_on_grid()
  end subroutine vsf_from_function

  subroutine calc_vsf_on_grid(iops)
    class(optical_properties) iops
    type(space_angle_grid) grid
    double precision th, ph, thp, php
    integer p, pp
    integer nomega
    double precision angle_diff
    double precision vsf_val
    double precision norm

    grid = iops%grid
    nomega = grid%angles%nomega

    ! Calculate cos VSF
    iops%vsf_cos = cos(iops%vsf_angles)

    ! Normalize cos VSF to 1/(2pi) on [-1, 1]
    iops%vsf_scat_coef = abs(trap_rule_uneven(iops%vsf_cos, iops%vsf_vals, iops%num_vsf))
    iops%vsf_vals(:) = iops%vsf_vals(:) / (2*pi * iops%vsf_scat_coef)

    ! write(*,*) 'norm = ', iops%vsf_scat_coef
    ! write(*,*) 'now: ', trap_rule_uneven(iops%vsf_cos, iops%vsf_vals, iops%num_vsf)
    ! write(*,*) 'cos: ', iops%vsf_cos
    ! write(*,*) 'vals: ', iops%vsf_vals

    do p=1, nomega
       th = grid%angles%theta_p(p)
       ph = grid%angles%phi_p(p)
       do  pp=1, nomega
          thp = grid%angles%theta_p(pp)
          php = grid%angles%phi_p(pp)
          ! TODO: Might be better to calculate average scattering
          ! from angular cell rather than only using center
          iops%vsf(p, pp) = iops%eval_vsf(angle_diff_3d(th,ph,thp,php))
       end do

       ! Normalize each row of VSF (midpoint rule)
       norm = sum(iops%vsf(p,:) * grid%angles%area_p(:))
       iops%vsf(p,:) = iops%vsf(p,:) / norm

       ! % / meter light scattered from cell pp into direction p.
       ! TODO: Could integrate VSF instead of just using value at center
       iops%vsf_integral(p, :) = iops%vsf(p, :) * grid%angles%area_p(:)
    end do

    ! Normalize VSF on unit sphere w.r.t. north pole
    !iops%vsf_scat_coef = sum(iops%vsf(1,:) * iops%grid%angles%area_p)
    !iops%vsf = iops%vsf / iops%vsf_scat_coef
    !iops%vsf_integral = iops%vsf_integral / iops%vsf_scat_coef
  end subroutine calc_vsf_on_grid

  subroutine iop_deinit(iops)
    class(optical_properties) iops
    deallocate(iops%vsf_angles)
    deallocate(iops%vsf_vals)
    deallocate(iops%vsf)
    deallocate(iops%vsf_integral)
    deallocate(iops%abs_water)
    deallocate(iops%scat_water)
    deallocate(iops%abs_grid)
  end subroutine iop_deinit

end module kelp_context
