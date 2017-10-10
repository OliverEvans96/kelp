module light_context
  use sag
  use rte_sparse_matrices
  !use hdf5
  implicit none

  type light_state
     double precision, dimension(:,:,:), allocatable :: irradiance
     double precision, dimension(:,:,:,:,:), allocatable :: radiance
     type(space_angle_grid) :: grid
     type(rte_mat) :: mat
   contains
     procedure :: init => light_init
     procedure :: calculate_radiance
     procedure :: calculate_irradiance
     procedure :: deinit => light_deinit
     !procedure :: to_hdf => light_to_hdf
  end type light_state

contains
  subroutine light_init(light, mat)
    class(light_state) light
    type(rte_mat) mat
    integer nx, ny, nz, ntheta, nphi

    light%mat = mat
    light%grid = mat%grid

    nx = light%grid%x%num
    ny = light%grid%y%num
    nz = light%grid%z%num
    ntheta = light%grid%theta%num
    nphi = light%grid%phi%num

    allocate(light%irradiance(nx, ny, nz))
    allocate(light%radiance(nx, ny, nz, ntheta, nphi))
  end subroutine light_init

  subroutine calculate_radiance(light)
    class(light_state) light
    integer i, j, k, l, m
    integer nx, ny, nz, ntheta, nphi
    integer index

    nx = light%grid%x%num
    ny = light%grid%y%num
    nz = light%grid%z%num
    ntheta = light%grid%theta%num
    nphi = light%grid%phi%num

    call light%mat%solve()

    index = 1

    ! Traverse solution vector in order
    ! so as to avoid calculating index
    do k=1, nz
       do j=1, ny
          do i=1, nx
             do m=1, nphi
                do l=1, ntheta
                   light%radiance(i,j,k,l,m) = light%mat%sol(index)
                   index = index + 1
                end do
             end do
          end do
       end do
    end do
  end subroutine calculate_radiance

  subroutine calculate_irradiance(light)
    class(light_state) light
    integer i, j, k
    integer nx, ny, nz

    do i=1, nx
       do j=1, ny
          do k=1, nz
             light%irradiance(i,j,k) = light%grid%integrate_angle_2d( &
                  light%radiance(i,j,k,:,:))
          end do
       end do
    end do

  end subroutine calculate_irradiance

!  subroutine light_to_hdf(light, radfile, irradfile)
!    class(light_state) light
!    character(len=*) radfile
!    character(len=*) irradfile
!
!    call hdf_write_radiance(radfile, light%radiance, light%grid)
!    call hdf_write_irradiance(irradfile, light%irradiance, light%grid)
!  end subroutine light_to_hdf

  subroutine light_deinit(light)
    class(light_state) light

    deallocate(light%irradiance)
    deallocate(light%radiance)
  end subroutine light_deinit
end module
