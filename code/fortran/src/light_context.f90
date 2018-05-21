module light_context
#include "lisf.h"
  use sag
  use rte_sparse_matrices
  !use hdf5
  implicit none

  type light_state
     double precision, dimension(:,:,:), allocatable :: irradiance
     double precision, dimension(:,:,:,:), allocatable :: radiance
     type(space_angle_grid) :: grid
     type(rte_mat) :: mat
   contains
     procedure :: init => light_init
     procedure :: init_grid => light_init_grid
     procedure :: calculate_radiance
     procedure :: calculate_irradiance
     procedure :: deinit => light_deinit
     !procedure :: to_hdf => light_to_hdf
  end type light_state

contains

  ! Init for use with mat
  subroutine light_init(light, mat)
    class(light_state) light
    type(rte_mat) mat
    integer nx, ny, nz, nomega

    light%mat = mat
    light%grid = mat%grid

    nx = light%grid%x%num
    ny = light%grid%y%num
    nz = light%grid%z%num
    nomega = light%grid%angles%nomega

    allocate(light%irradiance(nx, ny, nz))
    allocate(light%radiance(nx, ny, nz, nomega))
  end subroutine light_init

  ! Init for use without mat
  subroutine light_init_grid(light, grid)
    class(light_state) light
    type(space_angle_grid) grid
    integer nx, ny, nz, nomega

    light%grid = grid

    nx = light%grid%x%num
    ny = light%grid%y%num
    nz = light%grid%z%num
    nomega = light%grid%angles%nomega

    allocate(light%irradiance(nx, ny, nz))
    allocate(light%radiance(nx, ny, nz, nomega))
  end subroutine light_init_grid

  subroutine calculate_radiance(light)
    class(light_state) light
    integer i, j, k, p
    integer nx, ny, nz, nomega
    integer index

    nx = light%grid%x%num
    ny = light%grid%y%num
    nz = light%grid%z%num
    nomega = light%grid%angles%nomega

    ! call lis_vector_get_size(light%mat%x, ln, gn)

    ! write(*,*) 'ln =', ln
    ! write(*,*) 'gn =', gn

    index = 1

    ! Set initial guess from provided radiance
    ! Traverse solution vector in order
    ! so as to avoid calculating index
    do k=1, nz
       do i=1, nx
           do j=1, ny
             do p=1, nomega
                call lis_vector_set_value(LIS_INS_VALUE, index, &
                     light%radiance(i,j,k,p), light%mat%x, light%mat%ierr)
                if(light%mat%ierr .ne. 0) then
                   write(*,*) 'IG ERROR:', light%mat%ierr
                end if

                index = index + 1
             end do
          end do
       end do
    end do

    !call light%mat%initial_guess()

    ! Solve (LIS)
    call light%mat%solve()

    index = 1

    ! Extract solution
    do k=1, nz
       do i=1, nx
          do j=1, ny
             do p=1, nomega
                call lis_vector_get_value(light%mat%x, index, &
                     light%radiance(i,j,k,p), light%mat%ierr)
                if(light%mat%ierr .ne. 0) then
                   write(*,*) 'EXTRACT ERROR:', light%mat%ierr
                end if
                index = index + 1
             end do
          end do
       end do
    end do
  end subroutine calculate_radiance

  subroutine calculate_irradiance(light)
    class(light_state) light
    integer i, j, k
    integer nx, ny, nz
    double precision, dimension(light%grid%angles%nomega) :: tmp_rad

    nx = light%grid%x%num
    ny = light%grid%y%num
    nz = light%grid%z%num

    do i=1, nx
       do j=1, ny
          do k=1, nz
             ! Use temporary array to avoid creating one
             ! implicitly at every spatial grid point
             tmp_rad = light%radiance(i,j,k,:)
             light%irradiance(i,j,k) = &
                  light%grid%angles%integrate_points(tmp_rad)
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
