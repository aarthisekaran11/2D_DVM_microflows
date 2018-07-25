!-----------------------------------------------------------------------bl-
!--------------------------------------------------------------------------
! 
! DVM - a discrete velocity method for solving the Boltzmann equation
!
! Copyright (C) 2010,2011 The PECOS Development Team
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the Version 2 GNU General
! Public License as published by the Free Software Foundation.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this library; if not, write to the Free Software
! Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
! 02110-1301 USA
!
!-----------------------------------------------------------------------el-
! $Id: $
!--------------------------------------------------------------------------
module PhysicalGrid

  use ErrorCheck

  implicit none

  private

  double precision :: delta_x
  integer :: nspace, num_zones

  integer, allocatable, dimension(:) :: spatial_grid_reference

  public :: set_physical_grid
  public :: get_delta_x
  public :: get_nspace
  public :: get_num_zones
  public :: initialize_spatial_grid_reference
  public :: set_spatial_reference
  public :: get_spatial_reference

contains

  subroutine set_physical_grid( nspace_in, delta_x_in, num_zones_in )

    implicit none

    double precision, intent(in) :: delta_x_in
    integer, intent(in) :: nspace_in, num_zones_in

    nspace    = nspace_in
    delta_x   = delta_x_in
    num_zones = num_zones_in

    return
  end subroutine set_physical_grid

  subroutine get_nspace( nspace_out )

    implicit none

    integer :: nspace_out

    nspace_out = nspace

    return
  end subroutine get_nspace

  subroutine get_delta_x( delta_x_out )

    implicit none

    double precision :: delta_x_out

    delta_x_out = delta_x

    return
  end subroutine get_delta_x

  subroutine get_num_zones( num_zones_out )

    implicit none

    integer :: num_zones_out

    num_zones_out = num_zones

    return
  end subroutine get_num_zones

  subroutine initialize_spatial_grid_reference( )

    implicit none

    integer :: status

    allocate( spatial_grid_reference( 1:nspace ), STAT=status )
    call allocate_error_check( status, "spatial_grid_reference" )

    return
  end subroutine initialize_spatial_grid_reference

  subroutine set_spatial_reference( reference, start_values )

    implicit none

    integer, dimension(:), intent(in) :: start_values
    integer, dimension(:), intent(in) :: reference

    integer :: start, end_point
    integer :: i, j

    do i = 1, num_zones
       if( i .eq. num_zones )then
          end_point = nspace
       else
          end_point = start_values( i + 1 ) - 1
       end if
       start = start_values( i )

       do j = start, end_point
          spatial_grid_reference( j ) = reference( i )
       end do

    end do

    return
  end subroutine set_spatial_reference

  subroutine get_spatial_reference( space, reference )

    implicit none

    integer, intent(in) :: space
    integer :: reference

    reference = spatial_grid_reference( space )

    return
  end subroutine get_spatial_reference

end module PhysicalGrid
