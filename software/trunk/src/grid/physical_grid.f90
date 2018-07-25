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

  double precision :: delta_x, delta_y
  integer :: nx_space, ny_space

  public :: set_physical_grid
  public :: get_delta_x
  public :: get_nspace

contains

  subroutine set_physical_grid( nx_space_in, delta_x_in, ny_space_in, delta_y_in )

    implicit none

    double precision, intent(in) :: delta_x_in, delta_y_in
    integer, intent(in) :: nx_space_in, ny_space_in

    delta_x  = delta_x_in
    delta_y  = delta_y_in
    nx_space = nx_space_in
    ny_space = ny_space_in

    return
  end subroutine set_physical_grid

  subroutine get_nspace( nx_space_out, ny_space_out )

    implicit none

    integer :: nx_space_out, ny_space_out

    nx_space_out = nx_space
    ny_space_out = ny_space

    return
  end subroutine get_nspace

  subroutine get_delta_x( delta_x_out, delta_y_out )

    implicit none

    double precision :: delta_x_out, delta_y_out

    delta_x_out = delta_x
    delta_y_out = delta_y

    return
  end subroutine get_delta_x

end module PhysicalGrid
