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
module TimeStepping

  implicit none

  private

  double precision :: deltat

  integer :: num_time_steps

  public :: set_deltat
  public :: get_deltat
  public :: set_num_time_steps
  public :: get_num_time_steps

contains

  subroutine set_deltat( deltat_in )

    implicit none

    double precision, intent(in) :: deltat_in

    deltat = deltat_in

    return
  end subroutine set_deltat

  subroutine get_deltat( deltat_out )

    implicit none

    double precision, intent(out) :: deltat_out

    deltat_out = deltat

    return
  end subroutine get_deltat

  subroutine set_num_time_steps( num_time_steps_in )

    implicit none

    integer, intent(in) :: num_time_steps_in

    num_time_steps = num_time_steps_in

    return
  end subroutine set_num_time_steps

  subroutine get_num_time_steps( num_time_steps_out )

    implicit none

    integer, intent(out) :: num_time_steps_out

    num_time_steps_out = num_time_steps

    return
  end subroutine get_num_time_steps

end module TimeStepping
