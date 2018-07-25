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
program coll_com_vel_test
  
  use CollisionUtilities
  use Constants

  implicit none

  double precision, dimension(2), parameter :: x = (/-0.4d0, 1.2d0/)
  double precision, dimension(2), parameter :: y = (/-0.4d0, 1.2d0/)
  double precision, dimension(2), parameter :: z = (/-0.4d0, 1.2d0/)

  double precision, dimension(2), parameter :: mass = (/1.0d0, 0.1d0/)

  double precision, parameter :: exact_com_x = -0.254545454545454d0
  double precision, parameter :: exact_com_y = -0.254545454545454d0
  double precision, parameter :: exact_com_z = -0.254545454545454d0

  double precision, dimension(3) :: com

  double precision :: dp_error

  integer :: n

  integer :: return_flag

  return_flag = 0

  call compute_center_of_mass_velocity( x, y, z, mass, com )

  dp_error = com(1) - exact_com_x
  if( abs(dp_error) .gt. double_tol )then
     write(*,*) "Incorrect center of mass x location. Error = ", dp_error
     return_flag = 1
  end if
  
  dp_error = com(2) - exact_com_y
  if( abs(dp_error) .gt. double_tol )then
     write(*,*) "Incorrect center of mass y location. Error = ", dp_error
     return_flag = 1
  end if

  dp_error = com(3) - exact_com_z
  if( abs(dp_error) .gt. double_tol )then
     write(*,*) "Incorrect center of mass z location. Error = ", dp_error
     return_flag = 1
  end if

  call exit(return_flag)
  
end program coll_com_vel_test
