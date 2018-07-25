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
program interp_x_test

  implicit none

  type(VelocityGridType) :: vel_grid

  integer, parameter :: grid_ref = 1
  integer, parameter :: num_species = 1

  integer, parameter :: x_vel_grid_type = 0
  integer, parameter :: y_vel_grid_type = 0
  integer, parameter :: z_vel_grid_type = 0

  character(len=128), parameter :: filename = "./input_files/interp_test_grid.dat"

  type(MappingResultType) :: mapping

  double precision, parameter :: u = 0.1d0
  double precision, parameter :: v = 0.4d0
  double precision, parameter :: w = 0.2d0

  double precision, parameter :: f_off = 1.0d0

  integer, parameter :: o_i =  1
  integer, parameter :: o_j =  2
  integer, parameter :: o_k =  2
  integer, parameter :: a   =  1
  integer, parameter :: b   = -1
  integer, parameter :: c   = -1
  integer, parameter :: ext_map_flag = 1

  double precision, parameter :: f_o  = zero
  double precision, parameter :: f_ix = zero
  double precision, parameter :: f_iy = zero
  double precision, parameter :: f_iz = zero
  double precision, parameter :: f_ex = -0.275d0
  double precision, parameter :: f_ey = zero
  double precision, parameter :: f_ez = zero

  double precision :: dp_error
  integer :: i_error

  integer :: return_flag
  
  return_flag = 0

  call initialize_grid_input_arrays( num_species )

  do n = 1, num_species
     call set_vel_grid_type( x_vel_grid_type, y_vel_grid_type, z_vel_grid_type, n )
  end do

  call read_vel_grid( filename, num_species )

  do n = 1, num_species
     call create_velocity_grid( vel_grid, n, grid_ref )
  end do

  call destroy_file_read_arrays( num_species )
  call destroy_grid_input_arrays()

  call perform_remapping( mapping, u, v, w, f_off, vel_grid )
  
  i_error = mapping%stencil(1) - o_i
  if( i_error .ne. 0 )then
     write(*,*) "o_i failed test. Error = ", i_error
     return_flag = 1
  end if

  i_error = mapping%stencil(2) - o_j
  if( i_error .ne. 0 )then
     write(*,*) "o_j failed test. Error = ", i_error
     return_flag = 1
  end if

  i_error = mapping%stencil(3) - o_k
  if( i_error .ne. 0 )then
     write(*,*) "o_k failed test. Error = ", i_error
     return_flag = 1
  end if

  i_error = mapping%stencil(4) - a
  if( i_error .ne. 0 )then
     write(*,*) "a failed test. Error = ", i_error
     return_flag = 1
  end if

  i_error = mapping%stencil(5) - b
  if( i_error .ne. 0 )then
     write(*,*) "b failed test. Error = ", i_error
     return_flag = 1
  end if

  i_error = mapping%stencil(6) - c
  if( i_error .ne. 0 )then
     write(*,*) "c failed test. Error = ", i_error
     return_flag = 1
  end if

  i_error = mapping%stencil(7) - ext_map_flag
  if( i_error .ne. 0 )then
     write(*,*) "ext_map_flag failed test. Error = ", i_error
     return_flag = 1
  end if

  dp_error = mapping%stencil(1) - f_o
  if( abs(dp_error) .gt. double_tol )then
     write(*,*) "f_o failed test. Error = ", dp_error
     return_flag = 1
  end if

  dp_error = mapping%stencil(2) - f_ix
  if( abs(dp_error) .gt. double_tol )then
     write(*,*) "f_ix failed test. Error = ", dp_error
     return_flag = 1
  end if

  dp_error = mapping%stencil(3) - f_iy
  if( abs(dp_error) .gt. double_tol )then
     write(*,*) "f_iy failed test. Error = ", dp_error
     return_flag = 1
  end if

  dp_error = mapping%stencil(4) - f_iz
  if( abs(dp_error) .gt. double_tol )then
     write(*,*) "f_iz failed test. Error = ", dp_error
     return_flag = 1
  end if

  dp_error = mapping%stencil(5) - f_ex
  if( abs(dp_error) .gt. double_tol )then
     write(*,*) "f_ex failed test. Error = ", dp_error
     return_flag = 1
  end if

  dp_error = mapping%stencil(6) - f_ey
  if( abs(dp_error) .gt. double_tol )then
     write(*,*) "f_ey failed test. Error = ", dp_error
     return_flag = 1
  end if

  dp_error = mapping%stencil(7) - f_ez
  if( abs(dp_error) .gt. double_tol )then
     write(*,*) "f_ez failed test. Error = ", dp_error
     return_flag = 1
  end if

  call exit(return_flag)

end program interp_x_test
