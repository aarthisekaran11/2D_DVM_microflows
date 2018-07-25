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
program input_read_vel_grid_test

  use VelocityGrid
  use ReadInput
  use Constants
  use ErrorCheck

  implicit none

  type(VelocityGridType), allocatable, dimension(:) :: vel_grid

  integer, parameter :: grid_ref = 1

  integer, parameter :: num_species = 2

  integer, parameter :: x_vel_grid_type = 0
  integer, parameter :: y_vel_grid_type = 0
  integer, parameter :: z_vel_grid_type = 0

  character(len=128), parameter :: filename = "./input_files/coll_test_grid.dat"

  integer, parameter :: num_points_x = 3
  integer, parameter :: num_points_y = 3
  integer, parameter :: num_points_z = 3

  double precision, dimension(2), parameter :: xv_min = (/-0.4d0, -1.2d0/)
  double precision, dimension(2), parameter :: xv_max = (/ 0.4d0,  1.2d0/)
  double precision, dimension(2), parameter :: yv_min = (/-0.4d0, -1.2d0/)
  double precision, dimension(2), parameter :: yv_max = (/ 0.4d0,  1.2d0/)
  double precision, dimension(2), parameter :: zv_min = (/-0.4d0, -1.2d0/)
  double precision, dimension(2), parameter :: zv_max = (/ 0.4d0,  1.2d0/)

  integer :: n

  integer :: i_error
  double precision :: dp_error

  integer :: status
  integer :: return_flag

  return_flag = 0

  allocate( vel_grid(1:num_species), STAT=status )
  call allocate_error_check( status, "vel_grid" )
  
  call initialize_grid_input_arrays( num_species )

  do n = 1, num_species
     call set_vel_grid_type( x_vel_grid_type, y_vel_grid_type, z_vel_grid_type, n )
  end do

  call read_vel_grid( filename, num_species )

  do n = 1, num_species
     call create_velocity_grid( vel_grid(n), n, grid_ref )
  end do

  call destroy_file_read_arrays( num_species )
  call destroy_grid_input_arrays()

  do n = 1, num_species

     i_error = vel_grid(n)%num_points_x - num_points_x
     if( i_error .ne. 0 )then
        write(*,*) "Incorrect number of x points for species ", n
        write(*,*) "Error is: ", i_error
        return_flag = 1
     end if

     i_error = vel_grid(n)%num_points_y - num_points_y
     if( i_error .ne. 0 )then
        write(*,*) "Incorrect number of y points for species ", n
        write(*,*) "Error is: ", i_error
        return_flag = 1
     end if

     i_error = vel_grid(n)%num_points_z - num_points_z
     if( i_error .ne. 0 )then
        write(*,*) "Incorrect number of z points for species ", n
        write(*,*) "Error is: ", i_error
        return_flag = 1
     end if

     dp_error = vel_grid(n)%xv_min - xv_min(n)
     if( abs(dp_error) .gt. double_tol )then
        write(*,*) "Incorrect value of xv_min for species ", n
        write(*,*) "Error is: ", dp_error
        return_flag = 1
     end if

     dp_error = vel_grid(n)%xv_max - xv_max(n)
     if( abs(dp_error) .gt. double_tol )then
        write(*,*) "Incorrect value of xv_max for species ", n
        write(*,*) "Error is: ", dp_error
        return_flag = 1
     end if

     dp_error = vel_grid(n)%yv_min - yv_min(n)
     if( abs(dp_error) .gt. double_tol )then
        write(*,*) "Incorrect value of yv_min for species ", n
        write(*,*) "Error is: ", dp_error
        return_flag = 1
     end if

     dp_error = vel_grid(n)%yv_max - yv_max(n)
     if( abs(dp_error) .gt. double_tol )then
        write(*,*) "Incorrect value of yv_max for species ", n
        write(*,*) "Error is: ", dp_error
        return_flag = 1
     end if

     dp_error = vel_grid(n)%zv_min - zv_min(n)
     if( abs(dp_error) .gt. double_tol )then
        write(*,*) "Incorrect value of zv_min for species ", n
        write(*,*) "Error is: ", dp_error
        return_flag = 1
     end if

     dp_error = vel_grid(n)%zv_max - zv_max(n)
     if( abs(dp_error) .gt. double_tol )then
        write(*,*) "Incorrect value of zv_max for species ", n
        write(*,*) "Error is: ", dp_error
        return_flag = 1
     end if

  end do

  do n = 1, num_species
     call destroy_velocity_grid( vel_grid(n) )
  end do

  deallocate( vel_grid, STAT=status )
  call deallocate_error_check( status, "vel_grid" )

  call exit(return_flag)

end program input_read_vel_grid_test
