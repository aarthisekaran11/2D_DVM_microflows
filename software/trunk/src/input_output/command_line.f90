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
module CommandLineProcessor

  implicit none

  private

  integer, parameter :: max_num_command_line_args = 2
  integer, parameter :: num_optional_args = 1
  integer, parameter :: filename_arg = 1
  integer, parameter :: vel_grid_arg = 2

  character(len=128) :: input_filename, vel_grid_filename

  integer :: num_command_line_args

  public :: process_command_line_args
  public :: get_input_filename

contains

  subroutine process_command_line_args()

    implicit none

    integer :: length, status

    ! Intrinsic function to F2003.
    num_command_line_args = command_argument_count()

    ! Check to make sure we have the correct number command line args
    if( (num_command_line_args .gt. max_num_command_line_args) .or. &
        (num_command_line_args .lt. max_num_command_line_args - num_optional_args) )then
       write(*,*) "Error: Number of command line arguments incorrect."
       write(*,*) "Found ", num_command_line_args
       write(*,*) "Maximum ", max_num_command_line_args
       write(*,*) "Minimum ", max_num_command_line_args-num_optional_args
       write(*,*) "Format should be: dvm_exe <input file>"
       stop
    endif

    call get_command_argument( filename_arg, input_filename, length, status )
    if( status .gt. 0 )then
       write(*,*) "Error: get_command_argument failed."
       stop
    endif
    if( status .eq. -1 )then
       write(*,*) "Warning: Variable input_file_name was trucated during read."
    endif

    if( num_command_line_args .gt. 1 )then
       call get_command_argument( vel_grid_arg, vel_grid_filename, length, status )
       if( status .gt. 0 )then
          write(*,*) "Error: get_command_argument failed."
          stop
       endif
       if( status .eq. -1 )then
          write(*,*) "Warning: Variable input_file_name was trucated during read."
       endif
    end if

    return
  end subroutine process_command_line_args

  subroutine get_input_filename( input_filename_out, vel_grid_filename_out, num_command_line_args_out )

    implicit none

    character(len=128) :: input_filename_out, vel_grid_filename_out
    integer :: num_command_line_args_out

    input_filename_out        = input_filename
    vel_grid_filename_out     = vel_grid_filename
    num_command_line_args_out = num_command_line_args

    return
  end subroutine get_input_filename

end module CommandLineProcessor
