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
module ErrorCheck

  implicit none

contains

  subroutine allocate_error_check( error, string_in )
    
    implicit none

    integer, intent(in) :: error
    character(len=*), intent(in) :: string_in

    if( error .ne. 0 )then
       write(*,*) "Could not allocate variable ", string_in,". Stopping."
       stop
    endif

    return
  end subroutine allocate_error_check

  
  subroutine deallocate_error_check( error, string_in )

    implicit none

    integer, intent(in) :: error
    character(len=*), intent(in) :: string_in

    if( error .ne. 0 )then
       write(*,*) "Could not deallocate variable ", string_in,". Stopping."
       stop
    endif

    return
  end subroutine deallocate_error_check

end module ErrorCheck
