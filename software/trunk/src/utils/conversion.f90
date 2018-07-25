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
module Conversion

  implicit none

  private

  public :: int2logical
  public :: int2char

contains

  function int2logical(int_in) result(val)

    implicit none

    integer, intent(in) :: int_in
    logical :: val

    select case( int_in )
    case(0)
       val = .false.

    case(1)
       val = .true.

    case default
       write(*,"(a,i2,a)") "Error: Cannot convert integer value ", int_in, " to logical."
       stop

    end select

    return
  end function int2logical

  function int2char(int_in) result(val)

    implicit none

    integer, intent(in) :: int_in
    character(len=32) :: val

    write(val,"(I31)") int_in

    val = trim(adjustl(val))

    return
  end function int2char

end module Conversion
