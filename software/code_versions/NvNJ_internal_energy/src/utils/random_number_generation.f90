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
module RandomNumberGeneration

  implicit none

  private

  integer :: seed
  integer :: rand_method

  integer, parameter :: gfortran_rng = 0
  integer, parameter :: sfmt         = 1

  public :: set_seed
  public :: set_rand_method
  public :: initialize_rand
  public :: get_rand
  public :: destroy_rand

contains

  subroutine set_seed( seed_in )

    implicit none

    integer, intent(in) :: seed_in

    seed = seed_in

    return
  end subroutine set_seed
  
  subroutine set_rand_method( method_in )

    implicit none

    integer, intent(in) :: method_in

    rand_method = method_in

    return
  end subroutine set_rand_method

  subroutine initialize_rand()
    
    implicit none

    select case( rand_method )
    case( gfortran_rng )
       if( seed .eq. 0 ) call random_seed()
       if( seed .ne. 0 ) then
          write(*,*) "Error: Invalid value of seed: ", seed
          stop
       end if

    case( sfmt )

    case default
       write(*,*) "Error: Invalid value for rand_method: ", rand_method
       stop

    end select

    return
  end subroutine initialize_rand

  function get_rand() result( rand_value )

    implicit none

    double precision :: rand_value

    select case( rand_method )
    case( gfortran_rng )
       call random_number( rand_value )

    case( sfmt )

    case default
       write(*,*) "Error: Invalid value for rand_method: ", rand_method
       stop

    end select

    return
  end function get_rand

  subroutine destroy_rand()

    implicit none

    select case( rand_method )
    case( gfortran_rng )

    case( sfmt )

    case default
       write(*,*) "Error: Invalid value for rand_method: ", rand_method
       stop

    end select

    return
  end subroutine destroy_rand

end module RandomNumberGeneration
