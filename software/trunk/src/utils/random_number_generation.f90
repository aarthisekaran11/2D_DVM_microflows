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

  use iso_c_binding

  implicit none

  private

  integer :: seed
  integer :: rand_method

  integer, parameter :: gfortran_rng = 0
  integer, parameter :: sfmt         = 1
  integer, parameter :: gsl_rng      = 2

  ! Pointer to C-struct that houses GSL RNG type and state info
  type(c_ptr) :: gsl_rng_struct

  ! Pointer to C-struct that defines state struct
  type(c_ptr) :: gsl_rng_type

  interface

     function gsl_rng_env_setup() bind(c)
       import :: c_ptr
       implicit none
       type(c_ptr) :: gsl_rng_env_setup
     end function gsl_rng_env_setup

     subroutine get_gsl_rng_default(gsl_rng_type) bind(c)
       import :: c_ptr
       implicit none
       type(c_ptr), value :: gsl_rng_type
     end subroutine get_gsl_rng_default

     function gsl_rng_alloc(gsl_rng_type) bind(c)
       import :: c_ptr
       implicit none
       type(c_ptr), value :: gsl_rng_type
       type(c_ptr) :: gsl_rng_alloc
     end function gsl_rng_alloc

     function gsl_rng_uniform(gsl_rng_struct) bind(c)
       import :: c_ptr, c_double
       implicit none
       type(c_ptr), value :: gsl_rng_struct
       real(c_double) :: gsl_rng_uniform
     end function gsl_rng_uniform

     subroutine gsl_rng_free(gsl_rng_struct) bind(c)
       import :: c_ptr
       implicit none
       type(c_ptr), value :: gsl_rng_struct
     end subroutine gsl_rng_free

  end interface

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

    case( gsl_rng )
       gsl_rng_type = gsl_rng_env_setup()
       call get_gsl_rng_default( gsl_rng_type )
       gsl_rng_struct = gsl_rng_alloc( gsl_rng_type )

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

    case( gsl_rng )
       rand_value = gsl_rng_uniform(gsl_rng_struct)

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

    case( gsl_rng )
       call gsl_rng_free(gsl_rng_struct)

    case( sfmt )

    case default
       write(*,*) "Error: Invalid value for rand_method: ", rand_method
       stop

    end select

    return
  end subroutine destroy_rand

end module RandomNumberGeneration
