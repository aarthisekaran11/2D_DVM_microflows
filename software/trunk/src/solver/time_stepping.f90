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

  use Constants

  implicit none

  private

  double precision :: t_conv, t_coll
  double precision :: t_final, t_frac
  double precision :: deltat

  integer :: num_time_steps

  public :: set_times
  public :: get_t_final
  public :: calculate_t_conv
!!$  public :: calculate_deltat

  public :: set_deltat
  public :: get_deltat
  public :: set_num_time_steps
  public :: get_num_time_steps


contains

  subroutine set_times( t_final_in, t_frac_in )

    implicit none

    double precision, intent(in) :: t_final_in, t_frac_in

    t_final = t_final_in
    t_frac  = t_frac_in

    t_conv = bigdp

    return
  end subroutine set_times

  subroutine get_t_final( t_final_out )

    implicit none

    double precision :: t_final_out

    t_final_out = t_final

    return
  end subroutine get_t_final

  subroutine calculate_t_conv( dx, eta_min, eta_max )

    implicit none

    double precision, intent(in) :: dx
    double precision, intent(in) :: eta_min, eta_max

    if( dx / abs( eta_min ) .lt. t_conv ) t_conv = dx / abs( eta_min )
    if( dx / abs( eta_max ) .lt. t_conv ) t_conv = dx / abs( eta_max )

    write(*,*)"t convection"
    write(*,*)t_conv

    return
  end subroutine calculate_t_conv

!!$  subroutine calculate_deltat( properties, num_species, num_space )
!!$
!!$    use PhysicalProperties
!!$
!!$    implicit none
!!$
!!$    type(PropertiesType), dimension(:), intent(in) :: properties
!!$    integer, intent(in) :: num_species, num_space
!!$    
!!$    integer :: n, s
!!$
!!$    t_coll = bigdp
!!$
!!$    do n = 1, num_space
!!$       do s = 1, num_species
!!$
!!$          if( properties(n)%t_coll(s) .lt. t_coll ) t_coll = properties(n)%t_coll(s)
!!$
!!$       end do
!!$    end do
!!$    
!!$    write(*,*)"t collision"
!!$    write(*,*)t_coll

!!$    deltat = min( t_conv, t_coll )

!!$    write(*,*)"delta t"
!!$    write(*,*)deltat

!!$    return
!!$  end subroutine calculate_deltat

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
