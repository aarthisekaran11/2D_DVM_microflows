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
module Scaling

  use Constants
  use MathUtilities

  implicit none

  private

  integer :: csection_model, csection_coeff

  integer, parameter :: clarke = 1
  integer, parameter :: bird   = 2

  double precision :: c_ref, sigma_ref
  double precision :: L_scale, t_scale

  public :: set_cross_section
  public :: set_scaling_parameters
  public :: get_csection_model
  public :: calculate_cross_section

contains

  subroutine set_cross_section( csection_model_in, csection_coeff_in )

    implicit none

    integer, intent(in) :: csection_model_in, csection_coeff_in

    csection_model = csection_model_in
    csection_coeff = csection_coeff_in

    return
  end subroutine set_cross_section

  subroutine set_scaling_parameters( n_ref, temp_ref, m_ref, d_ref )

    implicit none

    double precision, intent(in) :: n_ref, temp_ref, m_ref, d_ref

    ! Reference speed
    c_ref = sqrt( 2.0d0*kb*temp_ref/m_ref ) ! [m/s]

    ! Reference Cross Section
    sigma_ref = pi*d_ref*d_ref

    ! Length Scale
    L_scale = one/(n_ref*sigma_ref)

    ! Time Scale
    t_scale = L_scale/c_ref

    return
  end subroutine set_scaling_parameters

  subroutine get_csection_model( csection_model_out )

    implicit none

    integer :: csection_model_out

    csection_model_out = csection_model

    return
  end subroutine get_csection_model
  
  subroutine calculate_cross_section( cross_section, m_red, omega, d )

    implicit none
    
    double precision, dimension(2), intent(in) :: d
    double precision, intent(in) :: m_red, omega

    double precision :: cross_section

    double precision :: gamma_fnct, sigma_vhs

    sigma_vhs = 0.25d0*( d(1) + d(2) )*( d(1) + d(2) )

    select case( csection_coeff )
    case( clarke )
       cross_section = sigma_vhs

    case( bird )
       gamma_fnct = gamma_calculator( 2.5d0 - omega )

       cross_section = sigma_vhs*m_red**( one_half - omega )/gamma_fnct

    case default
       write(*,*) "Error: Invalid value of csection_coeff: ", &
            csection_coeff
       stop

    end select
    
    return
  end subroutine calculate_cross_section

end module Scaling
