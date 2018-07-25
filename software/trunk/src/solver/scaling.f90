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

  integer, parameter :: clarke = 1
  integer, parameter :: bird   = 2
  integer, parameter :: vki    = 3

  integer :: csection_model, csection_coeff

  double precision :: Kn, Lc

  double precision :: c_ref, sigma_ref
  double precision :: L_scale, t_scale

  integer, parameter :: ref_unit = 50
  character(len=20), parameter :: ref_filename = 'reference_values.txt'

  public :: set_cross_section
  public :: get_csection_model
  public :: set_Lc
  public :: set_knudsen_number
  public :: get_knudsen_number
  public :: set_scaling_parameters
  public :: get_scaling_parameters
  public :: calculate_cross_section

contains

  subroutine set_cross_section( csection_model_in, csection_coeff_in )

    implicit none

    integer, intent(in) :: csection_model_in, csection_coeff_in

    csection_model = csection_model_in
    csection_coeff = csection_coeff_in

    return
  end subroutine set_cross_section

  subroutine get_csection_model( csection_model_out )

    implicit none

    integer :: csection_model_out

    csection_model_out = csection_model

    return
  end subroutine get_csection_model

  subroutine set_Lc( Lc_in )

    implicit none

    double precision, intent(in) :: Lc_in 

!!$    L_scale = Lc_in

    return
  end subroutine set_Lc

  subroutine set_knudsen_number( Kn_in )

    implicit none

    double precision, intent(in) :: Kn_in

    ! TODO: instead of setting Kn, why don't I set L and 
    ! TODO: calculate Kn, that seems more natural.

    Kn = Kn_in

    return
  end subroutine set_knudsen_number

  subroutine get_knudsen_number( Kn_out )

    implicit none
    
    double precision :: Kn_out

    Kn_out = Kn

    return
  end subroutine get_knudsen_number

  subroutine set_scaling_parameters( n_ref, temp_ref, m_ref, d_ref )

    implicit none

    double precision, intent(in) :: n_ref, temp_ref, m_ref, d_ref

    !TODO: option 1:
!!$    if( L_scale .lt. zero )then
!!$       Kn = on
!!$    else
!!$       Kn = one / ( L_scale * n_ref * sigma_ref )
!!$    end if
    !TODO: option 2:

    ! Reference speed
    c_ref = sqrt( two * kb * temp_ref / m_ref ) ! [m/s]

    ! Reference Cross Section
    sigma_ref = pi * d_ref * d_ref ! [m^2]

    ! Length Scale
    L_scale = one / ( Kn * n_ref * sigma_ref ) ! [m]

    ! Time Scale
    t_scale = L_scale / c_ref ! [s]

    call write_reference_values_to_file( n_ref, temp_ref, m_ref, &
         d_ref, c_ref, sigma_ref, Kn )

    return
  end subroutine set_scaling_parameters

  subroutine get_scaling_parameters( c_ref_out, sigma_ref_out, t_ref_out )
    
    implicit none

    double precision :: c_ref_out, sigma_ref_out, t_ref_out

    c_ref_out = c_ref
    sigma_ref_out = sigma_ref
    t_ref_out = t_scale

    return
  end subroutine get_scaling_parameters
  
  subroutine calculate_cross_section( cross_section, m_red, omega, d, viscosity, temp_visc )

    implicit none
    
    double precision, dimension(2), intent(in) :: d, viscosity, temp_visc
    double precision, intent(in) :: m_red, omega

    double precision :: cross_section

    double precision :: gamma_fnct, sigma_vhs

    sigma_vhs = 0.25d0 * ( d(1) + d(2) ) * ( d(1) + d(2) )

    select case( csection_coeff )
    case( clarke )
       cross_section = sigma_vhs

    case( bird )
       gamma_fnct = gamma_calculator( 2.5d0 - omega )

       cross_section = sigma_vhs * m_red**( one_half - omega ) / gamma_fnct

    case( vki )
       ! Only applicable to single species gases
       gamma_fnct = gamma_calculator( 4.5d0 - omega )

       ! 2*pi*15/16/sqrt(pi) ~= 3.32335097045
       cross_section = ( 3.32335097045d0 * m_red**( one - omega ) * temp_visc(1)**omega ) &
            / ( gamma_fnct * viscosity(1) )

    case default
       write(*,*) "Error: Invalid value of csection_coeff: ", &
            csection_coeff
       stop

    end select

    return
  end subroutine calculate_cross_section

  subroutine write_reference_values_to_file( dens, temp, mass, diam, c, sigma, Kn )

    implicit none

    double precision, intent(in) :: dens, temp, mass, diam, c, sigma, Kn
    double precision :: L, press, hf, e
    double precision :: t

    ! Calculate unknown ref values
    L      = one / ( Kn * dens * sigma )
    press  = two * kb * dens * temp
    hf     = mass * c * c * c
    e      = mass * c * c
    t      = L / c

    open( unit=ref_unit, file=ref_filename, status='unknown' )

    write( ref_unit, '(a)' ) 'Reference values:'
    write( ref_unit, '(a7,es16.8)' )  'Mass = ', mass
    write( ref_unit, '(a11,es16.8)' ) 'Diameter = ', diam
    write( ref_unit, '(a16,es16.8)' ) 'Cross-section = ', sigma
    write( ref_unit, '(a10,es16.8)' ) 'Density = ', dens
    write( ref_unit, '(a14,es16.8)' ) 'Temperature = ', temp
    write( ref_unit, '(a8,es16.8)' )  'Speed = ', c
    write( ref_unit, '(a11,es16.8)' ) 'Pressure = ', press
    write( ref_unit, '(a12,es16.8)' ) 'Heat Flux = ', hf
    write( ref_unit, '(a9,es16.8)' )  'Energy = ', e
    write( ref_unit, '(a22,es16.8)' ) 'Flow Knudsen number = ', Kn
    write( ref_unit, '(a30,es16.8)' ) 'Characteristic length scale = ', L
    write( ref_unit, '(a13,es16.8)' ) 'Time scale = ', t

    close( ref_unit )

    return
  end subroutine write_reference_values_to_file

end module Scaling
