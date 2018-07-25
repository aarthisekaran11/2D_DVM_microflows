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
module TransferProbabilities

  use Constants
  use SpeciesAndReferenceData

  implicit none

  private

  type ProbabilityType

     character(len=16) :: species_A, species_B

     double precision, dimension(:), allocatable :: rel_velocities
     double precision, dimension(:,:,:), allocatable :: rot_probability
     double precision, dimension(:,:,:), allocatable :: vib_probability

  end type ProbabilityType

  type(ProbabilityType), dimension(:), allocatable :: probabilities

  integer :: rot_method, vib_method
  integer, dimension(:), allocatable :: rot_equil, vib_equil

  integer, parameter :: file_read = 0

  integer, parameter :: rigid_rotor_LB    = 1
  integer, parameter :: nonrigid_rotor_LB = 2

  integer, parameter :: SHO_LB = 1
  integer, parameter :: AHO_LB = 2

  ! Subroutines
  public :: set_energy_methods
  public :: initialize_energy_exchange
  public :: compute_rot_level_change
  public :: compute_vib_level_change

contains

  subroutine set_energy_methods( rot_method_in, vib_method_in )

    implicit none

    integer, intent(in) :: rot_method_in, vib_method_in

    rot_method = rot_method_in
    vib_method = vib_method_in

    return
  end subroutine set_energy_methods

  subroutine initialize_energy_exchange( molecule )

    use DistFunc

    implicit none

    type(MoleculeType), dimension(:), intent(in) :: molecule
    integer :: length

    if( rot_method .eq. 0 .or. vib_method .eq. 0 )then
       length = num_species*num_species
       allocate( probabilities(1:length) )
    end if

    if( rot_method .eq. 0 ) call read_rot_energy_exchange_probabilities( molecule )

    if( vib_method .eq. 0 ) call read_vib_energy_exchange_probabilities( molecule )

    return
  end subroutine initialize_energy_exchange

  subroutine  read_rot_energy_exchange_probabilities( molecule )

    use DistFunc

    implicit none

    type(MoleculeType), dimension(:), intent(in) :: molecule

    character(len=16) :: species_A, species_B
    integer :: r_modes_A, r_modes_B

    integer :: g_length, r_levels

    integer :: n, m

    do n = 1, num_species
       do m = n, num_species
          species_A = molecule(n)%species_code
          species_B = molecule(m)%species_code

          r_modes_A = molecule(n)%rot_modes
          r_modes_B = molecule(m)%rot_modes

          if( r_modes_A .eq. 0 .and. r_modes_B .eq. 0 ) cycle

    end do

    return
  end subroutine read_rot_energy_exchange_probabilities

  subroutine compute_rot_level_change( Jp1, Jp2, d_E, J1, J2, g, omega, m_red, levels, &
       molecule1, molecule2, level_array1, level_array2 )

    implicit none

    type(MoleculeType), intent(in) :: molecule1, molecule2
    integer, intent(in) :: J1, J2
    integer, dimension(2), intent(in) :: levels
    double precision, intent(in) :: m_red, g, omega
    double precision, dimension(:), intent(in) :: level_array1, level_array2

    integer, dimension(2) :: modes, type
    double precision, dimension(2) :: mass, theta

    integer :: Jp1, Jp2
    double precision :: d_E

    double precision :: E_c, E_t

    mass(1)  = molecule1%mass
    theta(1) = molecule1%theta_r
    modes(1) = molecule1%rot_modes
    type(1)  = molecule1%molecule_type

    mass(2)  = molecule2%mass
    theta(2) = molecule2%theta_r
    modes(2) = molecule2%rot_modes
    type(2)  = molecule2%molecule_type

    select case( rot_method )
    case( file_read )


    case( rigid_rotor_LB )
       E_t = one_half*m_red*g*g
       E_c = E_t
       call rot_rigid_rotor_LB( Jp1, E_c, J1, omega, mass(1), theta(1), modes(1), &
            levels(1), level_array1, type(1) )
       call rot_rigid_rotor_LB( Jp2, E_c, J2, omega, mass(2), theta(2), modes(2), &
            levels(2), level_array2, type(2) )
       d_E = E_t - E_c

    case( nonrigid_rotor_LB )

    case default
       write(*,*) "Error: Invalid value of rot_method: ", rot_method
       stop

    end select

    return
  end subroutine compute_rot_level_change

  subroutine compute_vib_level_change( vp1, vp2, d_E, v1, v2, g, omega, m_red, levels, &
       molecule1, molecule2, level_array1, level_array2 )

    implicit none

    type(MoleculeType), intent(in) :: molecule1, molecule2
    integer, intent(in) :: v1, v2
    integer, dimension(2), intent(in) :: levels
    double precision, intent(in) :: m_red, g, omega
    double precision, dimension(:), intent(in) :: level_array1, level_array2

    integer, dimension(2) :: modes
    double precision, dimension(2) :: mass, theta, Y10, Y20

    integer :: vp1, vp2
    double precision :: d_E

    double precision :: E_c, E_t

    mass(1)  = molecule1%mass
    theta(1) = molecule1%theta_v
    modes(1) = molecule1%vib_modes
    Y10(1)   = molecule1%Y10
    Y20(1)   = molecule1%Y20

    mass(2)  = molecule2%mass
    theta(2) = molecule2%theta_v
    modes(2) = molecule2%vib_modes
    Y10(2)   = molecule2%Y10
    Y20(2)   = molecule2%Y20

    select case( vib_method )
    case( file_read )


    case( SHO_LB )
       E_t = one_half*m_red*g*g
       E_c = E_t
       call vib_SHO_LB( vp1, E_c, v1, omega, mass(1), theta(1), modes(1), levels(1), level_array1 )
       call vib_SHO_LB( vp2, E_c, v2, omega, mass(2), theta(2), modes(2), levels(2), level_array2 )
       d_E = E_t - E_c

    case( AHO_LB )
       E_t = one_half*m_red*g*g
       E_c = E_t
       call vib_AHO_LB( vp1, E_c, v1, omega, mass(1), Y10(1), Y20(1), modes(1), levels(1), level_array1 )
       call vib_AHO_LB( vp2, E_c, v2, omega, mass(2), Y10(2), Y20(2), modes(2), levels(2), level_array2 )
       d_E = E_t - E_c

    case default
       write(*,*) "Error: Invalid value of vib_method: ", vib_method
       stop

    end select

    return
  end subroutine compute_vib_level_change



end module TransferProbabilities
