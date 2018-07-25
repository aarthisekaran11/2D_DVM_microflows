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
module SpeciesAndReferenceData

  use Constants
  use ErrorCheck

  implicit none

  private

  integer, parameter :: constant       = 1
  integer, parameter :: temp_dependent = 2

  integer, parameter :: monatomic     = 0
  integer, parameter :: homonuclear   = 1
  integer, parameter :: heteronuclear = 2

  integer, parameter :: hard_sphere    = 1
  integer, parameter :: pseudo_maxwell = 2
  integer, parameter :: vhs            = 3
  integer, parameter :: vss            = 4

  type MoleculeType

     character(len=16) :: species_code

     integer :: rot_modes, vib_modes
     integer :: molecule_type
     
     double precision :: mass, diam, omega, alpha

     double precision :: theta_r, theta_v, theta_d

     integer :: rot_relax_form, vib_relax_form
     double precision :: cZr, cZv
     double precision :: Zr_inf, T_star
     double precision :: Zv_1, Zv_2
     double precision :: vib_cs, Av, Bv, coll_freq_cnst, p_Zvscale

     double precision :: even_spin, odd_spin

     double precision :: Y01, Y02, Y10, Y20

     ! coefficients for vibration temperature vs energy relation
     double precision, allocatable, dimension(:) :: c_vib
     double precision, allocatable, dimension(:) :: c_vib_low, c_vib_high
     double precision, allocatable, dimension(:) :: c_rot_low, c_rot_high
     double precision :: E_v_curve_fit_split
     double precision :: E_r_curve_fit_split

     ! coefficients for rotation temperature vs energy relation
     double precision, allocatable, dimension(:) :: c_rot
     
  end type MoleculeType

  character(len=16), allocatable, dimension(:) :: species_codes

  character(len=16) :: ref_species
  double precision :: temp_ref, dens_ref, diam_ref, mass_ref

  ! public data types
  public :: MoleculeType

  ! public subroutines
  public :: set_species_codes
  public :: set_reference_parameters
  public :: set_molecule_constants
  public :: get_reference_mol_data
  public :: get_species_rot_dof

contains

  subroutine get_species_rot_dof( species_code, rot_dof )

    implicit none

    double precision :: rot_dof
    
    character(len=16), intent(in) :: species_code

    select case( species_code )

       ! rotational degrees of freedom

    case( "He", "Ne", "Ar", "Xe", "N2_vib" )
       rot_dof = zero

    case( "H2", "N2", "O2", "N2_rot" )
       rot_dof = two

    case default
       write(*,*) "Error, invalid value of species_code. Given value is: ", species_code
       stop

    end select

    return
  end subroutine get_species_rot_dof

  subroutine set_species_codes( codes_in, num_species )

    implicit none

    character(len=16), dimension(:), intent(in) :: codes_in

    integer :: num_species, status

    allocate( species_codes(1:num_species), STAT=status )
    call allocate_error_check( status, "species_codes input" )

    species_codes = codes_in

    return
  end subroutine set_species_codes

  subroutine set_reference_parameters( temp_ref_in, dens_ref_in, ref_species_in )

    implicit none

    double precision, intent(in) :: temp_ref_in, dens_ref_in
    character(len=16), intent(in) :: ref_species_in
    
    temp_ref = temp_ref_in
    dens_ref = dens_ref_in
    ref_species = ref_species_in

    return
  end subroutine set_reference_parameters

  subroutine set_molecule_constants( molecule, num_species )
    
    use Scaling

    implicit none

    integer, intent(in) :: num_species

    type(MoleculeType), dimension(:) :: molecule

    type(MoleculeType) :: ref_molecule

    integer :: n

    ! Reference parameters
    call get_molecule_info( ref_molecule, ref_species )
    diam_ref = ref_molecule%diam
    mass_ref = ref_molecule%mass

    call set_scaling_parameters( dens_ref, temp_ref, mass_ref, diam_ref )

    ! Read molecular data for each species and scale by reference parameters
    do n = 1, num_species
       call get_molecule_info( molecule(n), species_codes(n) )

       ! Scaling and calculation of time independent constants
       molecule(n)%diam = molecule(n)%diam/diam_ref
       molecule(n)%mass = molecule(n)%mass/mass_ref

       if( molecule(n)%rot_modes .gt. 0 .or. molecule(n)%vib_modes .gt. 0 )then
          molecule(n)%theta_d = molecule(n)%theta_d/temp_ref
       end if

       if( molecule(n)%rot_modes .gt. 0 )then
          ! Relaxation rate coefficients
          molecule(n)%T_star  = molecule(n)%T_star/temp_ref

          ! Energy fraction coefficients
          molecule(n)%theta_r = molecule(n)%theta_r/temp_ref
          molecule(n)%Y01     = molecule(n)%Y01/(cm2J*two*kb*temp_ref)
          molecule(n)%Y02     = molecule(n)%Y02/(cm2J*two*kb*temp_ref)
       end if

       if( molecule(n)%vib_modes .gt. 0 )then
          ! Relaxation rate coefficients
          molecule(n)%Av      = &
               1.16d-3*molecule(n)%theta_v**(4.0d0/3.0d0)*(mass_ref/kg2amu)**(0.5d0)*temp_ref**(-one_third)
          molecule(n)%Bv      = -1.74d-5*molecule(n)%theta_v**(4.0d0/3.0d0)*(mass_ref/kg2amu)**(0.75d0)
          molecule(n)%vib_cs  = molecule(n)%vib_cs*angstrom*angstrom/diam_ref/diam_ref/pi

          ! Bird relaxation rate coefficients
          molecule(n)%Zv_1 = molecule(n)%Zv_1/( temp_ref**molecule(n)%omega )
          molecule(n)%Zv_2 = molecule(n)%Zv_2 * temp_ref**(-one_third)

          ! Energy fraction coefficients
          molecule(n)%theta_v = molecule(n)%theta_v/temp_ref
          molecule(n)%Y10     = molecule(n)%Y10/(cm2J*two*kb*temp_ref)
          molecule(n)%Y20     = molecule(n)%Y20/(cm2J*two*kb*temp_ref)
       end if

       ! Collision frequency constants - TODO: not correct
       molecule(n)%coll_freq_cnst = two*pi**(-0.5d0)*molecule(n)%diam*molecule(n)%diam*&
            (temp_ref/temp_diam)**(one_half-molecule(n)%omega)
       molecule(n)%p_Zvscale = Pa2atm*( pi*diam_ref*diam_ref*sqrt(two*kb*temp_ref/mass_ref) )/( kb*temp_ref )

    end do

    return
  end subroutine set_molecule_constants

  subroutine get_reference_mol_data( diam_ref_out, mass_ref_out, temp_ref_out, dens_ref_out )
    
    implicit none

    double precision :: diam_ref_out, mass_ref_out, temp_ref_out, dens_ref_out

    diam_ref_out = diam_ref
    mass_ref_out = mass_ref
    temp_ref_out = temp_ref
    dens_ref_out = dens_ref

    return
  end subroutine get_reference_mol_data

  subroutine get_molecule_info( molecule, species_code_in )

    use scaling

    implicit none

    character(len=16), intent(in) :: species_code_in

    type(MoleculeType) :: molecule

    integer :: csection_model

    molecule%species_code = species_code_in

    select case( species_code_in )
    case( "He" )
       molecule%diam  = 2.33e-10
       molecule%mass  = 6.65e-27
       molecule%omega = 0.66d0
       molecule%alpha = 1.31d0

       molecule%rot_modes     = 0
       molecule%vib_modes     = 0
       molecule%theta_d       = zero
       molecule%molecule_type = monatomic

    case( "Ne" )
       molecule%diam  = 2.77e-10
       molecule%mass  = 33.5e-27
       molecule%omega = 0.66d0
       molecule%alpha = 1.31d0

       molecule%rot_modes     = 0
       molecule%vib_modes     = 0
       molecule%theta_d       = zero
       molecule%molecule_type = monatomic

    case( "Ar" )
       molecule%diam  = 4.17e-10
       molecule%mass  = 66.3e-27
       molecule%omega = 0.81d0
       molecule%alpha = 1.40d0

       molecule%rot_modes     = 0
       molecule%vib_modes     = 0
       molecule%theta_d       = zero
       molecule%molecule_type = monatomic

    case( "Xe" )
       molecule%diam  = 5.74e-10
       molecule%mass  = 218.e-27
       molecule%omega = 0.85d0
       molecule%alpha = 1.44d0

       molecule%rot_modes     = 0
       molecule%vib_modes     = 0
       molecule%theta_d       = zero
       molecule%molecule_type = monatomic

    case( "H2" )
       molecule%diam  = 2.92e-10
       molecule%mass  = 3.34e-27
       molecule%omega = 0.67d0
       molecule%alpha = 1.35d0

       molecule%rot_modes     = 1
       molecule%vib_modes     = 1
       molecule%molecule_type = homonuclear

       molecule%theta_r = 85.4d0
       molecule%theta_v = 6159.0d0
       molecule%theta_d = 52000.0d0

       molecule%cZr  = 10.0d0
       molecule%cZv  = 10.0d0
       molecule%Zr_inf = one
       molecule%T_star = one
       molecule%Zv_1 = one
       molecule%Zv_2 = one
       molecule%vib_cs = one

       molecule%even_spin = 1.0d0
       molecule%odd_spin  = 3.0d0

       molecule%Y10    = 4401.0d0
       molecule%Y20    = -121.3d0
       molecule%Y01    = 60.8d0
       molecule%Y02    = -1.6d-2

    case( "N2" )
       molecule%diam  = 4.17e-10
       molecule%mass  = 46.5e-27
       molecule%omega = 0.74d0
       molecule%alpha = 1.36d0

       molecule%rot_modes     = 1
       molecule%vib_modes     = 1
       molecule%molecule_type = homonuclear

       molecule%theta_r = 2.88d0
       molecule%theta_v = 3371.0d0
       molecule%theta_d = 113500.0d0

       molecule%cZr  = 5.0d0
       molecule%cZv  = 1000.0d0
       molecule%Zr_inf = 15.7d0
       molecule%T_star = 80.0d0
       molecule%Zv_1 = 9.1d0
       molecule%Zv_2 = 220.0d0
       molecule%vib_cs = 1.4d0

       molecule%even_spin = 6.0d0
       molecule%odd_spin  = 3.0d0

       molecule%Y10    = 2359.0d0
       molecule%Y20    = -14.3d0
       molecule%Y01    = 2.01d0
       molecule%Y02    = -5.8d-6

    case( "N2_rot" )
       molecule%diam  = 4.17e-10
       molecule%mass  = 46.5e-27
       molecule%omega = 0.74d0
       molecule%alpha = 1.36d0

       molecule%rot_modes     = 1
       molecule%vib_modes     = 0
       molecule%molecule_type = homonuclear

       molecule%theta_r = 2.88d0
       molecule%theta_v = 3371.0d0
       molecule%theta_d = 113500.0d0

       molecule%cZr  = 5.0d0
       molecule%cZv  = 5.0d0
       molecule%Zr_inf = 15.7d0
       molecule%T_star = 80.0d0
       molecule%Zv_1 = 9.1d0
       molecule%Zv_2 = 220.0d0
       molecule%vib_cs = 1.4d0

       molecule%even_spin = 6.0d0
       molecule%odd_spin  = 3.0d0

       molecule%Y10    = 2359.0d0
       molecule%Y20    = -14.3d0
       molecule%Y01    = 2.01d0
       molecule%Y02    = -5.8d-6

    case( "N2_vib" )
       molecule%diam  = 4.17e-10
       molecule%mass  = 46.5e-27
       molecule%omega = 0.74d0
       molecule%alpha = 1.36d0

       molecule%rot_modes     = 0
       molecule%vib_modes     = 1
       molecule%molecule_type = homonuclear

       molecule%theta_r = 2.88d0
       molecule%theta_v = 3371.0d0
       molecule%theta_d = 113500.0d0

       molecule%cZr  = 5.0d0
       molecule%cZv  = 5.0d0
       molecule%Zr_inf = 15.7d0
       molecule%T_star = 80.0d0
       molecule%Zv_1 = 9.1d0
       molecule%Zv_2 = 220.0d0
       molecule%vib_cs = 1.4d0

       molecule%even_spin = 6.0d0
       molecule%odd_spin  = 3.0d0

       molecule%Y10    = 2359.0d0
       molecule%Y20    = -14.3d0
       molecule%Y01    = 2.01d0
       molecule%Y02    = -5.8d-6

    case( "O2" )
       molecule%diam  = 4.07e-10
       molecule%mass  = 53.12e-27
       molecule%omega = 0.77d0
       molecule%alpha = 1.40d0

       molecule%rot_modes     = 1
       molecule%vib_modes     = 1
       molecule%molecule_type = homonuclear

       molecule%theta_r = 2.07d0
       molecule%theta_v = 2256.0d0
       molecule%theta_d = 59500.0d0

       molecule%cZr  = 10.0d0
       molecule%cZv  = 10.0d0
       molecule%Zr_inf = 14.4d0
       molecule%T_star = 90.0d0
       molecule%Zv_1 = 56.5d0
       molecule%Zv_2 = 153.5d0
       molecule%vib_cs = 1.3d0

       molecule%even_spin = 0.0d0
       molecule%odd_spin  = 1.0d0

       molecule%Y10    = 1580.0d0
       molecule%Y20    = -12.0d0
       molecule%Y01    = 1.45d0
       molecule%Y02    = -4.8d-6

    case( "CO" )
       molecule%diam  = 4.19e-10
       molecule%mass  = 46.5e-27
       molecule%omega = 0.73d0
       molecule%alpha = 1.49d0

       molecule%rot_modes     = 1
       molecule%vib_modes     = 1
       molecule%molecule_type = heteronuclear

       molecule%theta_r = 2.77d0
       molecule%theta_v = 3103.0d0
       molecule%theta_d = 29700.0d0

       molecule%cZr  = 10.0d0
       molecule%cZv  = 10.0d0
       molecule%Zr_inf = 0.0d0
       molecule%T_star = 0.0d0
       molecule%Zv_1 = 37.7d0
       molecule%Zv_2 = 175.0d0
       molecule%vib_cs = 0.0d0

       molecule%even_spin = 0.0d0
       molecule%odd_spin  = 0.0d0

       molecule%Y10    = 2170.0d0
       molecule%Y20    = -13.3d0
       molecule%Y01    = 1.93d0
       molecule%Y02    = -6.1d-6

    case default
       write(*,*) "Error, invalid value of species_code. Given value is: ", species_code_in
       stop

    end select

    ! Adjust constants for collision model
    call get_csection_model( csection_model )

    select case( csection_model )
    case( hard_sphere )
       molecule%omega = 0.5d0
       molecule%alpha = 1.0d0

    case( pseudo_maxwell )
       molecule%omega = 1.0d0
       molecule%alpha = 1.0d0

    case( vhs )
       molecule%alpha = 1.0d0

    case( vss )
       ! nothing changed

    end select

    return
  end subroutine get_molecule_info

end module SpeciesAndReferenceData
