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
     
  end type MoleculeType

  character(len=16), allocatable, dimension(:) :: species_codes

  character(len=16) :: ref_species
  double precision :: temp_ref, dens_ref, diam_ref, mass_ref

  public :: MoleculeType

  public :: set_species_codes
  public :: set_reference_parameters
  public :: set_molecule_constants
  public :: get_reference_mol_data
  public :: get_species_rot_dof

!!$  public :: get_molecule_data
  public :: get_modes_data
!!$  public :: get_characteristic_temp_data
!!$  public :: get_molecular_type
!!$  public :: get_relaxation_rate_data

contains

!!$  subroutine get_molecule_data( species_code, mass, diam, omega, alpha )
!!$
!!$    implicit none
!!$
!!$    character(len=16), intent(in) :: species_code
!!$
!!$    double precision :: mass, diam, omega, alpha
!!$
!!$    select case( species_code )
!!$
!!$       ! diameter - [m]
!!$       ! mass     - [kg]
!!$       ! omega    - VHS coefficient
!!$       ! alpha    - VSS exponent
!!$
!!$    case( "He" )
!!$       diam  = 2.33e-10
!!$       mass  = 6.65e-27
!!$       omega = 0.66d0
!!$       alpha = 1.26d0
!!$
!!$    case( "Ne" )
!!$       diam  = 2.77e-10
!!$       mass  = 33.5e-27
!!$       omega = 0.66d0
!!$       alpha = 1.31d0
!!$
!!$    case( "Ar" )
!!$       diam  = 4.17e-10
!!$       mass  = 66.3e-27
!!$       omega = 0.81d0
!!$       alpha = 1.40d0
!!$
!!$    case( "Xe" )
!!$       diam  = 5.74e-10
!!$       mass  = 218.e-27
!!$       omega = 0.85d0
!!$       alpha = 1.44d0
!!$
!!$    case( "H2" )
!!$       
!!$
!!$    case( "N2", "N2_rot", "N2_vib" )
!!$       diam  = 4.17e-10
!!$       mass  = 46.5e-27
!!$       omega = 0.74d0
!!$       alpha = 1.36d0
!!$
!!$    case( "O2" )
!!$       diam  = 4.07e-10
!!$       mass  = 53.12e-27
!!$       omega = 0.77d0
!!$       alpha = 1.40d0
!!$
!!$    case default
!!$       write(*,*) "Error, invalid value of species_code. Given value is: ", species_code
!!$       stop
!!$
!!$    end select
!!$
!!$    return
!!$  end subroutine get_molecule_data

  subroutine get_modes_data( species_code, r_modes, v_modes )

    implicit none

    character(len=16), intent(in) :: species_code

    integer :: r_modes, v_modes

    select case( species_code )
       
       ! r_modes - number of rotational energy modes
       ! v_modes - number of vibrational energy modes

    case( "He", "Ne", "Ar", "Xe" )
       r_modes = 0
       v_modes = 0

    case( "H2", "N2", "O2" )
       r_modes = 1
       v_modes = 1

    case( "N2_rot" )
       r_modes = 1
       v_modes = 0

    case( "N2_vib" )
       r_modes = 0
       v_modes = 1

    case default
       write(*,*) "Error, invalid value of species_code. Given value is: ", species_code
       stop

    end select

    return
  end subroutine get_modes_data

!!$  subroutine get_characteristic_temp_data( species_code, theta_r, theta_v, theta_d )
!!$
!!$    implicit none
!!$
!!$    character(len=16), intent(in) :: species_code
!!$
!!$    double precision :: theta_r, theta_v, theta_d
!!$
!!$    select case( species_code )
!!$       
!!$       ! Characteristic temperatures - [K]
!!$
!!$    case( "He", "Ne", "Ar", "Xe" )
!!$       theta_r = zero
!!$       theta_v = zero
!!$       theta_d = zero
!!$
!!$    case( "N2", "N2_rot", "N2_vib" )
!!$       theta_r = 2.88d0
!!$       theta_v = 3390.0d0
!!$       theta_d = 113500.0d0
!!$
!!$    case( "O2" )
!!$       theta_r = 2.07d0
!!$       theta_v = 2256.0d0
!!$       theta_d = 59500.0d0
!!$
!!$    case default
!!$       write(*,*) "Error, invalid value of species_code. Given value is: ", species_code
!!$       stop
!!$
!!$    end select
!!$
!!$    return
!!$  end subroutine get_characteristic_temp_data

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

!!$  subroutine get_molecular_type( species_code, mol_type, even_spin, odd_spin )
!!$
!!$    implicit none
!!$    
!!$    character(len=16), intent(in) :: species_code
!!$
!!$    integer :: mol_type
!!$
!!$    double precision :: even_spin, odd_spin
!!$
!!$    select case( species_code )
!!$       
!!$       ! shape of molecules (homonuclear or heteronuclear)
!!$
!!$    case( "He", "Ne", "Ar", "Xe" )
!!$       mol_type = monatomic
!!$       
!!$    case( "N2", "N2_rot", "N2_vib" )
!!$       mol_type = homonuclear
!!$
!!$       even_spin = 6.0d0
!!$       odd_spin  = 3.0d0
!!$
!!$    case( "O2" )
!!$       mol_type = homonuclear
!!$
!!$       even_spin = 0.0d0
!!$       odd_spin  = 1.0d0
!!$
!!$    case default
!!$       write(*,*) "Error, invalid value of species_code. Given value is: ", species_code
!!$       stop
!!$
!!$    end select
!!$
!!$    return
!!$  end subroutine get_molecular_type

!!$  subroutine get_relaxation_rate_data( Cr, Cv1, Cv2, species_code )
!!$
!!$    implicit none
!!$
!!$    character(len=16), intent(in) :: species_code
!!$
!!$    double precision :: Cr, Cv1, Cv2
!!$
!!$    select case( species_code )
!!$    case( "He", "Ne", "Ar", "Xe" )
!!$       ! No rotational or vibrational energy
!!$
!!$    case( "N2", "N2_rot", "N2_vib" )
!!$       Cr  = 5.0d0
!!$       Cv1 = 10.0d0
!!$       Cv2 = 0.0d0
!!$
!!$    case( "O2" )
!!$       Cr  = 5.0d0
!!$       Cv1 = 1000.0d0
!!$       Cv2 = 0.0d0
!!$
!!$    case default
!!$       write(*,*) "Error, invalid value of species_code. Given value is: ", species_code
!!$       stop
!!$
!!$    end select
!!$
!!$    return
!!$  end subroutine get_relaxation_rate_data

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

    call get_molecule_info( ref_molecule, ref_species )
    diam_ref = ref_molecule%diam
    mass_ref = ref_molecule%mass

    call set_scaling_parameters( dens_ref, temp_ref, mass_ref, diam_ref )

    do n = 1, num_species
       call get_molecule_info( molecule(n), species_codes(n) )

       ! Scaling and calculation of time independent constants
       molecule(n)%diam = molecule(n)%diam/diam_ref
       molecule(n)%mass = molecule(n)%mass/mass_ref

       if( molecule(n)%rot_modes .gt. 0 )then
          molecule(n)%theta_r = molecule(n)%theta_r/temp_ref
          molecule(n)%Zr_inf  = molecule(n)%Zr_inf/temp_ref
          molecule(n)%Y01     = molecule(n)%Y01/(cm2J*two*kb*temp_ref)
          molecule(n)%Y02     = molecule(n)%Y02/(cm2J*two*kb*temp_ref)
       end if

       if( molecule(n)%vib_modes .gt. 0 )then
          molecule(n)%Av      = &
               1.16d-3*molecule(n)%theta_v**(4.0d0/3.0d0)*(mass_ref/kg2amu)**(0.5d0)*temp_ref**(-one_third)
          molecule(n)%Bv      = -1.74d-5*molecule(n)%theta_v**(4.0d0/3.0d0)*(mass_ref/kg2amu)**(0.75d0)
          molecule(n)%vib_cs  = molecule(n)%vib_cs*angstrom*angstrom/diam_ref/diam_ref/pi
          molecule(n)%theta_v = molecule(n)%theta_v/temp_ref
          molecule(n)%theta_d = molecule(n)%theta_d/temp_ref
          molecule(n)%Y10     = molecule(n)%Y10/(cm2J*two*kb*temp_ref)
          molecule(n)%Y20     = molecule(n)%Y20/(cm2J*two*kb*temp_ref)
       end if

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

    implicit none

    character(len=16), intent(in) :: species_code_in

    type(MoleculeType) :: molecule

    molecule%species_code = species_code_in

    select case( species_code_in )
    case( "He" )
       molecule%diam  = 2.33e-10
       molecule%mass  = 6.65e-27
       molecule%omega = 0.66d0
       molecule%alpha = 1.31d0

       molecule%rot_modes     = 0
       molecule%vib_modes     = 0
       molecule%molecule_type = monatomic

    case( "Ne" )
       molecule%diam  = 2.77e-10
       molecule%mass  = 33.5e-27
       molecule%omega = 0.66d0
       molecule%alpha = 1.31d0

       molecule%rot_modes     = 0
       molecule%vib_modes     = 0
       molecule%molecule_type = monatomic

    case( "Ar" )
       molecule%diam  = 4.17e-10
       molecule%mass  = 66.3e-27
       molecule%omega = 0.81d0
       molecule%alpha = 1.40d0

       molecule%rot_modes     = 0
       molecule%vib_modes     = 0
       molecule%molecule_type = monatomic

    case( "Xe" )
       molecule%diam  = 5.74e-10
       molecule%mass  = 218.e-27
       molecule%omega = 0.85d0
       molecule%alpha = 1.44d0

       molecule%rot_modes     = 0
       molecule%vib_modes     = 0
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

       molecule%cZr  = 10.0d0
       molecule%cZv  = 10.0d0
       molecule%Zr_inf = 80.0d0
       molecule%T_star = 15.7d0
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
       molecule%cZv  = 200.0d0
       molecule%Zr_inf = 80.0d0
       molecule%T_star = 15.7d0
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
       molecule%cZv  = 10.0d0
       molecule%Zr_inf = 80.0d0
       molecule%T_star = 15.7d0
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

       molecule%cZr  = 5.0d0
       molecule%cZv  = 200.0d0
       molecule%Zr_inf = 90.0d0
       molecule%T_star = 14.4d0
       molecule%Zv_1 = 56.5d0
       molecule%Zv_2 = 153.5d0
       molecule%vib_cs = 1.3d0

       molecule%even_spin = 0.0d0
       molecule%odd_spin  = 1.0d0

       molecule%Y10    = 1580.0d0
       molecule%Y20    = -12.0d0
       molecule%Y01    = 1.45d0
       molecule%Y02    = -4.8d-6

    case default
       write(*,*) "Error, invalid value of species_code. Given value is: ", species_code_in
       stop

    end select

    return
  end subroutine get_molecule_info

end module SpeciesAndReferenceData
