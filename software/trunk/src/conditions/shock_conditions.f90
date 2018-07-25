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
module ShockConditions

  use Constants
  use BoundaryConditions
  use DistFunc

  implicit none

  private

  type NormalShock
     double precision :: dens_ratio, T_ratio, P_ratio
     double precision :: u_down, M_down
     double precision :: u_shock, M_infty
  end type NormalShock

  integer, parameter :: specular_reflection = 1
  integer, parameter :: diffuse_reflection  = 2
  integer, parameter :: fixed_in_out_flow   = 3
  integer, parameter :: zero_gradient_flow  = 4

  public :: NormalShock

  public :: initialize_normal_shock_props
  public :: calculate_molecular_mass
  public :: calculate_specific_heat_ratio

contains

  subroutine initialize_normal_shock_props( shock_props, molecule )

    use BoundaryConditions
    use PhysicalGrid
    use MathUtilities
    use SpeciesAndReferenceData

    implicit none

    type(MoleculeType), dimension(:), intent(in) :: molecule

    type(NormalShock) :: shock_props

    double precision :: dens_up, T_up, u_up, P_up
    double precision :: mol_mass, gamma

    double precision :: u_s, a, b, c
    double precision :: a_up, a_down

    double precision :: M_down_local, M_infty_local, T_down_local, P_down_local

    integer :: LWall, RWall, BWall, TWall
    integer :: nx_space, ny_space

    ! Check nspace > 1
    call get_nspace( nx_space, ny_space )
    if( nx_space .eq. 1 ) return

    ! Get flags to test for wall at boundary
    call get_boundary_conditions_flags( LWall, RWall, BWall, TWall )

    ! Calculate upstream conditions, molecular mass, and specific heat ratio
    call calculate_upstream_conditions( dens_up, T_up, u_up, molecule )
    call calculate_molecular_mass( mol_mass, molecule )
    call calculate_specific_heat_ratio( gamma, molecule )
    
    P_up = dens_up*T_up!*one_half/mol_mass

    write(*,*) "Specific heat ratio  = ", gamma
    write(*,*) "Molecular mass       = ", mol_mass

    write(*,*)"Upstream Conditions: "
    write(*,'(a15,e16.10)') "density     = ",dens_up
    write(*,'(a15,e16.10)') "velocity    = ",u_up
    write(*,'(a15,e16.10)') "temperature = ",T_up
    write(*,'(a15,e16.10)') "pressure    = ",P_up
   
    a_up = sqrt( gamma*T_up/( two*mol_mass ) )
    
    if( LWall .eq. specular_reflection .or. RWall .eq. specular_reflection )then
       a = one

       b = one_half*( 3.0d0 - gamma )*u_up

       c = -one_half*( gamma + ( gamma - one )*u_up*u_up )

       u_s = abs( ( -b - sqrt( b*b - 4.0d0*a*c ) )/( two*a ) )
       shock_props%M_infty = ( abs( u_up ) + u_s )/a_up


    else
       if( u_up .gt. zero )then
          write(*,*) "Error. Flow going in the wrong direction, flow for stationary shock must be negative."
          stop
       end if

       u_s = zero
       shock_props%M_infty = abs( u_up )/a_up

    end if

    shock_props%u_shock = u_s

    M_infty_local = shock_props%M_infty

    write(*,*)"Mach number seen by shock = ",M_infty_local
    if( M_infty_local .lt. double_tol )then
       write(*,*)"Local Mach number is near zero. The flow isn't moving in x-direction"
       shock_props%T_ratio = one
       shock_props%P_ratio = one
       shock_props%dens_ratio = one
       shock_props%M_down = zero
       shock_props%u_down = zero
       return
    end if

    shock_props%T_ratio = ( one + two*gamma*( M_infty_local*M_infty_local - one )/( gamma + one ) )*&
         ( two + ( gamma - one )*M_infty_local*M_infty_local )/&
         ( ( gamma + one )*M_infty_local*M_infty_local )

    T_down_local = shock_props%T_ratio*T_up

    shock_props%P_ratio = ( one + two*gamma*( M_infty_local*M_infty_local - one )/( gamma + one ) )

    P_down_local = shock_props%P_ratio*P_up

    shock_props%dens_ratio = ( gamma + one )*M_infty_local*M_infty_local/&
         ( two + ( gamma - one )*M_infty_local*M_infty_local )

    shock_props%M_down = sqrt( ( one + one_half*( gamma - one )*M_infty_local*M_infty_local )/&
         ( gamma*M_infty_local*M_infty_local - one_half*( gamma - one ) ) )

    M_down_local = shock_props%M_down

    a_down = sqrt( gamma*T_down_local/( two*mol_mass ) )

    shock_props%u_down = dsgn(u_up)*a_down*M_down_local + u_s

    write(*,*)" "
    write(*,*)"Shock Jump Conditions: "
    write(*,*)"density ratio       = ",shock_props%dens_ratio
    write(*,*)"temperature ratio   = ",shock_props%T_ratio
    write(*,*)"pressure ratio      = ",shock_props%P_ratio
    write(*,*)"downstream velocity = ",shock_props%u_down
    write(*,*)"downstream Mach#    = ",M_down_local

    return
  end subroutine initialize_normal_shock_props

  subroutine calculate_upstream_conditions( dens_up, T_up, u_up, molecule )

    use SpeciesAndReferenceData

    implicit none

    double precision :: dens_up, T_up, u_up

    type(MoleculeType), dimension(:), intent(in) :: molecule

    double precision :: dens_tmp, T_tmp, u_tmp, v_tmp, w_tmp, tot_mass, mass

    integer :: n
    
    ! It is understood that for a shock, 
    ! every species velocity must be the same value at the right wall and
    ! every species temperature must be the same value at the right wall

    dens_up = zero
    T_up    = zero
    u_up    = zero
    tot_mass = zero

    do n = 1, num_species

       call get_right_wall_velocity( u_tmp, v_tmp, w_tmp, n )
       call get_right_wall_density( dens_tmp, n )
       call get_right_wall_temp( T_tmp, n )
       if( u_tmp .gt. zero )then
          call get_left_wall_velocity( u_tmp, v_tmp, w_tmp, n )
          call get_left_wall_density( dens_tmp, n )
          call get_left_wall_temp( T_tmp, n )
       end if
       mass = molecule(n)%mass

       dens_up = dens_up + dens_tmp
       T_up    = T_up + T_tmp*dens_tmp
       u_up = u_up + u_tmp*mass*dens_tmp

       tot_mass = tot_mass + dens_tmp*mass

    end do

    T_up = T_up/dens_up
    u_up = u_up/tot_mass

    return
  end subroutine calculate_upstream_conditions

  subroutine calculate_molecular_mass( mol_mass, molecule )

    use SpeciesAndReferenceData

    implicit none

    double precision :: mol_mass

    type(MoleculeType), dimension(:), intent(in) :: molecule

    double precision :: mass, dens, tot_dens

    integer :: n
    
    mol_mass = zero
    tot_dens = zero

    do n = 1, num_species

       call get_right_wall_density( dens, n )
       mass = molecule(n)%mass

       tot_dens = tot_dens + dens

       mol_mass = mol_mass + dens*mass

    end do

    mol_mass = mol_mass/tot_dens

    return
  end subroutine calculate_molecular_mass

  subroutine calculate_specific_heat_ratio( gamma, molecule )

    use SpeciesAndReferenceData

    implicit none
    
    double precision :: gamma

    type(MoleculeType), dimension(:), intent(in) :: molecule

    double precision :: vib_dof, temp_v, theta_v
    double precision :: mass, dens, dof, rot_dof, total_mass
    double precision :: cp_mix

    integer :: r_modes, v_modes

    character(len=16) :: species_code

    integer :: n

    cp_mix = zero
    total_mass = zero

    do n = 1, num_species

       r_modes = molecule(n)%rot_modes
       v_modes = molecule(n)%vib_modes

       call get_right_wall_density( dens, n )
       mass = molecule(n)%mass
       species_code = molecule(n)%species_code
       if( r_modes .gt. 0 )then
          call get_species_rot_dof( species_code, rot_dof )
       else
          rot_dof = 0.0d0
       end if
       
       if( v_modes .gt. 0 )then
          call get_right_wall_temp_vib( temp_v, n )
          theta_v = molecule(n)%theta_v
          vib_dof = ( two * theta_v / temp_v ) / ( exp( theta_v / temp_v ) - one )
       else
          vib_dof = 0.0d0
       end if

       total_mass = total_mass + mass*dens

       dof = rot_dof + 3.0d0 + vib_dof

       cp_mix = cp_mix + mass*dens*( dof*one_half + one )


    end do

    cp_mix = cp_mix/total_mass

    gamma = cp_mix/( cp_mix - one )

    return
  end subroutine calculate_specific_heat_ratio

end module ShockConditions
