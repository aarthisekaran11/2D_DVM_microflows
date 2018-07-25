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
module ReplenishingCollisions

  use Constants
  use ErrorCheck
  use CollisionUtilities

  implicit none

  private

  integer :: num_repl_locations

  integer, parameter :: x_only  = 1
  integer, parameter :: y_only  = 2
  integer, parameter :: z_only  = 4
  integer, parameter :: xy_only = 3
  integer, parameter :: xz_only = 5
  integer, parameter :: yz_only = 6
  integer, parameter :: xyz     = 7

  integer, parameter :: elastic       = 0
  integer, parameter :: rot_trans     = 1
  integer, parameter :: vib_trans     = 2
  integer, parameter :: rot_vib_trans = 3

  integer, parameter :: A = 1
  integer, parameter :: B = 2

  double precision, allocatable, dimension(:) :: fr1_post, fr2_post, fv1_post, fv2_post
  double precision, allocatable, dimension(:) :: fr3_post, fr4_post, fv3_post, fv4_post

  public :: set_num_repl_locations
  public :: create_internal_energy_replenish_arrays
  public :: destroy_internal_energy_replenish_arrays
  public :: replenish_collision

contains

  subroutine set_num_repl_locations( num_repl_locations_in )

    implicit none

    integer, intent(in) :: num_repl_locations_in

    num_repl_locations = num_repl_locations_in

    return
  end subroutine set_num_repl_locations

  subroutine create_internal_energy_replenish_arrays( r_modes, v_modes, r_levels, v_levels )

    implicit none

    integer, dimension(2), intent(in) :: r_modes, v_modes, r_levels, v_levels

    integer :: status

    if( r_modes(1) .gt. 0 .or. r_modes(2) .gt. 0 )then
       allocate( fr1_post( 1:r_levels(1) ), fr3_post( 1:r_levels(1) ), STAT=status )
       call allocate_error_check( status, "post collision fr1 and fr3" )

       allocate( fr2_post( 1:r_levels(2) ), fr4_post( 1:r_levels(2) ), STAT=status )
       call allocate_error_check( status, "post collision fr2 and fr4" )

    end if

    if( v_modes(1) .gt. 0 .or. v_modes(2) .gt. 0 )then
       allocate( fv1_post( 1:v_levels(1) ), fv3_post( 1:v_levels(1) ), STAT=status )
       call allocate_error_check( status, "post collision fv1 and fv3" )

       allocate( fv2_post( 1:v_levels(2) ), fv4_post( 1:v_levels(2) ), STAT=status )
       call allocate_error_check( status, "post collision fv2 and fv4" )

    end if

    return
  end subroutine create_internal_energy_replenish_arrays

  subroutine destroy_internal_energy_replenish_arrays( r_modes, v_modes )

    implicit none

    integer, dimension(2), intent(in) :: r_modes, v_modes

    integer :: status

    if( r_modes(1) .gt. 0 )then
       deallocate( fr1_post, fr3_post, STAT=status )
       call deallocate_error_check( status, "post collision fr1 and fr3" )
    end if

    if( r_modes(2) .gt. 0 )then
       deallocate( fr2_post, fr4_post, STAT=status )
       call deallocate_error_check( status, "post collision fr2 and fr4" )
    end if

    if( v_modes(1) .gt. 0 )then
       deallocate( fv1_post, fv3_post, STAT=status )
       call deallocate_error_check( status, "post collision fv1 and fv3" )
    end if

    if( v_modes(2) .gt. 0 )then
       deallocate( fv2_post, fv4_post, STAT=status )
       call deallocate_error_check( status, "post collision fv2 and fv4" )
    end if

    return
  end subroutine destroy_internal_energy_replenish_arrays

  subroutine replenish_collision( phi, depletion, kin_depl, fr1, fr2, fv1, fv2, r_level, v_level, i, j, k, &
       g, n, m, vel_grid, m_red, coll_type, exchange_mol, molecule )

    use DistFunc
    use SpeciesAndReferenceData
    use VelocityGrid
    use ReMapping

    implicit none

    type(DistFuncType), dimension(:) :: phi

    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    type(MoleculeType), dimension(:), intent(in) :: molecule
    double precision, dimension(:), intent(in) :: fr1, fr2, fv1, fv2
    double precision, intent(in) :: g
    double precision, dimension(2), intent(in) :: depletion, kin_depl
    integer, dimension(2), intent(in) :: i, j, k
    integer, dimension(4), intent(in) :: r_level, v_level
    integer, intent(in) :: coll_type, n, m, exchange_mol

    double precision, dimension(2) :: x, y, z
    double precision, dimension(2) :: mass, omega, alpha, theta_r, theta_v
    double precision, dimension(2) :: repl_frac, kin_repl
    double precision, dimension(2) :: er_prime
    double precision, dimension(3) :: v_com, xi, eta
    double precision :: delta_eps, g_prime, g_ratio
    double precision :: omega_AB, alpha_AB
    double precision :: m_red

    integer, dimension(2) :: r_modes, v_modes, r_levels, v_levels, level_p, mol_type
    integer :: repl

    type(MappingResultType), dimension(2) :: mapping

    logical, dimension(2) :: single_valued_rot

    ! Set grid variables
    x(1) = vel_grid(n)%x( i(1) ) ![m/s]
    y(1) = vel_grid(n)%y( j(1) ) ![m/s]
    z(1) = vel_grid(n)%z( k(1) ) ![m/s]

    x(2) = vel_grid(m)%x( i(2) ) ![m/s]
    y(2) = vel_grid(m)%y( j(2) ) ![m/s]
    z(2) = vel_grid(m)%z( k(2) ) ![m/s]

    ! Set species variables
    mass(1)  = molecule(n)%mass ![kg/#]
    omega(1) = molecule(n)%omega
    alpha(1) = molecule(n)%alpha

    mass(2)  = molecule(m)%mass ![kg/#]
    omega(2) = molecule(m)%omega
    alpha(2) = molecule(m)%alpha

    omega_AB = one_half*( omega(1) + omega(2) )
    alpha_AB = one_half*( alpha(1) + alpha(2) )

    ! Internal structure
    r_modes(1) = molecule(n)%rot_modes
    v_modes(1) = molecule(n)%vib_modes

    r_modes(2) = molecule(m)%rot_modes
    v_modes(2) = molecule(m)%vib_modes

    r_levels(1) = phi(n)%num_rot_levels
    v_levels(1) = phi(n)%num_vib_levels
    theta_r(1)  = molecule(n)%theta_r
    theta_v(1)  = molecule(n)%theta_v
    mol_type(1) = molecule(n)%molecule_type

    r_levels(2) = phi(m)%num_rot_levels
    v_levels(2) = phi(m)%num_vib_levels
    theta_r(2)  = molecule(m)%theta_r
    theta_v(2)  = molecule(m)%theta_v
    mol_type(2) = molecule(m)%molecule_type

    single_valued_rot(1) = phi(n)%single_valued_rot
    single_valued_rot(2) = phi(m)%single_valued_rot

    ! Calculate replenishment variables
    repl_frac = depletion/dble( num_repl_locations )
    kin_repl  = kin_depl/dble( num_repl_locations )

    !mass = kin_repl * mass
    !m_red = mass(1) * mass(2) / ( mass(1) + mass(2) )
    if( r_modes(1) .gt. 0 ) er_prime(1) = fr1(1)
    if( r_modes(2) .gt. 0 ) er_prime(2) = fr2(1)

    call compute_center_of_mass_velocity( x, y, z, mass, v_com )

    do repl = 1, num_repl_locations

       !TODO: I need a check for if the post-collision velocities are outside the velocity domain by
       !TODO: too much.

       select case( coll_type )
       case( elastic )
          delta_eps = zero
          g_prime = g
          g_ratio = one

          call find_post_collision_velocities( xi, eta, x, y, z, v_com, g_prime, g_ratio, alpha_AB, &
               vel_grid, n, m )

          call perform_remapping( mapping(1), xi(1), xi(2), xi(3), one, vel_grid(n) )
          call perform_remapping( mapping(2), eta(1), eta(2), eta(3), one, vel_grid(m) )

          call apply_phi_replenishment( phi(n), kin_repl(1), mapping(1) )
          call apply_phi_replenishment( phi(m), kin_repl(2), mapping(2) )

          call rot_replenishment_full( phi(n), repl_frac(1), mapping(1), fr1, r_modes(1) )
          call rot_replenishment_full( phi(m), repl_frac(2), mapping(2), fr2, r_modes(2) )
             
          call vib_replenishment_full( phi(n), repl_frac(1), mapping(1), fv1, v_modes(1) )
          call vib_replenishment_full( phi(m), repl_frac(2), mapping(2), fv2, v_modes(2) )

       case( rot_trans )
          if( exchange_mol .eq. A )then
             call compute_rot_level_change( level_p(1), delta_eps, r_level(1), &
                  g, omega_AB, m_red, mass(1), r_levels(1), molecule(n), phi(n)%rot_level, fr1 )

          else
             call compute_rot_level_change( level_p(2), delta_eps, r_level(2), &
                  g, omega_AB, m_red, mass(2), r_levels(2), molecule(m), phi(m)%rot_level, fr2 )

          end if

          call calculate_g_prime( g_prime, g, delta_eps, m_red )
          g_ratio = g_prime/g

          call find_post_collision_velocities( xi, eta, x, y, z, v_com, g_prime, g_ratio, alpha_AB, &
               vel_grid, n, m )

          call perform_remapping( mapping(1), xi(1), xi(2), xi(3), one, vel_grid(n) )
          call perform_remapping( mapping(2), eta(1), eta(2), eta(3), one, vel_grid(m) )

          call apply_phi_replenishment( phi(n), kin_repl(1), mapping(1) )
          call apply_phi_replenishment( phi(m), kin_repl(2), mapping(2) )

          if( exchange_mol .eq. A )then
             call rot_replenishment_full( phi(m), repl_frac(2), mapping(2), fr2, r_modes(2) )
             call vib_replenishment_full( phi(m), repl_frac(2), mapping(2), fv2, v_modes(2) )
!!$                call rot_replenishment_single( phi(m), kin_repl(2), mapping(2), r_level(2), r_modes(2) )
!!$                call vib_replenishment_single( phi(m), kin_repl(2), mapping(2), v_level(4), v_modes(2) )

             if( single_valued_rot(1) .eqv. .false. )then
                call rot_replenishment_single( phi(n), kin_repl(1), mapping(1), level_p(1), r_modes(1) )
                call vib_replenishment_single( phi(n), kin_repl(1), mapping(1), v_level(3), v_modes(1) )
             else
                call rot_replenishment_full( phi(n), repl_frac(1), mapping(1), fr1 + delta_eps/mass(1), r_modes(1) )
             end if

          else
             call rot_replenishment_full( phi(n), repl_frac(1), mapping(1), fr1, r_modes(1) )
             call vib_replenishment_full( phi(n), repl_frac(1), mapping(1), fv1, v_modes(1) )
!!$             call rot_replenishment_single( phi(n), kin_repl(1), mapping(1), r_level(1), r_modes(1) )
!!$             call vib_replenishment_single( phi(n), kin_repl(1), mapping(1), v_level(3), v_modes(1) )

             if( single_valued_rot(2) .eqv. .false. )then
                call rot_replenishment_single( phi(m), kin_repl(2), mapping(2), level_p(2), r_modes(2) )
                call vib_replenishment_single( phi(m), kin_repl(2), mapping(2), v_level(4), v_modes(2) )
             else
                call rot_replenishment_full( phi(m), repl_frac(2), mapping(2), fr2 + delta_eps/mass(2), r_modes(2) )
             end if

          end if
          
       case( vib_trans )
          if( exchange_mol .eq. A )then
             call compute_vib_level_change( level_p(1), delta_eps, v_level(1), &
                  g, omega_AB, m_red, mass(1), v_levels(1), molecule(n), phi(n)%vib_level )
    
          else
             call compute_vib_level_change( level_p(2), delta_eps, v_level(2), &
                  g, omega_AB, m_red, mass(2), v_levels(2), molecule(m), phi(m)%vib_level )

          end if

          call calculate_g_prime( g_prime, g, delta_eps, m_red )
          g_ratio = g_prime/g

          call find_post_collision_velocities( xi, eta, x, y, z, v_com, g_prime, g_ratio, alpha_AB, &
               vel_grid, n, m )

          call perform_remapping( mapping(1), xi(1), xi(2), xi(3), one, vel_grid(n) )
          call perform_remapping( mapping(2), eta(1), eta(2), eta(3), one, vel_grid(m) )

          call apply_phi_replenishment( phi(n), kin_repl(1), mapping(1) )
          call apply_phi_replenishment( phi(m), kin_repl(2), mapping(2) )

          if( exchange_mol .eq. A )then
             call vib_replenishment_single( phi(n), kin_repl(1), mapping(1), level_p(1), v_modes(1) )
             call rot_replenishment_single( phi(n), kin_repl(1), mapping(1), r_level(3), r_modes(1) )

!!$             call vib_replenishment_single( phi(m), kin_repl(2), mapping(2), v_level(2), v_modes(2) )
!!$             call rot_replenishment_single( phi(m), kin_repl(2), mapping(2), r_level(4), r_modes(2) )
!!$
             call vib_replenishment_full( phi(m), repl_frac(2), mapping(2), fv2, v_modes(2) )
             call rot_replenishment_full( phi(m), repl_frac(2), mapping(2), fr2, r_modes(2) )

          else
             call vib_replenishment_full( phi(n), repl_frac(1), mapping(1), fv1, v_modes(1) )
             call rot_replenishment_full( phi(n), repl_frac(1), mapping(1), fr1, r_modes(1) )
!!$             
!!$             call vib_replenishment_single( phi(n), kin_repl(1), mapping(1), v_level(1), v_modes(1) )
!!$             call rot_replenishment_single( phi(n), kin_repl(1), mapping(1), r_level(3), r_modes(1) )

             call vib_replenishment_single( phi(m), kin_repl(2), mapping(2), level_p(2), v_modes(2) )
             call rot_replenishment_single( phi(m), kin_repl(2), mapping(2), r_level(4), r_modes(2) )

          end if

       case default

       end select

    end do

    return
  end subroutine replenish_collision

  subroutine apply_phi_replenishment( phi, repl, mapping )

    use DistFunc
    use Remapping
    use TimeStepping

    implicit none

    type(DistFuncType) :: phi

    type(MappingResultType), intent(in) :: mapping
    double precision, intent(in) :: repl

    double precision :: repl_o, repl_ix, repl_iy, repl_iz, repl_ex, repl_ey, repl_ez

    integer :: i, j, k, a, b, c, ext_flag

    ! read in stencil locations
    i = mapping%stencil(1)
    j = mapping%stencil(2)
    k = mapping%stencil(3)
    a = mapping%stencil(4)
    b = mapping%stencil(5)
    c = mapping%stencil(6)
    ext_flag = mapping%stencil(7)

    ! read in replenishment values
    repl_o  = repl * mapping%mass(1)
    repl_ix = repl * mapping%mass(2)
    repl_iy = repl * mapping%mass(3)
    repl_iz = repl * mapping%mass(4)
    repl_ex = repl * mapping%mass(5)
    repl_ey = repl * mapping%mass(6)
    repl_ez = repl * mapping%mass(7)

    ! Origin and interior points
    phi%value( i, j, k )     = phi%value( i, j, k ) + repl_o
    phi%value( i + a, j, k ) = phi%value( i + a, j, k ) + repl_ix
    phi%value( i, j + b, k ) = phi%value( i, j + b, k ) + repl_iy
    phi%value( i, j, k + c ) = phi%value( i, j, k + c ) + repl_iz

    ! Exterior points
    select case( ext_flag )
    case( x_only )
       phi%value( i - a, j, k ) = phi%value( i - a, j, k ) + repl_ex

    case( y_only )
       phi%value( i, j - b, k ) = phi%value( i, j - b, k ) + repl_ey

    case( z_only )
       phi%value( i, j, k - c ) = phi%value( i, j, k - c ) + repl_ez

    case( xy_only )
       phi%value( i - a, j, k ) = phi%value( i - a, j, k ) + repl_ex
       phi%value( i, j - b, k ) = phi%value( i, j - b, k ) + repl_ey

    case( xz_only )
       phi%value( i - a, j, k ) = phi%value( i - a, j, k ) + repl_ex
       phi%value( i, j, k - c ) = phi%value( i, j, k - c ) + repl_ez

    case( yz_only )
       phi%value( i, j - b, k ) = phi%value( i, j - b, k ) + repl_ey
       phi%value( i, j, k - c ) = phi%value( i, j, k - c ) + repl_ez

    case( xyz )
       phi%value( i - a, j, k ) = phi%value( i - a, j, k ) + repl_ex
       phi%value( i, j - b, k ) = phi%value( i, j - b, k ) + repl_ey
       phi%value( i, j, k - c ) = phi%value( i, j, k - c ) + repl_ez

    case default
       !TODO: Consistent error handling
       write(*,*) "Error: Invalid value of ext_flag - ",&
            ext_flag
       stop

    end select

    return
  end subroutine apply_phi_replenishment

  subroutine rot_replenishment_full( phi, repl, mapping, fr, r_modes )

    use DistFunc
    use Remapping
    use TimeStepping

    implicit none

    type(DistFuncType) :: phi

    type(MappingResultType), intent(in) :: mapping
    double precision, intent(in) :: repl

    double precision, dimension(:), intent(in) :: fr
    integer, intent(in) :: r_modes

    double precision :: repl_o, repl_ix, repl_iy, repl_iz, repl_ex, repl_ey, repl_ez

    integer :: i, j, k, a, b, c, ext_flag

    !!! THIS IS THE SLOWEST PART OF THE CODE DUE TO ARRAY SIZES AND THE NUMBER OF CALLS

    if( r_modes .eq. 0 )return

    ! read in stencil locations
    i = mapping%stencil(1)
    j = mapping%stencil(2)
    k = mapping%stencil(3)
    a = mapping%stencil(4)
    b = mapping%stencil(5)
    c = mapping%stencil(6)
    ext_flag = mapping%stencil(7)

    ! read in replenishment values
    repl_o  = repl * mapping%mass(1)
    repl_ix = repl * mapping%mass(2)
    repl_iy = repl * mapping%mass(3)
    repl_iz = repl * mapping%mass(4)
    repl_ex = repl * mapping%mass(5)
    repl_ey = repl * mapping%mass(6)
    repl_ez = repl * mapping%mass(7)

    ! Origin and interior points
    phi%rot( :, i, j, k )     = phi%rot( :, i,j,k ) + repl_o*fr
    phi%rot( :, i + a, j, k ) = phi%rot( :, i + a, j, k ) + repl_ix*fr
    phi%rot( :, i, j + b, k ) = phi%rot( :, i, j + b, k ) + repl_iy*fr
    phi%rot( :, i, j, k + c ) = phi%rot( :, i, j, k + c ) + repl_iz*fr

    ! Exterior points
    select case( ext_flag )
    case( x_only )
       phi%rot( :, i - a, j, k ) = phi%rot( :, i - a, j, k ) + repl_ex*fr

    case( y_only )
       phi%rot( :, i, j - b, k ) = phi%rot( :, i, j - b, k ) + repl_ey*fr

    case( z_only )
       phi%rot( :, i, j, k - c ) = phi%rot( :, i, j, k - c ) + repl_ez*fr

    case( xy_only )
       phi%rot( :, i - a, j, k ) = phi%rot( :, i - a, j, k ) + repl_ex*fr
       phi%rot( :, i, j - b, k ) = phi%rot( :, i, j - b, k ) + repl_ey*fr

    case( xz_only )
       phi%rot( :, i - a, j, k ) = phi%rot( :, i - a, j, k ) + repl_ex*fr
       phi%rot( :, i, j, k - c ) = phi%rot( :, i, j, k - c ) + repl_ez*fr

    case( yz_only )
       phi%rot( :, i, j - b, k ) = phi%rot( :, i, j - b, k ) + repl_ey*fr
       phi%rot( :, i, j, k - c ) = phi%rot( :, i, j, k - c ) + repl_ez*fr

    case( xyz )
       phi%rot( :, i - a, j, k ) = phi%rot( :, i - a, j, k ) + repl_ex*fr
       phi%rot( :, i, j - b, k ) = phi%rot( :, i, j - b, k ) + repl_ey*fr
       phi%rot( :, i, j, k - c ) = phi%rot( :, i, j, k - c ) + repl_ez*fr


    case default
       !TODO: Consistent error handling
       write(*,*) "Error: Invalid value of ext_flag - ",&
            ext_flag
       stop

    end select

    return
  end subroutine rot_replenishment_full

  subroutine vib_replenishment_full( phi, repl, mapping, fv, v_modes )

    use DistFunc
    use Remapping

    implicit none

    type(DistFuncType) :: phi

    type(MappingResultType), intent(in) :: mapping
    double precision, intent(in) :: repl

    double precision, dimension(:), intent(in) :: fv
    integer, intent(in) :: v_modes

    double precision :: repl_o, repl_ix, repl_iy, repl_iz, repl_ex, repl_ey, repl_ez

    integer :: i, j, k, a, b, c, ext_flag

    if( v_modes .eq. 0 )return

    ! read in stencil locations
    i = mapping%stencil(1)
    j = mapping%stencil(2)
    k = mapping%stencil(3)
    a = mapping%stencil(4)
    b = mapping%stencil(5)
    c = mapping%stencil(6)
    ext_flag = mapping%stencil(7)

    ! read in replenishment values
    repl_o  = repl * mapping%mass(1)
    repl_ix = repl * mapping%mass(2)
    repl_iy = repl * mapping%mass(3)
    repl_iz = repl * mapping%mass(4)
    repl_ex = repl * mapping%mass(5)
    repl_ey = repl * mapping%mass(6)
    repl_ez = repl * mapping%mass(7)

    ! Origin and interior points
    phi%vib( :, i, j, k )     = phi%vib( :, i, j, k ) + repl_o*fv
    phi%vib( :, i + a, j, k ) = phi%vib( :, i + a, j, k ) + repl_ix*fv
    phi%vib( :, i, j + b, k ) = phi%vib( :, i, j + b, k ) + repl_iy*fv
    phi%vib( :, i, j, k + c ) = phi%vib( :, i, j, k + c ) + repl_iz*fv

    ! Exterior points
    select case( ext_flag )
    case( x_only )
       phi%vib( :, i - a, j, k ) = phi%vib( :, i - a, j, k ) + repl_ex*fv

    case( y_only )
       phi%vib( :, i, j - b, k ) = phi%vib( :, i, j - b, k ) + repl_ey*fv

    case( z_only )
       phi%vib( :, i, j, k - c ) = phi%vib( :, i, j, k - c ) + repl_ez*fv

    case( xy_only )
       phi%vib( :, i - a, j, k ) = phi%vib( :, i - a, j, k ) + repl_ex*fv
       phi%vib( :, i, j - b, k ) = phi%vib( :, i, j - b, k ) + repl_ey*fv

    case( xz_only )
       phi%vib( :, i - a, j, k ) = phi%vib( :, i - a, j, k ) + repl_ex*fv
       phi%vib( :, i, j, k - c ) = phi%vib( :, i, j, k - c ) + repl_ez*fv

    case( yz_only )
       phi%vib( :, i, j - b, k ) = phi%vib( :, i, j - b, k ) + repl_ey*fv
       phi%vib( :, i, j, k - c ) = phi%vib( :, i, j, k - c ) + repl_ez*fv

    case( xyz )
       phi%vib( :, i - a, j, k ) = phi%vib( :, i - a, j, k ) + repl_ex*fv
       phi%vib( :, i, j - b, k ) = phi%vib( :, i, j - b, k ) + repl_ey*fv
       phi%vib( :, i, j, k - c ) = phi%vib( :, i, j, k - c ) + repl_ez*fv

    case default
       !TODO: Consistent error handling
       write(*,*) "Error: Invalid value of ext_flag - ",&
            ext_flag
       stop

    end select

    return
  end subroutine vib_replenishment_full

  subroutine rot_replenishment_single( phi, repl, mapping, r_level, r_modes )

    use DistFunc
    use Remapping

    implicit none

    type(DistFuncType) :: phi

    type(MappingResultType), intent(in) :: mapping
    double precision, intent(in) :: repl

    integer, intent(in) :: r_level, r_modes

    double precision :: repl_o, repl_ix, repl_iy, repl_iz, repl_ex, repl_ey, repl_ez

    integer :: i, j, k, a, b, c, ext_flag

    if( r_modes .eq. 0 )return

    ! read in stencil locations
    i = mapping%stencil(1)
    j = mapping%stencil(2)
    k = mapping%stencil(3)
    a = mapping%stencil(4)
    b = mapping%stencil(5)
    c = mapping%stencil(6)
    ext_flag = mapping%stencil(7)

    ! read in replenishment values
    repl_o  = repl * mapping%mass(1)
    repl_ix = repl * mapping%mass(2)
    repl_iy = repl * mapping%mass(3)
    repl_iz = repl * mapping%mass(4)
    repl_ex = repl * mapping%mass(5)
    repl_ey = repl * mapping%mass(6)
    repl_ez = repl * mapping%mass(7)

    ! Origin and interior points
    phi%rot( r_level, i, j, k )     = phi%rot( r_level, i,j,k ) + repl_o
    phi%rot( r_level, i + a, j, k ) = phi%rot( r_level, i + a, j, k ) + repl_ix
    phi%rot( r_level, i, j + b, k ) = phi%rot( r_level, i, j + b, k ) + repl_iy
    phi%rot( r_level, i, j, k + c ) = phi%rot( r_level, i, j, k + c ) + repl_iz

    ! Exterior points
    select case( ext_flag )
    case( x_only )
       phi%rot( r_level, i - a, j, k ) = phi%rot( r_level, i - a, j, k ) + repl_ex

    case( y_only )
       phi%rot( r_level, i, j - b, k ) = phi%rot( r_level, i, j - b, k ) + repl_ey

    case( z_only )
       phi%rot( r_level, i, j, k - c ) = phi%rot( r_level, i, j, k - c ) + repl_ez

    case( xy_only )
       phi%rot( r_level, i - a, j, k ) = phi%rot( r_level, i - a, j, k ) + repl_ex
       phi%rot( r_level, i, j - b, k ) = phi%rot( r_level, i, j - b, k ) + repl_ey

    case( xz_only )
       phi%rot( r_level, i - a, j, k ) = phi%rot( r_level, i - a, j, k ) + repl_ex
       phi%rot( r_level, i, j, k - c ) = phi%rot( r_level, i, j, k - c ) + repl_ez

    case( yz_only )
       phi%rot( r_level, i, j - b, k ) = phi%rot( r_level, i, j - b, k ) + repl_ey
       phi%rot( r_level, i, j, k - c ) = phi%rot( r_level, i, j, k - c ) + repl_ez

    case( xyz )
       phi%rot( r_level, i - a, j, k ) = phi%rot( r_level, i - a, j, k ) + repl_ex
       phi%rot( r_level, i, j - b, k ) = phi%rot( r_level, i, j - b, k ) + repl_ey
       phi%rot( r_level, i, j, k - c ) = phi%rot( r_level, i, j, k - c ) + repl_ez

    case default
       !TODO: Consistent error handling
       write(*,*) "Error: Invalid value of ext_flag - ",&
            ext_flag
       stop

    end select

    return
  end subroutine rot_replenishment_single

  subroutine vib_replenishment_single( phi, repl, mapping, v_level, v_modes )

    use DistFunc
    use Remapping

    implicit none

    type(DistFuncType) :: phi

    type(MappingResultType), intent(in) :: mapping
    double precision, intent(in) :: repl

    integer, intent(in) :: v_level, v_modes

    double precision :: repl_o, repl_ix, repl_iy, repl_iz, repl_ex, repl_ey, repl_ez

    integer :: i, j, k, a, b, c, ext_flag

    if( v_modes .eq. 0 )return

    ! read in stencil locations
    i = mapping%stencil(1)
    j = mapping%stencil(2)
    k = mapping%stencil(3)
    a = mapping%stencil(4)
    b = mapping%stencil(5)
    c = mapping%stencil(6)
    ext_flag = mapping%stencil(7)

    ! read in replenishment values
    repl_o  = repl * mapping%mass(1)
    repl_ix = repl * mapping%mass(2)
    repl_iy = repl * mapping%mass(3)
    repl_iz = repl * mapping%mass(4)
    repl_ex = repl * mapping%mass(5)
    repl_ey = repl * mapping%mass(6)
    repl_ez = repl * mapping%mass(7)

    ! Origin and interior points
    phi%vib( v_level, i, j, k )     = phi%vib( v_level, i,j,k ) + repl_o
    phi%vib( v_level, i + a, j, k ) = phi%vib( v_level, i + a, j, k ) + repl_ix
    phi%vib( v_level, i, j + b, k ) = phi%vib( v_level, i, j + b, k ) + repl_iy
    phi%vib( v_level, i, j, k + c ) = phi%vib( v_level, i, j, k + c ) + repl_iz

    ! Exterior points
    select case( ext_flag )
    case( x_only )
       phi%vib( v_level, i - a, j, k ) = phi%vib( v_level, i - a, j, k ) + repl_ex

    case( y_only )
       phi%vib( v_level, i, j - b, k ) = phi%vib( v_level, i, j - b, k ) + repl_ey

    case( z_only )
       phi%vib( v_level, i, j, k - c ) = phi%vib( v_level, i, j, k - c ) + repl_ez

    case( xy_only )
       phi%vib( v_level, i - a, j, k ) = phi%vib( v_level, i - a, j, k ) + repl_ex
       phi%vib( v_level, i, j - b, k ) = phi%vib( v_level, i, j - b, k ) + repl_ey

    case( xz_only )
       phi%vib( v_level, i - a, j, k ) = phi%vib( v_level, i - a, j, k ) + repl_ex
       phi%vib( v_level, i, j, k - c ) = phi%vib( v_level, i, j, k - c ) + repl_ez

    case( yz_only )
       phi%vib( v_level, i, j - b, k ) = phi%vib( v_level, i, j - b, k ) + repl_ey
       phi%vib( v_level, i, j, k - c ) = phi%vib( v_level, i, j, k - c ) + repl_ez

    case( xyz )
       phi%vib( v_level, i - a, j, k ) = phi%vib( v_level, i - a, j, k ) + repl_ex
       phi%vib( v_level, i, j - b, k ) = phi%vib( v_level, i, j - b, k ) + repl_ey
       phi%vib( v_level, i, j, k - c ) = phi%vib( v_level, i, j, k - c ) + repl_ez

    case default
       !TODO: Consistent error handling
       write(*,*) "Error: Invalid value of ext_flag - ",&
            ext_flag
       stop

    end select

    return
  end subroutine vib_replenishment_single

  subroutine rot_replenishment_energy( phi, repl, mapping, fr, r_modes )

    use DistFunc
    use Remapping
    use TimeStepping
    use MathUtilities

    implicit none

    type(DistFuncType) :: phi

    type(MappingResultType), intent(in) :: mapping
    double precision, intent(in) :: repl

    double precision, intent(in) :: fr
    integer, intent(in) :: r_modes

    double precision :: repl_o, repl_ix, repl_iy, repl_iz, repl_ex, repl_ey, repl_ez

    integer :: i, j, k, a, b, c, ext_flag

    double precision :: pos_dens, neg_dens, neg_energy
    double precision :: e_r, er_o, er_ix, er_iy, er_iz, er_ex, er_ey, er_ez
    double precision :: energy
    logical :: flag_o, flag_ix, flag_iy, flag_iz, flag_ex, flag_ey, flag_ez

    double precision :: e_start, e_end

    if( r_modes .eq. 0 )return

    ! read in stencil locations
    i = mapping%stencil(1)
    j = mapping%stencil(2)
    k = mapping%stencil(3)
    a = mapping%stencil(4)
    b = mapping%stencil(5)
    c = mapping%stencil(6)
    ext_flag = mapping%stencil(7)

    ! read in replenishment values
    repl_o  = repl * mapping%mass(1)
    repl_ix = repl * mapping%mass(2)
    repl_iy = repl * mapping%mass(3)
    repl_iz = repl * mapping%mass(4)
    repl_ex = repl * mapping%mass(5)
    repl_ey = repl * mapping%mass(6)
    repl_ez = repl * mapping%mass(7)

    e_start = fr * repl
    e_end = zero

    ! Find negative flags
    flag_o =  .false.
    flag_ix = .false.
    flag_iy = .false.
    flag_iz = .false.
    flag_ex = .false.
    flag_ey = .false.
    flag_ez = .false.
    if( dsgn( repl_o ) .ne. dsgn( phi%value(i,j,k) ) ) flag_o = .true.
    if( dsgn( repl_ix ) .ne. dsgn( phi%value(i+a,j,k) ) ) flag_ix = .true.
    if( dsgn( repl_iy ) .ne. dsgn( phi%value(i,j+b,k) ) ) flag_iy = .true.
    if( dsgn( repl_iz ) .ne. dsgn( phi%value(i,j,k+c) ) ) flag_iz = .true.
    
    select case( ext_flag )
    case( x_only )
       if( dsgn( repl_ex ) .ne. dsgn( phi%value(i-a,j,k) ) ) flag_ex = .true.
    case( y_only )
       if( dsgn( repl_ey ) .ne. dsgn( phi%value(i,j-b,k) ) ) flag_ey = .true.
    case( z_only )
       if( dsgn( repl_ez ) .ne. dsgn( phi%value(i,j,k-c) ) ) flag_ez = .true.
    case( xy_only )
       if( dsgn( repl_ex ) .ne. dsgn( phi%value(i-a,j,k) ) ) flag_ex = .true.
       if( dsgn( repl_ey ) .ne. dsgn( phi%value(i,j-b,k) ) ) flag_ey = .true.
    case( xz_only )
       if( dsgn( repl_ex ) .ne. dsgn( phi%value(i-a,j,k) ) ) flag_ex = .true.
       if( dsgn( repl_ez ) .ne. dsgn( phi%value(i,j,k-c) ) ) flag_ez = .true.
    case( yz_only )
       if( dsgn( repl_ey ) .ne. dsgn( phi%value(i,j-b,k) ) ) flag_ey = .true.
       if( dsgn( repl_ez ) .ne. dsgn( phi%value(i,j,k-c) ) ) flag_ez = .true.
    case( xyz )
       if( dsgn( repl_ex ) .ne. dsgn( phi%value(i-a,j,k) ) ) flag_ex = .true.
       if( dsgn( repl_ey ) .ne. dsgn( phi%value(i,j-b,k) ) ) flag_ey = .true.
       if( dsgn( repl_ez ) .ne. dsgn( phi%value(i,j,k-c) ) ) flag_ez = .true.
    case default
    end select

    ! Perform negatives first
    pos_dens = zero
    neg_dens = zero
    neg_energy = zero

    if( flag_o )then
       if( abs( phi%value(i,j,k) ) .gt. double_tol )then
          neg_dens = neg_dens + repl_o
          er_o = phi%rot( 1, i,j,k ) / phi%value( i,j,k )
          neg_energy = neg_energy + repl_o * er_o
       else
          er_o = zero
       end if
    else
       pos_dens = pos_dens + repl_o
    end if

    if( flag_ix )then
       if( abs( phi%value(i,j,k) ) .gt. double_tol  )then
          neg_dens = neg_dens + repl_ix
          er_ix = phi%rot( 1, i+a,j,k ) / phi%value( i+a,j,k )
          neg_energy = neg_energy + repl_ix * er_ix
       else
          er_ix = zero
       end if
    else
       pos_dens = pos_dens + repl_ix
    end if

    if( flag_iy )then
       if( abs( phi%value(i,j+b,k) ) .gt. double_tol  )then
          neg_dens = neg_dens + repl_iy
          er_iy = phi%rot( 1, i,j+b,k ) / phi%value( i,j+b,k )
          neg_energy = neg_energy + repl_iy * er_iy
       else
          er_iy = zero
       end if
    else
       pos_dens = pos_dens + repl_iy
    end if

    if( flag_iz )then
       if( abs( phi%value(i,j,k+c) ) .gt. double_tol  )then
          neg_dens = neg_dens + repl_iz
          er_iz = phi%rot( 1, i,j,k+c ) / phi%value( i,j,k+c )
          neg_energy = neg_energy + repl_iz * er_iz
       else
          er_iz = zero
       end if
    else
       pos_dens = pos_dens + repl_iz
    end if

    if( flag_ex )then
       if( abs( phi%value(i-a,j,k) ) .gt. double_tol  )then
          neg_dens = neg_dens + repl_ex
          er_ex = phi%rot( 1, i-a,j,k ) / phi%value( i-a,j,k )
          neg_energy = neg_energy + repl_ex * er_ex
       else
          er_ex = zero
       end if
    else
       pos_dens = pos_dens + repl_ex
    end if

    if( flag_ey )then
       if( abs( phi%value(i,j-b,k) ) .gt. double_tol  )then
          neg_dens = neg_dens + repl_ey
          er_ey = phi%rot( 1, i,j-b,k ) / phi%value( i,j-b,k )
          neg_energy = neg_energy + repl_ey * er_ey
       else
          er_ey = zero
       end if
    else
       pos_dens = pos_dens + repl_ey
    end if

    if( flag_ez )then
       if( abs( phi%value(i,j,k-c) ) .gt. double_tol  )then
          neg_dens = neg_dens + repl_ez
          er_ez = phi%rot( 1, i,j,k-c ) / phi%value( i,j,k-c )
          neg_energy = neg_energy + repl_ez * er_ez
       else
          er_ez = zero
       end if
    else
       pos_dens = pos_dens + repl_ez
    end if
    if( abs( pos_dens ) .gt. double_tol )then
       e_r = ( fr * repl - neg_energy ) / pos_dens
    else
       write(*,*)"losing energy"
       write(*,*)repl, fr, fr * repl
       write(*,*)repl_o, repl_ix, repl_iy, repl_iz, repl_ex, repl_ey, repl_ez
       write(*,*)phi%value(i,j,k), phi%value(i+a,j,k), phi%value(i,j+b,k), phi%value(i,j,k+c), &
            phi%value(i-a,j,k), phi%value(i,j-b,k), phi%value(i,j,k-c)
       write(*,*)neg_dens, neg_energy
       write(*,*)i,j,k,a,b,c
       pos_dens = repl_o
       neg_energy = neg_energy - repl_o * er_o
       e_r = ( fr * repl - neg_energy ) / pos_dens
       flag_o = .false.
!!$       stop
    end if

    ! Origin and interior points
    energy = e_r
    if( flag_o ) energy = er_o
    phi%rot( :, i, j, k ) = phi%rot( :, i,j,k ) + repl_o * energy
    e_end = e_end + repl_o * energy

    energy = e_r
    if( flag_ix ) energy = er_ix
    phi%rot( :, i + a, j, k ) = phi%rot( :, i + a, j, k ) + repl_ix * energy
    e_end = e_end + repl_ix * energy
    
    energy = e_r
    if( flag_iy ) energy = er_iy
    phi%rot( :, i, j + b, k ) = phi%rot( :, i, j + b, k ) + repl_iy * energy
    e_end = e_end + repl_iy * energy

    energy = e_r
    if( flag_iz ) energy = er_iz
    phi%rot( :, i, j, k + c ) = phi%rot( :, i, j, k + c ) + repl_iz * energy
    e_end = e_end + repl_iz * energy

    ! Exterior points
    select case( ext_flag )
    case( x_only )
       energy = e_r
       if( flag_ex ) energy = er_ex
       phi%rot( :, i - a, j, k ) = phi%rot( :, i - a, j, k ) + repl_ex * energy
       e_end = e_end + repl_ex * energy

    case( y_only )
       energy = e_r
       if( flag_ey ) energy = er_ey
       phi%rot( :, i, j - b, k ) = phi%rot( :, i, j - b, k ) + repl_ey * energy
       e_end = e_end + repl_ey * energy

    case( z_only )
       energy = e_r
       if( flag_ez ) energy = er_ez
       phi%rot( :, i, j, k - c ) = phi%rot( :, i, j, k - c ) + repl_ez * energy
       e_end = e_end + repl_ez * energy

    case( xy_only )
       energy = e_r
       if( flag_ex ) energy = er_ex
       phi%rot( :, i - a, j, k ) = phi%rot( :, i - a, j, k ) + repl_ex * energy
       e_end = e_end + repl_ex * energy

       energy = e_r
       if( flag_ey ) energy = er_ey
       phi%rot( :, i, j - b, k ) = phi%rot( :, i, j - b, k ) + repl_ey * energy
       e_end = e_end + repl_ey * energy

    case( xz_only )
       energy = e_r
       if( flag_ex ) energy = er_ex
       phi%rot( :, i - a, j, k ) = phi%rot( :, i - a, j, k ) + repl_ex * energy
       e_end = e_end + repl_ex * energy

       energy = e_r
       if( flag_ez ) energy = er_ez
       phi%rot( :, i, j, k - c ) = phi%rot( :, i, j, k - c ) + repl_ez * energy
       e_end = e_end + repl_ez * energy

    case( yz_only )
       energy = e_r
       if( flag_ey ) energy = er_ey
       phi%rot( :, i, j - b, k ) = phi%rot( :, i, j - b, k ) + repl_ey * energy
       e_end = e_end + repl_ey * energy

       energy = e_r
       if( flag_ez ) energy = er_ez
       phi%rot( :, i, j, k - c ) = phi%rot( :, i, j, k - c ) + repl_ez * energy
       e_end = e_end + repl_ez * energy

    case( xyz )
       energy = e_r
       if( flag_ex ) energy = er_ex
       phi%rot( :, i - a, j, k ) = phi%rot( :, i - a, j, k ) + repl_ex * energy
       e_end = e_end + repl_ex * energy

       energy = e_r
       if( flag_ey ) energy = er_ey
       phi%rot( :, i, j - b, k ) = phi%rot( :, i, j - b, k ) + repl_ey * energy
       e_end = e_end + repl_ey * energy

       energy = e_r
       if( flag_ez ) energy = er_ez
       phi%rot( :, i, j, k - c ) = phi%rot( :, i, j, k - c ) + repl_ez * energy
       e_end = e_end + repl_ez * energy


    case default
       !TODO: Consistent error handling
       write(*,*) "Error: Invalid value of ext_flag - ",&
            ext_flag
       stop

    end select

!!$    write(*,*)e_start, e_end

    return
  end subroutine rot_replenishment_energy

  subroutine rot_replenishment_energy2( phi, mapping, energy, r_modes )

    use DistFunc
    use Remapping
    use TimeStepping
    use MathUtilities

    implicit none

    type(DistFuncType) :: phi

    type(MappingResultType), intent(in) :: mapping
    double precision, intent(in) :: energy
    integer, intent(in) :: r_modes

    double precision :: repl_o, repl_ix, repl_iy, repl_iz, repl_ex, repl_ey, repl_ez
    integer :: i, j, k, a, b, c, ext_flag

    if( r_modes .eq. 0 )return

    ! read in stencil locations
    i = mapping%stencil(1)
    j = mapping%stencil(2)
    k = mapping%stencil(3)
    a = mapping%stencil(4)
    b = mapping%stencil(5)
    c = mapping%stencil(6)
    ext_flag = mapping%stencil(7)

    ! read in replenishment values
    repl_o  = energy * mapping%mass(1)
    repl_ix = energy * mapping%mass(2)
    repl_iy = energy * mapping%mass(3)
    repl_iz = energy * mapping%mass(4)
    repl_ex = energy * mapping%mass(5)
    repl_ey = energy * mapping%mass(6)
    repl_ez = energy * mapping%mass(7)

    phi%rot( 1, i, j, k ) = phi%rot( 1, i,j,k ) + repl_o
    phi%rot( 1, i + a, j, k ) = phi%rot( 1, i + a, j, k ) + repl_ix
    phi%rot( 1, i, j + b, k ) = phi%rot( 1, i, j + b, k ) + repl_iy
    phi%rot( 1, i, j, k + c ) = phi%rot( 1, i, j, k + c ) + repl_iz

    ! Exterior points
    select case( ext_flag )
    case( x_only )
       phi%rot( 1, i - a, j, k ) = phi%rot( 1, i - a, j, k ) + repl_ex

    case( y_only )
       phi%rot( 1, i, j - b, k ) = phi%rot( 1, i, j - b, k ) + repl_ey

    case( z_only )
       phi%rot( 1, i, j, k - c ) = phi%rot( 1, i, j, k - c ) + repl_ez

    case( xy_only )
       phi%rot( 1, i - a, j, k ) = phi%rot( 1, i - a, j, k ) + repl_ex
       phi%rot( 1, i, j - b, k ) = phi%rot( 1, i, j - b, k ) + repl_ey

    case( xz_only )
       phi%rot( 1, i - a, j, k ) = phi%rot( 1, i - a, j, k ) + repl_ex
       phi%rot( 1, i, j, k - c ) = phi%rot( 1, i, j, k - c ) + repl_ez

    case( yz_only )
       phi%rot( 1, i, j - b, k ) = phi%rot( 1, i, j - b, k ) + repl_ey
       phi%rot( 1, i, j, k - c ) = phi%rot( 1, i, j, k - c ) + repl_ez

    case( xyz )
       phi%rot( 1, i - a, j, k ) = phi%rot( 1, i - a, j, k ) + repl_ex
       phi%rot( 1, i, j - b, k ) = phi%rot( 1, i, j - b, k ) + repl_ey
       phi%rot( 1, i, j, k - c ) = phi%rot( 1, i, j, k - c ) + repl_ez

    case default
       !TODO: Consistent error handling
       write(*,*) "Error: Invalid value of ext_flag - ",&
            ext_flag
       stop

    end select

    return
  end subroutine rot_replenishment_energy2

end module ReplenishingCollisions
