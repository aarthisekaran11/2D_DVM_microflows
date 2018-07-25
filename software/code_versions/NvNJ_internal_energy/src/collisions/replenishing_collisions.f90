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

  subroutine replenish_collision( phi, depletion, fr1, fr2, fv1, fv2, frv1, frv2, &
       level, i, j, k, g, n, m, vel_grid, m_red, coll_type, exchange_mol, molecule )

    use DistFunc
    use SpeciesAndReferenceData
    use VelocityGrid
    use ReMapping

    implicit none

    type(DistFuncType), dimension(:) :: phi

    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    type(MoleculeType), dimension(:), intent(in) :: molecule
    double precision, dimension(:,:), intent(in) :: frv1, frv2
    double precision, dimension(:), intent(in) :: fr1, fr2, fv1, fv2
    double precision, intent(in) :: g, depletion, m_red
    integer, dimension(2), intent(in) :: i, j, k, level
    integer, intent(in) :: coll_type, n, m, exchange_mol

    double precision, dimension(2) :: x, y, z
    double precision, dimension(2) :: mass, omega, alpha, theta_r, theta_v
    double precision, dimension(3) :: v_com, xi, eta
    double precision :: delta_eps, g_prime, g_ratio
    double precision :: omega_AB, alpha_AB
    double precision :: repl_frac, coll_repl

    integer, dimension(2) :: r_modes, v_modes, r_levels, v_levels, level_p, mol_type
    integer :: repl

    type(MappingResultType), dimension(2) :: mapping

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

    ! Calculate replenishment variables
    repl_frac = depletion/dble( num_repl_locations )

    call compute_center_of_mass_velocity( x, y, z, mass, v_com )

    do repl = 1, num_repl_locations

       select case( coll_type )
       case( elastic )
          delta_eps = zero
          g_prime = g
          g_ratio = one

          call find_post_collision_velocities( xi, eta, x, y, z, v_com, g_prime, g_ratio, alpha_AB )

          call perform_remapping( mapping(1), xi(1), xi(2), xi(3), repl_frac, vel_grid(n) )
          call perform_remapping( mapping(2), eta(1), eta(2), eta(3), repl_frac, vel_grid(m) )

          call apply_phi_replenishment( phi(n), mapping(1) )
          call apply_phi_replenishment( phi(m), mapping(2) )

          call energy_replenishment_full( phi(n), mapping(1), frv1, r_modes(1), v_modes(1) )
          call energy_replenishment_full( phi(m), mapping(2), frv2, r_modes(2), v_modes(2) )

       case( rot_trans )
          if( exchange_mol .eq. A )then
             call compute_rot_level_change( level_p(1), delta_eps, level(1), &
                  g, omega_AB, m_red, r_levels(1), molecule(n), phi(n)%rot_level )
          else
             call compute_rot_level_change( level_p(2), delta_eps, level(2), &
                  g, omega_AB, m_red, r_levels(2), molecule(m), phi(m)%rot_level )
          end if

          call calculate_g_prime( g_prime, g, delta_eps, m_red )
          g_ratio = g_prime/g

          call find_post_collision_velocities( xi, eta, x, y, z, v_com, g_prime, g_ratio, alpha_AB )

          call perform_remapping( mapping(1), xi(1), xi(2), xi(3), repl_frac, vel_grid(n) )
          call perform_remapping( mapping(2), eta(1), eta(2), eta(3), repl_frac, vel_grid(m) )
          
          call apply_phi_replenishment( phi(n), mapping(1) )
          call apply_phi_replenishment( phi(m), mapping(2) )

          if( exchange_mol .eq. A )then
             call rot_replenishment_single( phi(n), mapping(1), level_p(1), fv1, r_modes(1), v_modes(2) )
             call energy_replenishment_full( phi(m), mapping(2), frv2, r_modes(2), v_modes(2) )
          else
             call energy_replenishment_full( phi(n), mapping(1), frv1, r_modes(1), v_modes(1) )
             call rot_replenishment_single( phi(m), mapping(2), level_p(2), fv2, r_modes(2), v_modes(2) )
          end if

       case( vib_trans )
          if( exchange_mol .eq. A )then
             call compute_vib_level_change( level_p(1), delta_eps, level(1), &
                  g, omega_AB, m_red, v_levels(1), molecule(n), phi(n)%vib_level )
             level_p(2) = level(2)
          else
             call compute_vib_level_change( level_p(2), delta_eps, level(2), &
                  g, omega_AB, m_red, v_levels(2), molecule(m), phi(m)%vib_level )
             level_p(1) = level(1)
          end if

          call calculate_g_prime( g_prime, g, delta_eps, m_red )
          g_ratio = g_prime/g

          call find_post_collision_velocities( xi, eta, x, y, z, v_com, g_prime, g_ratio, alpha_AB )

          call perform_remapping( mapping(1), xi(1), xi(2), xi(3), repl_frac, vel_grid(n) )
          call perform_remapping( mapping(2), eta(1), eta(2), eta(3), repl_frac, vel_grid(m) )

          call apply_phi_replenishment( phi(n), mapping(1) )
          call apply_phi_replenishment( phi(m), mapping(2) )
          
          if( exchange_mol .eq. A )then
             call vib_replenishment_single( phi(n), mapping(1), level_p(1), fr1, r_modes(1), v_modes(1) )
             call energy_replenishment_full( phi(m), mapping(2), frv2, r_modes(2), v_modes(2) )
          else
             call energy_replenishment_full( phi(n), mapping(1), frv1, r_modes(1), v_modes(1) )
             call vib_replenishment_single( phi(m), mapping(2), level_p(2), fr2, r_modes(2), v_modes(2) )
          end if

       case default

       end select

    end do

    return
  end subroutine replenish_collision

  subroutine apply_phi_replenishment( phi, mapping )

    use DistFunc
    use Remapping
    use TimeStepping

    implicit none

    type(DistFuncType) :: phi

    type(MappingResultType), intent(in) :: mapping

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
    repl_o  = mapping%mass(1)
    repl_ix = mapping%mass(2)
    repl_iy = mapping%mass(3)
    repl_iz = mapping%mass(4)
    repl_ex = mapping%mass(5)
    repl_ey = mapping%mass(6)
    repl_ez = mapping%mass(7)


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

  subroutine energy_replenishment_full( phi, mapping, frv, r_modes, v_modes )

    use DistFunc
    use Remapping
    use TimeStepping

    implicit none

    type(DistFuncType) :: phi

    type(MappingResultType), intent(in) :: mapping

    double precision, dimension(:,:), intent(in) :: frv
    integer, intent(in) :: r_modes, v_modes

    double precision :: repl_o, repl_ix, repl_iy, repl_iz, repl_ex, repl_ey, repl_ez

    integer :: i, j, k, a, b, c, ext_flag

    if( r_modes .eq. 0 .and. v_modes .eq. 0 )return

    ! read in stencil locations
    i = mapping%stencil(1)
    j = mapping%stencil(2)
    k = mapping%stencil(3)
    a = mapping%stencil(4)
    b = mapping%stencil(5)
    c = mapping%stencil(6)
    ext_flag = mapping%stencil(7)

    ! read in replenishment values
    repl_o  = mapping%mass(1)
    repl_ix = mapping%mass(2)
    repl_iy = mapping%mass(3)
    repl_iz = mapping%mass(4)
    repl_ex = mapping%mass(5)
    repl_ey = mapping%mass(6)
    repl_ez = mapping%mass(7)

    ! Origin and interior points
    phi%int_energy( :, :, i, j, k )     = phi%int_energy( :, :, i,j,k ) + repl_o*frv
    phi%int_energy( :, :, i + a, j, k ) = phi%int_energy( :, :, i + a, j, k ) + repl_ix*frv
    phi%int_energy( :, :, i, j + b, k ) = phi%int_energy( :, :, i, j + b, k ) + repl_iy*frv
    phi%int_energy( :, :, i, j, k + c ) = phi%int_energy( :, :, i, j, k + c ) + repl_iz*frv


    ! Exterior points
    select case( ext_flag )
    case( x_only )
       phi%int_energy( :, :, i - a, j, k ) = phi%int_energy( :, :, i - a, j, k ) + repl_ex*frv

    case( y_only )
       phi%int_energy( :, :, i, j - b, k ) = phi%int_energy( :, :, i, j - b, k ) + repl_ey*frv

    case( z_only )
       phi%int_energy( :, :, i, j, k - c ) = phi%int_energy( :, :, i, j, k - c ) + repl_ez*frv

    case( xy_only )
       phi%int_energy( :, :, i - a, j, k ) = phi%int_energy( :, :, i - a, j, k ) + repl_ex*frv
       phi%int_energy( :, :, i, j - b, k ) = phi%int_energy( :, :, i, j - b, k ) + repl_ey*frv

    case( xz_only )
       phi%int_energy( :, :, i - a, j, k ) = phi%int_energy( :, :, i - a, j, k ) + repl_ex*frv
       phi%int_energy( :, :, i, j, k - c ) = phi%int_energy( :, :, i, j, k - c ) + repl_ez*frv

    case( yz_only )
       phi%int_energy( :, :, i, j - b, k ) = phi%int_energy( :, :, i, j - b, k ) + repl_ey*frv
       phi%int_energy( :, :, i, j, k - c ) = phi%int_energy( :, :, i, j, k - c ) + repl_ez*frv

    case( xyz )
       phi%int_energy( :, :, i - a, j, k ) = phi%int_energy( :, :, i - a, j, k ) + repl_ex*frv
       phi%int_energy( :, :, i, j - b, k ) = phi%int_energy( :, :, i, j - b, k ) + repl_ey*frv
       phi%int_energy( :, :, i, j, k - c ) = phi%int_energy( :, :, i, j, k - c ) + repl_ez*frv


    case default
       !TODO: Consistent error handling
       write(*,*) "Error: Invalid value of ext_flag - ",&
            ext_flag
       stop

    end select

    return
  end subroutine energy_replenishment_full

  subroutine rot_replenishment_single( phi, mapping, r_level, fv, r_modes, v_modes )

    use DistFunc
    use Remapping

    implicit none

    type(DistFuncType) :: phi

    type(MappingResultType), intent(in) :: mapping
    double precision, dimension(:), intent(in) :: fv
    integer, intent(in) :: r_level, r_modes, v_modes

    double precision :: repl_o, repl_ix, repl_iy, repl_iz, repl_ex, repl_ey, repl_ez

    integer :: i, j, k, a, b, c, ext_flag

    if( r_modes .eq. 0 .and. v_modes .eq. 0 )return

    ! read in stencil locations
    i = mapping%stencil(1)
    j = mapping%stencil(2)
    k = mapping%stencil(3)
    a = mapping%stencil(4)
    b = mapping%stencil(5)
    c = mapping%stencil(6)
    ext_flag = mapping%stencil(7)

    ! read in replenishment values
    repl_o  = mapping%mass(1)
    repl_ix = mapping%mass(2)
    repl_iy = mapping%mass(3)
    repl_iz = mapping%mass(4)
    repl_ex = mapping%mass(5)
    repl_ey = mapping%mass(6)
    repl_ez = mapping%mass(7)

    ! Origin and interior points
    phi%int_energy( r_level, :, i, j, k )     = phi%int_energy( r_level, :, i,j,k ) + repl_o*fv
    phi%int_energy( r_level, :, i + a, j, k ) = phi%int_energy( r_level, :, i + a, j, k ) + repl_ix*fv
    phi%int_energy( r_level, :, i, j + b, k ) = phi%int_energy( r_level, :, i, j + b, k ) + repl_iy*fv
    phi%int_energy( r_level, :, i, j, k + c ) = phi%int_energy( r_level, :, i, j, k + c ) + repl_iz*fv

    ! Exterior points
    select case( ext_flag )
    case( x_only )
       phi%int_energy( r_level, :, i - a, j, k ) = phi%int_energy( r_level, :, i - a, j, k ) + repl_ex*fv

    case( y_only )
       phi%int_energy( r_level, :, i, j - b, k ) = phi%int_energy( r_level, :, i, j - b, k ) + repl_ey*fv

    case( z_only )
       phi%int_energy( r_level, :, i, j, k - c ) = phi%int_energy( r_level, :, i, j, k - c ) + repl_ez*fv

    case( xy_only )
       phi%int_energy( r_level, :, i - a, j, k ) = phi%int_energy( r_level, :, i - a, j, k ) + repl_ex*fv
       phi%int_energy( r_level, :, i, j - b, k ) = phi%int_energy( r_level, :, i, j - b, k ) + repl_ey*fv

    case( xz_only )
       phi%int_energy( r_level, :, i - a, j, k ) = phi%int_energy( r_level, :, i - a, j, k ) + repl_ex*fv
       phi%int_energy( r_level, :, i, j, k - c ) = phi%int_energy( r_level, :, i, j, k - c ) + repl_ez*fv

    case( yz_only )
       phi%int_energy( r_level, :, i, j - b, k ) = phi%int_energy( r_level, :, i, j - b, k ) + repl_ey*fv
       phi%int_energy( r_level, :, i, j, k - c ) = phi%int_energy( r_level, :, i, j, k - c ) + repl_ez*fv

    case( xyz )
       phi%int_energy( r_level, :, i - a, j, k ) = phi%int_energy( r_level, :, i - a, j, k ) + repl_ex*fv
       phi%int_energy( r_level, :, i, j - b, k ) = phi%int_energy( r_level, :, i, j - b, k ) + repl_ey*fv
       phi%int_energy( r_level, :, i, j, k - c ) = phi%int_energy( r_level, :, i, j, k - c ) + repl_ez*fv

    case default
       !TODO: Consistent error handling
       write(*,*) "Error: Invalid value of ext_flag - ",&
            ext_flag
       stop

    end select

    return
  end subroutine rot_replenishment_single

  subroutine vib_replenishment_single( phi, mapping, v_level, fr, r_modes, v_modes )

    use DistFunc
    use Remapping

    implicit none

    type(DistFuncType) :: phi

    type(MappingResultType), intent(in) :: mapping
    double precision, dimension(:), intent(in) :: fr
    integer, intent(in) :: v_level, r_modes, v_modes

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
    repl_o  = mapping%mass(1)
    repl_ix = mapping%mass(2)
    repl_iy = mapping%mass(3)
    repl_iz = mapping%mass(4)
    repl_ex = mapping%mass(5)
    repl_ey = mapping%mass(6)
    repl_ez = mapping%mass(7)

    ! Origin and interior points
    phi%int_energy( :, v_level, i, j, k )     = phi%int_energy( :, v_level, i,j,k ) + repl_o*fr
    phi%int_energy( :, v_level, i + a, j, k ) = phi%int_energy( :, v_level, i + a, j, k ) + repl_ix*fr
    phi%int_energy( :, v_level, i, j + b, k ) = phi%int_energy( :, v_level, i, j + b, k ) + repl_iy*fr
    phi%int_energy( :, v_level, i, j, k + c ) = phi%int_energy( :, v_level, i, j, k + c ) + repl_iz*fr

    ! Exterior points
    select case( ext_flag )
    case( x_only )
       phi%int_energy( :, v_level, i - a, j, k ) = phi%int_energy( :, v_level, i - a, j, k ) + repl_ex*fr

    case( y_only )
       phi%int_energy( :, v_level, i, j - b, k ) = phi%int_energy( :, v_level, i, j - b, k ) + repl_ey*fr

    case( z_only )
       phi%int_energy( :, v_level, i, j, k - c ) = phi%int_energy( :, v_level, i, j, k - c ) + repl_ez*fr

    case( xy_only )
       phi%int_energy( :, v_level, i - a, j, k ) = phi%int_energy( :, v_level, i - a, j, k ) + repl_ex*fr
       phi%int_energy( :, v_level, i, j - b, k ) = phi%int_energy( :, v_level, i, j - b, k ) + repl_ey*fr

    case( xz_only )
       phi%int_energy( :, v_level, i - a, j, k ) = phi%int_energy( :, v_level, i - a, j, k ) + repl_ex*fr
       phi%int_energy( :, v_level, i, j, k - c ) = phi%int_energy( :, v_level, i, j, k - c ) + repl_ez*fr

    case( yz_only )
       phi%int_energy( :, v_level, i, j - b, k ) = phi%int_energy( :, v_level, i, j - b, k ) + repl_ey*fr
       phi%int_energy( :, v_level, i, j, k - c ) = phi%int_energy( :, v_level, i, j, k - c ) + repl_ez*fr

    case( xyz )
       phi%int_energy( :, v_level, i - a, j, k ) = phi%int_energy( :, v_level, i - a, j, k ) + repl_ex*fr
       phi%int_energy( :, v_level, i, j - b, k ) = phi%int_energy( :, v_level, i, j - b, k ) + repl_ey*fr
       phi%int_energy( :, v_level, i, j, k - c ) = phi%int_energy( :, v_level, i, j, k - c ) + repl_ez*fr

    case default
       !TODO: Consistent error handling
       write(*,*) "Error: Invalid value of ext_flag - ",&
            ext_flag
       stop

    end select

    return
  end subroutine vib_replenishment_single

end module ReplenishingCollisions
