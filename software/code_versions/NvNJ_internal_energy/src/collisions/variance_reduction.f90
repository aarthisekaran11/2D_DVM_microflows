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
module VarianceReduction

  use ErrorCheck
  use Constants
  use MathUtilities
  use CollisionUtilities
  use SpeciesAndReferenceData

  implicit none

  private

  type(CumulativeDF), allocatable, dimension(:,:,:) :: psi_eq
  type(CumulativeDF), allocatable, dimension(:,:) :: psi_neq

  integer, parameter :: elastic   = 0
  integer, parameter :: rot_trans = 1
  integer, parameter :: vib_trans = 2

  integer, parameter :: A = 1
  integer, parameter :: B = 2

  integer, parameter :: AA = 1
  integer, parameter :: AB = 2

  public :: create_variance_reduction
  public :: destroy_variance_reduction
  public :: compute_psi_vr
  public :: variance_reduction_kernel

contains

  subroutine create_variance_reduction( vel_grid, phi, molecule )

    use VelocityGrid
    use DistFunc 
    use PhysicalGrid

    implicit none

    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    type(MoleculeType), dimension(:), intent(in) :: molecule

    integer :: num_points, nspace, grid_ref
    integer :: r_modes, v_modes, r_levels, v_levels
    integer :: n, same, s
    integer :: status

    call get_nspace( nspace )

    allocate( psi_eq( 1:nspace, 1:num_species, 1:2 ), STAT=status )
    call allocate_error_check( status, "psi_eq(VR)" )

    allocate( psi_neq( 1:num_species, 1:2 ), STAT=status )
    call allocate_error_check( status, "psi_neq(VR)" )

    do same = 1, 2
       do s = 1, num_species

          call get_spatial_reference( 1, grid_ref )
          num_points = vel_grid( s, grid_ref )%num_points

          allocate( psi_neq(s,same)%cumul_df( 0:num_points ), STAT=status )
          call allocate_error_check( status, "psi_neq%cumul_df(VR)" )

          allocate( psi_neq(s,same)%sign( 0:num_points ), STAT=status )
          call allocate_error_check( status, "psi_neq%sign(VR)" )
          
          r_modes = molecule(s)%rot_modes
          v_modes = molecule(s)%vib_modes

          allocate( psi_neq(s,same)%cumul_rot_df( 0:num_points ), STAT=status )
          call allocate_error_check( status, "psi_neq%cumul_rot_df(VR)" )

          allocate( psi_neq(s,same)%cumul_vib_df( 0:num_points ), STAT=status )
          call allocate_error_check( status, "psi_neq%cumul_vib_df(VR)" )

          do n = 1, nspace

             call get_spatial_reference( n, grid_ref )
             num_points = vel_grid( s, grid_ref )%num_points

             allocate( psi_eq(n,s,same)%cumul_df( 0:num_points ), STAT=status )
             call allocate_error_check( status, "psi_eq%cumul_df(VR)" )

             allocate( psi_eq(n,s,same)%sign( 0:num_points ), STAT=status )
             call allocate_error_check( status, "psi_eq%sign(VR)" )     

             if( r_modes .gt. 0 )then
                r_levels = phi(s,n)%num_rot_levels
                allocate( psi_eq(n,s,same)%rot( 1:r_levels ), STAT=status )
                call allocate_error_check( status, "psi_eq%rot(VR)" )
             end if

             if( v_modes .gt. 0 )then
                v_levels = phi(s,n)%num_vib_levels
                allocate( psi_eq(n,s,same)%vib( 1:v_levels ), STAT=status )
                call allocate_error_check( status, "psi_eq%vib(VR)" )
             end if

          end do
       end do
    end do

    return
  end subroutine create_variance_reduction

  subroutine destroy_variance_reduction( molecule )

    use DistFunc
    use PhysicalGrid

    implicit none

    type(MoleculeType), dimension(:), intent(in) :: molecule  

    integer :: nspace, n, s, same
    integer :: r_modes, v_modes
    integer :: status

    call get_nspace( nspace )

    do same = 1, 2
       do s = 1, num_species

          deallocate( psi_neq(s,same)%cumul_df, STAT=status )
          call deallocate_error_check( status, "psi_neq%cumul_df(VR)" )

          deallocate( psi_neq(s,same)%sign, STAT=status )
          call deallocate_error_check( status, "psi_neq%sign(VR)" )

          r_modes = molecule(s)%rot_modes
          v_modes = molecule(s)%vib_modes

          deallocate( psi_neq(s,same)%cumul_rot_df, STAT=status )
          call deallocate_error_check( status, "psi_neq%cumul_rot_df(VR)" )

          deallocate( psi_neq(s,same)%cumul_vib_df, STAT=status )
          call deallocate_error_check( status, "psi_neq%cumul_vib_df(VR)" )

          do n = 1, nspace

             deallocate( psi_eq(n,s,same)%cumul_df, STAT=status )
             call deallocate_error_check( status, "psi_eq%cumul_df(VR)" )

             deallocate( psi_eq(n,s,same)%sign, STAT=status )
             call deallocate_error_check( status, "psi_eq%sign(VR)" )

             if( r_modes .gt. 0 )then
                deallocate( psi_eq(n,s,same)%rot, STAT=status )
                call deallocate_error_check( status, "psi_eq%rot(VR)" )
             end if

             if( v_modes .gt. 0 )then
                deallocate( psi_eq(n,s,same)%vib, STAT=status )
                call deallocate_error_check( status, "psi_eq%vib(VR)" )
             end if

          end do
       end do
    end do

    deallocate( psi_eq, STAT=status )
    call deallocate_error_check( status, "psi_eq(VR)" )

    deallocate( psi_neq, STAT=status )
    call deallocate_error_check( status, "psi_neq(VR)" )

    return
  end subroutine destroy_variance_reduction

  subroutine compute_psi_vr( phi, molecule, vel_grid, properties, new_eq_flag, energy_flag, same, s, n )

    use DistFunc
    use VelocityGrid
    use PhysicalProperties

    implicit none

    type(DistFuncType), intent(in) :: phi
    type(MoleculeType), intent(in) :: molecule
    type(VelocityGridType), intent(in) :: vel_grid
    type(PropertiesType), intent(in) :: properties
    integer, intent(in) :: same, s, n, energy_flag
    logical, intent(in) :: new_eq_flag

    double precision, allocatable, dimension(:) :: temporary

    double precision :: x, y, z
    double precision :: beta_x, beta_y, beta_z, beta3
    double precision :: dens, u, v, w, temp, mass
    double precision :: maxwell_coeff, phi_eq, phi_neq
    double precision :: energy_level, fraction

    integer :: num_points_x, num_points_y, num_points
    integer :: l, i, j, k, p
    integer :: i_min, j_min, k_min
    integer :: r_modes, v_modes, r_levels, v_levels
    integer :: mol_type

    ! Internal energy modes
    r_modes = molecule%rot_modes
    v_modes = molecule%vib_modes

    ! Grid properties
    i_min = vel_grid%i_min
    j_min = vel_grid%j_min
    k_min = vel_grid%k_min

    num_points_x = vel_grid%num_points_x
    num_points_y = vel_grid%num_points_y
    num_points   = vel_grid%num_points

    select case( same )
    case( AA )
       ! Individual species properties for when like molecules collide
       mass = molecule%mass
       dens = properties%dens(s)
       u    = properties%x_vel(s)
       v    = properties%y_vel(s)
       w    = properties%z_vel(s)
       temp = properties%species_temp(s)
       if( temp .lt. double_tol ) temp = properties%temp(s)

    case( AB )
       ! Mixture properties for when unlike molecules collide
       mass = molecule%mass
       dens = properties%dens(s)
       u    = properties%mix_x_vel
       v    = properties%mix_y_vel
       w    = properties%mix_z_vel
       temp = properties%mix_temp

    case default
       write(*,*) "Error: Invalid value of AA or AB: ", same
       stop

    end select

    ! Set initial array values
    psi_eq(n,s,same)%cumul_df(0)    = zero
    psi_neq(s,same)%cumul_df(0)     = zero

    psi_neq(s,same)%cumul_rot_df(0) = zero
    psi_neq(s,same)%cumul_vib_df(0) = zero

    psi_eq(n,s,same)%sign(0)  = one
    psi_neq(s,same)%sign(0)   = one
    
    if( new_eq_flag .eqv. .true. )then
       ! Coefficient for calculating the Maxwellian
       maxwell_coeff = dens*sqrt( mass*mass*mass )/sqrt( pi*pi*pi*temp*temp*temp )

       ! Calculate the rotational equilibrium distribution function
       if( r_modes .gt. 0 )then
          r_levels  = phi%num_rot_levels
          
          allocate( temporary(1:r_levels) )
          call compute_rot_distribution( psi_eq(n,s,same)%rot(:), temporary, molecule, temp, r_levels, s )
          deallocate( temporary )

       end if
       
       ! Calculate the vibrational equilibrium distribution function
       if( v_modes .gt. 0 )then
          v_levels = phi%num_vib_levels
          
          allocate( temporary(1:v_levels) )
          call compute_vib_distribution( psi_eq(n,s,same)%vib(:), temporary, molecule, temp, v_levels, s )
          deallocate( temporary )
          
       end if
    end if

    ! Loop over all the velocity points and calculate the equilibrium and deviation distributions
    do l = 1, num_points

       ! Find x,y,z coordinates of the velocity location
       call global2local_map( l-1, i_min, j_min, k_min, num_points_x, num_points_y, i, j, k )

       ! Either calculate a new equilibrium value or read the equilibrium value from the previous time 
       ! step's equilibrium distribution
       if( new_eq_flag .eqv. .true. )then
          x = vel_grid%x(i)
          y = vel_grid%y(j)
          z = vel_grid%z(k)

          beta_x = vel_grid%beta_x(i)
          beta_y = vel_grid%beta_y(j)
          beta_z = vel_grid%beta_z(k)
          beta3  = beta_x*beta_y*beta_z

          call compute_maxwellian( x, y, z, maxwell_coeff, mass, u, v, w, temp, phi_eq )
          phi_eq = phi_eq*beta3

          psi_eq(n,s,same)%cumul_df(l) = psi_eq(n,s,same)%cumul_df(l-1) + phi_eq
          psi_eq(n,s,same)%sign(l) = one

       else
          phi_eq = psi_eq(n,s,same)%cumul_df(l) - psi_eq(n,s,same)%cumul_df(l-1)

       end if

       phi_neq = phi%value(i,j,k) - phi_eq

       psi_neq(s,same)%cumul_df(l) = psi_neq(s,same)%cumul_df(l-1) + abs( phi_neq )
       psi_neq(s,same)%sign(l)     = dsgn( phi_neq )
       
       ! The rotational deviation density at each velocity location is the sum of the absolute magnitude 
       ! of each rotational level
       if( r_modes .gt. 0 )then
          phi_neq = zero
          do p = 1, r_levels
             phi_neq = phi_neq + abs( phi%rot(p,i,j,k) - psi_eq(n,s,same)%rot(p)*phi_eq )
          end do

          psi_neq(s,same)%cumul_rot_df(l) = psi_neq(s,same)%cumul_rot_df(l-1) + phi_neq

       else
          psi_neq(s,same)%cumul_rot_df(l) = zero
          
       end if

       ! The vibrational deviation density at each velocity location is the sum of the absolute magnitude 
       ! of each vibrational level
       if( v_modes .gt. 0 )then
          phi_neq = zero
          do p = 1, v_levels
             phi_neq = phi_neq + abs( phi%vib(p,i,j,k) - psi_eq(n,s,same)%vib(p)*phi_eq )
          end do

          psi_neq(s,same)%cumul_vib_df(l) = psi_neq(s,same)%cumul_vib_df(l-1) + phi_neq

       else
          psi_neq(s,same)%cumul_vib_df(l) = zero

       end if

    end do

    psi_neq(s,same)%neq_dens = psi_neq(s,same)%cumul_df(num_points)
    psi_neq(s,same)%rot_neq_dens = psi_neq(s,same)%cumul_rot_df(num_points)
    psi_neq(s,same)%vib_neq_dens = psi_neq(s,same)%cumul_vib_df(num_points)

    return
  end subroutine compute_psi_vr

  subroutine variance_reduction_kernel( phi, tot_colls, n, m, ns, coln_rms, molecule, vel_grid, properties )

    use DistFunc
    use VelocityGrid
    use PhysicalProperties
    use ReplenishingCollisions
    use TimeStepping
    use Scaling

    implicit none

    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(PropertiesType), intent(in) :: properties
    double precision, dimension(:,:), intent(in) :: coln_rms
    integer, intent(in) :: n, m, ns

    type(DistFuncType), dimension(:) :: phi
    integer :: tot_colls

    logical :: rot_flag, vib_flag

    double precision, dimension(:), allocatable :: frA, frB, fvA, fvB
    double precision, dimension(:), allocatable :: rotA, rotB, vibA, vibB

    double precision, dimension(5) :: depl_frac, depl_sign
    double precision, dimension(2) :: mass, diam, omega
    double precision, dimension(2) :: dens, neq_dens, temp
    double precision :: m_red, omega_AB, vhs_exponent, sigma
    double precision :: dt, factor_coeff
    double precision :: factor, depletion
    double precision :: g, g_sigma
    double precision :: phiA, phiB, phiA_eq, phiB_eq
    double precision :: r_sumA, r_sumB, v_sumA, v_sumB

    integer, dimension(4) :: rl, vl
    integer, dimension(2) :: r_modes, v_modes, r_levels, v_levels
    integer, dimension(2) :: i, j, k
    integer :: same
    integer :: coll, num_colls
    integer :: glA, glB
    integer :: status

    ! Initialize flags for energy activation
    rot_flag = .false.
    vib_flag = .false.

    ! Species specific constants
    mass(1)  = molecule(n)%mass
    diam(1)  = molecule(n)%diam
    omega(1) = molecule(n)%omega

    mass(2)  = molecule(m)%mass
    diam(2)  = molecule(m)%diam
    omega(2) = molecule(m)%omega

    ! Species combination constants
    m_red = mass(1)*mass(2)/( mass(1) + mass(2) )

    omega_AB = one_half*( omega(1) + omega(2) )
    call calculate_cross_section( sigma, m_red, omega_AB, diam )

    vhs_exponent = two - two*omega_AB

    ! Time step
    call get_deltat( dt )

    if( n .eq. m )then
       factor_coeff = dt*one_half
       same = AA
    else
       factor_coeff = dt
       same = AB
    end if

    ! Equilibrium species properties
    dens(1) = properties%dens(n)
    temp(1) = properties%tr_temp(n)

    dens(2) = properties%dens(m)
    temp(2) = properties%tr_temp(m)

    ! Internal Structure
    r_modes(1) = molecule(n)%rot_modes
    v_modes(1) = molecule(n)%vib_modes

    r_modes(2) = molecule(m)%rot_modes
    v_modes(2) = molecule(m)%vib_modes

    ! Nonequilibrium density
    neq_dens(1) = psi_neq(n,same)%neq_dens
    neq_dens(2) = psi_neq(m,same)%neq_dens

    ! Set levels and allocate energy arrays
    if( r_modes(1) .gt. 0 )then
       rot_flag = .true.
       r_levels(1) = phi(n)%num_rot_levels
       allocate( frA( 1:r_levels(1) ), rotA( 1:r_levels(1) ), STAT=status )
       call allocate_error_check( status, "frA, rotA (VR)" )
    end if

    if( r_modes(2) .gt. 0 )then
       rot_flag = .true.
       r_levels(2) = phi(m)%num_rot_levels
       allocate( frB( 1:r_levels(2) ), rotB( 1:r_levels(2) ), STAT=status )
       call allocate_error_check( status, "frB, rotB (VR)" )
    end if

    if( v_modes(1) .gt. 0 )then
       vib_flag = .true.
       v_levels(1) = phi(n)%num_vib_levels
       allocate( fvA( 1:v_levels(1) ), vibA( 1:v_levels(1) ), STAT=status )
       call allocate_error_check( status, "fvA, vibA (VR)" )
    end if

    if( v_modes(2) .gt. 0 )then
       vib_flag = .true.
       v_levels(2) = phi(m)%num_vib_levels
       allocate( fvB( 1:v_levels(2) ), vibB( 1:v_levels(2) ), STAT=status )
       call allocate_error_check( status, "fvB, vibB (VR)" )
    end if

    ! Separate depletion into elastic and inelastic parts (elastic, r-t A, r-t B, v-t A, v-t B)
    call split_depletion( depl_frac, m_red, properties, molecule, n, m )

    !======================================================================================================
    ! First collision integral: phi_eq*phi_noneq
    !                           A-B, not B-A
    !======================================================================================================
!!$    ! Calculate the number of collisions
!!$    call compute_num_coll_pairs( num_colls, dt, dens(1)/dens(1), neq_dens(2)/dens(2), &
!!$         temp, m_red, coln_rms, vel_grid, n, m )
!!$
!!$    ! Number of collisions can't be zero
!!$    if( num_colls .eq. 0 ) num_colls = 1
!!$    tot_colls = tot_colls + num_colls
!!$
!!$    ! Calculate factor for depletion
!!$    factor  = factor_coeff*dens(1)*neq_dens(2)/dble(num_colls)
!!$
!!$    do coll = 1, num_colls
!!$
!!$       ! Pick velocities to collide
!!$       call pick_collision_partner( i(1), j(1), k(1), glA, psi_eq(ns,n,same)%cumul_df, vel_grid(n) )
!!$       call pick_collision_partner( i(2), j(2), k(2), glB, psi_neq(m,same)%cumul_df, vel_grid(m) )
!!$
!!$       ! Calculate relative velocity between chosen collision partners
!!$       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )
!!$
!!$       ! Check for self collision
!!$       if( g .lt. double_tol ) cycle
!!$
!!$       ! Calculate g*sigma term
!!$       g_sigma = sigma*g**vhs_exponent
!!$
!!$       ! Temporary arrays for phi
!!$       phiA = phi(n)%value( i(1), j(1), k(1) )
!!$       phiB = phi(m)%value( i(2), j(2), k(2) )
!!$
!!$       ! Elastic sign
!!$       depl_sign(1) = psi_eq(ns,n,same)%sign( glA )*psi_neq(m,same)%sign( glB )
!!$
!!$       ! Set internal energy arrays and pick levels for energy transfer collisions
!!$       if( r_modes(1) .gt. 0 )then
!!$          rotA = phi(n)%rot( :, i(1), j(1), k(1) )
!!$          call pick_level( rl(1), rotA, one, r_levels(1) )          
!!$          depl_sign(2) = dsgn( phiB )*dsgn( phi(n)%rot( rl(1), i(1), j(1), k(1) ) )
!!$          call energy_fraction_array( frA, r_sumA, rotA, r_levels(1), dsgn(phiA) )
!!$       end if
!!$
!!$       if( r_modes(2) .gt. 0 )then
!!$          rotB = phi(m)%rot( :, i(2), j(2), k(2) )
!!$          call pick_level( rl(2), rotB, one, r_levels(2) )
!!$          depl_sign(3) = dsgn( phiA )*dsgn( phi(m)%rot( rl(2), i(2), j(2), k(2) ) )
!!$          call energy_fraction_array( frB, r_sumB, rotB, r_levels(2), dsgn(phiB) )
!!$       end if
!!$
!!$       if( v_modes(1) .gt. 0 )then
!!$          vibA = phi(n)%vib( :, i(1), j(1), k(1) )
!!$          call pick_level( vl(1), vibA, one, v_levels(1) )
!!$          depl_sign(4) = dsgn( phiB )*dsgn( phi(n)%vib( vl(1), i(1), j(1), k(1) ) )
!!$          call energy_fraction_array( fvA, v_sumA, vibA, v_levels(1), dsgn(phiA) )
!!$       end if
!!$
!!$       if( v_modes(2) .gt. 0 )then
!!$          vibB = phi(m)%vib( :, i(2), j(2), k(2) )
!!$          call pick_level( vl(2), vibB, one, v_levels(2) )
!!$          depl_sign(5) = dsgn( phiA )*dsgn( phi(m)%vib( vl(2), i(2), j(2), k(2) ) )
!!$          call energy_fraction_array( fvB, v_sumB, vibB, v_levels(2), dsgn(phiB) )
!!$       end if
!!$
!!$       !**ELASTIC*********************************************************************************************
!!$       depletion = depl_sign(1)*depl_frac(1)*factor*g_sigma
!!$
!!$       call vr_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), 0, &
!!$            frA, fvA, r_modes(1), v_modes(1) )
!!$       call vr_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), 0, &
!!$            frB, fvB, r_modes(2), v_modes(2) )
!!$       call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
!!$            g, n, m, vel_grid, m_red, elastic, A, molecule )
!!$    end do
!!$
!!$    !**ROTATIONAL - TRANSLATIONAL**************************************************************************
!!$    !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
!!$
!!$    !======================================================================================================
!!$    ! Second collision integral: phi_noneq*phi_eq
!!$    !                            B-A, not A-B
!!$    !======================================================================================================
!!$    ! Calculate the number of collisions
!!$    call compute_num_coll_pairs( num_colls, dt, neq_dens(1)/dens(1), dens(2)/dens(2), &
!!$         temp, m_red, coln_rms, vel_grid, n, m )
!!$
!!$    ! Number of collisions can't be zero
!!$    if( num_colls .eq. 0 ) num_colls = 1
!!$    tot_colls = tot_colls + num_colls
!!$
!!$    ! Calculate factor for depletion
!!$    factor  = factor_coeff*neq_dens(1)*dens(2)/dble(num_colls)
!!$
!!$    do coll = 1, num_colls
!!$
!!$       ! Pick velocities to collide
!!$       call pick_collision_partner( i(1), j(1), k(1), glA, psi_neq(n,same)%cumul_df, vel_grid(n) )
!!$       call pick_collision_partner( i(2), j(2), k(2), glB, psi_eq(ns,m,same)%cumul_df, vel_grid(m) )
!!$
!!$       ! Calculate relative velocity between chosen collision partners
!!$       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )
!!$
!!$       ! Check for self collision
!!$       if( g .lt. double_tol ) cycle
!!$
!!$       ! Calculate g*sigma term
!!$       g_sigma = sigma*g**vhs_exponent
!!$
!!$       ! Temporary arrays for phi
!!$       phiA = phi(n)%value( i(1), j(1), k(1) )
!!$       phiB = phi(m)%value( i(2), j(2), k(2) )
!!$
!!$       ! Elastic sign
!!$       depl_sign(1) = psi_neq(n,same)%sign( glA )*psi_eq(ns,m,same)%sign( glB )
!!$
!!$       ! Set internal energy arrays and pick levels for energy transfer collisions
!!$       if( r_modes(1) .gt. 0 )then
!!$          rotA = phi(n)%rot( :, i(1), j(1), k(1) )
!!$          call pick_level( rl(1), rotA, one, r_levels(1) )          
!!$          depl_sign(2) = dsgn( phiB )*dsgn( phi(n)%rot( rl(1), i(1), j(1), k(1) ) )
!!$          call energy_fraction_array( frA, r_sumA, rotA, r_levels(1), dsgn(phiA) )
!!$       end if
!!$
!!$       if( r_modes(2) .gt. 0 )then
!!$          rotB = phi(m)%rot( :, i(2), j(2), k(2) )
!!$          call pick_level( rl(2), rotB, one, r_levels(2) )
!!$          depl_sign(3) = dsgn( phiA )*dsgn( phi(m)%rot( rl(2), i(2), j(2), k(2) ) )
!!$          call energy_fraction_array( frB, r_sumB, rotB, r_levels(2), dsgn(phiB) )
!!$       end if
!!$
!!$       if( v_modes(1) .gt. 0 )then
!!$          vibA = phi(n)%vib( :, i(1), j(1), k(1) )
!!$          call pick_level( vl(1), vibA, one, v_levels(1) )
!!$          depl_sign(4) = dsgn( phiB )*dsgn( phi(n)%vib( vl(1), i(1), j(1), k(1) ) )
!!$          call energy_fraction_array( fvA, v_sumA, vibA, v_levels(1), dsgn(phiA) )
!!$       end if
!!$
!!$       if( v_modes(2) .gt. 0 )then
!!$          vibB = phi(m)%vib( :, i(2), j(2), k(2) )
!!$          call pick_level( vl(2), vibB, one, v_levels(2) )
!!$          depl_sign(5) = dsgn( phiA )*dsgn( phi(m)%vib( vl(2), i(2), j(2), k(2) ) )
!!$          call energy_fraction_array( fvB, v_sumB, vibB, v_levels(2), dsgn(phiB) )
!!$       end if
!!$
!!$       !**ELASTIC*********************************************************************************************
!!$       depletion = depl_sign(1)*depl_frac(1)*factor*g_sigma
!!$
!!$       call vr_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), 0, &
!!$            frA, fvA, r_modes(1), v_modes(1) )
!!$       call vr_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), 0, &
!!$            frB, fvB, r_modes(2), v_modes(2) )
!!$       call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
!!$            g, n, m, vel_grid, m_red, elastic, A, molecule )
!!$    end do
!!$
!!$    !**ROTATIONAL - TRANSLATIONAL**************************************************************************
!!$    !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
!!$
!!$    !======================================================================================================
!!$    ! Third collision integral: phi_noneq*phi_noneq
!!$    !                           A-B and B-A
!!$    !======================================================================================================
!!$    ! Calculate the number of collisions
!!$    call compute_num_coll_pairs( num_colls, dt, neq_dens(1)/dens(1), neq_dens(2)/dens(2), &
!!$         temp, m_red, coln_rms, vel_grid, n, m )
!!$
!!$    ! Number of collisions can't be zero
!!$    if( num_colls .eq. 0 ) num_colls = 1
!!$    tot_colls = tot_colls + num_colls
!!$
!!$    ! Calculate factor for depletion
!!$    factor  = factor_coeff*neq_dens(1)*neq_dens(2)/dble(num_colls)
!!$
!!$    do coll = 1, num_colls
!!$
!!$       ! Pick velocities to collide
!!$       call pick_collision_partner( i(1), j(1), k(1), glA, psi_neq(n,same)%cumul_df, vel_grid(n) )
!!$       call pick_collision_partner( i(2), j(2), k(2), glB, psi_neq(m,same)%cumul_df, vel_grid(m) )
!!$
!!$       ! Calculate relative velocity between chosen collision partners
!!$       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )
!!$
!!$       ! Check for self collision
!!$       if( g .lt. double_tol ) cycle
!!$
!!$       ! Calculate g*sigma term
!!$       g_sigma = sigma*g**vhs_exponent
!!$
!!$       ! Temporary arrays for phi
!!$       phiA = phi(n)%value( i(1), j(1), k(1) )
!!$       phiB = phi(m)%value( i(2), j(2), k(2) )
!!$
!!$       ! Elastic sign
!!$       depl_sign(1) = psi_neq(n,same)%sign( glA )*psi_neq(m,same)%sign( glB )
!!$
!!$       ! Set internal energy arrays and pick levels for energy transfer collisions
!!$       if( r_modes(1) .gt. 0 )then
!!$          rotA = phi(n)%rot( :, i(1), j(1), k(1) )
!!$          call pick_level( rl(1), rotA, one, r_levels(1) )          
!!$          depl_sign(2) = dsgn( phiB )*dsgn( phi(n)%rot( rl(1), i(1), j(1), k(1) ) )
!!$          call energy_fraction_array( frA, r_sumA, rotA, r_levels(1), dsgn(phiA) )
!!$       end if
!!$
!!$       if( r_modes(2) .gt. 0 )then
!!$          rotB = phi(m)%rot( :, i(2), j(2), k(2) )
!!$          call pick_level( rl(2), rotB, one, r_levels(2) )
!!$          depl_sign(3) = dsgn( phiA )*dsgn( phi(m)%rot( rl(2), i(2), j(2), k(2) ) )
!!$          call energy_fraction_array( frB, r_sumB, rotB, r_levels(2), dsgn(phiB) )
!!$       end if
!!$
!!$       if( v_modes(1) .gt. 0 )then
!!$          vibA = phi(n)%vib( :, i(1), j(1), k(1) )
!!$          call pick_level( vl(1), vibA, one, v_levels(1) )
!!$          depl_sign(4) = dsgn( phiB )*dsgn( phi(n)%vib( vl(1), i(1), j(1), k(1) ) )
!!$          call energy_fraction_array( fvA, v_sumA, vibA, v_levels(1), dsgn(phiA) )
!!$       end if
!!$
!!$       if( v_modes(2) .gt. 0 )then
!!$          vibB = phi(m)%vib( :, i(2), j(2), k(2) )
!!$          call pick_level( vl(2), vibB, one, v_levels(2) )
!!$          depl_sign(5) = dsgn( phiA )*dsgn( phi(m)%vib( vl(2), i(2), j(2), k(2) ) )
!!$          call energy_fraction_array( fvB, v_sumB, vibB, v_levels(2), dsgn(phiB) )
!!$       end if
!!$
!!$       !**ELASTIC*********************************************************************************************
!!$       depletion = depl_sign(1)*depl_frac(1)*factor*g_sigma
!!$
!!$       call vr_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), 0, &
!!$            frA, fvA, r_modes(1), v_modes(1) )
!!$       call vr_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), 0, &
!!$            frB, fvB, r_modes(2), v_modes(2) )
!!$       call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
!!$            g, n, m, vel_grid, m_red, elastic, A, molecule )
!!$    end do
!!$
!!$    !**ROTATIONAL - TRANSLATIONAL**************************************************************************
!!$    !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
!!$
!!$    if( r_modes(1) .gt. 0 ) deallocate( frA, rotA )
!!$    if( r_modes(2) .gt. 0 ) deallocate( frB, rotB )
!!$    if( v_modes(1) .gt. 0 ) deallocate( fvA, vibA )
!!$    if( v_modes(2) .gt. 0 ) deallocate( fvB, vibB )

!!$    !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
!!$    if( vib_flag )then
!!$       ! Species Properties
!!$       neq_dens(1) = psi_neq(n,same)%vib_neq_dens
!!$       neq_dens(2) = psi_neq(m,same)%vib_neq_dens
!!$
!!$       ! Calculate total number of collisions
!!$       call compute_num_coll_pairs( num_colls, dt, dens(1)/dens(1), neq_dens(2)/dens(2), &
!!$            temp, m_red, coln_rms, vel_grid, n, m )
!!$       if( num_colls .eq. 0 ) num_colls = 1
!!$       tot_colls = tot_colls + num_colls
!!$
!!$       ! Calculate factor for depletion
!!$       factor = factor_coeff*dens(1)*neq_dens(2)/dble(num_colls)
!!$
!!$       ! Loop over all collisions
!!$       do coll = 1, num_colls
!!$          ! Pick velocities to collide
!!$          call pick_collision_partner( i(1), j(1), k(1), glA, psi_eq(ns,n,same)%cumul_df, vel_grid(n) )
!!$          call pick_collision_partner( i(2), j(2), k(2), glB, psi_neq(m,same)%cumul_vib_df, vel_grid(m) )
!!$
!!$          ! Calculate relative velocity between chosen collision partners
!!$          call compute_relative_velocity( i, j, k, g, vel_grid, n, m )
!!$
!!$          ! Check for self collision
!!$          if( g .lt. double_tol*vel_grid(n)%beta3_min .and. g .lt. double_tol*vel_grid(m)%beta3_min ) cycle
!!$
!!$          ! Calculate g*sigma term
!!$          g_sigma = sigma*g**vhs_exponent
!!$
!!$          ! Temporary arrays for phi
!!$          phiA = phi(n)%value( i(1), j(1), k(1) )
!!$          phiB = phi(m)%value( i(2), j(2), k(2) )
!!$
!!$          phiA_eq = psi_eq(ns,n,same)%cumul_df( glA ) - psi_eq(ns,n,same)%cumul_df( glA - 1 )
!!$          phiB_eq = psi_eq(ns,n,same)%cumul_df( glB ) - psi_eq(ns,n,same)%cumul_df( glB - 1 )
!!$
!!$          ! Find depletion energy levels
!!$          if( v_modes(1) .gt. 0 )then
!!$             vibA = psi_eq(ns,n,same)%vib*phiA_eq
!!$             call pick_level( vl(1), vibA, one, v_levels(1) )
!!$             depl_sign = one
!!$
!!$             if( r_modes(1) .gt. 0 )then
!!$                rotA = phi(n)%rot( :, i(1), j(1), k(1) )
!!$                call energy_fraction_array( frA, r_sumA, rotA, r_levels(1), dsgn(phiA) )
!!$             end if
!!$          end if
!!$
!!$          if( v_modes(2) .gt. 0 )then
!!$             vibB = phi(m)%vib( :, i(2), j(2), k(2) ) - psi_eq(ns,m,same)%vib*phiB_eq
!!$             call pick_level( vl(2), vibB, one, v_levels(2) )
!!$             depl_sign = &
!!$                  dsgn( phi(m)%vib( vl(2), i(2), j(2), k(2) ) - psi_eq(ns,m,same)%vib( vl(2) )*phiB_eq )
!!$
!!$             if( r_modes(2) .gt. 0 )then
!!$                rotB = phi(m)%rot( :, i(2), j(2), k(2) )
!!$                call energy_fraction_array( frB, r_sumB, rotB, r_levels(2), dsgn(phiB) )
!!$             end if
!!$          end if
!!$
!!$          ! Calculate depletion
!!$          depletion = depl_sign*depl_frac(3)*factor*g_sigma
!!$
!!$          ! deplete
!!$          call vr_deplete( vib_trans, phi(n), depletion, i(1), j(1), k(1), vl(1), &
!!$               frA, fvA, r_modes(1), v_modes(1) )
!!$          call vr_deplete( vib_trans, phi(m), depletion, i(2), j(2), k(2), vl(2), &
!!$               frB, fvB, r_modes(2), v_modes(2) )
!!$
!!$          ! replenish
!!$          call replenish_collision( phi, depletion, frA, frB, fvA, fvB, vl, i, j, k, &
!!$               g, n, m, vel_grid, m_red, vib_trans, molecule )
!!$
!!$       end do
!!$    end if
!!$
!!$    !**ROTATIONAL - TRANSLATIONAL**************************************************************************
!!$    if( rot_flag )then
!!$       ! Species Properties
!!$       neq_dens(1) = psi_neq(n,same)%rot_neq_dens
!!$       neq_dens(2) = psi_neq(m,same)%rot_neq_dens
!!$
!!$       ! Calculate total number of collisions
!!$       call compute_num_coll_pairs( num_colls, dt, dens(1)/dens(1), neq_dens(2)/dens(2), &
!!$            temp, m_red, coln_rms, vel_grid, n, m )
!!$       if( num_colls .eq. 0 ) num_colls = 1
!!$       tot_colls = tot_colls + num_colls
!!$
!!$       ! Calculate factor for depletion
!!$       factor = factor_coeff*dens(1)*neq_dens(2)/dble(num_colls)
!!$
!!$       do coll = 1, num_colls
!!$
!!$          ! Pick velocities to collide
!!$          call pick_collision_partner( i(1), j(1), k(1), glA, psi_eq(ns,n,same)%cumul_df, vel_grid(n) )
!!$          call pick_collision_partner( i(2), j(2), k(2), glB, psi_neq(m,same)%cumul_rot_df, vel_grid(m) )
!!$
!!$          ! Calculate relative velocity between chosen collision partners
!!$          call compute_relative_velocity( i, j, k, g, vel_grid, n, m )
!!$
!!$          ! Check for self collision
!!$          if( g .lt. double_tol*vel_grid(n)%beta3_min .and. g .lt. double_tol*vel_grid(m)%beta3_min ) cycle
!!$
!!$          ! Calculate g*sigma term
!!$          g_sigma = sigma*g**vhs_exponent
!!$
!!$          ! Temporary arrays for phi
!!$          phiA = phi(n)%value( i(1), j(1), k(1) )
!!$          phiB = phi(m)%value( i(2), j(2), k(2) )
!!$
!!$          phiA_eq = psi_eq(ns,n,same)%cumul_df( glA ) - psi_eq(ns,n,same)%cumul_df( glA - 1 )
!!$          phiB_eq = psi_eq(ns,n,same)%cumul_df( glB ) - psi_eq(ns,n,same)%cumul_df( glB - 1 )
!!$
!!$          ! Find depletion energy levels
!!$          if( r_modes(1) .gt. 0 )then
!!$             rotA = psi_eq(ns,n,same)%rot*phiA_eq
!!$             call pick_level( rl(1), rotA, one, r_levels(1) )
!!$             depl_sign = one
!!$
!!$             if( v_modes(1) .gt. 0 )then
!!$                vibA = phi(n)%vib( :, i(1), j(1), k(1) )
!!$                call energy_fraction_array( fvA, v_sumA, vibA, v_levels(1), dsgn(phiA) )
!!$             end if
!!$          end if
!!$
!!$          if( r_modes(2) .gt. 0 )then
!!$             rotB = phi(m)%rot( :, i(2), j(2), k(2) ) - psi_eq(ns,m,same)%rot*phiB_eq
!!$             call pick_level( rl(2), rotB, one, r_levels(2) )
!!$             depl_sign = &
!!$                  dsgn( phi(m)%rot( rl(2), i(2), j(2), k(2) ) - psi_eq(ns,m,same)%rot( rl(2) )*phiB_eq )
!!$
!!$             if( v_modes(2) .gt. 0 )then
!!$                vibB = phi(m)%vib( :, i(2), j(2), k(2) )
!!$                call energy_fraction_array( fvB, v_sumB, vibB, v_levels(2), dsgn(phiB) )
!!$             end if
!!$          end if
!!$
!!$          ! Calculate depletion
!!$          depletion = depl_sign*depl_frac(2)*factor*g_sigma
!!$
!!$          ! deplete
!!$          call vr_deplete( rot_trans, phi(n), depletion, i(1), j(1), k(1), rl(1), &
!!$               frA, fvA, r_modes(1), v_modes(1) )
!!$          call vr_deplete( rot_trans, phi(m), depletion, i(2), j(2), k(2), rl(2), &
!!$               frB, fvB, r_modes(2), v_modes(2) )
!!$
!!$          ! replenish
!!$          call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
!!$               g, n, m, vel_grid, m_red, rot_trans, molecule )
!!$
!!$       end do
!!$    end if
!!$
!!$    !**ELASTIC*********************************************************************************************
!!$
!!$    ! Species Properties
!!$    neq_dens(1) = psi_neq(n,same)%neq_dens
!!$    neq_dens(2) = psi_neq(m,same)%neq_dens
!!$
!!$    ! Calculate total number of collisions
!!$    call compute_num_coll_pairs( num_colls, dt, dens(1)/dens(1), neq_dens(2)/dens(2), &
!!$         temp, m_red, coln_rms, vel_grid, n, m )
!!$    if( num_colls .eq. 0 ) num_colls = 1
!!$    tot_colls = tot_colls + num_colls
!!$
!!$    ! Calculate factor for depletion
!!$    factor = factor_coeff*dens(1)*neq_dens(2)/dble(num_colls)
!!$
!!$    do coll = 1, num_colls
!!$       ! Pick velocities to collide
!!$       call pick_collision_partner( i(1), j(1), k(1), glA, psi_eq(ns,n,same)%cumul_df, vel_grid(n) )
!!$       call pick_collision_partner( i(2), j(2), k(2), glB, psi_neq(m,same)%cumul_df, vel_grid(m) )
!!$
!!$       depl_sign = psi_eq(ns,n,same)%sign( glA )
!!$       depl_sign = depl_sign*psi_neq(m,same)%sign( glB )
!!$
!!$       ! Calculate relative velocity between chosen collision partners
!!$       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )
!!$
!!$       ! Check for self collision
!!$       if( g .lt. double_tol*vel_grid(n)%beta3_min .and. g .lt. double_tol*vel_grid(m)%beta3_min ) cycle
!!$
!!$       ! Calculate g*sigma term
!!$       g_sigma = sigma*g**vhs_exponent
!!$
!!$       ! Temporary arrays for phi
!!$       phiA = phi(n)%value( i(1), j(1), k(1) )
!!$       phiB = phi(m)%value( i(2), j(2), k(2) )
!!$
!!$       ! Find depletion energy levels
!!$       if( v_modes(1) .gt. 0 )then
!!$          vibA = phi(n)%vib( :, i(1), j(1), k(1) )
!!$          call energy_fraction_array( fvA, v_sumA, vibA, v_levels(1), dsgn(phiA) )
!!$       end if
!!$
!!$       if( v_modes(2) .gt. 0 )then
!!$          vibB = phi(m)%vib( :, i(2), j(2), k(2) )
!!$          call energy_fraction_array( fvB, v_sumB, vibB, v_levels(2), dsgn(phiB) )
!!$       end if
!!$
!!$       if( r_modes(1) .gt. 0 )then
!!$          rotA = phi(n)%rot( :, i(1), j(1), k(1) )
!!$          call energy_fraction_array( frA, r_sumA, rotA, r_levels(1), dsgn(phiA) )
!!$       end if
!!$
!!$       if( r_modes(2) .gt. 0 )then
!!$          rotB = phi(m)%rot( :, i(2), j(2), k(2) )
!!$          call energy_fraction_array( frB, r_sumB, rotB, r_levels(2), dsgn(phiB) )
!!$       end if
!!$
!!$       ! Calculate depletion
!!$       depletion = depl_sign*depl_frac(1)*factor*g_sigma
!!$
!!$       ! deplete
!!$       call vr_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), 0, &
!!$            frA, fvA, r_modes(1), v_modes(1) )
!!$       call vr_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), 0, &
!!$            frB, fvB, r_modes(2), v_modes(2) )
!!$
!!$       ! replenish
!!$       call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
!!$            g, n, m, vel_grid, m_red, elastic, molecule )
!!$
!!$    end do
!!$
!!$
!!$    !======================================================================================================
!!$    ! Second collision integral: phi_noneq*phi_eq
!!$    !                            B-A, not A-B
!!$    !======================================================================================================
!!$    !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
!!$    if( vib_flag )then
!!$
!!$       ! Species Properties
!!$       neq_dens(1) = psi_neq(n,same)%vib_neq_dens
!!$       neq_dens(2) = psi_neq(m,same)%vib_neq_dens
!!$
!!$       ! Calculate total number of collisions
!!$       call compute_num_coll_pairs( num_colls, dt, neq_dens(1)/dens(1), dens(2)/dens(2), &
!!$            temp, m_red, coln_rms, vel_grid, n, m )
!!$       if( num_colls .eq. 0 ) num_colls = 1
!!$       tot_colls = tot_colls + num_colls
!!$
!!$       ! Calculate factor for depletion
!!$       factor = factor_coeff*neq_dens(1)*dens(2)/dble(num_colls)
!!$
!!$       ! Loop over all collisions
!!$       do coll = 1, num_colls
!!$          ! Pick velocities to collide
!!$          call pick_collision_partner( i(1), j(1), k(1), glA, psi_neq(n,same)%cumul_vib_df, vel_grid(n) )
!!$          call pick_collision_partner( i(2), j(2), k(2), glB, psi_eq(ns,m,same)%cumul_df, vel_grid(m) )
!!$
!!$          ! Calculate relative velocity between chosen collision partners
!!$          call compute_relative_velocity( i, j, k, g, vel_grid, n, m )
!!$
!!$          ! Check for self collision
!!$          if( g .lt. double_tol*vel_grid(n)%beta3_min .and. g .lt. double_tol*vel_grid(m)%beta3_min ) cycle
!!$
!!$          ! Calculate g*sigma term
!!$          g_sigma = sigma*g**vhs_exponent
!!$
!!$          ! Temporary arrays for phi
!!$          phiA = phi(n)%value( i(1), j(1), k(1) )
!!$          phiB = phi(m)%value( i(2), j(2), k(2) )
!!$
!!$          phiA_eq = psi_eq(ns,n,same)%cumul_df( glA ) - psi_eq(ns,n,same)%cumul_df( glA - 1 )
!!$          phiB_eq = psi_eq(ns,n,same)%cumul_df( glB ) - psi_eq(ns,n,same)%cumul_df( glB - 1 )
!!$
!!$          ! Find depletion energy levels
!!$          if( v_modes(1) .gt. 0 )then
!!$             vibA = phi(n)%vib( :, i(1), j(1), k(1) ) - psi_eq(ns,n,same)%vib*phiA_eq
!!$             call pick_level( vl(1), vibA, one, v_levels(1) )
!!$             depl_sign = &
!!$                  dsgn( phi(n)%vib( vl(1), i(1), j(1), k(1) ) - psi_eq(ns,n,same)%vib( vl(1) )*phiA_eq )
!!$
!!$             if( r_modes(1) .gt. 0 )then
!!$                rotA = phi(n)%rot( :, i(1), j(1), k(1) )
!!$                call energy_fraction_array( frA, r_sumA, rotA, r_levels(1), dsgn(phiA) )
!!$             end if
!!$          end if
!!$
!!$          if( v_modes(2) .gt. 0 )then
!!$             vibB = psi_eq(ns,m,same)%vib*phiB_eq
!!$             call pick_level( vl(2), vibB, one, v_levels(2) )
!!$             depl_sign = depl_sign*one
!!$
!!$             if( r_modes(2) .gt. 0 )then
!!$                rotB = phi(m)%rot( :, i(2), j(2), k(2) )
!!$                call energy_fraction_array( frB, r_sumB, rotB, r_levels(2), dsgn(phiB) )
!!$             end if
!!$          end if
!!$
!!$          ! Calculate depletion
!!$          depletion = depl_sign*depl_frac(3)*factor*g_sigma
!!$
!!$          ! deplete
!!$          call vr_deplete( vib_trans, phi(n), depletion, i(1), j(1), k(1), vl(1), &
!!$               frA, fvA, r_modes(1), v_modes(1) )
!!$          call vr_deplete( vib_trans, phi(m), depletion, i(2), j(2), k(2), vl(2), &
!!$               frB, fvB, r_modes(2), v_modes(2) )
!!$
!!$          ! replenish
!!$          call replenish_collision( phi, depletion, frA, frB, fvA, fvB, vl, i, j, k, &
!!$               g, n, m, vel_grid, m_red, vib_trans, molecule )
!!$
!!$       end do
!!$    end if
!!$
!!$    !**ROTATIONAL - TRANSLATIONAL**************************************************************************
!!$    if( rot_flag )then
!!$
!!$       ! Species Properties
!!$       neq_dens(1) = psi_neq(n,same)%rot_neq_dens
!!$       neq_dens(2) = psi_neq(m,same)%rot_neq_dens
!!$
!!$       ! Calculate total number of collisions
!!$       call compute_num_coll_pairs( num_colls, dt, neq_dens(1)/dens(1), dens(2)/dens(2), &
!!$            temp, m_red, coln_rms, vel_grid, n, m )
!!$       if( num_colls .eq. 0 ) num_colls = 1
!!$       tot_colls = tot_colls + num_colls
!!$
!!$       ! Calculate factor for depletion
!!$       factor = factor_coeff*neq_dens(1)*dens(2)/dble(num_colls)
!!$
!!$       do coll = 1, num_colls
!!$          ! Pick velocities to collide
!!$          call pick_collision_partner( i(1), j(1), k(1), glA, psi_neq(n,same)%cumul_rot_df, vel_grid(n) )
!!$          call pick_collision_partner( i(2), j(2), k(2), glB, psi_eq(ns,m,same)%cumul_df, vel_grid(m) )
!!$
!!$          ! Calculate relative velocity between chosen collision partners
!!$          call compute_relative_velocity( i, j, k, g, vel_grid, n, m )
!!$
!!$          ! Check for self collision
!!$          if( g .lt. double_tol*vel_grid(n)%beta3_min .and. g .lt. double_tol*vel_grid(m)%beta3_min ) cycle
!!$
!!$          ! Calculate g*sigma term
!!$          g_sigma = sigma*g**vhs_exponent
!!$
!!$          ! Temporary arrays for phi
!!$          phiA = phi(n)%value( i(1), j(1), k(1) )
!!$          phiB = phi(m)%value( i(2), j(2), k(2) )
!!$
!!$          phiA_eq = psi_eq(ns,n,same)%cumul_df( glA ) - psi_eq(ns,n,same)%cumul_df( glA - 1 )
!!$          phiB_eq = psi_eq(ns,n,same)%cumul_df( glB ) - psi_eq(ns,n,same)%cumul_df( glB - 1 )
!!$
!!$          ! Find depletion energy levels
!!$          if( r_modes(1) .gt. 0 )then
!!$             rotA = phi(n)%rot( :, i(1), j(1), k(1) ) - psi_eq(ns,n,same)%rot*phiA_eq
!!$             call pick_level( rl(1), rotA, one, r_levels(1) )
!!$             depl_sign = &
!!$                  dsgn( phi(n)%rot( rl(1), i(1), j(1), k(1) ) - psi_eq(ns,n,same)%rot( rl(1) )*phiA_eq )
!!$
!!$             if( v_modes(1) .gt. 0 )then
!!$                vibA = phi(n)%vib( :, i(1), j(1), k(1) )
!!$                call energy_fraction_array( fvA, v_sumA, vibA, v_levels(1), dsgn(phiA) )
!!$             end if
!!$          end if
!!$
!!$          if( r_modes(2) .gt. 0 )then
!!$             rotB = psi_eq(ns,m,same)%rot*phiB_eq
!!$             call pick_level( rl(2), rotB, one, r_levels(2) )
!!$             depl_sign = depl_sign*one
!!$
!!$             if( v_modes(2) .gt. 0 )then
!!$                vibB = phi(m)%vib( :, i(2), j(2), k(2) )
!!$                call energy_fraction_array( fvB, v_sumB, vibB, v_levels(2), dsgn(phiB) )
!!$             end if
!!$          end if
!!$
!!$          ! Calculate depletion
!!$          depletion = depl_sign*depl_frac(2)*factor*g_sigma
!!$
!!$          ! deplete
!!$          call vr_deplete( rot_trans, phi(n), depletion, i(1), j(1), k(1), rl(1), &
!!$               frA, fvA, r_modes(1), v_modes(1) )
!!$          call vr_deplete( rot_trans, phi(m), depletion, i(2), j(2), k(2), rl(2), &
!!$               frB, fvB, r_modes(2), v_modes(2) )
!!$
!!$          ! replenish
!!$          call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
!!$               g, n, m, vel_grid, m_red, rot_trans, molecule )
!!$
!!$       end do
!!$    end if
!!$
!!$    !**ELASTIC*********************************************************************************************
!!$
!!$    ! Species Properties
!!$    neq_dens(1) = psi_neq(n,same)%neq_dens
!!$    neq_dens(2) = psi_neq(m,same)%neq_dens
!!$
!!$    ! Calculate total number of collisions
!!$    call compute_num_coll_pairs( num_colls, dt, neq_dens(1)/dens(1), dens(2)/dens(2), &
!!$         temp, m_red, coln_rms, vel_grid, n, m )
!!$    if( num_colls .eq. 0 ) num_colls = 1
!!$    tot_colls = tot_colls + num_colls
!!$
!!$    ! Calculate factor for depletion
!!$    factor = factor_coeff*neq_dens(1)*dens(2)/dble(num_colls)
!!$
!!$    do coll = 1, num_colls
!!$       ! Pick velocities to collide
!!$       call pick_collision_partner( i(1), j(1), k(1), glA, psi_neq(n,same)%cumul_df, vel_grid(n) )
!!$       call pick_collision_partner( i(2), j(2), k(2), glB, psi_eq(ns,m,same)%cumul_df, vel_grid(m) )
!!$
!!$       depl_sign = psi_neq(n,same)%sign( glA )
!!$       depl_sign = depl_sign*psi_eq(ns,m,same)%sign( glB )
!!$
!!$       ! Calculate relative velocity between chosen collision partners
!!$       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )
!!$
!!$       ! Check for self collision
!!$       if( g .lt. double_tol*vel_grid(n)%beta3_min .and. g .lt. double_tol*vel_grid(m)%beta3_min ) cycle
!!$
!!$       ! Calculate g*sigma term
!!$       g_sigma = sigma*g**vhs_exponent
!!$
!!$       ! Temporary arrays for phi
!!$       phiA = phi(n)%value( i(1), j(1), k(1) )
!!$       phiB = phi(m)%value( i(2), j(2), k(2) )
!!$
!!$       ! Find depletion energy levels
!!$       if( v_modes(1) .gt. 0 )then
!!$          vibA = phi(n)%vib( :, i(1), j(1), k(1) )
!!$          call energy_fraction_array( fvA, v_sumA, vibA, v_levels(1), dsgn(phiA) )
!!$       end if
!!$
!!$       if( v_modes(2) .gt. 0 )then
!!$          vibB = phi(m)%vib( :, i(2), j(2), k(2) )
!!$          call energy_fraction_array( fvB, v_sumB, vibB, v_levels(2), dsgn(phiB) )
!!$       end if
!!$
!!$       if( r_modes(1) .gt. 0 )then
!!$          rotA = phi(n)%rot( :, i(1), j(1), k(1) )
!!$          call energy_fraction_array( frA, r_sumA, rotA, r_levels(1), dsgn(phiA) )
!!$       end if
!!$
!!$       if( r_modes(2) .gt. 0 )then
!!$          rotB = phi(m)%rot( :, i(2), j(2), k(2) )
!!$          call energy_fraction_array( frB, r_sumB, rotB, r_levels(2), dsgn(phiB) )
!!$       end if
!!$
!!$       ! Calculate depletion
!!$       depletion = depl_sign*depl_frac(1)*factor*g_sigma
!!$
!!$       ! deplete
!!$       call vr_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), 0, &
!!$            frA, fvA, r_modes(1), v_modes(1) )
!!$       call vr_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), 0, &
!!$            frB, fvB, r_modes(2), v_modes(2) )
!!$
!!$       ! replenish
!!$       call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
!!$            g, n, m, vel_grid, m_red, elastic, molecule )
!!$
!!$    end do
!!$
!!$    !======================================================================================================
!!$    ! Third collision integral: phi_noneq*phi_noneq
!!$    !                           A-B and B-A
!!$    !======================================================================================================
!!$    !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
!!$    if( vib_flag )then
!!$
!!$       ! Species Properties
!!$       neq_dens(1) = psi_neq(n,same)%vib_neq_dens
!!$       neq_dens(2) = psi_neq(m,same)%vib_neq_dens
!!$
!!$       ! Calculate total number of collisions
!!$       call compute_num_coll_pairs( num_colls, dt, neq_dens(1)/dens(1), neq_dens(2)/dens(2), &
!!$            temp, m_red, coln_rms, vel_grid, n, m )
!!$       if( num_colls .eq. 0 ) num_colls = 1
!!$       tot_colls = tot_colls + num_colls
!!$
!!$       ! Calculate factor for depletion
!!$       factor = factor_coeff*neq_dens(1)*neq_dens(2)/dble(num_colls)
!!$
!!$       ! Loop over all collisions
!!$       do coll = 1, num_colls
!!$          ! Pick velocities to collide
!!$          call pick_collision_partner( i(1), j(1), k(1), glA, psi_neq(n,same)%cumul_vib_df, vel_grid(n) )
!!$          call pick_collision_partner( i(2), j(2), k(2), glB, psi_neq(m,same)%cumul_vib_df, vel_grid(m) )
!!$
!!$          ! Calculate relative velocity between chosen collision partners
!!$          call compute_relative_velocity( i, j, k, g, vel_grid, n, m )
!!$
!!$          ! Check for self collision
!!$          if( g .lt. double_tol*vel_grid(n)%beta3_min .and. g .lt. double_tol*vel_grid(m)%beta3_min ) cycle
!!$
!!$          ! Calculate g*sigma term
!!$          g_sigma = sigma*g**vhs_exponent
!!$
!!$          ! Temporary arrays for phi
!!$          phiA = phi(n)%value( i(1), j(1), k(1) )
!!$          phiB = phi(m)%value( i(2), j(2), k(2) )
!!$
!!$          phiA_eq = psi_eq(ns,n,same)%cumul_df( glA ) - psi_eq(ns,n,same)%cumul_df( glA - 1 )
!!$          phiB_eq = psi_eq(ns,n,same)%cumul_df( glB ) - psi_eq(ns,n,same)%cumul_df( glB - 1 )
!!$
!!$          ! Find depletion energy levels
!!$          if( v_modes(1) .gt. 0 )then
!!$             vibA = phi(n)%vib( :, i(1), j(1), k(1) ) - psi_eq(ns,n,same)%vib*phiA_eq
!!$             call pick_level( vl(1), vibA, one, v_levels(1) )
!!$             depl_sign = &
!!$                  dsgn( phi(n)%vib( vl(1), i(1), j(1), k(1) ) - psi_eq(ns,n,same)%vib( vl(1) )*phiA_eq )
!!$
!!$             if( r_modes(1) .gt. 0 )then
!!$                rotA = phi(n)%rot( :, i(1), j(1), k(1) )
!!$                call energy_fraction_array( frA, r_sumA, rotA, r_levels(1), dsgn(phiA) )
!!$             end if
!!$          end if
!!$
!!$          if( v_modes(2) .gt. 0 )then
!!$             vibB = phi(m)%vib( :, i(2), j(2), k(2) ) - psi_eq(ns,m,same)%vib*phiB_eq
!!$             call pick_level( vl(2), vibB, one, v_levels(2) )
!!$             depl_sign = depl_sign*&
!!$                  dsgn( phi(m)%vib( vl(2), i(2), j(2), k(2) ) - psi_eq(ns,m,same)%vib( vl(2) )*phiB_eq )
!!$
!!$             if( r_modes(2) .gt. 0 )then
!!$                rotB = phi(m)%rot( :, i(2), j(2), k(2) )
!!$                call energy_fraction_array( frB, r_sumB, rotB, r_levels(2), dsgn(phiB) )
!!$             end if
!!$          end if
!!$
!!$          ! Calculate depletion
!!$          depletion = depl_sign*depl_frac(3)*factor*g_sigma
!!$
!!$          ! deplete
!!$          call vr_deplete( vib_trans, phi(n), depletion, i(1), j(1), k(1), vl(1), &
!!$               frA, fvA, r_modes(1), v_modes(1) )
!!$          call vr_deplete( vib_trans, phi(m), depletion, i(2), j(2), k(2), vl(2), &
!!$               frB, fvB, r_modes(2), v_modes(2) )
!!$
!!$          ! replenish
!!$          call replenish_collision( phi, depletion, frA, frB, fvA, fvB, vl, i, j, k, &
!!$               g, n, m, vel_grid, m_red, vib_trans, molecule )
!!$
!!$       end do
!!$    end if
!!$
!!$    !**ROTATIONAL - TRANSLATIONAL**************************************************************************
!!$    if( rot_flag )then
!!$
!!$       ! Species Properties
!!$       neq_dens(1) = psi_neq(n,same)%rot_neq_dens
!!$       neq_dens(2) = psi_neq(m,same)%rot_neq_dens
!!$
!!$       ! Calculate total number of collisions
!!$       call compute_num_coll_pairs( num_colls, dt, neq_dens(1)/dens(1), neq_dens(2)/dens(2), &
!!$            temp, m_red, coln_rms, vel_grid, n, m )
!!$       if( num_colls .eq. 0 ) num_colls = 1
!!$       tot_colls = tot_colls + num_colls
!!$
!!$       ! Calculate factor for depletion
!!$       factor = factor_coeff*neq_dens(1)*neq_dens(2)/dble(num_colls)
!!$
!!$       do coll = 1, num_colls
!!$          ! Pick velocities to collide
!!$          call pick_collision_partner( i(1), j(1), k(1), glA, psi_neq(n,same)%cumul_rot_df, vel_grid(n) )
!!$          call pick_collision_partner( i(2), j(2), k(2), glB, psi_neq(m,same)%cumul_rot_df, vel_grid(m) )
!!$
!!$          ! Calculate relative velocity between chosen collision partners
!!$          call compute_relative_velocity( i, j, k, g, vel_grid, n, m )
!!$
!!$          ! Check for self collision
!!$          if( g .lt. double_tol*vel_grid(n)%beta3_min .and. g .lt. double_tol*vel_grid(m)%beta3_min ) cycle
!!$
!!$          ! Calculate g*sigma term
!!$          g_sigma = sigma*g**vhs_exponent
!!$
!!$          ! Temporary arrays for phi
!!$          phiA = phi(n)%value( i(1), j(1), k(1) )
!!$          phiB = phi(m)%value( i(2), j(2), k(2) )
!!$
!!$          phiA_eq = psi_eq(ns,n,same)%cumul_df( glA ) - psi_eq(ns,n,same)%cumul_df( glA - 1 )
!!$          phiB_eq = psi_eq(ns,n,same)%cumul_df( glB ) - psi_eq(ns,n,same)%cumul_df( glB - 1 )
!!$
!!$          ! Find depletion energy levels
!!$          if( r_modes(1) .gt. 0 )then
!!$             rotA = phi(n)%rot( :, i(1), j(1), k(1) ) - psi_eq(ns,n,same)%rot*phiA_eq
!!$             call pick_level( rl(1), rotA, one, r_levels(1) )
!!$             depl_sign = &
!!$                  dsgn( phi(n)%rot( rl(1), i(1), j(1), k(1) ) - psi_eq(ns,n,same)%rot( rl(1) )*phiA_eq )
!!$
!!$             if( v_modes(1) .gt. 0 )then
!!$                vibA = phi(n)%vib( :, i(1), j(1), k(1) )
!!$                call energy_fraction_array( fvA, v_sumA, vibA, v_levels(1), dsgn(phiA) )
!!$             end if
!!$          end if
!!$
!!$          if( r_modes(2) .gt. 0 )then
!!$             rotB = phi(m)%rot( :, i(2), j(2), k(2) ) - psi_eq(ns,m,same)%rot*phiB_eq
!!$             call pick_level( rl(2), rotB, one, r_levels(2) )
!!$             depl_sign = depl_sign*&
!!$                  dsgn( phi(m)%rot( rl(2), i(2), j(2), k(2) ) - psi_eq(ns,m,same)%rot( rl(2) )*phiB_eq )
!!$
!!$             if( v_modes(2) .gt. 0 )then
!!$                vibB = phi(m)%vib( :, i(2), j(2), k(2) )
!!$                call energy_fraction_array( fvB, v_sumB, vibB, v_levels(2), dsgn(phiB) )
!!$             end if
!!$          end if
!!$
!!$          ! Calculate depletion
!!$          depletion = depl_sign*depl_frac(2)*factor*g_sigma
!!$
!!$          ! deplete
!!$          call vr_deplete( rot_trans, phi(n), depletion, i(1), j(1), k(1), rl(1), &
!!$               frA, fvA, r_modes(1), v_modes(1) )
!!$          call vr_deplete( rot_trans, phi(m), depletion, i(2), j(2), k(2), rl(2), &
!!$               frB, fvB, r_modes(2), v_modes(2) )
!!$
!!$          ! replenish
!!$          call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
!!$               g, n, m, vel_grid, m_red, rot_trans, molecule )
!!$
!!$       end do
!!$    end if
!!$
!!$    !**ELASTIC*********************************************************************************************
!!$
!!$    ! Species Properties
!!$    neq_dens(1) = psi_neq(n,same)%neq_dens
!!$    neq_dens(2) = psi_neq(m,same)%neq_dens
!!$
!!$    ! Calculate total number of collisions
!!$    call compute_num_coll_pairs( num_colls, dt, neq_dens(1)/dens(1), neq_dens(2)/dens(2), &
!!$         temp, m_red, coln_rms, vel_grid, n, m )
!!$    if( num_colls .eq. 0 ) num_colls = 1
!!$    tot_colls = tot_colls + num_colls
!!$
!!$    ! Calculate factor for depletion
!!$    factor = factor_coeff*neq_dens(1)*neq_dens(2)/dble(num_colls)
!!$
!!$    do coll = 1, num_colls
!!$       ! Pick velocities to collide
!!$       call pick_collision_partner( i(1), j(1), k(1), glA, psi_neq(n,same)%cumul_df, vel_grid(n) )
!!$       call pick_collision_partner( i(2), j(2), k(2), glB, psi_neq(m,same)%cumul_df, vel_grid(m) )
!!$
!!$       depl_sign = psi_neq(n,same)%sign( glA )
!!$       depl_sign = depl_sign*psi_neq(m,same)%sign( glB )
!!$
!!$       ! Calculate relative velocity between chosen collision partners
!!$       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )
!!$
!!$       ! Check for self collision
!!$       if( g .lt. double_tol*vel_grid(n)%beta3_min .and. g .lt. double_tol*vel_grid(m)%beta3_min ) cycle
!!$
!!$       ! Calculate g*sigma term
!!$       g_sigma = sigma*g**vhs_exponent
!!$
!!$       ! Temporary arrays for phi
!!$       phiA = phi(n)%value( i(1), j(1), k(1) )
!!$       phiB = phi(m)%value( i(2), j(2), k(2) )
!!$
!!$       ! Find depletion energy levels
!!$       if( v_modes(1) .gt. 0 )then
!!$          vibA = phi(n)%vib( :, i(1), j(1), k(1) )
!!$          call energy_fraction_array( fvA, v_sumA, vibA, v_levels(1), dsgn(phiA) )
!!$       end if
!!$
!!$       if( v_modes(2) .gt. 0 )then
!!$          vibB = phi(m)%vib( :, i(2), j(2), k(2) )
!!$          call energy_fraction_array( fvB, v_sumB, vibB, v_levels(2), dsgn(phiB) )
!!$       end if
!!$
!!$       if( r_modes(1) .gt. 0 )then
!!$          rotA = phi(n)%rot( :, i(1), j(1), k(1) )
!!$          call energy_fraction_array( frA, r_sumA, rotA, r_levels(1), dsgn(phiA) )
!!$       end if
!!$       
!!$       if( r_modes(2) .gt. 0 )then
!!$          rotB = phi(m)%rot( :, i(2), j(2), k(2) )
!!$          call energy_fraction_array( frB, r_sumB, rotB, r_levels(2), dsgn(phiB) )
!!$       end if
!!$
!!$       ! Calculate depletion
!!$       depletion = depl_sign*depl_frac(1)*factor*g_sigma
!!$
!!$       ! deplete
!!$       call vr_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), 0, &
!!$            frA, fvA, r_modes(1), v_modes(1) )
!!$       call vr_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), 0, &
!!$            frB, fvB, r_modes(2), v_modes(2) )
!!$
!!$       ! replenish
!!$       call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
!!$            g, n, m, vel_grid, m_red, elastic, molecule )
!!$
!!$    end do
!!$
!!$    if( r_modes(1) .gt. 0 ) deallocate( frA, rotA )
!!$    if( r_modes(2) .gt. 0 ) deallocate( frB, rotB )
!!$    if( v_modes(1) .gt. 0 ) deallocate( fvA, vibA )
!!$    if( v_modes(2) .gt. 0 ) deallocate( fvB, vibB )

    return
  end subroutine variance_reduction_kernel
  
  subroutine vr_deplete( coll_type, phi, depletion, i, j, k, level, fr, fv, r_modes, v_modes )

    use DistFunc

    implicit none

    double precision, dimension(:), intent(in) :: fr, fv
    double precision, intent(in) :: depletion

    integer, intent(in) :: coll_type, level, i, j, k, r_modes, v_modes

    type(DistFuncType) :: phi

    select case( coll_type )
    case( vib_trans )
       if( v_modes .gt. 0 )then
          phi%vib( level, i, j, k ) = phi%vib( level, i, j, k ) - depletion
       end if
       if( r_modes .gt. 0 )then
          phi%rot( :, i, j, k ) = phi%rot( :, i, j, k ) - depletion*fr
       end if
       phi%value( i, j, k ) = phi%value( i, j, k ) - depletion

    case( rot_trans )
       if( r_modes .gt. 0 )then
          phi%rot( level, i, j, k ) = phi%rot( level, i, j, k ) - depletion
       end if
       if( v_modes .gt. 0 )then
          phi%vib( :, i, j, k ) = phi%vib( :, i, j, k ) - depletion*fv
       end if
       phi%value( i, j, k ) = phi%value( i, j, k ) - depletion

    case( elastic )
       if( v_modes .gt. 0 )then
          phi%vib( :, i, j, k ) = phi%vib( :, i, j, k ) - depletion*fv
       end if
       if( r_modes .gt. 0 )then
          phi%rot( :, i, j, k ) = phi%rot( :, i, j, k ) - depletion*fr
       end if
       phi%value( i, j, k ) = phi%value( i, j, k ) - depletion

    case default

    end select

    return
  end subroutine vr_deplete

end module VarianceReduction
