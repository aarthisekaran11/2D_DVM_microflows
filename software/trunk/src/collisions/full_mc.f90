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
module FullMonteCarlo

  use ErrorCheck
  use CollisionUtilities
  use DistFunc
  use Constants
  use SpeciesAndReferenceData

  implicit none

  private

  integer, parameter :: elastic   = 0
  integer, parameter :: rot_trans = 1
  integer, parameter :: vib_trans = 2

  integer, parameter :: A = 1
  integer, parameter :: B = 2

  type(CumulativeDF), allocatable, dimension(:) :: cdf

  public :: create_full_monte_carlo
  public :: destroy_full_monte_carlo
  public :: compute_cdf_fmc
  public :: full_monte_carlo_kernel

contains

  subroutine create_full_monte_carlo( vel_grid )

    use VelocityGrid
    use PhysicalGrid

    implicit none

    type(VelocityGridType), dimension(:), intent(in) :: vel_grid

    integer :: num_points, n
    integer :: status

    allocate( cdf( 1:num_species ), STAT=status )
    call allocate_error_check( status, "cdf(FMC)" )

    do n = 1, num_species

       num_points = vel_grid(n)%num_points

       allocate( cdf(n)%cumul_df( 0:num_points ), STAT=status )
       call allocate_error_check( status, "cdf%cumul_df(FMC)" )

       allocate( cdf(n)%sign( 0:num_points ), STAT=status )
       call allocate_error_check( status, "cdf%sign(FMC)" )

    end do

    return
  end subroutine create_full_monte_carlo

  subroutine destroy_full_monte_carlo( molecule )

    implicit none

    type(MoleculeType), dimension(:), intent(in) :: molecule

    integer :: status, n
    integer :: r_modes, v_modes

    do n = 1, num_species
       r_modes = molecule(n)%rot_modes
       v_modes = molecule(n)%vib_modes

       deallocate( cdf(n)%cumul_df, STAT=status )
       call deallocate_error_check( status, "cdf%cumul_df(FMC)" )

       deallocate( cdf(n)%sign, STAT=status )
       call deallocate_error_check( status, "cdf%sign(FMC)" )

    end do

    deallocate( cdf, STAT=status )
    call deallocate_error_check( status, "cdf(FMC)" )

    return
  end subroutine destroy_full_monte_carlo

  subroutine compute_cdf_fmc( phi, vel_grid, n )

    use DistFunc
    use VelocityGrid
    use MathUtilities

    implicit none

    type(DistFuncType), intent(in) :: phi
    type(VelocityGridType), intent(in) :: vel_grid
    integer, intent(in) :: n

    double precision :: neg_dens

    double precision :: div_Nxy, div_Nx
    integer :: num_points_x, num_points_y, num_points
    integer :: l, i, j, k
    integer :: i_min, j_min, k_min

    integer :: ref
    double precision :: max

    ! Grid properties
    i_min = vel_grid%i_min
    j_min = vel_grid%j_min
    k_min = vel_grid%k_min

    num_points_x = vel_grid%num_points_x
    num_points_y = vel_grid%num_points_y
    num_points   = vel_grid%num_points
    div_Nxy = vel_grid%div_Nxy
    div_Nx  = vel_grid%div_Nx

    ! Set initial array values
    cdf(n)%cumul_df(0) = zero
    cdf(n)%sign(0)     = one

    neg_dens = zero

    ! Loop over all the velocity points and calculate the cumulative distribution
    do l = 1, num_points

       ! Find x,y,z coordinates of the velocity location
       call global2local_map( l-1, i_min, j_min, k_min, num_points_x, num_points_y, div_Nxy, div_Nx, i, j, k )

       ! Magnitude of the distribution function is used
       cdf(n)%cumul_df(l) = cdf(n)%cumul_df(l-1) + abs( phi%value(i,j,k) )
       cdf(n)%sign(l)     = dsgn( phi%value(i,j,k) )

       ! Test for negative density
       ! TODO: if I know dens and abs_dens, can I just calculate this outsid the loop?
!!$       if( cdf(n)%sign(l) .lt. zero ) neg_dens = neg_dens + phi%value(i,j,k)

    end do

    l = 1
    max = 0.001d0 * cdf(n)%cumul_df(num_points)
    do ref = 1, 1000
       cdf(n)%search_ref(ref,1) = l-1
       do while( dble(ref)*max .gt. cdf(n)%cumul_df(l) )
          l = l + 1
          if( l .gt. num_points )exit
       end do
       cdf(n)%search_ref(ref,2) = l
    end do
    cdf(n)%search_ref(1000,2) = num_points

    cdf(n)%neg_dens = neg_dens

    return
  end subroutine compute_cdf_fmc

  subroutine full_monte_carlo_kernel( phi, phi_old, tot_colls, n, m, coln_rms, molecule, vel_grid, properties )

    use VelocityGrid
    use PhysicalProperties
    use ReplenishingCollisions
    use TimeStepping
    use MathUtilities
    use Scaling
    use RandomNumberGeneration

    implicit none

    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    type(PropertiesType), intent(in) :: properties
    type(DistFuncType), dimension(:), intent(in) :: phi_old
    double precision, dimension(:,:), intent(in) :: coln_rms
    integer, intent(in) :: n, m

    type(DistFuncType), dimension(:) :: phi
    integer :: tot_colls

    integer :: num_points

    double precision, allocatable, dimension(:) :: rotA, rotB, vibA, vibB
    double precision, allocatable, dimension(:) :: frA, frB, fvA, fvB
    double precision, dimension(5) :: depl_frac, depl_sign
    double precision, dimension(2) :: dens, temp, neg_dens, abs_dens
    double precision, dimension(2) :: mass, diam, omega
    double precision, dimension(2) :: viscosity, temp_visc
    double precision, dimension(2) :: depletion, kin_depl
    double precision :: omega_AB, sigma, vhs_exponent
    double precision :: dt, m_red
    double precision :: Kn
    double precision :: g, g_sigma
    double precision :: factor!, depletion
    double precision :: r_sumA, r_sumB, v_sumA, v_sumB
    double precision :: deplA, deplB, depl_adj
    double precision :: phiA, phiB

    integer, dimension(2) :: r_modes, v_modes, r_levels, v_levels
    integer, dimension(2) :: i, j, k
    integer, dimension(4) :: rl, vl
    integer :: num_colls, coll
    integer :: glA, glB

    integer :: status

    logical :: rot_flag, vib_flag
    logical, dimension(2) :: single_valued_rot
    double precision, dimension(2) :: dE

    ! Initialize flags for energy activation
    rot_flag = .false.
    vib_flag = .false.

    single_valued_rot(1) = phi(n)%single_valued_rot
    single_valued_rot(2) = phi(m)%single_valued_rot

    ! Species specific properties
    num_points  = vel_grid(n)%num_points
    dens(1)     = properties%dens(n)
    temp(1)     = properties%tr_temp(n)
    neg_dens(1) = one_half * ( dens(1) - cdf(n)%cumul_df(num_points) )

    num_points  = vel_grid(m)%num_points
    dens(2)     = properties%dens(m)
    temp(2)     = properties%tr_temp(m)
    neg_dens(2) = one_half * ( dens(2) - cdf(m)%cumul_df(num_points) )

    ! Species specific constants
    mass(1)  = molecule(n)%mass
    diam(1)  = molecule(n)%diam
    omega(1) = molecule(n)%omega
    viscosity(1) = molecule(n)%viscosity
    temp_visc(1) = molecule(n)%T_viscosity

    mass(2)  = molecule(m)%mass
    diam(2)  = molecule(m)%diam
    omega(2) = molecule(m)%omega
    viscosity(2) = molecule(m)%viscosity
    temp_visc(2) = molecule(m)%T_viscosity

    dE(1) = properties%rot_energy(n) / properties%dens(n)
    dE(2) = properties%rot_energy(m) / properties%dens(m)

    ! Density used for depletion density and 
    ! number of collisions accounts for negative mass
    abs_dens(1) = dens(1) - two*neg_dens(1)
    abs_dens(2) = dens(2) - two*neg_dens(2)

    ! Species combination constants
    m_red = mass(1)*mass(2)/( mass(1) + mass(2) )

    call get_knudsen_number( Kn )

    omega_AB = one_half*( omega(1) + omega(2) )
    call calculate_cross_section( sigma, m_red, omega_AB, diam, viscosity, temp_visc )

    vhs_exponent = two - two*omega_AB

    ! Time step
    call get_deltat( dt )

    ! Internal structure
    r_modes(1) = molecule(n)%rot_modes
    v_modes(1) = molecule(n)%vib_modes

    r_modes(2) = molecule(m)%rot_modes
    v_modes(2) = molecule(m)%vib_modes

    ! Set levels and allocate energy arrays for diatomic/polyatomic molecules
    if( r_modes(1) .gt. 0 )then
       rot_flag = .true.
       r_levels(1) = phi(n)%num_rot_levels
       allocate( frA( 1:r_levels(1) ), rotA( 1:r_levels(1) ), STAT=status )
       call allocate_error_check( status, "frA, rotA (fmc)" )
    end if

    if( r_modes(2) .gt. 0 )then
       rot_flag = .true.
       r_levels(2) = phi(m)%num_rot_levels
       allocate( frB( 1:r_levels(2) ), rotB( 1:r_levels(2) ), STAT=status )
       call allocate_error_check( status, "frB, rotB (fmc)" )
    end if

    if( v_modes(1) .gt. 0 )then
       vib_flag = .true.
       v_levels(1) = phi(n)%num_vib_levels
       allocate( fvA( 1:v_levels(1) ), vibA( 1:v_levels(1) ), STAT=status )
       call allocate_error_check( status, "fvA, vibA (fmc)" )
    end if

    if( v_modes(2) .gt. 0 )then
       vib_flag = .true.
       v_levels(2) = phi(m)%num_vib_levels
       allocate( fvB( 1:v_levels(2) ), vibB( 1:v_levels(2) ), STAT=status )
       call allocate_error_check( status, "fvB, vibB (fmc)" )
    end if

    ! Separate depletion into elastic and inelastic parts (elastic, RT-A, RT-B, VT-A, VT-B)
    call split_depletion( depl_frac, m_red, properties, molecule, n, m )

    ! Calculate the number of collisions
    call compute_num_coll_pairs( num_colls, dt, Kn, abs_dens(1)/dens(1), abs_dens(2)/dens(2), &
         temp, mass, coln_rms, vel_grid, n, m )

    ! Calculate factor for depletion
    factor = ( dt / Kn ) * sigma * abs_dens(1) * abs_dens(2) / dble( num_colls )
    if( n .eq. m ) factor = one_half*factor

    do coll = 1, num_colls
       ! Pick velocities to collide
       call pick_collision_partner( i(1), j(1), k(1), glA, cdf(n)%cumul_df, cdf(n)%search_ref, vel_grid(n) )
       call pick_collision_partner( i(2), j(2), k(2), glB, cdf(m)%cumul_df, cdf(m)%search_ref, vel_grid(m) )

       ! Relative velocity between chosen collision velocities
       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )

       ! Check for self collision - if true, move to next collision
       if( g .lt. double_tol ) cycle

       ! calculate g*sigma term
       g_sigma = g**vhs_exponent
      
       ! Number density at colliding velocities before collision process begins
       phiA = phi_old(n)%value( i(1), j(1), k(1) )
       phiB = phi_old(m)%value( i(2), j(2), k(2) )

       ! Colliding partner - A
       if( r_modes(1) .gt. 0 .and. v_modes(1) .gt. 0 )then
          vibA = phi_old(n)%vib( :, i(1), j(1), k(1) )
          call absolute_normalized_array( fvA, vibA )

          rotA = phi_old(n)%rot( :, i(1), j(1), k(1) )
          if( single_valued_rot(1) .eqv. .false. )then
             call absolute_normalized_array( frA, rotA )

             r_sumA = sum( frA )
             v_sumA = sum( fvA )
             frA    = frA * abs( v_sumA )
             fvA    = fvA * abs( r_sumA )
             
             call pick_level( rl(1), frA, r_levels(1) )
             call pick_level( vl(3), fvA, v_levels(1) )

             call pick_level( vl(1), fvA, v_levels(1) )
             call pick_level( rl(3), frA, r_levels(1) )
             
             depl_sign(2) = dsgn( frA( rl(1) ) )
             depl_sign(4) = dsgn( fvA( vl(1) ) )
          else
             if( factor * g_sigma .gt. 1.d0 * abs( phiA ) )then
!!$                frA = rotA / ( factor * g_sigma )
                frA = dsgn( rotA(1) ) * dsgn( phiA ) * dE(1)
             else
                frA = rotA / phiA
             end if
!!$             frA = dE(1)

             r_sumA = dsgn( phiA )
             v_sumA = sum( fvA )

             rl(1) = 1
             call pick_level( vl(1), fvA, v_levels(1) )

             depl_sign(2) = dsgn( phiA )
             depl_sign(4) = dsgn( fvA( vl(1) ) )
          end if
         
       else if( r_modes(1) .gt. 0 )then
          rotA = phi_old(n)%rot( :, i(1), j(1), k(1) )
          if( single_valued_rot(1) .eqv. .false. )then
             call absolute_normalized_array( frA, rotA )
             
             r_sumA = (sum( frA ))
             v_sumA = one
             
             call pick_level( rl(1), frA, r_levels(1) )

             depl_sign(2) = dsgn( frA( rl(1) ) )
          else
             if( factor * g_sigma .gt. 1.d0 * abs( phiA ) )then
!!$                frA = rotA / ( factor * g_sigma )
                frA = dsgn( rotA(1) ) * dsgn( phiA ) * dE(1)
             else
                frA = rotA / phiA
             end if
!!$             frA = dE(1)

             r_sumA = dsgn( phiA )
             v_sumA = one

             rl(1) = 1

             depl_sign(2) = dsgn( phiA )
          end if


       else if( v_modes(1) .gt. 0 )then
          vibA = phi_old(n)%vib( :, i(1), j(1), k(1) )
          call absolute_normalized_array( fvA, vibA )

          v_sumA = ( sum( fvA ) )
          r_sumA = dsgn( v_sumA )

          call pick_level( vl(1), fvA, v_levels(1) )

          depl_sign(4) = dsgn( fvA( vl(1) ) )

       else
          r_sumA = dsgn( phiA )
          v_sumA = dsgn( phiA )

       end if

       ! Colliding partner - B
       if( r_modes(2) .gt. 0 .and. v_modes(2) .gt. 0 )then
          vibB = phi_old(m)%vib( :, i(2), j(2), k(2) )
          call absolute_normalized_array( fvB, vibB )

          rotB = phi_old(m)%rot( :, i(2), j(2), k(2) )
          if( single_valued_rot(2) .eqv. .false. )then
             call absolute_normalized_array( frB, rotB )
             
             r_sumB = sum( frB )
             v_sumB = sum( fvB )
             frB    = frB * abs( v_sumB )
             fvB    = fvB * abs( r_sumB )
             
             call pick_level( rl(2), frB, r_levels(2) )
             call pick_level( vl(4), fvB, v_levels(2) )

             call pick_level( vl(2), fvB, v_levels(2) )
             call pick_level( rl(4), frB, r_levels(2) )
             
             depl_sign(3) = dsgn( frB( rl(2) ) )
             depl_sign(5) = dsgn( fvB( vl(2) ) )
          else
             if( factor * g_sigma .gt. 1.d0 * abs( phiB ) )then
!!$                frB = rotB / ( factor * g_sigma )
                frB = dsgn( rotB(1) ) * dsgn( phiB ) * dE(2)
             else
                frB = rotB / phiB
             end if
!!$             frB = dE(2)

             r_sumB = dsgn( phiB )
             v_sumB = sum( fvB )

             rl(2) = 1
             call pick_level( vl(2), fvB, v_levels(2) )

             depl_sign(3) = dsgn( phiB )
             depl_sign(5) = dsgn( fvB( vl(2) ) )
          end if

       else if( r_modes(2) .gt. 0 )then
          rotB = phi_old(m)%rot( :, i(2), j(2), k(2) )
          if( single_valued_rot(2) .eqv. .false. )then
             call absolute_normalized_array( frB, rotB )
             
             r_sumB = ( sum( frB ) )
             v_sumB = one
             
             call pick_level( rl(2), frB, r_levels(2) )

             depl_sign(3) = dsgn( frB( rl(2) ) )
          else
             if( factor * g_sigma .gt. 1.d0 * abs( phiB ) )then
!!$                frB = rotB / ( factor * g_sigma )
                frB = dsgn( rotB(1) ) * dsgn( phiB ) * dE(2)
             else
                frB = rotB / phiB
             end if
!!$             frB = dE(2)

             r_sumB = one
             v_sumB = one

             rl(2) = 1

             depl_sign(3) = dsgn( phiB )
          end if

       else if( v_modes(2) .gt. 0 )then
          vibB = phi_old(m)%vib( :, i(2), j(2), k(2) )
          call absolute_normalized_array( fvB, vibB )

          v_sumB = ( sum( fvB ) )
          r_sumB = dsgn( v_sumB )

          call pick_level( vl(2), fvB, v_levels(2) )

          depl_sign(5) = dsgn( fvB( vl(2) ) )

       else
          r_sumB = dsgn( phiB )
          v_sumB = dsgn( phiB )

       end if

       !==ELASTIC=============================================================================================!
       deplA    = r_sumA * abs(v_sumA)
       deplB    = r_sumB * abs(v_sumB)
       depl_adj = deplA * deplB

       depletion(1) = deplB    * depl_frac(1) * factor * g_sigma
       depletion(2) = deplA    * depl_frac(1) * factor * g_sigma
       kin_depl(1)  = depl_adj * depl_frac(1) * factor * g_sigma
       kin_depl(2)  = depl_adj * depl_frac(1) * factor * g_sigma

       call fmc_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
            i(1), j(1), k(1), 0, 0, frA, fvA, r_modes(1), v_modes(1) )
       call fmc_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
            i(2), j(2), k(2), 0, 0, frB, fvB, r_modes(2), v_modes(2) )

       call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
            rl, vl, i, j, k, g, n, m, vel_grid, m_red, elastic, A, molecule )


       !==ROTATIONAL - TRANSLATIONAL==========================================================================!
       if( rot_flag )then
          if( r_modes(1) .gt. 0 )then
             deplA    = depl_sign(2)
             deplB    = r_sumB * abs(v_sumB)
             depl_adj = deplA * deplB
             
             depletion(1) = depl_adj * depl_frac(2) * factor * g_sigma
             depletion(2) = deplA    * depl_frac(2) * factor * g_sigma
             kin_depl(1)  = depl_adj * depl_frac(2) * factor * g_sigma
             kin_depl(2)  = depl_adj * depl_frac(2) * factor * g_sigma

             call fmc_deplete( rot_trans, phi(n), depletion(1), kin_depl(1), i(1), j(1), k(1), rl(1), vl(3), &
                  frA, fvA, r_modes(1), v_modes(1) )
             call fmc_deplete( elastic, phi(m), depletion(2), kin_depl(2), i(2), j(2), k(2), rl(2), vl(4), &
                  frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
                  rl, vl, i, j, k, g, n, m, vel_grid, m_red, rot_trans, A, molecule )

          end if


          if( r_modes(2) .gt. 0 )then
             deplA = r_sumA * abs(v_sumA)
             deplB = depl_sign(3)
             depl_adj = deplA * deplB

             depletion(1) = deplB    * depl_frac(3) * factor * g_sigma
             depletion(2) = deplA    * depl_frac(3) * factor * g_sigma * deplB
             kin_depl(1)  = depl_adj * depl_frac(3) * factor * g_sigma
             kin_depl(2)  = depl_adj * depl_frac(3) * factor * g_sigma

             call fmc_deplete( elastic, phi(n), depletion(1), kin_depl(1), i(1), j(1), k(1), rl(1), vl(3), &
                  frA, fvA, r_modes(1), v_modes(1) )
             call fmc_deplete( rot_trans, phi(m), depletion(2), kin_depl(2), i(2), j(2), k(2), rl(2), vl(4), &
                  frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
                  rl, vl, i, j, k, g, n, m, vel_grid, m_red, rot_trans, B, molecule )

          end if
       end if


       !==VIBRATIONAL - TRANSLATIONAL=========================================================================!
       if( vib_flag )then
          if( v_modes(1) .gt. 0 )then
             deplA = depl_sign(4)
             deplB = r_sumB * abs(v_sumB)
             depl_adj = deplA * deplB

             depletion(1) = deplB    * depl_frac(4) * factor * g_sigma * deplA
             depletion(2) = deplA    * depl_frac(4) * factor * g_sigma
             kin_depl(1)  = depl_adj * depl_frac(4) * factor * g_sigma
             kin_depl(2)  = depl_adj * depl_frac(4) * factor * g_sigma

             call fmc_deplete( vib_trans, phi(n), depletion(1), kin_depl(1), i(1), j(1), k(1), rl(3), vl(1), &
                  frA, fvA, r_modes(1), v_modes(1) )
             call fmc_deplete( elastic, phi(m), depletion(2), kin_depl(2), i(2), j(2), k(2), rl(4), vl(2), &
                  frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
                  rl, vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, A, molecule )

          end if

          if( v_modes(2) .gt. 0 )then
             deplA = r_sumA * abs(v_sumA)
             deplB = depl_sign(5)
             depl_adj = deplA * deplB

             depletion(1) = deplB    * depl_frac(5) * factor * g_sigma
             depletion(2) = deplA    * depl_frac(5) * factor * g_sigma * deplB
             kin_depl(1)  = depl_adj * depl_frac(5) * factor * g_sigma
             kin_depl(2)  = depl_adj * depl_frac(5) * factor * g_sigma

             call fmc_deplete( elastic, phi(n), depletion(1), kin_depl(1), i(1), j(1), k(1), rl(3), vl(1), &
                  frA, fvA, r_modes(1), v_modes(1) )
             call fmc_deplete( vib_trans, phi(m), depletion(2), kin_depl(2), i(2), j(2), k(2), rl(4), vl(2), &
                  frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
                  rl, vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, B, molecule )

          end if
       end if

       !==VIBRATIONAL - ROTATIONAL - TRANSLATIONAL============================================================!
       ! I do not currently calculate these energy exchanges since this process is somewhat unknown

    end do

    tot_colls = num_colls

    ! deallocate energy arrays
    if( r_modes(1) .gt. 0 ) deallocate( frA, rotA )
    if( r_modes(2) .gt. 0 ) deallocate( frB, rotB )
    if( v_modes(1) .gt. 0 ) deallocate( fvA, vibA )
    if( v_modes(2) .gt. 0 ) deallocate( fvB, vibB )

    return
  end subroutine full_monte_carlo_kernel

  subroutine fmc_deplete( coll_type, phi, depletion, kin_depl, i, j, k, &
       r_level, v_level, fr, fv, r_modes, v_modes )

    implicit none

    double precision, dimension(:), intent(in) :: fr, fv
    double precision, intent(in) :: depletion, kin_depl

    integer, intent(in) :: coll_type, r_level, v_level, i, j, k, r_modes, v_modes

    type(DistFuncType) :: phi

    select case( coll_type )
    case( vib_trans )
       if( v_modes .gt. 0 ) phi%vib( v_level, i, j, k ) = phi%vib( v_level, i, j, k ) - kin_depl
       if( r_modes .gt. 0 ) phi%rot( r_level, i, j, k ) = phi%rot( r_level, i, j, k ) - kin_depl
       phi%value( i, j, k ) = phi%value( i, j, k ) - kin_depl

    case( rot_trans )
       if( r_modes .gt. 0 )then
          if( phi%single_valued_rot .eqv. .false. )then
             phi%rot( r_level, i, j, k ) = phi%rot( r_level, i, j, k ) - kin_depl
          else
             phi%rot( :, i, j, k ) = phi%rot( :, i, j, k ) - depletion*fr
          end if
       end if
       if( v_modes .gt. 0 ) phi%vib( v_level, i, j, k ) = phi%vib( v_level, i, j, k ) - kin_depl
       phi%value( i, j, k ) = phi%value( i, j, k ) - kin_depl

    case( elastic )
       if( v_modes .gt. 0 ) phi%vib( :, i, j, k ) = phi%vib( :, i, j, k ) - depletion*fv
       if( r_modes .gt. 0 ) phi%rot( :, i, j, k ) = phi%rot( :, i, j, k ) - depletion*fr
       phi%value( i, j, k ) = phi%value( i, j, k ) - kin_depl

    case default

    end select

    return
  end subroutine fmc_deplete

  subroutine test_level_depletion( depletion, extra_dens, phi_val, frac )

    implicit none

    double precision, intent(in) :: phi_val, frac

    double precision :: depletion, extra_dens

    if( phi_val .gt. zero .and. phi_val - frac*depletion .lt. zero )then
       extra_dens = extra_dens + depletion - phi_val/frac
       depletion  = phi_val/frac
    end if

    return
  end subroutine test_level_depletion

end module FullMonteCarlo
