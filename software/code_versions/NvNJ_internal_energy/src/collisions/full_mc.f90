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

  type(CumulativeDF), allocatable, dimension(:) :: psi

  public :: create_full_monte_carlo
  public :: destroy_full_monte_carlo
  public :: compute_psi_fmc
  public :: full_monte_carlo_kernel

contains

  subroutine create_full_monte_carlo( vel_grid, molecule )

    use VelocityGrid
    use PhysicalGrid
    
    implicit none

    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    type(MoleculeType), dimension(:), intent(in) :: molecule

    integer :: num_points, grid_ref, n
    integer :: r_modes, v_modes
    integer :: status

    allocate( psi( 1:num_species ), STAT=status )
    call allocate_error_check( status, "psi(FMC)" )

    do n = 1, num_species
       r_modes = molecule(n)%rot_modes
       v_modes = molecule(n)%vib_modes

       call get_spatial_reference( 1, grid_ref )
       num_points = vel_grid(n,grid_ref)%num_points

       allocate( psi(n)%cumul_df( 0:num_points ), STAT=status )
       call allocate_error_check( status, "psi%cumul_df(FMC)" )
       
       allocate( psi(n)%sign( 0:num_points ), STAT=status )
       call allocate_error_check( status, "psi%sign(FMC)" )
          
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

       deallocate( psi(n)%cumul_df, STAT=status )
       call deallocate_error_check( status, "psi%cumul_df(FMC)" )
       
       deallocate( psi(n)%sign, STAT=status )
       call deallocate_error_check( status, "psi%sign(FMC)" )

    end do

    deallocate( psi, STAT=status )
    call deallocate_error_check( status, "psi(FMC)" )

    return
  end subroutine destroy_full_monte_carlo
  
  subroutine compute_psi_fmc( phi, molecule, vel_grid, n )

    use DistFunc
    use VelocityGrid
    use MathUtilities

    implicit none

    type(DistFuncType), intent(in) :: phi
    type(MoleculeType), intent(in) :: molecule
    type(VelocityGridType), intent(in) :: vel_grid
    integer, intent(in) :: n

    double precision :: neg_dens, magnitude

    integer :: num_points_x, num_points_y, num_points
    integer :: l, i, j, k, m
    integer :: i_min, j_min, k_min
    integer :: r_modes, v_modes, r_levels, v_levels

    ! Grid properties
    i_min = vel_grid%i_min
    j_min = vel_grid%j_min
    k_min = vel_grid%k_min

    num_points_x = vel_grid%num_points_x
    num_points_y = vel_grid%num_points_y
    num_points   = vel_grid%num_points

    ! Set initial array values
    psi(n)%cumul_df(0) = zero
    psi(n)%sign(0)     = one

    neg_dens = zero

    ! Loop over all the velocity points and calculate the cumulative distribution
    do l = 1, num_points

       ! Find x,y,z coordinates of the velocity location
       call global2local_map( l-1, i_min, j_min, k_min, num_points_x, num_points_y, i, j, k )

       psi(n)%cumul_df(l) = psi(n)%cumul_df(l-1) + abs( phi%value(i,j,k) )
       psi(n)%sign(l) = dsgn( phi%value(i,j,k) )

       ! Test for negative density
       if( psi(n)%sign(l) .lt. zero ) neg_dens = neg_dens + phi%value(i,j,k)

    end do

    psi(n)%neg_dens = neg_dens

    return
  end subroutine compute_psi_fmc

  subroutine full_monte_carlo_kernel( phi, tot_colls, n, m, coln_rms, molecule, vel_grid, properties )

    use VelocityGrid
    use PhysicalProperties
    use ReplenishingCollisions
    use TimeStepping
    use MathUtilities
    use Scaling

    implicit none

    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    type(PropertiesType), intent(in) :: properties
    double precision, dimension(:,:), intent(in) :: coln_rms
    integer, intent(in) :: n, m

    type(DistFuncType), dimension(:) :: phi
    integer :: tot_colls

    double precision, allocatable, dimension(:) :: rotA, rotB, vibA, vibB
    double precision, allocatable, dimension(:) :: frA, frB, fvA, fvB
    double precision, allocatable, dimension(:,:) :: rotvibA, rotvibB, frvA, frvB
    double precision, dimension(5) :: depl_frac, depl_sign
    double precision, dimension(2) :: dens, temp, neg_dens, abs_dens
    double precision, dimension(2) :: mass, diam, omega
    double precision :: omega_AB, sigma, vhs_exponent
    double precision :: dt, m_red
    double precision :: g, g_sigma
    double precision :: depletion, factor
    double precision :: r_sumA, r_sumB, v_sumA, v_sumB
    double precision :: phiA, phiB, rot_sgn, vib_sgn
    double precision :: extra_dens

    integer, dimension(2) :: r_modes, v_modes, r_levels, v_levels
    integer, dimension(2) :: i, j, k
    integer, dimension(4) :: rl, vl
    integer :: num_colls, coll
    integer :: glA, glB

    integer :: status

    logical :: rot_flag, vib_flag

    ! Initialize flags for energy activation
    rot_flag = .false.
    vib_flag = .false.

    ! Species specific properties
    dens(1)     = properties%dens(n)
    temp(1)     = properties%temp(n)
    neg_dens(1) = psi(n)%neg_dens

    dens(2)     = properties%dens(m)
    temp(2)     = properties%temp(m)
    neg_dens(2) = psi(m)%neg_dens

    ! Species specific constants
    mass(1)  = molecule(n)%mass
    diam(1)  = molecule(n)%diam
    omega(1) = molecule(n)%omega
    
    mass(2)  = molecule(m)%mass
    diam(2)  = molecule(m)%diam
    omega(2) = molecule(m)%omega

    abs_dens(1) = dens(1) - two*neg_dens(1)
    abs_dens(2) = dens(2) - two*neg_dens(2)

    ! Species combination constants
    m_red = mass(1)*mass(2)/( mass(1) + mass(2) )

    omega_AB = one_half*( omega(1) + omega(2) )
    call calculate_cross_section( sigma, m_red, omega_AB, diam )

    vhs_exponent = two - two*omega_AB
    
    ! Time step
    call get_deltat( dt )

    ! Internal structure
    r_modes(1) = molecule(n)%rot_modes
    v_modes(1) = molecule(n)%vib_modes

    r_modes(2) = molecule(m)%rot_modes
    v_modes(2) = molecule(m)%vib_modes

    ! Set levels and allocate energy arrays
    if( r_modes(1) .gt. 0 .and. v_modes(1) .gt. 0 )then
       rot_flag = .true.
       vib_flag = .true.
       r_levels(1) = phi(n)%num_rot_levels
       v_levels(1) = phi(n)%num_vib_levels

       allocate( frA( 1:r_levels(1) ), rotA( 1:r_levels(1) ), STAT=status )
       call allocate_error_check( status, "frA, rotA (fmc)" )

       allocate( fvA( 1:v_levels(1) ), vibA( 1:v_levels(1) ), STAT=status )
       call allocate_error_check( status, "fvA, vibA (fmc)" )

       allocate( frvA( 1:r_levels(1), 1:v_levels(1) ), rotvibA( 1:r_levels(1), 1:v_levels(1) ), STAT=status )
       call allocate_error_check( status, "frvA, rotvibA (fmc)" )
    end if

    if( r_modes(2) .gt. 0 .and. v_modes(2) .gt. 0 )then
       rot_flag = .true.
       vib_flag = .true.
       r_levels(2) = phi(m)%num_rot_levels
       v_levels(2) = phi(m)%num_vib_levels

       allocate( frB( 1:r_levels(2) ), rotB( 1:r_levels(2) ), STAT=status )
       call allocate_error_check( status, "frB, rotB (fmc)" )

       allocate( fvB( 1:v_levels(2) ), vibB( 1:v_levels(2) ), STAT=status )
       call allocate_error_check( status, "fvB, vibB (fmc)" )

       allocate( frvB( 1:r_levels(2), 1:v_levels(2) ), rotvibB( 1:r_levels(2), 1:v_levels(2) ), STAT=status )
       call allocate_error_check( status, "frvB, rotvibB (fmc)" )
    end if

    ! Separate depletion into elastic and inelastic parts (elastic, r-t A, r-t B, v-t A, v-t B)
    call split_depletion( depl_frac, m_red, properties, molecule, n, m )

    ! Calculate the number of collisions
    call compute_num_coll_pairs( num_colls, dt, abs_dens(1), abs_dens(2), &
         temp, m_red, coln_rms, vel_grid, n, m )
    tot_colls = num_colls

    ! Calculate factor for depletion
    factor = dt*abs_dens(1)*abs_dens(2)/dble( num_colls )
    if( n .eq. m ) factor = one_half*factor

    extra_dens = zero

    do coll = 1, num_colls

       ! ISSUE: if each collision type is separate, I need to keep a consitant sum of extra mass without
       ! ISSUE: messing up rates.
       if( extra_dens .gt. factor )then
          extra_dens = extra_dens - factor
          num_colls = num_colls + 1
       end if

       ! Pick velocities to collide
       call pick_collision_partner( i(1), j(1), k(1), glA, psi(n)%cumul_df, vel_grid(n) )
       call pick_collision_partner( i(2), j(2), k(2), glB, psi(m)%cumul_df, vel_grid(m) )

       ! Relative velocity between chosen collision velocities
       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )

       ! Check for self collision
       if( g .lt. double_tol ) cycle

       ! calculate g*sigma term
       g_sigma = sigma*g**vhs_exponent

       ! Temporary variable for number density at colliding velocities
       phiA = phi(n)%value( i(1), j(1), k(1) )
       phiB = phi(m)%value( i(2), j(2), k(2) )

       ! Elastic sign
       depl_sign(1) = dsgn( phiA )*dsgn( phiB )

       if( r_modes(1) .gt. 0 .and. v_modes(1) .gt. 0 )then
          rotA = phi(n)%rot(:,i(1),j(1),k(1))
          vibA = phi(n)%vib(:,i(1),j(1),k(1))
          call pick_level( rl(1), rotA, one, r_levels(1) )
          call pick_level( vl(1), vibA, one, v_levels(1) )

          rotA = phi(n)%int_energy(:,vl(1),i(1),j(1),k(1))
          vibA = phi(n)%int_energy(rl(1),:,i(1),j(1),k(1))
          rotvibA = phi(n)%int_energy(:,:,i(1),j(1),k(1))

          rot_sgn = dsgn( sum(vibA) )
          vib_sgn = dsgn( sum(rotA) )

          depl_sign(2) = dsgn( phiB )*rot_sgn
          depl_sign(4) = dsgn( phiB )*vib_sgn

          call energy_fraction_array_1D( frA, r_sumA, rotA, r_levels(1), vib_sgn )
          call energy_fraction_array_1D( fvA, v_sumA, vibA, v_levels(1), rot_sgn )
          call energy_fraction_array_2D( frvA, phiA, rotvibA, r_levels(1), v_levels(1), dsgn(phiA) )
       end if

       if( r_modes(2) .gt. 0 .and. v_modes(2) .gt. 0 )then
          rotB = phi(m)%rot(:,i(2),j(2),k(2))
          vibB = phi(m)%vib(:,i(2),j(2),k(2))
          call pick_level( rl(2), rotB, one, r_levels(2) )
          call pick_level( vl(2), vibB, one, v_levels(2) )

          rotB = phi(m)%int_energy(:,vl(2),i(2),j(2),k(2))
          vibB = phi(m)%int_energy(rl(2),:,i(2),j(2),k(2))
          rotvibB = phi(m)%int_energy(:,:,i(2),j(2),k(2))

          rot_sgn = dsgn( sum(vibB) )
          vib_sgn = dsgn( sum(rotB) )

          depl_sign(3) = dsgn( phiA )*rot_sgn
          depl_sign(5) = dsgn( phiA )*vib_sgn

          call energy_fraction_array_1D( frB, r_sumB, rotB, r_levels(2), vib_sgn )
          call energy_fraction_array_1D( fvB, v_sumB, vibB, v_levels(2), rot_sgn )
          call energy_fraction_array_2D( frvB, phiB, rotvibB, r_levels(2), v_levels(2), dsgn(phiB) )
       end if

       !==ELASTIC=============================================================================================!
       depletion = depl_sign(1)*depl_frac(1)*factor*g_sigma

       call fmc_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), 0, &
            frA, fvA, frvA, r_modes(1), v_modes(1), r_levels(1), v_levels(1) )
       call fmc_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), 0, &
            frB, fvB, frvB, r_modes(2), v_modes(2), r_levels(2), v_levels(2) )

       call replenish_collision( phi, depletion, frA, frB, fvA, fvB, frvA, frvB, &
            rl, i, j, k, g, n, m, vel_grid, m_red, elastic, A, molecule )


       !==ROTATIONAL - TRANSLATIONAL==========================================================================!
       if( rot_flag )then
          if( r_modes(1) .gt. 0 )then 
             depletion = depl_sign(2)*depl_frac(2)*factor*g_sigma
             
             call fmc_deplete( rot_trans, phi(n), depletion, i(1), j(1), k(1), rl(1), &
                  frA, fvA, frvA, r_modes(1), v_modes(1), r_levels(1), v_levels(1) )
             call fmc_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), rl(2), &
                  frB, fvB, frvB, r_modes(2), v_modes(2), r_levels(2), v_levels(2) )

             call replenish_collision( phi, depletion, frA, frB, fvA, fvB, frvA, frvB, &
                  rl, i, j, k, g, n, m, vel_grid, m_red, rot_trans, A, molecule )

          end if


          if( r_modes(2) .gt. 0 )then 
             depletion = depl_sign(3)*depl_frac(3)*factor*g_sigma

             call fmc_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), rl(1), &
                  frA, fvA, frvA, r_modes(1), v_modes(1), r_levels(1), v_levels(1) )
             call fmc_deplete( rot_trans, phi(m), depletion, i(2), j(2), k(2), rl(2), &
                  frB, fvB, frvB, r_modes(2), v_modes(2), r_levels(2), v_levels(2) )

             call replenish_collision( phi, depletion, frA, frB, fvA, fvB, frvA, frvB, & 
                  rl, i, j, k, g, n, m, vel_grid, m_red, rot_trans, B, molecule )

          end if
       end if


       !==VIBRATIONAL - TRANSLATIONAL=========================================================================!
       if( vib_flag )then
          if( v_modes(1) .gt. 0 )then
             depletion = depl_sign(4)*depl_frac(4)*factor*g_sigma

             call fmc_deplete( vib_trans, phi(n), depletion, i(1), j(1), k(1), vl(1), &
                  frA, fvA, frvA, r_modes(1), v_modes(1), r_levels(1), v_levels(1) )
             call fmc_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), vl(2), &
                  frB, fvB, frvB, r_modes(2), v_modes(2), r_levels(2), v_levels(2) )

             call replenish_collision( phi, depletion, frA, frB, fvA, fvB, frvA, frvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, A, molecule )

          end if

          if( v_modes(2) .gt. 0 )then
             depletion = depl_sign(5)*depl_frac(5)*factor*g_sigma

             call fmc_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), vl(1), &
                  frA, fvA, frvA, r_modes(1), v_modes(1), r_levels(1), v_levels(1) )
             call fmc_deplete( vib_trans, phi(m), depletion, i(2), j(2), k(2), vl(2), &
                  frB, fvB, frvB, r_modes(2), v_modes(2), r_levels(2), v_levels(2) )

             call replenish_collision( phi, depletion, frA, frB, fvA, fvB, frvA, frvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, B, molecule )

          end if
       end if

       !==VIBRATIONAL - ROTATIONAL - TRANSLATIONAL============================================================!
       ! I do not currently calculate these energy exchanges since this process is somewhat unknown

    end do

    ! deallocate energy arrays
    if( r_modes(1) .gt. 0 ) deallocate( frA, rotA )
    if( r_modes(2) .gt. 0 ) deallocate( frB, rotB )
    if( v_modes(1) .gt. 0 ) deallocate( fvA, vibA )
    if( v_modes(2) .gt. 0 ) deallocate( fvB, vibB )

    return
  end subroutine full_monte_carlo_kernel

  subroutine fmc_deplete( coll_type, phi, depletion, i, j, k, level, fr, fv, frv, &
       r_modes, v_modes, r_levels, v_levels )

    implicit none

    double precision, dimension(:,:), intent(in) :: frv
    double precision, dimension(:), intent(in) :: fr, fv
    double precision, intent(in) :: depletion
    
    integer, intent(in) :: coll_type, level, i, j, k
    integer, intent(in) :: r_modes, v_modes, r_levels, v_levels

    type(DistFuncType) :: phi

    if( v_modes .gt. 0 .and. r_modes .gt. 0 )then
       select case( coll_type )
       case( vib_trans )
          phi%int_energy( 1:r_levels, level, i, j, k ) = &
               phi%int_energy( 1:r_levels, level, i, j, k ) - depletion*fr

       case( rot_trans )
          phi%int_energy( level, 1:v_levels, i, j, k ) = &
               phi%int_energy( level, 1:v_levels, i, j, k ) - depletion*fv

       case( elastic )
          phi%int_energy( 1:r_levels, 1:v_levels, i, j, k ) = &
               phi%int_energy( 1:r_levels, 1:v_levels, i, j, k ) - depletion*frv

       case default

       end select
    end if

    phi%value( i, j, k ) = phi%value( i, j, k ) - depletion

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
