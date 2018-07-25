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

  type(CumulativeDF), allocatable, dimension(:,:,:) :: cdf_eq
  type(CumulativeDF), allocatable, dimension(:,:)   :: cdf_neq

  integer, parameter :: elastic   = 0
  integer, parameter :: rot_trans = 1
  integer, parameter :: vib_trans = 2

  integer, parameter :: A = 1
  integer, parameter :: B = 2

  integer, parameter :: AA = 1
  integer, parameter :: AB = 2

  double precision :: cdf_tolerance

  public :: set_cdf_tolerance
  public :: create_variance_reduction
  public :: destroy_variance_reduction
  public :: compute_cdf_vr
  public :: variance_reduction_kernel

contains

  subroutine set_cdf_tolerance( cdf_tolerance_in )
    
    implicit none

    double precision, intent(in) :: cdf_tolerance_in

    cdf_tolerance = cdf_tolerance_in

    return
  end subroutine set_cdf_tolerance

  subroutine create_variance_reduction( vel_grid, phi, molecule )

    use VelocityGrid
    use DistFunc 
    use PhysicalGrid

    implicit none

    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(MoleculeType), dimension(:), intent(in) :: molecule

    integer :: num_points, nx_space, ny_space, grid_ref
    integer :: r_modes, v_modes, r_levels, v_levels
    integer :: n, same, s
    integer :: status

    call get_nspace( nx_space, ny_space )

    allocate( cdf_eq( 1:nx_space, 1:num_species, 1:2 ), STAT=status )
    call allocate_error_check( status, "cdf_eq(VR)" )

    allocate( cdf_neq( 1:num_species, 1:2 ), STAT=status )
    call allocate_error_check( status, "cdf_neq(VR)" )

    do same = 1, 2
       do s = 1, num_species

          call get_spatial_reference( 1, grid_ref )
          num_points = vel_grid( s, grid_ref )%num_points

          allocate( cdf_neq(s,same)%cumul_df( 0:num_points ), STAT=status )
          call allocate_error_check( status, "cdf_neq%cumul_df(VR)" )

          allocate( cdf_neq(s,same)%sign( 0:num_points ), STAT=status )
          call allocate_error_check( status, "cdf_neq%sign(VR)" )
          
          r_modes = molecule(s)%rot_modes
          v_modes = molecule(s)%vib_modes

          allocate( cdf_neq(s,same)%kin_df( 0:num_points ), STAT=status )
          call allocate_error_check( status, "cdf_neq%kin_df(VR)" )

          allocate( cdf_neq(s,same)%rot_df( 0:num_points ), STAT=status )
          call allocate_error_check( status, "cdf_neq%rot_df(VR)" )

          allocate( cdf_neq(s,same)%vib_df( 0:num_points ), STAT=status )
          call allocate_error_check( status, "cdf_neq%vib_df(VR)" )

          allocate( cdf_neq(s,same)%coeff( 0:num_points ), STAT=status )
          call allocate_error_check( status, "cdf_neq%coeff" )
          allocate( cdf_neq(s,same)%coeff2( 0:num_points ), STAT=status )
          call allocate_error_check( status, "cdf_neq%coeff2" )

          do n = 1, nx_space

             call get_spatial_reference( n, grid_ref )
             num_points = vel_grid( s, grid_ref )%num_points

             allocate( cdf_eq(n,s,same)%cumul_df( 0:num_points ), STAT=status )
             call allocate_error_check( status, "cdf_eq%cumul_df(VR)" )

             allocate( cdf_eq(n,s,same)%sign( 0:num_points ), STAT=status )
             call allocate_error_check( status, "cdf_eq%sign(VR)" )     
             
             cdf_eq(n,s,same)%dens = zero
             cdf_eq(n,s,same)%u    = zero
             cdf_eq(n,s,same)%v    = zero
             cdf_eq(n,s,same)%w    = zero
             cdf_eq(n,s,same)%temp = zero

             if( r_modes .gt. 0 )then
                r_levels = phi(s,n)%num_rot_levels
                allocate( cdf_eq(n,s,same)%rot_eq( 1:r_levels ), STAT=status )
                call allocate_error_check( status, "cdf_eq%rot(VR)" )
             end if

             if( v_modes .gt. 0 )then
                v_levels = phi(s,n)%num_vib_levels
                allocate( cdf_eq(n,s,same)%vib_eq( 1:v_levels ), STAT=status )
                call allocate_error_check( status, "cdf_eq%vib(VR)" )
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

    integer :: nx_space, ny_space, n, s, same
    integer :: r_modes, v_modes
    integer :: status

    call get_nspace( nx_space, ny_space )

    do same = 1, 2
       do s = 1, num_species

          deallocate( cdf_neq(s,same)%cumul_df, STAT=status )
          call deallocate_error_check( status, "cdf_neq%cumul_df(VR)" )

          deallocate( cdf_neq(s,same)%sign, STAT=status )
          call deallocate_error_check( status, "cdf_neq%sign(VR)" )

          r_modes = molecule(s)%rot_modes
          v_modes = molecule(s)%vib_modes

          deallocate( cdf_neq(s,same)%kin_df, STAT=status )
          call deallocate_error_check( status, "cdf_neq%kin_df(VR)" )

          deallocate( cdf_neq(s,same)%rot_df, STAT=status )
          call deallocate_error_check( status, "cdf_neq%rot_df(VR)" )

          deallocate( cdf_neq(s,same)%vib_df, STAT=status )
          call deallocate_error_check( status, "cdf_neq%vib_df(VR)" )

          do n = 1, nx_space

             deallocate( cdf_eq(n,s,same)%cumul_df, STAT=status )
             call deallocate_error_check( status, "cdf_eq%cumul_df(VR)" )

             deallocate( cdf_eq(n,s,same)%sign, STAT=status )
             call deallocate_error_check( status, "cdf_eq%sign(VR)" )

             if( r_modes .gt. 0 )then
                deallocate( cdf_eq(n,s,same)%rot_eq, STAT=status )
                call deallocate_error_check( status, "cdf_eq%rot(VR)" )
             end if

             if( v_modes .gt. 0 )then
                deallocate( cdf_eq(n,s,same)%vib_eq, STAT=status )
                call deallocate_error_check( status, "cdf_eq%vib(VR)" )
             end if

          end do
       end do
    end do

    deallocate( cdf_eq, STAT=status )
    call deallocate_error_check( status, "cdf_eq(VR)" )

    deallocate( cdf_neq, STAT=status )
    call deallocate_error_check( status, "cdf_neq(VR)" )

    return
  end subroutine destroy_variance_reduction

  subroutine compute_cdf_vr( phi, molecule, vel_grid, properties, new_eq_flag, same, s, n )

    use DistFunc
    use VelocityGrid
    use PhysicalProperties

    implicit none

    type(DistFuncType), intent(in) :: phi
    type(MoleculeType), intent(in) :: molecule
    type(VelocityGridType), intent(in) :: vel_grid
    type(PropertiesType), intent(in) :: properties
    integer, intent(in) :: same, s, n
    logical :: new_eq_flag

    double precision, allocatable, dimension(:) :: temporary

    double precision :: x, y, z
    double precision :: beta_x, beta_y, beta_z, beta3
    double precision :: dens, u, v, w, temp, mass
    double precision :: maxwell_coeff, phi_eq, phi_neq, rot_neq, vib_neq, tot_neq
    double precision :: coeff
    double precision :: energy_level, fraction

    double precision :: dens_ratio

    integer :: num_points_x, num_points_y, num_points
    integer :: l, i, j, k, p
    integer :: i_min, j_min, k_min
    integer :: r_modes, v_modes, r_levels, v_levels
    integer :: mol_type

    integer :: ref
    double precision :: max

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

     maxwell_coeff = dens * sqrt( mass * mass * mass / ( pi * pi * pi * temp * temp * temp ) )

     new_eq_flag = .false.
     if( abs( dens - cdf_eq(n,s,same)%dens ) .gt. cdf_tolerance*dens .or. &
          abs( u - cdf_eq(n,s,same)%u ) .gt. cdf_tolerance*u .or. &
          abs( temp - cdf_eq(n,s,same)%temp ) .gt. cdf_tolerance*temp )then
        new_eq_flag = .true.
     end if

    ! Set initial array values
    cdf_eq(n,s,same)%cumul_df(0) = zero
    cdf_neq(s,same)%cumul_df(0)  = zero

    cdf_neq(s,same)%kin_df(0) = zero
    cdf_neq(s,same)%rot_df(0) = zero
    cdf_neq(s,same)%vib_df(0) = zero

    cdf_eq(n,s,same)%sign(0)  = one
    cdf_neq(s,same)%sign(0)   = one

    ! Compute quantities which are independent of velocity
    if( new_eq_flag )then
       ! Coefficient for calculating the Maxwellian
       

       do l = 1, num_points
          
          ! Find x,y,z coordinates of the velocity location
          call global2local_map( l-1, i_min, j_min, k_min, num_points_x, num_points_y, i, j, k )
          
          ! Velocity and grid spacing
          x = vel_grid%x(i)
          y = vel_grid%y(j)
          z = vel_grid%z(k)

          beta_x = vel_grid%beta_x(i)
          beta_y = vel_grid%beta_y(j)
          beta_z = vel_grid%beta_z(k)
          beta3  = beta_x * beta_y * beta_z

          call compute_maxwellian( x, y, z, maxwell_coeff, mass, u, v, w, temp, phi_eq )
          phi_eq = phi_eq * beta3

          cdf_eq(n,s,same)%cumul_df(l) = cdf_eq(n,s,same)%cumul_df(l-1) + phi_eq
          cdf_eq(n,s,same)%sign(l)     = one

       end do

       ! Adjust equilibrium density for discretiazation error
       dens_ratio = dens / cdf_eq(n,s,same)%cumul_df(num_points)

       cdf_eq(n,s,same)%cumul_df(:) = dens_ratio * cdf_eq(n,s,same)%cumul_df(:)

       cdf_eq(n,s,same)%dens = dens
       cdf_eq(n,s,same)%u    = u
       cdf_eq(n,s,same)%temp = temp

       ! Calculate the rotational equilibrium distribution function
       if( r_modes .gt. 0 )then
          r_levels  = phi%num_rot_levels
          allocate( temporary(1:r_levels) ) ! Second term in subroutine call must be an array
          call compute_rot_distribution( cdf_eq(n,s,same)%rot_eq(:), temporary, molecule, temp, r_levels, s )
          deallocate( temporary )
       end if
       
       ! Calculate the vibrational equilibrium distribution function
       if( v_modes .gt. 0 )then
          v_levels = phi%num_vib_levels
          allocate( temporary(1:v_levels) ) ! Second term in subroutine call must be an array
          call compute_vib_distribution( cdf_eq(n,s,same)%vib_eq(:), temporary, molecule, temp, v_levels, s )
          deallocate( temporary )
       end if

    end if

    ! Loop over all the velocity points and calculate the equilibrium and deviation distributions
    do l = 1, num_points

       ! Find x,y,z coordinates of the velocity location
       call global2local_map( l-1, i_min, j_min, k_min, num_points_x, num_points_y, i, j, k )

       phi_eq  = cdf_eq(n,s,same)%cumul_df(l) - cdf_eq(n,s,same)%cumul_df(l-1)
       phi_neq = phi%value(i,j,k) - phi_eq

       ! Construct cumulative distribution function and record associated sign values
       cdf_neq(s,same)%kin_df(l) = phi_neq ! NOT stored as a cumulative function
       cdf_neq(s,same)%sign(l)   = dsgn( phi_neq )

       ! The rotational deviation density at each velocity location is the sum of the absolute magnitude 
       ! of each rotational level
       if( r_modes .gt. 0 )then
          rot_neq = sum( abs( phi%rot(:,i,j,k) - cdf_eq(n,s,same)%rot_eq(:) * phi_eq ) )
       else
          rot_neq = zero
       end if

       cdf_neq(s,same)%rot_df(l) = rot_neq ! NOT stored as a cumulative function

       ! The vibrational deviation density at each velocity location is the sum of the absolute magnitude 
       ! of each vibrational level
       if( v_modes .gt. 0 )then
          vib_neq = sum( abs( phi%vib(:,i,j,k) - cdf_eq(n,s,same)%vib_eq(:) * phi_eq ) )
       else
          vib_neq = zero
       end if

       cdf_neq(s,same)%vib_df(l) = vib_neq ! NOT stored as a cumulative function
       
       ! Calculate the total non-equilibrium density and set coefficients for later?
       if( r_modes .gt. 0 .and. v_modes .gt. 0 )then
          tot_neq = ( phi_eq * rot_neq + phi_eq * vib_neq + rot_neq * vib_neq ) &
               / sqrt( ( phi_eq + rot_neq ) * ( phi_eq + vib_neq ) )
          
          if( tot_neq .gt. zero )then
             cdf_neq(s,same)%coeff(l)     = one / ( phi_eq * rot_neq + phi_eq * vib_neq + rot_neq * vib_neq )
             cdf_neq(s,same)%coeff2(l)    = one
          else
             cdf_neq(s,same)%coeff(l)     = zero
             cdf_neq(s,same)%coeff2(l)    = zero
          end if

          ! Phi equilibrium is reduced because some is used in the deviation distribution
          phi_eq = phi_eq * phi_eq / sqrt( ( phi_eq + rot_neq ) * ( phi_eq + vib_neq ) )

       else if( r_modes .gt. 0 )then
          tot_neq = rot_neq
          cdf_neq(s,same)%coeff(l)     = one
          cdf_neq(s,same)%coeff2(l)    = zero

       else if( v_modes .gt. 0 )then
          tot_neq = vib_neq
          cdf_neq(s,same)%coeff(l)     = one
          cdf_neq(s,same)%coeff2(l)    = zero

       else
          tot_neq = abs( phi_neq )

       end if

       ! Cumulative distribution function
       cdf_neq(s,same)%cumul_df(l)  = cdf_neq(s,same)%cumul_df(l-1)  + tot_neq
       cdf_eq(n,s,same)%cumul_df(l) = cdf_eq(n,s,same)%cumul_df(l-1) + phi_eq

    end do

    ! Deviation density - n_dev
    cdf_neq(s,same)%neq_dens = cdf_neq(s,same)%cumul_df(num_points)

    ! Equilibrium density - n_eq
    cdf_eq(n,s,same)%eq_dens = cdf_eq(n,s,same)%cumul_df(num_points)

    l = 1
    max = 0.001d0 * cdf_eq(n,s,same)%eq_dens
    do ref = 1, 1000
       cdf_eq(n,s,same)%search_ref(ref,1) = l-1
       do while( dble(ref)*max .gt. cdf_eq(n,s,same)%cumul_df(l) )
          l = l + 1
          if( l .gt. num_points )exit
       end do
       cdf_eq(n,s,same)%search_ref(ref,2) = l
    end do
    cdf_eq(n,s,same)%search_ref(1000,2) = num_points

    l = 1
    max = 0.001d0 * cdf_neq(s,same)%neq_dens
    do ref = 1, 1000
       cdf_neq(s,same)%search_ref(ref,1) = l-1
       do while( dble(ref)*max .gt. cdf_neq(s,same)%cumul_df(l) )
          l = l + 1
          if( l .gt. num_points )exit
       end do
       cdf_neq(s,same)%search_ref(ref,2) = l
    end do
    cdf_neq(s,same)%search_ref(1000,2) = num_points

    return
  end subroutine compute_cdf_vr

  subroutine variance_reduction_kernel( phi, phi_old, tot_colls, n, m, ns, coln_rms, &
       molecule, vel_grid, properties, ntime )

    use DistFunc
    use VelocityGrid
    use PhysicalProperties
    use ReplenishingCollisions
    use TimeStepping
    use Scaling
    use RandomNumberGeneration

    implicit none

    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    type(DistFuncType), dimension(:), intent(in) :: phi_old
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(PropertiesType), intent(in) :: properties
    double precision, dimension(:,:), intent(in) :: coln_rms
    integer, intent(in) :: n, m, ns, ntime

    type(DistFuncType), dimension(:) :: phi
    integer :: tot_colls

    logical :: rot_flag, vib_flag

    double precision, dimension(:), allocatable :: frA, frB, fvA, fvB
    double precision, dimension(:), allocatable :: rotA, rotB, vibA, vibB

    double precision, dimension(:), allocatable :: rotA_pos, rotA_neg
    double precision, dimension(:), allocatable :: rotB_pos, rotB_neg
    double precision, dimension(:), allocatable :: vibA_pos, vibA_neg
    double precision, dimension(:), allocatable :: vibB_pos, vibB_neg

    double precision :: coeff
    double precision :: sum_pos, sum_neg
    double precision :: rot_sum_pos, rot_sum_neg, vib_sum_pos, vib_sum_neg
    double precision :: rot_sum, vib_sum
    double precision :: sumA, sumB
    double precision :: sum_rtA, sum_rtB, sum_vtA, sum_vtB
    double precision :: deplA, deplB, depl_adj

    double precision, dimension(5) :: depl_frac, depl_sign
    double precision, dimension(2) :: mass, diam, omega
    double precision, dimension(2) :: dens, eq_dens, neq_dens, temp, kin_neq_dens, kin_eq_dens
    double precision, dimension(2) :: depletion, kin_depl
    double precision :: m_red, omega_AB, vhs_exponent, sigma
    double precision :: dt, factor_coeff
    double precision :: factor!, depletion
    double precision :: g, g_sigma

    double precision :: densA, densB

    double precision :: phiA, phiB, phiA_eq, phiB_eq, phiA_neq, phiB_neq
    double precision :: phiA_kneq, phiB_kneq
    double precision :: phiA_rneq, phiB_rneq, phiA_vneq, phiB_vneq
    double precision :: signA, signB

    double precision :: sign_test, array_signA, array_signB
    double precision :: A_rot, A_vib, A_rv, B_rot, B_vib, B_rv
    double precision :: A_coeff, A_coeff2, B_coeff, B_coeff2
    double precision :: A_r, A_v, B_r, B_v
    double precision :: r_sumA, r_sumB, v_sumA, v_sumB

    integer, dimension(4) :: rl, vl
    integer, dimension(2) :: r_modes, v_modes, r_levels, v_levels
    integer, dimension(2) :: i, j, k
    integer :: same
    integer :: coll, num_colls
    integer :: glA, glB
    integer :: status

    integer :: num_pointsA, num_pointsB, test
    double precision :: affected, total_depl

    num_pointsA = vel_grid(n)%num_points
    num_pointsB = vel_grid(m)%num_points

    ! Initialize flags for energy activation
    rot_flag = .false.
    vib_flag = .false.

    tot_colls = 0

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

    ! Factor coefficient accounts for double counting if both 
    ! collision partners are drawn from the same species
    if( n .eq. m )then
       factor_coeff = dt*one_half
       same = AA
    else
       factor_coeff = dt
       same = AB
    end if

    ! Equilibrium species properties
    dens(1) = properties%dens(n)
    temp(1) = properties%species_temp(n)
!!$
    dens(2) = properties%dens(m)
    temp(2) = properties%species_temp(m)

    ! Internal Structure
    r_modes(1) = molecule(n)%rot_modes
    v_modes(1) = molecule(n)%vib_modes

    r_modes(2) = molecule(m)%rot_modes
    v_modes(2) = molecule(m)%vib_modes

    ! Nonequilibrium density
    eq_dens(1) = cdf_eq(ns,n,same)%eq_dens
    eq_dens(2) = cdf_eq(ns,m,same)%eq_dens
    neq_dens(1) = cdf_neq(n,same)%neq_dens
    neq_dens(2) = cdf_neq(m,same)%neq_dens

    kin_neq_dens(1) = sum( abs( cdf_neq(n,same)%kin_df ) )
    kin_neq_dens(2) = sum( abs( cdf_neq(m,same)%kin_df ) )
    kin_eq_dens(1)  = properties%dens(n)
    kin_eq_dens(2)  = properties%dens(m)

!!$    write(*,*)"Non-equilibrium density = ", neq_dens
!!$
!!$    write(*,*)"total deviation: ",sum( cdf_neq(n,same)%kin_df(:) )

    ! Set levels and allocate energy arrays
    if( r_modes(1) .gt. 0 )then
       rot_flag = .true.
       r_levels(1) = phi(n)%num_rot_levels
       allocate( frA( 1:r_levels(1) ), rotA( 1:r_levels(1) ), STAT=status )
       call allocate_error_check( status, "frA, rotA (VR)" )
       allocate( rotA_pos( 1:r_levels(1) ), rotA_neg( 1:r_levels(1) ), STAT=status )
    end if

    if( r_modes(2) .gt. 0 )then
       rot_flag = .true.
       r_levels(2) = phi(m)%num_rot_levels
       allocate( frB( 1:r_levels(2) ), rotB( 1:r_levels(2) ), STAT=status )
       call allocate_error_check( status, "frB, rotB (VR)" )
       allocate( rotB_pos( 1:r_levels(2) ), rotB_neg( 1:r_levels(2) ), STAT=status )
    end if

    if( v_modes(1) .gt. 0 )then
       vib_flag = .true.
       v_levels(1) = phi(n)%num_vib_levels
       allocate( fvA( 1:v_levels(1) ), vibA( 1:v_levels(1) ), STAT=status )
       call allocate_error_check( status, "fvA, vibA (VR)" )
       allocate( vibA_pos( 1:v_levels(1) ), vibA_neg( 1:v_levels(1) ), STAT=status )
    end if

    if( v_modes(2) .gt. 0 )then
       vib_flag = .true.
       v_levels(2) = phi(m)%num_vib_levels
       allocate( fvB( 1:v_levels(2) ), vibB( 1:v_levels(2) ), STAT=status )
       call allocate_error_check( status, "fvB, vibB (VR)" )
       allocate( vibB_pos( 1:v_levels(2) ), vibB_neg( 1:v_levels(2) ), STAT=status )
    end if

    ! Separate depletion into elastic and inelastic parts (elastic, r-t A, r-t B, v-t A, v-t B)
    call split_depletion( depl_frac, m_red, properties, molecule, n, m )

    !======================================================================================================
    ! First collision integral: phi_eq*phi_noneq
    !                           A-B, not B-A
    !======================================================================================================
    ! Calculate the number of collisions
    ! TODO: the maximum number of collisions is bounded by 2(n^2)
    densA = eq_dens(1)/dens(1)
    densB = neq_dens(2)/dens(2)
    call compute_num_coll_pairs( num_colls, dt, densA, densB, &
         temp, mass, coln_rms, vel_grid, n, m )

    ! Number of collisions can't be zero
    if( num_colls .eq. 0 ) num_colls = 1
    tot_colls = tot_colls + num_colls

    ! Calculate factor for depletion
    factor = factor_coeff * sigma * eq_dens(1) * neq_dens(2) / dble(num_colls)

    do coll = 1, num_colls
       ! Pick velocities to collide
       call pick_collision_partner( i(1), j(1), k(1), glA, cdf_eq(ns,n,same)%cumul_df, &
            cdf_eq(ns,n,same)%search_ref, vel_grid(n) )
       call pick_collision_partner( i(2), j(2), k(2), glB, cdf_neq(m,same)%cumul_df, &
            cdf_neq(m,same)%search_ref, vel_grid(m) )

       ! Calculate relative velocity between chosen collision partners
       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )

       ! Check for self collision
       if( g .lt. double_tol ) cycle

       ! Calculate g*sigma term
       g_sigma = g**vhs_exponent

       ! Distribution function values needed for calculations:
       phiA = phi_old(n)%value( i(1), j(1), k(1) ) 
       phiB = phi_old(m)%value( i(2), j(2), k(2) ) 

       phiB_eq   = cdf_eq(ns,m,same)%cumul_df( glB ) - cdf_eq(ns,m,same)%cumul_df( glB-1 )

       phiB_kneq = cdf_neq(m,same)%kin_df( glB )
       phiB_rneq = cdf_neq(m,same)%rot_df( glB )
       phiB_vneq = cdf_neq(m,same)%vib_df( glB )

       signB     = cdf_neq(m,same)%sign( glB )

       ! Equilibrium collision partner
       sumA = one

       if( r_modes(1) .gt. 0 )then
          frA = cdf_eq(ns,n,same)%rot_eq(:)

          sumA    = sum( frA )
          sum_vtA = sumA

          call pick_level( rl(1), frA, one, r_levels(1) )
          depl_sign(2) = one
       end if

       if( v_modes(1) .gt. 0 )then
          fvA = cdf_eq(ns,n,same)%vib_eq(:)

          sumA    = sum( fvA )
          sum_rtA = sumA

          call pick_level( vl(1), fvA, one, v_levels(1) )
          depl_sign(4) = one
       end if

       ! Deviation collision partner
       coeff = cdf_neq(m,same)%coeff( glB )

       if( r_modes(2) .gt. 0 .and. v_modes(2) .gt. 0 )then
          frB = ( phi_old(m)%rot(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%rot_eq * phiB_eq ) / phiB_rneq
          fvB = ( phi_old(m)%vib(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%vib_eq * phiB_eq ) / phiB_vneq
          
          ! Elastic array - used for when particle does not exchange energy
          frB  = &
               ( phiB_eq * phiB_rneq * coeff )   * frB + &
               ( phiB_eq * phiB_vneq * coeff )   * ( phiB_kneq / phiB_vneq ) * cdf_eq(ns,m,same)%rot_eq + &
               ( phiB_vneq * phiB_rneq * coeff ) * abs( phiB_kneq / phiB_vneq ) * frB
          fvB  = &
               ( phiB_eq * phiB_rneq * coeff )   * ( phiB_kneq / phiB_rneq ) * cdf_eq(ns,m,same)%vib_eq + &
               ( phiB_eq * phiB_vneq * coeff )   * fvB + &
               ( phiB_vneq * phiB_rneq * coeff ) * abs( phiB_kneq / phiB_rneq ) * fvB

          call pick_level( rl(2), frB, one, r_levels(2) )
          call pick_level( vl(2), fvB, one, v_levels(2) )

          depl_sign(3) = dsgn( frB( rl(2) ) )
          depl_sign(5) = dsgn( fvB( vl(2) ) )

          ! Inelastic array - used for when particle undergoes r-t or v-t exchange
          rotB = &
               ( phiB_eq * phiB_vneq * coeff )       * depl_sign(5) * cdf_eq(ns,m,same)%rot_eq + &
               ( one - phiB_eq * phiB_vneq * coeff ) * frB
          vibB = &
               ( phiB_eq * phiB_rneq * coeff )       * depl_sign(3) * cdf_eq(ns,m,same)%vib_eq + &
               ( one - phiB_eq * phiB_rneq * coeff ) * fvB

          sumB    = sum( fvB )
          sum_vtB = sum( rotB )
          sum_rtB = sum( vibB )

       else if( r_modes(2) .gt. 0 )then
          frB = ( phi_old(m)%rot(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%rot_eq * phiB_eq ) / phiB_rneq

          sumB    = sum( frB )
          sum_vtB = sumB

          call pick_level( rl(2), frB, one, r_levels(2) )

          depl_sign(3) = dsgn( frB( rl(2) ) )

       else if( v_modes(2) .gt. 0 )then
          fvB = ( phi_old(m)%vib(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%vib_eq * phiB_eq ) / phiB_vneq

          sumB    = sum( fvB )
          sum_rtB = sumB

          call pick_level( vl(2), fvB, one, v_levels(2) )

          depl_sign(5) = dsgn( fvB( vl(2) ) )

       else
          sumB = signB

       end if

       !**ELASTIC*********************************************************************************************
       deplA = sumB
       deplB = sumA

       depletion(1) = deplA        * depl_frac(1) * factor * g_sigma
       depletion(2) = deplB        * depl_frac(1) * factor * g_sigma
       kin_depl(1)  = deplA * sumA * depl_frac(1) * factor * g_sigma
       kin_depl(2)  = deplB * sumB * depl_frac(1) * factor * g_sigma

       call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
            i(1), j(1), k(1), 0, frA, fvA, r_modes(1), v_modes(1) )
       call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
            i(2), j(2), k(2), 0, frB, fvB, r_modes(2), v_modes(2) )
       call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
            rl, i, j, k, g, n, m, vel_grid, m_red, elastic, A, molecule )

       !**ROTATIONAL - TRANSLATIONAL**************************************************************************
       if( rot_flag )then ! TODO: this call to rot_flag seems excessive....
          if( r_modes(1) .gt. 0 )then
             deplA = sumB
             deplB = sum_rtA

             depletion(1) = deplA           * depl_frac(2) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(2) * factor * g_sigma
             kin_depl(1)  = deplA * sum_rtA * depl_frac(2) * factor * g_sigma
             kin_depl(2)  = deplB * sumB    * depl_frac(2) * factor * g_sigma

             call vr_deplete( rot_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, vibA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, vibA, fvB, &
                  rl, i, j, k, g, n, m, vel_grid, m_red, rot_trans, A, molecule )

          end if

          if( r_modes(2) .gt. 0 )then
             deplA = sum_rtB
             deplB = sumA

             depletion(1) = deplA           * depl_frac(3) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(3) * factor * g_sigma
             kin_depl(1)  = deplA * sumA    * depl_frac(3) * factor * g_sigma
             kin_depl(2)  = deplB * sum_rtB * depl_frac(3) * factor * g_sigma

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( rot_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, vibB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, vibB, rl, i, j, k, &
                  g, n, m, vel_grid, m_red, rot_trans, B, molecule )

          end if
       end if

       !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
       if( vib_flag )then
          if( v_modes(1) .gt. 0 )then
             deplA = sumB
             deplB = sum_vtA

             depletion(1) = deplA           * depl_frac(4) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(4) * factor * g_sigma
             kin_depl(1)  = deplA * sum_vtA * depl_frac(4) * factor * g_sigma
             kin_depl(2)  = deplB * sumB    * depl_frac(4) * factor * g_sigma

             call vr_deplete( vib_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), rotA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, rotA, frB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, A, molecule )

          end if

          if( v_modes(2) .gt. 0 )then             
             deplA = sum_vtB
             deplB = sumA

             depletion(1) = deplA           * depl_frac(5) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(5) * factor * g_sigma
             kin_depl(1)  = deplA * sumA    * depl_frac(5) * factor * g_sigma
             kin_depl(2)  = deplB * sum_vtB * depl_frac(5) * factor * g_sigma

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( vib_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), rotB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, rotB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, B, molecule )

          end if
       end if

    end do


    !======================================================================================================
    ! Second collision integral: phi_noneq*phi_eq
    !                            B-A, not A-B
    !======================================================================================================
    ! Calculate the number of collisions
    densA = neq_dens(1)/dens(1)
    densB = eq_dens(2)/dens(2)
    call compute_num_coll_pairs( num_colls, dt, densA, densB, &
         temp, mass, coln_rms, vel_grid, n, m )
    
    ! Number of collisions can't be zero
    if( num_colls .eq. 0 ) num_colls = 1
    tot_colls = tot_colls + num_colls

    ! Calculate factor for depletion
    factor = factor_coeff * sigma * neq_dens(1) * eq_dens(2) / dble(num_colls)

    do coll = 1, num_colls
       ! Pick velocities to collide
       call pick_collision_partner( i(1), j(1), k(1), glA, cdf_neq(n,same)%cumul_df, &
            cdf_neq(n,same)%search_ref, vel_grid(n) )
       call pick_collision_partner( i(2), j(2), k(2), glB, cdf_eq(ns,m,same)%cumul_df, &
            cdf_eq(ns,m,same)%search_ref, vel_grid(m) )

       ! Calculate relative velocity between chosen collision partners
       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )

       ! Check for self collision
       if( g .lt. double_tol ) cycle

       ! Calculate g*sigma term
       g_sigma = g**vhs_exponent

       ! Temporary arrays for phi
       phiA = abs( phi_old(n)%value( i(1), j(1), k(1) ) )
       phiB = abs( phi_old(m)%value( i(2), j(2), k(2) ) )

       phiA_eq  = cdf_eq(ns,n,same)%cumul_df( glA ) - cdf_eq(ns,n,same)%cumul_df( glA-1 )

       phiA_kneq  = cdf_neq(n,same)%kin_df( glA )
       phiA_rneq = cdf_neq(n,same)%rot_df( glA )
       phiA_vneq = cdf_neq(n,same)%vib_df( glA )

       signA = cdf_neq(n,same)%sign( glA )

       depl_sign = one
       sumA = one
       sumB = one
       sum_rtA = one
       sum_vtA = one
       sum_rtB = one
       sum_vtB = one


       ! Deviation collision partner
       coeff = cdf_neq(n,same)%coeff( glA )

       if( r_modes(1) .gt. 0 .and. v_modes(1) .gt. 0 )then
          frA = ( phi_old(m)%rot(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%rot_eq * phiA_eq ) / phiA_rneq
          fvA = ( phi_old(m)%vib(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%vib_eq * phiA_eq ) / phiA_vneq
          
          ! Elastic array - used for when particle does not exchange energy
          frA  = &
               ( phiA_eq * phiA_rneq * coeff )   * frA + &
               ( phiA_eq * phiA_vneq * coeff )   * ( phiA_kneq / phiA_vneq ) * cdf_eq(ns,n,same)%rot_eq + &
               ( phiA_vneq * phiA_rneq * coeff ) * abs( phiA_kneq / phiA_vneq ) * frA
          fvA  = &
               ( phiA_eq * phiA_rneq * coeff )   * ( phiA_kneq / phiA_rneq ) * cdf_eq(ns,n,same)%vib_eq + &
               ( phiA_eq * phiA_vneq * coeff )   * fvA + &
               ( phiA_vneq * phiA_rneq * coeff ) * abs( phiA_kneq / phiA_rneq ) * fvA

          call pick_level( rl(1), frA, one, r_levels(1) )
          call pick_level( vl(1), fvA, one, v_levels(1) )

          depl_sign(2) = dsgn( frA( rl(1) ) )
          depl_sign(4) = dsgn( fvA( vl(1) ) )

          ! Inelastic array - used for when particle undergoes r-t or v-t exchange
          rotA = &
               ( phiA_eq * phiA_vneq * coeff )       * depl_sign(5) * cdf_eq(ns,n,same)%rot_eq + &
               ( one - phiA_eq * phiA_vneq * coeff ) * frA
          vibA = &
               ( phiA_eq * phiA_rneq * coeff )       * depl_sign(3) * cdf_eq(ns,n,same)%vib_eq + &
               ( one - phiA_eq * phiA_rneq * coeff ) * fvA

          sumA    = sum( fvA )
          sum_vtA = sum( rotA )
          sum_rtA = sum( vibA )

       else if( r_modes(1) .gt. 0 )then
          frA = ( phi_old(m)%rot(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%rot_eq * phiA_eq ) / phiA_rneq

          sumA    = sum( frA )
          sum_vtA = sumA

          call pick_level( rl(1), frA, one, r_levels(1) )

          depl_sign(2) = dsgn( frA( rl(1) ) )

       else if( v_modes(1) .gt. 0 )then
          fvA = ( phi_old(m)%vib(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%vib_eq * phiA_eq ) / phiA_vneq

          sumA    = sum( fvA )
          sum_rtA = sumA

          call pick_level( vl(1), fvA, one, v_levels(1) )

          depl_sign(4) = dsgn( fvA( vl(1) ) )

       else
          sumA = signA

       end if

       ! Equilibrium collision partner
       sumB = one

       if( r_modes(2) .gt. 0 )then
          frB = cdf_eq(ns,m,same)%rot_eq(:)

          sumB    = sum( frB )
          sum_vtB = sumB

          call pick_level( rl(2), frB, one, r_levels(2) )
          depl_sign(3) = one
       end if

       if( v_modes(2) .gt. 0 )then
          fvB = cdf_eq(ns,m,same)%vib_eq(:)

          sumB    = sum( fvB )
          sum_rtB = sumB

          call pick_level( vl(2), fvB, one, v_levels(2) )
          depl_sign(5) = one
       end if


       !**ELASTIC*********************************************************************************************
       deplA = sumB
       deplB = sumA

       depletion(1) = deplA        * depl_frac(1) * factor * g_sigma
       depletion(2) = deplB        * depl_frac(1) * factor * g_sigma
       kin_depl(1)  = deplA * sumA * depl_frac(1) * factor * g_sigma
       kin_depl(2)  = deplB * sumB * depl_frac(1) * factor * g_sigma

       call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
            i(1), j(1), k(1), 0, frA, fvA, r_modes(1), v_modes(1) )
       call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
            i(2), j(2), k(2), 0, frB, fvB, r_modes(2), v_modes(2) )
       call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
            rl, i, j, k, g, n, m, vel_grid, m_red, elastic, A, molecule )

       !**ROTATIONAL - TRANSLATIONAL**************************************************************************
       if( rot_flag )then ! TODO: this call to rot_flag seems excessive....
          if( r_modes(1) .gt. 0 )then
             deplA = sumB
             deplB = sum_rtA

             depletion(1) = deplA           * depl_frac(2) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(2) * factor * g_sigma
             kin_depl(1)  = deplA * sum_rtA * depl_frac(2) * factor * g_sigma
             kin_depl(2)  = deplB * sumB    * depl_frac(2) * factor * g_sigma

             call vr_deplete( rot_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, vibA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, vibA, fvB, &
                  rl, i, j, k, g, n, m, vel_grid, m_red, rot_trans, A, molecule )

          end if


          if( r_modes(2) .gt. 0 )then
             deplA = sum_rtB
             deplB = sumA

             depletion(1) = deplA           * depl_frac(3) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(3) * factor * g_sigma
             kin_depl(1)  = deplA * sumA    * depl_frac(3) * factor * g_sigma
             kin_depl(2)  = deplB * sum_rtB * depl_frac(3) * factor * g_sigma

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( rot_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, vibB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, vibB, rl, i, j, k, &
                  g, n, m, vel_grid, m_red, rot_trans, B, molecule )

          end if
       end if

       !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
       if( vib_flag )then
          if( v_modes(1) .gt. 0 )then
             deplA = sumB
             deplB = sum_vtA

             depletion(1) = deplA           * depl_frac(4) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(4) * factor * g_sigma
             kin_depl(1)  = deplA * sum_vtA * depl_frac(4) * factor * g_sigma
             kin_depl(2)  = deplB * sumB    * depl_frac(4) * factor * g_sigma

             call vr_deplete( vib_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), rotA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, rotA, frB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, A, molecule )

          end if

          if( v_modes(2) .gt. 0 )then            
             deplA = sum_vtB
             deplB = sumA

             depletion(1) = deplA           * depl_frac(5) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(5) * factor * g_sigma
             kin_depl(1)  = deplA * sumA    * depl_frac(5) * factor * g_sigma
             kin_depl(2)  = deplB * sum_vtB * depl_frac(5) * factor * g_sigma

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( vib_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), rotB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, rotB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, B, molecule )

          end if
       end if

    end do


    !======================================================================================================
    ! Third collision integral: phi_noneq*phi_noneq
    !                           A-B and B-A
    !======================================================================================================
    ! Calculate the number of collisions
    densA = neq_dens(1)/dens(1)
    densB = neq_dens(2)/dens(2)
    call compute_num_coll_pairs( num_colls, dt, densA, densB, &
         temp, mass, coln_rms, vel_grid, n, m )

    ! Number of collisions can't be zero
    if( num_colls .eq. 0 ) num_colls = 1
    tot_colls = tot_colls + num_colls

    ! Calculate factor for depletion
    factor = factor_coeff * sigma * neq_dens(1) * neq_dens(2) / dble(num_colls)

    do coll = 1, num_colls

       ! Pick velocities to collide
       call pick_collision_partner( i(1), j(1), k(1), glA, cdf_neq(n,same)%cumul_df, &
            cdf_neq(n,same)%search_ref, vel_grid(n) )
       call pick_collision_partner( i(2), j(2), k(2), glB, cdf_neq(m,same)%cumul_df, &
            cdf_neq(m,same)%search_ref, vel_grid(m) )

       ! Calculate relative velocity between chosen collision partners
       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )

       ! Check for self collision
       if( g .lt. double_tol ) cycle

       ! Calculate g*sigma term
       g_sigma = g**vhs_exponent

       ! Temporary arrays for phi
       phiA = abs( phi_old(n)%value( i(1), j(1), k(1) ) )
       phiB = abs( phi_old(m)%value( i(2), j(2), k(2) ) )

       phiA_eq  = cdf_eq(ns,n,same)%cumul_df( glA ) - cdf_eq(ns,n,same)%cumul_df( glA-1 )
       phiB_eq  = cdf_eq(ns,m,same)%cumul_df( glB ) - cdf_eq(ns,m,same)%cumul_df( glB-1 )

       phiA_kneq  = cdf_neq(n,same)%kin_df( glA )
       phiA_rneq = cdf_neq(n,same)%rot_df( glA )
       phiA_vneq = cdf_neq(n,same)%vib_df( glA )

       phiB_kneq = cdf_neq(m,same)%kin_df( glB )
       phiB_rneq = cdf_neq(m,same)%rot_df( glB )
       phiB_vneq = cdf_neq(m,same)%vib_df( glB )

       signA = cdf_neq(n,same)%sign( glA )
       signB = cdf_neq(m,same)%sign( glB )

       depl_sign = one
       sumA = one
       sumB = one
       sum_rtA = one
       sum_vtA = one
       sum_rtB = one
       sum_vtB = one
       
       ! Deviation collision partner - A
       coeff = cdf_neq(n,same)%coeff( glA )

       if( r_modes(1) .gt. 0 .and. v_modes(1) .gt. 0 )then
          frA = ( phi_old(m)%rot(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%rot_eq * phiA_eq ) / phiA_rneq
          fvA = ( phi_old(m)%vib(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%vib_eq * phiA_eq ) / phiA_vneq
          
          ! Elastic array - used for when particle does not exchange energy
          frA  = &
               ( phiA_eq * phiA_rneq * coeff )   * frA + &
               ( phiA_eq * phiA_vneq * coeff )   * ( phiA_kneq / phiA_vneq ) * cdf_eq(ns,n,same)%rot_eq + &
               ( phiA_vneq * phiA_rneq * coeff ) * abs( phiA_kneq / phiA_vneq ) * frA
          fvA  = &
               ( phiA_eq * phiA_rneq * coeff )   * ( phiA_kneq / phiA_rneq ) * cdf_eq(ns,n,same)%vib_eq + &
               ( phiA_eq * phiA_vneq * coeff )   * fvA + &
               ( phiA_vneq * phiA_rneq * coeff ) * abs( phiA_kneq / phiA_rneq ) * fvA

          call pick_level( rl(1), frA, one, r_levels(1) )
          call pick_level( vl(1), fvA, one, v_levels(1) )

          depl_sign(2) = dsgn( frA( rl(1) ) )
          depl_sign(4) = dsgn( fvA( vl(1) ) )

          ! Inelastic array - used for when particle undergoes r-t or v-t exchange
          rotA = &
               ( phiA_eq * phiA_vneq * coeff )       * depl_sign(5) * cdf_eq(ns,n,same)%rot_eq + &
               ( one - phiA_eq * phiA_vneq * coeff ) * frA
          vibA = &
               ( phiA_eq * phiA_rneq * coeff )       * depl_sign(3) * cdf_eq(ns,n,same)%vib_eq + &
               ( one - phiA_eq * phiA_rneq * coeff ) * fvA

          sumA    = sum( fvA )
          sum_vtA = sum( rotA )
          sum_rtA = sum( vibA )

       else if( r_modes(1) .gt. 0 )then
          frA = ( phi_old(m)%rot(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%rot_eq * phiA_eq ) / phiA_rneq

          sumA    = sum( frA )
          sum_vtA = sumA

          call pick_level( rl(1), frA, one, r_levels(1) )

          depl_sign(2) = dsgn( frA( rl(1) ) )

       else if( v_modes(1) .gt. 0 )then
          fvA = ( phi_old(m)%vib(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%vib_eq * phiA_eq ) / phiA_vneq

          sumA    = sum( fvA )
          sum_rtA = sumA

          call pick_level( vl(1), fvA, one, v_levels(1) )

          depl_sign(4) = dsgn( fvA( vl(1) ) )

       else
          sumA = signA

       end if

       ! Deviation collision partner - B
       coeff = cdf_neq(m,same)%coeff( glB )

       if( r_modes(2) .gt. 0 .and. v_modes(2) .gt. 0 )then
          frB = ( phi_old(m)%rot(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%rot_eq * phiB_eq ) / phiB_rneq
          fvB = ( phi_old(m)%vib(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%vib_eq * phiB_eq ) / phiB_vneq
          
          ! Elastic array - used for when particle does not exchange energy
          frB  = &
               ( phiB_eq * phiB_rneq * coeff )   * frB + &
               ( phiB_eq * phiB_vneq * coeff )   * ( phiB_kneq / phiB_vneq ) * cdf_eq(ns,m,same)%rot_eq + &
               ( phiB_vneq * phiB_rneq * coeff ) * abs( phiB_kneq / phiB_vneq ) * frB
          fvB  = &
               ( phiB_eq * phiB_rneq * coeff )   * ( phiB_kneq / phiB_rneq ) * cdf_eq(ns,m,same)%vib_eq + &
               ( phiB_eq * phiB_vneq * coeff )   * fvB + &
               ( phiB_vneq * phiB_rneq * coeff ) * abs( phiB_kneq / phiB_rneq ) * fvB

          call pick_level( rl(2), frB, one, r_levels(2) )
          call pick_level( vl(2), fvB, one, v_levels(2) )

          depl_sign(3) = dsgn( frB( rl(2) ) )
          depl_sign(5) = dsgn( fvB( vl(2) ) )

          ! Inelastic array - used for when particle undergoes r-t or v-t exchange
          rotB = &
               ( phiB_eq * phiB_vneq * coeff )       * depl_sign(5) * cdf_eq(ns,m,same)%rot_eq + &
               ( one - phiB_eq * phiB_vneq * coeff ) * frB
          vibB = &
               ( phiB_eq * phiB_rneq * coeff )       * depl_sign(3) * cdf_eq(ns,m,same)%vib_eq + &
               ( one - phiB_eq * phiB_rneq * coeff ) * fvB

          sumB    = sum( fvB )
          sum_vtB = sum( rotB )
          sum_rtB = sum( vibB )

       else if( r_modes(2) .gt. 0 )then
          frB = ( phi_old(m)%rot(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%rot_eq * phiB_eq ) / phiB_rneq

          sumB    = sum( frB )
          sum_vtB = sumB

          call pick_level( rl(2), frB, one, r_levels(2) )

          depl_sign(3) = dsgn( frB( rl(2) ) )

       else if( v_modes(2) .gt. 0 )then
          fvB = ( phi_old(m)%vib(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%vib_eq * phiB_eq ) / phiB_vneq

          sumB    = sum( fvB )
          sum_rtB = sumB

          call pick_level( vl(2), fvB, one, v_levels(2) )

          depl_sign(5) = dsgn( fvB( vl(2) ) )

       else
          sumB = signB

       end if

       depl_sign(1) = dsgn( sumA ) * dsgn( sumB )
       depl_adj     = sumA * sumB


       !**ELASTIC*********************************************************************************************
       deplA = sumB
       deplB = sumA

       depletion(1) = deplA        * depl_frac(1) * factor * g_sigma
       depletion(2) = deplB        * depl_frac(1) * factor * g_sigma
       kin_depl(1)  = deplA * sumA * depl_frac(1) * factor * g_sigma
       kin_depl(2)  = deplB * sumB * depl_frac(1) * factor * g_sigma

       call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
            i(1), j(1), k(1), 0, frA, fvA, r_modes(1), v_modes(1) )
       call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
            i(2), j(2), k(2), 0, frB, fvB, r_modes(2), v_modes(2) )
       call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
            rl, i, j, k, g, n, m, vel_grid, m_red, elastic, A, molecule )

       !**ROTATIONAL - TRANSLATIONAL**************************************************************************
       if( rot_flag )then ! TODO: this call to rot_flag seems excessive....
          if( r_modes(1) .gt. 0 )then
             deplA = sumB
             deplB = sum_rtA

             depletion(1) = deplA           * depl_frac(2) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(2) * factor * g_sigma
             kin_depl(1)  = deplA * sum_rtA * depl_frac(2) * factor * g_sigma
             kin_depl(2)  = deplB * sumB    * depl_frac(2) * factor * g_sigma

             call vr_deplete( rot_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, vibA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, vibA, fvB, &
                  rl, i, j, k, g, n, m, vel_grid, m_red, rot_trans, A, molecule )

          end if


          if( r_modes(2) .gt. 0 )then
             deplA = sum_rtB
             deplB = sumA

             depletion(1) = deplA           * depl_frac(3) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(3) * factor * g_sigma
             kin_depl(1)  = deplA * sumA    * depl_frac(3) * factor * g_sigma
             kin_depl(2)  = deplB * sum_rtB * depl_frac(3) * factor * g_sigma

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( rot_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, vibB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, vibB, rl, i, j, k, &
                  g, n, m, vel_grid, m_red, rot_trans, B, molecule )

          end if
       end if

       !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
       if( vib_flag )then
          if( v_modes(1) .gt. 0 )then
             deplA = sumB
             deplB = sum_vtA

             depletion(1) = deplA           * depl_frac(4) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(4) * factor * g_sigma
             kin_depl(1)  = deplA * sum_vtA * depl_frac(4) * factor * g_sigma
             kin_depl(2)  = deplB * sumB    * depl_frac(4) * factor * g_sigma

             call vr_deplete( vib_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), rotA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, rotA, frB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, A, molecule )

          end if

          if( v_modes(2) .gt. 0 )then            
             deplA = sum_vtB
             deplB = sumA

             depletion(1) = deplA           * depl_frac(5) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(5) * factor * g_sigma
             kin_depl(1)  = deplA * sumA    * depl_frac(5) * factor * g_sigma
             kin_depl(2)  = deplB * sum_vtB * depl_frac(5) * factor * g_sigma

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( vib_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), rotB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, rotB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, B, molecule )

          end if
       end if

    end do


    ! deallocate energy arrays
    if( r_modes(1) .gt. 0 ) deallocate( frA, rotA )
    if( r_modes(2) .gt. 0 ) deallocate( frB, rotB )
    if( v_modes(1) .gt. 0 ) deallocate( fvA, vibA )
    if( v_modes(2) .gt. 0 ) deallocate( fvB, vibB )

    if( r_modes(1) .gt. 0 ) deallocate( rotA_pos, rotA_neg )
    if( r_modes(2) .gt. 0 ) deallocate( rotB_pos, rotB_neg )
    if( v_modes(1) .gt. 0 ) deallocate( vibA_pos, vibA_neg )
    if( v_modes(2) .gt. 0 ) deallocate( vibB_pos, vibB_neg )

    return
  end subroutine variance_reduction_kernel

  subroutine vr_deplete( coll_type, phi, depletion, kin_depl, i, j, k, level, fr, fv, r_modes, v_modes )

    use DistFunc

    implicit none

    double precision, dimension(:), intent(in) :: fr, fv
    double precision, intent(in) :: depletion, kin_depl

    integer, intent(in) :: coll_type, level, i, j, k, r_modes, v_modes

    type(DistFuncType) :: phi

    select case( coll_type )
    case( vib_trans )
       if( v_modes .gt. 0 )then
          phi%vib( level, i, j, k ) = phi%vib( level, i, j, k ) - kin_depl
       end if
       if( r_modes .gt. 0 )then
          phi%rot( :, i, j, k ) = phi%rot( :, i, j, k ) - depletion*fr
       end if
       phi%value( i, j, k ) = phi%value( i, j, k ) - kin_depl

    case( rot_trans )
       if( r_modes .gt. 0 )then
          phi%rot( level, i, j, k ) = phi%rot( level, i, j, k ) - kin_depl
       end if
       if( v_modes .gt. 0 )then
          phi%vib( :, i, j, k ) = phi%vib( :, i, j, k ) - depletion*fv
       end if
       phi%value( i, j, k ) = phi%value( i, j, k ) - kin_depl

    case( elastic )
       if( v_modes .gt. 0 )then
          phi%vib( :, i, j, k ) = phi%vib( :, i, j, k ) - depletion*fv
       end if
       if( r_modes .gt. 0 )then
          phi%rot( :, i, j, k ) = phi%rot( :, i, j, k ) - depletion*fr
       end if
       phi%value( i, j, k ) = phi%value( i, j, k ) - kin_depl

    case default

    end select

    return
  end subroutine vr_deplete

end module VarianceReduction
