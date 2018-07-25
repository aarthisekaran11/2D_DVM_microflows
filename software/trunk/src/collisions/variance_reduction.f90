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

  type(CumulativeDF), allocatable, dimension(:,:,:,:) :: cdf_eq
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

    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    type(DistFuncType), dimension(:,:,:), intent(in) :: phi
    type(MoleculeType), dimension(:), intent(in) :: molecule

    integer :: num_points, nx_space, ny_space
    integer :: r_modes, v_modes, r_levels, v_levels
    integer :: nx, ny, same, s
    integer :: status

    call get_nspace( nx_space, ny_space )

    allocate( cdf_eq( 1:nx_space, 1:ny_space, 1:num_species, 1:2 ), STAT=status )
    call allocate_error_check( status, "cdf_eq(VR)" )

    allocate( cdf_neq( 1:num_species, 1:2 ), STAT=status )
    call allocate_error_check( status, "cdf_neq(VR)" )

    do same = 1, 2
       do s = 1, num_species

          num_points = vel_grid(s)%num_points

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

          do nx = 1, nx_space
             do ny = 1, ny_space

                allocate( cdf_eq(nx,ny,s,same)%cumul_df( 0:num_points ), STAT=status )
                call allocate_error_check( status, "cdf_eq%cumul_df(VR)" )

                allocate( cdf_eq(nx,ny,s,same)%eq_df( 0:num_points ), STAT=status )
                call allocate_error_check( status, "cdf_eq%eq_df(VR)" )

                allocate( cdf_eq(nx,ny,s,same)%sign( 0:num_points ), STAT=status )
                call allocate_error_check( status, "cdf_eq%sign(VR)" )     
             
                cdf_eq(nx,ny,s,same)%dens = zero
                cdf_eq(nx,ny,s,same)%u    = zero
                cdf_eq(nx,ny,s,same)%v    = zero
                cdf_eq(nx,ny,s,same)%w    = zero
                cdf_eq(nx,ny,s,same)%temp = zero

                if( r_modes .gt. 0 )then
                   r_levels = phi(s,nx,ny)%num_rot_levels
                   allocate( cdf_eq(nx,ny,s,same)%rot_eq( 1:r_levels ), STAT=status )
                   call allocate_error_check( status, "cdf_eq%rot(VR)" )
                end if

                if( v_modes .gt. 0 )then
                   v_levels = phi(s,nx,ny)%num_vib_levels
                   allocate( cdf_eq(nx,ny,s,same)%vib_eq( 1:v_levels ), STAT=status )
                   call allocate_error_check( status, "cdf_eq%vib(VR)" )
                end if

             end do
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

    integer :: nx_space, ny_space, nx, ny, s, same
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

          do nx = 1, nx_space
             do ny = 1, ny_space

                deallocate( cdf_eq(nx,ny,s,same)%cumul_df, STAT=status )
                call deallocate_error_check( status, "cdf_eq%cumul_df(VR)" )

                deallocate( cdf_eq(nx,ny,s,same)%eq_df, STAT=status )
                call deallocate_error_check( status, "cdf_eq%eq_df(VR)" )

                deallocate( cdf_eq(nx,ny,s,same)%sign, STAT=status )
                call deallocate_error_check( status, "cdf_eq%sign(VR)" )

                if( r_modes .gt. 0 )then
                   deallocate( cdf_eq(nx,ny,s,same)%rot_eq, STAT=status )
                   call deallocate_error_check( status, "cdf_eq%rot(VR)" )
                end if

                if( v_modes .gt. 0 )then
                   deallocate( cdf_eq(nx,ny,s,same)%vib_eq, STAT=status )
                   call deallocate_error_check( status, "cdf_eq%vib(VR)" )
                end if

             end do
          end do
       end do
    end do

    deallocate( cdf_eq, STAT=status )
    call deallocate_error_check( status, "cdf_eq(VR)" )

    deallocate( cdf_neq, STAT=status )
    call deallocate_error_check( status, "cdf_neq(VR)" )

    return
  end subroutine destroy_variance_reduction

  subroutine compute_cdf_vr( phi, molecule, vel_grid, properties, new_eq_flag, same, s, nx, ny )

    use DistFunc
    use VelocityGrid
    use PhysicalProperties

    implicit none

    type(DistFuncType), intent(in) :: phi
    type(MoleculeType), intent(in) :: molecule
    type(VelocityGridType), intent(in) :: vel_grid
    type(PropertiesType), intent(in) :: properties
    integer, intent(in) :: same, s, nx, ny
    logical :: new_eq_flag

    double precision, allocatable, dimension(:) :: temporary

    double precision :: x, y, z
    double precision :: beta_x, beta_y, beta_z, beta3
    double precision :: dens, u, v, w, temp, mass
    double precision :: maxwell_coeff, phi_eq, phi_neq, rot_neq, vib_neq, tot_neq

    double precision :: dens_ratio

    double precision :: div_Nxy, div_Nx
    integer :: num_points_x, num_points_y, num_points
    integer :: l, i, j, k
    integer :: i_min, j_min, k_min
    integer :: r_modes, v_modes, r_levels, v_levels

    integer :: ref
    double precision :: max

    double precision :: eq_magnitude, Lt, phit

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
    div_Nxy = vel_grid%div_Nxy
    div_Nx  = vel_grid%div_Nx

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

     if( dens .gt. zero )then
        maxwell_coeff = dens * sqrt( mass * mass * mass / ( pi * pi * pi * temp * temp * temp ) )
     else
        maxwell_coeff = zero
     end if

     new_eq_flag = .false.
     if( abs( dens - cdf_eq(nx,ny,s,same)%dens ) .gt. cdf_tolerance*dens .or. &
          abs( u - cdf_eq(nx,ny,s,same)%u ) .gt. cdf_tolerance*u .or. &
          abs( v - cdf_eq(nx,ny,s,same)%v ) .gt. cdf_tolerance*v .or. &
          abs( temp - cdf_eq(nx,ny,s,same)%temp ) .gt. cdf_tolerance*temp )then
        new_eq_flag = .true.
     end if

    ! Set initial array values
    cdf_eq(nx,ny,s,same)%cumul_df(0) = zero
    cdf_eq(nx,ny,s,same)%eq_df(0) = zero
    cdf_neq(s,same)%cumul_df(0)  = zero

    cdf_neq(s,same)%kin_df(0) = zero
    cdf_neq(s,same)%rot_df(0) = zero
    cdf_neq(s,same)%vib_df(0) = zero

    cdf_eq(nx,ny,s,same)%sign(0)  = one
    cdf_neq(s,same)%sign(0)   = one

    ! Compute quantities which are independent of velocity
    if( new_eq_flag )then
       ! Coefficient for calculating the Maxwellian
       

       do l = 1, num_points
          
          ! Find x,y,z coordinates of the velocity location
          call global2local_map( l-1, i_min, j_min, k_min, num_points_x, num_points_y, &
               div_Nxy, div_Nx, i, j, k )

          ! Velocity and grid spacing
          x = vel_grid%x(i)
          y = vel_grid%y(j)
          z = vel_grid%z(k)

          beta_x = vel_grid%beta_x(i)
          beta_y = vel_grid%beta_y(j)
          beta_z = vel_grid%beta_z(k)
          beta3  = beta_x * beta_y * beta_z

          if( temp .gt. zero )then
             call compute_maxwellian( x, y, z, maxwell_coeff, mass, u, v, w, temp, phi_eq )
          else
             phi_eq = zero
          end if
          phi_eq = phi_eq * beta3

          cdf_eq(nx,ny,s,same)%eq_df(l) = cdf_eq(nx,ny,s,same)%eq_df(l-1) + phi_eq
          cdf_eq(nx,ny,s,same)%sign(l)  = one

       end do

       ! Adjust equilibrium density for discretiazation error
       !dens_ratio = dens / cdf_eq(nx,ny,s,same)%cumul_df(num_points)

       !cdf_eq(nx,ny,s,same)%cumul_df(0:num_points) = dens_ratio * cdf_eq(nx,ny,s,same)%cumul_df(0:num_points)

       cdf_eq(nx,ny,s,same)%dens = dens
       cdf_eq(nx,ny,s,same)%u    = u
       cdf_eq(nx,ny,s,same)%v    = v
       cdf_eq(nx,ny,s,same)%temp = temp

       ! Calculate the rotational equilibrium distribution function
       if( r_modes .gt. 0 )then
          r_levels  = phi%num_rot_levels
          allocate( temporary(1:r_levels) ) ! Second term in subroutine call must be an array
          call compute_rot_distribution( cdf_eq(nx,ny,s,same)%rot_eq, temporary, molecule, temp, r_levels, s )
          deallocate( temporary )
       end if
       
       ! Calculate the vibrational equilibrium distribution function
       if( v_modes .gt. 0 )then
          v_levels = phi%num_vib_levels
          allocate( temporary(1:v_levels) ) ! Second term in subroutine call must be an array
          call compute_vib_distribution( cdf_eq(nx,ny,s,same)%vib_eq, temporary, molecule, temp, v_levels, s )
          deallocate( temporary )
       end if

    end if

    ! Loop over all the velocity points and calculate the equilibrium and deviation distributions
    do l = 1, num_points

       ! Find x,y,z coordinates of the velocity location
       call global2local_map( l-1, i_min, j_min, k_min, num_points_x, num_points_y, &
            div_Nxy, div_Nx, i, j, k )

       phi_eq  = cdf_eq(nx,ny,s,same)%eq_df(l) - cdf_eq(nx,ny,s,same)%eq_df(l-1)
       phi_neq = phi%value(i,j,k) - phi_eq

       ! Construct cumulative distribution function and record associated sign values
       cdf_neq(s,same)%kin_df(l) = phi_neq ! NOT stored as a cumulative function
       cdf_neq(s,same)%sign(l)   = dsgn( phi_neq )

       ! The internal energy deviation density at each velocity location is the sum of the absolute magnitude 
       ! of each rotational level
       if( r_modes .gt. 0 .and. v_modes .gt. 0 )then
          eq_magnitude = phi_eq
          rot_neq = sum( abs( phi%rot(:,i,j,k) - cdf_eq(nx,ny,s,same)%rot_eq(:) * eq_magnitude ) )
          vib_neq = sum( abs( phi%vib(:,i,j,k) - cdf_eq(nx,ny,s,same)%vib_eq(:) * eq_magnitude ) )

       else if( r_modes .gt. 0 )then
          rot_neq = sum( abs( phi%rot(:,i,j,k) - cdf_eq(nx,ny,s,same)%rot_eq(:) * phi_eq ) )
          vib_neq = zero

       else if( v_modes .gt. 0 )then
          rot_neq = 0
          vib_neq = sum( abs( phi%vib(:,i,j,k) - cdf_eq(nx,ny,s,same)%vib_eq(:) * phi_eq ) )

       else
          rot_neq = zero
          vib_neq = zero

       end if

       cdf_neq(s,same)%rot_df(l) = rot_neq ! NOT stored as a cumulative function
       cdf_neq(s,same)%vib_df(l) = vib_neq ! NOT stored as a cumulative function

       
       ! Calculate the total non-equilibrium density and set coefficients for later?
       if( r_modes .gt. 0 .and. v_modes .gt. 0 )then
          tot_neq = rot_neq + vib_neq
          
          if( tot_neq .gt. zero )then
             cdf_neq(s,same)%coeff(l) = one / ( rot_neq + vib_neq )
          else
             cdf_neq(s,same)%coeff(l) = zero
          end if

       else if( r_modes .gt. 0 )then
          tot_neq = rot_neq
          cdf_neq(s,same)%coeff(l) = one

       else if( v_modes .gt. 0 )then
          tot_neq = vib_neq
          cdf_neq(s,same)%coeff(l) = one

       else
          tot_neq = abs( phi_neq )

       end if

       ! Cumulative distribution function
       cdf_neq(s,same)%cumul_df(l)  = cdf_neq(s,same)%cumul_df(l-1)  + tot_neq
       cdf_eq(nx,ny,s,same)%cumul_df(l) = cdf_eq(nx,ny,s,same)%cumul_df(l-1) + phi_eq

    end do

    ! Deviation density - n_dev
    cdf_neq(s,same)%neq_dens = cdf_neq(s,same)%cumul_df(num_points)

    ! Equilibrium density - n_eq
    cdf_eq(nx,ny,s,same)%eq_dens = cdf_eq(nx,ny,s,same)%cumul_df(num_points)

    l = 1
    max = 0.001d0 * cdf_eq(nx,ny,s,same)%eq_dens
    do ref = 1, 1000
       cdf_eq(nx,ny,s,same)%search_ref(ref,1) = l-1
       do while( dble(ref)*max .gt. cdf_eq(nx,ny,s,same)%cumul_df(l) )
          l = l + 1
          if( l .gt. num_points )exit
       end do
       cdf_eq(nx,ny,s,same)%search_ref(ref,2) = l
    end do
    cdf_eq(nx,ny,s,same)%search_ref(1000,2) = num_points

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

subroutine variance_reduction_kernel( phi, phi_old, tot_colls, n, m, nx, ny, coln_rms, &
       molecule, vel_grid, properties )

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
    integer, intent(in) :: n, m, nx, ny

    type(DistFuncType), dimension(:) :: phi
    integer :: tot_colls

    logical :: rot_flag, vib_flag

    double precision, dimension(:), allocatable :: frA, frB, fvA, fvB

    double precision :: coeff
    double precision :: sumA, sumB

    double precision, dimension(5) :: depl_frac, depl_sign
    double precision, dimension(2) :: mass, diam, omega
    double precision, dimension(2) :: viscosity, temp_visc
    double precision, dimension(2) :: temp
    double precision, dimension(2) :: dens, eq_dens, neq_dens, kin_neq_dens, kin_eq_dens, nc_dens_neq
    double precision, dimension(2) :: depletion, kin_depl
    double precision :: m_red, omega_AB, vhs_exponent, sigma
    double precision :: dt, factor_coeff
    double precision :: factor
    double precision :: g, g_sigma
    double precision :: Kn

    double precision :: densA, densB

    double precision :: phiA, phiB, phiA_eq, phiB_eq, phiA_neq, phiB_neq
    double precision :: phiA_kneq, phiB_kneq
    double precision :: phiA_rneq, phiB_rneq, phiA_vneq, phiB_vneq
    double precision :: signA, signB

    integer, dimension(4) :: rl, vl
    integer, dimension(2) :: r_modes, v_modes, r_levels, v_levels
    integer, dimension(2) :: i, j, k
    integer :: same
    integer :: coll, num_colls
    integer :: glA, glB
    integer :: status

    integer :: num_pointsA, num_pointsB

    double precision :: eq_magA, eq_magB
    double precision :: neq_magA, neq_magB

    double precision :: P1, P2, P3, P4, Lt

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
    viscosity(1) = molecule(n)%viscosity
    temp_visc(1) = molecule(n)%T_viscosity

    mass(2)  = molecule(m)%mass
    diam(2)  = molecule(m)%diam
    omega(2) = molecule(m)%omega
    viscosity(2) = molecule(m)%viscosity
    temp_visc(2) = molecule(m)%T_viscosity

    ! Species combination constants
    m_red = mass(1)*mass(2)/( mass(1) + mass(2) )

    omega_AB = one_half*( omega(1) + omega(2) )
    call calculate_cross_section( sigma, m_red, omega_AB, diam, viscosity, temp_visc )

    vhs_exponent = two - two*omega_AB

    ! Knudsen number and time step
    call get_knudsen_number( Kn )
    call get_deltat( dt )

    ! Factor coefficient accounts for double counting if both 
    ! collision partners are drawn from the same species
    if( n .eq. m )then
       factor_coeff = dt * one_half / Kn
       same = AA
    else
       factor_coeff = dt / Kn
       same = AB
    end if

    ! Equilibrium species properties
    dens(1) = properties%dens(n)
    temp(1) = properties%species_temp(n)

    dens(2) = properties%dens(m)
    temp(2) = properties%species_temp(m)

    ! Internal Structure
    r_modes(1) = molecule(n)%rot_modes
    v_modes(1) = molecule(n)%vib_modes

    r_modes(2) = molecule(m)%rot_modes
    v_modes(2) = molecule(m)%vib_modes

    ! Nonequilibrium density
    eq_dens(1) = cdf_eq(nx,ny,n,same)%eq_dens
    eq_dens(2) = cdf_eq(nx,ny,m,same)%eq_dens
    neq_dens(1) = cdf_neq(n,same)%neq_dens
    neq_dens(2) = cdf_neq(m,same)%neq_dens

    kin_neq_dens(1) = sum( abs( cdf_neq(n,same)%kin_df ) )
    kin_neq_dens(2) = sum( abs( cdf_neq(m,same)%kin_df ) )
    kin_eq_dens(1)  = properties%dens(n)
    kin_eq_dens(2)  = properties%dens(m)

    nc_dens_neq(1) = neq_dens(1)
    nc_dens_neq(2) = neq_dens(2)
!!$    write(*,*)neq_dens(1), neq_dens(2), dens(1), dens(2), kin_neq_dens(1), kin_neq_dens(2)
    ! Set levels and allocate energy arrays
    if( r_modes(1) .gt. 0 )then
       rot_flag = .true.
       r_levels(1) = phi(n)%num_rot_levels
       allocate( frA( 1:r_levels(1) ), STAT=status )
       call allocate_error_check( status, "frA(VR)" )
    end if

    if( r_modes(2) .gt. 0 )then
       rot_flag = .true.
       r_levels(2) = phi(m)%num_rot_levels
       allocate( frB( 1:r_levels(2) ), STAT=status )
       call allocate_error_check( status, "frB(VR)" )
    end if

    if( v_modes(1) .gt. 0 )then
       vib_flag = .true.
       v_levels(1) = phi(n)%num_vib_levels
       allocate( fvA( 1:v_levels(1) ), STAT=status )
       call allocate_error_check( status, "fvA(VR)" )
    end if

    if( v_modes(2) .gt. 0 )then
       vib_flag = .true.
       v_levels(2) = phi(m)%num_vib_levels
       allocate( fvB( 1:v_levels(2) ), STAT=status )
       call allocate_error_check( status, "fvB(VR)" )
    end if

    ! Separate depletion into elastic and inelastic parts (elastic, r-t A, r-t B, v-t A, v-t B)
    call split_depletion( depl_frac, m_red, properties, molecule, n, m )
    !===============================================================================================
    ! First collision integral: phi_eq*phi_noneq
    !                           A-B, not B-A
    !===============================================================================================
    ! Calculate the number of collisions
    ! TODO: the maximum number of collisions is bounded by 2(n^2)
    densA = one
    densB = one_half * neq_dens(2) / dens(2)
    call compute_num_coll_pairs( num_colls, dt, Kn, densA, densB, &
         temp, mass, coln_rms, vel_grid, n, m )
!!$    write(*,*)num_colls
    ! Number of collisions can't be zero
    if( num_colls .gt. 0 )then
       ! Calculate factor for depletion
       factor = factor_coeff * sigma * eq_dens(1) * neq_dens(2) / dble(num_colls)
    else
       factor = zero
    end if

    tot_colls = tot_colls + num_colls

    do coll = 1, num_colls
       ! Pick velocities to collide
       call pick_collision_partner( i(1), j(1), k(1), glA, cdf_eq(nx,ny,n,same)%cumul_df, &
            cdf_eq(nx,ny,n,same)%search_ref, vel_grid(n) )
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

       phiB_eq = cdf_eq(nx,ny,m,same)%eq_df( glB ) - cdf_eq(nx,ny,m,same)%eq_df( glB-1 )

       phiB_kneq = cdf_neq(m,same)%kin_df( glB )
       
       if( r_modes(2) .gt. 0 ) phiB_rneq = cdf_neq(m,same)%rot_df( glB )
       if( v_modes(2) .gt. 0 ) phiB_vneq = cdf_neq(m,same)%vib_df( glB )

       phiB_neq = phiB_rneq + phiB_vneq

       signB     = cdf_neq(m,same)%sign( glB )

       ! Equilibrium collision partner
       if( r_modes(1) .gt. 0 .and. v_modes(1) .gt. 0 )then
          frA = cdf_eq(nx,ny,n,same)%rot_eq(:)
          fvA = cdf_eq(nx,ny,n,same)%vib_eq(:)

          sumA = one

          call pick_level( rl(1), frA, r_levels(1) )
          call pick_level( vl(3), fvA, v_levels(1) )

          call pick_level( vl(1), fvA, v_levels(1) )
          call pick_level( rl(3), frA, r_levels(1) )

          depl_sign(2) = one
          depl_sign(4) = one

       else if( r_modes(1) .gt. 0 )then
          frA = cdf_eq(nx,ny,n,same)%rot_eq(:)

          sumA = one

          call pick_level( rl(1), frA, r_levels(1) )

          depl_sign(2) = one

       else if( v_modes(1) .gt. 0 )then
          fvA = cdf_eq(nx,ny,n,same)%vib_eq(:)

          sumA = one

          call pick_level( vl(1), fvA, v_levels(1) )

          depl_sign(4) = one

       else
          sumA = one

       end if

       ! New Idea:
       ! Deviation collision partner:
       !    frB = P(rd,vE)*sum(fvE)*frd + P(rE,vd)*sum(fvd)*frE + P(rd,vd)*sum(fvd)*frd
       !    fvB = P(rd,vE)*sum(frd)*fvE + P(rE,vd)*sum(frE)*fvd + P(rd,vd)*sum(frd)*fvd
       !    rlevel - pick from |frB|, sign[frB(rlevel)]
       !    vlevel - pick from |fvB|, sign[fvB(vlevel)]
       !
       !    Elastic: dPhi*frB and dPhi*fvB, dPhi_B = dPhi*sum(frB) = dPhi*sum(fvB)
       !    R-T:     dPhi_B = dPhi*( P(rd,vE)*sum(fvE) + P(rE,vd)*sum(fvd) + P(rd,vd)*sum(fvd) )
       !    V-T:     dPhi_B = dPhi*( P(rd,vE)*sum(frd) + P(rE,vd)*sum(frE) + P(rd,vd)*sum(frd) )
       !

       ! Deviation collision partner
       coeff = cdf_neq(m,same)%coeff( glB )

       if( r_modes(2) .gt. 0 .and. v_modes(2) .gt. 0 )then
          frB = ( phi_old(m)%rot(:,i(2),j(2),k(2)) - cdf_eq(nx,ny,m,same)%rot_eq * phiB_eq ) / phiB_rneq
          fvB = ( phi_old(m)%vib(:,i(2),j(2),k(2)) - cdf_eq(nx,ny,m,same)%vib_eq * phiB_eq ) / phiB_vneq

          sumB = phiB_kneq / phiB_neq

          call pick_level( rl(2), frB, r_levels(2) )
          call pick_level( vl(4), fvB, v_levels(2) )

          call pick_level( vl(2), fvB, v_levels(2) )
          call pick_level( rl(4), frB, r_levels(2) )

          frB = frB * phiB_rneq / phiB_neq
          fvB = fvB * phiB_vneq / phiB_neq

          depl_sign(3) = ( phiB_rneq / phiB_neq ) * dsgn( frB( rl(2) ) )
          depl_sign(5) = ( phiB_vneq / phiB_neq ) * dsgn( fvB( vl(2) ) )

       else if( r_modes(2) .gt. 0 )then
          frB = ( phi_old(m)%rot(:,i(2),j(2),k(2)) - cdf_eq(nx,ny,m,same)%rot_eq * phiB_eq ) / phiB_rneq

          call pick_level( rl(2), frB, r_levels(2) )

          depl_sign(3) = dsgn( frB( rl(2) ) )

          sumB = sum( frB )

       else if( v_modes(2) .gt. 0 )then
          fvB = ( phi_old(m)%vib(:,i(2),j(2),k(2)) - cdf_eq(nx,ny,m,same)%vib_eq * phiB_eq ) / phiB_vneq

          call pick_level( vl(2), fvB, v_levels(2) )

          depl_sign(5) = dsgn( fvB( vl(2) ) )

          sumB = sum( fvB )

       else
          sumB = signB

       end if

       call vr_depletion_routine( sumA, sumB, depl_frac, depl_sign, &
            phi, phi_old, frA, frB, fvA, fvB, rl, vl, i, j, k, g, n, m, vel_grid, m_red, molecule, &
            r_modes, v_modes, factor, g_sigma, rot_flag, vib_flag )

    end do

    !===============================================================================================
    ! Second collision integral: phi_noneq*phi_eq
    !                            B-A, not A-B
    !===============================================================================================
    ! Calculate the number of collisions
    densA = one_half * neq_dens(1) / dens(1)
    densB = one
    call compute_num_coll_pairs( num_colls, dt, Kn, densA, densB, &
         temp, mass, coln_rms, vel_grid, n, m )
!!$    write(*,*)num_colls
    ! Number of collisions can't be zero
    if( num_colls .gt. 0 )then
       ! Calculate factor for depletion
       factor = factor_coeff * sigma * neq_dens(1) * eq_dens(2) / dble(num_colls)
    else
       factor = zero
    end if

    tot_colls = tot_colls + num_colls

    do coll = 1, num_colls
       ! Pick velocities to collide
       call pick_collision_partner( i(1), j(1), k(1), glA, cdf_neq(n,same)%cumul_df, &
            cdf_neq(n,same)%search_ref, vel_grid(n) )
       call pick_collision_partner( i(2), j(2), k(2), glB, cdf_eq(nx,ny,m,same)%cumul_df, &
            cdf_eq(nx,ny,m,same)%search_ref, vel_grid(m) )

       ! Calculate relative velocity between chosen collision partners
       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )

       ! Check for self collision
       if( g .lt. double_tol ) cycle

       ! Calculate g*sigma term
       g_sigma = g**vhs_exponent

       ! Temporary arrays for phi
       phiA = abs( phi_old(n)%value( i(1), j(1), k(1) ) )
       phiB = abs( phi_old(m)%value( i(2), j(2), k(2) ) )

       phiA_eq = cdf_eq(nx,ny,n,same)%eq_df( glA ) - cdf_eq(nx,ny,n,same)%eq_df( glA-1 )

       phiA_kneq = cdf_neq(n,same)%kin_df( glA )
       if( r_modes(1) .gt. 0 ) phiA_rneq = cdf_neq(n,same)%rot_df( glA )
       if( v_modes(1) .gt. 0 ) phiA_vneq = cdf_neq(n,same)%vib_df( glA )

       phiA_neq = phiA_rneq + phiA_vneq

       signA = cdf_neq(n,same)%sign( glA )

       depl_sign = one
       sumA = one
       sumB = one

       ! Deviation collision partner
       coeff = cdf_neq(n,same)%coeff( glA )

       if( r_modes(1) .gt. 0 .and. v_modes(1) .gt. 0 )then
          frA = ( phi_old(n)%rot(:,i(1),j(1),k(1)) - cdf_eq(nx,ny,n,same)%rot_eq * phiA_eq ) / phiA_rneq
          fvA = ( phi_old(n)%vib(:,i(1),j(1),k(1)) - cdf_eq(nx,ny,n,same)%vib_eq * phiA_eq ) / phiA_vneq

          sumA = phiA_kneq / phiA_neq

          call pick_level( rl(1), frA, r_levels(1) )
          call pick_level( vl(3), fvA, v_levels(1) )

          call pick_level( vl(1), fvA, v_levels(1) )
          call pick_level( rl(3), frA, r_levels(1) )

          frA = frA * phiA_rneq / phiA_neq
          fvA = fvA * phiA_vneq / phiA_neq

          depl_sign(2) = ( phiA_rneq / phiA_neq ) * dsgn( frA( rl(1) ) )
          depl_sign(4) = ( phiA_vneq / phiA_neq ) * dsgn( fvA( vl(1) ) )

       else if( r_modes(1) .gt. 0 )then
          frA = ( phi_old(n)%rot(:,i(1),j(1),k(1)) - cdf_eq(nx,ny,n,same)%rot_eq * phiA_eq ) / phiA_rneq

          call pick_level( rl(1), frA, r_levels(1) )

          depl_sign(2) = dsgn( frA( rl(1) ) )

          sumA = sum( frA )

       else if( v_modes(1) .gt. 0 )then
          fvA = ( phi_old(n)%vib(:,i(1),j(1),k(1)) - cdf_eq(nx,ny,n,same)%vib_eq * phiA_eq ) / phiA_vneq

          call pick_level( vl(1), fvA, v_levels(1) )

          depl_sign(4) = dsgn( fvA( vl(1) ) )

          sumA = sum( fvA )

       else
          sumA = signA

       end if

       ! Equilibrium collision partner
       sumB = one

       if( r_modes(2) .gt. 0 .and. v_modes(2) .gt. 0 )then
          frB = cdf_eq(nx,ny,m,same)%rot_eq(:)
          fvB = cdf_eq(nx,ny,m,same)%vib_eq(:)

          sumB = one

          call pick_level( rl(2), frB, r_levels(2) )
          call pick_level( vl(4), fvB, v_levels(2) )

          call pick_level( vl(2), fvB, v_levels(2) )
          call pick_level( rl(4), frB, r_levels(2) )

          depl_sign(3) = one
          depl_sign(5) = one

       else if( r_modes(2) .gt. 0 )then
          frB = cdf_eq(nx,ny,m,same)%rot_eq(:)

          sumB = one

          call pick_level( rl(2), frB, r_levels(2) )

          depl_sign(3) = one

       else if( v_modes(2) .gt. 0 )then
          fvB = cdf_eq(nx,ny,m,same)%vib_eq(:)

          sumB = one

          call pick_level( vl(2), fvB, v_levels(2) )

          depl_sign(5) = one

       else
          sumB = one

       end if

       call vr_depletion_routine( sumA, sumB, depl_frac, depl_sign, &
            phi, phi_old, frA, frB, fvA, fvB, rl, vl, i, j, k, g, n, m, vel_grid, m_red, molecule, &
            r_modes, v_modes, factor, g_sigma, rot_flag, vib_flag )

    end do

    !================================================================================================
    ! Third collision integral: phi_noneq*phi_noneq
    !                           A-B and B-A
    !================================================================================================
    ! Calculate the number of collisions
    densA = one_half * neq_dens(1) / dens(1)
    densB = one_half * neq_dens(2) / dens(2)
    call compute_num_coll_pairs( num_colls, dt, Kn, densA, densB, &
         temp, mass, coln_rms, vel_grid, n, m )
!!$    write(*,*)num_colls
    ! Number of collisions can't be zero
    if( num_colls .gt. 0 )then
       ! Calculate factor for depletion
       factor = factor_coeff * sigma * neq_dens(1) * neq_dens(2) / dble(num_colls)
    else
       factor = zero
    end if
       
    tot_colls = tot_colls + num_colls

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

       phiA_eq = cdf_eq(nx,ny,n,same)%eq_df( glA ) - cdf_eq(nx,ny,n,same)%eq_df( glA-1 )
       phiB_eq = cdf_eq(nx,ny,m,same)%eq_df( glB ) - cdf_eq(nx,ny,m,same)%eq_df( glB-1 )

       phiA_kneq = cdf_neq(n,same)%kin_df( glA )
       if( r_modes(1) .gt. 0 ) phiA_rneq = cdf_neq(n,same)%rot_df( glA )
       if( v_modes(1) .gt. 0 ) phiA_vneq = cdf_neq(n,same)%vib_df( glA )

       phiB_kneq = cdf_neq(m,same)%kin_df( glB )
       if( r_modes(2) .gt. 0 ) phiB_rneq = cdf_neq(m,same)%rot_df( glB )
       if( v_modes(2) .gt. 0 ) phiB_vneq = cdf_neq(m,same)%vib_df( glB )

       phiA_neq = phiA_rneq + phiA_vneq
       phiB_neq = phiB_rneq + phiB_vneq

       signA = cdf_neq(n,same)%sign( glA )
       signB = cdf_neq(m,same)%sign( glB )
       
       depl_sign = one
       sumA = one
       sumB = one
       
        ! Deviation collision partner - A
       coeff = cdf_neq(n,same)%coeff( glA )

       if( r_modes(1) .gt. 0 .and. v_modes(1) .gt. 0 )then
          frA = ( phi_old(n)%rot(:,i(1),j(1),k(1)) - cdf_eq(nx,ny,n,same)%rot_eq * phiA_eq ) / phiA_rneq
          fvA = ( phi_old(n)%vib(:,i(1),j(1),k(1)) - cdf_eq(nx,ny,n,same)%vib_eq * phiA_eq ) / phiA_vneq
          
          sumA = phiA_kneq / phiA_neq

          call pick_level( rl(1), frA, r_levels(1) )
          call pick_level( vl(3), fvA, v_levels(1) )

          call pick_level( vl(1), fvA, v_levels(1) )
          call pick_level( rl(3), frA, r_levels(1) )

          frA = frA * phiA_rneq / phiA_neq
          fvA = fvA * phiA_vneq / phiA_neq

          depl_sign(2) = ( phiA_rneq / phiA_neq ) * dsgn( frA( rl(1) ) )
          depl_sign(4) = ( phiA_vneq / phiA_neq ) * dsgn( fvA( vl(1) ) )

       else if( r_modes(1) .gt. 0 )then
          frA = ( phi_old(n)%rot(:,i(1),j(1),k(1)) - cdf_eq(nx,ny,n,same)%rot_eq * phiA_eq ) / phiA_rneq

          call pick_level( rl(1), frA, r_levels(1) )

          depl_sign(2) = dsgn( frA( rl(1) ) )

          sumA = sum( frA )

       else if( v_modes(1) .gt. 0 )then
          fvA = ( phi_old(n)%vib(:,i(1),j(1),k(1)) - cdf_eq(nx,ny,n,same)%vib_eq * phiA_eq ) / phiA_vneq

          call pick_level( vl(1), fvA, v_levels(1) )

          depl_sign(4) = dsgn( fvA( vl(1) ) )

          sumA = sum( fvA )

       else
          sumA = signA

       end if

       ! Deviation collision partner - B
       coeff = cdf_neq(m,same)%coeff( glB )

       if( r_modes(2) .gt. 0 .and. v_modes(2) .gt. 0 )then
          frB = ( phi_old(m)%rot(:,i(2),j(2),k(2)) - cdf_eq(nx,ny,m,same)%rot_eq * phiB_eq ) / phiB_rneq
          fvB = ( phi_old(m)%vib(:,i(2),j(2),k(2)) - cdf_eq(nx,ny,m,same)%vib_eq * phiB_eq ) / phiB_vneq

          sumB = phiB_kneq / phiB_neq

          call pick_level( rl(2), frB, r_levels(2) )
          call pick_level( vl(4), fvB, v_levels(2) )

          call pick_level( vl(2), fvB, v_levels(2) )
          call pick_level( rl(4), frB, r_levels(2) )

          frB = frB * phiB_rneq / phiB_neq
          fvB = fvB * phiB_vneq / phiB_neq

          depl_sign(3) = ( phiB_rneq / phiB_neq ) * dsgn( frB( rl(2) ) )
          depl_sign(5) = ( phiB_vneq / phiB_neq ) * dsgn( fvB( vl(2) ) )

       else if( r_modes(2) .gt. 0 )then
          frB = ( phi_old(m)%rot(:,i(2),j(2),k(2)) - cdf_eq(nx,ny,m,same)%rot_eq * phiB_eq ) / phiB_rneq

          call pick_level( rl(2), frB, r_levels(2) )

          depl_sign(3) = dsgn( frB( rl(2) ) )

          sumB = sum( frB )

       else if( v_modes(2) .gt. 0 )then
          fvB = ( phi_old(m)%vib(:,i(2),j(2),k(2)) - cdf_eq(nx,ny,m,same)%vib_eq * phiB_eq ) / phiB_vneq

          call pick_level( vl(2), fvB, v_levels(2) )

          depl_sign(5) = dsgn( fvB( vl(2) ) )

          sumB = sum( fvB )

       else
          sumB = signB

       end if

       call vr_depletion_routine( sumA, sumB, depl_frac, depl_sign, &
            phi, phi_old, frA, frB, fvA, fvB, rl, vl, i, j, k, g, n, m, vel_grid, m_red, molecule, &
            r_modes, v_modes, factor, g_sigma, rot_flag, vib_flag )

    end do

    ! deallocate energy arrays
    if( r_modes(1) .gt. 0 ) deallocate( frA )
    if( r_modes(2) .gt. 0 ) deallocate( frB )
    if( v_modes(1) .gt. 0 ) deallocate( fvA )
    if( v_modes(2) .gt. 0 ) deallocate( fvB )

    return
  end subroutine variance_reduction_kernel

  subroutine vr_depletion_routine( sumA, sumB, depl_frac, depl_sign, &
       phi, phi_old, frA, frB, fvA, fvB, rl, vl, i, j, k, g, n, m, vel_grid, m_red, molecule, &
       r_modes, v_modes, factor, g_sigma, rot_flag, vib_flag )

    use DistFunc
    use VelocityGrid
    use ReplenishingCollisions
    use PhysicalProperties

    implicit none

    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    type(DistFuncType), dimension(:), intent(in) :: phi_old
    type(MoleculeType), dimension(:), intent(in) :: molecule

    logical, intent(in) :: rot_flag, vib_flag

    integer, intent(in) :: n, m
    integer, dimension(2), intent(in) :: i, j, k, r_modes, v_modes
    integer, dimension(4), intent(in) :: rl, vl

    double precision, intent(in) :: sumA, sumB
    double precision, intent(in) :: g, m_red, factor, g_sigma
    double precision, dimension(5), intent(in) :: depl_frac, depl_sign
    double precision, dimension(:), intent(in) :: frA, frB, fvA, fvB

    type(DistFuncType), dimension(:) :: phi

    double precision, dimension(2) :: depletion, kin_depl

    !**ELASTIC*********************************************************************************************
    depletion(1) = sumB * depl_frac(1) * factor * g_sigma
    depletion(2) = sumA * depl_frac(1) * factor * g_sigma
    kin_depl(1)  = sumA * depletion(1)
    kin_depl(2)  = sumB * depletion(2)

    call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
         i(1), j(1), k(1), 0, 0, frA, fvA, r_modes(1), v_modes(1) )
    call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
         i(2), j(2), k(2), 0, 0, frB, fvB, r_modes(2), v_modes(2) )
    call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
         rl, vl, i, j, k, g, n, m, vel_grid, m_red, elastic, A, molecule )

    !**ROTATIONAL - TRANSLATIONAL**************************************************************************
    if( rot_flag )then ! TODO: this call to rot_flag seems excessive....
       if( r_modes(1) .gt. 0 )then
          depletion(1) = sumB * depl_sign(2) * depl_frac(2) * factor * g_sigma
          depletion(2) = one  * depl_sign(2) * depl_frac(2) * factor * g_sigma
          kin_depl(1)  = one  * depletion(1)
          kin_depl(2)  = sumB * depletion(2)

          call vr_deplete( rot_trans, phi(n), depletion(1), kin_depl(1), &
               i(1), j(1), k(1), rl(1), vl(3), frA, fvA, r_modes(1), v_modes(1) )
          call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
               i(2), j(2), k(2), rl(2), vl(4), frB, fvB, r_modes(2), v_modes(2) )
          call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
               rl, vl, i, j, k, g, n, m, vel_grid, m_red, rot_trans, A, molecule )

       end if

       if( r_modes(2) .gt. 0 )then
          depletion(1) = one  * depl_sign(3) * depl_frac(3) * factor * g_sigma
          depletion(2) = sumA * depl_sign(3) * depl_frac(3) * factor * g_sigma
          kin_depl(1)  = sumA * depletion(1)
          kin_depl(2)  = one  * depletion(2)

          call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
               i(1), j(1), k(1), rl(1), vl(3), frA, fvA, r_modes(1), v_modes(1) )
          call vr_deplete( rot_trans, phi(m), depletion(2), kin_depl(2), &
               i(2), j(2), k(2), rl(2), vl(4), frB, fvB, r_modes(2), v_modes(2) )
          call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
               rl, vl, i, j, k, g, n, m, vel_grid, m_red, rot_trans, B, molecule )

       end if
    end if

    !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
    if( vib_flag )then
       if( v_modes(1) .gt. 0 )then
          depletion(1) = sumB * depl_sign(4) * depl_frac(4) * factor * g_sigma
          depletion(2) = one  * depl_sign(4) * depl_frac(4) * factor * g_sigma
          kin_depl(1)  = one  * depletion(1)
          kin_depl(2)  = sumB * depletion(2)

          call vr_deplete( vib_trans, phi(n), depletion(1), kin_depl(1), &
               i(1), j(1), k(1), rl(3), vl(1), frA, fvA, r_modes(1), v_modes(1) )
          call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
               i(2), j(2), k(2), rl(4), vl(2), frB, fvB, r_modes(2), v_modes(2) )
          call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
               rl, vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, A, molecule )

       end if

       if( v_modes(2) .gt. 0 )then
          depletion(1) = one  * depl_sign(5) * depl_frac(5) * factor * g_sigma
          depletion(2) = sumA * depl_sign(5) * depl_frac(5) * factor * g_sigma
          kin_depl(1)  = sumA * depletion(1)
          kin_depl(2)  = one  * depletion(2)

          call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
               i(1), j(1), k(1), rl(3), vl(1), frA, fvA, r_modes(1), v_modes(1) )
          call vr_deplete( vib_trans, phi(m), depletion(2), kin_depl(2), &
               i(2), j(2), k(2), rl(4), vl(2), frB, fvB, r_modes(2), v_modes(2) )
          call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
               rl, vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, B, molecule )

       end if
    end if

    return
  end subroutine vr_depletion_routine

  subroutine vr_deplete( coll_type, phi, depletion, kin_depl, i, j, k, r_level, v_level, fr, fv, r_modes, v_modes )

    use DistFunc

    implicit none

    double precision, dimension(:), intent(in) :: fr, fv
    double precision, intent(in) :: depletion, kin_depl

    integer, intent(in) :: coll_type, r_level, v_level, i, j, k, r_modes, v_modes

    type(DistFuncType) :: phi

    select case( coll_type )
    case( vib_trans )
       if( v_modes .gt. 0 )then
          phi%vib( v_level, i, j, k ) = phi%vib( v_level, i, j, k ) - kin_depl
       end if
       if( r_modes .gt. 0 )then
          phi%rot( r_level, i, j, k ) = phi%rot( r_level, i, j, k ) - kin_depl
       end if
       phi%value( i, j, k ) = phi%value( i, j, k ) - kin_depl

    case( rot_trans )
       if( r_modes .gt. 0 )then
          phi%rot( r_level, i, j, k ) = phi%rot( r_level, i, j, k ) - kin_depl
       end if
       if( v_modes .gt. 0 )then
          phi%vib( v_level, i, j, k ) = phi%vib( v_level, i, j, k ) - kin_depl
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
