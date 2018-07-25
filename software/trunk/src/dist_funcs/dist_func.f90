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
module DistFunc

  use ErrorCheck
  use Constants
  use SpeciesAndReferenceData

  implicit none

  private

  ! Derived type specification here
  type DistFuncType

     ! Density values for each node in 3D velocity space
     double precision, allocatable, dimension(:,:,:) :: value

     ! Number of energy levels
     integer :: num_rot_levels, num_vib_levels

     ! Energy level values
     double precision, allocatable, dimension(:) :: rot_level
     double precision, allocatable, dimension(:) :: vib_level

     ! Energy values for each energy level at each node in 3D velocity spaces
     double precision, allocatable, dimension(:,:,:,:) :: rot
     double precision, allocatable, dimension(:,:,:,:) :: vib

     ! Single valued rotation
     logical :: single_valued_rot

  end type DistFuncType

  integer :: num_species

  integer, dimension(:), allocatable :: rot_func, vib_func

  character(len=16), allocatable, dimension(:) :: species_codes
  integer, allocatable, dimension(:) :: rot_levels, vib_levels

  integer, parameter :: hard_sphere    = 1
  integer, parameter :: pseudo_maxwell = 2
  integer, parameter :: vhs            = 3
  integer, parameter :: vss            = 4

  integer, parameter :: homonuclear   = 1
  integer, parameter :: heteronuclear = 2

  integer, parameter :: rigid_rotor    = 1
  integer, parameter :: nonrigid_rotor = 2
  integer, parameter :: single_valued_rot = 3
  integer, parameter :: reduced_rot    = 4

  integer, parameter :: sho = 1
  integer, parameter :: aho = 2

  public :: DistFuncType

  ! Declare public data here
  public :: num_species

  ! Declare public subroutines/functions here
  public :: initialize_df_input_arrays
  public :: destroy_df_input_arrays
  public :: set_species_codes
  public :: set_num_energy_levels
  public :: get_num_energy_levels
  public :: create_dist_func
  public :: destroy_dist_func
  public :: perform_dist_func_update
  public :: compute_exact_gaussian
  public :: compute_maxwellian
  public :: compute_bkw
  public :: set_initial_energy_df_conditions
  public :: get_int_energy_df_function
  public :: compute_rot_distribution
  public :: compute_vib_distribution

contains

  subroutine initialize_df_input_arrays( )

    implicit none

    integer :: status

    allocate( rot_levels(1:num_species), STAT=status )
    call allocate_error_check( status, "rot_levels input" )

    allocate( vib_levels(1:num_species), STAT=status )
    call allocate_error_check( status, "vib_levels input" )

    return
  end subroutine initialize_df_input_arrays

  subroutine destroy_df_input_arrays( )

    implicit none

    integer :: status

    deallocate( rot_levels, STAT=status )
    call deallocate_error_check( status, "rot_levels input" )

    deallocate( vib_levels, STAT=status )
    call deallocate_error_check( status, "vib_levels input" )

    return
  end subroutine destroy_df_input_arrays

  subroutine set_num_energy_levels( rot_levels_in, vib_levels_in, species )

    implicit none

    integer, intent(in) :: rot_levels_in, vib_levels_in, species   

    rot_levels(species) = rot_levels_in
    vib_levels(species) = vib_levels_in

    return
  end subroutine set_num_energy_levels

  subroutine get_num_energy_levels( rot_levels_out, vib_levels_out, species )

    implicit none

    integer, intent(in) :: species
    integer :: rot_levels_out, vib_levels_out

    rot_levels_out = rot_levels(species)
    vib_levels_out = vib_levels(species)

    return
  end subroutine get_num_energy_levels

  subroutine create_dist_func( phi, vel_grid, molecule, species )

    use VelocityGrid
    use Scaling

    implicit none

    type(VelocityGridType), intent(in) :: vel_grid
    type(MoleculeType), intent(in) :: molecule
    integer, intent(in) :: species

    type(DistFuncType) :: phi

    integer :: status
    integer :: i_min, i_max, j_min, j_max, k_min, k_max

    integer :: rot_modes, vib_modes

    rot_modes = molecule%rot_modes
    vib_modes = molecule%vib_modes

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    ! Density distribution function
    allocate( phi%value(i_min:i_max, j_min:j_max, k_min:k_max), STAT=status )
    call allocate_error_check( status, "phi%value" )

    ! Internal Energy distribution function
    if( rot_modes .gt. 0 )then
       phi%num_rot_levels = rot_levels( species )

       allocate( phi%rot_level( 1:rot_levels(species) ), STAT=status )
       call allocate_error_check( status, "phi%rot_level" )

       allocate( phi%rot( 1:rot_levels(species), i_min:i_max, j_min:j_max, k_min:k_max ), STAT=status )
       call allocate_error_check( status, "phi%rot" )

       phi%single_valued_rot = .false.
       if( rot_func(species) .eq. single_valued_rot ) phi%single_valued_rot = .true.

    end if

    if( vib_modes .gt. 0 )then
       phi%num_vib_levels = vib_levels( species )

       allocate( phi%vib_level( 1:vib_levels(species) ), STAT=status )
       call allocate_error_check( status, "phi%vib_level" )

       allocate( phi%vib( 1:vib_levels(species), i_min:i_max, j_min:j_max, k_min:k_max ), STAT=status )
       call allocate_error_check( status, "phi%vib" )

    end if

    return
  end subroutine create_dist_func

  subroutine destroy_dist_func( phi, molecule )

    implicit none

    type(MoleculeType), intent(in) :: molecule

    type(DistFuncType) :: phi

    integer :: status

    deallocate( phi%value, STAT=status )
    call deallocate_error_check( status, "phi%value" )

    if( molecule%rot_modes .gt. 0 )then
       deallocate( phi%rot_level, STAT=status )
       call deallocate_error_check( status, "phi%rot_level" )
       deallocate( phi%rot, STAT=status )
       call deallocate_error_check( status, "phi%rot" )
    end if

    if( molecule%vib_modes .gt. 0 )then
       deallocate( phi%vib_level, STAT=status )
       call deallocate_error_check( status, "phi%vib_level" )
       deallocate( phi%vib, STAT=status )
       call deallocate_error_check( status, "phi%vib" )
    end if

    return
  end subroutine destroy_dist_func

  subroutine perform_dist_func_update( phi, delta_phi, molecule, vel_grid )

    use VelocityGrid

    implicit none

    type(DistFuncType), intent(in) :: delta_phi
    type(VelocityGridType), intent(in) :: vel_grid
    type(MoleculeType), intent(in) :: molecule

    type(DistFuncType) :: phi

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: i, j, k, l

    integer :: r_modes, v_modes, r_levels, v_levels
    double precision :: ksum, rsum, vsum ! TODO: check for tdf = sum(rdf) = sum(vdf)

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    r_modes  = molecule%rot_modes
    v_modes  = molecule%vib_modes

    r_levels = phi%num_rot_levels
    v_levels = phi%num_vib_levels

    ksum = zero
    rsum = zero
    vsum = zero

    do k = k_min,k_max
       do j = j_min,j_max
          do i = i_min,i_max

             phi%value(i,j,k) = delta_phi%value(i,j,k)
             ksum = delta_phi%value(i,j,k)

             if( r_modes .gt. 0 )then
                do l = 1, r_levels
                   phi%rot(l,i,j,k) = delta_phi%rot(l,i,j,k)
                end do
                !rsum = sum( delta_phi%rot(:,i,j,k) )
             end if

             if( v_modes .gt. 0 )then
                do l = 1, v_levels
                   phi%vib(l,i,j,k) = delta_phi%vib(l,i,j,k)
                end do
                !vsum = sum( delta_phi%vib(:,i,j,k) )
             end if

!!$             if( abs( ksum - rsum ) .gt. double_tol .or. abs( ksum - vsum ) .gt. double_tol )then
!!$                write(*,*)"Inconsistant mass: ",ksum,rsum,vsum
!!$             end if

          end do
       end do
    end do

    return
  end subroutine perform_dist_func_update

  subroutine compute_maxwellian( x, y, z, maxwell_coeff, mass, u, v, w, temp, maxwell )

    implicit none

    double precision, intent(in) :: x, y, z, maxwell_coeff, mass
    double precision, intent(in) :: u, v, w, temp
    double precision :: maxwell

    double precision :: Cxsq, Cysq, Czsq, Csq

    ! maxwell_coeff = n * ( m / ( pi * T ) )^(3/2)

    ! This is (x - u)^2. Doing FLOP optimization since this will get called alot.
    Cxsq = ( x - u )
    Cxsq = Cxsq * Cxsq

    Cysq = ( y - v )
    Cysq = Cysq * Cysq

    Czsq = ( z - w )
    Czsq = Czsq * Czsq

    Csq = Cxsq + Cysq + Czsq

    maxwell = maxwell_coeff * exp( -Csq * mass / temp ) !TODO: exp is slow, can I speed this up?

    return
  end subroutine compute_maxwellian

  subroutine compute_exact_gaussian( x1, x, x2, y1, y, y2, z1, z, z2, C1, C2, u, v, w, df_value )

    implicit none

    double precision, intent(in) :: x1, x, x2, y1, y, y2, z1, z, z2
    double precision, intent(in) :: C1, C2, u, v, w

    double precision :: df_value

    double precision :: xl, xr, yl, yr, zl, zr
    double precision :: C_x, C_y, C_z

    xl = sqrt(C2) * ( x1 + one_half * ( x - x1 ) - u )
    xr = sqrt(C2) * ( x  + one_half * ( x2 - x ) - u )
    yl = sqrt(C2) * ( y1 + one_half * ( y - y1 ) - v )
    yr = sqrt(C2) * ( y  + one_half * ( y2 - y ) - v )
    zl = sqrt(C2) * ( z1 + one_half * ( z - z1 ) - w )
    zr = sqrt(C2) * ( z  + one_half * ( z2 - z ) - w )

    C_x = erf( xr ) - erf( xl )
    C_y = erf( yl ) - erf( yr )
    C_z = erf( zl ) - erf( zr )
    
    df_value = 0.125d0 * C1 * ( pi / C2 )**(1.5d0) * C_x * C_y * C_z

    return
  end subroutine compute_exact_gaussian

  subroutine compute_bkw( x, y, z, density, temp, mass, u, v, w, xk, norm, bkw )

    implicit none

    !It is expected that time has been scaled appropriatelyx
    double precision, intent(in) :: x, y, z, density, u, v, w
    double precision, intent(in) :: mass, temp, xk, norm
    double precision :: bkw

    double precision :: Cxsq, Cysq, Czsq, Csq

    !xk = 1.d0 - 0.4d0 * exp( -scaled_time / 6.0d0 )

    !FLOP optimization
    !opt = mass / ( pi * xk * temp )

    !norm = density * sqrt( opt * opt * opt ) / ( 2.d0 * xk )

    Cxsq = ( x - u ) * ( x - u )
    Cysq = ( y - v ) * ( y - v )
    Czsq = ( z - w ) * ( z - w )
    Csq = Cxsq + Cysq + Czsq

    bkw = norm * ( 5.d0 * xk - 3.d0 + 2.d0 * ( 1.d0 - xk ) * Csq * mass / ( xk * temp ) ) * &
         exp( -Csq * mass / ( xk * temp ) )

    return
  end subroutine compute_bkw

  subroutine set_initial_energy_df_conditions( rot_func_in, vib_func_in )

    implicit none

    integer, dimension(:), intent(in) :: rot_func_in, vib_func_in

    integer :: status

    allocate( rot_func(1:num_species), STAT=status )
    call allocate_error_check( status, "rot_func" )

    allocate( vib_func(1:num_species), STAT=status )
    call allocate_error_check( status, "vib_func" )

    rot_func = rot_func_in
    vib_func = vib_func_in

    return
  end subroutine set_initial_energy_df_conditions

  subroutine get_int_energy_df_function( rot_func_out, vib_func_out, species )

    implicit none

    integer, intent(in) :: species

    integer :: rot_func_out, vib_func_out

    rot_func_out = rot_func( species )
    vib_func_out = vib_func( species )

    return
  end subroutine get_int_energy_df_function

  subroutine compute_rot_distribution( rot_df, rot_levels, molecule, temp, levels, species )

    implicit none

    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: temp
    integer, intent(in) :: levels, species

    double precision, dimension(:) :: rot_df, rot_levels

    double precision :: fraction, energy_level

    integer :: l

    do l = 1, levels
       select case( rot_func( species ) )
       case( rigid_rotor )
          call compute_rigid_rotor( fraction, energy_level, l, levels, temp, molecule )

       case( nonrigid_rotor )
          call compute_nonrigid_rotor( fraction, energy_level, l, levels, temp, molecule )

       case( reduced_rot )
          call compute_reduced_rot( fraction, energy_level, l, levels, temp, molecule )

       case( single_valued_rot )
          call compute_single_valued_rot( fraction, energy_level, temp, molecule )
          write(*,*)fraction, temp
       case default
          write(*,*) "Error: Invalid value of rot_func. Value is: ", rot_func( species )
          stop

       end select

       rot_levels(l) = energy_level
       rot_df(l)     = fraction
    end do

    return
  end subroutine compute_rot_distribution

  subroutine compute_vib_distribution( vib_df, vib_levels, molecule, temp, levels, species )

    implicit none

    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: temp
    integer, intent(in) :: levels, species

    double precision, dimension(:) :: vib_df, vib_levels

    double precision :: fraction, energy_level

    integer :: l

    do l = 1, levels
       select case( vib_func( species ) )
       case( sho )
          call compute_sho( fraction, energy_level, l, levels, temp, molecule )

       case( aho )
          call compute_aho( fraction, energy_level, l, levels, temp, molecule )

       case default
          write(*,*) "Error: Invalid value of vib_func. Value is: ", vib_func( species )
          stop

       end select

       vib_levels(l) = energy_level
       vib_df(l)     = fraction
    end do

    return
  end subroutine compute_vib_distribution

  subroutine compute_rigid_rotor( fraction, energy_level, l, levels, temp, molecule )

    implicit none

    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: temp
    integer, intent(in) :: l, levels

    double precision :: fraction, energy_level

    double precision :: theta, q_even, q_odd
    double precision :: denominator, J

    integer :: i

    q_even = molecule%even_spin
    q_odd  = molecule%odd_spin
    theta  = molecule%theta_r

    denominator = zero

    do i = 1, levels
       J = dble(i-1)
       energy_level = J * ( J + one ) * theta
       ! NOTE: scaled energy level is (1/2)*J*(J + 1 )*theta but the (1/2) is dropped since it is
       !       always multiplied by 2 in the exponential: exp( -2*e/T )

       if( energy_level/temp .gt. 685 )cycle ! Compiler flag "overflow" doesnt like larger values

       if( mod((i-1),2) .eq. 0 )then
          denominator = denominator + q_even * ( two * J + one ) * exp( -energy_level / temp )
       else
          denominator = denominator + q_odd * ( two * J + one ) * exp( -energy_level / temp )
       end if
    end do

    J = dble(l-1)
    energy_level = one_half * J * ( J + one ) * theta
    
    if( two*energy_level/temp .gt. 685 )then ! Compiler flag "overflow" doesnt like larger values
       fraction = zero

    else
       if( mod((l-1),2) .eq. 0 )then
          fraction = q_even * ( two * J + one ) * exp( -two * energy_level / temp ) / denominator
       else
          fraction = q_odd * ( two * J + one ) * exp( -two * energy_level / temp ) / denominator
       end if

    end if

    return
  end subroutine compute_rigid_rotor

  subroutine compute_nonrigid_rotor( fraction, energy_level, l, levels, temp, molecule )

    implicit none

    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: temp
    integer, intent(in) :: l, levels

    double precision :: fraction, energy_level

    double precision :: q_even, q_odd, Y01, Y02
    double precision :: denominator, J
    
    integer :: i

    denominator = zero

    Y01 = molecule%Y01
    Y02 = molecule%Y02

    q_even = molecule%even_spin
    q_odd  = molecule%odd_spin

    do i = 1, levels
       J = dble(i-1)
       energy_level = Y01 * J * (J + one) + Y02 * J * J * ( J + one ) * ( J + one )

       if( mod((i-1),2) .eq. 0 )then
          denominator = denominator + q_even * ( two * J + one ) * exp( -two * energy_level / temp )
       else
          denominator = denominator + q_odd * ( two * J + one ) * exp( -two * energy_level / temp )
       end if

    end do

    J = dble(l-1)
    energy_level = Y01 * J * ( J + one ) + Y02 * J * J * ( J + one ) * ( J + one )

    if( mod((l-1),2) .eq. 0 )then
       fraction = q_even * ( two * J + one ) * exp( -two * energy_level / temp ) / denominator
    else
       fraction = q_odd * ( two * J + one ) * exp( -two * energy_level / temp ) / denominator
    end if

    return
  end subroutine compute_nonrigid_rotor

  subroutine compute_reduced_rigid_rotor( fraction, energy_level, l, levels, temp, molecule )

    implicit none

    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: temp
    integer, intent(in) :: l, levels

    double precision :: fraction, energy_level

    double precision :: theta, q_even, q_odd
    double precision :: denominator, J

    integer :: i

    !Issue: nuclear spin would need to be ignored?
    !Issue: or the lowest odd level would become the skip number (5 in this case)

    q_even = molecule%even_spin
    q_odd  = molecule%odd_spin
    theta  = molecule%theta_r

    denominator = zero

    do i = 1, levels
       J = 5.0d0*dble(i-1)
       energy_level = J * ( J + one ) * theta
       ! NOTE: scaled energy level is (1/2)*J*(J + 1 )*theta but the (1/2) is dropped since it is
       !       always multiplied by 2 in the exponential: exp( -2*e/T )

       if( energy_level/temp .gt. 685 )cycle ! Compiler flag "overflow" doesnt like larger values

       if( mod((i-1),2) .eq. 0 )then
          denominator = denominator + q_even * ( two * J + one ) * exp( -energy_level / temp )
       else
          denominator = denominator + q_odd * ( two * J + one ) * exp( -energy_level / temp )
       end if
    end do

    J = 5.0d0*dble(l-1)
    energy_level = one_half * J * ( J + one ) * theta
    
    if( two*energy_level/temp .gt. 685 )then ! Compiler flag "overflow" doesnt like larger values
       fraction = zero

    else
       if( mod((l-1),2) .eq. 0 )then
          fraction = q_even * ( two * J + one ) * exp( -two * energy_level / temp ) / denominator
       else
          fraction = q_odd * ( two * J + one ) * exp( -two * energy_level / temp ) / denominator
       end if

    end if

    return
  end subroutine compute_reduced_rigid_rotor

  subroutine compute_reduced_rot( fraction, energy_level, l, levels, temp, molecule )
    
    implicit none

    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: temp
    integer, intent(in) :: l, levels

    double precision :: fraction, energy_level
    double precision :: coeff, denominator

    integer :: J, i

    ! Do I do this? Or do I just skip levels?

    coeff = 0.1d0 * one_half * temp
    denominator = 0.0d0

    do i = 1, levels
       J = dble(i-1)
       energy_level = coeff * J
       denominator = denominator + exp( -two * energy_level / temp )
    end do

    J = dble(l-1)
    energy_level = coeff * J
    fraction = exp( -two * energy_level / temp ) / denominator

  end subroutine compute_reduced_rot

  subroutine compute_SHO( fraction, energy_level, l, levels, temp, molecule )

    implicit none

    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: temp
    integer, intent(in) :: l, levels

    double precision :: fraction, energy_level
    double precision :: theta

    double precision :: denominator
    double precision :: v
    integer :: i

    theta = molecule%theta_v

    denominator = zero

    do i = 1, levels
       v = dble(i-1)
       energy_level = v * theta
       ! NOTE: scaled energy level is (1/2)*v*theta but the (1/2) is dropped since it is
       !       always multiplied by 2 in the exponential: exp( -2*e/T )
       denominator = denominator + exp( -energy_level / temp )
    end do

    v = dble(l-1)
    energy_level = one_half * v * theta

    fraction = exp( -two * energy_level / temp ) / denominator

    return
  end subroutine compute_SHO

  subroutine compute_AHO( fraction, energy_level, l, levels, temp, molecule )

    implicit none
    
    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: temp
    integer, intent(in) :: l, levels

    double precision :: fraction, energy_level

    double precision :: Y10, Y20
    double precision :: denominator, v

    integer :: i

    denominator = zero

    Y10 = molecule%Y10
    Y20 = molecule%Y20

    do i = 1, levels
       v = dble(i-1)
       energy_level = Y10 * v + Y20 * v * ( v + one )

       denominator = denominator + exp( -two * energy_level / temp )
    end do

    v = dble(l-1)
    energy_level = Y10 * v + Y20 * v * ( v + one )

    fraction = exp( -two * energy_level / temp ) / denominator

    return
  end subroutine compute_AHO

  ! TODO: Unused but left in code in case further work requires it
  subroutine discretized_vib_df( fraction, energy_level, l, low, high, temp, molecule )

    implicit none

    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: temp
    integer, intent(in) :: l, low, high

    double precision :: fraction, energy_level

    double precision :: v, Y10, Y20

    integer :: i

    fraction = zero

    Y10 = molecule%Y10
    Y20 = molecule%Y20

    do i = low, high
       v = dble(i-1)
       energy_level = Y10*v + Y20*v*( v + one )
       fraction = fraction + exp( -two*energy_level/temp)
    end do

    v = dble(l-1)
    energy_level = Y10*v + Y20*v*( v + one )
    
    return
  end subroutine discretized_vib_df

  subroutine compute_single_valued_rot( fraction, energy_level, temp, molecule )

    implicit none
    
    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: temp
    
    double precision :: fraction, energy_level

    double precision :: mass, rot_dof

    energy_level = one

    mass = molecule%mass
    call get_species_rot_dof( molecule%species_code, rot_dof )

    fraction = 0.25d0 * temp * rot_dof / mass

    return
  end subroutine compute_single_valued_rot

end module DistFunc
