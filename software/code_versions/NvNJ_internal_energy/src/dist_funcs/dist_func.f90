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

     ! Number of energy levels
     integer :: num_rot_levels, num_vib_levels

     ! Density values for each node in 3D velocity space
     double precision, allocatable, dimension(:,:,:) :: value

     double precision, allocatable, dimension(:) :: rot_level, vib_level
     double precision, allocatable, dimension(:,:,:,:) :: rot, vib

     ! Energy level values - columns are rotational states and rows are vibrational states
     double precision, allocatable, dimension(:,:) :: int_energy_level

     ! Energy fractions - (rot,vib,i,j,k)
     double precision, allocatable, dimension(:,:,:,:,:) :: int_energy

  end type DistFuncType
  
  integer :: num_species

  integer, dimension(:), allocatable :: init_rot_df, init_vib_df

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
  public :: compute_maxwellian
  public :: compute_bkw
  public :: set_initial_energy_df_conditions
  public :: compute_rot_distribution
  public :: compute_vib_distribution
  public :: compute_rot_vib_distribution

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
    integer, intent(in) :: species

    type(DistFuncType) :: phi
    type(MoleculeType) :: molecule

    integer :: status
    integer :: i_min, i_max, j_min, j_max, k_min, k_max

    integer :: rot_modes, vib_modes

    integer :: csection_model

    double precision :: mass_ref, diam_ref, temp_ref, dens_ref

    call get_reference_mol_data( diam_ref, mass_ref, temp_ref, dens_ref )
    
    ! Density distribution function

    ! TODO: consider moving somewhere else
    call get_csection_model( csection_model )
    
    select case( csection_model )
    case( hard_sphere )
       molecule%omega = one_half
       molecule%alpha = one

    case( pseudo_maxwell )
       molecule%omega = one
       molecule%alpha = one

    case( vhs, vss )
       ! Leave as initialized

    case default
       write(*,*) "Error: Invalid value of csection_model: ", csection_model
       stop

    end select

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    rot_modes = molecule%rot_modes
    vib_modes = molecule%vib_modes

    allocate( phi%value(i_min:i_max, j_min:j_max, k_min:k_max), STAT=status )
    call allocate_error_check( status, "phi%value" )

    ! Internal Energy distribution function
    if( ( rot_modes .gt. 0 .and. vib_modes .eq. 0 ) .or. ( rot_modes .eq. 0 .and. vib_modes .gt. 0 ) )then
       write(*,*)"ERROR: combined method only works with both active rot_modes and vib_modes"
       stop
    end if

    if( rot_modes .gt. 0 .and. vib_modes .gt. 0 )then
       phi%num_rot_levels = rot_levels( species )
       phi%num_vib_levels = vib_levels( species )

       allocate( phi%rot_level( 1:rot_levels(species) ), STAT=status )
       call allocate_error_check( status, "phi%rot_level" )
       allocate( phi%vib_level( 1:vib_levels(species) ), STAT=status )
       call allocate_error_check( status, "phi%vib_level" )
       allocate( phi%rot( 1:rot_levels(species), i_min:i_max, j_min:j_max, k_min:k_max ), STAT=status )
       call allocate_error_check( status, "phi%rot" )
       allocate( phi%vib( 1:vib_levels(species), i_min:i_max, j_min:j_max, k_min:k_max ), STAT=status )
       call allocate_error_check( status, "phi%vib" )

       allocate( phi%int_energy_level( 1:rot_levels(species), 1:vib_levels(species) ), STAT=status )
       call allocate_error_check( status, "phi%int_energy_level" )
       allocate( phi%int_energy( 1:rot_levels(species), 1:vib_levels(species), &
            i_min:i_max, j_min:j_max, k_min:k_max ), STAT=status )
       call allocate_error_check( status, "phi%int_energy" )
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

    if( molecule%rot_modes .gt. 0 .and. molecule%vib_modes .gt. 0 )then
       deallocate( phi%rot_level, STAT=status )
       call deallocate_error_check( status, "phi%rot_level" )
       deallocate( phi%vib_level, STAT=status )
       call deallocate_error_check( status, "phi%vib_level" )
       deallocate( phi%rot, STAT=status )
       call deallocate_error_check( status, "phi%rot" )
       deallocate( phi%vib, STAT=status )
       call deallocate_error_check( status, "phi%vib" )

       deallocate( phi%int_energy_level, STAT=status )
       call deallocate_error_check( status, "phi%int_energy_level" )
       deallocate( phi%int_energy, STAT=status )
       call deallocate_error_check( status, "phi%int_energy" )
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
    integer :: i, j, k, lr, lv
    
    integer :: r_modes, v_modes, r_levels, v_levels

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

    do k = k_min,k_max
       do j = j_min,j_max
          do i = i_min,i_max

             phi%value(i,j,k) = delta_phi%value(i,j,k)

             if( r_modes .gt. 0 .and. v_modes .gt. 0 )then
                phi%rot(:,i,j,k) = zero
                phi%vib(:,i,j,k) = zero
                do lr = 1, r_levels
                   do lv = 1, v_levels
                      phi%int_energy(lr,lv,i,j,k) = delta_phi%int_energy(lr,lv,i,j,k)
                      phi%rot(lr,i,j,k) = phi%rot(lr,i,j,k) + delta_phi%int_energy(lr,lv,i,j,k)
                      phi%vib(lv,i,j,k) = phi%vib(lv,i,j,k) + delta_phi%int_energy(lr,lv,i,j,k)
                   end do
                end do
             end if

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

    ! This is (x - u)^2. Doing FLOP optimization since this will get called alot.
    Cxsq = (x - u)
    Cxsq = Cxsq*Cxsq

    Cysq = (y - v)
    Cysq = Cysq*Cysq

    Czsq = (z - w)
    Czsq = Czsq*Czsq

    Csq = Cxsq + Cysq + Czsq

    maxwell = maxwell_coeff*exp(-Csq*mass/temp) !TODO: exp is slow, can I speed this up?

    return
  end subroutine compute_maxwellian

  subroutine compute_bkw( x, y, z, density, temp, mass, x_vel, scaled_time, bkw )

    implicit none

    !It is expected that time has been scaled appropriatelyx
    double precision, intent(in) :: x, y, z, density, x_vel
    double precision, intent(in) :: mass, temp, scaled_time
    double precision :: bkw

    double precision :: norm, Cxsq, Cysq, Czsq, Csq
    double precision :: xk
    double precision :: opt

    xk = 1.d0 - 0.4d0*exp(-scaled_time/6.d0)

    !FLOP optimization
    opt = mass/(pi*xk*temp)

    norm = density*sqrt(opt*opt*opt)/( 2.d0*xk )

    Cxsq = ( x - x_vel )*( x - x_vel )
    Cysq = y*y
    Czsq = z*z
    Csq = Cxsq + Cysq + Czsq

    bkw = norm*( 5.d0*xk - 3.d0 + 2.d0*( 1.d0 - xk )*Csq*mass/( xk*temp ) )*&
         exp( -Csq*mass/( xk*temp ) )
    
    return
  end subroutine compute_bkw

  subroutine set_initial_energy_df_conditions( init_rot_df_in, init_vib_df_in )
  
    implicit none

    integer, dimension(:), intent(in) :: init_rot_df_in, init_vib_df_in

    integer :: status

    allocate( init_rot_df(1:num_species), STAT=status )
    call allocate_error_check( status, "init_rot_df" )
    
    allocate( init_vib_df(1:num_species), STAT=status )
    call allocate_error_check( status, "init_vib_df" )

    init_rot_df = init_rot_df_in
    init_vib_df = init_vib_df_in

    return
  end subroutine set_initial_energy_df_conditions

  subroutine compute_rot_distribution( rot_df, rot_levels, molecule, temp, levels, species )

    implicit none

    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: temp
    integer, intent(in) :: levels, species

    double precision, dimension(:) :: rot_df, rot_levels

    double precision :: fraction, energy_level

    integer :: l

    do l = 1, levels
       select case( init_rot_df(species) )
       case( rigid_rotor )
          call compute_rigid_rotor( fraction, energy_level, l, levels, temp, molecule )

       case( nonrigid_rotor )
          call compute_nonrigid_rotor( fraction, energy_level, l, levels, temp, molecule )

       case default
          write(*,*) "Error: Invalid value of init_rot_df. Value is: ", init_rot_df(species)
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
       select case( init_vib_df(species) )
       case( sho )
          call compute_sho( fraction, energy_level, l, levels, temp, molecule )

       case( aho )
          call compute_aho( fraction, energy_level, l, levels, temp, molecule )

       case default
          write(*,*) "Error: Invalid value of init_vib_df. Value is: ", init_vib_df(species)
          stop

       end select

       vib_levels(l) = energy_level
       vib_df(l)     = fraction
    end do

    return
  end subroutine compute_vib_distribution

  subroutine compute_rot_vib_distribution( int_energy_df, int_energy_levels, molecule, &
       temp_r, temp_v, r_levels, v_levels, species )

    implicit none

    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: temp_r, temp_v
    integer, intent(in) :: r_levels, v_levels, species
    
    double precision, dimension(:,:) :: int_energy_df, int_energy_levels

    double precision :: fractionr, fractionv, energy_levelr, energy_levelv

    integer :: lr, lv

    do lr = 1, r_levels
       select case( init_rot_df(species) )
       case( rigid_rotor )
          call compute_rigid_rotor( fractionr, energy_levelr, lr, r_levels, temp_r, molecule )

       case( nonrigid_rotor )
          call compute_nonrigid_rotor( fractionr, energy_levelr, lr, r_levels, temp_r, molecule )

       case default
          write(*,*) "Error: Invalid value of init_rot_df. Value is: ", init_rot_df(species)
          stop

       end select

       do lv = 1, v_levels
          select case( init_vib_df(species) )
          case( sho )
             call compute_sho( fractionv, energy_levelv, lv, v_levels, temp_v, molecule )

          case( aho )
             call compute_aho( fractionv, energy_levelv, lv, v_levels, temp_v, molecule )

          case default
             write(*,*) "Error: Invalid value of init_vib_df. Value is: ", init_vib_df(species)

          end select

          int_energy_levels(lr,lv) = energy_levelr + energy_levelv
          int_energy_df(lr,lv) = fractionr*fractionv
          
       end do
    end do

    return
  end subroutine compute_rot_vib_distribution

  subroutine compute_rigid_rotor( fraction, energy_level, l, levels, temp, molecule )

    implicit none

    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: temp
    integer, intent(in) :: l, levels

    double precision :: fraction, energy_level

    double precision :: theta, q_even, q_odd
    double precision :: denominator, J, factor

    integer :: i

    q_even = molecule%even_spin
    q_odd  = molecule%odd_spin
    theta  = molecule%theta_r

    denominator = zero

    if( temp .gt. ten*theta )then
       denominator = (q_even + q_odd)*temp/(two*theta)

    else
       do i = 1, levels
          J = dble(i-1)
          energy_level = J*( J + one )*theta
          ! NOTE: scaled energy level is (1/2)*J*(J + 1 )*theta but the (1/2) is dropped since it is
          !       always multiplied by 2 in the exponential: exp( -2*e/T )

          if( energy_level/temp .gt. 685 )cycle

          if( mod((i-1),2) .eq. 0 )then
             denominator = denominator + q_even*( two*J + one )*exp( -energy_level/temp )
          else
             denominator = denominator + q_odd*( two*J + one )*exp( -energy_level/temp )
          end if
          !denominator = denominator + ( two*J + one )*exp( -two*energy_level/temp )
       end do

    end if

    J = dble(l-1)
    energy_level = one_half*J*( J + one )*theta

    !fraction = ( two*J + one )*exp( -two*energy_level/temp )*denominator
    
    if( two*energy_level/temp .gt. 685 )then
       fraction = zero

    else
       if( mod((l-1),2) .eq. 0 )then
          fraction = q_even*( two*J + one )*exp( -two*energy_level/temp )/denominator
       else
          fraction = q_odd*( two*J + one )*exp( -two*energy_level/temp )/denominator
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
       energy_level = Y01*J*(J + one) + Y02*J*J*( J + one )*( J + one )

       if( mod((i-1),2) .eq. 0 )then
          denominator = denominator + q_even*( two*J + one )*exp( -energy_level/temp )
       else
          denominator = denominator + q_odd*( two*J + one )*exp( -energy_level/temp )
       end if

    end do

    J = dble(l-1)
    energy_level = Y01*J*(J + one) + Y02*J*J*( J + one )*( J + one )

    if( mod((i-1),2) .eq. 0 )then
       fraction = q_even*( two*J + one )*exp( -two*energy_level/temp )/denominator
    else
       fraction = q_odd*( two*J + one )*exp( -two*energy_level/temp )/denominator
    end if

    return
  end subroutine compute_nonrigid_rotor

  subroutine compute_SHO( fraction, energy_level, l, levels, temp, molecule )

    implicit none

    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: temp
    integer, intent(in) :: l, levels

    double precision :: fraction, energy_level
    double precision :: theta, part_func

!!$    double precision :: denominator
    double precision :: v

    theta = molecule%theta_v
    part_func = one/( one - exp( -theta/temp ) )

!!$    integer :: i

!!$    denominator = zero
!!$
!!$    do i = 1, levels
!!$       v = dble(i-1)
!!$       energy_level = one_half*v*theta
!!$       denominator = denominator + exp( -two*energy_level/temp )
!!$    end do

    v = dble(l-1)
    energy_level = one_half*v*theta

    fraction = exp( -two*energy_level/temp )/part_func!denominator

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
       energy_level = Y10*v + Y20*v*( v + one )

       denominator = denominator + exp( -two*energy_level/temp )
    end do

    v = dble(l-1)
    energy_level = Y10*v + Y20*v*( v + one )

    fraction = exp( -two*energy_level/temp )/denominator

    return
  end subroutine compute_AHO

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

end module DistFunc
