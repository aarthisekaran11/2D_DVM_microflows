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
module PhysicalProperties

  use DistFunc
  use VelocityGrid
  use SpeciesAndReferenceData
  use ErrorCheck
  use Constants

  implicit none

  private
  
  type PropertiesType

     ! Mixture Properties
     double precision :: mix_dens, mix_x_vel, mix_y_vel, mix_z_vel
     double precision :: mix_avg_speed, mix_pressure
     double precision :: mix_tr_temp, mix_rot_temp, mix_vib_temp, mix_temp
     double precision :: mix_tr_energy, mix_rot_energy, mix_vib_energy, mix_energy

     ! Species Properties
     double precision, allocatable, dimension(:) :: dens, x_vel, y_vel, z_vel
     double precision, allocatable, dimension(:) :: avg_speed, pressure, coll_freq, non_vib_dof
     double precision, allocatable, dimension(:) :: tr_temp, rot_temp, vib_temp, temp

     ! Energy
     double precision, allocatable, dimension(:) :: tr_energy, rot_energy, vib_energy, total_energy

     ! Degrees of Freedom and Relaxation Rates
     double precision, allocatable, dimension(:) :: vib_dof, rot_dof

     ! VR Species Properties
     double precision, allocatable, dimension(:) :: species_temp

  end type PropertiesType

  integer, parameter :: x_dir = 1
  integer, parameter :: y_dir = 2
  integer, parameter :: z_dir = 3

  public :: PropertiesType
     
  public :: create_physical_properties
  public :: destroy_physical_properties
  public :: calculate_physical_properties

  ! Output property functions
  public :: compute_speed
  public :: compute_mach_number
  public :: compute_directional_temperature
  public :: compute_entropy
  public :: compute_heat_flux
  public :: compute_shear_stress
  public :: compute_moment

contains

  subroutine create_physical_properties( properties )

    implicit none

    type(PropertiesType) :: properties
    integer :: status

    allocate( properties%dens( num_species ), STAT=status )
    call allocate_error_check( status, "properties%dens" )

    allocate( properties%x_vel( num_species ), STAT=status )
    call allocate_error_check( status, "properties%x_vel" )

    allocate( properties%y_vel( num_species ), STAT=status )
    call allocate_error_check( status, "properties%y_vel" )

    allocate( properties%z_vel( num_species ), STAT=status )
    call allocate_error_check( status, "properties%z_vel" )
    
    allocate( properties%tr_temp( num_species ), STAT=status )
    call allocate_error_check( status, "properties%tr_temp" )

    allocate( properties%rot_temp( num_species ), STAT=status )
    call allocate_error_check( status, "properties%rot_temp" )
    
    allocate( properties%vib_temp( num_species ), STAT=status )
    call allocate_error_check( status, "properties%vib_temp" )

    allocate( properties%temp( num_species ), STAT=status )
    call allocate_error_check( status, "properties%temp" )

    allocate( properties%species_temp( num_species ), STAT=status )
    call allocate_error_check( status, "properties%species_temp" )

    allocate( properties%tr_energy( num_species ), STAT=status )
    call allocate_error_check( status, "properties%tr_energy" )

    allocate( properties%rot_energy( num_species ), STAT=status )
    call allocate_error_check( status, "properties%rot_energy" )
    
    allocate( properties%vib_energy( num_species ), STAT=status )
    call allocate_error_check( status, "properties%vib_energy" )

    allocate( properties%total_energy( num_species ), STAT=status )
    call allocate_error_check( status, "properties%total_energy" )

    allocate( properties%rot_dof( num_species ), STAT=status )
    call allocate_error_check( status, "properties%rot_dof" )

    allocate( properties%vib_dof( num_species ), STAT=status )
    call allocate_error_check( status, "properties%vib_dof" )

    allocate( properties%avg_speed( num_species ), STAT=status )
    call allocate_error_check( status, "properties%avg_speed" )

    allocate( properties%pressure( num_species ), STAT=status )
    call allocate_error_check( status, "properties%pressure" )

    allocate( properties%coll_freq( num_species ), STAT=status )
    call allocate_error_check( status, "properties%coll_freq" )

    allocate( properties%non_vib_dof( num_species ), STAT=status )
    call allocate_error_check( status, "properties%non_vib_dof" )

    return
  end subroutine create_physical_properties

  subroutine destroy_physical_properties( properties )

    implicit none

    type(PropertiesType) :: properties
    integer :: status

    deallocate( properties%dens, STAT=status )
    call deallocate_error_check( status, "properties%dens" )

    deallocate( properties%x_vel, STAT=status )
    call deallocate_error_check( status, "properties%x_vel" )

    deallocate( properties%y_vel, STAT=status )
    call deallocate_error_check( status, "properties%y_vel" )

    deallocate( properties%z_vel, STAT=status )
    call deallocate_error_check( status, "properties%z_vel" )
    
    deallocate( properties%tr_temp, STAT=status )
    call deallocate_error_check( status, "properties%tr_temp" )

    deallocate( properties%rot_temp, STAT=status )
    call deallocate_error_check( status, "properties%rot_temp" )
    
    deallocate( properties%vib_temp, STAT=status )
    call deallocate_error_check( status, "properties%vib_temp" )

    deallocate( properties%temp, STAT=status )
    call deallocate_error_check( status, "properties%temp" )

    deallocate( properties%species_temp, STAT=status )
    call deallocate_error_check( status, "properties%species_temp" )

    deallocate( properties%tr_energy, STAT=status )
    call deallocate_error_check( status, "properties%tr_energy" )

    deallocate( properties%rot_energy, STAT=status )
    call deallocate_error_check( status, "properties%rot_energy" )
    
    deallocate( properties%vib_energy, STAT=status )
    call deallocate_error_check( status, "properties%vib_energy" )

    deallocate( properties%total_energy, STAT=status )
    call deallocate_error_check( status, "properties%total_energy" )

    deallocate( properties%rot_dof, STAT=status )
    call deallocate_error_check( status, "properties%rot_dof" )

    deallocate( properties%vib_dof, STAT=status )
    call deallocate_error_check( status, "properties%vib_dof" )

    deallocate( properties%avg_speed, STAT=status )
    call deallocate_error_check( status, "properties%avg_speed" )

    deallocate( properties%pressure, STAT=status )
    call deallocate_error_check( status, "properties%pressure" )

    deallocate( properties%coll_freq, STAT=status )
    call deallocate_error_check( status, "properties%coll_freq" )

    deallocate( properties%non_vib_dof, STAT=status )
    call deallocate_error_check( status, "properties%non_vib_dof" )

    return
  end subroutine destroy_physical_properties

  subroutine calculate_physical_properties( props, phi, molecule, vel_grid )

    implicit none

    type(PropertiesType) :: props

    type(DistFuncType), dimension(:), intent(in) :: phi
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(VelocityGridType), dimension(:), intent(in) :: vel_grid

    integer :: n

    double precision :: x_vel, y_vel, z_vel, avg_speed, coll_freq
    double precision :: dens, u, v, w
    double precision :: temp, tr_energy
    double precision :: macro_dens
    double precision :: rot_dof, rot_energy, rot_temp
    double precision :: vib_dof, vib_energy, vib_temp
    double precision :: total_energy

    ! Calculate individual species density, velocity, and temperature (species referenced)
    do n = 1, num_species

       call compute_density( dens, phi(n), vel_grid(n) )
       call compute_x_velocity( u, dens, phi(n), vel_grid(n) )
       call compute_y_velocity( v, dens, phi(n), vel_grid(n) )
       call compute_z_velocity( w, dens, phi(n), vel_grid(n) )
       call compute_avg_speed( avg_speed, dens, phi(n), vel_grid(n) )
       call compute_temperature( temp, dens, u, v, w, phi(n), molecule(n), vel_grid(n) )
       call compute_coll_freq( coll_freq, dens, temp, molecule(n) )

       props%dens(n)  = dens
       props%x_vel(n) = u
       props%y_vel(n) = v
       props%z_vel(n) = w
       props%species_temp(n) = temp

       props%pressure(n) = dens*temp

       props%avg_speed(n) = avg_speed
       props%coll_freq(n) = coll_freq

    end do

    ! Macro density is the mass density of the mixture
    call compute_macro_density( macro_dens, props%dens, molecule )

    ! Calculate mixture density and velocities
    call compute_total_density( dens, props%dens )
    call compute_total_x_vel( x_vel, props%dens, props%x_vel, macro_dens, molecule )
    call compute_total_y_vel( y_vel, props%dens, props%y_vel, macro_dens, molecule )
    call compute_total_z_vel( z_vel, props%dens, props%z_vel, macro_dens, molecule )
    call compute_total_avg_speed( avg_speed, props%dens, props%avg_speed, macro_dens, molecule )

    props%mix_dens  = dens
    props%mix_x_vel = x_vel
    props%mix_y_vel = y_vel
    props%mix_z_vel = z_vel
    props%mix_avg_speed = avg_speed

    ! Calculate individual species temperature (mixture referenced)
    do n = 1, num_species

       call compute_temperature( temp, props%dens(n), x_vel, y_vel, z_vel, phi(n), molecule(n), vel_grid(n) )
       tr_energy = 0.75d0*props%dens(n)*temp

       props%tr_temp(n) = temp
       props%tr_energy(n)  = tr_energy

    end do

    ! Calculate mixture temperature
    call compute_total_temperature( temp, dens, props%dens, props%tr_temp )
    props%mix_tr_temp = temp
    
    props%mix_pressure = dens*temp

    ! Calculate internal energy related individual species properties
    do n = 1, num_species
       
       call compute_rotational_energy( rot_energy, phi(n), molecule(n), vel_grid(n) )
       call get_species_rot_dof( molecule(n)%species_code, rot_dof )
       call compute_rotational_temperature( rot_temp, rot_energy, rot_dof, molecule(n) )

       call compute_vibrational_energy( vib_energy, phi(n), molecule(n), vel_grid(n) )
       call compute_vibrational_temperature( vib_temp, props%dens(n), vib_energy, molecule(n) )
       call compute_vib_dof( vib_dof, vib_temp, molecule(n) )
!!$       call compute_vibrational_temperature2( vib_temp, props%dens(n), phi(n), molecule(n), vel_grid(n) )
!!$       call compute_vib_dof2( vib_dof, vib_temp, vib_energy, molecule(n) )


       props%rot_energy(n) = rot_energy
       props%rot_temp(n)   = rot_temp
       props%rot_dof(n)    = rot_dof
       
       props%vib_energy(n) = vib_energy
       props%vib_temp(n)   = vib_temp
       props%vib_dof(n)    = vib_dof

       props%non_vib_dof(n) = ( 3.0d0 + rot_dof )/( 3.0d0 + rot_dof + vib_dof )

       props%total_energy(n) = tr_energy + rot_energy + vib_energy

       call compute_overall_temperature( temp, props%tr_energy(n), rot_energy, vib_energy, rot_dof, vib_dof )
       props%temp(n) = temp

       tr_energy = 0.75d0*props%dens(n)*props%species_temp(n)
       call compute_overall_temperature( temp, tr_energy, rot_energy, vib_energy, rot_dof, vib_dof )
       props%species_temp(n) = temp

    end do

    call compute_total_energy( rot_energy, props%rot_energy )
    call compute_total_internal_temperature( rot_temp, props%dens, props%rot_temp, props%rot_dof )
    call compute_total_dof()

    call compute_total_energy( vib_energy, props%vib_energy )
    call compute_total_internal_temperature( vib_temp, props%dens, props%vib_temp, props%vib_dof )
    call compute_total_dof()

    call compute_overall_temperature( temp, 0.75d0*props%mix_tr_temp*props%mix_dens, &
         rot_energy, vib_energy, rot_dof, vib_dof )

    props%mix_rot_energy = rot_energy
    props%mix_rot_temp   = rot_temp

    props%mix_vib_energy = vib_energy
    props%mix_vib_temp   = vib_temp

    props%mix_temp       = temp

    tr_energy = 0.75d0*props%mix_dens*props%mix_tr_temp
    props%mix_tr_energy = tr_energy
    total_energy = tr_energy + rot_energy + vib_energy
    props%mix_energy = total_energy

    return
  end subroutine calculate_physical_properties

  ! Single species properties
  subroutine compute_density( dens, phi, vel_grid )

    implicit none

    type(DistFuncType), intent(in) :: phi
    type(VelocityGridType), intent(in) :: vel_grid
    double precision :: dens

    integer :: i,j,k

    integer :: i_min,i_max,j_min,j_max,k_min,k_max

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    dens = 0.d0
    
    do k = k_min,k_max
       do j = j_min,j_max
          do i = i_min,i_max

             dens = dens + phi%value(i,j,k)

          enddo
       enddo
    enddo

    return
  end subroutine compute_density

  subroutine compute_x_velocity( u, density, phi, vel_grid )

    implicit none

    type(DistFuncType), intent(in) :: phi
    type(VelocityGridType), intent(in) :: vel_grid
    double precision, intent(in) :: density
    double precision :: u

    double precision :: x
    integer :: i,j,k

    integer :: i_min,i_max,j_min,j_max,k_min,k_max

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    u = 0.d0

    do k = k_min,k_max
       do j = j_min,j_max
          do i = i_min,i_max

             x = vel_grid%x(i)

             u = u + phi%value(i,j,k)*x

          enddo
       enddo
    enddo

    u = u/density

    return
  end subroutine compute_x_velocity

  subroutine compute_y_velocity( v, density, phi, vel_grid )

    implicit none

    type(DistFuncType), intent(in) :: phi
    type(VelocityGridType), intent(in) :: vel_grid
    double precision, intent(in) :: density
    double precision :: v

    double precision :: y
    integer :: i,j,k

    integer :: i_min,i_max,j_min,j_max,k_min,k_max

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    v = 0.d0

    do k = k_min,k_max
       do j = j_min,j_max
          y = vel_grid%y(j)

          do i = i_min,i_max

             v = v + phi%value(i,j,k)*y

          enddo
       enddo
    enddo

    v = v/density

    return
  end subroutine compute_y_velocity

  subroutine compute_z_velocity( w, density, phi, vel_grid )    

    implicit none

    type(DistFuncType), intent(in) :: phi
    type(VelocityGridType), intent(in) :: vel_grid
    double precision, intent(in) :: density
    double precision :: w

    double precision :: z
    integer :: i,j,k

    integer :: i_min,i_max,j_min,j_max,k_min,k_max

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    w = 0.d0

    do k = k_min,k_max
       z = vel_grid%z(k)

       do j = j_min,j_max
          do i = i_min,i_max

             w = w + phi%value(i,j,k)*z

          enddo
       enddo
    enddo

    w = w/density

    return
  end subroutine compute_z_velocity

  subroutine compute_temperature( temp, dens, u, v, w, phi, molecule, vel_grid )  

    implicit none

    type(DistFuncType), intent(in) :: phi
    type(MoleculeType), intent(in) :: molecule
    type(VelocityGridType), intent(in) :: vel_grid
    double precision, intent(in) :: dens, u, v, w
    double precision :: temp

    integer :: i, j, k
    double precision :: x, y, z
    double precision :: Cxsq, Cysq, Czsq, Csq

    integer :: i_min,i_max,j_min,j_max,k_min,k_max
    double precision :: mass

    mass = molecule%mass

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    temp = zero

    do k = k_min, k_max
       z = vel_grid%z(k)
       Czsq = (z - w)
       Czsq = Czsq*Czsq
       
       do j = j_min, j_max
          y = vel_grid%y(j)
          Cysq = (y - v)
          Cysq = Cysq*Cysq

          do i = i_min, i_max
             x = vel_grid%x(i)
             Cxsq = (x - u)
             Cxsq = Cxsq*Cxsq

             Csq = Cxsq + Cysq + Czsq

             temp = temp + phi%value(i,j,k)*Csq
             
          enddo
       enddo
    enddo

    temp = 2.d0*mass*temp/(3.d0*dens)

    return
  end subroutine compute_temperature

  subroutine compute_rotational_energy( rot_energy, phi, molecule, vel_grid )

    implicit none

    double precision :: rot_energy

    type(DistFuncType), intent(in) :: phi
    type(MoleculeType), intent(in) :: molecule
    type(VelocityGridType), intent(in) :: vel_grid

    integer :: l, i, j, k

    integer :: r_modes, r_levels
    integer :: i_min, i_max, j_min, j_max, k_min, k_max

    double precision :: mass, fraction, energy_level

    r_modes = molecule%rot_modes
    rot_energy = zero

    if( r_modes .eq. 0 )return

    r_levels = phi%num_rot_levels

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    mass = molecule%mass

    do k = k_min,k_max
       do j = j_min,j_max
          do i = i_min,i_max

             do l = 1, r_modes*r_levels

                energy_level = phi%rot_level(l)
                fraction = phi%rot(l,i,j,k)
                rot_energy = rot_energy + fraction*energy_level

             end do

          end do
       end do
    end do

    rot_energy = rot_energy*mass

    return
  end subroutine compute_rotational_energy

  subroutine compute_rotational_temperature( T_rot, rot_energy, rot_dof, molecule )

    implicit none

    double precision :: T_rot

    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: rot_energy, rot_dof

    T_rot = zero

    if( molecule%rot_modes .eq. 0 )return

    T_rot = (4.0d0/rot_dof)*rot_energy ! J/m^3*m^3 - J (boltz const scales to Kelvin)   

    return
  end subroutine compute_rotational_temperature

subroutine compute_vibrational_energy( vib_energy, phi, molecule, vel_grid )
  
    implicit none

    double precision :: vib_energy

    type(DistFuncType), intent(in) :: phi
    type(MoleculeType), intent(in) :: molecule
    type(VelocityGridType), intent(in) :: vel_grid

    integer :: v_modes, v_levels
    integer :: l, i, j, k
    integer :: i_min, i_max, j_min, j_max, k_min, k_max

    double precision :: mass, fraction, energy_level

    v_modes = molecule%vib_modes
    v_levels = phi%num_vib_levels
    vib_energy = zero

    if( v_modes .eq. 0 )return

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    mass = molecule%mass

    do k = k_min,k_max
       do j = j_min,j_max
          do i = i_min,i_max

             do l = 1, v_modes*v_levels

                energy_level = phi%vib_level(l)
                fraction = phi%vib(l,i,j,k)
                vib_energy = vib_energy + fraction*energy_level

             end do

          end do
       end do
    end do

    vib_energy = vib_energy*mass
    
    return
  end subroutine compute_vibrational_energy

  subroutine compute_vibrational_temperature2( T_vib, dens, phi, molecule, vel_grid )

    implicit none

    type(DistFuncType), intent(in) :: phi
    type(MoleculeType), intent(in) :: molecule
    type(VelocityGridType), intent(in) :: vel_grid
    double precision, intent(in) :: dens

    double precision :: T_vib

    double precision :: dens0, dens1, energy_level

    integer :: v_modes
    integer :: i, j, k
    integer :: i_min, i_max, j_min, j_max, k_min, k_max

    v_modes = molecule%vib_modes
    T_vib   = zero

    if( v_modes .eq. 0 )return

    dens0 = zero
    dens1 = zero

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    do k = k_min,k_max
       do j = j_min,j_max
          do i = i_min,i_max

             dens0 = dens0 + phi%vib(1,i,j,k)
             dens1 = dens1 + phi%vib(2,i,j,k)

          end do
       end do
    end do

    energy_level = phi%vib_level(2)

    T_vib = two*energy_level/( log( dens0/dens1 ) )
    
    return
  end subroutine compute_vibrational_temperature2

  subroutine compute_vib_dof2( vib_dof, T_vib, vib_energy, molecule )
    
    implicit none

    double precision :: vib_dof

    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: T_vib
    double precision, intent(in) :: vib_energy

    ! TODO: The current vibrational dof formula is particular to the simple harmonic oscillator

    vib_dof = zero

    if( molecule%vib_modes .eq. 0 )return

    vib_dof = four*vib_energy/T_vib

    return
  end subroutine compute_vib_dof2

  subroutine compute_vibrational_temperature( T_vib, dens, vib_energy, molecule )

    implicit none

    double precision :: T_vib

    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: dens, vib_energy
    
    double precision :: theta_v

    ! TODO: The current vibrational temperature formula is particular to the simple harmonic oscillator

    T_vib = zero

    if( molecule%vib_modes .eq. 0 )return

    theta_v = molecule%theta_v

    T_vib = theta_v/log( one_half*theta_v*dens/vib_energy + one )

    return
  end subroutine compute_vibrational_temperature

  subroutine compute_vib_dof( vib_dof, T_vib, molecule )
    
    implicit none

    double precision :: vib_dof

    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: T_vib

    double precision :: theta_v

    ! TODO: The current vibrational dof formula is particular to the simple harmonic oscillator

    vib_dof = zero

    if( molecule%vib_modes .eq. 0 )return

    theta_v = molecule%theta_v

    vib_dof = ( two*theta_v/T_vib )/( exp( theta_v/T_vib ) - one )

    return
  end subroutine compute_vib_dof

  subroutine compute_overall_temperature( temp_ov, trans, rot, vib, rot_dof, vib_dof )

    implicit none

    double precision :: temp_ov

    double precision, intent(in) :: trans, rot, vib
    double precision, intent(in) :: vib_dof, rot_dof

    double precision :: numerator, denominator

!!$    numerator = 3.0d0*trans
!!$    denominator = 3.0d0
!!$
!!$    numerator = numerator + rot_dof*rot
!!$    denominator = denominator + rot_dof
!!$
!!$    numerator = numerator + vib_dof*vib
!!$    denominator = denominator + vib_dof
!!$
!!$    temp_ov = numerator/denominator

    numerator = 4.0d0*( trans + vib + rot )
    denominator = 3.0d0 + rot_dof + vib_dof

    temp_ov = numerator/denominator

    return
  end subroutine compute_overall_temperature

  subroutine compute_avg_speed( avg_speed, dens, phi, vel_grid )

    implicit none

    type(DistFuncType), intent(in) :: phi
    type(VelocityGridType), intent(in) :: vel_grid
    double precision, intent(in) :: dens
    
    double precision :: avg_speed

    integer :: i, j, k
    double precision :: x, y, z
    double precision :: Cxsq, Cysq, Czsq, Csq

    integer :: i_min,i_max,j_min,j_max,k_min,k_max

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    avg_speed = zero

    do k = k_min, k_max
       z = vel_grid%z(k)
       Czsq = z*z
       
       do j = j_min, j_max
          y = vel_grid%y(j)
          Cysq = y*y

          do i = i_min, i_max
             x = vel_grid%x(i)
             Cxsq = x*x

             Csq = Cxsq + Cysq + Czsq

             avg_speed = avg_speed + phi%value(i,j,k)*sqrt(Csq)
             
          enddo
       enddo
    enddo

    avg_speed = avg_speed/dens
    
    return
  end subroutine compute_avg_speed

  subroutine compute_coll_freq( coll_freq, dens, temp, molecule )

    implicit none

    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: dens, temp

    double precision :: coll_freq

    double precision :: omega, coll_freq_cnst

    omega = molecule%omega
    coll_freq_cnst = molecule%coll_freq_cnst

    coll_freq = coll_freq_cnst*dens*temp**( one - omega )

    return
  end subroutine compute_coll_freq

  ! Mixture Properties
  subroutine compute_total_density( tot_dens, dens )

    implicit none

    double precision, dimension(:), intent(in) :: dens

    integer :: n
    double precision :: tot_dens, species_dens

    tot_dens = 0.0d0

    do n = 1, num_species

       species_dens = dens( n )
       tot_dens = tot_dens + species_dens

    end do

    return
  end subroutine compute_total_density

  subroutine compute_macro_density( macro_dens, dens, molecule )

    implicit none

    type(MoleculeType), dimension(:), intent(in) :: molecule
    double precision, dimension(:), intent(in) :: dens

    double precision :: macro_dens

    integer :: n
    double precision :: species_dens, species_mass

    macro_dens = 0.0d0

    do n = 1, num_species

       species_dens = dens( n )
       species_mass = molecule(n)%mass
       macro_dens = macro_dens + species_mass*species_dens

    end do

    return
  end subroutine compute_macro_density

  subroutine compute_total_x_vel( tot_xvel, dens, xvel, macro_dens, molecule )

    implicit none

    type(MoleculeType), dimension(:), intent(in) :: molecule
    double precision, dimension(:), intent(in) :: dens, xvel
    double precision, intent(in) :: macro_dens

    double precision :: tot_xvel
    double precision :: species_dens, species_xvel
    double precision :: mass

    integer :: n

    tot_xvel = zero

    do n = 1, num_species

       species_dens = dens( n )
       species_xvel = xvel( n )
       mass = molecule(n)%mass
       tot_xvel = tot_xvel + species_dens*mass*species_xvel

    end do

    tot_xvel = tot_xvel/macro_dens

    return
  end subroutine compute_total_x_vel

  subroutine compute_total_y_vel( tot_yvel, dens, yvel, macro_dens, molecule )

    implicit none

    type(MoleculeType), dimension(:), intent(in) :: molecule
    double precision, dimension(:), intent(in) :: dens, yvel
    double precision, intent(in) :: macro_dens

    double precision :: tot_yvel
    double precision :: species_dens, species_yvel
    double precision :: mass

    integer :: n

    tot_yvel = 0.0d0

    do n = 1, num_species

       species_dens = dens( n )
       species_yvel = yvel( n )
       mass = molecule(n)%mass
       tot_yvel = tot_yvel + species_dens*mass*species_yvel

    end do

    tot_yvel = tot_yvel/macro_dens

    return
  end subroutine compute_total_y_vel

  subroutine compute_total_z_vel( tot_zvel, dens, zvel, macro_dens, molecule )

    implicit none

    type(MoleculeType), dimension(:), intent(in) :: molecule
    double precision, dimension(:), intent(in) :: dens, zvel
    double precision, intent(in) :: macro_dens

    double precision :: tot_zvel
    double precision :: species_dens, species_zvel
    double precision :: mass

    integer :: n

    tot_zvel = zero

    do n = 1, num_species

       species_dens = dens( n )
       species_zvel = zvel( n )
       mass = molecule(n)%mass
       tot_zvel = tot_zvel + species_dens*mass*species_zvel

    end do

    tot_zvel = tot_zvel/macro_dens

    return
  end subroutine compute_total_z_vel

  subroutine compute_total_temperature( tot_temp, tot_dens, dens, temp )

    implicit none

    double precision, dimension(:), intent(in) :: dens, temp
    double precision, intent(in) :: tot_dens
    double precision :: tot_temp, species_temp
    double precision :: species_dens

    integer :: n

    tot_temp = zero

    do n = 1, num_species

       species_dens = dens( n )
       species_temp = temp( n )*species_dens
       tot_temp = tot_temp + species_temp

    end do

    tot_temp = tot_temp/tot_dens

    return
  end subroutine compute_total_temperature

  subroutine compute_total_avg_speed( tot_avg_speed, dens, avg_speed, macro_dens, molecule )

    implicit none

    type(MoleculeType), dimension(:), intent(in) :: molecule
    double precision, dimension(:), intent(in) :: dens, avg_speed
    double precision, intent(in) :: macro_dens

    double precision :: tot_avg_speed
    double precision :: mass

    integer :: n

    tot_avg_speed = zero
    
    do n = 1, num_species

       mass = molecule(n)%mass
       tot_avg_speed = tot_avg_speed + mass*dens(n)*avg_speed(n)

    end do

    tot_avg_speed = tot_avg_speed/macro_dens

    return
  end subroutine compute_total_avg_speed
       

  subroutine compute_total_dof()

    return
  end subroutine compute_total_dof

  subroutine compute_total_energy( mix_energy, energy )

    implicit none

    double precision :: mix_energy

    double precision, dimension(:), intent(in) :: energy

    integer :: n 

    mix_energy = zero

    do n = 1, num_species

       mix_energy = mix_energy + energy(n)

    end do

    return
  end subroutine compute_total_energy

  subroutine compute_total_internal_temperature( temp, species_dens, species_temp, dof )

    implicit none

    double precision, dimension(:), intent(in) :: species_dens, species_temp, dof

    double precision :: temp

    double precision :: dens

    integer :: n 

    temp = zero
    dens = zero

    do n = 1, num_species

       temp = temp + species_dens(n)*dof(n)*species_temp(n)
       dens = dens + species_dens(n)*dof(n)

    end do

    if( dens .lt. double_tol )then
       temp = zero
    else
       temp = temp/dens
    end if

    return
  end subroutine compute_total_internal_temperature

  ! Properties NOT part of the properties structure - Only needed on output
  subroutine compute_speed( speed, u, v, w )

    implicit none

    double precision, intent(in) :: u, v, w
    
    double precision :: speed

    speed = sqrt( u*u + v*v + w*w )

    return
  end subroutine compute_speed

  subroutine compute_mach_number( mach, speed, temp )

    implicit none

    double precision, intent(in) :: speed, temp

    double precision :: mach

    double precision :: speed_of_sound, gamma

    gamma = zero

    speed_of_sound = sqrt( one_half*gamma*temp )

    mach = speed/speed_of_sound

    return
  end subroutine compute_mach_number

  subroutine compute_directional_temperature( temp, dens, vel, direction, phi, molecule, vel_grid )

    implicit none

    type(DistFuncType), intent(in) :: phi
    type(MoleculeType), intent(in) :: molecule
    type(VelocityGridType), intent(in) :: vel_grid
    double precision, intent(in) :: dens, vel
    integer, intent(in) :: direction

    double precision :: temp

    integer :: i, j, k
    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    double precision :: C, Csq, mass

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    mass = molecule%mass

    temp = zero

    do k = k_min, k_max
       do j = j_min, j_max
          do i = i_min, i_max

             select case( direction )
             case( x_dir )
                C = vel_grid%x(i)

             case( y_dir )
                C = vel_grid%y(j)

             case( z_dir )
                C = vel_grid%z(k)

             case default
                write(*,*) "Error: Invalid value of direction: ", direction
                stop

             end select

             Csq = ( C - vel )*( C - vel )

             temp = temp + phi%value(i,j,k)*Csq

          end do
       end do
    end do

    temp = two_thirds*mass*temp/dens

    return
  end subroutine compute_directional_temperature

  subroutine compute_entropy( entropy, dens, phi, vel_grid )

    implicit none

    type(DistFuncType), intent(in) :: phi
    type(VelocityGridType), intent(in) :: vel_grid
    double precision, intent(in) :: dens

    double precision :: entropy
    
    integer :: i, j, k
    integer :: i_min, i_max, j_min, j_max, k_min, k_max

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    entropy = zero

    do k = k_min,k_max
       do j = j_min,j_max
          do i = i_min,i_max

             if( phi%value(i,j,k) .gt. zero )then
                entropy = entropy - phi%value(i,j,k)*log( phi%value(i,j,k) )
             endif

          enddo
       enddo
    enddo

    entropy = entropy/dens

    return
  end subroutine compute_entropy

  subroutine compute_heat_flux( q, u, v, w, direction, phi, molecule, vel_grid )

    implicit none

    type(DistFuncType), intent(in) :: phi
    type(MoleculeType), intent(in) :: molecule
    type(VelocityGridType), intent(in) :: vel_grid
    double precision, intent(in) :: u, v, w
    integer, intent(in) :: direction

    double precision :: q

    integer :: i, j, k
    integer :: i_min, i_max, j_min, j_max, k_min, k_max

    double precision :: x, y, z, Cxsq, Cysq, Czsq, Csq
    double precision :: mass

    mass = molecule%mass

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    q = zero

    do k = k_min, k_max
       z = vel_grid%z(k)
       Czsq = ( z - w )*( z - w )

       do j = j_min, j_max
          y = vel_grid%y(j)
          Cysq = ( y - v )*( y - v )

          do i = i_min, i_max
             x = vel_grid%x(i)
             Cxsq = ( x - u )*( x - u )

             Csq = Cxsq + Cysq + Czsq

             select case( direction )
             case( x_dir )
                q = q + ( x - u )*Csq*phi%value(i,j,k)

             case( y_dir )
                q = q + ( y - v )*Csq*phi%value(i,j,k)

             case( z_dir )
                q = q + ( z - w )*Csq*phi%value(i,j,k)

             case default
                write(*,*) "Error: Invalid value of direction: ", direction
                stop

             end select

          end do
       end do
    end do

    q = q*mass

    return
  end subroutine compute_heat_flux

  subroutine compute_shear_stress( tau, dens, u, v, w, temp, direction1, direction2, phi, vel_grid )

    implicit none

    type(DistFuncType), intent(in) :: phi
    type(VelocityGridType), intent(in) :: vel_grid
    double precision, intent(in) :: dens, u, v, w, temp
    integer, intent(in) :: direction1, direction2

    double precision :: tau

    integer :: i, j, k
    integer :: i_min, i_max, j_min, j_max, k_min, k_max

    double precision :: x, y, z
    double precision :: vel1, vel2

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    tau = zero

    do k = k_min, k_max
       z = vel_grid%z(k)

       do j = j_min, j_max
          y = vel_grid%y(j)

          do i = i_min, i_max
             x = vel_grid%x(i)

             select case( direction1 )
             case( x_dir )
                vel1 = x - u
             case( y_dir )
                vel1 = y - v
             case( z_dir )
                vel1 = z - w
             case default
                write(*,*) "Error: Invalid direction 1 value: ", direction1
                stop
             end select

             select case( direction2 )
             case( x_dir )
                vel2 = x - u
             case( y_dir )
                vel2 = y - v
             case( z_dir )
                vel2 = z - w
             case default
                write(*,*) "Error: Invalid direction 2 value: ", direction2
                stop
             end select

             tau = tau - vel1*vel2*phi%value(i,j,k)

          end do
       end do
    end do

    tau = two*tau

    if( direction1 .eq. direction2 )then
       tau = tau + one_half*dens*temp
    end if      

    return
  end subroutine compute_shear_stress

  subroutine compute_moment( moment, moment_flag, dens, u, v, w, temp, phi, molecule, vel_grid )

    implicit none

    type(DistFuncType), intent(in) :: phi
    type(MoleculeType), intent(in) :: molecule
    type(VelocityGridType), intent(in) :: vel_grid
    double precision, intent(in) :: dens, u, v, w, temp
    integer, intent(in) :: moment_flag

    double precision :: moment

    double precision :: mass
    double precision :: d_moment_flag
    double precision :: nfact, mom_MB, moment_change
    double precision :: x_vel, y_vel, z_vel, Cxsq, Cysq, Czsq, C

    integer :: i, j, k, in
    integer :: i_min, i_max, j_min, j_max, k_min, k_max

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    mass = molecule%mass

    d_moment_flag = dble( moment_flag )

    nfact = one
    in = 1
    do
       moment_change = dble( moment_flag + in )

       if( moment_flag + in .ge. 1 )then
          nfact = nfact*moment_change
          in = in - 2

       else
          exit

       end if

    end do

    mom_MB = nfact*temp**( 1.5d0 + d_moment_flag )/( two**( one_half*d_moment_flag ) )

    moment = zero

    do k = k_min, k_max
       z_vel = vel_grid%z(k)
       Czsq  = ( z_vel - w )
       Czsq  = Czsq*Czsq

       do j = j_min, j_max
          y_vel = vel_grid%y(j)
          Cysq  = ( y_vel - v )
          Cysq  = Cysq*Cysq

          do i = i_min, i_max
             x_vel = vel_grid%x(i)
             Cxsq  = ( x_vel - u )
             Cxsq  = Cxsq*Cxsq

             C = sqrt( Cxsq + Cysq + Czsq )

             moment = moment + C**( d_moment_flag )*phi%value(i,j,k)

          end do
       end do
    end do

    moment = mass**( one_half*d_moment_flag )*moment/( dens*mom_MB )

    return
  end subroutine compute_moment

end module PhysicalProperties
