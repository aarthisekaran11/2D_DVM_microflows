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
module BoundaryConditions

  use ErrorCheck
  use Constants

  implicit none

  private

  integer :: boundary_condition_LW, boundary_condition_RW

  double precision, allocatable, dimension(:) :: density_LW, u_LW, temp_LW
  double precision, allocatable, dimension(:) :: density_RW, u_RW, temp_RW
  double precision, allocatable, dimension(:) :: temp_rot_LW, temp_rot_RW
  double precision, allocatable, dimension(:) :: temp_vib_LW, temp_vib_RW

  integer, parameter :: specular_reflection = 1
  integer, parameter :: diffuse_reflection  = 2
  integer, parameter :: fixed_in_out_flow   = 3
  integer, parameter :: zero_gradient_flow  = 4

  public :: initialize_property_arrays
  public :: set_boundary_conditions_flags
  public :: get_boundary_conditions_flags
  public :: set_boundary_properties
  public :: set_temp_rot
  public :: set_temp_vib
  public :: get_right_wall_temp
  public :: get_right_wall_density
  public :: get_right_wall_velocity
  public :: get_right_wall_temp_rot
  public :: get_right_wall_temp_vib
  public :: get_left_wall_temp
  public :: get_left_wall_density
  public :: get_left_wall_velocity
  public :: get_left_wall_temp_rot
  public :: get_left_wall_temp_vib
  public :: apply_boundary_conditions

contains

  subroutine set_boundary_conditions_flags( left_wall_bc_in, right_wall_bc_in )

    implicit none

    integer, intent(in) :: left_wall_bc_in, right_wall_bc_in

    boundary_condition_LW = left_wall_bc_in
    boundary_condition_RW = right_wall_bc_in

    return
  end subroutine set_boundary_conditions_flags

  subroutine get_boundary_conditions_flags( left_wall_bc_out, right_wall_bc_out )

    implicit none

    integer :: left_wall_bc_out, right_wall_bc_out

    left_wall_bc_out  = boundary_condition_LW
    right_wall_bc_out = boundary_condition_RW

    return
  end subroutine get_boundary_conditions_flags

  subroutine initialize_property_arrays( )

    use DistFunc, only : num_species

    implicit none

    integer :: status

    allocate( density_LW( num_species ), STAT=status )
    call allocate_error_check( status, "density_LW" )
    allocate( u_LW( num_species ), STAT=status )
    call allocate_error_check( status, "u_LW" )
    allocate( temp_LW( num_species ), STAT=status )
    call allocate_error_check( status, "temp_LW" )
    allocate( temp_rot_LW( num_species ), STAT=status )
    call allocate_error_check( status, "temp_rot_LW" )
    allocate( temp_vib_LW( num_species ), STAT=status )
    call allocate_error_check( status, "temp_vib_LW" )

    allocate( density_RW( num_species ), STAT=status )
    call allocate_error_check( status, "density_RW" )
    allocate( u_RW( num_species ), STAT=status )
    call allocate_error_check( status, "u_RW" )
    allocate( temp_RW( num_species ), STAT=status )
    call allocate_error_check( status, "temp_RW" )
    allocate( temp_rot_RW( num_species ), STAT=status )
    call allocate_error_check( status, "temp_rot_RW" )
    allocate( temp_vib_RW( num_species ), STAT=status )
    call allocate_error_check( status, "temp_vib_RW" )

    return
  end subroutine initialize_property_arrays
  
  subroutine set_boundary_properties( density_LW_in, u_LW_in, temp_LW_in, &
       density_RW_in, u_RW_in, temp_RW_in, species )

    implicit none

    double precision, intent(in) :: density_LW_in, u_LW_in, temp_LW_in
    double precision, intent(in) :: density_RW_in, u_RW_in, temp_RW_in
    integer, intent(in) :: species

    density_LW( species ) = density_LW_in
    u_LW( species ) = u_LW_in
    temp_LW( species ) = temp_LW_in

    density_RW( species ) = density_RW_in  
    u_RW( species ) = u_RW_in
    temp_RW( species ) = temp_RW_in

    return
  end subroutine set_boundary_properties

  subroutine set_temp_rot( temp_rot_LW_in, temp_rot_RW_in, species )

    implicit none

    double precision, intent(in) :: temp_rot_LW_in, temp_rot_RW_in
    integer, intent(in) :: species

    temp_rot_LW( species ) = temp_rot_LW_in
    temp_rot_RW( species ) = temp_rot_RW_in

    return
  end subroutine set_temp_rot

  subroutine set_temp_vib( temp_vib_LW_in, temp_vib_RW_in, species )

    implicit none

    double precision, intent(in) :: temp_vib_LW_in, temp_vib_RW_in
    integer, intent(in) :: species

    temp_vib_LW( species ) = temp_vib_LW_in
    temp_vib_RW( species ) = temp_vib_RW_in

    return
  end subroutine set_temp_vib

  subroutine get_right_wall_temp( temp_RW_set, species )

    implicit none

    double precision :: temp_RW_set
    integer, intent(in) :: species

    temp_RW_set = temp_RW( species )
    return
  end subroutine get_right_wall_temp

  subroutine get_left_wall_temp( temp_LW_set, species )

    implicit none

    double precision :: temp_LW_set
    integer, intent(in) :: species

    temp_LW_set = temp_LW( species )
    return
  end subroutine get_left_wall_temp

  subroutine get_right_wall_density( density_RW_set, species )

    implicit none

    double precision :: density_RW_set
    integer, intent(in) :: species

    density_RW_set = density_RW( species )
    return
  end subroutine get_right_wall_density

  subroutine get_left_wall_density( density_LW_set, species )

    implicit none

    double precision :: density_LW_set
    integer, intent(in) :: species

    density_LW_set = density_LW( species )
    return
  end subroutine get_left_wall_density

  subroutine get_right_wall_velocity( u_RW_set, species )

    implicit none

    double precision :: u_RW_set
    integer, intent(in) :: species

    u_RW_set = u_RW( species )
    return
  end subroutine get_right_wall_velocity

  subroutine get_left_wall_velocity( u_LW_set, species )

    implicit none

    double precision :: u_LW_set
    integer, intent(in) :: species

    u_LW_set = u_LW( species )
    return
  end subroutine get_left_wall_velocity

  subroutine get_right_wall_temp_rot( temp_rot_RW_set, species )

    implicit none

    double precision :: temp_rot_RW_set
    integer, intent(in) :: species

    temp_rot_RW_set = temp_rot_RW( species )
    return
  end subroutine get_right_wall_temp_rot

  subroutine get_left_wall_temp_rot( temp_rot_LW_set, species )

    implicit none

    double precision :: temp_rot_LW_set
    integer, intent(in) :: species

    temp_rot_LW_set = temp_rot_LW( species )
    return
  end subroutine get_left_wall_temp_rot

  subroutine get_right_wall_temp_vib( temp_vib_RW_set, species )

    implicit none

    double precision :: temp_vib_RW_set
    integer, intent(in) :: species

    temp_vib_RW_set = temp_vib_RW( species )

    return
  end subroutine get_right_wall_temp_vib

  subroutine get_left_wall_temp_vib( temp_vib_LW_set, species )

    implicit none

    double precision :: temp_vib_LW_set
    integer, intent(in) :: species

    temp_vib_LW_set = temp_vib_LW( species )

    return
  end subroutine get_left_wall_temp_vib

  subroutine apply_boundary_conditions( phi, phi_conv, molecule, vel_grid, cfl, species )

    use DistFunc
    use VelocityGrid
    use PhysicalGrid
    use SpeciesAndReferenceData

    implicit none

    type(DistFuncType), dimension(:), intent(in) :: phi
    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    type(MoleculeType), intent(in) :: molecule

    type(DistFuncType), dimension(:) :: phi_conv

    double precision, intent(in) :: cfl
    integer, intent(in) :: species

    integer :: nspace
    integer :: ns, index_1, index_2, gr_1, gr_2

    call get_nspace( nspace )
    
    select case( boundary_condition_LW )
    case( specular_reflection )
       ns = 1

       call get_spatial_reference( ns, gr_1 )
       call get_spatial_reference( ns+1, gr_2 )
       if( gr_1 .ne. gr_2 )then
          write(*,*) "Grid types at left boundary elements must be the same."
          stop
       end if

       call apply_specular_wall( phi, phi_conv, molecule, vel_grid(gr_1), cfl, ns )

    case( diffuse_reflection )

    case( fixed_in_out_flow )
       ns = 1

       call get_spatial_reference( ns, gr_1 )

       call apply_fixed_in_out_flow( phi_conv, ns, species, molecule, vel_grid(gr_1) )

    case( zero_gradient_flow )
       index_1 = 1
       index_2 = 2

       call get_spatial_reference( index_1, gr_1 )
       call get_spatial_reference( index_2, gr_2 )
       if( gr_1 .ne. gr_2 )then
          write(*,*)"Grid types at left boundary elements must be the same."
          stop
       end if

       call apply_zero_gradient_flow( phi_conv, molecule, vel_grid(gr_1), index_1, index_2 )

    case default
       write(*,*) "Error: Invalid value of boundary_condition_LW - :", &
            boundary_condition_LW
    end select

    select case( boundary_condition_RW )
    case( specular_reflection )
       ns = nspace

       call get_spatial_reference( nspace, gr_1 )
       call get_spatial_reference( nspace-1, gr_2 )
       if( gr_1 .ne. gr_2 )then
          write(*,*)"Grid types at right boundary elements must be the same."
          stop
       end if

       call apply_specular_wall( phi, phi_conv, molecule, vel_grid(gr_1), cfl, ns )

    case( diffuse_reflection )

    case( fixed_in_out_flow )
       ns = nspace

       call get_spatial_reference( ns, gr_1 )

       call apply_fixed_in_out_flow( phi_conv, ns, species, molecule, vel_grid(gr_1) )

    case( zero_gradient_flow )
       index_1 = nspace
       index_2 = nspace-1

       call get_spatial_reference( index_1, gr_1 )
       call get_spatial_reference( index_2, gr_2 )
       if( gr_1 .ne. gr_2 )then
          write(*,*)"grid types at right boundary elements must be the same"
          stop
       end if

       call apply_zero_gradient_flow( phi_conv, molecule, vel_grid(gr_1), index_1, index_2 )

    case default
       write(*,*) "Error: Invalid value of boundary_condition_RW - :", &
            boundary_condition_RW
    end select

    return
  end subroutine apply_boundary_conditions

  subroutine apply_specular_wall( phi, phi_conv, molecule, vel_grid, cfl, ns )

    use DistFunc
    use VelocityGrid
    use PhysicalGrid
    use SpeciesAndReferenceData

    implicit none

    type(DistFuncType), dimension(:), intent(in) :: phi
    type(VelocityGridType), intent(in) :: vel_grid
    type(MoleculeType), intent(in) :: molecule

    type(DistFuncType), dimension(:) :: phi_conv

    double precision, intent(in) :: cfl

    integer, intent(in) :: ns

    double precision :: factor, x_vel, opp_vel

    integer :: r_modes, v_modes

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: r_levels, v_levels
    integer :: nspace
    integer :: i, j, k, l

    call get_nspace( nspace )

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    r_modes  = molecule%rot_modes
    v_modes  = molecule%vib_modes
    r_levels = phi(ns)%num_rot_levels
    v_levels = phi(ns)%num_vib_levels

    do k = k_min, k_max
       do j = j_min, j_max
          do i = i_min, i_max

             x_vel = vel_grid%x(i)
             factor = cfl*x_vel

             if( ns .eq. 1 )then
                if( x_vel .lt. zero )then

                   phi_conv(ns)%value(i,j,k) = phi(ns)%value(i,j,k) - &
                        factor*( phi(ns+1)%value(i,j,k) - phi(ns)%value(i,j,k) )

                   if( r_modes .gt. 0 )then
                      do l = 1, r_levels
                         phi_conv(ns)%rot(l,i,j,k) = phi(ns)%rot(l,i,j,k) - &
                              factor*( phi(ns+1)%rot(l,i,j,k) - phi(ns)%rot(l,i,j,k) )
                      end do
                   end if

                   if( r_modes .gt. 0 )then
                      do l = 1, v_levels
                         phi_conv(ns)%vib(l,i,j,k) = phi(ns)%vib(l,i,j,k) - &
                              factor*( phi(ns+1)%vib(l,i,j,k) - phi(ns)%vib(l,i,j,k) )
                      end do
                   end if
                   
                   opp_vel = -x_vel

                   ! TODO: flip each velocity and map back onto the grid.
                   ! TODO: uniform grid with x_vel = 0 at a node would map exactly as a reflection.

                else

                   ! TODO: reduce the current density how? phi_conv = phi_conv + phi - cfl*phi?

                end if
             end if

             if( ns .eq. nspace )then
                if( x_vel .gt. zero )then

                   phi_conv(ns)%value(i,j,k) = phi(ns)%value(i,j,k) - &
                        factor*( phi(ns)%value(i,j,k) - phi(ns-1)%value(i,j,k) )

                   if( r_modes .gt. 0 )then
                      do l = 1, r_levels
                         phi_conv(ns)%rot(l,i,j,k) = phi(ns)%rot(l,i,j,k) - &
                              factor*( phi(ns)%rot(l,i,j,k) - phi(ns-1)%rot(l,i,j,k) )
                      end do
                   end if

                   if( r_modes .gt. 0 )then
                      do l = 1, v_levels
                         phi_conv(ns)%vib(l,i,j,k) = phi(ns)%vib(l,i,j,k) - &
                              factor*( phi(ns)%vib(l,i,j,k) - phi(ns-1)%vib(l,i,j,k) )
                      end do
                   end if

                   opp_vel = -x_vel

                   ! TODO: flip each velocity and map back onto the grid.
                   ! TODO: uniform grid with x_vel = 0 at a node would map exactly as a reflection.

                else

                   ! TODO: reduce the current density how? phi_conv = phi_conv + phi - cfl*phi?

                end if
             end if

          end do
       end do
    end do
                
    return
  end subroutine apply_specular_wall

  subroutine apply_fixed_in_out_flow( phi_conv, ns, species, molecule, vel_grid )

    use DistFunc
    use VelocityGrid
    use PhysicalGrid
    use SpeciesAndReferenceData

    implicit none

    type(DistFuncType), dimension(:) :: phi_conv

    type(VelocityGridType), intent(in) :: vel_grid
    type(MoleculeType), intent(in) :: molecule

    integer, intent(in) :: ns, species

    double precision :: dens, u, v, w, temp, temp_rot, temp_vib, mass
    double precision :: maxwell_value, maxwell_coeff

    double precision :: energy_level, fraction
    double precision :: vib_part_func
    double precision :: theta_r, theta_v
    double precision, dimension(:), allocatable :: rot_levels, rot_df, vib_levels, vib_df

    double precision :: x, y, z
    double precision :: beta_x, beta_y, beta_z, beta3

    double precision :: even_spin, odd_spin

    integer :: r_modes, v_modes
    integer :: nspace
    integer :: mol_shape

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: r_levels, v_levels
    integer :: i, j, k, l
    
    call get_nspace( nspace )

    ! Get internal energy structure
    r_modes  = molecule%rot_modes
    v_modes  = molecule%vib_modes

    mol_shape = molecule%molecule_type
    even_spin = molecule%even_spin
    odd_spin  = molecule%odd_spin

    ! Get boundary conditions
    if( ns .eq. 1 )then
       call get_left_wall_density( dens, species )
       call get_left_wall_temp( temp, species )
       call get_left_wall_velocity( u, species )
       v = zero
       w = zero
       if( r_modes .gt. 0 ) call get_left_wall_temp_rot( temp_rot, species )
       if( v_modes .gt. 0 ) call get_left_wall_temp_vib( temp_vib, species )
       
    else if( ns .eq. nspace )then
       call get_right_wall_density( dens, species )
       call get_right_wall_temp( temp, species )
       call get_right_wall_velocity( u, species )
       v = zero
       w = zero
       if( r_modes .gt. 0 ) call get_right_wall_temp_rot( temp_rot, species )
       if( v_modes .gt. 0 ) call get_right_wall_temp_vib( temp_vib, species )

    else
       write(*,*) "Error. ns is not on a boundary. Given value is: ",ns
       stop

    end if

    if( r_modes .gt. 0 )then
       r_levels = phi_conv(ns)%num_rot_levels
       allocate( rot_levels(1:r_levels) )
       allocate( rot_df(1:r_levels) )
       call compute_rot_distribution( rot_df, rot_levels, molecule, temp_rot, r_levels, species )
    end if

    if( v_modes .gt. 0 )then
       v_levels = phi_conv(ns)%num_vib_levels
       allocate( vib_levels(1:v_levels) )
       allocate( vib_df(1:v_levels) )
       call compute_vib_distribution( vib_df, vib_levels, molecule, temp_vib, v_levels, species )
    end if
    
    theta_r = molecule%theta_r
    theta_v = molecule%theta_v

    mass = molecule%mass

    maxwell_coeff = dens*sqrt( mass*mass*mass )/sqrt( pi*pi*pi*temp*temp*temp )
    if( v_modes .gt. 0 ) vib_part_func = one/( one - exp( -theta_v/temp_vib ) )

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    do k = k_min, k_max
       z = vel_grid%z(k)
       beta_z = vel_grid%beta_z(k)

       do j = j_min,j_max
          y = vel_grid%y(j)
          beta_y = vel_grid%beta_y(j)

          do i = i_min,i_max
             x = vel_grid%x(i)
             beta_x = vel_grid%beta_x(i)
             beta3 = beta_x*beta_y*beta_z

             call compute_maxwellian( x, y, z, maxwell_coeff, mass, u, v, w, temp, maxwell_value )

             phi_conv(ns)%value(i,j,k) = maxwell_value*beta3

             if( r_modes .gt. 0 )then
                phi_conv(ns)%rot(:,i,j,k) = maxwell_value*beta3*rot_df
!!$                do l = 1, r_levels
!!$                   call compute_rigid_rotor( fraction, energy_level, l, r_levels, temp_rot, molecule )
!!$                   phi_conv(ns)%rot(l,i,j,k) = fraction*maxwell_value*beta3
!!$                end do
             end if

             if( v_modes .gt. 0 )then
                phi_conv(ns)%vib(:,i,j,k) = maxwell_value*beta3*vib_df
!!$                do l = 1, v_levels
!!$                   call compute_SHO( fraction, energy_level, l, v_levels, temp_vib, molecule )
!!$                   phi_conv(ns)%vib(l,i,j,k) = fraction*maxwell_value*beta3
!!$                end do
             end if

          end do
       end do
    end do

    return
  end subroutine apply_fixed_in_out_flow

  subroutine apply_zero_gradient_flow( phi_conv, molecule, vel_grid, index_1, index_2 )

    use DistFunc
    use VelocityGrid
    use SpeciesAndReferenceData

    implicit none

    type(DistFuncType), dimension(:) :: phi_conv

    type(VelocityGridType), intent(in) :: vel_grid
    type(MoleculeType), intent(in) :: molecule
    
    integer, intent(in) :: index_1, index_2

    integer :: r_modes, v_modes

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: r_levels, v_levels
    integer :: i, j, k, l
    
    ! Get internal energy structure
    r_modes  = molecule%rot_modes
    v_modes  = molecule%vib_modes

    r_levels = phi_conv(index_1)%num_rot_levels
    v_levels = phi_conv(index_1)%num_vib_levels

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    do k = k_min, k_max
       do j = j_min, j_max
          do i = i_min, i_max

             phi_conv(index_1)%value(i,j,k) = phi_conv(index_2)%value(i,j,k)

             if( r_modes .gt. 0 )then
                do l = 1, r_levels
                   phi_conv(index_1)%rot(l,i,j,k) = phi_conv(index_2)%rot(l,i,j,k)
                end do
             end if

             if( v_modes .gt. 0 )then
                do l = 1, v_levels
                   phi_conv(index_1)%vib(l,i,j,k) = phi_conv(index_2)%vib(l,i,j,k)
                end do
             end if

          end do
       end do
    end do

    return
  end subroutine apply_zero_gradient_flow

end module BoundaryConditions
