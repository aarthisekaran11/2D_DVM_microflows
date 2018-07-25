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

  double precision :: acc_coeff

  double precision, allocatable, dimension(:) :: density_LW, u_LW, v_LW, w_LW, temp_LW
  double precision, allocatable, dimension(:) :: density_RW, u_RW, v_RW, w_RW, temp_RW
  double precision, allocatable, dimension(:) :: temp_rot_LW, temp_rot_RW
  double precision, allocatable, dimension(:) :: temp_vib_LW, temp_vib_RW

  integer, parameter :: specular_reflection = 1
  integer, parameter :: diffuse_reflection  = 2
  integer, parameter :: fixed_in_out_flow   = 3
  integer, parameter :: zero_gradient_flow  = 4

  integer, parameter :: x_only  = 1
  integer, parameter :: y_only  = 2
  integer, parameter :: z_only  = 4
  integer, parameter :: xy_only = 3
  integer, parameter :: xz_only = 5
  integer, parameter :: yz_only = 6
  integer, parameter :: xyz     = 7

  public :: initialize_property_arrays
  public :: set_boundary_conditions_flags
  public :: get_boundary_conditions_flags
  public :: set_accomodation_coefficient
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

  subroutine set_accomodation_coefficient( acc_coeff_in )
    
    implicit none

    double precision, intent(in) :: acc_coeff_in

    acc_coeff = acc_coeff_in

    return
  end subroutine set_accomodation_coefficient

  subroutine initialize_property_arrays( )

    use DistFunc, only : num_species

    implicit none

    integer :: status

    allocate( density_LW( num_species ), STAT=status )
    call allocate_error_check( status, "density_LW" )
    allocate( u_LW( num_species ), STAT=status )
    call allocate_error_check( status, "u_LW" )
    allocate( v_LW( num_species ), STAT=status )
    call allocate_error_check( status, "v_LW" )
    allocate( w_LW( num_species ), STAT=status )
    call allocate_error_check( status, "w_LW" )
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
    allocate( v_RW( num_species ), STAT=status )
    call allocate_error_check( status, "v_RW" )
    allocate( w_RW( num_species ), STAT=status )
    call allocate_error_check( status, "w_RW" )
    allocate( temp_RW( num_species ), STAT=status )
    call allocate_error_check( status, "temp_RW" )
    allocate( temp_rot_RW( num_species ), STAT=status )
    call allocate_error_check( status, "temp_rot_RW" )
    allocate( temp_vib_RW( num_species ), STAT=status )
    call allocate_error_check( status, "temp_vib_RW" )

    return
  end subroutine initialize_property_arrays
  
  subroutine set_boundary_properties( density_LW_in, u_LW_in, v_LW_in, w_LW_in, temp_LW_in, &
       density_RW_in, u_RW_in, v_RW_in, w_RW_in, temp_RW_in, species )

    implicit none

    double precision, intent(in) :: density_LW_in, u_LW_in, v_LW_in, w_LW_in, temp_LW_in
    double precision, intent(in) :: density_RW_in, u_RW_in, v_RW_in, w_RW_in, temp_RW_in
    integer, intent(in) :: species

    density_LW( species ) = density_LW_in
    u_LW( species )       = u_LW_in
    v_LW( species )       = v_LW_in
    w_LW( species )       = w_LW_in
    temp_LW( species )    = temp_LW_in

    density_RW( species ) = density_RW_in  
    u_RW( species )       = u_RW_in
    v_RW( species )       = v_RW_in
    w_RW( species )       = w_RW_in
    temp_RW( species )    = temp_RW_in

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

  subroutine get_right_wall_velocity( u_RW_set, v_RW_set, w_RW_set, species )

    implicit none

    double precision :: u_RW_set, v_RW_set, w_RW_set
    integer, intent(in) :: species

    u_RW_set = u_RW( species )
    v_RW_set = v_RW( species )
    w_RW_set = w_RW( species )
    return
  end subroutine get_right_wall_velocity

  subroutine get_left_wall_velocity( u_LW_set, v_LW_set, w_LW_set, species )

    implicit none

    double precision :: u_LW_set, v_LW_set, w_LW_set
    integer, intent(in) :: species

    u_LW_set = u_LW( species )
    v_LW_set = v_LW( species )
    w_LW_set = w_LW( species )
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

    integer :: nx_space, ny_space
    integer :: nx, index_1, index_2, gr_1, gr_2

    call get_nspace( nx_space, ny_space )
    
    select case( boundary_condition_LW )
    case( specular_reflection )
       nx = 1

       call get_spatial_reference( nx, gr_1 )
       call get_spatial_reference( nx+1, gr_2 )
       if( gr_1 .ne. gr_2 )then
          write(*,*) "Grid types at left boundary elements must be the same."
          stop
       end if

       call apply_specular_wall( phi, phi_conv, molecule, vel_grid(gr_1), cfl, nx )

    case( diffuse_reflection )
       nx = 1

       call get_spatial_reference( nx, gr_1 )
       call get_spatial_reference( nx+1, gr_2 )
       if( gr_1 .ne. gr_2 )then
          write(*,*) "Grid types at left boundary elements must be the same."
          stop
       end if

       call apply_diffuse_reflection( phi, phi_conv, species, molecule, vel_grid(gr_1), cfl, nx )


    case( fixed_in_out_flow )
       nx = 1

       call get_spatial_reference( nx, gr_1 )

       call apply_fixed_in_out_flow( phi_conv, nx, species, molecule, vel_grid(gr_1) )

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
       nx = nx_space

       call get_spatial_reference( nx_space, gr_1 )
       call get_spatial_reference( nx_space-1, gr_2 )
       if( gr_1 .ne. gr_2 )then
          write(*,*)"Grid types at right boundary elements must be the same."
          stop
       end if

       call apply_specular_wall( phi, phi_conv, molecule, vel_grid(gr_1), cfl, nx )

    case( diffuse_reflection )
       nx = nx_space

       call get_spatial_reference( nx_space, gr_1 )
       call get_spatial_reference( nx_space-1, gr_2 )
       if( gr_1 .ne. gr_2 )then
          write(*,*)"Grid types at right boundary elements must be the same."
          stop
       end if

       call apply_diffuse_reflection( phi, phi_conv, species, molecule, vel_grid(gr_1), cfl, nx )

    case( fixed_in_out_flow )
       nx = nx_space

       call get_spatial_reference( nx, gr_1 )

       call apply_fixed_in_out_flow( phi_conv, nx, species, molecule, vel_grid(gr_1) )

    case( zero_gradient_flow )
       index_1 = nx_space
       index_2 = nx_space-1

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

  subroutine apply_specular_wall( phi, phi_conv, molecule, vel_grid, cfl, nx )

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

    integer, intent(in) :: nx

    double precision :: factor, x_vel, opp_vel

    integer :: r_modes, v_modes

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: r_levels, v_levels
    integer :: nx_space, ny_space
    integer :: i, j, k, l

    call get_nspace( nx_space, ny_space )

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    r_modes  = molecule%rot_modes
    v_modes  = molecule%vib_modes
    r_levels = phi(nx)%num_rot_levels
    v_levels = phi(nx)%num_vib_levels

    do k = k_min, k_max
       do j = j_min, j_max
          do i = i_min, i_max

             x_vel = vel_grid%x(i)
             factor = cfl*x_vel

             if( nx .eq. 1 )then
                if( x_vel .lt. zero )then

                   phi_conv(nx)%value(i,j,k) = phi(nx)%value(i,j,k) - &
                        factor*( phi(nx+1)%value(i,j,k) - phi(nx)%value(i,j,k) )

                   if( r_modes .gt. 0 )then
                      do l = 1, r_levels
                         phi_conv(nx)%rot(l,i,j,k) = phi(nx)%rot(l,i,j,k) - &
                              factor*( phi(nx+1)%rot(l,i,j,k) - phi(nx)%rot(l,i,j,k) )
                      end do
                   end if

                   if( r_modes .gt. 0 )then
                      do l = 1, v_levels
                         phi_conv(nx)%vib(l,i,j,k) = phi(nx)%vib(l,i,j,k) - &
                              factor*( phi(nx+1)%vib(l,i,j,k) - phi(nx)%vib(l,i,j,k) )
                      end do
                   end if
                   
                   opp_vel = -x_vel

                   ! TODO: flip each velocity and map back onto the grid.
                   ! TODO: uniform grid with x_vel = 0 at a node would map exactly as a reflection.

                else

                   ! TODO: reduce the current density how? phi_conv = phi_conv + phi - cfl*phi?

                end if
             end if

             if( nx .eq. nx_space )then
                if( x_vel .gt. zero )then

                   phi_conv(nx)%value(i,j,k) = phi(nx)%value(i,j,k) - &
                        factor*( phi(nx)%value(i,j,k) - phi(nx-1)%value(i,j,k) )

                   if( r_modes .gt. 0 )then
                      do l = 1, r_levels
                         phi_conv(nx)%rot(l,i,j,k) = phi(nx)%rot(l,i,j,k) - &
                              factor*( phi(nx)%rot(l,i,j,k) - phi(nx-1)%rot(l,i,j,k) )
                      end do
                   end if

                   if( r_modes .gt. 0 )then
                      do l = 1, v_levels
                         phi_conv(nx)%vib(l,i,j,k) = phi(nx)%vib(l,i,j,k) - &
                              factor*( phi(nx)%vib(l,i,j,k) - phi(nx-1)%vib(l,i,j,k) )
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

  subroutine apply_diffuse_reflection( phi, phi_conv, species, molecule, vel_grid, cfl, nx )

    use DistFunc
    use VelocityGrid
    use PhysicalGrid
    use SpeciesAndReferenceData
    use Remapping

    implicit none

    type(DistFuncType), dimension(:), intent(in) :: phi
    type(VelocityGridType), intent(in) :: vel_grid
    type(MoleculeType), intent(in) :: molecule

    double precision, intent(in) :: cfl

    integer, intent(in) :: nx, species

    type(DistFuncType), dimension(:) :: phi_conv

    double precision, dimension(:), allocatable :: rot_levels, rot_df, vib_levels, vib_df

    integer :: nx_space, ny_space
    integer :: r_modes, v_modes
    integer :: r_levels, v_levels

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: i, j, k, l

    double precision :: mass
    double precision :: dens, kin_temp, u, v, w, rot_temp, vib_temp
    double precision :: maxwell_coeff
    double precision :: computed_dens, maxwell_value

    double precision :: x, y, z
    double precision :: beta_x, beta_y, beta_z, beta3

    integer :: index_1, index_2

    double precision :: opp_vel, factor

    type(MappingResultType) :: mapping

    call get_nspace( nx_space, ny_space )

    r_modes  = molecule%rot_modes
    v_modes  = molecule%vib_modes

    ! phi_conv must be zeroed first for diffuse reflection
    phi_conv(nx)%value = zero


    ! Get boundary properties
    mass = molecule%mass

    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    ! Wall is normal to the x-direction in this version of DVM
    v = zero
    w = zero

    if( nx .eq. 1 )then
       call get_left_wall_temp( kin_temp, species )
       call get_left_wall_velocity( u, v, w, species )
       i_min = vel_grid%i_min
       i_max = vel_grid%i_neg
       
    else if( nx .eq. nx_space )then
       call get_right_wall_temp( kin_temp, species )
       call get_right_wall_velocity( u, v, w, species )
       i_min = vel_grid%i_pos
       i_max = vel_grid%i_max

    else
       write(*,*) "Error. nx is not on a boundary. Given value is: ",nx
       stop

    end if
    
    ! Rotational and vibrational temperatures are set to the wall temperature
    if( r_modes .gt. 0 ) rot_temp = kin_temp
    if( v_modes .gt. 0 ) vib_temp = kin_temp

    ! Calculate flux into the wall
    dens = zero
    do k = k_min, k_max
       do j = j_min,j_max
          do i = i_min,i_max
             dens = dens + phi(nx)%value(i,j,k)
          end do
       end do
    end do
    dens = acc_coeff * dens

    ! Calculate post-reflection half-Maxwellian  - we assume vibrational and rotational temperatures remain 
    !                                              constant but the energy becomes proportionally distributed 
    !                                              throughout the outgoing flux
    if( r_modes .gt. 0 )then
       r_levels = phi_conv(nx)%num_rot_levels
       allocate( rot_levels(1:r_levels), rot_df(1:r_levels) )
       call compute_rot_distribution( rot_df, rot_levels, molecule, rot_temp, r_levels, species )
    end if

    if( v_modes .gt. 0 )then
       v_levels = phi_conv(nx)%num_vib_levels
       allocate( vib_levels(1:v_levels), vib_df(1:v_levels) )
       call compute_vib_distribution( vib_df, vib_levels, molecule, vib_temp, v_levels, species )
    end if

    ! Determine start and end indices in x direction
    if( nx .eq. 1 )then
       i_min = vel_grid%i_pos
       i_max = vel_grid%i_max
    else
       i_min = vel_grid%i_min
       i_max = vel_grid%i_neg
    end if

    maxwell_coeff = sqrt( mass * mass * mass / ( pi * pi * pi * kin_temp * kin_temp * kin_temp ) )
    computed_dens = zero

    ! Diffuse reflection of a distribution function is a half Maxwellian at the wall temperature and velocity 
    ! with density equal to the flux of density into the wall
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

             call compute_maxwellian( x, y, z, maxwell_coeff, mass, u, v, w, kin_temp, maxwell_value )

             phi_conv(nx)%value(i,j,k) = maxwell_value * beta3
             computed_dens = computed_dens + maxwell_value * beta3

             if( r_modes .gt. 0 ) phi_conv(nx)%rot(:,i,j,k) = maxwell_value*beta3*rot_df
             if( v_modes .gt. 0 ) phi_conv(nx)%vib(:,i,j,k) = maxwell_value*beta3*vib_df

          end do
       end do
    end do

    phi_conv(nx)%value = ( dens / computed_dens ) * phi_conv(nx)%value

    ! Specular reflection and 1st order upwind finite difference
    if( nx .eq. 1 )then
       i_min = vel_grid%i_min
       i_max = vel_grid%i_neg
       index_1 = 1
       index_2 = 2
    else
       i_min = vel_grid%i_pos
       i_max = vel_grid%i_max
       index_1 = nx_space
       index_2 = nx_space - 1
    end if

    do k = k_min, k_max
       z = vel_grid%z(k)

       do j = j_min,j_max
          y = vel_grid%y(j)

          do i = i_min,i_max
             x = vel_grid%x(i)
             factor = cfl * x

             ! Finite difference of distribution moving towards the wall
             phi_conv(index_1)%value(i,j,k) = phi(index_1)%value(i,j,k) - &
                  abs(factor)*( phi(index_1)%value(i,j,k) - phi(index_2)%value(i,j,k) )

             if( r_modes .gt. 0 )then
                do l = 1, r_levels
                   phi_conv(index_1)%rot(l,i,j,k) = phi(index_1)%rot(l,i,j,k) - &
                        abs(factor)*( phi(index_1)%rot(l,i,j,k) - phi(index_2)%rot(l,i,j,k) )
                end do
             end if

             if( r_modes .gt. 0 )then
                do l = 1, v_levels
                   phi_conv(index_1)%vib(l,i,j,k) = phi(index_1)%vib(l,i,j,k) - &
                        abs(factor)*( phi(index_1)%vib(l,i,j,k) - phi(index_2)%vib(l,i,j,k) )
                end do
             end if

             ! Specular reflection
             opp_vel = -x
             
             dens = ( one - acc_coeff ) * phi(nx)%value(i,j,k)
             if( r_modes .gt. 0 ) rot_df = ( one - acc_coeff ) * phi(nx)%rot(:,i,j,k)
             if( v_modes .gt. 0 ) vib_df = ( one - acc_coeff ) * phi(nx)%vib(:,i,j,k)

             call perform_remapping( mapping, opp_vel, y, z, one, vel_grid )
    
             call apply_bc_remapping( phi_conv(nx), mapping, dens )
             if( r_modes .gt. 0 ) call apply_bc_rot_remapping( phi_conv(nx), mapping, rot_df )
             if( v_modes .gt. 0 ) call apply_bc_vib_remapping( phi_conv(nx), mapping, vib_df )

          end do
       end do
    end do

    if( r_modes .gt. 0 ) deallocate( rot_levels, rot_df )
    if( v_modes .gt. 0 ) deallocate( vib_levels, vib_df )

    return
  end subroutine apply_diffuse_reflection

  subroutine apply_fixed_in_out_flow( phi_conv, nx, species, molecule, vel_grid )

    use DistFunc
    use VelocityGrid
    use PhysicalGrid
    use SpeciesAndReferenceData

    implicit none

    type(DistFuncType), dimension(:) :: phi_conv

    type(VelocityGridType), intent(in) :: vel_grid
    type(MoleculeType), intent(in) :: molecule

    integer, intent(in) :: nx, species

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
    integer :: nx_space, ny_space
    integer :: mol_shape

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: r_levels, v_levels
    integer :: i, j, k, l
    
    call get_nspace( nx_space, ny_space )

    ! Get internal energy structure
    r_modes  = molecule%rot_modes
    v_modes  = molecule%vib_modes

    mol_shape = molecule%molecule_type
    even_spin = molecule%even_spin
    odd_spin  = molecule%odd_spin

    ! Get boundary conditions
    if( nx .eq. 1 )then
       call get_left_wall_density( dens, species )
       call get_left_wall_temp( temp, species )
       call get_left_wall_velocity( u, v, w, species )
       v = zero
       w = zero
       if( r_modes .gt. 0 ) call get_left_wall_temp_rot( temp_rot, species )
       if( v_modes .gt. 0 ) call get_left_wall_temp_vib( temp_vib, species )
       
    else if( nx .eq. nx_space )then
       call get_right_wall_density( dens, species )
       call get_right_wall_temp( temp, species )
       call get_right_wall_velocity( u, v, w, species )
       v = zero
       w = zero
       if( r_modes .gt. 0 ) call get_right_wall_temp_rot( temp_rot, species )
       if( v_modes .gt. 0 ) call get_right_wall_temp_vib( temp_vib, species )

    else
       write(*,*) "Error. nx is not on a boundary. Given value is: ",nx
       stop

    end if

    if( r_modes .gt. 0 )then
       r_levels = phi_conv(nx)%num_rot_levels
       allocate( rot_levels(1:r_levels) )
       allocate( rot_df(1:r_levels) )
       call compute_rot_distribution( rot_df, rot_levels, molecule, temp_rot, r_levels, species )
    end if

    if( v_modes .gt. 0 )then
       v_levels = phi_conv(nx)%num_vib_levels
       allocate( vib_levels(1:v_levels) )
       allocate( vib_df(1:v_levels) )
       call compute_vib_distribution( vib_df, vib_levels, molecule, temp_vib, v_levels, species )
    end if
    
    theta_r = molecule%theta_r
    theta_v = molecule%theta_v

    mass = molecule%mass

    maxwell_coeff = dens * sqrt( mass * mass * mass / ( pi * pi * pi * temp * temp * temp ) )
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

             phi_conv(nx)%value(i,j,k) = maxwell_value*beta3

             if( r_modes .gt. 0 ) phi_conv(nx)%rot(:,i,j,k) = maxwell_value*beta3*rot_df

             if( v_modes .gt. 0 ) phi_conv(nx)%vib(:,i,j,k) = maxwell_value*beta3*vib_df

          end do
       end do
    end do

    if( r_modes .gt. 0 ) deallocate( rot_levels, rot_df )
    if( v_modes .gt. 0 ) deallocate( vib_levels, vib_df )

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

  subroutine apply_bc_remapping( phi, mapping, dens )

    use DistFunc
    use Remapping

    implicit none

    type(MappingResultType), intent(in) :: mapping
    double precision, intent(in) :: dens

    type(DistFuncType) :: phi

    double precision :: repl_o, repl_ix, repl_iy, repl_iz, repl_ex, repl_ey, repl_ez
    double precision :: repl_xyz

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
    phi%value( i, j, k )     = phi%value( i, j, k ) + dens * repl_o
    phi%value( i + a, j, k ) = phi%value( i + a, j, k ) + dens * repl_ix
    phi%value( i, j + b, k ) = phi%value( i, j + b, k ) + dens * repl_iy
    phi%value( i, j, k + c ) = phi%value( i, j, k + c ) + dens * repl_iz

    ! Exterior points
    select case( ext_flag )
    case( x_only )
       phi%value( i - a, j, k ) = phi%value( i - a, j, k ) + dens * repl_ex

    case( y_only )
       phi%value( i, j - b, k ) = phi%value( i, j - b, k ) + dens * repl_ey

    case( z_only )
       phi%value( i, j, k - c ) = phi%value( i, j, k - c ) + dens * repl_ez

    case( xy_only )
       phi%value( i - a, j, k ) = phi%value( i - a, j, k ) + dens * repl_ex
       phi%value( i, j - b, k ) = phi%value( i, j - b, k ) + dens * repl_ey

    case( xz_only )
       phi%value( i - a, j, k ) = phi%value( i - a, j, k ) + dens * repl_ex
       phi%value( i, j, k - c ) = phi%value( i, j, k - c ) + dens * repl_ez

    case( yz_only )
       phi%value( i, j - b, k ) = phi%value( i, j - b, k ) + dens * repl_ey
       phi%value( i, j, k - c ) = phi%value( i, j, k - c ) + dens * repl_ez

    case( xyz )
       phi%value( i - a, j, k ) = phi%value( i - a, j, k ) + dens * repl_ex
       phi%value( i, j - b, k ) = phi%value( i, j - b, k ) + dens * repl_ey
       phi%value( i, j, k - c ) = phi%value( i, j, k - c ) + dens * repl_ez

    case default
       !TODO: Consistent error handling
       write(*,*) "Error: Invalid value of ext_flag - ",&
            ext_flag
       stop

    end select

    return
  end subroutine apply_bc_remapping

  subroutine apply_bc_rot_remapping( phi, mapping, rot_df )
    
    use DistFunc
    use Remapping

    implicit none

    type(MappingResultType), intent(in) :: mapping
    double precision, dimension(:), intent(in) :: rot_df

    type(DistFuncType) :: phi

    double precision :: repl_o, repl_ix, repl_iy, repl_iz, repl_ex, repl_ey, repl_ez
    double precision :: repl_xyz

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
    phi%rot( :, i, j, k )     = phi%rot( :, i,j,k ) + repl_o * rot_df
    phi%rot( :, i + a, j, k ) = phi%rot( :, i + a, j, k ) + repl_ix * rot_df
    phi%rot( :, i, j + b, k ) = phi%rot( :, i, j + b, k ) + repl_iy * rot_df
    phi%rot( :, i, j, k + c ) = phi%rot( :, i, j, k + c ) + repl_iz * rot_df

    ! Exterior points
    select case( ext_flag )
    case( x_only )
       phi%rot( :, i - a, j, k ) = phi%rot( :, i - a, j, k ) + repl_ex * rot_df

    case( y_only )
       phi%rot( :, i, j - b, k ) = phi%rot( :, i, j - b, k ) + repl_ey * rot_df

    case( z_only )
       phi%rot( :, i, j, k - c ) = phi%rot( :, i, j, k - c ) + repl_ez * rot_df

    case( xy_only )
       phi%rot( :, i - a, j, k ) = phi%rot( :, i - a, j, k ) + repl_ex * rot_df
       phi%rot( :, i, j - b, k ) = phi%rot( :, i, j - b, k ) + repl_ey * rot_df

    case( xz_only )
       phi%rot( :, i - a, j, k ) = phi%rot( :, i - a, j, k ) + repl_ex * rot_df
       phi%rot( :, i, j, k - c ) = phi%rot( :, i, j, k - c ) + repl_ez * rot_df

    case( yz_only )
       phi%rot( :, i, j - b, k ) = phi%rot( :, i, j - b, k ) + repl_ey * rot_df
       phi%rot( :, i, j, k - c ) = phi%rot( :, i, j, k - c ) + repl_ez * rot_df

    case( xyz )
       phi%rot( :, i - a, j, k ) = phi%rot( :, i - a, j, k ) + repl_ex * rot_df
       phi%rot( :, i, j - b, k ) = phi%rot( :, i, j - b, k ) + repl_ey * rot_df
       phi%rot( :, i, j, k - c ) = phi%rot( :, i, j, k - c ) + repl_ez * rot_df


    case default
       !TODO: Consistent error handling
       write(*,*) "Error: Invalid value of ext_flag - ",&
            ext_flag
       stop

    end select

    return
  end subroutine apply_bc_rot_remapping

  subroutine apply_bc_vib_remapping( phi, mapping, vib_df )

    use DistFunc
    use Remapping

    implicit none

    type(MappingResultType), intent(in) :: mapping
    double precision, dimension(:), intent(in) :: vib_df

    type(DistFuncType) :: phi

    double precision :: repl_o, repl_ix, repl_iy, repl_iz, repl_ex, repl_ey, repl_ez
    double precision :: repl_xyz

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
    phi%vib( :, i, j, k )     = phi%vib( :, i,j,k ) + repl_o * vib_df
    phi%vib( :, i + a, j, k ) = phi%vib( :, i + a, j, k ) + repl_ix * vib_df
    phi%vib( :, i, j + b, k ) = phi%vib( :, i, j + b, k ) + repl_iy * vib_df
    phi%vib( :, i, j, k + c ) = phi%vib( :, i, j, k + c ) + repl_iz * vib_df

    ! Exterior points
    select case( ext_flag )
    case( x_only )
       phi%vib( :, i - a, j, k ) = phi%vib( :, i - a, j, k ) + repl_ex * vib_df

    case( y_only )
       phi%vib( :, i, j - b, k ) = phi%vib( :, i, j - b, k ) + repl_ey * vib_df

    case( z_only )
       phi%vib( :, i, j, k - c ) = phi%vib( :, i, j, k - c ) + repl_ez * vib_df

    case( xy_only )
       phi%vib( :, i - a, j, k ) = phi%vib( :, i - a, j, k ) + repl_ex * vib_df
       phi%vib( :, i, j - b, k ) = phi%vib( :, i, j - b, k ) + repl_ey * vib_df

    case( xz_only )
       phi%vib( :, i - a, j, k ) = phi%vib( :, i - a, j, k ) + repl_ex * vib_df
       phi%vib( :, i, j, k - c ) = phi%vib( :, i, j, k - c ) + repl_ez * vib_df

    case( yz_only )
       phi%vib( :, i, j - b, k ) = phi%vib( :, i, j - b, k ) + repl_ey * vib_df
       phi%vib( :, i, j, k - c ) = phi%vib( :, i, j, k - c ) + repl_ez * vib_df

    case( xyz )
       phi%vib( :, i - a, j, k ) = phi%vib( :, i - a, j, k ) + repl_ex * vib_df
       phi%vib( :, i, j - b, k ) = phi%vib( :, i, j - b, k ) + repl_ey * vib_df
       phi%vib( :, i, j, k - c ) = phi%vib( :, i, j, k - c ) + repl_ez * vib_df


    case default
       !TODO: Consistent error handling
       write(*,*) "Error: Invalid value of ext_flag - ",&
            ext_flag
       stop

    end select

    return
  end subroutine apply_bc_vib_remapping

end module BoundaryConditions
