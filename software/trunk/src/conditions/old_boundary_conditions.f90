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
  ! Added use DistFunc for boundary nodes
  use DistFunc

  implicit none

  private

  integer :: boundary_condition_LW, boundary_condition_RW
  integer :: boundary_condition_BW, boundary_condition_TW

  ! TODO: change to a species dependent variable
  double precision :: acc_coeff
  
  ! TODO: adding evaporation/absorption boundary condition
  logical, allocatable, dimension(:) :: evaporation_flag
  double precision, allocatable, dimension(:) :: evaporation_rate

  double precision, allocatable, dimension(:) :: density_LW, u_LW, v_LW, w_LW, temp_LW
  double precision, allocatable, dimension(:) :: density_RW, u_RW, v_RW, w_RW, temp_RW
  double precision, allocatable, dimension(:) :: density_BW, u_BW, v_BW, w_BW, temp_BW
  double precision, allocatable, dimension(:) :: density_TW, u_TW, v_TW, w_TW, temp_TW
  double precision, allocatable, dimension(:) :: temp_rot_LW, temp_rot_RW, temp_rot_BW, temp_rot_TW
  double precision, allocatable, dimension(:) :: temp_vib_LW, temp_vib_RW, temp_vib_BW, temp_vib_TW

  ! ADDING BOUNDARY NODES HERE
  type(DistFuncType), dimension(4) :: boundary_nodes

  integer, parameter :: wall_reflection     = 1
  integer, parameter :: fixed_in_out_flow   = 2
  integer, parameter :: zero_gradient_flow  = 3

  integer, parameter :: x_only  = 1
  integer, parameter :: y_only  = 2
  integer, parameter :: z_only  = 4
  integer, parameter :: xy_only = 3
  integer, parameter :: xz_only = 5
  integer, parameter :: yz_only = 6
  integer, parameter :: xyz     = 7

  integer, parameter :: left   = 1
  integer, parameter :: right  = 2
  integer, parameter :: bottom = 3
  integer, parameter :: top    = 4

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
  public :: get_top_wall_temp
  public :: get_top_wall_density
  public :: get_top_wall_velocity
  public :: get_top_wall_temp_rot
  public :: get_top_wall_temp_vib
  public :: get_bottom_wall_temp
  public :: get_bottom_wall_density
  public :: get_bottom_wall_velocity
  public :: get_bottom_wall_temp_rot
  public :: get_bottom_wall_temp_vib
  public :: apply_boundary_conditions

contains

  subroutine set_boundary_conditions_flags( left_wall_bc_in, right_wall_bc_in, &
       bottom_wall_bc_in, top_wall_bc_in )

    implicit none

    integer, intent(in) :: left_wall_bc_in, right_wall_bc_in
    integer, intent(in) :: bottom_wall_bc_in, top_wall_bc_in

    boundary_condition_LW = left_wall_bc_in
    boundary_condition_RW = right_wall_bc_in
    boundary_condition_BW = bottom_wall_bc_in
    boundary_condition_TW = top_wall_bc_in

    return
  end subroutine set_boundary_conditions_flags

  subroutine get_boundary_conditions_flags( left_wall_bc_out, right_wall_bc_out, &
       bottom_wall_bc_out, top_wall_bc_out )

    implicit none

    integer :: left_wall_bc_out, right_wall_bc_out
    integer :: bottom_wall_bc_out, top_wall_bc_out

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

    allocate( density_BW( num_species ), STAT=status )
    call allocate_error_check( status, "density_BW" )
    allocate( u_BW( num_species ), STAT=status )
    call allocate_error_check( status, "u_BW" )
    allocate( v_BW( num_species ), STAT=status )
    call allocate_error_check( status, "v_BW" )
    allocate( w_BW( num_species ), STAT=status )
    call allocate_error_check( status, "w_BW" )
    allocate( temp_BW( num_species ), STAT=status )
    call allocate_error_check( status, "temp_BW" )
    allocate( temp_rot_BW( num_species ), STAT=status )
    call allocate_error_check( status, "temp_rot_BW" )
    allocate( temp_vib_BW( num_species ), STAT=status )
    call allocate_error_check( status, "temp_vib_BW" )

    allocate( density_TW( num_species ), STAT=status )
    call allocate_error_check( status, "density_TW" )
    allocate( u_TW( num_species ), STAT=status )
    call allocate_error_check( status, "u_TW" )
    allocate( v_TW( num_species ), STAT=status )
    call allocate_error_check( status, "v_TW" )
    allocate( w_TW( num_species ), STAT=status )
    call allocate_error_check( status, "w_TW" )
    allocate( temp_TW( num_species ), STAT=status )
    call allocate_error_check( status, "temp_TW" )
    allocate( temp_rot_TW( num_species ), STAT=status )
    call allocate_error_check( status, "temp_rot_TW" )
    allocate( temp_vib_TW( num_species ), STAT=status )
    call allocate_error_check( status, "temp_vib_TW" )

    return
  end subroutine initialize_property_arrays
  
  subroutine set_boundary_properties( density_LW_in, u_LW_in, v_LW_in, w_LW_in, temp_LW_in, &
       density_RW_in, u_RW_in, v_RW_in, w_RW_in, temp_RW_in, &
       density_BW_in, u_BW_in, v_BW_in, w_BW_in, temp_BW_in, &
       density_TW_in, u_TW_in, v_TW_in, w_TW_in, temp_TW_in, species )

    implicit none

    double precision, intent(in) :: density_LW_in, u_LW_in, v_LW_in, w_LW_in, temp_LW_in
    double precision, intent(in) :: density_RW_in, u_RW_in, v_RW_in, w_RW_in, temp_RW_in
    double precision, intent(in) :: density_BW_in, u_BW_in, v_BW_in, w_BW_in, temp_BW_in
    double precision, intent(in) :: density_TW_in, u_TW_in, v_TW_in, w_TW_in, temp_TW_in
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

    density_BW( species ) = density_BW_in
    u_BW( species )       = u_BW_in
    v_BW( species )       = v_BW_in
    w_BW( species )       = w_BW_in
    temp_BW( species )    = temp_BW_in

    density_TW( species ) = density_TW_in  
    u_TW( species )       = u_TW_in
    v_TW( species )       = v_TW_in
    w_TW( species )       = w_TW_in
    temp_TW( species )    = temp_TW_in
    
    return
  end subroutine set_boundary_properties

  subroutine set_temp_rot( temp_rot_LW_in, temp_rot_RW_in, &
       temp_rot_BW_in, temp_rot_TW_in, species )

    implicit none

    double precision, intent(in) :: temp_rot_LW_in, temp_rot_RW_in
    double precision, intent(in) :: temp_rot_BW_in, temp_rot_TW_in
    integer, intent(in) :: species

    temp_rot_LW( species ) = temp_rot_LW_in
    temp_rot_RW( species ) = temp_rot_RW_in
    temp_rot_BW( species ) = temp_rot_BW_in
    temp_rot_TW( species ) = temp_rot_TW_in

    return
  end subroutine set_temp_rot

  subroutine set_temp_vib( temp_vib_LW_in, temp_vib_RW_in, &
       temp_vib_BW_in, temp_vib_TW_in,species )

    implicit none

    double precision, intent(in) :: temp_vib_LW_in, temp_vib_RW_in
    double precision, intent(in) :: temp_vib_BW_in, temp_vib_TW_in
    integer, intent(in) :: species

    temp_vib_LW( species ) = temp_vib_LW_in
    temp_vib_RW( species ) = temp_vib_RW_in
    temp_vib_BW( species ) = temp_vib_BW_in
    temp_vib_TW( species ) = temp_vib_TW_in

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

  subroutine get_top_wall_temp( temp_TW_set, species )

    implicit none

    double precision :: temp_TW_set
    integer, intent(in) :: species

    temp_TW_set = temp_TW( species )
    return
  end subroutine get_top_wall_temp

  subroutine get_bottom_wall_temp( temp_BW_set, species )

    implicit none

    double precision :: temp_BW_set
    integer, intent(in) :: species

    temp_BW_set = temp_BW( species )
    return
  end subroutine get_bottom_wall_temp

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

  subroutine get_top_wall_density( density_TW_set, species )

    implicit none

    double precision :: density_TW_set
    integer, intent(in) :: species

    density_TW_set = density_TW( species )
    return
  end subroutine get_top_wall_density

  subroutine get_bottom_wall_density( density_BW_set, species )

    implicit none

    double precision :: density_BW_set
    integer, intent(in) :: species

    density_BW_set = density_BW( species )
    return
  end subroutine get_bottom_wall_density

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

  subroutine get_top_wall_velocity( u_TW_set, v_TW_set, w_TW_set, species )

    implicit none

    double precision :: u_TW_set, v_TW_set, w_TW_set
    integer, intent(in) :: species

    u_TW_set = u_TW( species )
    v_TW_set = v_TW( species )
    w_TW_set = w_TW( species )
    return
  end subroutine get_top_wall_velocity

  subroutine get_bottom_wall_velocity( u_BW_set, v_BW_set, w_BW_set, species )

    implicit none

    double precision :: u_BW_set, v_BW_set, w_BW_set
    integer, intent(in) :: species

    u_BW_set = u_BW( species )
    v_BW_set = v_BW( species )
    w_BW_set = w_BW( species )
    return
  end subroutine get_bottom_wall_velocity

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

  subroutine get_top_wall_temp_rot( temp_rot_TW_set, species )

    implicit none

    double precision :: temp_rot_TW_set
    integer, intent(in) :: species

    temp_rot_TW_set = temp_rot_TW( species )
    return
  end subroutine get_top_wall_temp_rot

  subroutine get_bottom_wall_temp_rot( temp_rot_BW_set, species )

    implicit none

    double precision :: temp_rot_BW_set
    integer, intent(in) :: species

    temp_rot_BW_set = temp_rot_BW( species )
    return
  end subroutine get_bottom_wall_temp_rot

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

  subroutine get_top_wall_temp_vib( temp_vib_TW_set, species )

    implicit none

    double precision :: temp_vib_TW_set
    integer, intent(in) :: species

    temp_vib_TW_set = temp_vib_TW( species )

    return
  end subroutine get_top_wall_temp_vib

  subroutine get_bottom_wall_temp_vib( temp_vib_BW_set, species )

    implicit none

    double precision :: temp_vib_BW_set
    integer, intent(in) :: species

    temp_vib_BW_set = temp_vib_BW( species )

    return
  end subroutine get_bottom_wall_temp_vib

  subroutine apply_boundary_conditions( phi, phi_conv, molecule, vel_grid, species )

    use DistFunc
    use VelocityGrid
    use PhysicalGrid
    use SpeciesAndReferenceData

    implicit none

    type(VelocityGridType), intent(in) :: vel_grid
    type(MoleculeType), intent(in) :: molecule

    type(DistFuncType), dimension(:,:) :: phi, phi_conv

    !double precision, intent(in) :: cfl
    integer, intent(in) :: species

    integer :: nx_space, ny_space
    integer :: nx, ny, index_1, index_2

    call get_nspace( nx_space, ny_space )

    do nx = 1, nx_space
       ! Bottom wall
       ny = 1
       call apply_wall_reflection( phi, phi_conv, species, molecule, vel_grid, nx, ny )
       ! Top wall
       ny = ny_space
       call apply_wall_reflection( phi, phi_conv, species, molecule, vel_grid, nx, ny )
    end do
    do ny = 2, ny_space-1
       ! Left wall
       nx = 1
       call apply_wall_reflection( phi, phi_conv, species, molecule, vel_grid, nx, ny )
       ! Right wall
       nx = nx_space
       call apply_wall_reflection( phi, phi_conv, species, molecule, vel_grid, nx, ny )
    end do
    
!    select case( boundary_condition_LW )
!    case( wall_reflection )
!       nx = 1
!       call apply_wall_reflection( phi, phi_conv, species, molecule, vel_grid, nx )
!
!    case( fixed_in_out_flow )
!       nx = 1
!       call apply_fixed_in_out_flow( phi_conv, nx, species, molecule, vel_grid )
!
!    case( zero_gradient_flow )
!       index_1 = 1
!       index_2 = 2
!       call apply_zero_gradient_flow( phi_conv, molecule, vel_grid, index_1, index_2 )
!
!    case default
!       write(*,*) "Error: Invalid value of boundary_condition_LW - :", &
!            boundary_condition_LW
!    end select
!
!    select case( boundary_condition_RW )
!    case( wall_reflection )
!       nx = nx_space
!       call apply_wall_reflection( phi, phi_conv, species, molecule, vel_grid, nx )
!
!    case( fixed_in_out_flow )
!       nx = nx_space
!       call apply_fixed_in_out_flow( phi_conv, nx, species, molecule, vel_grid )
!
!    case( zero_gradient_flow )
!       index_1 = nx_space
!       index_2 = nx_space-1
!       call apply_zero_gradient_flow( phi_conv, molecule, vel_grid, index_1, index_2 )
!
!    case default
!       write(*,*) "Error: Invalid value of boundary_condition_RW - :", &
!            boundary_condition_RW
!    end select

    return
  end subroutine apply_boundary_conditions

  subroutine apply_wall_reflection( phi, phi_conv, species, molecule, vel_grid, nx, ny )
    
    use DistFunc
    use VelocityGrid
    use SpeciesAndReferenceData
    use PhysicalGrid
    use TimeStepping
    use Remapping

    implicit none

    type(VelocityGridType), intent(in) :: vel_grid
    type(MoleculeType), intent(in) :: molecule

    integer, intent(in) :: nx, ny, species

    type(DistFuncType), dimension(:,:) :: phi, phi_conv

    double precision, dimension(:), allocatable :: rot_levels, rot_df, vib_levels, vib_df

    double precision :: evaporation
    double precision :: x_factor, y_factor

    integer :: nx_space, ny_space
    integer :: index, wall_index
    integer :: x_in_min, x_in_max, x_out_min, x_out_max
    integer :: y_in_min, y_in_max, y_out_min, y_out_max
    integer :: i_zero, j_zero, i_min, i_max, j_min, j_max, k_min, k_max
    integer :: i, j, k, i_ref, j_ref
    integer :: x_side, y_side
    logical :: i_zero_flag = .false.
    logical :: j_zero_flag = .false.

    integer :: r_modes, v_modes, r_levels, v_levels

    double precision :: dt, dx, dy

    double precision :: mass, kin_temp, u, v, w
    double precision :: dens_in, dens_out, mom_in, mom_out
    double precision :: flux_in, flux_out
    double precision :: force_on_wall
    double precision :: maxwell_coeff, maxwell_value
    double precision :: dens, computed_dens

    type(MappingResultType) :: mapping
    double precision :: x, y, z
    double precision :: beta_x, beta_y, beta_z, beta3
    double precision :: x_cfl, y_cfl

    ! Molecular properties
    r_modes = molecule%rot_modes
    v_modes = molecule%vib_modes

    mass = molecule%mass

    ! The convection distribution must be set to zero initially for this boundary
    phi_conv(nx,ny)%value = zero
    if( r_modes .gt. zero ) phi_conv(nx,ny)%rot = zero
    if( v_modes .gt. zero ) phi_conv(nx,ny)%vib = zero

    ! Initialize the boundary nodes and set them to zero
    do i = 1, 4
       call create_dist_func( boundary_nodes(i), vel_grid, molecule, species )
       boundary_nodes(i)%value = zero
       if( r_modes .gt. zero ) boundary_nodes(i)%rot = zero
       if( v_modes .gt. zero ) boundary_nodes(i)%vib = zero
    end do

    ! Time step
    call get_deltat( dt )

    ! Physical Grid
    call get_delta_x( dx, dy )
    call get_nspace( nx_space, ny_space )

    ! Get velocity domain
    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    i_zero  = vel_grid%i_zero
    j_zero  = vel_grid%j_zero
    

    ! Check for existance of zero velocity point
    if( i_zero .ne. vel_grid%i_min .and. i_zero .ne. vel_grid%i_neg .and. &
         i_zero .ne. vel_grid%i_pos .and. i_zero .ne. vel_grid%i_max )then
       i_zero_flag = .true.
    end if

    if( j_zero .ne. vel_grid%j_min .and. j_zero .ne. vel_grid%j_neg .and. &
         j_zero .ne. vel_grid%j_pos .and. j_zero .ne. vel_grid%j_max )then
       j_zero_flag = .true.
    end if

!    else
!       write(*,*) "Error: nx is not on a boundary. Given values is: ", nx
!       stop
!
!    end if

    ! CFL condition factor
    x_factor = dt/dx
    y_factor = dt/dy

    ! Evaporation flux from the wall [#/m^3/s]
    evaporation = zero

!    ! Calculate properties into the wall
!    dens_in = zero
!    mom_in  = zero
!    flux_in = zero
!    do k = k_min, k_max
!       do j = j_min, j_max
!          do i = in_min, in_max
!             x = vel_grid%x(i)
!             cfl = factor * x
!             dens_in = dens_in + phi(nx)%value(i,j,k)
!             mom_in  = mom_in + x * phi(nx)%value(i,j,k)
!             flux_in = flux_in + abs(x) * acc_coeff * phi(nx)%value(i,j,k) 
!          end do
!       end do
!    end do
!
!    ! DIFFUSE REFLECTION
!    !---------------------------------------------------------------------------------------------------------
!    if( acc_coeff .gt. zero .or. evaporation .ne. zero )then
!       ! Rotational and vibrational temperatures are set to the wall temperature
!       if( r_modes .gt. 0 )then
!          r_levels = phi(nx)%num_rot_levels
!          allocate( rot_levels(1:r_levels), rot_df(1:r_levels) )
!          call compute_rot_distribution( rot_df, rot_levels, molecule, kin_temp, r_levels, species )
!       end if
!       
!       if( v_modes .gt. 0 )then
!          v_levels = phi(nx)%num_vib_levels
!          allocate( vib_levels(1:v_levels), vib_df(1:v_levels) )
!          call compute_vib_distribution( vib_df, vib_levels, molecule, kin_temp, v_levels, species )
!       end if
!       
!       ! Diffuse Reflection
!       dens = acc_coeff * dens_in! + dt * evaporation
!       maxwell_coeff = sqrt( mass * mass * mass / ( pi * pi * pi * kin_temp * kin_temp * kin_temp ) )
!       computed_dens = zero
!       flux_out = zero
!
!       do k = k_min, k_max
!          z = vel_grid%z(k)
!          beta_z = vel_grid%beta_z(k)
!          
!          do j = j_min, j_max
!             y = vel_grid%y(j)
!             beta_y = vel_grid%beta_y(j)
!
!             do i = out_min, out_max
!                x = vel_grid%x(i)
!                beta_x = vel_grid%beta_x(i)
!                beta3 = beta_x*beta_y*beta_z
!
!                ! TODO: should only have to do this once
!                call compute_maxwellian( x, y, z, maxwell_coeff, mass, u, v, w, kin_temp, maxwell_value )
!
!                phi_conv(nx)%value(i,j,k) = maxwell_value * beta3
!
!                flux_out = flux_out + abs(x) * maxwell_value * beta3
!                computed_dens = computed_dens + phi_conv(nx)%value(i,j,k)
!
!                if( r_modes .gt. 0 ) phi_conv(nx)%rot(:,i,j,k) = maxwell_value * beta3 * rot_df
!                if( v_modes .gt. 0 ) phi_conv(nx)%vib(:,i,j,k) = maxwell_value * beta3 * vib_df
!
!             end do
!          end do
!       end do
!
!       phi_conv(nx)%value = ( flux_in / flux_out ) * phi_conv(nx)%value
!       if( r_modes .gt. 0 ) phi_conv(nx)%rot = ( flux_in / flux_out ) * phi_conv(nx)%rot
!       if( v_modes .gt. 0 ) phi_conv(nx)%vib = ( flux_in / flux_out ) * phi_conv(nx)%vib
!
!       
!
!    end if

    ! Left wall
    !=====================================================================================================
    if( nx .eq. 1 )then
!       print *, "Left boundary node gets wall props"
       ! Get wall properties
       call get_left_wall_temp( kin_temp, species )
       call get_left_wall_velocity( u, v, w, species )
       ! Calculate properties into the wall
       dens_in = zero
       mom_in  = zero
       flux_in = zero
       do k = k_min, k_max
          do j = j_min, j_max
             do i = vel_grid%i_min, vel_grid%i_neg
                x = vel_grid%x(i)
                dens_in = dens_in + phi(nx,ny)%value(i,j,k)
                mom_in  = mom_in + x * phi(nx,ny)%value(i,j,k)
                flux_in = flux_in + abs(x) * acc_coeff * phi(nx,ny)%value(i,j,k) 
             end do
          end do
       end do

       ! DIFFUSE REFLECTION
       !---------------------------------------------------------------------------------------------------------
       if( acc_coeff .gt. zero .or. evaporation .ne. zero )then
          ! Rotational and vibrational temperatures are set to the wall temperature
          if( r_modes .gt. 0 )then
             r_levels = phi(nx,ny)%num_rot_levels
             allocate( rot_levels(1:r_levels), rot_df(1:r_levels) )
             call compute_rot_distribution( rot_df, rot_levels, molecule, kin_temp, r_levels, species )
          end if
          
          if( v_modes .gt. 0 )then
             v_levels = phi(nx,ny)%num_vib_levels
             allocate( vib_levels(1:v_levels), vib_df(1:v_levels) )
             call compute_vib_distribution( vib_df, vib_levels, molecule, kin_temp, v_levels, species )
          end if
          
          ! Diffuse Reflection
          dens = acc_coeff * dens_in
          maxwell_coeff = sqrt( mass * mass * mass / ( pi * pi * pi * kin_temp * kin_temp * kin_temp ) )
          computed_dens = zero
          flux_out = zero

          do k = k_min, k_max
             z = vel_grid%z(k)
             beta_z = vel_grid%beta_z(k)
             
             do j = j_min, j_max
                y = vel_grid%y(j)
                beta_y = vel_grid%beta_y(j)

                do i = vel_grid%i_pos, vel_grid%i_max
                   x = vel_grid%x(i)
                   beta_x = vel_grid%beta_x(i)
                   beta3 = beta_x*beta_y*beta_z

                   ! TODO: should only have to do this once
                   call compute_maxwellian( x, y, z, maxwell_coeff, mass, u, v, w, kin_temp, maxwell_value )

                   boundary_nodes(1)%value(i,j,k) = maxwell_value * beta3

                   flux_out = flux_out + abs(x) * maxwell_value * beta3
                   computed_dens = computed_dens + boundary_nodes(1)%value(i,j,k)

                   if( r_modes .gt. 0 ) boundary_nodes(1)%rot(:,i,j,k) = maxwell_value * beta3 * rot_df
                   if( v_modes .gt. 0 ) boundary_nodes(1)%vib(:,i,j,k) = maxwell_value * beta3 * vib_df

                end do
             end do
          end do

          boundary_nodes(1)%value = ( flux_in / flux_out ) * boundary_nodes(1)%value
          if( r_modes .gt. 0 ) boundary_nodes(1)%rot = ( flux_in / flux_out ) * boundary_nodes(1)%rot
          if( v_modes .gt. 0 ) boundary_nodes(1)%vib = ( flux_in / flux_out ) * boundary_nodes(1)%vib

       end if

       ! SPECULAR REFLECTION
       !---------------------------------------------------------------------------------------------------------
       if( acc_coeff .lt. one )then
          do k = k_min, k_max
             z = vel_grid%z(k)
   
             do j = j_min, j_max
                y = vel_grid%y(j)
   
                do i = vel_grid%i_min, vel_grid%i_neg
                   x = -vel_grid%x(i)
                   x_cfl = x_factor * x
   
                   dens = ( one - acc_coeff ) * phi(nx,ny)%value(i,j,k)
   !                print *, "left",i,j,k,dens
   
                   if( r_modes .gt. 0 ) rot_df = ( one - acc_coeff ) * phi(nx,ny)%rot(:,i,j,k)
                   if( v_modes .gt. 0 ) vib_df = ( one - acc_coeff ) * phi(nx,ny)%vib(:,i,j,k)
   
                   call perform_remapping( mapping, x, y, z, one, vel_grid )
   
                   call apply_bc_remapping( boundary_nodes(1), mapping, dens )
                   if( r_modes .gt. 0 ) call apply_bc_rot_remapping( boundary_nodes(1), mapping, rot_df )
                   if( v_modes .gt. 0 ) call apply_bc_vib_remapping( boundary_nodes(1), mapping, vib_df )
   
                end do
             end do
          end do
       end if
       ! Test
!       print *, "Left wall boundary node"
!       do k = k_min, k_max
!          do j = j_min, j_max
!             do i = vel_grid%i_min, vel_grid%i_neg
!                ! Set mirror index of i assuming i_min=0
!                i_ref = vel_grid%i_max-i
!                print *, i, j, k, phi(nx,ny)%value(i,j,k) - boundary_nodes(1)%value(i_ref,j,k)
!             end do
!          end do
!       end do
    else
!       print *, "Left boundary node gets value of node to the left (not at left wall)"
       ! Not at left wall so boundary node gets value of cell to the left
       boundary_nodes(1)%value=phi(nx-1,ny)%value
       if( i_zero_flag )then
          boundary_nodes(1)%value(i_zero,:,:) = 0
       end if
       ! Test
!       print *, "TOP"
!       do k = k_min, k_max
!          do j = j_min, j_max
!             do i = vel_grid%i_min, vel_grid%i_neg
!                print *, i, j, k, phi(nx-1,ny)%value(i,j,k) - boundary_nodes(1)%value(i,j,k)
!             end do
!             do i = vel_grid%i_pos, vel_grid%i_max
!                print *, i, j, k, phi(nx-1,ny)%value(i,j,k) - boundary_nodes(1)%value(i,j,k)
!             end do
!          end do
!       end do
    end if

    ! Right wall
    !=========================================================================================================
    if( nx .eq. nx_space )then
 !      print *, "Right boundary node gets wall props"
       ! Get wall properties
       call get_right_wall_temp( kin_temp, species )
       call get_right_wall_velocity( u, v, w, species )
       ! Calculate properties into the wall
       dens_in = zero
       mom_in  = zero
       flux_in = zero
       do k = k_min, k_max
          do j = j_min, j_max
             do i = vel_grid%i_pos, vel_grid%i_max
                x = vel_grid%x(i)
                dens_in = dens_in + phi(nx,ny)%value(i,j,k)
                mom_in  = mom_in + x * phi(nx,ny)%value(i,j,k)
                flux_in = flux_in + abs(x) * acc_coeff * phi(nx,ny)%value(i,j,k) 
             end do
          end do
       end do

       ! DIFFUSE REFLECTION
       !---------------------------------------------------------------------------------------------------------
       if( acc_coeff .gt. zero .or. evaporation .ne. zero )then
          ! Rotational and vibrational temperatures are set to the wall temperature
          if( r_modes .gt. 0 )then
             r_levels = phi(nx,ny)%num_rot_levels
             allocate( rot_levels(1:r_levels), rot_df(1:r_levels) )
             call compute_rot_distribution( rot_df, rot_levels, molecule, kin_temp, r_levels, species )
          end if
          
          if( v_modes .gt. 0 )then
             v_levels = phi(nx,ny)%num_vib_levels
             allocate( vib_levels(1:v_levels), vib_df(1:v_levels) )
             call compute_vib_distribution( vib_df, vib_levels, molecule, kin_temp, v_levels, species )
          end if
          
          ! Diffuse Reflection
          dens = acc_coeff * dens_in
          maxwell_coeff = sqrt( mass * mass * mass / ( pi * pi * pi * kin_temp * kin_temp * kin_temp ) )
          computed_dens = zero
          flux_out = zero

          do k = k_min, k_max
             z = vel_grid%z(k)
             beta_z = vel_grid%beta_z(k)
             
             do j = j_min, j_max
                y = vel_grid%y(j)
                beta_y = vel_grid%beta_y(j)

                do i = vel_grid%i_min, vel_grid%i_neg
                   x = vel_grid%x(i)
                   beta_x = vel_grid%beta_x(i)
                   beta3 = beta_x*beta_y*beta_z

                   ! TODO: should only have to do this once
                   call compute_maxwellian( x, y, z, maxwell_coeff, mass, u, v, w, kin_temp, maxwell_value )

                   boundary_nodes(2)%value(i,j,k) = maxwell_value * beta3

                   flux_out = flux_out + abs(x) * maxwell_value * beta3
                   computed_dens = computed_dens + boundary_nodes(2)%value(i,j,k)

                   if( r_modes .gt. 0 ) boundary_nodes(2)%rot(:,i,j,k) = maxwell_value * beta3 * rot_df
                   if( v_modes .gt. 0 ) boundary_nodes(2)%vib(:,i,j,k) = maxwell_value * beta3 * vib_df

                end do
             end do
          end do

          boundary_nodes(2)%value = ( flux_in / flux_out ) * boundary_nodes(2)%value
          if( r_modes .gt. 0 ) boundary_nodes(2)%rot = ( flux_in / flux_out ) * boundary_nodes(2)%rot
          if( v_modes .gt. 0 ) boundary_nodes(2)%vib = ( flux_in / flux_out ) * boundary_nodes(2)%vib

       end if

       ! SPECULAR REFLECTION
       !---------------------------------------------------------------------------------------------------------
       if( acc_coeff .lt. one )then
          do k = k_min, k_max
             z = vel_grid%z(k)
   
             do j = j_min, j_max
                y = vel_grid%y(j)
   
                do i = vel_grid%i_pos, vel_grid%i_max
                   x = -vel_grid%x(i)
                   x_cfl = x_factor * x
   
                   dens = ( one - acc_coeff ) * phi(nx,ny)%value(i,j,k)
   
                   if( r_modes .gt. 0 ) rot_df = ( one - acc_coeff ) * phi(nx,ny)%rot(:,i,j,k)
                   if( v_modes .gt. 0 ) vib_df = ( one - acc_coeff ) * phi(nx,ny)%vib(:,i,j,k)
   
                   call perform_remapping( mapping, x, y, z, one, vel_grid )
   
                   call apply_bc_remapping( boundary_nodes(2), mapping, dens )
                   if( r_modes .gt. 0 ) call apply_bc_rot_remapping( boundary_nodes(2), mapping, rot_df )
                   if( v_modes .gt. 0 ) call apply_bc_vib_remapping( boundary_nodes(2), mapping, vib_df )
   
                end do
             end do
          end do
       end if
       ! Test
!       print *, "Right wall boundary node"
!       do k = k_min, k_max
!          do j = j_min, j_max
!             do i = vel_grid%i_pos, vel_grid%i_max
!                ! Set mirror index of i assuming i_min=0
!                i_ref = vel_grid%i_max-i
!                print *, i, j, k, phi(nx,ny)%value(i,j,k) - boundary_nodes(2)%value(i_ref,j,k)
!             end do
!          end do
!       end do
    else
!       print *, "Right boundary node gets value of node to the right (not at right wall)"
       ! Not at right wall so boundary node gets value of cell to the left
       boundary_nodes(2)%value=phi(nx+1,ny)%value
       if( i_zero_flag )then
          boundary_nodes(2)%value(i_zero,:,:) = 0
       end if
       ! Test
!       print *, "RIGHT"
!       do k = k_min, k_max
!          do j = j_min, j_max
!             do i = vel_grid%i_min, vel_grid%i_neg
!                print *, i, j, k, phi(nx+1,ny)%value(i,j,k) - boundary_nodes(2)%value(i,j,k)
!             end do
!             do i = vel_grid%i_pos, vel_grid%i_max
!                print *, i, j, k, phi(nx+1,ny)%value(i,j,k) - boundary_nodes(2)%value(i,j,k)
!             end do
!          end do
!       end do
    end if

    ! Bottom wall
    !=========================================================================================================
    if( ny .eq. 1 )then
!       print *, "Bottom boundary node gets wall props"
       ! Get wall properties
!       call get_left_wall_temp( kin_temp, species )
!       call get_left_wall_velocity( u, v, w, species )
       call get_bottom_wall_temp( kin_temp, species )
       call get_bottom_wall_velocity( u, v, w, species )
       ! Calculate properties into the wall
       dens_in = zero
       mom_in  = zero
       flux_in = zero
       do k = k_min, k_max
          do i = i_min, i_max
             do j = vel_grid%j_min, vel_grid%j_neg
                y = vel_grid%y(j)
                dens_in = dens_in + phi(nx,ny)%value(i,j,k)
                mom_in  = mom_in + y * phi(nx,ny)%value(i,j,k)
                flux_in = flux_in + abs(y) * acc_coeff * phi(nx,ny)%value(i,j,k) 
             end do
          end do
       end do

       ! DIFFUSE REFLECTION
       !---------------------------------------------------------------------------------------------------------
       if( acc_coeff .gt. zero .or. evaporation .ne. zero )then
          ! Rotational and vibrational temperatures are set to the wall temperature
          if( r_modes .gt. 0 )then
             r_levels = phi(nx,ny)%num_rot_levels
             allocate( rot_levels(1:r_levels), rot_df(1:r_levels) )
             call compute_rot_distribution( rot_df, rot_levels, molecule, kin_temp, r_levels, species )
          end if
          
          if( v_modes .gt. 0 )then
             v_levels = phi(nx,ny)%num_vib_levels
             allocate( vib_levels(1:v_levels), vib_df(1:v_levels) )
             call compute_vib_distribution( vib_df, vib_levels, molecule, kin_temp, v_levels, species )
          end if
          
          ! Diffuse Reflection
          dens = acc_coeff * dens_in
          maxwell_coeff = sqrt( mass * mass * mass / ( pi * pi * pi * kin_temp * kin_temp * kin_temp ) )
          computed_dens = zero
          flux_out = zero

          do k = k_min, k_max
             z = vel_grid%z(k)
             beta_z = vel_grid%beta_z(k)
             
             do i = i_min, i_max
                x = vel_grid%x(i)
                beta_x = vel_grid%beta_x(i)

                do j = vel_grid%j_pos, vel_grid%j_max
                   y = vel_grid%y(j)
                   beta_y = vel_grid%beta_y(j)
                   beta3 = beta_x*beta_y*beta_z

                   ! TODO: should only have to do this once
                   call compute_maxwellian( x, y, z, maxwell_coeff, mass, u, v, w, kin_temp, maxwell_value )

                   boundary_nodes(3)%value(i,j,k) = maxwell_value * beta3

                   flux_out = flux_out + abs(y) * maxwell_value * beta3
                   computed_dens = computed_dens + boundary_nodes(3)%value(i,j,k)

                   if( r_modes .gt. 0 ) boundary_nodes(3)%rot(:,i,j,k) = maxwell_value * beta3 * rot_df
                   if( v_modes .gt. 0 ) boundary_nodes(3)%vib(:,i,j,k) = maxwell_value * beta3 * vib_df

                end do
             end do
          end do

          boundary_nodes(3)%value = ( flux_in / flux_out ) * boundary_nodes(3)%value
          if( r_modes .gt. 0 ) boundary_nodes(3)%rot = ( flux_in / flux_out ) * boundary_nodes(3)%rot
          if( v_modes .gt. 0 ) boundary_nodes(3)%vib = ( flux_in / flux_out ) * boundary_nodes(3)%vib

       end if

       ! SPECULAR REFLECTION
       !---------------------------------------------------------------------------------------------------------
       if( acc_coeff .lt. one )then
          do k = k_min, k_max
             z = vel_grid%z(k)
   
             do j = vel_grid%j_min, vel_grid%j_neg
                y = -vel_grid%y(j)
   
                do i = i_min, i_max
                   x = vel_grid%x(i)
   
                   dens = ( one - acc_coeff ) * phi(nx,ny)%value(i,j,k)
   
                   if( r_modes .gt. 0 ) rot_df = ( one - acc_coeff ) * phi(nx,ny)%rot(:,i,j,k)
                   if( v_modes .gt. 0 ) vib_df = ( one - acc_coeff ) * phi(nx,ny)%vib(:,i,j,k)
   
                   call perform_remapping( mapping, x, y, z, one, vel_grid )
   
                   call apply_bc_remapping( boundary_nodes(3), mapping, dens )
                   if( r_modes .gt. 0 ) call apply_bc_rot_remapping( boundary_nodes(3), mapping, rot_df )
                   if( v_modes .gt. 0 ) call apply_bc_vib_remapping( boundary_nodes(3), mapping, vib_df )
   
                end do
             end do
          end do
       end if
       ! Test
!       print *, "Bottom wall boundary node"
!       do k = k_min, k_max
!          do i = i_min, i_max
!             do j = vel_grid%j_min, vel_grid%j_neg
!                ! Set mirror index of j assuming j_min=0
!                j_ref = vel_grid%j_max-j
!                print *, i, j, k, phi(nx,ny)%value(i,j,k) - boundary_nodes(3)%value(i,j_ref,k)
!             end do
!          end do
!       end do
    else
!       print *, "Bottom boundary node gets value of node below (not at bottom wall)"
       ! Not at bottom wall so boundary node gets value of cell below
       boundary_nodes(3)%value=phi(nx,ny-1)%value
       if( j_zero_flag )then
          boundary_nodes(3)%value(:,j_zero,:) = 0
       end if
       ! Test
!       print *, "BOTTOM"
!       do k = k_min, k_max
!          do i = i_min, i_max
!             do j = vel_grid%j_min, vel_grid%j_neg
!                print *, i, j, k, phi(nx,ny-1)%value(i,j,k) - boundary_nodes(3)%value(i,j,k)
!             end do
!             do j = vel_grid%j_pos, vel_grid%j_max
!                print *, i, j, k, phi(nx,ny-1)%value(i,j,k) - boundary_nodes(3)%value(i,j,k)
!             end do
!          end do
!       end do
    end if


    ! Top wall
    if( ny .eq. ny_space )then
!       print *, "Top boundary node gets wall props"
       ! Get wall properties
!       call get_right_wall_temp( kin_temp, species )
!       call get_right_wall_velocity( u, v, w, species )
       call get_top_wall_temp( kin_temp, species )
       call get_top_wall_velocity( u, v, w, species )
       ! Calculate properties into the wall
       dens_in = zero
       mom_in  = zero
       flux_in = zero
       do k = k_min, k_max
          do i = i_min, i_max
             do j = vel_grid%j_pos, vel_grid%j_max
                y = vel_grid%y(j)
                dens_in = dens_in + phi(nx,ny)%value(i,j,k)
                mom_in  = mom_in + y * phi(nx,ny)%value(i,j,k)
                flux_in = flux_in + abs(y) * acc_coeff * phi(nx,ny)%value(i,j,k) 
             end do
          end do
       end do

       ! DIFFUSE REFLECTION
       !---------------------------------------------------------------------------------------------------------
       if( acc_coeff .gt. zero .or. evaporation .ne. zero )then
          ! Rotational and vibrational temperatures are set to the wall temperature
          if( r_modes .gt. 0 )then
             r_levels = phi(nx,ny)%num_rot_levels
             allocate( rot_levels(1:r_levels), rot_df(1:r_levels) )
             call compute_rot_distribution( rot_df, rot_levels, molecule, kin_temp, r_levels, species )
          end if
          
          if( v_modes .gt. 0 )then
             v_levels = phi(nx,ny)%num_vib_levels
             allocate( vib_levels(1:v_levels), vib_df(1:v_levels) )
             call compute_vib_distribution( vib_df, vib_levels, molecule, kin_temp, v_levels, species )
          end if
          
          ! Diffuse Reflection
          dens = acc_coeff * dens_in
          maxwell_coeff = sqrt( mass * mass * mass / ( pi * pi * pi * kin_temp * kin_temp * kin_temp ) )
          computed_dens = zero
          flux_out = zero

          do k = k_min, k_max
             z = vel_grid%z(k)
             beta_z = vel_grid%beta_z(k)
             
             do i = i_min, i_max
                x = vel_grid%x(i)
                beta_x = vel_grid%beta_x(i)

                do j = vel_grid%j_min, vel_grid%j_neg
                   y = vel_grid%y(j)
                   beta_y = vel_grid%beta_y(j)
                   beta3 = beta_x*beta_y*beta_z

                   ! TODO: should only have to do this once
                   call compute_maxwellian( x, y, z, maxwell_coeff, mass, u, v, w, kin_temp, maxwell_value )

                   boundary_nodes(4)%value(i,j,k) = maxwell_value * beta3

                   flux_out = flux_out + abs(y) * maxwell_value * beta3
                   computed_dens = computed_dens + boundary_nodes(4)%value(i,j,k)

                   if( r_modes .gt. 0 ) boundary_nodes(4)%rot(:,i,j,k) = maxwell_value * beta3 * rot_df
                   if( v_modes .gt. 0 ) boundary_nodes(4)%vib(:,i,j,k) = maxwell_value * beta3 * vib_df

                end do
             end do
          end do

          boundary_nodes(4)%value = ( flux_in / flux_out ) * boundary_nodes(4)%value
          if( r_modes .gt. 0 ) boundary_nodes(4)%rot = ( flux_in / flux_out ) * boundary_nodes(4)%rot
          if( v_modes .gt. 0 ) boundary_nodes(4)%vib = ( flux_in / flux_out ) * boundary_nodes(4)%vib

       end if

       ! SPECULAR REFLECTION
       !---------------------------------------------------------------------------------------------------------
       if( acc_coeff .lt. one )then
          do k = k_min, k_max
             z = vel_grid%z(k)
   
             do j = vel_grid%j_pos, vel_grid%j_max
                y = -vel_grid%y(j)
   
                do i = i_min, i_max
                   x = vel_grid%x(i)
   
                   dens = ( one - acc_coeff ) * phi(nx,ny)%value(i,j,k)
   
                   if( r_modes .gt. 0 ) rot_df = ( one - acc_coeff ) * phi(nx,ny)%rot(:,i,j,k)
                   if( v_modes .gt. 0 ) vib_df = ( one - acc_coeff ) * phi(nx,ny)%vib(:,i,j,k)
                   call perform_remapping( mapping, x, y, z, one, vel_grid )
                   
                   call apply_bc_remapping( boundary_nodes(4), mapping, dens )
                   if( r_modes .gt. 0 ) call apply_bc_rot_remapping( boundary_nodes(4), mapping, rot_df )
                   if( v_modes .gt. 0 ) call apply_bc_vib_remapping( boundary_nodes(4), mapping, vib_df )
   
                end do
             end do
          end do
       end if
       ! Test
!       print *, "Top wall boundary node"
!       do k = k_min, k_max
!          do i = i_min, i_max
!             do j = vel_grid%j_pos, vel_grid%j_max
!                ! Set mirror index of i assuming i_min=0
!                j_ref = vel_grid%j_max-j
!                print *, i, j, k, phi(nx,ny)%value(i,j,k) - boundary_nodes(4)%value(i,j,k)
!             end do
!          end do
!       end do
    else
!       print *, "Top boundary node gets value of node above (not at top wall)"
       ! Not at top wall so boundary node gets value of cell above
       boundary_nodes(4)%value=phi(nx,ny+1)%value
       if( j_zero_flag )then
          boundary_nodes(4)%value(:,j_zero,:) = 0
       end if
       ! Test
!       print *, "TOP"
!       do k = k_min, k_max
!          do i = i_min, i_max
!             do j = vel_grid%j_min, vel_grid%j_neg
!                print *, i, j, k, phi(nx,ny+1)%value(i,j,k) - boundary_nodes(4)%value(i,j,k)
!             end do
!             do j = vel_grid%j_pos, vel_grid%j_max
!                print *, i, j, k, phi(nx,ny+1)%value(i,j,k) - boundary_nodes(4)%value(i,j,k)
!             end do
!          end do
!       end do
    end if

!    ! Set Phi and Calculate properties off the wall
!    dens_out = zero
!    mom_out  = zero
!    do k = k_min, k_max
!       do j = j_min, j_max
!          do i = out_min, out_max
!             x = vel_grid%x(i)
!             dens_out = dens_out + boundary_nodes(wall_index)%value(i,j,k)
!             mom_out  = mom_out + x * boundary_nodes(wall_index)%value(i,j,k)
!          end do
!       end do
!    end do
!
!    ! Print information to file
!    force_on_wall = mass * ( mom_in - mom_out ) / dt
!
!    if( nx .eq. 1 )then
!       open( unit = 10, file = "left_wall.dat", status="old", position="append" )
!       write(10,fmt="(5e24.14)") dens_in, dens_out, mom_in, mom_out, force_on_wall
!       close( 10 )
!    end if
!    if( nx .eq. nx_space )then
!       open( unit = 11, file = "right_wall.dat", status="old", position="append" )
!       write(11,fmt="(5e24.14)") dens_in, dens_out, mom_in, mom_out, force_on_wall
!       close( 11 )
!    end if

    ! First-Order Upwind FD Effects
    !TODO: Check how interpolation to other velocity nodes works
    !TODO: ADD INTERNAL ENERGY
    !--------------------------------------------------------------------------------------------------------

    do k = k_min, k_max
       z = vel_grid%z(k)

       do j = vel_grid%j_min, vel_grid%j_neg
          y = vel_grid%y(j)
          y_cfl = y_factor * y

          do i = vel_grid%i_min, vel_grid%i_neg
             x = vel_grid%x(i)
             x_cfl = x_factor * x

             phi_conv(nx,ny)%value(i,j,k) = phi(nx,ny)%value(i,j,k) &
                  + abs(x_cfl) * ( boundary_nodes(2)%value(i,j,k) - phi(nx,ny)%value(i,j,k) ) &
                  + abs(y_cfl) * ( boundary_nodes(4)%value(i,j,k) - phi(nx,ny)%value(i,j,k) )

          end do

          if( i_zero_flag )then
             phi_conv(nx,ny)%value(i_zero,j,k) = phi(nx,ny)%value(i_zero,j,k) &
                  + boundary_nodes(1)%value(i_zero,j,k) + boundary_nodes(2)%value(i_zero,j,k) &
                  + abs(y_cfl) * ( boundary_nodes(4)%value(i_zero,j,k) - phi(nx,ny)%value(i_zero,j,k) )
          end if

          do i = vel_grid%i_pos, vel_grid%i_max
             x = vel_grid%x(i)
             x_cfl = x_factor * x

             phi_conv(nx,ny)%value(i,j,k) = phi(nx,ny)%value(i,j,k) &
                  + abs(x_cfl) * ( boundary_nodes(1)%value(i,j,k) - phi(nx,ny)%value(i,j,k) ) &
                  + abs(y_cfl) * ( boundary_nodes(4)%value(i,j,k) - phi(nx,ny)%value(i,j,k) )

          end do

       end do

       if( j_zero_flag )then
          do i = vel_grid%i_min, vel_grid%i_neg
             x = vel_grid%x(i)
             x_cfl = x_factor * x

             phi_conv(nx,ny)%value(i,j_zero,k) = phi(nx,ny)%value(i,j_zero,k) &
                  + abs(x_cfl) * ( boundary_nodes(2)%value(i,j_zero,k) - phi(nx,ny)%value(i,j_zero,k) ) &
                  +  boundary_nodes(3)%value(i,j_zero,k) + boundary_nodes(4)%value(i,j_zero,k)

          end do

          if( i_zero_flag )then
             phi_conv(nx,ny)%value(i_zero,j_zero,k) = phi(nx,ny)%value(i_zero,j_zero,k) &
                  + boundary_nodes(1)%value(i_zero,j_zero,k) + boundary_nodes(2)%value(i_zero,j_zero,k) &
                  + boundary_nodes(3)%value(i_zero,j_zero,k) + boundary_nodes(4)%value(i_zero,j_zero,k)

          end if

          do i = vel_grid%i_pos, vel_grid%i_max
             x = vel_grid%x(i)
             x_cfl = x_factor * x

             phi_conv(nx,ny)%value(i,j_zero,k) = phi(nx,ny)%value(i,j_zero,k) &
                  + abs(x_cfl) * ( boundary_nodes(1)%value(i,j_zero,k) - phi(nx,ny)%value(i,j_zero,k) ) &
                  + boundary_nodes(3)%value(i,j_zero,k) + boundary_nodes(4)%value(i,j_zero,k)

          end do
       
       end if

       do j = vel_grid%j_pos, vel_grid%j_max
          y = vel_grid%y(j)
          y_cfl = y_factor * y

          do i = vel_grid%i_min, vel_grid%i_neg
             x = vel_grid%x(i)
             x_cfl = x_factor * x

             phi_conv(nx,ny)%value(i,j,k) = phi(nx,ny)%value(i,j,k) &
                  + abs(x_cfl) * ( boundary_nodes(2)%value(i,j,k) - phi(nx,ny)%value(i,j,k) ) &
                  + abs(y_cfl) * ( boundary_nodes(3)%value(i,j,k) - phi(nx,ny)%value(i,j,k) )

          end do

          if( i_zero_flag )then
             phi_conv(nx,ny)%value(i_zero,j,k) = phi(nx,ny)%value(i_zero,j,k) &
                  + boundary_nodes(1)%value(i_zero,j_zero,k) + boundary_nodes(2)%value(i_zero,j_zero,k) &
                  + abs(y_cfl) * ( boundary_nodes(3)%value(i_zero,j,k) - phi(nx,ny)%value(i_zero,j,k) )
          end if

          do i = vel_grid%i_pos, vel_grid%i_max
             x = vel_grid%x(i)
             x_cfl = x_factor * x

             phi_conv(nx,ny)%value(i,j,k) = phi(nx,ny)%value(i,j,k) &
                  + abs(x_cfl) * ( boundary_nodes(1)%value(i,j,k) - phi(nx,ny)%value(i,j,k) ) &
                  + abs(y_cfl) * ( boundary_nodes(3)%value(i,j,k) - phi(nx,ny)%value(i,j,k) )

          end do
       
       end do

    end do

    ! Deallocate energy arrays
    if( r_modes .gt. 0 ) deallocate( rot_levels, rot_df )
    if( v_modes .gt. 0 ) deallocate( vib_levels, vib_df )

    ! Deallocate boundary nodes
    call destroy_dist_func( boundary_nodes(1), molecule )
    call destroy_dist_func( boundary_nodes(2), molecule )
    call destroy_dist_func( boundary_nodes(3), molecule )
    call destroy_dist_func( boundary_nodes(4), molecule )

    return
  end subroutine apply_wall_reflection
          
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
    double precision :: computed_dens, maxwell_value, maxwell_coeff

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
    integer :: i, j, k
    
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
    computed_dens = zero
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
             computed_dens = computed_dens + maxwell_value * beta3

             if( r_modes .gt. 0 ) phi_conv(nx)%rot(:,i,j,k) = maxwell_value*beta3*rot_df
             if( v_modes .gt. 0 ) phi_conv(nx)%vib(:,i,j,k) = maxwell_value*beta3*vib_df

          end do
       end do
    end do

    phi_conv(nx)%value = ( dens / computed_dens ) * phi_conv(nx)%value

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
