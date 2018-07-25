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
module InitialConditions

  use Constants
  use ErrorCheck

  implicit none

  private

  integer, dimension(:), allocatable :: init_vel_df, init_rot_df, init_vib_df
  integer, dimension(:), allocatable :: init_domain, shock_thickness, shock_middle
  double precision, dimension(:), allocatable :: x_delta_loc, y_delta_loc, z_delta_loc

  integer, parameter :: left_wall        = 1
  integer, parameter :: right_wall       = 2
  integer, parameter :: left_right_1D    = 3
  integer, parameter :: rankine_hugoniot = 4

  integer, parameter :: maxwell        = 1
  integer, parameter :: bkw            = 2
  integer, parameter :: delta_function = 3

  integer, parameter :: rigid_rotor    = 1
  integer, parameter :: nonrigid_rotor = 2

  integer, parameter :: sho = 1
  integer, parameter :: aho = 2

  double precision :: calculated_dens, dens_ratio

  public :: set_initial_df_conditions
  public :: set_delta_function_location
  public :: set_spatial_domain_conditions
  public :: get_spatial_domain_conditions
  public :: initialize_dist_func
  public :: generate_bkw

contains

  subroutine set_initial_df_conditions( init_vel_df_in )!, init_rot_df_in, init_vib_df_in )

    use DistFunc

    implicit none

    integer, dimension(:), intent(in) :: init_vel_df_in!, init_rot_df_in, init_vib_df_in

    integer :: status

    allocate( init_vel_df(1:num_species), STAT=status )
    call allocate_error_check( status, "init_vel_df" )

!!$    allocate( init_rot_df(1:num_species), STAT=status )
!!$    call allocate_error_check( status, "init_rot_df" )
!!$    
!!$    allocate( init_vib_df(1:num_species), STAT=status )
!!$    call allocate_error_check( status, "init_vib_df" )

    init_vel_df = init_vel_df_in
!!$    init_rot_df = init_rot_df_in
!!$    init_vib_df = init_vib_df_in

    return
  end subroutine set_initial_df_conditions

  subroutine set_delta_function_location( x_delta_loc_in, y_delta_loc_in, z_delta_loc_in )

    use DistFunc

    implicit none

    double precision, dimension(:), intent(in) :: x_delta_loc_in, y_delta_loc_in, z_delta_loc_in

    integer :: status

    allocate( x_delta_loc(1:num_species), STAT=status )
    call allocate_error_check( status, "x_delta_loc" )

    allocate( y_delta_loc(1:num_species), STAT=status )
    call allocate_error_check( status, "y_delta_loc" )

    allocate( z_delta_loc(1:num_species), STAT=status )
    call allocate_error_check( status, "z_delta_loc" )

    x_delta_loc = x_delta_loc_in
    y_delta_loc = y_delta_loc_in
    z_delta_loc = z_delta_loc_in

    return
  end subroutine set_delta_function_location

  subroutine set_spatial_domain_conditions( init_domain_in, shock_thickness_in, shock_middle_in )

    use DistFunc

    implicit none

    integer, dimension(:) :: init_domain_in, shock_thickness_in, shock_middle_in

    integer :: status

    allocate( init_domain(1:num_species), STAT=status )
    call allocate_error_check( status, "init_domain" )

    allocate( shock_thickness(1:num_species), STAT=status )
    call allocate_error_check( status, "shock_thickness" )

    allocate( shock_middle(1:num_species), STAT=status )
    call allocate_error_check( status, "shock_middle" )

    init_domain     = init_domain_in
    shock_thickness = shock_thickness_in
    shock_middle    = shock_middle_in

    return
  end subroutine set_spatial_domain_conditions

  subroutine get_spatial_domain_conditions( init_domain_out, species )
    
    implicit none

    integer, intent(in) :: species

    integer :: init_domain_out

    init_domain_out = init_domain( species )

    return
  end subroutine get_spatial_domain_conditions

  subroutine initialize_dist_func( phi, vel_grid, molecule, shock_props, species )

    use DistFunc
    use VelocityGrid
    use PhysicalGrid
    use ShockConditions
    use BoundaryConditions
    use SpeciesAndReferenceData

    implicit none

    type(DistFuncType), dimension(:) :: phi
    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    type(MoleculeType), intent(in) :: molecule
    type(NormalShock), intent(in) :: shock_props
    integer, intent(in) :: species

    double precision :: dens_LW, u_LW, v_LW, w_LW, temp_LW
    double precision :: dens_RW, u_RW, v_RW, w_RW, temp_RW
    double precision :: mass

    double precision :: dens_local, temp_local
    double precision :: u_local, v_local, w_local
    double precision :: dens_old, u_old
    double precision :: dens_ratio
    double precision :: opt, maxwell_coeff, bkw_xk, bkw_norm

    double precision :: temp_rot_local, temp_vib_local
    double precision :: temp_rot_LW, temp_rot_RW, temp_vib_LW, temp_vib_RW
    double precision :: theta_r, theta_v
    double precision :: vib_part_func

    double precision :: x, y, z
    double precision :: beta_x, beta_y, beta_z, beta

    double precision :: df_value, energy_level, fraction
    double precision, dimension(:), allocatable :: rot_df, rot_levels, vib_df, vib_levels

    double precision :: even_spin, odd_spin

    integer :: nx, grid_ref, dens_midpoint, temp_midpoint
    integer :: shock_thick, half_shock, shock_mid
    integer :: i, j, k, i_min, i_max, j_min, j_max, k_min, k_max
    integer :: nx_space, ny_space

    integer :: mol_shape
    integer :: r_modes, v_modes
    integer :: l, r_levels, v_levels

    logical :: properties_reset

    properties_reset = .false.

    call get_nspace( nx_space, ny_space )

    ! Read species constants
    r_modes = molecule%rot_modes
    v_modes = molecule%vib_modes

    mass    = molecule%mass
    theta_r = molecule%theta_r
    theta_v = molecule%theta_v

    mol_shape = molecule%molecule_type
    even_spin = molecule%even_spin
    odd_spin  = molecule%odd_spin

    ! Get boundary condition properties we need
    call get_left_wall_density( dens_LW, species )
    call get_right_wall_density( dens_RW, species )
    call get_left_wall_velocity( u_LW, v_LW, w_LW, species )
    call get_right_wall_velocity( u_RW, v_RW, w_RW, species )
    call get_left_wall_temp( temp_LW, species )
    call get_right_wall_temp( temp_RW, species )

    if( r_modes .gt. 0 )then
       call get_right_wall_temp_rot( temp_rot_RW, species )
       call get_left_wall_temp_rot( temp_rot_LW, species )
    end if

    if( v_modes .gt. 0 )then
       call get_right_wall_temp_vib( temp_vib_RW, species )
       call get_left_wall_temp_vib( temp_vib_LW, species )
    end if

    ! Set Rankine-Hugoniot conditions at left wall if rankine hugoniot
    ! TODO: hard coded for left wall to be downstream
    if( init_domain(species) .eq. rankine_hugoniot )then
       shock_thick = shock_thickness(species)
       shock_mid   = shock_middle(species)

       dens_LW    = shock_props%dens_ratio*dens_RW
       u_LW       = shock_props%u_down
       temp_LW    = shock_props%T_ratio*temp_RW

       ! Update boundary properties to shock properties
       call set_boundary_properties( dens_LW, u_LW, v_LW, w_LW, temp_LW, &
            dens_RW, u_RW, v_RW, w_RW, temp_RW, species )

       ! Internal energy boundaries
       if( r_modes .gt. 0 )then
          temp_rot_LW = temp_LW
          call set_temp_rot( temp_rot_LW, temp_rot_RW, species )
       end if

       if( v_modes .gt. 0 )then
          temp_vib_LW = temp_LW
          call set_temp_vib( temp_vib_LW, temp_vib_RW, species )
       end if

    end if

    ! Set initial local macroscopic properties
    select case( init_domain(species) )
    case( left_wall, left_right_1D, rankine_hugoniot )
       dens_local = dens_LW
       u_local    = u_LW
       v_local    = v_LW
       w_local    = w_LW
       temp_local = temp_LW

       if( r_modes .gt. 0 ) temp_rot_local = temp_rot_LW
       if( v_modes .gt. 0 ) temp_vib_local = temp_vib_LW

    case( right_wall )
       dens_local = dens_RW
       u_local    = u_RW
       v_local    = v_RW
       w_local    = w_RW
       temp_local = temp_RW

       if( r_modes .gt. 0 ) temp_rot_local = temp_rot_RW
       if( v_modes .gt. 0 ) temp_vib_local = temp_vib_RW

    case default
       write(*,*) "Error: Invalid value of init_domain. Given value is: ", init_domain(species)
       stop
    end select

    ! Set shock midpoint, and temperature change midpoint (leads shock)
    half_shock    = ceiling( one_half * shock_thick )
    dens_midpoint = shock_mid
    temp_midpoint = ceiling( shock_mid + one_half * half_shock )

    ! Loop over every spatial node
!!$    do ny = 1, ny_space
    do nx = 1, nx_space

       ! Internal energy array allocation
       if( r_modes .gt. 0 )then
          r_levels = phi(nx)%num_rot_levels
          allocate( rot_df(1:r_levels) )
          allocate( rot_levels(1:r_levels) )
       end if

       if( v_modes .gt. 0 )then
          v_levels = phi(nx)%num_vib_levels
          allocate( vib_df(1:v_levels) )
          allocate( vib_levels(1:v_levels) )
       end if

       ! If steady shock, determine macroscopic properties
       ! One of several predefined shock profile equations are used for initialization of shock shape
       if( init_domain(species) .eq. rankine_hugoniot )then
          if( properties_reset .eqv. .false. )then ! Haven't reached right wall region yet

             dens_old = dens_local
             u_old    = u_local



             if( nx .gt. dens_midpoint + half_shock )then
                dens_local = dens_RW
                u_local    = u_RW

             else if( nx .gt. dens_midpoint )then
                call initial_shock_profile( dens_local, dens_LW, dens_RW, half_shock, dens_midpoint, nx )
                dens_ratio = ( dens_local - dens_RW )/( dens_LW - dens_RW )
                u_local = u_RW + ( u_LW - u_RW )*dens_ratio

             else if( nx .ge. dens_midpoint - half_shock )then
                call initial_shock_profile( dens_local, dens_RW, dens_LW, half_shock, nx, dens_midpoint )
                dens_ratio = ( dens_local - dens_RW )/( dens_LW - dens_RW )
                u_local = u_RW + ( u_LW - u_RW )*dens_ratio

             else
                dens_local = dens_LW
                u_local    = u_LW

             end if

             if( nx .gt. temp_midpoint + half_shock )then
                temp_local = temp_RW
                properties_reset = .true.

                if( r_modes .gt. 0 ) temp_rot_local = temp_rot_RW
                if( v_modes .gt. 0 ) temp_vib_local = temp_vib_RW


             else if( nx .gt. temp_midpoint )then
                call initial_shock_profile( temp_local, temp_LW, temp_RW, half_shock, temp_midpoint, nx )

                if( r_modes .gt. 0 ) temp_rot_local = temp_local
                if( v_modes .gt. 0 ) temp_vib_local = temp_local

             else if( nx .ge. temp_midpoint - half_shock )then
                call initial_shock_profile( temp_local, temp_RW, temp_LW, half_shock, nx, temp_midpoint )

                if( r_modes .gt. 0 ) temp_rot_local = temp_local
                if( v_modes .gt. 0 ) temp_vib_local = temp_local

             else
                temp_local = temp_LW

                if( r_modes .gt. 0 ) temp_rot_local = temp_local
                if( v_modes .gt. 0 ) temp_rot_local = temp_local

             end if

          end if
       end if

       ! Velocity independent coefficient used for calculating initial Maxwellian and BKW values
       opt           = mass / ( pi * temp_local ) ! FLOP optimization
       maxwell_coeff = dens_local * sqrt( opt * opt * opt )
       bkw_xk        = one - 0.4d0 * exp( -zero / 6.0d0 ) ! hardcode time = 0
       bkw_norm      = one_half * dens_local * sqrt( opt * opt * opt ) * bkw_xk**( -2.5d0 )

       ! Calculate initial rotational and vibrational distributions
       if( r_modes .gt. 0 )then
          call compute_rot_distribution( rot_df, rot_levels, molecule, temp_rot_local, r_levels, species )
          phi(nx)%rot_level = rot_levels
       end if

       if( v_modes .gt. 0 )then
          vib_part_func = one/( one - exp( -theta_v/temp_vib_local ) )
          call compute_vib_distribution( vib_df, vib_levels, molecule, temp_vib_local, v_levels, species )
          phi(nx)%vib_level = vib_levels
       end if

       ! Index limits
       call get_spatial_reference( nx, grid_ref )
       i_min = vel_grid(grid_ref)%i_min
       i_max = vel_grid(grid_ref)%i_max
       j_min = vel_grid(grid_ref)%j_min
       j_max = vel_grid(grid_ref)%j_max
       k_min = vel_grid(grid_ref)%k_min
       k_max = vel_grid(grid_ref)%k_max

       ! Discrete distribution density - used to remove discretization error
       calculated_dens = zero

       do k = k_min, k_max
          z = vel_grid(grid_ref)%z(k)
          beta_z = vel_grid(grid_ref)%beta_z(k)

          do j = j_min, j_max
             y = vel_grid(grid_ref)%y(j)
             beta_y = vel_grid(grid_ref)%beta_y(j)

             do i = i_min,i_max
                x = vel_grid(grid_ref)%x(i)
                beta_x = vel_grid(grid_ref)%beta_x(i)
                beta = beta_x*beta_y*beta_z

                select case( init_vel_df(species) )
                case( maxwell )
                   call compute_maxwellian( x, y, z, maxwell_coeff, mass, &
                        u_local, v_local, w_local, temp_local, df_value )

                case( bkw )
                   call compute_bkw( x, y, z, dens_local, temp_local, &
                        mass, u_local, v_local, w_local, bkw_xk, bkw_norm, df_value )

                case( delta_function )
                   ! only one delta function per species allowed
                   if(  abs( x_delta_loc(species) - x ) .lt. double_tol .and. &
                        abs( y_delta_loc(species) - y ) .lt. double_tol .and. &
                        abs( z_delta_loc(species) - z ) .lt. double_tol         )then
                      df_value = dens_local / beta
                   else
                      df_value = zero
                   end if

                case default
                   write(*,*) "Error: Invalid value of init_dist_func_flag. Value is: ", init_vel_df(species)
                   stop

                end select

                phi(nx)%value(i,j,k) = df_value * beta

                if( r_modes .gt. 0 )then
                   phi(nx)%rot(:,i,j,k) = df_value * beta * rot_df
                end if

                if( v_modes .gt. 0 )then
                   phi(nx)%vib(:,i,j,k) = df_value * beta * vib_df
                end if

                ! Track total density in discretized distribution
                calculated_dens = calculated_dens + phi(nx)%value(i,j,k)

             end do
          end do
       end do

       ! At this point, the discrete distribution density is adjusted to match the input density
       dens_ratio = dens_local / calculated_dens

       phi(nx)%value = dens_ratio * phi(nx)%value
       if( r_modes .gt. 0 ) phi(nx)%rot = dens_ratio * phi(nx)%rot
       if( v_modes .gt. 0 ) phi(nx)%vib = dens_ratio * phi(nx)%vib


       ! Internal energy array deallocation
       if( r_modes .gt. 0 )then
          deallocate( rot_df )
          deallocate( rot_levels )
       end if

       if( v_modes .gt. 0 )then
          deallocate( vib_df )
          deallocate( vib_levels )
       end if

    end do
!!$    end do

    return
  end subroutine initialize_dist_func

  subroutine initial_shock_profile( local, LW, RW, thick, midpoint, nx )

    implicit none

    double precision :: local

    double precision, intent(in) :: LW, RW
    integer, intent(in) :: thick, midpoint, nx

    double precision :: exponent, coeff, local_x

    ! The current equation used for shock profiles is a linear fit
    exponent = one_half * ( LW - RW ) / dble( thick * thick )
    coeff    = RW
    local_x  = dble( midpoint - nx + thick )

    local = coeff + exponent * local_x * local_x

    return
  end subroutine initial_shock_profile

  subroutine generate_bkw( bkw, mass, dens, temp, x_vel, vel_grid, time )

    use DistFunc
    use VelocityGrid

    implicit none

    double precision, intent(in) :: time, dens, temp, x_vel, mass
    type(VelocityGridType), intent(in) :: vel_grid

    type(DistFuncType) :: bkw

    double precision :: df_value
    double precision :: opt, bkw_xk, bkw_norm
    double precision :: x, y, z, beta_x, beta_y, beta_z, beta3
    double precision :: calc_dens

    integer :: i, j, k
    integer :: i_min, i_max, j_min, j_max, k_min, k_max

    opt       = mass / ( pi * temp ) ! FLOP optimization
    bkw_xk    = one - 0.4d0 * exp( -time / 6.0d0 )
    bkw_norm  = one_half * dens * sqrt( opt * opt * opt ) * bkw_xk**( -2.5d0 )

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    calc_dens = zero

    do k = k_min,k_max
       z = vel_grid%z(k)
       beta_z = vel_grid%beta_z(k)

       do j = j_min,j_max
          y = vel_grid%y(j)
          beta_y = vel_grid%beta_y(j)

          do i = i_min,i_max
             x = vel_grid%x(i)
             beta_x = vel_grid%beta_x(i)
             beta3 = beta_x*beta_y*beta_z

             ! hardcoded v=0, w=0
             call compute_bkw( x, y, z, dens, temp, mass, x_vel, zero, zero, bkw_xk, bkw_norm, df_value )

             bkw%value(i,j,k) = df_value * beta3
             calc_dens = calc_dens + df_value * beta3
             
          end do
       end do
    end do

    bkw%value = ( dens / calc_dens ) * bkw%value

    return
  end subroutine generate_bkw

end module InitialConditions

