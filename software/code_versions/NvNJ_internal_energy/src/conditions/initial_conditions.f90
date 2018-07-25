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
  integer, dimension(:), allocatable :: init_domain, shock_thickness
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

  public :: set_initial_df_conditions
  public :: set_delta_function_location
  public :: set_spatial_domain_conditions
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

  subroutine set_spatial_domain_conditions( init_domain_in, shock_thickness_in )

    use DistFunc
    
    implicit none

    integer, dimension(:) :: init_domain_in, shock_thickness_in

    integer :: status

    allocate( init_domain(1:num_species), STAT=status )
    call allocate_error_check( status, "init_domain" )

    allocate( shock_thickness(1:num_species), STAT=status )
    call allocate_error_check( status, "shock_thickness" )

    init_domain     = init_domain_in
    shock_thickness = shock_thickness_in

    return
  end subroutine set_spatial_domain_conditions

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

    double precision :: dens_LW, u_LW, temp_LW
    double precision :: dens_RW, u_RW, temp_RW
    double precision :: mass

    double precision :: dens_local, u_local, temp_local
    double precision :: dens_old, u_old
    double precision :: dens_ratio, maxwell_coeff

    double precision :: temp_rot_local, temp_vib_local
    double precision :: temp_rot_LW, temp_rot_RW, temp_vib_LW, temp_vib_RW
    double precision :: theta_r, theta_v
    double precision :: vib_part_func

    double precision :: x, y, z
    double precision :: beta_x, beta_y, beta_z, beta

    double precision :: df_value, energy_level, fraction

    double precision, dimension(:,:), allocatable :: int_energy_df, int_energy_levels

    double precision :: even_spin, odd_spin

    integer :: ns, grid_ref, dens_midpoint, temp_midpoint
    integer :: shock_thick, half_shock
    integer :: i, j, k, i_min, i_max, j_min, j_max, k_min, k_max
    integer :: nspace

    integer :: mol_shape
    integer :: r_modes, v_modes
    integer :: l, r_levels, v_levels

    logical :: properties_reset

!!$!>>>>>>>>>>>>>>>>>>>
    double precision :: temp_vib_temp
!!$!>>>>>>>>>>>>>>>>>>>

    properties_reset = .false.

    call get_nspace( nspace )

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
    call get_left_wall_velocity( u_LW, species )
    call get_right_wall_velocity( u_RW, species )
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

    if( init_domain(species) .eq. rankine_hugoniot )then
       shock_thick = shock_thickness(species)

       dens_LW    = shock_props%dens_ratio*dens_RW
       u_LW       = shock_props%u_down
       temp_LW    = shock_props%T_ratio*temp_RW

       ! Update boundary properties to shock properties
       call set_boundary_properties( dens_LW, u_LW, temp_LW, dens_RW, u_RW, temp_RW, species )
       
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
       temp_local = temp_LW

       if( r_modes .gt. 0 ) temp_rot_local = temp_rot_LW
       if( v_modes .gt. 0 ) temp_vib_local = temp_vib_LW

    case( right_wall )
       dens_local = dens_RW
       u_local    = u_RW
       temp_local = temp_RW

       if( r_modes .gt. 0 ) temp_rot_local = temp_rot_RW
       if( v_modes .gt. 0 ) temp_vib_local = temp_vib_RW

    case default
       write(*,*) "Error: Invalid value of init_domain. Given value is: ", init_domain(species)
       stop
    end select

    half_shock = floor( one_half*shock_thick )
    dens_midpoint = floor( one_half*nspace )
    temp_midpoint = floor( one_half*nspace + one_half*half_shock )
    
    do ns = 1, nspace

       if( r_modes .gt. 0 .and. v_modes .gt. 0 )then
          r_levels = phi(ns)%num_rot_levels
          v_levels = phi(ns)%num_vib_levels
          allocate( int_energy_df(1:r_levels,1:v_levels) )
          allocate( int_energy_levels(1:r_levels,1:v_levels) )
       end if

       ! If steady shock, determine macroscopic properties
       if( init_domain(species) .eq. rankine_hugoniot )then
          if( properties_reset .eqv. .false. )then ! Haven't reached right wall region yet

             dens_old = dens_local
             u_old    = u_local

             if( ns .ge. dens_midpoint + half_shock )then
                dens_local = dens_RW
                u_local    = u_RW

             else if( ns .gt. dens_midpoint )then
                call initial_shock_profile( dens_local, dens_LW, dens_RW, half_shock, dens_midpoint, ns )
                dens_ratio = ( dens_local - dens_RW )/( dens_LW - dens_RW )
                u_local = u_RW + ( u_LW - u_RW )*dens_ratio

                write(*,*) dens_local, dens_old, dens_ratio, u_old, u_local

             else if( ns .ge. dens_midpoint - half_shock )then
                call initial_shock_profile( dens_local, dens_RW, dens_LW, half_shock, ns, dens_midpoint )
                dens_ratio = ( dens_local - dens_RW )/( dens_LW - dens_RW )
                u_local = u_RW + ( u_LW - u_RW )*dens_ratio

             else
                dens_local = dens_LW
                u_local    = u_LW

             end if

             if( ns .ge. temp_midpoint + half_shock )then
                temp_local = temp_RW
                properties_reset = .true.

                if( r_modes .gt. 0 ) temp_rot_local = temp_rot_RW
                if( v_modes .gt. 0 ) temp_vib_local = temp_vib_RW

             else if( ns .gt. temp_midpoint )then
                call initial_shock_profile( temp_local, temp_LW, temp_RW, half_shock, temp_midpoint, ns )
                
                if( r_modes .gt. 0 ) temp_rot_local = temp_local
                if( v_modes .gt. 0 ) temp_vib_local = temp_local

             else if( ns .ge. temp_midpoint - half_shock )then
                call initial_shock_profile( temp_local, temp_RW, temp_LW, half_shock, ns, temp_midpoint )

                if( r_modes .gt. 0 ) temp_rot_local = temp_local
                if( v_modes .gt. 0 ) temp_vib_local = temp_local

             else
                temp_local = temp_LW

                if( r_modes .gt. 0 ) temp_rot_local = temp_local
                if( v_modes .gt. 0 ) temp_rot_local = temp_local

             end if

          end if
       end if

       maxwell_coeff = dens_local*sqrt( mass*mass*mass )/sqrt( pi*pi*pi*temp_local*temp_local*temp_local )
       
       if( r_modes .gt. 0 .and. v_modes .gt. 0 )then
          call compute_rot_vib_distribution( int_energy_df, int_energy_levels, molecule, &
               temp_rot_local, temp_vib_local, r_levels, v_levels, species )
          phi(ns)%int_energy_level = int_energy_levels

          call compute_rot_distribution( phi(ns)%rot(:,1,1,1), phi(ns)%rot_level, molecule, &
               temp_rot_local, r_levels, species )
          call compute_vib_distribution( phi(ns)%vib(:,1,1,1), phi(ns)%vib_level, molecule, &
               temp_vib_local, v_levels, species )

       end if

       call get_spatial_reference( ns, grid_ref )
       i_min = vel_grid(grid_ref)%i_min
       i_max = vel_grid(grid_ref)%i_max
       j_min = vel_grid(grid_ref)%j_min
       j_max = vel_grid(grid_ref)%j_max
       k_min = vel_grid(grid_ref)%k_min
       k_max = vel_grid(grid_ref)%k_max

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
                   ! hardcoded v, w = 0
                   call compute_maxwellian( x, y, z, maxwell_coeff, mass, &
                        u_local, zero, zero, temp_local, df_value )

                case( bkw )
                   ! hardcoded time, v, w = 0, time = 0
                   call compute_bkw( x, y, z, dens_local, temp_local, &
                        mass, u_local, zero, df_value )

                case( delta_function )
                   ! only one delta function per species allowed
                   if( abs( x_delta_loc(species) - x ) .lt. double_tol .and. &
                       abs( y_delta_loc(species) - y ) .lt. double_tol .and. &
                       abs( z_delta_loc(species) - z ) .lt. double_tol         )then
                      df_value = dens_local/beta
                   else
                      df_value = zero
                   end if

                case default
                   write(*,*) "Error: Invalid value of init_dist_func_flag. Value is: ", init_vel_df(species)
                   stop
                   
                end select

                phi(ns)%value(i,j,k) = df_value*beta

                if( r_modes .gt. 0 .and. v_modes .gt. 0 )then
                   phi(ns)%int_energy(:,:,i,j,k) = df_value*beta*int_energy_df
                   do l = 1, r_levels
                      phi(ns)%rot(l,i,j,k) = sum( phi(ns)%int_energy(l,1:v_levels,i,j,k) )
                   end do
                   do l = 1, v_levels
                      phi(ns)%vib(l,i,j,k) = sum( phi(ns)%int_energy(1:r_levels,l,i,j,k) )
                   end do
                end if

             end do
          end do
       end do

    end do

    return
  end subroutine initialize_dist_func

  subroutine initial_shock_profile( local, LW, RW, thick, midpoint, ns )

    implicit none

    double precision :: local

    double precision, intent(in) :: LW, RW
    integer, intent(in) :: thick, midpoint, ns

    double precision :: exponent, coeff, local_x

    exponent = one_half*( LW - RW )/dble( thick )/dble( thick )
    coeff    = RW
    local_x  = dble( midpoint - ns + thick )

    local = coeff + exponent*local_x*local_x

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
    double precision :: x, y, z, beta_x, beta_y, beta_z, beta3
    
    integer :: i, j, k
    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

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

             call compute_bkw( x, y, z, dens, temp, mass, x_vel, time, df_value )

             bkw%value(i,j,k) = df_value*beta3

          end do
       end do
    end do

    return
  end subroutine generate_bkw

end module InitialConditions

