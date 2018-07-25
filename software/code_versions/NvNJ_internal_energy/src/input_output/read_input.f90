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
module ReadInput

  use DistFunc
  use VelocityGrid
  use PhysicalGrid
  use PhysicalProperties
  use InitialConditions
  use BoundaryConditions
  use CollisionIntegral
  use CollisionUtilities
  use ReplenishingCollisions
  use Convection
  use TimeStepping
  use SpeciesAndReferenceData
  use Scaling
  use RandomNumberGeneration
  use ErrorCheck
  use ScreenOutput
  use Visualization
  use Restart

  implicit none

  private

  public :: read_input
  public :: read_vel_grid

contains

  subroutine read_input()

    use CommandLineProcessor
    
    implicit none

    character(len=128) :: input_filename, vel_grid_filename
    integer :: num_args

    call get_input_filename( input_filename, vel_grid_filename, num_args )

    if( num_args .eq. 1 ) call read_input_grvy( input_filename )
    if( num_args .eq. 2 ) call read_input_grvy( input_filename, vel_grid_filename )

    return
  end subroutine read_input

  subroutine read_input_grvy( input_filename, vel_grid_filename )

    use grvy

    implicit none

    character(len=128), intent(in) :: input_filename
    character(len=128), optional, intent(in) :: vel_grid_filename

    integer :: status, n

    double precision :: n_ref, T_ref
    character(len=16) :: species_ref

    integer :: num_rot_levels, num_vib_levels
    character(len=16), dimension(:), allocatable :: species_codes

    integer :: csection_method, csection_coeff

    double precision :: delta_x
    integer :: nspace, num_zones
    integer, dimension(:), allocatable :: reference, start_values

    integer :: x_grid_type, y_grid_type, z_grid_type
    integer :: num_points_x, num_points_y, num_points_z
    double precision :: xv_min, xv_max, yv_min, yv_max, zv_min, zv_max

    double precision :: deltat
    integer :: ntstep 

    integer, allocatable, dimension(:) :: InitC, shock_thickness
    integer, allocatable, dimension(:) :: InitDF, InitRot, InitVib
    integer :: LWallBC, RWallBC
    double precision, allocatable, dimension(:) :: x_delta_loc, y_delta_loc, z_delta_loc
    double precision :: dens_LW, x_vel_LW, tr_temp_LW, rot_temp_LW, vib_temp_LW
    double precision :: dens_RW, x_vel_RW, tr_temp_RW, rot_temp_RW, vib_temp_RW

    integer :: iColl, NRepl
    double precision, allocatable, dimension(:) :: ColnRMS

    integer :: iConv

    integer :: verbosity, display_freq, screen_output_loc

    integer :: recalc_equil_limit, RotColl, VibColl
    integer :: Z_rot_method, Z_vib_method

    integer :: enable_vel_df_vis, enable_rot_df_vis, enable_vib_df_vis
    integer :: enable_rot_levels_vis, enable_vib_levels_vis
    integer :: enable_properties_vis, enable_collisions_vis
    integer :: enable_bkw_analytic_vis
    character(len=128) :: df_vis_filename, rot_levels_vis_filename, vib_levels_vis_filename
    character(len=128) :: properties_vis_filename, collisions_vis_filename
    character(len=128) :: rotvib_levels_vis_filename
    integer :: df_vis_dump_freq, rot_levels_vis_dump_freq, vib_levels_vis_dump_freq
    integer :: properties_vis_dump_freq, collisions_vis_dump_freq
    integer :: rotvib_levels_vis_dump_freq
    integer :: num_output_locations
    integer, allocatable, dimension(:) :: output_locations
    integer :: df_output_dimension
    integer :: write_dens, write_x_vel, write_y_vel, write_z_vel, write_speed, write_mach
    integer :: write_tr_temp, write_rot_temp, write_vib_temp, write_tot_temp
    integer :: write_temp_x, write_temp_y, write_temp_z
    integer :: write_rot_dof, write_vib_dof
    integer :: write_tr_energy, write_rot_energy, write_vib_energy, write_tot_energy
    integer :: write_entropy, write_heat_flux_x, write_heat_flux_y, write_heat_flux_z
    integer :: write_shear_stress
    integer :: write_moment
    integer, allocatable, dimension(:) :: moment_numbers

    integer :: write_rst, write_freq, read_rst
    character(len=128) :: write_restart_filename, read_restart_filename

    integer :: rand_method, seed

    call grvy_input_fopen( trim(adjustl(input_filename)), status )

    !===============================================================================================
    ! Read reference information
    !===============================================================================================
    call grvy_input_fread_double( 'reference/n_ref', n_ref, status )
    call grvy_input_fread_double( 'reference/T_ref', T_ref, status )
    call grvy_input_fread_char( 'reference/species_ref', species_ref, status )

    call set_reference_parameters( T_ref, n_ref, species_ref )

    !===============================================================================================
    ! Read species data
    !===============================================================================================
    call grvy_input_fread_int( 'species/num_species', num_species, status )
    
    allocate( species_codes(1:num_species), STAT=status )
    call allocate_error_check( status, "species_codes input" )

    do n = 1, num_species
       call grvy_input_fread_char_ivec( "species/species_codes", species_codes(n), n, status )
    end do

    call set_species_codes( species_codes, num_species )

    call initialize_df_input_arrays()
    do n = 1, num_species
       call grvy_input_fread_int_ivec( 'species/num_rot_levels', num_rot_levels, n, status )
       call grvy_input_fread_int_ivec( 'species/num_vib_levels', num_vib_levels, n, status )

       call set_num_energy_levels( num_rot_levels, num_vib_levels, n )
    end do

    !===============================================================================================
    ! Read cross section model data
    !===============================================================================================
    call grvy_input_fread_int( 'csection_model/csection_method', csection_method, status )
    call grvy_input_fread_int( 'csection_model/csection_coeff', csection_coeff, status )

    call set_cross_section( csection_method, csection_coeff )

    !===============================================================================================
    ! Read physical grid definitions
    !===============================================================================================
    call grvy_input_fread_double( 'phy_grid/delta_x', delta_x, status )
    call grvy_input_fread_int( 'phy_grid/nspace', nspace, status )
    call grvy_input_fread_int( 'phy_grid/num_zones', num_zones, status )

    call set_physical_grid( nspace, delta_x, num_zones )

    call initialize_spatial_grid_reference()
    allocate( reference(1:num_zones), start_values(1:num_zones), STAT=status )
    call allocate_error_check( status, "reference, start_values" )

    call grvy_input_fread_int_vec( 'phy_grid/reference', reference, num_zones, status )
    call grvy_input_fread_int_vec( 'phy_grid/start_values', start_values, num_zones, status )

    call set_spatial_reference( reference, start_values )

    deallocate( reference, start_values, STAT=status )
    call deallocate_error_check( status, "reference, start_values" )

    !===============================================================================================
    ! Read velocity grid definitions
    !===============================================================================================
    call initialize_grid_input_arrays( num_species )
    
    do n = 1, num_species
       call grvy_input_fread_int_ivec( 'vel_grid/x_vel_grid_type', x_grid_type, n, status )
       call grvy_input_fread_int_ivec( 'vel_grid/y_vel_grid_type', y_grid_type, n, status )
       call grvy_input_fread_int_ivec( 'vel_grid/z_vel_grid_type', z_grid_type, n, status )
       
       call set_vel_grid_type( x_grid_type, y_grid_type, z_grid_type, n )
    end do

    if( present(vel_grid_filename) )then
       call read_vel_grid( vel_grid_filename, num_species )

    else
       do n = 1, num_species
          call grvy_input_fread_int_ivec( 'vel_grid/num_points_x', num_points_x, n, status )
          call grvy_input_fread_int_ivec( 'vel_grid/num_points_y', num_points_y, n, status )
          call grvy_input_fread_int_ivec( 'vel_grid/num_points_z', num_points_z, n, status )
          call set_vel_grid_num_points( num_points_x, num_points_y, num_points_z, n )

          call grvy_input_fread_double_ivec( 'vel_grid/xv_min', xv_min, n, status )
          call grvy_input_fread_double_ivec( 'vel_grid/xv_max', xv_max, n, status )
          call grvy_input_fread_double_ivec( 'vel_grid/yv_min', yv_min, n, status )
          call grvy_input_fread_double_ivec( 'vel_grid/yv_max', yv_max, n, status )
          call grvy_input_fread_double_ivec( 'vel_grid/zv_min', zv_min, n, status )
          call grvy_input_fread_double_ivec( 'vel_grid/zv_max', zv_max, n, status )
          call set_vel_grid_bounds( xv_min, xv_max, yv_min, yv_max, zv_min, zv_max, n )
       end do

    end if

    !===============================================================================================
    ! Read time-stepping parameters
    !===============================================================================================
    call grvy_input_fread_double( 'time_stepping/deltat', deltat, status )
    call grvy_input_fread_int( 'time_stepping/ntstep', ntstep, status )

    call set_deltat( deltat )
    call set_num_time_steps( ntstep )

    !===============================================================================================
    ! Read boundary conditions
    !===============================================================================================
    allocate( InitC(1:num_species), shock_thickness(1:num_species), STAT=status )
    call allocate_error_check( status, "InitC, shock_thickness" )

    call grvy_input_fread_int_vec( 'boundary/InitC', InitC, num_species, status )
    call grvy_input_fread_int_vec( 'boundary/shock_thickness', shock_thickness, num_species, status )
    call set_spatial_domain_conditions( InitC, shock_thickness )

    allocate( InitDF(1:num_species), InitRot(1:num_species), InitVib(1:num_species), STAT=status )
    call allocate_error_check( status, "InitDF, InitRot, InitVib" )

    call grvy_input_fread_int_vec( 'boundary/InitDF',  InitDF, num_species, status )
    call grvy_input_fread_int_vec( 'boundary/InitRot', InitRot, num_species, status )
    call grvy_input_fread_int_vec( 'boundary/InitVib', InitVib, num_species, status )
    call set_initial_df_conditions( InitDF )!, InitRot, InitVib )
    call set_initial_energy_df_conditions( InitRot, InitVib )

    allocate( x_delta_loc(1:num_species), y_delta_loc(1:num_species), z_delta_loc(1:num_species), STAT=status )
    call allocate_error_check( status, "delta_loc" )

    call grvy_input_fread_double_vec( 'boundary/x_delta_loc', x_delta_loc, num_species, status )
    call grvy_input_fread_double_vec( 'boundary/y_delta_loc', y_delta_loc, num_species, status )
    call grvy_input_fread_double_vec( 'boundary/z_delta_loc', z_delta_loc, num_species, status )
    call set_delta_function_location( x_delta_loc, y_delta_loc, z_delta_loc )

    call grvy_input_fread_int( 'boundary/LWallBC', LWallBC, status )
    call grvy_input_fread_int( 'boundary/RWallBC', RWallBC, status )
    call set_boundary_conditions_flags( LWallBC, RWallBC )

    call initialize_property_arrays()

    do n = 1, num_species
       call grvy_input_fread_double_ivec( 'boundary/dens_LW', dens_LW, n, status )
       call grvy_input_fread_double_ivec( 'boundary/x_vel_LW', x_vel_LW, n, status )
       call grvy_input_fread_double_ivec( 'boundary/tr_temp_LW', tr_temp_LW, n, status )
       call grvy_input_fread_double_ivec( 'boundary/rot_temp_LW', rot_temp_LW, n, status )
       call grvy_input_fread_double_ivec( 'boundary/vib_temp_LW', vib_temp_LW, n, status )

       call grvy_input_fread_double_ivec( 'boundary/dens_RW', dens_RW, n, status )
       call grvy_input_fread_double_ivec( 'boundary/x_vel_RW', x_vel_RW, n, status )
       call grvy_input_fread_double_ivec( 'boundary/tr_temp_RW', tr_temp_RW, n, status )
       call grvy_input_fread_double_ivec( 'boundary/rot_temp_RW', rot_temp_RW, n, status )
       call grvy_input_fread_double_ivec( 'boundary/vib_temp_RW', vib_temp_RW, n, status )

       call set_boundary_properties( dens_LW, x_vel_LW, tr_temp_LW, dens_RW, x_vel_RW, tr_temp_RW, n )
       call set_temp_rot( rot_temp_LW, rot_temp_RW, n )
       call set_temp_vib( vib_temp_LW, vib_temp_RW, n )
    end do

    deallocate( InitC, shock_thickness, STAT=status )
    call deallocate_error_check( status, "InitC, shock_thickness" )

    deallocate( InitDF, InitRot, InitVib, STAT=status )
    call deallocate_error_check( status, "InitDF, InitRot, InitVib" )

    deallocate( x_delta_loc, y_delta_loc, z_delta_loc, STAT=status )
    call deallocate_error_check( status, "delta_loc" )

    !===============================================================================================
    ! Read collision properties
    !===============================================================================================
    call grvy_input_fread_int( 'collision/iColl', iColl, status )
    call set_collision_integral_solver( iColl )

    call grvy_input_fread_int( 'collision/recalc_equil_limit', recalc_equil_limit, status )
    call set_recalculate_equil_limit( recalc_equil_limit )

    call grvy_input_fread_int( 'collision/RotColl', RotColl, status )
    call grvy_input_fread_int( 'collision/VibColl', VibColl, status )
    call set_energy_methods( RotColl, VibColl )

    call grvy_input_fread_int( 'collision/NRepl', NRepl, status )
    call set_num_repl_locations( NRepl )

    allocate( ColnRMS( 1:(num_species*num_species) ), STAT=status )
    call allocate_error_check( status, "ColnRMS (input)" )

    call grvy_input_fread_double_vec( 'collision/ColnRMS', ColnRMS, num_species*num_species, status )
    call initialize_coln_rms( ColnRMS )

    call grvy_input_fread_int( 'collision/Z_rot_method', Z_rot_method, status )
    call grvy_input_fread_int( 'collision/Z_vib_method', Z_vib_method, status )
    call set_relaxation_rate_method( Z_rot_method, Z_vib_method )

    deallocate( ColnRMS, STAT=status )
    call deallocate_error_check( status, "ColnRMS (input)" )

    !===============================================================================================
    ! Read Convection information
    !===============================================================================================
    call grvy_input_fread_int( 'convection/iConv', iConv, status )
    call set_convection_update( iConv )

    !===============================================================================================
    ! Read Screen Output information
    !===============================================================================================
    call grvy_input_fread_int( 'screen_out/verbosity', verbosity, status )
    call grvy_input_fread_int( 'screen_out/display_freq', display_freq, status )
    call grvy_input_fread_int( 'screen_out/screen_output_loc', screen_output_loc, status )

    call set_verbosity( verbosity )
    call set_display_printing_freq( display_freq )
    call set_screen_output_location( screen_output_loc )

    !===============================================================================================
    ! Read Visualization information
    !===============================================================================================
    call grvy_input_fread_int( 'vis/enable_vel_df_vis', enable_vel_df_vis, status )
    call grvy_input_fread_int( 'vis/enable_rot_df_vis', enable_rot_df_vis, status )
    call grvy_input_fread_int( 'vis/enable_vib_df_vis', enable_vib_df_vis, status )
    call grvy_input_fread_int( 'vis/enable_rot_levels_vis', enable_rot_levels_vis, status )
    call grvy_input_fread_int( 'vis/enable_vib_levels_vis', enable_vib_levels_vis, status )
    call grvy_input_fread_int( 'vis/enable_properties_vis', enable_properties_vis, status )
    call grvy_input_fread_int( 'vis/enable_collisions_vis', enable_collisions_vis, status )
    call enable_vis( enable_vel_df_vis, enable_rot_df_vis, enable_vib_df_vis, &
         enable_rot_levels_vis, enable_vib_levels_vis, enable_properties_vis, enable_collisions_vis )

    call grvy_input_fread_int( 'vis/enable_bkw_analytic_vis', enable_bkw_analytic_vis, status )
    call enable_bkw_analytic( enable_bkw_analytic_vis )

    call grvy_input_fread_char( 'vis/df_vis_filename', df_vis_filename, status )
    call grvy_input_fread_char( 'vis/rot_levels_vis_filename', rot_levels_vis_filename, status )
    call grvy_input_fread_char( 'vis/vib_levels_vis_filename', vib_levels_vis_filename, status )
    call grvy_input_fread_char( 'vis/properties_vis_filename', properties_vis_filename, status )
    call grvy_input_fread_char( 'vis/collisions_vis_filename', collisions_vis_filename, status )
    call grvy_input_fread_char( 'vis/rotvib_levels_vis_filename', rotvib_levels_vis_filename, status )
    call set_vis_filename( df_vis_filename, rot_levels_vis_filename, vib_levels_vis_filename, &
         properties_vis_filename, collisions_vis_filename, rotvib_levels_vis_filename )

    call grvy_input_fread_int( 'vis/df_vis_dump_freq', df_vis_dump_freq, status )
    call grvy_input_fread_int( 'vis/rot_levels_vis_dump_freq', rot_levels_vis_dump_freq, status )
    call grvy_input_fread_int( 'vis/vib_levels_vis_dump_freq', vib_levels_vis_dump_freq, status )
    call grvy_input_fread_int( 'vis/properties_vis_dump_freq', properties_vis_dump_freq, status )
    call grvy_input_fread_int( 'vis/collisions_vis_dump_freq', collisions_vis_dump_freq, status )
    call grvy_input_fread_int( 'vis/rotvib_levels_vis_dump_freq', rotvib_levels_vis_dump_freq, status )
    call set_vis_dump_freq( df_vis_dump_freq, rot_levels_vis_dump_freq, vib_levels_vis_dump_freq, &
         properties_vis_dump_freq, collisions_vis_dump_freq, rotvib_levels_vis_dump_freq )

    call grvy_input_fread_int( 'vis/num_output_locations', num_output_locations, status )

    allocate( output_locations(1:num_output_locations), STAT=status )
    call allocate_error_check( status, "output_locations (input)" )

    call grvy_input_fread_int_vec( 'vis/output_locations', output_locations, num_output_locations, status )
    call set_output_locations( num_output_locations, output_locations )

    call grvy_input_fread_int( 'vis/df_output_dimension', df_output_dimension, status )
    call set_df_output_dimension( df_output_dimension )

    call grvy_input_fread_int("vis/write_dens", write_dens, status )
    
    call grvy_input_fread_int("vis/write_x_vel", write_x_vel, status )
    call grvy_input_fread_int("vis/write_y_vel", write_y_vel, status )
    call grvy_input_fread_int("vis/write_z_vel", write_z_vel, status )
    call grvy_input_fread_int('vis/write_speed', write_speed, status )
    call grvy_input_fread_int("vis/write_mach", write_mach, status )

    call grvy_input_fread_int("vis/write_tr_temp", write_tr_temp, status )
    call grvy_input_fread_int("vis/write_rot_temp", write_rot_temp, status )
    call grvy_input_fread_int("vis/write_vib_temp", write_vib_temp, status )
    call grvy_input_fread_int("vis/write_tot_temp", write_tot_temp, status )
    call grvy_input_fread_int("vis/write_temp_x", write_temp_x, status )
    call grvy_input_fread_int("vis/write_temp_y", write_temp_y, status )
    call grvy_input_fread_int("vis/write_temp_z", write_temp_z, status )

    call grvy_input_fread_int('vis/write_tr_energy', write_tr_energy, status )
    call grvy_input_fread_int("vis/write_rot_energy", write_rot_energy, status )
    call grvy_input_fread_int("vis/write_vib_energy", write_vib_energy, status )
    call grvy_input_fread_int("vis/write_tot_energy", write_tot_energy, status )

    call grvy_input_fread_int("vis/write_rot_dof", write_rot_dof, status )
    call grvy_input_fread_int("vis/write_vib_dof", write_vib_dof, status )

    call grvy_input_fread_int("vis/write_entropy", write_entropy, status )
    call grvy_input_fread_int("vis/write_heat_flux_x", write_heat_flux_x, status )
    call grvy_input_fread_int("vis/write_heat_flux_y", write_heat_flux_y, status )
    call grvy_input_fread_int("vis/write_heat_flux_z", write_heat_flux_z, status )
    call grvy_input_fread_int("vis/write_shear_stress", write_shear_stress, status )

    call set_properties_output_vars( write_dens, &
         write_x_vel, write_y_vel, write_z_vel, write_speed, write_mach, &
         write_tr_temp, write_rot_temp, write_vib_temp, write_tot_temp, &
         write_temp_x, write_temp_y, write_temp_z, &
         write_rot_dof, write_vib_dof, &
         write_tr_energy, write_rot_energy, write_vib_energy, write_tot_energy, &
         write_entropy, write_heat_flux_x, write_heat_flux_y, write_heat_flux_z, write_shear_stress )
    
    call grvy_input_fread_int( 'vis/write_moment', write_moment, status )
    if( write_moment .gt. 0 )then
       allocate( moment_numbers(1:write_moment), STAT=status )
       call allocate_error_check( status, "moment_numbers (input)" )
       call grvy_input_fread_int_vec( 'vis/moment_numbers', moment_numbers, write_moment, status )
    end if

    call set_moment_output( write_moment, moment_numbers )

    deallocate( output_locations, STAT=status )
    call deallocate_error_check( status, "output_locations (input)" )

    if( write_moment .gt. 0 )then
       deallocate( moment_numbers, STAT=status )
       call deallocate_error_check( status, "moment_numbers (input)" )
    end if

    !===============================================================================================
    ! Read restart file info
    !===============================================================================================
    call grvy_input_fread_int( 'restart/write_rst', write_rst, status )
    call grvy_input_fread_int( 'restart/write_freq', write_freq, status )
    call grvy_input_fread_int( 'restart/read_rst', read_rst, status )
    call set_restart_flags( write_rst, write_freq, read_rst )
    
    call grvy_input_fread_char( 'restart/write_restart_filename', write_restart_filename, status )
    call grvy_input_fread_char( 'restart/read_restart_filename', read_restart_filename, status )
    call set_restart_filename( write_restart_filename, read_restart_filename )

    !===============================================================================================
    ! Read random number generation info
    !===============================================================================================
    call grvy_input_fread_int( 'rng/rand_method', rand_method, status )
    call set_rand_method( rand_method )

    call grvy_input_fread_int( 'rng/seed', seed, status )
    call set_seed( seed )

    return
  end subroutine read_input_grvy

  subroutine read_vel_grid( vel_grid_filename, num_species )

    implicit none

    character(len=128), intent(in) :: vel_grid_filename
    integer, intent(in) :: num_species
    integer, parameter :: vel_grid_file_unit = 1

    integer :: num_points_x, num_points_y, num_points_z
    integer :: i, n, status

    double precision, allocatable, dimension(:) :: x, y, z
 
    open( unit=vel_grid_file_unit, file=trim(adjustl(vel_grid_filename)), status='old' )

    do n = 1, num_species

       call initialize_file_read_arrays( n )

       read(vel_grid_file_unit,*) num_points_x
       read(vel_grid_file_unit,*) num_points_y
       read(vel_grid_file_unit,*) num_points_z
       
       call set_vel_grid_num_points( num_points_x, num_points_y, num_points_z, n )

       allocate( x(1:num_points_x), y(1:num_points_y), z(1:num_points_z), STAT=status )
       call allocate_error_check( status, "x, y, z (input)" )

       do i = 1, num_points_x
          read(vel_grid_file_unit,*) x(i)
       end do
       
       do i = 1, num_points_y
          read(vel_grid_file_unit,*) y(i)
       end do
       
       do i = 1, num_points_z
          read(vel_grid_file_unit,*) z(i)
       end do

       call set_velocity_vectors( x, y, z, n )

       deallocate( x, y, z, STAT=status )
       call deallocate_error_check( status, "x, y, z (input)" )

    end do

    close( vel_grid_file_unit )

    return
  end subroutine read_vel_grid

end module ReadInput
