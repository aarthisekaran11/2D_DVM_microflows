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
  use VarianceReduction
  use input
  use strings

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

!    use grvy

    implicit none

    character(len=128), intent(in) :: input_filename
    character(len=128), optional, intent(in) :: vel_grid_filename
    character(len=40),dimension(:,:),allocatable :: args
    integer, parameter :: inp_line = 500
    integer :: inp_line_true

    integer :: status, n

    double precision :: Kn

    double precision :: n_ref, T_ref
    character(len=16) :: species_ref

    integer,dimension(:), allocatable :: num_rot_levels, num_vib_levels
    character(len=16), dimension(:), allocatable :: species_codes

    integer :: csection_method, csection_coeff

    double precision :: Lc
    double precision :: delta_x, delta_y
    integer :: nx_space, ny_space

    integer, dimension(:), allocatable :: x_grid_type, y_grid_type, z_grid_type
    integer, dimension(:), allocatable :: num_points_x, num_points_y, num_points_z
    double precision,dimension(:), allocatable :: xv_min, xv_max, yv_min, yv_max, zv_min, zv_max

    double precision :: t_final, t_frac
    double precision :: deltat
    integer :: ntstep 

    integer, allocatable, dimension(:) :: InitC, shock_thickness, spatial_middle
    integer, allocatable, dimension(:) :: InitDF, InitRot, InitVib
    integer :: LWallBC, RWallBC, BWallBC, TWallBC
    double precision :: acc_coeff
    double precision, allocatable, dimension(:) :: x_delta_loc, y_delta_loc, z_delta_loc
    double precision, dimension(:), allocatable :: dens_LW, u_LW, v_LW, w_LW, tr_temp_LW, rot_temp_LW, vib_temp_LW
    double precision, dimension(:), allocatable :: dens_RW, u_RW, v_RW, w_RW, tr_temp_RW, rot_temp_RW, vib_temp_RW
    double precision, dimension(:), allocatable :: dens_BW, u_BW, v_BW, w_BW, tr_temp_BW, rot_temp_BW, vib_temp_BW
    double precision, dimension(:), allocatable :: dens_TW, u_TW, v_TW, w_TW, tr_temp_TW, rot_temp_TW, vib_temp_TW
    double precision :: dens_init, u_init, v_init, w_init, tr_temp_init, rot_temp_init, vib_temp_init

    integer :: iColl, NRepl
    double precision, allocatable, dimension(:) :: ColnRMS

    integer :: iConv

    integer :: verbosity, display_freq!, screen_output_loc
    integer, dimension(2) :: screen_output_loc

    integer :: RotColl, VibColl
    double precision :: cdf_tolerance
    integer :: Z_rot_method, Z_vib_method

    integer :: enable_vel_df_vis, enable_rot_df_vis, enable_vib_df_vis
    integer :: enable_rot_levels_vis, enable_vib_levels_vis
    integer :: enable_properties_vis, enable_collisions_vis
    integer :: enable_speed_vis, enable_speed_rot_vis, enable_speed_vib_vis
    integer :: enable_bkw_analytic_vis
    character(len=128) :: df_vis_filename, rot_levels_vis_filename, vib_levels_vis_filename
    character(len=128) :: properties_vis_filename, collisions_vis_filename
    character(len=128) :: speed_vis_filename
    integer :: df_vis_dump_freq, rot_levels_vis_dump_freq, vib_levels_vis_dump_freq
    integer :: properties_vis_dump_freq, collisions_vis_dump_freq, speed_vis_dump_freq
    integer :: wall_heat_flux_dump_freq
    integer :: num_output_locations
    integer, allocatable, dimension(:) :: x_output_locations, y_output_locations
    integer :: df_output_dimension
    integer :: write_dens, write_x_vel, write_y_vel, write_z_vel, write_speed, write_mach
    integer :: write_tr_temp, write_rot_temp, write_vib_temp, write_tot_temp
    integer :: write_temp_x, write_temp_y, write_temp_z, write_pressure
    integer :: write_rot_dof, write_vib_dof
    integer :: write_tr_energy, write_rot_energy, write_vib_energy, write_tot_energy
    integer :: write_entropy, write_heat_flux_x, write_heat_flux_y, write_heat_flux_z
    integer :: write_kin_temp_relax, write_rot_temp_relax, write_vib_temp_relax
    integer :: write_stress, write_moment
    integer, allocatable, dimension(:) :: stress_directions, moment_numbers

    integer :: write_rst, write_freq, read_rst
    character(len=128) :: write_restart_filename, read_restart_filename

    integer :: rand_method, seed


    allocate(args(2,inp_line),STAT=status
    call allocate_error_check(status,'args')
!    call grvy_input_fopen( trim(adjustl(input_filename)), status )
     call file_read(inp_line,trim(adjustl(input_filename)),args,inp_line_true)
    !===============================================================================================
    ! Read reference information
    !===============================================================================================
    call variable_read(args,'Kn',Kn,inp_line_true)

    call set_knudsen_number( Kn )

    call variable_read(args,'Lc',Lc,inp_line_true)


    call set_Lc( Lc )

    call variable_read(args,'n_ref',n_ref,inp_line_true)
    call variable_read(args,'T_ref',T_ref,inp_line_true)
    call variable_read(args,'species_ref',species_ref,inp_line_true)



    call set_reference_parameters( T_ref, n_ref, species_ref )

    !===============================================================================================
    ! Read species data
    !===============================================================================================
    call variable_read(args,'num_species',num_species,inp_line_true)

    allocate( species_codes(1:num_species), STAT=status )
    call allocate_error_check( status, "species_codes input" )


    call variable_read_vec(args,'species_codes',species_codes,inp_line_true,num_species)

    call set_species_codes( species_codes, num_species )

    call initialize_df_input_arrays()

    allocate(num_rot_levels(1:num_species),STAT=status)
    call allocate_error_check(status,'num_rot_levels')    

    allocate(num_vib_levels(1:num_species),STAT=status)
    call allocate_error_check(status,'num_vib_levels')    


    call variable_read_vec(args,'num_rot_levels',num_rot_levels,inp_line_true,num_species)
    call variable_read_vec(args,'num_vib_levels',num_vib_levels,inp_line_true,num_species)
    do n = 1, num_species
       call set_num_energy_levels( num_rot_levels(n), num_vib_levels(n), n )
    end do

    !===============================================================================================
    ! Read cross section model data
    !===============================================================================================
    call variable_read(args,'csection_method',csection_method,inp_line_true)
    call variable_read(args,'csection_coeff',csection_coeff,inp_line_true)
    call set_cross_section( csection_method, csection_coeff )

    !===============================================================================================
    ! Read physical grid definitions
    !===============================================================================================
    call variable_read(args,'delta_x',delta_x,inp_line_true)
    call variable_read(args,'nx_space',nx_space,inp_line_true)
    call variable_read(args,'delta_y',delta_y,inp_line_true)
    call variable_read(args,'ny_space',ny_space,inp_line_true)



    call set_physical_grid( nx_space, delta_x, ny_space, delta_y )

    !===============================================================================================
    ! Read velocity grid definitions
    !===============================================================================================
    call initialize_grid_input_arrays( num_species )

    allocate(x_grid_type(1:num_species),STAT=status)
    call allocate_error_check(status,'x_vel_grid_tpye')   

    allocate(y_grid_type(1:num_species),STAT=status)
    call allocate_error_check(status,'y_vel_grid_tpye')   

    allocate(z_grid_type(1:num_species),STAT=status)
    call allocate_error_check(status,'z_vel_grid_tpye')   
    
    call variable_read_vec(args,'x_vel_grid_type',x_grid_type,inp_line_true,num_species)    
    call variable_read_vec(args,'y_vel_grid_type',y_grid_type,inp_line_true,num_species)    
    call variable_read_vec(args,'z_vel_grid_type',z_grid_type,inp_line_true,num_species)    



    do n = 1, num_species 
       call set_vel_grid_type( x_grid_type(n), y_grid_type(n), z_grid_type(n), n )
    end do

    deallocate(x_grid_type,y_grid_type,z_grid_type,STAT=status)
    call deallocate_error_check(status,'grid_type')

    if( present(vel_grid_filename) )then
       call read_vel_grid( vel_grid_filename, num_species )

    else

    allocate(num_points_x(1:num_species),STAT=status)
    call allocate_error_check(status,'num_points_x')   

    allocate(num_points_y(1:num_species),STAT=status)
    call allocate_error_check(status,'num_points_y')   

    allocate(num_points_z(1:num_species),STAT=status)
    call allocate_error_check(status,'num_points_z')   

    call variable_read_vec(args,'num_points_x',num_points_x,inp_line_true,num_species)
    call variable_read_vec(args,'num_points_y',num_points_y,inp_line_true,num_species)
    call variable_read_vec(args,'num_points_z',num_points_z,inp_line_true,num_species)

    allocate(xv_min(1:num_species),STAT=status)
    call allocate_error_check(status,'xv_min')   

    allocate(xv_max(1:num_species),STAT=status)
    call allocate_error_check(status,'xv_max')   

    allocate(yv_min(1:num_species),STAT=status)
    call allocate_error_check(status,'yv_min')   

    allocate(yv_max(1:num_species),STAT=status)
    call allocate_error_check(status,'yv_max')   

    allocate(zv_min(1:num_species),STAT=status)
    call allocate_error_check(status,'zv_min')   

    allocate(zv_max(1:num_species),STAT=status)
    call allocate_error_check(status,'zv_max')   

    call variable_read_vec(args,'xv_min',xv_min,inp_line_true,num_species)
    call variable_read_vec(args,'xv_max',xv_max,inp_line_true,num_species)
    call variable_read_vec(args,'yv_min',yv_min,inp_line_true,num_species)
    call variable_read_vec(args,'yv_max',yv_max,inp_line_true,num_species)
    call variable_read_vec(args,'zv_min',zv_min,inp_line_true,num_species)
    call variable_read_vec(args,'zv_max',zv_max,inp_line_true,num_species)



       do n = 1, num_species
          call set_vel_grid_num_points( num_points_x(n), num_points_y(n), num_points_z(n), n )
!
          call set_vel_grid_bounds( xv_min(n), xv_max(n), yv_min(n), yv_max(n), zv_min(n), zv_max(n), n )
       end do

    end if


    deallocate(num_points_x,num_points_y,num_points_z,STAT=status)
    call deallocate_error_check(status,'num_points')

    deallocate(xv_min,xv_max,yv_min,yv_max,zv_min,zv_max,STAT=status)
    call deallocate_error_check(status,'v_min/max')
    !===============================================================================================
    ! Read time-stepping parameters
    !===============================================================================================
    call variable_read(args,'deltat',deltat,inp_line_true)
    call variable_read(args,'ntstep',ntstep,inp_line_true)
!    call grvy_input_fread_double( 'time_stepping/deltat', deltat, status )
!    call grvy_input_fread_int( 'time_stepping/ntstep', ntstep, status )

    call set_deltat( deltat )
    call set_num_time_steps( ntstep )

    call variable_read(args,'t_final',t_final,inp_line_true)
    call variable_read(args,'t_final',t_final,inp_line_true)
!    call grvy_input_fread_double( 'time_stepping/t_final', t_final, status )
!    call grvy_input_fread_double( 'time_stepping/t_frac', t_frac, status )

    call set_times( t_final, t_frac )

    !===============================================================================================
    ! Read boundary conditions
    !===============================================================================================
    allocate( InitC(1:num_species), shock_thickness(1:num_species), spatial_middle(1:num_species), STAT=status )
    call allocate_error_check( status, "InitC, shock_thickness" )

    call variable_read_vec(args,'InitC',InitC,inp_line_true,num_species)
    call variable_read_vec(args,'shock_thickness',shock_thickness,inp_line_true,num_species)
    call variable_read_vec(args,'spatial_middle',spatial_middle,inp_line_true,num_species)


    call set_spatial_domain_conditions( InitC, shock_thickness, spatial_middle )

    allocate( InitDF(1:num_species), InitRot(1:num_species), InitVib(1:num_species), STAT=status )
    call allocate_error_check( status, "InitDF, InitRot, InitVib" )


    call variable_read_vec(args,'InitDF',InitDF,inp_line_true,num_species)
    call variable_read_vec(args,'InitRot',InitRot,inp_line_true,num_species)
    call variable_read_vec(args,'InitVib',InitVib,inp_line_true,num_species)


    call set_initial_df_conditions( InitDF )!, InitRot, InitVib )
    call set_initial_energy_df_conditions( InitRot, InitVib )

    allocate( x_delta_loc(1:num_species), y_delta_loc(1:num_species), z_delta_loc(1:num_species), STAT=status )
    call allocate_error_check( status, "delta_loc" )

    call variable_read_vec(args,'x_delta_loc',x_delta_loc,inp_line_true,num_species)
    call variable_read_vec(args,'y_delta_loc',y_delta_loc,inp_line_true,num_species)
    call variable_read_vec(args,'z_delta_loc',z_delta_loc,inp_line_true,num_species)

    call set_delta_function_location( x_delta_loc, y_delta_loc, z_delta_loc )

    call variable_read(args,'LWallBC',LWallBC,inp_line_true)
    call variable_read(args,'RWallBC',RWallBC,inp_line_true)
    call variable_read(args,'BWallBC',BWallBC,inp_line_true)
    call variable_read(args,'TWallBC',TWallBC,inp_line_true)

    call set_boundary_conditions_flags( LWallBC, RWallBC, BWallBC, TWallBC )

    call variable_read(args,'acc_coeff',acc_coeff,inp_line_true)
!    call grvy_input_fread_double( 'boundary/acc_coeff', acc_coeff, status )
    call set_accomodation_coefficient( acc_coeff )

    call initialize_property_arrays()
    call initialize_domain_init_arrays()

    allocate(dens_LW(1:num_species),u_LW(1:num_species),v_LW(1:num_species),w_LW(1:num_species),&
    tr_temp_LW(1:num_species),rot_temp_LW(1:num_species),vib_temp_LW(1:num_species),STAT=status)

    call allocate_error_check(status,'dens_LW,u_LW,v_LW,w_LW,tr_temp_LW,rot_temp_LW,vib_temp_LW')
   
    allocate(dens_RW(1:num_species),u_RW(1:num_species),v_RW(1:num_species),w_RW(1:num_species),&
    tr_temp_RW(1:num_species),rot_temp_RW(1:num_species),vib_temp_RW(1:num_species),STAT=status)

    call allocate_error_check(status,'dens_RW,u_RW,v_RW,w_RW,tr_temp_RW,rot_temp_RW,vib_temp_RW')

    allocate(dens_BW(1:num_species),u_BW(1:num_species),v_BW(1:num_species),w_BW(1:num_species),&
    tr_temp_BW(1:num_species),rot_temp_BW(1:num_species),vib_temp_BW(1:num_species),STAT=status)

    call allocate_error_check(status,'dens_BW,u_BW,v_BW,w_BW,tr_temp_BW,rot_temp_BW,vib_temp_BW')

    allocate(dens_TW(1:num_species),u_TW(1:num_species),v_TW(1:num_species),w_TW(1:num_species),&
    tr_temp_TW(1:num_species),rot_temp_TW(1:num_species),vib_temp_TW(1:num_species),STAT=status)

    call allocate_error_check(status,'dens_TW,u_TW,v_TW,w_TW,tr_temp_TW,rot_temp_TW,vib_temp_TW')



    call variable_read_vec(args,'dens_LW',dens_LW,inp_line_true,num_species)
    call variable_read_vec(args,'x_vel_LW',u_LW,inp_line_true,num_species)
    call variable_read_vec(args,'y_vel_LW',v_LW,inp_line_true,num_species)
    call variable_read_vec(args,'z_vel_LW',w_LW,inp_line_true,num_species)
    call variable_read_vec(args,'tr_temp_LW',tr_temp_LW,inp_line_true,num_species)
    call variable_read_vec(args,'rot_temp_LW',rot_temp_LW,inp_line_true,num_species)
    call variable_read_vec(args,'vib_temp_LW',vib_temp_LW,inp_line_true,num_species)

    call variable_read_vec(args,'dens_RW',dens_RW,inp_line_true,num_species)
    call variable_read_vec(args,'x_vel_RW',u_RW,inp_line_true,num_species)
    call variable_read_vec(args,'y_vel_RW',v_RW,inp_line_true,num_species)
    call variable_read_vec(args,'z_vel_RW',w_RW,inp_line_true,num_species)
    call variable_read_vec(args,'tr_temp_RW',tr_temp_RW,inp_line_true,num_species)
    call variable_read_vec(args,'rot_temp_RW',rot_temp_RW,inp_line_true,num_species)
    call variable_read_vec(args,'vib_temp_RW',vib_temp_RW,inp_line_true,num_species)

    call variable_read_vec(args,'dens_BW',dens_BW,inp_line_true,num_species)
    call variable_read_vec(args,'x_vel_BW',u_BW,inp_line_true,num_species)
    call variable_read_vec(args,'y_vel_BW',v_BW,inp_line_true,num_species)
    call variable_read_vec(args,'z_vel_BW',w_BW,inp_line_true,num_species)
    call variable_read_vec(args,'tr_temp_BW',tr_temp_BW,inp_line_true,num_species)
    call variable_read_vec(args,'rot_temp_BW',rot_temp_BW,inp_line_true,num_species)
    call variable_read_vec(args,'vib_temp_BW',vib_temp_BW,inp_line_true,num_species)

    call variable_read_vec(args,'dens_TW',dens_TW,inp_line_true,num_species)
    call variable_read_vec(args,'x_vel_TW',u_TW,inp_line_true,num_species)
    call variable_read_vec(args,'y_vel_TW',v_TW,inp_line_true,num_species)
    call variable_read_vec(args,'z_vel_TW',w_TW,inp_line_true,num_species)
    call variable_read_vec(args,'tr_temp_TW',tr_temp_TW,inp_line_true,num_species)
    call variable_read_vec(args,'rot_temp_TW',rot_temp_TW,inp_line_true,num_species)
    call variable_read_vec(args,'vib_temp_TW',vib_temp_TW,inp_line_true,num_species)

    call variable_read_vec(args,'dens_init',dens_init,inp_line_true,num_species)
    call variable_read_vec(args,'x_vel_init',u_init,inp_line_true,num_species)
    call variable_read_vec(args,'y_vel_init',v_init,inp_line_true,num_species)
    call variable_read_vec(args,'z_vel_init',w_init,inp_line_true,num_species)
    call variable_read_vec(args,'tr_temp_init',tr_temp_init,inp_line_true,num_species)
    call variable_read_vec(args,'rot_temp_init',rot_temp_init,inp_line_true,num_species)
    call variable_read_vec(args,'vib_temp_init',vib_temp_init,inp_line_true,num_species)


    do n = 1, num_species 
       call set_domain_init_properties( dens_init, u_init, v_init, w_init, &
            tr_temp_init, rot_temp_init, vib_temp_init, n)

       call set_boundary_properties( dens_LW, u_LW, v_LW, w_LW, tr_temp_LW, &
            dens_RW, u_RW, v_RW, w_RW, tr_temp_RW, dens_BW, u_BW, v_BW, w_BW, &
            tr_temp_BW, dens_TW, u_TW, v_TW, w_TW, tr_temp_TW, n )
       call set_temp_rot( rot_temp_LW, rot_temp_RW, rot_temp_BW, rot_temp_TW, n )
       call set_temp_vib( vib_temp_LW, vib_temp_RW, vib_temp_BW, vib_temp_TW, n )
    end do

    deallocate( InitC, shock_thickness, STAT=status )
    call deallocate_error_check( status, "InitC, shock_thickness" )

    deallocate( InitDF, InitRot, InitVib, STAT=status )
    call deallocate_error_check( status, "InitDF, InitRot, InitVib" )

    deallocate( x_delta_loc, y_delta_loc, z_delta_loc, STAT=status )
    call deallocate_error_check( status, "delta_loc" )
    deallocate(dens_LW,u_LW,v_LW,w_LW,tr_temp_LW,rot_temp_LW,vib_temp_LW,STAT=status)
    call deallocate_error_check(status,'LW array')

    deallocate(dens_RW,u_RW,v_RW,w_RW,tr_temp_RW,rot_temp_RW,vib_temp_RW,STAT=status)
    call deallocate_error_check(status,'RW array')

    deallocate(dens_BW,u_BW,v_BW,w_BW,tr_temp_BW,rot_temp_BW,vib_temp_BW,STAT=status)
    call deallocate_error_check(status,'BW array')

    deallocate(dens_TW,u_TW,v_TW,w_TW,tr_temp_TW,rot_temp_TW,vib_temp_TW,STAT=status)
    call deallocate_error_check(status,'TW array')


    !===============================================================================================
    ! Read collision properties
    !===============================================================================================
    call variable_read(args,'iColl',iColl,inp_line_true)

    call set_collision_integral_solver( iColl )


!!$    call grvy_input_fread_int( 'collision/coll_limit', coll_limit, status )
!!$    call set_coll_limit( coll_limit )

    call variable_read(args,'cdf_tolerance',cdf_tolerance,inp_line_true)

    call set_cdf_tolerance( cdf_tolerance )

    call variable_read(args,'RotColl',RotColl,inp_line_true)
    call variable_read(args,'VibColl',VibColl,inp_line_true)

    call set_energy_methods( RotColl, VibColl )

    call variable_read(args,'NRepl',NRepl,inp_line_true)

    call set_num_repl_locations( NRepl )


    allocate( ColnRMS( 1:(num_species*num_species) ), STAT=status )
    call allocate_error_check( status, "ColnRMS (input)" )

    call variable_read_vec(args,'ColnRMS',ColnRMS,inp_line_true,num_species)

    call initialize_coln_rms( ColnRMS )


    call variable_read(args,'Z_rot_method',Z_rot_method,inp_line_true)
    call variable_read(args,'Z_vib_method',Z_vib_method,inp_line_true)

    call set_relaxation_rate_method( Z_rot_method, Z_vib_method )

    deallocate( ColnRMS, STAT=status )
    call deallocate_error_check( status, "ColnRMS (input)" )

    !===============================================================================================
    ! Read Convection information
    !===============================================================================================

    call variable_read(args,'iConv',iConv,inp_line_true)

    call set_convection_update( iConv )

    !===============================================================================================
    ! Read Screen Output information
    !===============================================================================================

    call variable_read(args,'verbosity',verbosity,inp_line_true)
    call variable_read(args,'display_freq',display_freq,inp_line_true)
    call variable_read_vec(args,'screen_output_loc',screen_output_loc,inp_line_true,num_species)
    !reasonably confident the above is the correct way to do it, will have to
    !see
!    call grvy_input_fread_int( 'screen_out/verbosity', verbosity, status )
!    call grvy_input_fread_int( 'screen_out/display_freq', display_freq, status )
!    call grvy_input_fread_int( 'screen_out/screen_output_loc', screen_output_loc, status )
!    call grvy_input_fread_int_ivec( 'screen_out/screen_output_loc', screen_output_loc(1), 1, status )
!    call grvy_input_fread_int_ivec( 'screen_out/screen_output_loc', screen_output_loc(2), 2, status )

    call set_verbosity( verbosity )
    call set_display_printing_freq( display_freq )
    call set_screen_output_location( screen_output_loc )

    !===============================================================================================
    ! Read Visualization information
    !===============================================================================================
    call variable_read(args,'enable_vel_df_vis',enable_vel_df_vis,inp_line_true)
    call variable_read(args,'enable_rot_df_vis',enable_rot_df_vis,inp_line_true)
    call variable_read(args,'enable_vib_df_vis',enable_vib_df_vis,inp_line_true)
    call variable_read(args,'enable_rot_levels_vis',enable_rot_levels_vis,inp_line_true)
    call variable_read(args,'enable_vib_levels_vis',enable_vib_levels_vis,inp_line_true)
    call variable_read(args,'enable_properties_vis',enable_properties_vis,inp_line_true)
    call variable_read(args,'enable_collisions_vis',enable_collisions_vis,inp_line_true)
    call variable_read(args,'enable_speed_vis',enable_speed_vis,inp_line_true)
    call variable_read(args,'enable_speed_rot_vis',enable_speed_rot_vis,inp_line_true)
    call variable_read(args,'enable_speed_vib_vis',enable_speed_vib_vis,inp_line_true)

!    call grvy_input_fread_int( 'vis/enable_vel_df_vis', enable_vel_df_vis, status )
!    call grvy_input_fread_int( 'vis/enable_rot_df_vis', enable_rot_df_vis, status )
!    call grvy_input_fread_int( 'vis/enable_vib_df_vis', enable_vib_df_vis, status )
!    call grvy_input_fread_int( 'vis/enable_rot_levels_vis', enable_rot_levels_vis, status )
!    call grvy_input_fread_int( 'vis/enable_vib_levels_vis', enable_vib_levels_vis, status )
!    call grvy_input_fread_int( 'vis/enable_properties_vis', enable_properties_vis, status )
!    call grvy_input_fread_int( 'vis/enable_collisions_vis', enable_collisions_vis, status )
!    call grvy_input_fread_int( 'vis/enable_speed_vis', enable_speed_vis, status )
!    call grvy_input_fread_int( 'vis/enable_speed_rot_vis', enable_speed_rot_vis, status )
!    call grvy_input_fread_int( 'vis/enable_speed_vib_vis', enable_speed_vib_vis, status )

    call enable_vis( enable_vel_df_vis, enable_rot_df_vis, enable_vib_df_vis, &
         enable_rot_levels_vis, enable_vib_levels_vis, enable_properties_vis, enable_collisions_vis, &
         enable_speed_vis, enable_speed_rot_vis, enable_speed_vib_vis )

    call variable_read(args,'enable_bkw_analytic_vis',enable_bkw_analytic_vis,inp_line_true)

!    call grvy_input_fread_int( 'vis/enable_bkw_analytic_vis', enable_bkw_analytic_vis, status )
    call enable_bkw_analytic( enable_bkw_analytic_vis )

    call variable_read(args,'df_vis_filename',df_vis_filename,inp_line_true)
    call variable_read(args,'rot_levels_vis_filename',rot_levels_vis_filename,inp_line_true)
    call variable_read(args,'vib_levels_vis_filename',vib_levels_vis_filename,inp_line_true)
    call variable_read(args,'properties_vis_filename',properties_vis_filename,inp_line_true)
    call variable_read(args,'collisions_vis_filename',collisions_vis_filename,inp_line_true)
    call variable_read(args,'speed_vis_filename',speed_vis_filename,inp_line_true)
!    call grvy_input_fread_char( 'vis/df_vis_filename', df_vis_filename, status )
!    call grvy_input_fread_char( 'vis/rot_levels_vis_filename', rot_levels_vis_filename, status )
!    call grvy_input_fread_char( 'vis/vib_levels_vis_filename', vib_levels_vis_filename, status )
!    call grvy_input_fread_char( 'vis/properties_vis_filename', properties_vis_filename, status )
!    call grvy_input_fread_char( 'vis/collisions_vis_filename', collisions_vis_filename, status )
!    call grvy_input_fread_char( 'vis/speed_vis_filename', speed_vis_filename, status )
    call set_vis_filename( df_vis_filename, rot_levels_vis_filename, vib_levels_vis_filename, &
         properties_vis_filename, collisions_vis_filename, speed_vis_filename )

    call variable_read(args,'df_vis_dump_freq',df_vis_dump_freq,inp_line_true)
    call variable_read(args,'rot_levels_vis_dump_freq',rot_levels_vis_dump_freq,inp_line_true)
    call variable_read(args,'vib_levels_vis_dump_freq',vib_levels_vis_dump_freq,inp_line_true)
    call variable_read(args,'properties_vis_dump_freq',properties_vis_dump_freq,inp_line_true)
    call variable_read(args,'collisions_vis_dump_freq',collisions_vis_dump_freq,inp_line_true)
    call variable_read(args,'speed_vis_dump_freq',speed_vis_dump_freq,inp_line_true)
!    call grvy_input_fread_int( 'vis/df_vis_dump_freq', df_vis_dump_freq, status )
!    call grvy_input_fread_int( 'vis/rot_levels_vis_dump_freq', rot_levels_vis_dump_freq, status )
!    call grvy_input_fread_int( 'vis/vib_levels_vis_dump_freq', vib_levels_vis_dump_freq, status )
!    call grvy_input_fread_int( 'vis/properties_vis_dump_freq', properties_vis_dump_freq, status )
!    call grvy_input_fread_int( 'vis/collisions_vis_dump_freq', collisions_vis_dump_freq, status )
!    call grvy_input_fread_int( 'vis/speed_vis_dump_freq', speed_vis_dump_freq, status )
    call variable_read(args,'wall_heat_flux_dump_freq',wall_heat_flux_dump_freq,inp_line_true)
!    call grvy_input_fread_int( 'vis/wall_heat_flux_dump_freq', wall_heat_flux_dump_freq, status )
    call set_vis_dump_freq( df_vis_dump_freq, rot_levels_vis_dump_freq, vib_levels_vis_dump_freq, &
         properties_vis_dump_freq, collisions_vis_dump_freq, speed_vis_dump_freq )
    call set_wall_heat_flux_dump_freq( wall_heat_flux_dump_freq )

    call variable_read(args,'num_output_locations',num_output_locations,inp_line_true)
!    call grvy_input_fread_int( 'vis/num_output_locations', num_output_locations, status )
!    call grvy_input_fread_int_ivec( 'vis/num_rot_levels', num_rot_levels, n, status )

    allocate( x_output_locations(1:num_output_locations), STAT=status )
    call allocate_error_check( status, "x_output_locations (input)" )

    allocate( y_output_locations(1:num_output_locations), STAT=status )
    call allocate_error_check( status, "y_output_locations (input)" )
    call variable_read_vec(args,'x_output_locations',x_output_locations,inp_line_true,num_output_locations')
    call variable_read_vec(args,'y_output_locations',y_output_locations,inp_line_true,num_output_locations')

!    call grvy_input_fread_int_vec( 'vis/x_output_locations', x_output_locations, num_output_locations, status )
!    call grvy_input_fread_int_vec( 'vis/y_output_locations', y_output_locations, num_output_locations, status )

    call set_output_locations( num_output_locations, x_output_locations, y_output_locations )

    call variable_read(args,'df_output_dimension',df_output_dimension,inp_line_true)
!    call grvy_input_fread_int( 'vis/df_output_dimension', df_output_dimension, status )
    call set_df_output_dimension( df_output_dimension )

    call variable_read(args,'write_dens',write_dens,inp_line_true)
!    call grvy_input_fread_int("vis/write_dens", write_dens, status )
    
    call variable_read(args,'write_x_vel',write_x_vel,inp_line_true)
    call variable_read(args,'write_y_vel',write_y_vel,inp_line_true)
    call variable_read(args,'write_z_vel',write_z_vel,inp_line_true)
    call variable_read(args,'write_speed',write_speed,inp_line_true)
    call variable_read(args,'write_mach',write_mach,inp_line_true)
!    call grvy_input_fread_int("vis/write_x_vel", write_x_vel, status )
!    call grvy_input_fread_int("vis/write_y_vel", write_y_vel, status )
!    call grvy_input_fread_int("vis/write_z_vel", write_z_vel, status )
!    call grvy_input_fread_int('vis/write_speed', write_speed, status )
!    call grvy_input_fread_int("vis/write_mach", write_mach, status )

    call variable_read(args,'write_tr_temp',write_tr_temp,inp_line_true)
    call variable_read(args,'write_rot_temp',write_rot_temp,inp_line_true)
    call variable_read(args,'write_vib_temp',write_vib_temp,inp_line_true)
    call variable_read(args,'write_tot_temp',write_tot_temp,inp_line_true)
    call variable_read(args,'write_temp_x',write_temp_x,inp_line_true)
    call variable_read(args,'write_temp_y',write_temp_y,inp_line_true)
    call variable_read(args,'write_temp_z',write_temp_z,inp_line_true)
!    call grvy_input_fread_int("vis/write_tr_temp", write_tr_temp, status )
!    call grvy_input_fread_int("vis/write_rot_temp", write_rot_temp, status )
!    call grvy_input_fread_int("vis/write_vib_temp", write_vib_temp, status )
!    call grvy_input_fread_int("vis/write_tot_temp", write_tot_temp, status )
!    call grvy_input_fread_int("vis/write_temp_x", write_temp_x, status )
!    call grvy_input_fread_int("vis/write_temp_y", write_temp_y, status )
!    call grvy_input_fread_int("vis/write_temp_z", write_temp_z, status )

    call variable_read(args,'write_pressure',write_pressure,inp_line_true)
!    call grvy_input_fread_int("vis/write_pressure", write_pressure, status )

    call variable_read(args,'write_tr_energy',write_tr_energy,inp_line_true)
    call variable_read(args,'write_rot_energy',write_rot_energy,inp_line_true)
    call variable_read(args,'write_vib_energy',write_vib_energy,inp_line_true)
    call variable_read(args,'write_tot_energy',write_tot_energy,inp_line_true)
!    call grvy_input_fread_int('vis/write_tr_energy', write_tr_energy, status )
!    call grvy_input_fread_int("vis/write_rot_energy", write_rot_energy, status )
!    call grvy_input_fread_int("vis/write_vib_energy", write_vib_energy, status )
!    call grvy_input_fread_int("vis/write_tot_energy", write_tot_energy, status )

    call variable_read(args,'write_rot_dof',write_rot_dof,inp_line_true)
    call variable_read(args,'write_vib_dof',write_vib_dof,inp_line_true)
!    call grvy_input_fread_int("vis/write_rot_dof", write_rot_dof, status )
!    call grvy_input_fread_int("vis/write_vib_dof", write_vib_dof, status )

    call variable_read(args,'write_entropy',write_entropy,inp_line_true)
    call variable_read(args,'write_heat_flux_x',write_heat_flux_x,inp_line_true)
    call variable_read(args,'write_heat_flux_y',write_heat_flux_y,inp_line_true)
    call variable_read(args,'write_heat_flux_z',write_heat_flux_z,inp_line_true)
!    call grvy_input_fread_int("vis/write_entropy", write_entropy, status )
!    call grvy_input_fread_int("vis/write_heat_flux_x", write_heat_flux_x, status )
!    call grvy_input_fread_int("vis/write_heat_flux_y", write_heat_flux_y, status )
!    call grvy_input_fread_int("vis/write_heat_flux_z", write_heat_flux_z, status )

    call variable_read(args,'write_kin_temp_relax',write_kin_temp_relax,inp_line_true)
    call variable_read(args,'write_rot_temp_relax',write_rot_temp_relax,inp_line_true)
    call variable_read(args,'write_vib_temp_relax',write_vib_temp_relax,inp_line_true)
!    call grvy_input_fread_int("vis/write_kin_temp_relax", write_kin_temp_relax, status )
!    call grvy_input_fread_int("vis/write_rot_temp_relax", write_rot_temp_relax, status )
!    call grvy_input_fread_int("vis/write_vib_temp_relax", write_vib_temp_relax, status )

    call set_properties_output_vars( write_dens, &
         write_x_vel, write_y_vel, write_z_vel, write_speed, write_mach, &
         write_tr_temp, write_rot_temp, write_vib_temp, write_tot_temp, &
         write_temp_x, write_temp_y, write_temp_z, write_pressure, &
         write_rot_dof, write_vib_dof, &
         write_tr_energy, write_rot_energy, write_vib_energy, write_tot_energy, &
         write_entropy, write_heat_flux_x, write_heat_flux_y, write_heat_flux_z, &
         write_kin_temp_relax, write_rot_temp_relax, write_vib_temp_relax )
    
    call variable_read(args,'write_stress',write_stress,inp_line_true)
!    call grvy_input_fread_int("vis/write_stress", write_stress, status )
    if( write_stress .gt. 0 )then
       allocate( stress_directions(1:write_stress), STAT=status )
       call allocate_error_check( status, "stress_directions (input)" )
       call variable_read_vec(args,'stress_directions',stress_directions,inp_line_true,write_stress)
!       call grvy_input_fread_int_vec( "vis/stress_directions", stress_directions, write_stress, status )
    end if
    call set_stress_output( write_stress, stress_directions )

    call variable_read(args,'write_moment',write_moment,inp_line_true)
!    call grvy_input_fread_int( 'vis/write_moment', write_moment, status )
    if( write_moment .gt. 0 )then
       allocate( moment_numbers(1:write_moment), STAT=status )
       call allocate_error_check( status, "moment_numbers (input)" )
       call variable_read_vec(args,'moment_numbers',moment_numbers,inp_line_true,write_moment)
!       call grvy_input_fread_int_vec( 'vis/moment_numbers', moment_numbers, write_moment, status )
    end if
    call set_moment_output( write_moment, moment_numbers )

    deallocate( x_output_locations, STAT=status )
    call deallocate_error_check( status, "x_output_locations (input)" )

    deallocate( y_output_locations, STAT=status )
    call deallocate_error_check( status, "y_output_locations (input)" )

    if( write_stress .gt. 0 )then
       deallocate( stress_directions, STAT=status )
       call deallocate_error_check( status, "stress_directions" )
    end if

    if( write_moment .gt. 0 )then
       deallocate( moment_numbers, STAT=status )
       call deallocate_error_check( status, "moment_numbers (input)" )
    end if

    !===============================================================================================
    ! Read restart file info
    !===============================================================================================
    call variable_read(args,'write_rst',write_rst,inp_line_true)
    call variable_read(args,'write_freq',write_freq,inp_line_true)
    call variable_read(args,'read_rst',read_rst,inp_line_true)
!    call grvy_input_fread_int( 'restart/write_rst', write_rst, status )
!    call grvy_input_fread_int( 'restart/write_freq', write_freq, status )
!    call grvy_input_fread_int( 'restart/read_rst', read_rst, status )
    call set_restart_flags( write_rst, write_freq, read_rst )
    
    call variable_read(args,'write_restart_filename',write_restart_filename,inp_line_true)
    call variable_read(args,'read_restart_filename',read_restart_filename,inp_line_true)
!    call grvy_input_fread_char( 'restart/write_restart_filename', write_restart_filename, status )
!    call grvy_input_fread_char( 'restart/read_restart_filename', read_restart_filename, status )
    call set_restart_filename( write_restart_filename, read_restart_filename )

    !===============================================================================================
    ! Read random number generation info
    !===============================================================================================
    call variable_read(args,'rand_method',rand_method,inp_line_true)
!    call grvy_input_fread_int( 'rng/rand_method', rand_method, status )
    call set_rand_method( rand_method )

    call variable_read(args,'seed',seed,inp_line_true)
!    call grvy_input_fread_int( 'rng/seed', seed, status )
    call set_seed( seed )

    deallocate(args,STAT=status)
    call deallocate_error_check(status,'args')

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

       read(vel_grid_file_unit,*) num_points_x
       read(vel_grid_file_unit,*) num_points_y
       read(vel_grid_file_unit,*) num_points_z
       
       call set_vel_grid_num_points( num_points_x, num_points_y, num_points_z, n )

       call initialize_file_read_arrays( n )

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
