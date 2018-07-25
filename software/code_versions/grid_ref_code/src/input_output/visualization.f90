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
module Visualization

  use DistFunc
  use SpeciesAndReferenceData
  use Constants
  use ErrorCheck
  use Conversion

  implicit none

  private

  type SpeedType
     integer :: num_speeds
     double precision, allocatable, dimension(:) :: speeds
     integer, allocatable, dimension(:) :: i
     integer, allocatable, dimension(:) :: j
     integer, allocatable, dimension(:) :: k
  end type SpeedType

  ! File units
  integer, parameter :: base_unit = 100
  integer :: df_vis_unit
  integer :: rot_levels_vis_unit
  integer :: vib_levels_vis_unit
  integer :: properties_vis_unit
  integer :: collisions_vis_unit
  integer :: bkw_analytic_vis_unit
  integer :: speed_vis_unit
  integer :: speed_rot_vis_unit, speed_vib_vis_unit

  ! File input variables
  logical :: enable_vel_df_vis
  logical :: enable_rot_df_vis
  logical :: enable_vib_df_vis
  logical :: enable_df_vis
  logical :: enable_rot_levels_vis
  logical :: enable_vib_levels_vis
  logical :: enable_properties_vis
  logical :: enable_collisions_vis
  logical :: enable_bkw_analytic_vis
  logical :: enable_speed_vis
  logical :: enable_speed_rot_vis
  logical :: enable_speed_vib_vis

  character(len=128) :: df_vis_filename
  character(len=128) :: rot_levels_vis_filename
  character(len=128) :: vib_levels_vis_filename
  character(len=128) :: properties_vis_filename
  character(len=128) :: collisions_vis_filename
  character(len=128) :: speed_vis_filename

  integer :: df_vis_dump_freq
  integer :: rot_levels_vis_dump_freq
  integer :: vib_levels_vis_dump_freq
  integer :: properties_vis_dump_freq
  integer :: collisions_vis_dump_freq
  integer :: speed_vis_dump_freq

  integer :: num_output_locations
  integer :: df_output_dimension
  integer, allocatable, dimension(:) :: output_locations

  ! Output properites
  integer, parameter :: number_properties = 28
  logical, dimension(1:number_properties) :: properties_to_output
  character(len=128), dimension(1:number_properties) :: property_names
      
  integer ::  write_moment
  integer, allocatable, dimension(:) :: moment_numbers

  double precision, allocatable, dimension(:,:) :: total_data, total_moments
  double precision, allocatable, dimension(:,:,:) :: species_data, species_moments

  ! Analytic BKW
  type(DistFuncType), dimension(1,1) :: bkw

  ! Speed DF
  type(SpeedType), allocatable, dimension(:) :: sorted_speed

  ! Dimension parameters
  integer, parameter :: oneD   = 1
  integer, parameter :: twoD   = 2
  integer, parameter :: threeD = 3

  ! Direction parameters
  integer, parameter :: x_dir = 1
  integer, parameter :: y_dir = 2
  integer, parameter :: z_dir = 3

  public :: enable_vis
  public :: enable_bkw_analytic
  public :: set_vis_filename
  public :: set_vis_dump_freq
  public :: set_output_locations
  public :: set_df_output_dimension
  public :: set_properties_output_vars
  public :: set_moment_output
  public :: initialize_visualization
  public :: destroy_visualization
  public :: write_vis_data

contains

  !===================================================================================================
  ! Setter Subroutines
  !===================================================================================================

  subroutine enable_vis( enable_vel_df_vis_in, enable_rot_df_vis_in, enable_vib_df_vis_in, &
       enable_rot_levels_vis_in, enable_vib_levels_vis_in, enable_properties_vis_in, &
       enable_collisions_vis_in, &
       enable_speed_vis_in, enable_speed_rot_vis_in, enable_speed_vib_vis_in )

    implicit none

    integer, intent(in) :: enable_vel_df_vis_in, enable_rot_df_vis_in, enable_vib_df_vis_in
    integer, intent(in) :: enable_rot_levels_vis_in, enable_vib_levels_vis_in
    integer, intent(in) :: enable_properties_vis_in, enable_collisions_vis_in
    integer, intent(in) :: enable_speed_vis_in, enable_speed_rot_vis_in, enable_speed_vib_vis_in

    enable_vel_df_vis = int2logical( enable_vel_df_vis_in )
    enable_rot_df_vis = int2logical( enable_rot_df_vis_in )
    enable_vib_df_vis = int2logical( enable_vib_df_vis_in )

    if( enable_vel_df_vis .or. enable_rot_df_vis .or. enable_vib_df_vis ) enable_df_vis = .true.

    enable_rot_levels_vis = int2logical( enable_rot_levels_vis_in )
    enable_vib_levels_vis = int2logical( enable_vib_levels_vis_in )
    enable_properties_vis = int2logical( enable_properties_vis_in )
    enable_collisions_vis = int2logical( enable_collisions_vis_in )

    enable_speed_vis     = int2logical( enable_speed_vis_in )
    enable_speed_rot_vis = int2logical( enable_speed_rot_vis_in )
    enable_speed_vib_vis = int2logical( enable_speed_vib_vis_in )

    return
  end subroutine enable_vis

  subroutine enable_bkw_analytic( enable_bkw_analytic_vis_in )
    
    implicit none
    
    integer, intent(in) :: enable_bkw_analytic_vis_in

    enable_bkw_analytic_vis = int2logical( enable_bkw_analytic_vis_in )

    return
  end subroutine enable_bkw_analytic

  subroutine set_vis_filename( df_vis_filename_in, rot_levels_vis_filename_in, vib_levels_vis_filename_in, &
       properties_vis_filename_in, collisions_vis_filename_in, speed_vis_filename_in )

    implicit none

    character(len=128), intent(in) :: df_vis_filename_in
    character(len=128), intent(in) :: rot_levels_vis_filename_in, vib_levels_vis_filename_in
    character(len=128), intent(in) :: properties_vis_filename_in, collisions_vis_filename_in
    character(len=128), intent(in) :: speed_vis_filename_in

    df_vis_filename         = df_vis_filename_in
    rot_levels_vis_filename = rot_levels_vis_filename_in
    vib_levels_vis_filename = vib_levels_vis_filename_in
    properties_vis_filename = properties_vis_filename_in
    collisions_vis_filename = collisions_vis_filename_in
    speed_vis_filename      = speed_vis_filename_in

    return
  end subroutine set_vis_filename
  
  subroutine set_vis_dump_freq( df_vis_dump_freq_in, rot_levels_vis_dump_freq_in, &
       vib_levels_vis_dump_freq_in, properties_vis_dump_freq_in, collisions_vis_dump_freq_in, &
       speed_vis_dump_freq_in )

    implicit none

    integer, intent(in) :: df_vis_dump_freq_in
    integer, intent(in) :: rot_levels_vis_dump_freq_in, vib_levels_vis_dump_freq_in
    integer, intent(in) :: properties_vis_dump_freq_in, collisions_vis_dump_freq_in
    integer, intent(in) :: speed_vis_dump_freq_in

    df_vis_dump_freq         = df_vis_dump_freq_in
    rot_levels_vis_dump_freq = rot_levels_vis_dump_freq_in
    vib_levels_vis_dump_freq = vib_levels_vis_dump_freq_in
    properties_vis_dump_freq = properties_vis_dump_freq_in
    collisions_vis_dump_freq = collisions_vis_dump_freq_in
    speed_vis_dump_freq      = speed_vis_dump_freq_in

    return
  end subroutine set_vis_dump_freq

  subroutine set_output_locations( num_output_locations_in, output_locations_in )

    implicit none

    integer, intent(in) :: num_output_locations_in
    integer, dimension(:), intent(in) :: output_locations_in

    integer :: status

    num_output_locations = num_output_locations_in

    allocate( output_locations(1:num_output_locations), STAT=status )
    call allocate_error_check( status, "output_locations" )

    output_locations = output_locations_in

    return
  end subroutine set_output_locations

  subroutine set_df_output_dimension( df_output_dimension_in )
    
    implicit none

    integer, intent(in) :: df_output_dimension_in

    df_output_dimension = df_output_dimension_in

    return
  end subroutine set_df_output_dimension

  subroutine set_properties_output_vars( write_dens_in, &
       write_x_vel_in, write_y_vel_in, write_z_vel_in, write_speed_in, write_mach_in, &
       write_tr_temp_in, write_rot_temp_in, write_vib_temp_in, write_tot_temp_in, &
       write_temp_x_in, write_temp_y_in, write_temp_z_in, write_pressure_in, &
       write_rot_dof_in, write_vib_dof_in, &
       write_tr_energy_in, write_rot_energy_in, write_vib_energy_in, write_tot_energy_in, &
       write_entropy_in, write_heat_flux_x_in, write_heat_flux_y_in, write_heat_flux_z_in, &
       write_shear_stress_in, &
       write_kin_temp_relax_in, write_rot_temp_relax_in, write_vib_temp_relax_in )

    implicit none

    integer, intent(in) :: write_dens_in, write_x_vel_in, write_y_vel_in, write_z_vel_in
    integer, intent(in) :: write_speed_in, write_mach_in
    integer, intent(in) :: write_tr_temp_in, write_rot_temp_in, write_vib_temp_in
    integer, intent(in) :: write_tot_temp_in
    integer, intent(in) :: write_temp_x_in, write_temp_y_in, write_temp_z_in
    integer, intent(in) :: write_pressure_in
    integer, intent(in) :: write_rot_dof_in, write_vib_dof_in
    integer, intent(in) :: write_tot_energy_in, write_rot_energy_in, write_vib_energy_in, write_tr_energy_in
    integer, intent(in) :: write_entropy_in, write_heat_flux_x_in, write_heat_flux_y_in, write_heat_flux_z_in
    integer, intent(in) :: write_shear_stress_in
    integer, intent(in) :: write_kin_temp_relax_in, write_rot_temp_relax_in, write_vib_temp_relax_in


    properties_to_output(1)  = int2logical( write_dens_in )

    properties_to_output(2)  = int2logical( write_x_vel_in )
    properties_to_output(3)  = int2logical( write_y_vel_in )
    properties_to_output(4)  = int2logical( write_z_vel_in )
    properties_to_output(5)  = int2logical( write_speed_in )
    properties_to_output(6)  = int2logical( write_mach_in )
    
    properties_to_output(7)  = int2logical( write_tr_temp_in )
    properties_to_output(8)  = int2logical( write_rot_temp_in )
    properties_to_output(9)  = int2logical( write_vib_temp_in )
    properties_to_output(10) = int2logical( write_tot_temp_in )

    properties_to_output(11) = int2logical( write_temp_x_in ) 
    properties_to_output(12) = int2logical( write_temp_y_in )
    properties_to_output(13) = int2logical( write_temp_z_in )

    properties_to_output(14) = int2logical( write_tr_energy_in )
    properties_to_output(15) = int2logical( write_rot_energy_in )
    properties_to_output(16) = int2logical( write_vib_energy_in )
    properties_to_output(17) = int2logical( write_tot_energy_in )

    properties_to_output(18) = int2logical( write_entropy_in )
    properties_to_output(19) = int2logical( write_heat_flux_x_in )
    properties_to_output(20) = int2logical( write_heat_flux_y_in )
    properties_to_output(21) = int2logical( write_heat_flux_z_in )
    properties_to_output(22) = int2logical( write_shear_stress_in )

    properties_to_output(23) = int2logical( write_rot_dof_in )
    properties_to_output(24) = int2logical( write_vib_dof_in )

    ! These variable are for analytic solutions to relaxations
    properties_to_output(27) = int2logical( write_kin_temp_relax_in )
    properties_to_output(25) = int2logical( write_rot_temp_relax_in )
    properties_to_output(26) = int2logical( write_vib_temp_relax_in )

    properties_to_output(28) = int2logical( write_pressure_in )
    
    
    property_names(1)  = '"Density"'
    property_names(2)  = '"X-Velocity"'
    property_names(3)  = '"Y-Velocity"'
    property_names(4)  = '"Z-Velocity"'
    property_names(5)  = '"Speed"'
    property_names(6)  = '"Mach Number"'
    property_names(7)  = '"Trans. Temp."'
    property_names(8)  = '"Rot. Temp."'
    property_names(9)  = '"Vib. Temp."'
    property_names(10) = '"Total Temp."'
    property_names(11) = '"Temp. X"'
    property_names(12) = '"Temp. Y"'
    property_names(13) = '"Temp. Z"'
    property_names(14) = '"Trans. Energy"'
    property_names(15) = '"Rot. Energy"'
    property_names(16) = '"Vib. Energy"'
    property_names(17) = '"Total Energy"'
    property_names(18) = '"Entropy"'
    property_names(19) = '"Heat Flux X"'
    property_names(20) = '"Heat Flux Y"'
    property_names(21) = '"Heat Flux Z"'
    property_names(22) = '"Shear Stress"'
    property_names(23) = '"Rot DOF"'
    property_names(24) = '"Vib DOF"'
    property_names(27) = '"Trans. Temp Analytic"'
    property_names(25) = '"Rot. Temp Analytic"'
    property_names(26) = '"Vib. Temp Analytic"'
    property_names(28) = '"Pressure"'

    return
  end subroutine set_properties_output_vars

  subroutine set_moment_output( write_moment_in, moment_numbers_in )

    implicit none

    integer, intent(in) :: write_moment_in
    integer, dimension(:), intent(in) :: moment_numbers_in

    integer :: status

    write_moment = write_moment_in

    if( write_moment .gt. 0 )then
       allocate( moment_numbers(1:write_moment), STAT=status )
       call allocate_error_check( status, "moment_numbers" )

       moment_numbers = moment_numbers_in
    end if

    return
  end subroutine set_moment_output

  !===================================================================================================
  ! Initialize and Destroy
  !===================================================================================================
  subroutine initialize_visualization( vel_grid, molecule, end_time )

    use VelocityGrid
    use PhysicalGrid
    use SpeciesAndReferenceData

    implicit none

    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    type(MoleculeType), dimension(:), intent(in) :: molecule
    integer, intent(in) :: end_time
    integer :: status, nx_space, ny_space

    call get_nspace( nx_space, ny_space )
    
    if( enable_bkw_analytic_vis )then

       if( num_species .gt. 1 .or. nx_space .gt. 1 )then
          write(*,*) 'Error: Analytic solution to bkw function only available for single species, 0D.'
          stop
       end if

       call create_dist_func( bkw(1,1), vel_grid(1,1), molecule(1), num_species )

    end if

    ! Properties data structures
    allocate( total_data( 1:nx_space, 1:number_properties ), STAT=status )
    call allocate_error_check( status, "total_data" )

    allocate( species_data( 1:num_species, 1:nx_space, 1:number_properties ), STAT=status )
    call allocate_error_check( status, "species_data" )

    if( write_moment .gt. 0 )then
       allocate( total_moments( 1:nx_space, 1:write_moment ), STAT=status )
       call allocate_error_check( status, "total_moments" )

       allocate( species_moments( 1:num_species, 1:nx_space, 1:write_moment ), STAT=status )
       call allocate_error_check( status, "species_moments" )
    end if

    if( enable_speed_vis .or. enable_speed_rot_vis .or. enable_speed_vib_vis )then
       allocate( sorted_speed( 1:num_species ), STAT=status )
       call allocate_error_check( status, "sorted_speed" )

       call sort_speed( vel_grid(:,1) )
    end if

    df_vis_unit           = base_unit
    rot_levels_vis_unit   = df_vis_unit + num_species*num_output_locations
    vib_levels_vis_unit   = rot_levels_vis_unit + num_species*num_output_locations
    properties_vis_unit   = vib_levels_vis_unit + num_species*num_output_locations
    collisions_vis_unit   = properties_vis_unit + num_species + 1
    speed_vis_unit        = collisions_vis_unit + 1
    speed_rot_vis_unit    = speed_vis_unit + num_species*num_output_locations
    speed_vib_vis_unit    = speed_rot_vis_unit + num_species*num_output_locations
    bkw_analytic_vis_unit = speed_vib_vis_unit + num_species*num_output_locations

    call write_headers(end_time)

    return
  end subroutine initialize_visualization

  subroutine destroy_visualization()

    implicit none

    integer :: status

    ! Properties data structures
    deallocate( total_data, STAT=status )
    call deallocate_error_check( status, "total_data" )

    deallocate( species_data, STAT=status )
    call deallocate_error_check( status, "species_data" )

    if( enable_speed_vis .or. enable_speed_rot_vis .or. enable_speed_vib_vis )then
       deallocate( sorted_speed, STAT=status )
       call deallocate_error_check( status, "sorted_speed" )
    end if

    return
  end subroutine destroy_visualization

  subroutine write_headers( end_time )

    use PhysicalGrid

    implicit none

    integer, intent(in) :: end_time
    character(len=128) :: filename
    integer :: file_unit, n, i, loc, level
    integer :: rot_levels, vib_levels
    integer :: nx_space, ny_space
    integer :: print_times

    call get_nspace( nx_space, ny_space )

    ! Distribution function output header
    if( enable_df_vis )then

       do i = 1, num_species

          call get_num_energy_levels( rot_levels, vib_levels, i )

          do loc = 1, num_output_locations

             file_unit = df_vis_unit + (i-1)*num_output_locations + (loc-1)
             filename = trim(adjustl(df_vis_filename))//"_"//&
                  trim(adjustl(int2char(i)))//"_"//&
                  trim(adjustl(int2char(output_locations(loc))))//".dat"

             open( unit=file_unit, file=filename, status="unknown" )

             write( file_unit, advance="yes", fmt="(a)" ) 'title = "Distribution Function"'

             select case( df_output_dimension )
             case( oneD )
                write( file_unit, advance="no", fmt="(a)" )&
                     'variables = "X-Velocity", "Beta3", "Time"'

             case( twoD )
                write( file_unit, advance="no", fmt="(a)" )&
                     'variables = "X-Velocity", "Y-Velocity", "Beta3", "Time"'

             case( threeD )
                write( file_unit, advance="no", fmt="(a)" )&
                     'variables = "X-Velocity", "Y-Velocity", "Z-Velocity", "Beta3", "Time"'

             case default
                write(*,*) "Error: invalid df_output_dimension: ", df_output_dimension
                stop

             end select

             if( enable_vel_df_vis )then
                write( file_unit, advance="no", fmt="(a)" ) ', "Vel DistFunc", "Deviation"'
             end if

             if( enable_rot_df_vis )then
                do level = 1, rot_levels
                   write( file_unit, advance="no", fmt="(a7,i3,a1)" ) ', "Rot ', level,'"'
                end do
                write( file_unit, advance="no", fmt="(a)" ) ', "Rot Energy"'
             end if

             if( enable_vib_df_vis )then
                do level = 1, vib_levels
                   write( file_unit, advance="no", fmt="(a7,i3,a1)" ) ', "Vib ', level,'"'
                end do
                write( file_unit, advance="no", fmt="(a)" ) ', "Vib Energy"'
             end if

             write( file_unit, advance="yes", fmt="(1x)" )


             close( file_unit )

          end do
       end do

    end if

    if( enable_rot_levels_vis )then

       do i = 1, num_species
          do loc = 1, num_output_locations

             file_unit = rot_levels_vis_unit + (i-1)*num_output_locations + (loc-1)
             filename = trim(adjustl(rot_levels_vis_filename))//"_"//&
                  trim(adjustl(int2char(i)))//"_"//&
                  trim(adjustl(int2char(output_locations(loc))))//".dat"

             open( unit=file_unit, file=filename, status="unknown" )

             write( file_unit, advance="yes", fmt="(a)" ) 'title = "Rotational Levels"'
             write( file_unit, advance="yes", fmt="(a)" )&
                  'variables = "Time", "Level Number", "Level Energy", "Density Fraction", "Energy", &
                  "Spin", "Degeneracy"'

             close( file_unit )

             !>>>>>>>>>>
             open( unit=77+loc, file="vel_dependent_rot_levels_"//trim(adjustl(int2char(loc)))//".dat", &
                  status="unknown" )
             write( 77+loc, advance="yes", fmt="(a)" ) 'title = "Vel Rotational Levels"'
             write( 77+loc, advance="yes", fmt="(a)" )&
                  'variables = "Time", "Level Number", "Level Energy", "Density Fraction", "Energy", &
                  "Spin", "Degeneracy"'
             close( 77+loc )
             !<<<<<<<<<<

          end do
       end do

    end if

    if( enable_vib_levels_vis )then

       do i = 1, num_species
          do loc = 1, num_output_locations

             file_unit = vib_levels_vis_unit + (i-1)*num_output_locations + (loc-1)
             filename = trim(adjustl(vib_levels_vis_filename))//"_"//&
                  trim(adjustl(int2char(i)))//"_"//&
                  trim(adjustl(int2char(output_locations(loc))))//".dat"

             open( unit=file_unit, file=filename, status="unknown" )

             write( file_unit, advance="yes", fmt="(a)" ) 'title = "Vibational Levels"'
             write( file_unit, advance="yes", fmt="(a)" )&
                  'variables = "Time", "Level Number", "Level Energy", "Density Fraction", "Energy"'

             close( file_unit )

          end do
       end do

    end if

    if( enable_properties_vis )then

       file_unit = properties_vis_unit
       filename = trim(adjustl(properties_vis_filename))//".dat"

       open( unit=file_unit, file=filename, status="unknown" )

       write( file_unit, advance="yes", fmt="(a)" ) 'title = "Macroscopic Properties"'
       write( file_unit, advance="no", fmt="(a)" ) 'variables = "Time", "Location"'

       do i = 1, number_properties
          if( properties_to_output(i) )&
               write( file_unit, advance="no", fmt="(a3,a14)" ) ", ",property_names(i)
       end do

       if( write_moment .gt. 0 )then
          do i = 1, write_moment
             write( file_unit, advance="no", fmt="(a10,i2,a1)" ) ', "Moment ',moment_numbers(i),'"'
          end do
       end if

       write( file_unit, advance="yes", fmt="(1x)" )

       if( nx_space .eq. 1 )then
          if( mod(end_time,properties_vis_dump_freq) .eq. 0 )then
             print_times = end_time/properties_vis_dump_freq + 1
          else
             print_times = end_time/properties_vis_dump_freq + 2
          end if
          write( file_unit, advance="yes", fmt="(a7,i4)" ) "zone i=", print_times
       end if

       close( file_unit )

       do n = 1, num_species

          file_unit = properties_vis_unit + n
          filename = trim(adjustl(properties_vis_filename))//"_"//trim(adjustl(int2char(n)))//".dat"

          open( unit=file_unit, file=filename, status="unknown" )

          write( file_unit, advance="yes", fmt="(a)" ) 'title = "Macroscopic Properties"'
          write( file_unit, advance="no", fmt="(a)" ) 'variables = "Time", "Location"'

          do i = 1, number_properties
             if( properties_to_output(i) )&
                  write( file_unit, advance="no", fmt="(a2,a)" ) ", ",trim(adjustl(property_names(i)))
          end do

          if( write_moment .gt. 0 )then
             do i = 1, write_moment
                write( file_unit, advance="no", fmt="(a10,i2,a1)" ) ', "Moment ',moment_numbers(i),'"'
             end do
          end if

          write( file_unit, advance="yes", fmt="(1x)" )

          if( nx_space .eq. 1 )then
             if( mod(end_time,properties_vis_dump_freq) .eq. 0 )then
                print_times = end_time/properties_vis_dump_freq + 1
             else
                print_times = end_time/properties_vis_dump_freq + 2
             end if
             write( file_unit, advance="yes", fmt="(a7,i4)" ) "zone i=", print_times
          end if

          close( file_unit )

       end do

    end if

    if( enable_collisions_vis )then

       file_unit = collisions_vis_unit
       filename = trim(adjustl(collisions_vis_filename))//".dat"

       open( unit=file_unit, file=filename, status="unknown" )

       write( file_unit, advance="yes", fmt="(a)" ) 'title = "Number of Collisions"'
       write( file_unit, advance="no", fmt="(a)" ) 'variables = "Time", "Location"'

       do n = 1, num_species
          do i = n, num_species
             write( file_unit, advance="no", fmt="(a13,i2,a1,i2,a1)" )', "Coll Type ',n,'-',i,'"'
          end do
       end do

       write( file_unit, advance="yes", fmt="(a)" )', "Total Collisions"'

       if( nx_space .eq. 1 )then
          if( mod(end_time,collisions_vis_dump_freq) .eq. 0 )then
             print_times = end_time/collisions_vis_dump_freq + 1
          else
             print_times = end_time/collisions_vis_dump_freq + 2
          end if
          write( file_unit, advance="yes", fmt="(a7,i4)" ) "zone i=", print_times
       end if

       close( file_unit )

    end if

    if( enable_bkw_analytic_vis .and. enable_properties_vis )then

       file_unit = bkw_analytic_vis_unit
       filename = trim(adjustl(properties_vis_filename))//"_bkw.dat"

       open( unit=file_unit, file=filename, status="unknown" )

       open( unit=file_unit, file=filename, status="unknown" )

       write( file_unit, advance="yes", fmt="(a)" ) 'title = "Macroscopic Properties"'
       write( file_unit, advance="no", fmt="(a)" ) 'variables = "Time", "Location"'

       do i = 1, number_properties
          if( properties_to_output(i) )&
               write( file_unit, advance="no", fmt="(a3,a14)" ) ", ",property_names(i)
       end do

       if( write_moment .gt. 0 )then
          do i = 1, write_moment
             write( file_unit, advance="no", fmt="(a10,i2,a1)" ) ', "Moment ',moment_numbers(i),'"'
          end do
       end if

       write( file_unit, advance="yes", fmt="(1x)" )
       
       if( nx_space .eq. 1 )then
          if( mod(end_time,properties_vis_dump_freq) .eq. 0 )then
             print_times = end_time/properties_vis_dump_freq + 1
          else
             print_times = end_time/properties_vis_dump_freq + 2
          end if
          write( file_unit, advance="yes", fmt="(a7,i4)" ) "zone i=", print_times
       end if

       close( file_unit )

    end if

    if( enable_speed_vis )then
       do i = 1, num_species
          do loc = 1, num_output_locations

             file_unit = speed_vis_unit + (i-1)*num_output_locations + (loc-1)
             filename = trim(adjustl(speed_vis_filename))//"_"//&
                  trim(adjustl(int2char(i)))//"_"//&
                  trim(adjustl(int2char(output_locations(loc))))//".dat"

             open( unit=file_unit, file=filename, status="unknown" )

             write( file_unit, advance="yes", fmt="(a)" ) 'title = "Speed Distribution"'
             write( file_unit, advance="yes", fmt="(a)" ) 'variables = "Speed", "DF", "Deviation DF"'

          end do
       end do

    end if

    if( enable_speed_rot_vis )then
       do i = 1, num_species
          do loc = 1, num_output_locations

             file_unit = speed_rot_vis_unit + (i-1)*num_output_locations + (loc-1)
             filename = trim(adjustl(speed_vis_filename))//"_"//&
                  trim(adjustl(int2char(i)))//"_"//&
                  trim(adjustl(int2char(output_locations(loc))))//"_rot.dat"

             open( unit=file_unit, file=filename, status="unknown" )

             write( file_unit, advance="yes", fmt="(a)" ) 'title = "Speed Distribution with Rotational Levels"'
             write( file_unit, advance="yes", fmt="(a)" )&
                  'variables = "Speed", "Level", "Rot DF", "Normalized Rot DF", "Vel DF", "Deviation DF", &
                  "Total Deviation"'

          end do
       end do

    end if

    if( enable_speed_vib_vis )then
       do i = 1, num_species
          do loc = 1, num_output_locations

             file_unit = speed_vib_vis_unit + (i-1)*num_output_locations + (loc-1)
             filename = trim(adjustl(speed_vis_filename))//"_"//&
                  trim(adjustl(int2char(i)))//"_"//&
                  trim(adjustl(int2char(output_locations(loc))))//"_vib.dat"

             open( unit=file_unit, file=filename, status="unknown" )

             write( file_unit, advance="yes", fmt="(a)" ) 'title = "Speed Distribution with Vibrational Levels"'
             write( file_unit, advance="yes", fmt="(a)" )&
                  'variables = "Speed", "Level", "Vib DF", "Normalized Vib DF", "Vel DF", "Deviation DF", &
                  "Total Deviation"'

          end do
       end do

    end if

    return
  end subroutine write_headers

  !===================================================================================================
  !  Output Data
  !===================================================================================================
  subroutine write_vis_data( time, ntime, end_time, phi, molecule, vel_grid, properties, num_collisions )

    use VelocityGrid
    use PhysicalProperties
    use InitialConditions

    implicit none

    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(VelocityGridType),dimension(:,:), intent(in) :: vel_grid
    type(PropertiesType), dimension(:), intent(in) :: properties
    double precision, intent(in) :: time
    integer, intent(in) :: ntime, end_time
    integer, dimension(:,:,:), intent(in) :: num_collisions

    double precision :: dens, temp, x_vel, mass
    
    if( enable_df_vis )then
       if( mod(ntime,end_time) .eq. 0 .or. mod(ntime,df_vis_dump_freq) .eq. 0 )&
            call write_df( time, phi, molecule, properties, vel_grid )
!!$       if( ntime .gt. 65 )&
!!$            call write_df( time, phi, molecule, properties, vel_grid )
    end if

    if( enable_rot_levels_vis )then
       if( mod(ntime,end_time) .eq. 0 .or. mod(ntime,rot_levels_vis_dump_freq) .eq. 0  )&
            call write_rot_levels( time, phi, molecule, properties, vel_grid )
    end if

    if( enable_vib_levels_vis )then
       if( mod(ntime,end_time) .eq. 0 .or. mod(ntime,vib_levels_vis_dump_freq) .eq. 0  )&
            call write_vib_levels( time, phi, molecule, properties, vel_grid )
    end if

    if( enable_properties_vis )then
       if( mod(ntime,end_time) .eq. 0 .or. mod(ntime,properties_vis_dump_freq) .eq. 0  )&
            call write_properties( time, phi, molecule, vel_grid, properties )
    end if

    if( enable_collisions_vis )then
       if( mod(ntime,end_time) .eq. 0 .or. mod(ntime,collisions_vis_dump_freq) .eq. 0  )&
            call write_collisions( time, num_collisions )
    end if

    if( enable_bkw_analytic_vis .and. enable_properties_vis )then
       if( mod(ntime,end_time) .eq. 0 .or. mod(ntime,properties_vis_dump_freq) .eq. 0 )then

          dens  = properties(1)%dens(1)
          temp  = properties(1)%tr_temp(1)
          x_vel = properties(1)%x_vel(1)
          mass  = molecule(1)%mass

          call generate_bkw( bkw(1,1), mass, dens, temp, x_vel, vel_grid(1,1), time )

          call write_bkw_properties( time, molecule, vel_grid, properties )

       end if
    end if

    if( enable_speed_vis )then
       if( mod(ntime,end_time) .eq. 0 .or. mod(ntime,speed_vis_dump_freq) .eq. 0 )&   
            call write_speed_df( phi, vel_grid, properties, molecule )
    end if

    if( enable_speed_rot_vis )then
       if( mod(ntime,end_time) .eq. 0 .or. mod(ntime,speed_vis_dump_freq) .eq. 0 )&
            call write_speed_rot_df( phi, vel_grid, properties, molecule )
!!$       if( ntime .gt. 65 )&
!!$            call write_speed_rot_df( phi, vel_grid, properties, molecule )
    end if

    if( enable_speed_vib_vis )then
       if( mod(ntime,end_time) .eq. 0 .or. mod(ntime,speed_vis_dump_freq) .eq. 0 )&
            call write_speed_vib_df( phi, vel_grid, properties, molecule )
!!$       if( ntime .gt. 145 )&
!!$            call write_speed_vib_df( phi, vel_grid, molecule )
    end if

    return
  end subroutine write_vis_data

  subroutine write_df( time, phi, molecule, properties, vel_grid )

    use VelocityGrid
    use PhysicalGrid
    use PhysicalProperties

    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(PropertiesType), dimension(:), intent(in) :: properties
    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    double precision, intent(in) :: time

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: r_modes, v_modes, r_levels, v_levels
    integer :: num_points_x, num_points_y, num_points_z
    integer :: i_zero, j_zero, k_zero
    integer :: n, l, loc, i, j, k
    integer :: grid_ref, dfloc

    double precision :: x, y, z
    double precision :: beta_x, beta_y, beta_z, beta3
    double precision :: energy
    double precision :: mass, dens, temp, u, v, w, coeff, df_eq

    integer :: file_unit
    character(len=128) :: filename
    
    do n = 1, num_species

       r_modes = molecule(n)%rot_modes
       v_modes = molecule(n)%vib_modes

       do loc = 1, num_output_locations

          dfloc = output_locations(loc)

          file_unit = df_vis_unit + (n-1)*num_output_locations + (loc-1)
          filename = trim(adjustl(df_vis_filename))//"_"//&
               trim(adjustl(int2char(n)))//"_"//&
               trim(adjustl(int2char(output_locations(loc))))//".dat"

          open( unit=file_unit, file=filename, status="old", position="append" )

          call get_spatial_reference( output_locations(loc), grid_ref )
          
          ! Equilibrium calculation
          mass  = molecule(n)%mass
          dens  = properties(dfloc)%dens(n)
          u     = properties(dfloc)%x_vel(n)
          v     = properties(dfloc)%y_vel(n)
          w     = properties(dfloc)%z_vel(n)
          temp  = properties(dfloc)%temp(n)
          coeff = dens * sqrt( mass * mass * mass ) / sqrt( pi * pi * pi * temp * temp * temp )

          i_min = vel_grid(n,grid_ref)%i_min
          i_max = vel_grid(n,grid_ref)%i_max
          j_min = vel_grid(n,grid_ref)%j_min
          j_max = vel_grid(n,grid_ref)%j_max
          k_min = vel_grid(n,grid_ref)%k_min
          k_max = vel_grid(n,grid_ref)%k_max

          r_levels = phi(n,dfloc)%num_rot_levels
          v_levels = phi(n,dfloc)%num_vib_levels

          num_points_x = vel_grid(n,grid_ref)%num_points_x
          num_points_y = vel_grid(n,grid_ref)%num_points_y
          num_points_z = vel_grid(n,grid_ref)%num_points_z

          do i = i_min, i_max
             x = vel_grid(n,grid_ref)%x(i)
             if( abs( x ) .lt. double_tol )then
                i_zero = i
                exit
             end if
          end do

          do j = j_min, j_max
             y = vel_grid(n,grid_ref)%y(j)
             if( abs( y ) .lt. double_tol )then
                j_zero = j
                exit
             end if
          end do

          do k = k_min, k_max
             z = vel_grid(n,grid_ref)%z(k)
             if( abs( z ) .lt. double_tol )then
                k_zero = k
                exit
             end if
          end do

          select case( df_output_dimension )
          case( oneD )
             write( file_unit, advance="yes", fmt="(a7,i4)" )&
                  "zone i=",num_points_x
             j = j_zero
             k = k_zero
             beta_y = vel_grid(n,grid_ref)%beta_y(j)
             beta_z = vel_grid(n,grid_ref)%beta_z(k)
             do i = i_min, i_max
                x = vel_grid(n,grid_ref)%x(i)
                beta_x = vel_grid(n,grid_ref)%beta_x(i)
                beta3 = beta_x*beta_y*beta_z
                call compute_maxwellian( x, zero, zero, coeff, mass, u, v, w, temp, df_eq )
                df_eq = df_eq * beta3

                write( file_unit, advance="no", fmt="(3e24.16)" ) x, beta3, time

                if( enable_vel_df_vis )then
                   write( file_unit, advance="no", fmt="(2e24.16)" ) &
                        phi(n,dfloc)%value(i,j,k), phi(n,dfloc)%value(i,j,k) - df_eq
                end if

                if( enable_rot_df_vis )then
                   energy = zero
                   do l = 1, r_levels
                      if( r_modes .gt. 0 )then
                         write( file_unit, advance="no", fmt="(e24.16)" )&
                              phi(n,dfloc)%rot(l,i,j,k)
                         energy = energy + phi(n,dfloc)%rot(l,i,j,k)*phi(n,dfloc)%rot_level(l)
                      else
                         write( file_unit, advance="no", fmt="(e24.16)" ) zero
                      end if
                   end do
                   write( file_unit, advance="no", fmt="(e24.16)" ) energy
                end if

                if( enable_vib_df_vis )then
                   energy = zero
                   do l = 1, v_levels
                      if( v_modes .gt. 0 )then
                         write( file_unit, advance="no", fmt="(e24.16)" )&
                              phi(n,dfloc)%vib(l,i,j,k)
                         energy = energy + phi(n,dfloc)%vib(l,i,j,k)*phi(n,dfloc)%vib_level(l)
                      else
                         write( file_unit, advance="no", fmt="(e24.16)" ) zero
                      end if
                   end do
                   write( file_unit, advance="no", fmt="(e24.16)" ) energy
                end if

                write( file_unit, advance="yes", fmt="(1x)" )

             end do

          case( twoD )
             write( file_unit, advance="yes", fmt="(a7,i3,a4,i3)" )&
                  "zone i=",num_points_x,", j=",num_points_y

             k = k_zero
             beta_z = vel_grid(n,grid_ref)%beta_z(k)
             do j = j_min, j_max
                y = vel_grid(n,grid_ref)%y(j)
                beta_y = vel_grid(n,grid_ref)%beta_y(j)
                do i = i_min, i_max
                   x = vel_grid(n,grid_ref)%x(i)
                   beta_x = vel_grid(n,grid_ref)%beta_x(i)
                   beta3 = beta_x*beta_y*beta_z
                   call compute_maxwellian( x, y, zero, coeff, mass, u, v, w, temp, df_eq )
                   df_eq = df_eq * beta3

                   write( file_unit, advance="no", fmt="(4e24.16)" ) x, y, beta3, time

                   if( enable_vel_df_vis )then
                      write( file_unit, advance="no", fmt="(2e24.16)" ) &
                           phi(n,dfloc)%value(i,j,k), phi(n,dfloc)%value(i,j,k) - df_eq
                   end if

                   if( enable_rot_df_vis )then
                      energy = zero
                      do l = 1, r_levels
                         if( r_modes .gt. 0 )then
                            write( file_unit, advance="no", fmt="(e24.16)" )&
                                 phi(n,dfloc)%rot(l,i,j,k)
                            energy = energy + phi(n,dfloc)%rot(l,i,j,k)*phi(n,dfloc)%rot_level(l)
                         else
                            write( file_unit, advance="no", fmt="(e24.16)" ) zero
                         end if
                      end do
                      write( file_unit, advance="no", fmt="(e24.16)" ) energy
                   end if

                   if( enable_vib_df_vis )then
                      energy = zero
                      do l = 1, v_levels
                         if( v_modes .gt. 0 )then
                            write( file_unit, advance="no", fmt="(e24.16)" )&
                                 phi(n,dfloc)%vib(l,i,j,k)
                            energy = energy + phi(n,dfloc)%vib(l,i,j,k)*phi(n,dfloc)%vib_level(l)
                         else
                            write( file_unit, advance="no", fmt="(e24.16)" ) zero
                         end if
                      end do
                      write( file_unit, advance="no", fmt="(e24.16)" ) energy
                   end if

                   write( file_unit, advance="yes", fmt="(1x)" )

                end do
             end do

          case( threeD )
             write( file_unit, advance="yes", fmt="(a7,i4,a9,i4,a9,i4)" )&
                  "zone i=",num_points_x,", j=",num_points_y,", k=",num_points_z

             do k = k_min, k_max
                z = vel_grid(n,grid_ref)%z(k)
                beta_z = vel_grid(n,grid_ref)%beta_z(k)
                do j = j_min, j_max
                   y = vel_grid(n,grid_ref)%y(j)
                   beta_y = vel_grid(n,grid_ref)%beta_y(j)
                   do i = i_min, i_max
                      x = vel_grid(n,grid_ref)%x(i)
                      beta_x = vel_grid(n,grid_ref)%beta_x(i)
                      beta3 = beta_x*beta_y*beta_z
                      call compute_maxwellian( x, y, z, coeff, mass, u, v, w, temp, df_eq )
                      df_eq = df_eq * beta3

                      write( file_unit, advance="no", fmt="(5e24.16)" ) x, y, z, beta3, time

                      if( enable_vel_df_vis )then
                         write( file_unit, advance="no", fmt="(2e24.16)" ) &
                              phi(n,dfloc)%value(i,j,k), phi(n,dfloc)%value(i,j,k) - df_eq
                      end if

                      if( enable_rot_df_vis )then
                         energy = zero
                         do l = 1, r_levels
                            if( r_modes .gt. 0 )then
                               write( file_unit, advance="no", fmt="(e24.16)" )&
                                    phi(n,dfloc)%rot(l,i,j,k)
                               energy = energy + phi(n,dfloc)%rot(l,i,j,k)*phi(n,dfloc)%rot_level(l)
                            else
                               write( file_unit, advance="no", fmt="(e24.16)" ) zero
                            end if
                         end do
                         write( file_unit, advance="no", fmt="(e24.16)" ) energy
                      end if

                      if( enable_vib_df_vis )then
                         energy = zero
                         do l = 1, v_levels
                            if( v_modes .gt. 0 )then
                               write( file_unit, advance="no", fmt="(e24.16)" )&
                                    phi(n,dfloc)%vib(l,i,j,k)
                               energy = energy + phi(n,dfloc)%vib(l,i,j,k)*phi(n,dfloc)%vib_level(l)
                            else
                               write( file_unit, advance="no", fmt="(e24.16)" ) zero
                            end if
                         end do
                         write( file_unit, advance="no", fmt="(e24.16)" ) energy
                      end if

                      write( file_unit, advance="yes", fmt="(1x)" )

                   end do
                end do
             end do

          case default
             write(*,*) "Error: invalid df_output_dimension: ", df_output_dimension
             stop

          end select

          close( file_unit )

       end do
    end do

    return
  end subroutine write_df

  subroutine write_rot_levels( time, phi, molecule, props, vel_grid )

    use VelocityGrid
    use PhysicalGrid
    use PhysicalProperties

    implicit none

    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(PropertiesType), dimension(:), intent(in) :: props
    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    double precision, intent(in) :: time

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: modes, levels
    integer :: n, loc, i, j, k, l
    integer :: grid_ref, dfloc

    double precision :: energy_level, density_fraction
    double precision :: spin, degeneracy, density

    integer :: file_unit
    character(len=128) :: filename
    
    do n = 1, num_species

       modes = molecule(n)%rot_modes

       do loc = 1, num_output_locations

          dfloc = output_locations(loc)

          file_unit = rot_levels_vis_unit + (n-1)*num_output_locations + (loc-1)
          filename = trim(adjustl(rot_levels_vis_filename))//"_"//&
               trim(adjustl(int2char(n)))//"_"//&
               trim(adjustl(int2char(output_locations(loc))))//".dat"

          open( unit=file_unit, file=filename, status="old", position="append" )
          !>>>>>>>>>>
!!$          open( unit=77+loc, file="vel_dependent_rot_levels_"//trim(adjustl(int2char(loc)))//".dat", &
!!$                  status="old", position="append" )
          !<<<<<<<<<<

          call get_spatial_reference( output_locations(loc), grid_ref )

          i_min = vel_grid(n,grid_ref)%i_min
          i_max = vel_grid(n,grid_ref)%i_max
          j_min = vel_grid(n,grid_ref)%j_min
          j_max = vel_grid(n,grid_ref)%j_max
          k_min = vel_grid(n,grid_ref)%k_min
          k_max = vel_grid(n,grid_ref)%k_max

          levels = phi(n,dfloc)%num_rot_levels

          if( modes .eq. 0 ) cycle

          write( file_unit, advance="yes", fmt="(a7,i3)" )"zone i=",levels
          !>>>>>>>>>>
!!$          write( 77+loc, advance="yes", fmt="(a7,i3)" )"zone i=",levels
!!$          i = i_min + 6
!!$          j = j_min + 15
!!$          k = k_min + 15
!!$          density = sum(phi(n,dfloc)%rot(:,i,j,k))
!!$          do l = 1, levels
!!$             energy_level = phi(n,dfloc)%rot_level(l)
!!$             density_fraction = zero
!!$             if( mod(l-1,2) .eq. 0 )then
!!$                spin = molecule(n)%even_spin
!!$             else
!!$                spin = molecule(n)%odd_spin
!!$             end if
!!$             degeneracy = two * dble(l-1) + one
!!$
!!$             write( 77+loc, advance="yes", fmt="(e24.16,i5,5e24.16)" )&
!!$                  time, (l-1), energy_level, phi(n,dfloc)%rot(l,i,j,k)/density, &
!!$                  phi(n,dfloc)%rot(l,i,j,k)*energy_level, spin, degeneracy
!!$          end do
!!$
!!$          write( 77+loc, advance="yes", fmt="(a7,i3)" )"zone i=",levels
!!$          i = i_min + 6
!!$          j = j_min + 16
!!$          k = k_min + 15
!!$          density = sum(phi(n,dfloc)%rot(:,i,j,k))
!!$          do l = 1, levels
!!$             energy_level = phi(n,dfloc)%rot_level(l)
!!$             density_fraction = zero
!!$             if( mod(l-1,2) .eq. 0 )then
!!$                spin = molecule(n)%even_spin
!!$             else
!!$                spin = molecule(n)%odd_spin
!!$             end if
!!$             degeneracy = two * dble(l-1) + one
!!$
!!$             write( 77+loc, advance="yes", fmt="(e24.16,i5,5e24.16)" )&
!!$                  time, (l-1), energy_level, phi(n,dfloc)%rot(l,i,j,k)/density, &
!!$                  phi(n,dfloc)%rot(l,i,j,k)*energy_level, spin, degeneracy
!!$          end do
!!$
!!$          write( 77+loc, advance="yes", fmt="(a7,i3)" )"zone i=",levels
!!$          i = i_min + 6
!!$          j = j_min + 17
!!$          k = k_min + 15
!!$          density = sum(phi(n,dfloc)%rot(:,i,j,k))
!!$          do l = 1, levels
!!$             energy_level = phi(n,dfloc)%rot_level(l)
!!$             density_fraction = zero
!!$             if( mod(l-1,2) .eq. 0 )then
!!$                spin = molecule(n)%even_spin
!!$             else
!!$                spin = molecule(n)%odd_spin
!!$             end if
!!$             degeneracy = two * dble(l-1) + one
!!$
!!$             write( 77+loc, advance="yes", fmt="(e24.16,i5,5e24.16)" )&
!!$                  time, (l-1), energy_level, phi(n,dfloc)%rot(l,i,j,k)/density, &
!!$                  phi(n,dfloc)%rot(l,i,j,k)*energy_level, spin, degeneracy
!!$          end do
!!$
!!$          write( 77+loc, advance="yes", fmt="(a7,i3)" )"zone i=",levels
!!$          i = i_min + 10
!!$          j = j_min + 15
!!$          k = k_min + 15
!!$          density = sum(phi(n,dfloc)%rot(:,i,j,k))
!!$          do l = 1, levels
!!$             energy_level = phi(n,dfloc)%rot_level(l)
!!$             density_fraction = zero
!!$             if( mod(l-1,2) .eq. 0 )then
!!$                spin = molecule(n)%even_spin
!!$             else
!!$                spin = molecule(n)%odd_spin
!!$             end if
!!$             degeneracy = two * dble(l-1) + one
!!$
!!$             write( 77+loc, advance="yes", fmt="(e24.16,i5,5e24.16)" )&
!!$                  time, (l-1), energy_level, phi(n,dfloc)%rot(l,i,j,k)/density, &
!!$                  phi(n,dfloc)%rot(l,i,j,k)*energy_level, spin, degeneracy
!!$          end do
!!$
!!$          write( 77+loc, advance="yes", fmt="(a7,i3)" )"zone i=",levels
!!$          i = i_min + 15
!!$          j = j_min + 15
!!$          k = k_min + 15
!!$          density = sum(phi(n,dfloc)%rot(:,i,j,k))
!!$          do l = 1, levels
!!$             energy_level = phi(n,dfloc)%rot_level(l)
!!$             density_fraction = zero
!!$             if( mod(l-1,2) .eq. 0 )then
!!$                spin = molecule(n)%even_spin
!!$             else
!!$                spin = molecule(n)%odd_spin
!!$             end if
!!$             degeneracy = two * dble(l-1) + one
!!$
!!$             write( 77+loc, advance="yes", fmt="(e24.16,i5,5e24.16)" )&
!!$                  time, (l-1), energy_level, phi(n,dfloc)%rot(l,i,j,k)/density, &
!!$                  phi(n,dfloc)%rot(l,i,j,k)*energy_level, spin, degeneracy
!!$          end do
!!$
!!$          write( 77+loc, advance="yes", fmt="(a7,i3)" )"zone i=",levels
!!$          i = i_min + 15
!!$          j = j_min + 17
!!$          k = k_min + 15
!!$          density = sum(phi(n,dfloc)%rot(:,i,j,k))
!!$          do l = 1, levels
!!$             energy_level = phi(n,dfloc)%rot_level(l)
!!$             density_fraction = zero
!!$             if( mod(l-1,2) .eq. 0 )then
!!$                spin = molecule(n)%even_spin
!!$             else
!!$                spin = molecule(n)%odd_spin
!!$             end if
!!$             degeneracy = two * dble(l-1) + one
!!$
!!$             write( 77+loc, advance="yes", fmt="(e24.16,i5,5e24.16)" )&
!!$                  time, (l-1), energy_level, phi(n,dfloc)%rot(l,i,j,k)/density, &
!!$                  phi(n,dfloc)%rot(l,i,j,k)*energy_level, spin, degeneracy
!!$          end do
!!$
!!$          write( 77+loc, advance="yes", fmt="(a7,i3)" )"zone i=",levels
!!$          i = i_min + 15
!!$          j = j_min + 19
!!$          k = k_min + 15
!!$          density = sum(phi(n,dfloc)%rot(:,i,j,k))
!!$          do l = 1, levels
!!$             energy_level = phi(n,dfloc)%rot_level(l)
!!$             density_fraction = zero
!!$             if( mod(l-1,2) .eq. 0 )then
!!$                spin = molecule(n)%even_spin
!!$             else
!!$                spin = molecule(n)%odd_spin
!!$             end if
!!$             degeneracy = two * dble(l-1) + one
!!$
!!$             write( 77+loc, advance="yes", fmt="(e24.16,i5,5e24.16)" )&
!!$                  time, (l-1), energy_level, phi(n,dfloc)%rot(l,i,j,k)/density, &
!!$                  phi(n,dfloc)%rot(l,i,j,k)*energy_level, spin, degeneracy
!!$          end do
!!$
!!$          write( 77+loc, advance="yes", fmt="(a7,i3)" )"zone i=",levels
!!$          i = i_min + 18
!!$          j = j_min + 15
!!$          k = k_min + 15
!!$          density = sum(phi(n,dfloc)%rot(:,i,j,k))
!!$          do l = 1, levels
!!$             energy_level = phi(n,dfloc)%rot_level(l)
!!$             density_fraction = zero
!!$             if( mod(l-1,2) .eq. 0 )then
!!$                spin = molecule(n)%even_spin
!!$             else
!!$                spin = molecule(n)%odd_spin
!!$             end if
!!$             degeneracy = two * dble(l-1) + one
!!$
!!$             write( 77+loc, advance="yes", fmt="(e24.16,i5,5e24.16)" )&
!!$                  time, (l-1), energy_level, phi(n,dfloc)%rot(l,i,j,k)/density, &
!!$                  phi(n,dfloc)%rot(l,i,j,k)*energy_level, spin, degeneracy
!!$          end do
!!$
!!$          write( 77+loc, advance="yes", fmt="(a7,i3)" )"zone i=",levels
!!$          i = i_min + 18
!!$          j = j_min + 17
!!$          k = k_min + 15
!!$          density = sum(phi(n,dfloc)%rot(:,i,j,k))
!!$          do l = 1, levels
!!$             energy_level = phi(n,dfloc)%rot_level(l)
!!$             density_fraction = zero
!!$             if( mod(l-1,2) .eq. 0 )then
!!$                spin = molecule(n)%even_spin
!!$             else
!!$                spin = molecule(n)%odd_spin
!!$             end if
!!$             degeneracy = two * dble(l-1) + one
!!$
!!$             write( 77+loc, advance="yes", fmt="(e24.16,i5,5e24.16)" )&
!!$                  time, (l-1), energy_level, phi(n,dfloc)%rot(l,i,j,k)/density, &
!!$                  phi(n,dfloc)%rot(l,i,j,k)*energy_level, spin, degeneracy
!!$          end do
!!$
!!$          write( 77+loc, advance="yes", fmt="(a7,i3)" )"zone i=",levels
!!$          i = i_min + 4
!!$          j = j_min + 15
!!$          k = k_min + 15
!!$          density = sum(phi(n,dfloc)%rot(:,i,j,k))
!!$          do l = 1, levels
!!$             energy_level = phi(n,dfloc)%rot_level(l)
!!$             density_fraction = zero
!!$             if( mod(l-1,2) .eq. 0 )then
!!$                spin = molecule(n)%even_spin
!!$             else
!!$                spin = molecule(n)%odd_spin
!!$             end if
!!$             degeneracy = two * dble(l-1) + one
!!$
!!$             write( 77+loc, advance="yes", fmt="(e24.16,i5,5e24.16)" )&
!!$                  time, (l-1), energy_level, phi(n,dfloc)%rot(l,i,j,k)/density, &
!!$                  phi(n,dfloc)%rot(l,i,j,k)*energy_level, spin, degeneracy
!!$          end do
!!$
!!$          close( 77+loc )

          !<<<<<<<<<<

          density = props(loc)%dens(n)

          do l = 1, levels
             energy_level = phi(n,dfloc)%rot_level(l)
             density_fraction = zero
             if( mod(l-1,2) .eq. 0 )then
                spin = molecule(n)%even_spin
             else
                spin = molecule(n)%odd_spin
             end if
             degeneracy = two * dble(l-1) + one

             do k = k_min, k_max
                do j = j_min, j_max
                   do i = i_min, i_max
                      density_fraction = density_fraction + phi(n,dfloc)%rot(l,i,j,k)
                   end do
                end do
             end do
             
             write( file_unit, advance="yes", fmt="(e24.16,i5,5e24.16)" )&
                  time, (l-1), energy_level, density_fraction/density, density_fraction*energy_level, &
                  spin, degeneracy

          end do

          close( file_unit )

       end do
    end do

    return
  end subroutine write_rot_levels

  subroutine write_vib_levels( time, phi, molecule, props, vel_grid )

    use VelocityGrid
    use PhysicalGrid
    use PhysicalProperties

    implicit none

    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(PropertiesType), dimension(:), intent(in) :: props
    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    double precision, intent(in) :: time

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: modes, levels
    integer :: n, loc, i, j, k, l
    integer :: grid_ref, dfloc

    double precision :: energy_level, density_fraction, density

    integer :: file_unit
    character(len=128) :: filename


    do n = 1, num_species

       modes = molecule(n)%vib_modes

       do loc = 1, num_output_locations

          dfloc = output_locations(loc)

          file_unit = vib_levels_vis_unit + (n-1)*num_output_locations + (loc-1)
          filename = trim(adjustl(vib_levels_vis_filename))//"_"//&
               trim(adjustl(int2char(n)))//"_"//&
               trim(adjustl(int2char(output_locations(loc))))//".dat"

          open( unit=file_unit, file=filename, status="old", position="append" )

          call get_spatial_reference( output_locations(loc), grid_ref )

          i_min = vel_grid(n,grid_ref)%i_min
          i_max = vel_grid(n,grid_ref)%i_max
          j_min = vel_grid(n,grid_ref)%j_min
          j_max = vel_grid(n,grid_ref)%j_max
          k_min = vel_grid(n,grid_ref)%k_min
          k_max = vel_grid(n,grid_ref)%k_max

          levels = phi(n,dfloc)%num_vib_levels

          if( modes .eq. 0 ) cycle

          write( file_unit, advance="yes", fmt="(a7,i3)" )"zone i=",levels

          density = props(loc)%dens(n)

          do l = 1, levels
             energy_level = phi(n,dfloc)%vib_level(l)
             density_fraction = zero

             do k = k_min, k_max
                do j = j_min, j_max
                   do i = i_min, i_max
                      density_fraction = density_fraction + phi(n,dfloc)%vib(l,i,j,k)
                   end do
                end do
             end do

             write( file_unit, advance="yes", fmt="(e24.16,i3,3e24.16)" )&
                  time, (l-1), energy_level, density_fraction/density, density_fraction*energy_level
          end do

          close( file_unit )

       end do
    end do

    return
  end subroutine write_vib_levels

  subroutine write_properties( time, phi, molecule, vel_grid, properties )

    use VelocityGrid
    use PhysicalProperties
    use PhysicalGrid

    implicit none

    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    type(PropertiesType), dimension(:), intent(in) :: properties
    double precision, intent(in) :: time

    integer :: nx_space, ny_space, np, ns, prop, i

    double precision :: delta_x, delta_y

    integer :: file_unit
    character(len=128) :: filename

    call get_nspace( nx_space, ny_space )
    call get_delta_x( delta_x, delta_y )

    call get_data_to_output( phi, molecule, vel_grid, properties, time )

    file_unit = properties_vis_unit
    filename = trim(adjustl(properties_vis_filename))//".dat"

    open( unit=file_unit, file=filename, status="old", position="append" )

    if( nx_space .gt. 1 )then
       write( file_unit, advance="yes", fmt="(a7,i4)" ) "zone i=",nx_space
    end if

    do np = 1, nx_space
       write( file_unit, advance="no", fmt="(2e24.15)" ) time, dble(np-1)*delta_x

       do prop = 1, number_properties
          if( properties_to_output(prop) ) &
               write( file_unit, advance="no", fmt="(e24.16)" ) total_data(np,prop)
       end do

       if( write_moment .gt. 0 )then
          do i = 1, write_moment
             write( file_unit, advance="no", fmt="(e24.16)" ) total_moments(np,i)
          end do
       end if

       write( file_unit, advance="yes", fmt="(1x)" )
    end do
    
    close( file_unit )

    do ns = 1, num_species

       file_unit = properties_vis_unit + ns
       filename = trim(adjustl(properties_vis_filename))//"_"//trim(adjustl(int2char(ns)))//".dat"

       open( unit=file_unit, file=filename, status="old", position="append" )

       if( nx_space .gt. 1 )then
          write( file_unit, advance="yes", fmt="(a7,i4)" ) "zone i=",nx_space
       end if

       do np = 1, nx_space
          write( file_unit, advance="no", fmt="(2e24.15)" ) time, dble(np-1)*delta_x

          do prop = 1, number_properties
             if( properties_to_output(prop) ) write( file_unit, advance="no", fmt="(e24.16)" )&
                  species_data(ns,np,prop)
          end do

          if( write_moment .gt. 0 )then
             do i = 1, write_moment
                write( file_unit, advance="no", fmt="(e24.16)" ) species_moments(ns,np,i)
             end do
          end if

          write( file_unit, advance="yes", fmt="(1x)" )
       end do

       close( file_unit )

    end do

    return
  end subroutine write_properties

  subroutine write_bkw_properties( time, molecule, vel_grid, properties )

    use VelocityGrid
    use PhysicalProperties
    use PhysicalGrid

    implicit none

    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(PropertiesType), dimension(:), intent(in) :: properties
    double precision, intent(in) :: time

    integer :: nx_space, ny_space, prop, i

    double precision :: delta_x, delta_y

    integer :: file_unit
    character(len=128) :: filename

    call get_nspace( nx_space, ny_space )
    call get_delta_x( delta_x, delta_y )

    call get_data_to_output( bkw, molecule, vel_grid, properties, time )

    file_unit = bkw_analytic_vis_unit
    filename = trim(adjustl(properties_vis_filename))//"_bkw.dat"

    open( unit=file_unit, file=filename, status="old", position="append" )

    if( nx_space .gt. 1 )then
       write( file_unit, advance="yes", fmt="(a7,i4)" ) "zone i=",nx_space
    end if

    write( file_unit, advance="no", fmt="(2e24.15)" ) time, zero*delta_x

    do prop = 1, number_properties
       if( properties_to_output(prop) ) &
            write( file_unit, advance="no", fmt="(e24.16)" ) species_data(1,1,prop)
    end do

    if( write_moment .gt. 0 )then
       do i = 1, write_moment
          write( file_unit, advance="no", fmt="(e24.16)" ) species_moments(1,1,i)
       end do
    end if

    write( file_unit, advance="yes", fmt="(1x)" )

    close( file_unit )

    return
  end subroutine write_bkw_properties

  subroutine write_collisions( time, num_collisions )

    use VelocityGrid
    use PhysicalGrid

    implicit none

    integer, dimension(:,:,:), intent(in) :: num_collisions
    double precision, intent(in) :: time

    integer :: n, m
    integer :: np, nx_space, ny_space
    integer :: tot_collisions
    double precision :: delta_x, delta_y

    integer :: file_unit
    character(len=128) :: filename

    call get_nspace( nx_space, ny_space )
    call get_delta_x( delta_x, delta_y )

    file_unit = collisions_vis_unit
    filename = trim(adjustl(collisions_vis_filename))//".dat"

    open( unit=file_unit, file=filename, status="old", position="append" )

    if( nx_space .gt. 1 )then
       write( file_unit, advance="yes", fmt="(a7,i4)" ) "zone i=", nx_space
    end if

    do np = 1, nx_space
       write( file_unit, advance="no", fmt="(2e24.16)" ) time, dble(np)*delta_x

       tot_collisions = 0
       do n = 1, num_species
          do m = n, num_species
             write( file_unit, advance="no", fmt="(i12)" ) num_collisions(np,n,m)
             tot_collisions = tot_collisions + num_collisions(np,n,m)
          end do
       end do

       write( file_unit, advance="yes", fmt="(i12)" ) tot_collisions
    end do

    close( file_unit )

    return
  end subroutine write_collisions

  subroutine write_speed_df( phi, vel_grid, properties, molecule )

    use VelocityGrid
    use PhysicalGrid
    use PhysicalProperties

    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(PropertiesType), dimension(:), intent(in) :: properties

    double precision :: df, dev, df_eq

    double precision :: mass, dens, temp, u, v, w, coeff
    double precision :: x, y, z
    double precision :: beta_x, beta_y, beta_z, beta3

    integer :: i, j, k
    integer :: num_points, num_speeds
    integer :: n, loc, l
    integer :: grid_ref, dfloc

    integer :: file_unit
    character(len=128) :: filename

    do n = 1, num_species
       do loc = 1, num_output_locations

          dfloc = output_locations(loc)

          file_unit = speed_vis_unit + (n-1)*num_output_locations + (loc-1)
          filename = trim(adjustl(speed_vis_filename))//"_"//&
               trim(adjustl(int2char(n)))//"_"//&
               trim(adjustl(int2char(output_locations(loc))))//".dat"

          open( unit=file_unit, file=filename, status="old", position="append" )

          call get_spatial_reference( output_locations(loc), grid_ref )
          num_points = vel_grid(n,grid_ref)%num_points
          num_speeds = sorted_speed(n)%num_speeds

          write( file_unit, advance="yes", fmt="(a7,i4)" ) "zone i=",num_speeds

          ! Equilibrium calculation
          mass  = molecule(n)%mass
          dens  = properties(dfloc)%dens(n)
          u     = properties(dfloc)%x_vel(n)
          v     = properties(dfloc)%y_vel(n)
          w     = properties(dfloc)%z_vel(n)
          temp  = properties(dfloc)%temp(n)
          coeff = dens * sqrt( mass * mass * mass ) / sqrt( pi * pi * pi * temp * temp * temp )

          i = sorted_speed(n)%i(1)
          j = sorted_speed(n)%j(1)
          k = sorted_speed(n)%k(1)

          beta_x = vel_grid(n,grid_ref)%beta_x(i)
          beta_y = vel_grid(n,grid_ref)%beta_y(j)
          beta_z = vel_grid(n,grid_ref)%beta_z(k)
          beta3  = beta_x * beta_y * beta_z

          x = vel_grid(n,grid_ref)%x(i)
          y = vel_grid(n,grid_ref)%y(j)
          z = vel_grid(n,grid_ref)%z(k)

          call compute_maxwellian( x, y, z, coeff, mass, u, v, w, temp, df_eq )
          df_eq = df_eq * beta3

          df = phi(n,dfloc)%value(i,j,k)
          dev = df - df_eq

          do l = 2, num_points

             i = sorted_speed(n)%i(l)
             j = sorted_speed(n)%j(l)
             k = sorted_speed(n)%k(l)

             beta_x = vel_grid(n,grid_ref)%beta_x(i)
             beta_y = vel_grid(n,grid_ref)%beta_y(j)
             beta_z = vel_grid(n,grid_ref)%beta_z(k)
             beta3  = beta_x * beta_y * beta_z
             
             x = vel_grid(n,grid_ref)%x(i)
             y = vel_grid(n,grid_ref)%y(j)
             z = vel_grid(n,grid_ref)%z(k)

             call compute_maxwellian( x, y, z, coeff, mass, u, v, w, temp, df_eq )
             df_eq = df_eq * beta3

             if( sorted_speed(n)%speeds(l) .ne. sorted_speed(n)%speeds(l-1) )then
                write( file_unit, advance="yes", fmt="(3e24.16)" ) sorted_speed(n)%speeds(l-1), df, dev
                df = phi(n,dfloc)%value(i,j,k)
                dev = df - df_eq
             else
                df = df + phi(n,dfloc)%value(i,j,k)
                dev = dev + phi(n,dfloc)%value(i,j,k) - df_eq
             end if

          end do

          if( sorted_speed(n)%speeds(num_points) .eq. sorted_speed(n)%speeds(num_points-1) )&
               write( file_unit, advance="yes", fmt="(3e24.16)" ) sorted_speed(n)%speeds(num_points), df, dev

          close( file_unit )

       end do
    end do

    return
  end subroutine write_speed_df

  subroutine write_speed_rot_df( phi, vel_grid, properties, molecule )

    use VelocityGrid
    use PhysicalGrid
    use PhysicalProperties

    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(PropertiesType), dimension(:), intent(in) :: properties

    double precision :: mass, dens, temp, u, v, w, coeff
    double precision :: x, y, z
    double precision :: beta_x, beta_y, beta_z, beta3

    double precision :: df, df_eq, dev, total_dev
    double precision, allocatable, dimension(:) :: rot_df, rot_eq, rot_dev, temporary

    integer :: i, j, k
    integer :: num_points, num_speeds, rot_modes, rot_levels
    integer :: n, loc, l, m
    integer :: grid_ref, dfloc

    double precision :: delta_speed, base_speed

    integer :: status
    integer :: file_unit
    character(len=128) :: filename

    do n = 1, num_species

       rot_modes = molecule(n)%rot_modes

       if( rot_modes .eq. 0 ) cycle

       do loc = 1, num_output_locations

          dfloc = output_locations(loc)

          file_unit = speed_rot_vis_unit + (n-1)*num_output_locations + (loc-1)
          filename = trim(adjustl(speed_vis_filename))//"_"//&
               trim(adjustl(int2char(n)))//"_"//&
               trim(adjustl(int2char(dfloc)))//"_rot.dat"

          open( unit=file_unit, file=filename, status="old", position="append" )

          call get_spatial_reference( output_locations(loc), grid_ref )
          num_points = vel_grid(n,grid_ref)%num_points

          num_speeds = sorted_speed(n)%num_speeds
          rot_levels = phi(n,dfloc)%num_rot_levels

          allocate( rot_df(1:rot_levels), STAT=status )
          call allocate_error_check( status, "rot levels for speed output" )

          allocate( rot_eq(1:rot_levels), STAT=status )
          call allocate_error_check( status, "equilibrium rot levels for speed output" )

          allocate( rot_dev(1:rot_levels), STAT=status )
          call allocate_error_check( status, "deviation rot levels for speed output" )

          write( file_unit, advance="yes", fmt="(a7,i4,a4,i4)" )&
               "zone i=",rot_levels,", j=",num_speeds
          
          ! Equilibrium calculation
          mass  = molecule(n)%mass
          dens  = properties(dfloc)%dens(n)
          u     = properties(dfloc)%x_vel(n)
          v     = properties(dfloc)%y_vel(n)
          w     = properties(dfloc)%z_vel(n)
          temp  = properties(dfloc)%temp(n)
          coeff = dens * sqrt( mass * mass * mass ) / sqrt( pi * pi * pi * temp * temp * temp )

          allocate( temporary(1:rot_levels) )
          call compute_rot_distribution( rot_eq, temporary, molecule(n), temp, rot_levels, n )
          deallocate( temporary )

          i = sorted_speed(n)%i(1)
          j = sorted_speed(n)%j(1)
          k = sorted_speed(n)%k(1)

          beta_x = vel_grid(n,grid_ref)%beta_x(i)
          beta_y = vel_grid(n,grid_ref)%beta_y(j)
          beta_z = vel_grid(n,grid_ref)%beta_z(k)
          beta3  = beta_x * beta_y * beta_z

          x = vel_grid(n,grid_ref)%x(i)
          y = vel_grid(n,grid_ref)%y(j)
          z = vel_grid(n,grid_ref)%z(k)

          call compute_maxwellian( x, y, z, coeff, mass, u, v, w, temp, df_eq )
          df_eq = df_eq * beta3
          
          df = phi(n,dfloc)%value(i,j,k)
          rot_df = phi(n,dfloc)%rot(:,i,j,k)

          dev = df - df_eq
          rot_dev = rot_df - rot_eq * df_eq
          total_dev = sum(abs(rot_dev))
          
          delta_speed = double_tol * vel_grid(n,grid_ref)%beta3_min**(0.5d0)
          base_speed = sorted_speed(n)%speeds(1)

          do l = 2, num_points

             i = sorted_speed(n)%i(l)
             j = sorted_speed(n)%j(l)
             k = sorted_speed(n)%k(l)

             beta_x = vel_grid(n,grid_ref)%beta_x(i)
             beta_y = vel_grid(n,grid_ref)%beta_y(j)
             beta_z = vel_grid(n,grid_ref)%beta_z(k)
             beta3  = beta_x * beta_y * beta_z
             
             x = vel_grid(n,grid_ref)%x(i)
             y = vel_grid(n,grid_ref)%y(j)
             z = vel_grid(n,grid_ref)%z(k)

             call compute_maxwellian( x, y, z, coeff, mass, u, v, w, temp, df_eq )
             df_eq = df_eq * beta3

             if( sorted_speed(n)%speeds(l) .ge. base_speed + delta_speed )then
                do m = 1, rot_levels
                   write( file_unit, advance="yes", fmt="(7e24.16)" )&
                        sorted_speed(n)%speeds(l-1), dble(m), &
                        rot_df(m), rot_df(m)/df, df, &
                        rot_dev(m), total_dev
                end do

                df = phi(n,dfloc)%value(i,j,k)
                rot_df = phi(n,dfloc)%rot(:,i,j,k)

                dev = df - df_eq
                rot_dev = rot_df - rot_eq * df_eq
                total_dev = sum(abs(rot_dev))
                base_speed = sorted_speed(n)%speeds(l)
             else
                df = df + phi(n,dfloc)%value(i,j,k)
                rot_df = rot_df + phi(n,dfloc)%rot(:,i,j,k)
                dev = dev + phi(n,dfloc)%value(i,j,k) - df_eq
                rot_dev = rot_dev + phi(n,dfloc)%rot(:,i,j,k) - rot_eq * df_eq
                total_dev = total_dev + sum(abs(phi(n,dfloc)%rot(:,i,j,k) - rot_eq * df_eq))
             end if

          end do

          if( sorted_speed(n)%speeds(num_points) .lt. base_speed + delta_speed )then
             do m = 1, rot_levels
                write( file_unit, advance="yes", fmt="(7e24.16)" )&
                     sorted_speed(n)%speeds(l-1), dble(m), &
                     rot_df(m), rot_df(m)/df, df, &
                     rot_dev(m), total_dev
             end do
          end if

          close( file_unit )

          deallocate( rot_df, STAT=status )
          call deallocate_error_check( status, "rot levels for speed output" )

          deallocate( rot_eq, STAT=status )
          call deallocate_error_check( status, "equilibrium rot levels for speed output" )

          deallocate( rot_dev, STAT=status )
          call deallocate_error_check( status, "deviation rot levels for speed output" )

       end do
    end do

    return
  end subroutine write_speed_rot_df

  subroutine write_speed_vib_df( phi, vel_grid, properties, molecule )

    use VelocityGrid
    use PhysicalGrid
    use PhysicalProperties

    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(PropertiesType), dimension(:), intent(in) :: properties

    double precision :: mass, dens, temp, u, v, w, coeff
    double precision :: x, y, z
    double precision :: beta_x, beta_y, beta_z, beta3

    double precision :: df, df_eq, dev, total_dev
    double precision, allocatable, dimension(:) :: vib_df, vib_eq, vib_dev, temporary

    double precision :: base_speed, delta_speed

    integer :: i, j, k
    integer :: num_points, num_speeds, vib_modes, vib_levels
    integer :: n, loc, l, m
    integer :: grid_ref, dfloc

    integer :: status
    integer :: file_unit
    character(len=128) :: filename

    do n = 1, num_species

       vib_modes = molecule(n)%vib_modes

       if( vib_modes .eq. 0 ) cycle

       do loc = 1, num_output_locations

          dfloc = output_locations(loc)

          file_unit = speed_vib_vis_unit + (n-1)*num_output_locations + (loc-1)
          filename = trim(adjustl(speed_vis_filename))//"_"//&
               trim(adjustl(int2char(n)))//"_"//&
               trim(adjustl(int2char(dfloc)))//"_vib.dat"

          open( unit=file_unit, file=filename, status="old", position="append" )

          call get_spatial_reference( output_locations(loc), grid_ref )
          num_points = vel_grid(n,grid_ref)%num_points

          num_speeds = sorted_speed(n)%num_speeds
          vib_levels = phi(n,dfloc)%num_vib_levels

          allocate( vib_df(1:vib_levels), STAT=status )
          call allocate_error_check( status, "vib levels for speed output" )

          allocate( vib_eq(1:vib_levels), STAT=status )
          call allocate_error_check( status, "equilibrium vib levels for speed output" )

          allocate( vib_dev(1:vib_levels), STAT=status )
          call allocate_error_check( status, "deviation vib levels for speed output" )

          write( file_unit, advance="yes", fmt="(a7,i4,a4,i4)" )&
               "zone i=",vib_levels,", j=",num_speeds

          ! Equilibrium calculation
          mass  = molecule(n)%mass
          dens  = properties(dfloc)%dens(n)
          u     = properties(dfloc)%x_vel(n)
          v     = properties(dfloc)%y_vel(n)
          w     = properties(dfloc)%z_vel(n)
          temp  = properties(dfloc)%temp(n)
          coeff = dens * sqrt( mass * mass * mass ) / sqrt( pi * pi * pi * temp * temp * temp )
          
          allocate( temporary(1:vib_levels) )
          call compute_vib_distribution( vib_eq, temporary, molecule(n), temp, vib_levels, n )
          deallocate( temporary )

          i = sorted_speed(n)%i(1)
          j = sorted_speed(n)%j(1)
          k = sorted_speed(n)%k(1)

          beta_x = vel_grid(n,grid_ref)%beta_x(i)
          beta_y = vel_grid(n,grid_ref)%beta_y(j)
          beta_z = vel_grid(n,grid_ref)%beta_z(k)
          beta3  = beta_x * beta_y * beta_z

          x = vel_grid(n,grid_ref)%x(i)
          y = vel_grid(n,grid_ref)%y(j)
          z = vel_grid(n,grid_ref)%z(k)

          call compute_maxwellian( x, y, z, coeff, mass, u, v, w, temp, df_eq )
          df_eq = df_eq * beta3

          df = phi(n,dfloc)%value(i,j,k)
          vib_df = phi(n,dfloc)%vib(:,i,j,k)

          dev = df - df_eq
          vib_dev = vib_df - vib_eq * df_eq
          total_dev = sum(abs(vib_dev))

          delta_speed = double_tol * vel_grid(n,grid_ref)%beta3_min**(0.5d0)
          base_speed = sorted_speed(n)%speeds(1)

          do l = 2, num_points
             
             i = sorted_speed(n)%i(l)
             j = sorted_speed(n)%j(l)
             k = sorted_speed(n)%k(l)

             beta_x = vel_grid(n,grid_ref)%beta_x(i)
             beta_y = vel_grid(n,grid_ref)%beta_y(j)
             beta_z = vel_grid(n,grid_ref)%beta_z(k)
             beta3  = beta_x * beta_y * beta_z

             x = vel_grid(n,grid_ref)%x(i)
             y = vel_grid(n,grid_ref)%y(j)
             z = vel_grid(n,grid_ref)%z(k)

             call compute_maxwellian( x, y, z, coeff, mass, u, v, w, temp, df_eq )
             df_eq = df_eq * beta3

             if( sorted_speed(n)%speeds(l) .ge. base_speed + delta_speed )then
                do m = 1, vib_levels
                   write( file_unit, advance="yes", fmt="(7e24.16)" )&
                        sorted_speed(n)%speeds(l-1), dble(m), &
                        vib_df(m), vib_df(m)/df, df, &
                        vib_dev(m), total_dev
                end do

                df = phi(n,dfloc)%value(i,j,k)
                vib_df = phi(n,dfloc)%vib(:,i,j,k)

                dev = df - df_eq
                vib_dev = vib_df - vib_eq * df_eq
                total_dev = sum(abs(vib_dev))
                base_speed = sorted_speed(n)%speeds(l)
             else
                df = df + phi(n,dfloc)%value(i,j,k)
                vib_df = vib_df + phi(n,dfloc)%vib(:,i,j,k)
                dev = dev + phi(n,dfloc)%value(i,j,k) - df_eq
                vib_dev = vib_dev + phi(n,dfloc)%vib(:,i,j,k) - vib_eq * df_eq
                total_dev = total_dev + sum(abs(phi(n,dfloc)%vib(:,i,j,k) - vib_eq * df_eq))
             end if

!!$             if( sorted_speed(n)%speeds(l) .ne. sorted_speed(n)%speeds(l-1) )then
!!$                do m = 1, vib_levels
!!$                   write( file_unit, advance="yes", fmt="(5e24.16)" )&
!!$                        sorted_speed(n)%speeds(l-1), dble(m), &
!!$                        vib_df(m), vib_df(m)/df, df
!!$                end do
!!$                df = phi(n,dfloc)%value(i,j,k)
!!$                vib_df = phi(n,dfloc)%vib(:,i,j,k)
!!$             else
!!$                df = df + phi(n,dfloc)%value(i,j,k)
!!$                vib_df = vib_df + phi(n,dfloc)%vib(:,i,j,k)
!!$             end if

          end do
          
          if( sorted_speed(n)%speeds(num_points) .lt. base_speed + delta_speed )then
             do m = 1, vib_levels
                write( file_unit, advance="yes", fmt="(7e24.16)" )&
                     sorted_speed(n)%speeds(l-1), dble(m), &
                     vib_df(m), vib_df(m)/df, df, &
                     vib_dev(m), total_dev
             end do
          end if
!!$          if( sorted_speed(n)%speeds(num_points) .eq. sorted_speed(n)%speeds(num_points-1) )then
!!$             do m = 1, vib_levels
!!$                write( file_unit, advance="yes", fmt="(4e24.16)" )&
!!$                     dble(m), sorted_speed(n)%speeds(num_points), &
!!$                     vib_df(m), vib_df(m)/df, df
!!$             end do
!!$          end if

          close( file_unit )

          deallocate( vib_df, STAT=status )
          call deallocate_error_check( status, "vib levels for speed output" )

          deallocate( vib_eq, STAT=status )
          call deallocate_error_check( status, "equilibrium vib levels for speed output" )

          deallocate( vib_dev, STAT=status )
          call deallocate_error_check( status, "deviation vib levels for speed output" )

       end do
    end do

    return
  end subroutine write_speed_vib_df

  !===================================================================================================
  !  Calculate Data to Output
  !===================================================================================================
  subroutine get_data_to_output( phi, molecule, vel_grid, props, time )

    use VelocityGrid
    use PhysicalProperties
    use PhysicalGrid
    use ShockConditions
    use BoundaryConditions

    implicit none

    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    type(PropertiesType), dimension(:), intent(in) :: props
    double precision, intent(in) :: time

    integer :: np, ns, i, grid_ref
    integer :: nx_space, ny_space

    double precision :: prop_value, gamma, mol_mass
    double precision :: T_eq, T_rot, T_vib, T_kin

    ! np - spatial element
    ! ns - species number

    call get_nspace( nx_space, ny_space )

    do np = 1, nx_space

       call get_spatial_reference( np, grid_ref )

       ! Density
       if( properties_to_output(1) )then
          do ns = 1, num_species
             species_data(ns,np,1) = props(np)%dens(ns)
          end do
          total_data(np,1) = props(np)%mix_dens
       end if

       ! X-Velocity
       if( properties_to_output(2) )then
          do ns = 1, num_species
             species_data(ns,np,2) = props(np)%x_vel(ns)
          end do
          total_data(np,2) = props(np)%mix_x_vel
       end if

       ! Y-Velocity
       if( properties_to_output(3) )then
          do ns = 1, num_species
             species_data(ns,np,3) = props(np)%y_vel(ns)
          end do
          total_data(np,3) = props(np)%mix_y_vel
       end if

       ! Z-Velocity
       if( properties_to_output(4) )then
          do ns = 1, num_species
             species_data(ns,np,4) = props(np)%z_vel(ns)
          end do
          total_data(np,4) = props(np)%mix_z_vel
       end if
       
       ! Speed
       if( properties_to_output(5) )then
          do ns = 1, num_species
             call compute_speed( prop_value, props(np)%x_vel(ns), props(np)%y_vel(ns), props(np)%z_vel(ns) )
             species_data(ns,np,5) = prop_value
          end do
          call compute_speed( prop_value, props(np)%mix_x_vel, props(np)%mix_y_vel, props(np)%mix_z_vel )
          total_data(np,5) = prop_value
       end if

       ! Mach Number
       if( properties_to_output(6) )then
          call calculate_specific_heat_ratio( gamma, molecule )
          call calculate_molecular_mass( mol_mass, molecule )
          do ns = 1, num_species
             prop_value = props(np)%x_vel(ns) / sqrt( one_half * gamma * props(np)%tr_temp(ns) / mol_mass )
             species_data(ns,np,6) = prop_value
          end do
          prop_value = props(np)%mix_x_vel/sqrt( one_half * gamma * props(np)%mix_tr_temp / mol_mass )
          total_data(np,6) = prop_value
       end if

       ! Translational Temperature
       if( properties_to_output(7) )then
          do ns = 1, num_species
             species_data(ns,np,7) = props(np)%tr_temp(ns)
          end do
          total_data(np,7) = props(np)%mix_tr_temp
       end if

       ! Rotational Temperature
       if( properties_to_output(8) )then
          do ns = 1, num_species
             species_data(ns,np,8) = props(np)%rot_temp(ns)
          end do
          total_data(np,8) = props(np)%mix_rot_temp
       end if

       ! Vibrational Temperature
       if( properties_to_output(9) )then
          do ns = 1, num_species
             species_data(ns,np,9) = props(np)%vib_temp(ns)
          end do
          total_data(np,9) = props(np)%mix_vib_temp
       end if

       ! Total Temperature
       if( properties_to_output(10) )then
          do ns = 1, num_species
             species_data(ns,np,10) = props(np)%temp(ns)
          end do
          total_data(np,10) = props(np)%mix_temp
       end if

       ! Temperature X
       if( properties_to_output(11) )then
          do ns = 1, num_species
             call compute_directional_temperature( prop_value, props(np)%dens(ns), &
                  props(np)%x_vel(ns), x_dir, phi(ns,np), molecule(ns), vel_grid(ns,grid_ref) )
             species_data(ns,np,11) = prop_value
          end do
          total_data(np,11) = zero
       end if

       ! Temperature Y
       if( properties_to_output(12) )then
          do ns = 1, num_species
            call compute_directional_temperature( prop_value, props(np)%dens(ns), &
                  props(np)%y_vel(ns), y_dir, phi(ns,np), molecule(ns), vel_grid(ns,grid_ref) )
             species_data(ns,np,12) = prop_value
          end do
          total_data(np,12) = zero
       end if

       ! Temperature X
       if( properties_to_output(13) )then
          do ns = 1, num_species
            call compute_directional_temperature( prop_value, props(np)%dens(ns), &
                  props(np)%z_vel(ns), z_dir, phi(ns,np), molecule(ns), vel_grid(ns,grid_ref) )
             species_data(ns,np,13) = prop_value
          end do
          total_data(np,13) = zero
       end if

       ! Translational Energy
       if( properties_to_output(14) )then
          do ns = 1, num_species
             species_data(ns,np,14) = props(np)%tr_energy(ns)
          end do
          total_data(np,14) = props(np)%mix_tr_energy
       end if
       
       ! Rotational Energy
       if( properties_to_output(15) )then
          do ns = 1, num_species
             species_data(ns,np,15) = props(np)%rot_energy(ns)
          end do
          total_data(np,15) = props(np)%mix_rot_energy
       end if

       ! Vibrational Energy
       if( properties_to_output(16) )then
          do ns = 1, num_species
             species_data(ns,np,16) = props(np)%vib_energy(ns)
          end do
          total_data(np,16) = props(np)%mix_vib_energy
       end if

       ! Total Energy
       if( properties_to_output(17) )then
          do ns = 1, num_species
             species_data(ns,np,17) = props(np)%total_energy(ns)
          end do
          total_data(np,17) = props(np)%mix_energy
       end if

       ! Entropy
       if( properties_to_output(18) )then
          do ns = 1, num_species
             call compute_entropy( prop_value, props(np)%dens(ns), molecule(ns)%mass, &
                  phi(ns,np), vel_grid(ns,grid_ref) )
             species_data(ns,np,18) = prop_value
          end do
          total_data(np,18) = prop_value
       end if

       ! Heat Flux X
       if( properties_to_output(19) )then
          do ns = 1, num_species
             call compute_heat_flux( prop_value, props(np)%x_vel(ns), props(np)%y_vel(ns), props(np)%z_vel(ns), &
                  props(np)%dens(ns), x_dir, phi(ns,np), molecule(ns), vel_grid(ns,grid_ref) )
             species_data(ns,np,19) = prop_value
          end do
          total_data(np,19) = zero
       end if

       ! Heat Flux Y
       if( properties_to_output(20) )then
          do ns = 1, num_species
             call compute_heat_flux( prop_value, props(np)%x_vel(ns), props(np)%y_vel(ns), props(np)%z_vel(ns), &
                  props(np)%dens(ns), y_dir, phi(ns,np), molecule(ns), vel_grid(ns,grid_ref) )
             species_data(ns,np,20) = prop_value
          end do
          total_data(np,20) = zero
       end if
       
       ! Heat Flux Z
       if( properties_to_output(21) )then
          do ns = 1, num_species
             call compute_heat_flux( prop_value, props(np)%x_vel(ns), props(np)%y_vel(ns), props(np)%z_vel(ns), &
                  props(np)%dens(ns), z_dir, phi(ns,np), molecule(ns), vel_grid(ns,grid_ref) )
             species_data(ns,np,21) = prop_value
          end do
          total_data(np,21) = zero
       end if

       ! Shear Stress
       if( properties_to_output(22) )then
          do ns = 1, num_species
             species_data(ns,np,22) = zero
          end do
          total_data(np,22) = zero
       end if

       ! Moments
       if( write_moment .gt. 0 )then
          do i = 1, write_moment
             do ns = 1, num_species
                call compute_moment( prop_value, moment_numbers(i), props(np)%dens(ns), &
                     props(np)%x_vel(ns), props(np)%y_vel(ns), props(np)%z_vel(ns), props(np)%tr_temp(ns), &
                     phi(ns,np), molecule(ns), vel_grid(ns,grid_ref) )
                species_moments(ns,np,i) = prop_value
             end do
             total_moments(np,i) = zero
          end do
       end if

       ! Rotational Degrees of Freedom
       if( properties_to_output(23) )then
          do ns = 1, num_species
             species_data(ns,np,23) = props(np)%rot_dof(ns)
          end do
          total_data(np,23) = zero
       end if

       ! Vibrational Degrees of Freedom
       if( properties_to_output(24) )then
          do ns = 1, num_species
             species_data(ns,np,24) = props(np)%vib_dof(ns)
          end do
          total_data(np,24) = zero
       end if

       ! Analytic relaxations
       if( properties_to_output(25) )then
          if( abs( time ) .lt. double_tol )then
             do ns = 1, num_species
                species_data(ns,np,25) = props(np)%rot_temp(ns)
             end do
          else
             do ns = 1, num_species
!!$                species_data(ns,np,26) = species_data(ns,np,26) + &
!!$                     0.05d0 * ( species_data(ns,np,25) - species_data(ns,np,26) ) / &
!!$                     ( molecule(ns)%cZr * 3.0d0 / ( 3.0d0 + props(np)%rot_dof(ns) ) )
                T_eq = 4.0d0 * ( props(np)%tr_energy(ns) + props(np)%rot_energy(ns) )/&
                     ( props(np)%dens(ns) * ( 3.0d0 + props(np)%rot_dof(ns) ) )!props(np)%temp(ns)
                call get_right_wall_temp_rot( T_rot, ns )
                species_data(ns,np,25) = T_eq - ( T_eq - T_rot ) * &
                     exp( -2.0d0 * time/( molecule(ns)%cZr ) )
             end do
          end if
          total_data(np,25) = zero
       end if

       if( properties_to_output(26) )then
          if( abs( time ) .lt. double_tol )then
             do ns = 1, num_species
                species_data(ns,np,26) = props(np)%vib_temp(ns)
             end do
          else
             do ns = 1, num_species
!!$                species_data(ns,np,27) = species_data(ns,np,27) + &
!!$                     0.05d0 * ( species_data(ns,np,25) - species_data(ns,np,27) ) / &
!!$                     ( molecule(ns)%cZv * 3.0d0 / ( 3.0d0 + props(np)%vib_dof(ns) ) )
                T_eq = 4.0d0 * ( props(np)%tr_energy(ns) + props(np)%vib_energy(ns) )/&
                     ( props(np)%dens(ns) * ( 3.0d0 + props(np)%vib_dof(ns) ) )!props(np)%temp(ns) 
                call get_right_wall_temp_vib( T_vib, ns )
                species_data(ns,np,26) = T_eq - ( T_eq - T_vib ) * &
                     exp( -2.0d0 * time/( molecule(ns)%cZv ) )
             end do
          end if
          total_data(np,26) = zero
       end if


       if( properties_to_output(27) )then
          if( abs( time ) .lt. double_tol )then
             do ns = 1, num_species
                species_data(ns,np,27) = props(np)%tr_temp(ns)
             end do
          else
             do ns = 1, num_species
                species_data(ns,np,27) = ( 4.0d0 / 3.0d0 ) * ( props(np)%total_energy(ns) - &
                     ( props(np)%rot_dof(ns)/4.0d0 ) * species_data(ns,np,25) - &
                     ( props(np)%vib_dof(ns)/4.0d0 ) * species_data(ns,np,26)    )
!!$                species_data(ns,np,25) = ( 1.0d0 / 3.0d0 ) * ( 3.0d0 * species_data(ns,np,25) - &
!!$                     0.05d0 * props(np)%rot_dof(ns) * ( species_data(ns,np,25) - species_data(ns,np,26) ) / &
!!$                     ( molecule(ns)%cZr * 3.0d0 / ( 3.0d0 + props(np)%rot_dof(ns) ) ) - &
!!$                     0.05d0 * props(np)%vib_dof(ns) * ( species_data(ns,np,25) - species_data(ns,np,27) ) / &
!!$                     ( molecule(ns)%cZv * 3.0d0 / ( 3.0d0 + props(np)%vib_dof(ns) ) ) )
             end do
          end if
          total_data(np,27) = zero
       end if

       ! Pressure
       if( properties_to_output(28) )then
          do ns = 1, num_species
             prop_value = props(np)%dens(ns) * props(np)%tr_temp(ns)
             species_data(ns,np,28) = prop_value
          end do
          prop_value = props(np)%mix_dens * props(np)%mix_tr_temp
          total_data(np,28) = prop_value
       end if

    end do

    return
  end subroutine get_data_to_output

  subroutine sort_speed( vel_grid )

    use SortAndSearch
    use VelocityGrid

    implicit none

    type(VelocityGridType), dimension(:), intent(in) :: vel_grid

    integer :: count
    integer :: n, num_points
    integer :: i, j, k
    integer :: i_min, i_max, j_min, j_max, k_min, k_max

    double precision :: x, y, z

    double precision :: base_speed, delta_speed

    do n = 1, num_species
       
       num_points = vel_grid(n)%num_points
       
       i_min = vel_grid(n)%i_min
       i_max = vel_grid(n)%i_max
       j_min = vel_grid(n)%j_min
       j_max = vel_grid(n)%j_max
       k_min = vel_grid(n)%k_min
       k_max = vel_grid(n)%k_max

       allocate( sorted_speed(n)%speeds( 1:num_points ) )
       allocate( sorted_speed(n)%i( 1:num_points ) )
       allocate( sorted_speed(n)%j( 1:num_points ) )
       allocate( sorted_speed(n)%k( 1:num_points ) )
       
       count = 0

       do i = i_min, i_max
          x = vel_grid(n)%x(i)

          do j = j_min, j_max
             y = vel_grid(n)%y(j)

             do k = k_min, k_max
                z = vel_grid(n)%z(k)

                count = count+1
                
                sorted_speed(n)%i(count) = i
                sorted_speed(n)%j(count) = j
                sorted_speed(n)%k(count) = k

                sorted_speed(n)%speeds(count) = sqrt( x*x + y*y + z*z )

             end do
          end do
       end do

       call ssort( sorted_speed(n)%speeds, sorted_speed(n)%i, sorted_speed(n)%j, sorted_speed(n)%k, num_points )

       delta_speed = double_tol * vel_grid(n)%beta3_min**(0.5d0)

       count = 0
       base_speed = sorted_speed(n)%speeds(1)
       do i = 2, num_points
          if( sorted_speed(n)%speeds(i) .lt. base_speed + delta_speed )then
             cycle
          else
             count = count + 1
             base_speed = sorted_speed(n)%speeds(i)
          end if
!!$          if( sorted_speed(n)%speeds(i) .eq. sorted_speed(n)%speeds(i-1) )then
!!$             cycle
!!$          else
!!$             count = count + 1
!!$          end if
       end do

       if( sorted_speed(n)%speeds(num_points) .lt. base_speed + delta_speed )then
          count = count + 1
       end if
!!$       if( sorted_speed(n)%speeds(num_points) .eq. sorted_speed(n)%speeds(num_points-1) )then
!!$          count = count + 1
!!$       end if

       sorted_speed(n)%num_speeds = count

    end do

    return
  end subroutine sort_speed

end module Visualization
