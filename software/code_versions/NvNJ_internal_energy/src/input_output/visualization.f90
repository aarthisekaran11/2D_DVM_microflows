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

  ! File units
  integer, parameter :: base_unit = 100
  integer :: df_vis_unit
  integer :: rot_levels_vis_unit
  integer :: vib_levels_vis_unit
  integer :: properties_vis_unit
  integer :: collisions_vis_unit
  integer :: bkw_analytic_vis_unit
  integer :: rotvib_levels_vis_unit

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

  character(len=128) :: df_vis_filename
  character(len=128) :: rot_levels_vis_filename
  character(len=128) :: vib_levels_vis_filename
  character(len=128) :: properties_vis_filename
  character(len=128) :: collisions_vis_filename
  character(len=128) :: rotvib_levels_vis_filename

  integer :: df_vis_dump_freq
  integer :: rot_levels_vis_dump_freq
  integer :: vib_levels_vis_dump_freq
  integer :: properties_vis_dump_freq
  integer :: collisions_vis_dump_freq
  integer :: rotvib_levels_vis_dump_freq

  integer :: num_output_locations
  integer :: df_output_dimension
  integer, allocatable, dimension(:) :: output_locations

  ! Output properites
  integer, parameter :: number_properties = 24
  logical, dimension(1:number_properties) :: properties_to_output
  character(len=128), dimension(1:number_properties) :: property_names
      
  integer ::  write_moment
  integer, allocatable, dimension(:) :: moment_numbers

  double precision, allocatable, dimension(:,:) :: total_data, total_moments
  double precision, allocatable, dimension(:,:,:) :: species_data, species_moments

  ! Analytic BKW
  type(DistFuncType), dimension(1,1) :: bkw

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
       enable_collisions_vis_in )

    implicit none

    integer, intent(in) :: enable_vel_df_vis_in, enable_rot_df_vis_in, enable_vib_df_vis_in
    integer, intent(in) :: enable_rot_levels_vis_in, enable_vib_levels_vis_in
    integer, intent(in) :: enable_properties_vis_in, enable_collisions_vis_in

    enable_vel_df_vis = int2logical( enable_vel_df_vis_in )
    enable_rot_df_vis = int2logical( enable_rot_df_vis_in )
    enable_vib_df_vis = int2logical( enable_vib_df_vis_in )

    if( enable_vel_df_vis .or. enable_rot_df_vis .or. enable_vib_df_vis ) enable_df_vis = .true.

    enable_rot_levels_vis = int2logical( enable_rot_levels_vis_in )
    enable_vib_levels_vis = int2logical( enable_vib_levels_vis_in )
    enable_properties_vis = int2logical( enable_properties_vis_in )
    enable_collisions_vis = int2logical( enable_collisions_vis_in )

    return
  end subroutine enable_vis

  subroutine enable_bkw_analytic( enable_bkw_analytic_vis_in )
    
    implicit none
    
    integer, intent(in) :: enable_bkw_analytic_vis_in

    enable_bkw_analytic_vis = int2logical( enable_bkw_analytic_vis_in )

    return
  end subroutine enable_bkw_analytic

  subroutine set_vis_filename( df_vis_filename_in, rot_levels_vis_filename_in, vib_levels_vis_filename_in, &
       properties_vis_filename_in, collisions_vis_filename_in, rotvib_levels_vis_filename_in )

    implicit none

    character(len=128), intent(in) :: df_vis_filename_in
    character(len=128), intent(in) :: rot_levels_vis_filename_in, vib_levels_vis_filename_in
    character(len=128), intent(in) :: properties_vis_filename_in, collisions_vis_filename_in
    character(len=128), intent(in) :: rotvib_levels_vis_filename_in

    df_vis_filename         = df_vis_filename_in
    rot_levels_vis_filename = rot_levels_vis_filename_in
    vib_levels_vis_filename = vib_levels_vis_filename_in
    properties_vis_filename = properties_vis_filename_in
    collisions_vis_filename = collisions_vis_filename_in
    rotvib_levels_vis_filename = rotvib_levels_vis_filename_in

    return
  end subroutine set_vis_filename
  
  subroutine set_vis_dump_freq( df_vis_dump_freq_in, rot_levels_vis_dump_freq_in, &
       vib_levels_vis_dump_freq_in, properties_vis_dump_freq_in, collisions_vis_dump_freq_in, &
       rotvib_levels_vis_dump_freq_in )

    implicit none

    integer, intent(in) :: df_vis_dump_freq_in
    integer, intent(in) :: rot_levels_vis_dump_freq_in, vib_levels_vis_dump_freq_in
    integer, intent(in) :: properties_vis_dump_freq_in, collisions_vis_dump_freq_in
    integer, intent(in) :: rotvib_levels_vis_dump_freq_in

    df_vis_dump_freq         = df_vis_dump_freq_in
    rot_levels_vis_dump_freq = rot_levels_vis_dump_freq_in
    vib_levels_vis_dump_freq = vib_levels_vis_dump_freq_in
    properties_vis_dump_freq = properties_vis_dump_freq_in
    collisions_vis_dump_freq = collisions_vis_dump_freq_in
    rotvib_levels_vis_dump_freq = rotvib_levels_vis_dump_freq_in

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
       write_temp_x_in, write_temp_y_in, write_temp_z_in, &
       write_rot_dof_in, write_vib_dof_in, &
       write_tr_energy_in, write_rot_energy_in, write_vib_energy_in, write_tot_energy_in, &
       write_entropy_in, write_heat_flux_x_in, write_heat_flux_y_in, write_heat_flux_z_in, &
       write_shear_stress_in )

    implicit none

    integer, intent(in) :: write_dens_in, write_x_vel_in, write_y_vel_in, write_z_vel_in
    integer, intent(in) :: write_speed_in, write_mach_in
    integer, intent(in) :: write_tr_temp_in, write_rot_temp_in, write_vib_temp_in
    integer, intent(in) :: write_tot_temp_in
    integer, intent(in) :: write_temp_x_in, write_temp_y_in, write_temp_z_in
    integer, intent(in) :: write_rot_dof_in, write_vib_dof_in
    integer, intent(in) :: write_tot_energy_in, write_rot_energy_in, write_vib_energy_in, write_tr_energy_in
    integer, intent(in) :: write_entropy_in, write_heat_flux_x_in, write_heat_flux_y_in, write_heat_flux_z_in
    integer, intent(in) :: write_shear_stress_in


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
  subroutine initialize_visualization( vel_grid, molecule )

    use VelocityGrid
    use PhysicalGrid
    use SpeciesAndReferenceData

    implicit none

    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    type(MoleculeType), dimension(:), intent(in) :: molecule
    integer :: status, nspace

    call get_nspace( nspace )
    
    if( enable_bkw_analytic_vis )then

       if( num_species .gt. 1 .or. nspace .gt. 1 )then
          write(*,*) 'Error: Analytic solution to bkw function only available for single species, 0D.'
          stop
       end if

       call create_dist_func( bkw(1,1), vel_grid(1,1), molecule(1), num_species )

    end if

    ! Properties data structures
    allocate( total_data( 1:nspace, 1:number_properties ), STAT=status )
    call allocate_error_check( status, "total_data" )

    allocate( species_data( 1:num_species, 1:nspace, 1:number_properties ), STAT=status )
    call allocate_error_check( status, "species_data" )

    if( write_moment .gt. 0 )then
       allocate( total_moments( 1:nspace, 1:write_moment ), STAT=status )
       call allocate_error_check( status, "total_moments" )

       allocate( species_moments( 1:num_species, 1:nspace, 1:write_moment ), STAT=status )
       call allocate_error_check( status, "species_moments" )
    end if

    df_vis_unit            = base_unit
    rot_levels_vis_unit    = df_vis_unit + num_species*num_output_locations
    vib_levels_vis_unit    = rot_levels_vis_unit + num_species*num_output_locations
    rotvib_levels_vis_unit = vib_levels_vis_unit + num_species*num_output_locations
    properties_vis_unit    = rotvib_levels_vis_unit + num_species*num_output_locations
    collisions_vis_unit    = properties_vis_unit + num_species + 1
    bkw_analytic_vis_unit  = collisions_vis_unit + 1

    call write_headers()

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

    return
  end subroutine destroy_visualization

  subroutine write_headers()

    implicit none

    character(len=128) :: filename
    integer :: file_unit, n, i, loc, level
    integer :: rot_levels, vib_levels

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
                write( file_unit, advance="no", fmt="(a)" ) ', "Vel DistFunc"'
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
                  'variables = "Time", "Level Number", "Level Energy", "Density Fraction", "Energy"'

             close( file_unit )

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

    if( enable_rot_levels_vis .and. enable_vib_levels_vis )then

       do i = 1, num_species
          do loc = 1, num_output_locations

             file_unit = rotvib_levels_vis_unit + (i-1)*num_output_locations + (loc-1)
             filename = trim(adjustl(rotvib_levels_vis_filename))//"_"//&
                  trim(adjustl(int2char(i)))//"_"//&
                  trim(adjustl(int2char(output_locations(loc))))//".dat"

             open( unit=file_unit, file=filename, status="unknown" )

             write( file_unit, advance="yes", fmt="(a)" ) 'title = "Internal Energy Levels"'
             write( file_unit, advance="yes", fmt="(a)" )&
                  'variables = "Time", "Vib Level Number", "Rot Level Number", &
                  "Level Energy", "Density Fraction", "Energy"'

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

       close( file_unit )

    end if

    return
  end subroutine write_headers

  !===================================================================================================
  !  Output Data
  !===================================================================================================
  subroutine write_vis_data( time, ntime, phi, molecule, vel_grid, properties, num_collisions )

    use VelocityGrid
    use PhysicalProperties
    use InitialConditions

    implicit none

    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(VelocityGridType),dimension(:,:), intent(in) :: vel_grid
    type(PropertiesType), dimension(:), intent(in) :: properties
    double precision, intent(in) :: time
    integer, intent(in) :: ntime
    integer, dimension(:,:,:), intent(in) :: num_collisions

    double precision :: dens, temp, x_vel, mass

    if( enable_df_vis .and. mod(ntime,df_vis_dump_freq) .eq. 0 )&
         call write_df( time, phi, molecule, vel_grid )

    if( enable_rot_levels_vis .and. mod(ntime,rot_levels_vis_dump_freq) .eq. 0  )&
         call write_rot_levels( time, phi, molecule, vel_grid )

    if( enable_vib_levels_vis .and. mod(ntime,vib_levels_vis_dump_freq) .eq. 0  )&
         call write_vib_levels( time, phi, molecule, vel_grid )

    if( enable_rot_levels_vis .and. enable_vib_levels_vis .and. &
         mod(ntime,rotvib_levels_vis_dump_freq) .eq. 0 )&
         call write_rotvib_levels( time, phi, molecule, vel_grid )

    if( enable_properties_vis .and. mod(ntime,properties_vis_dump_freq) .eq. 0  )&
         call write_properties( time, phi, molecule, vel_grid, properties )

    if( enable_collisions_vis .and. mod(ntime,collisions_vis_dump_freq) .eq. 0  )&
         call write_collisions( time, num_collisions )
    
    if( enable_bkw_analytic_vis .and. enable_properties_vis .and. &
         mod(ntime,properties_vis_dump_freq) .eq. 0 )then

       dens  = properties(1)%dens(1)
       temp  = properties(1)%temp(1)
       x_vel = properties(1)%x_vel(1)
       mass  = molecule(1)%mass

       call generate_bkw( bkw(1,1), mass, dens, temp, x_vel, vel_grid(1,1), time )

       call write_bkw_properties( time, molecule, vel_grid, properties )

    end if

    return
  end subroutine write_vis_data

  subroutine write_df( time, phi, molecule, vel_grid )

    use VelocityGrid
    use PhysicalGrid

    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(MoleculeType), dimension(:), intent(in) :: molecule
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

                write( file_unit, advance="no", fmt="(3e24.16)" ) x, beta3, time

                if( enable_vel_df_vis )then
                   write( file_unit, advance="no", fmt="(e24.16)" ) phi(n,dfloc)%value(i,j,k)
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

                   write( file_unit, advance="no", fmt="(4e24.16)" ) x, y, beta3, time

                   if( enable_vel_df_vis )then
                      write( file_unit, advance="no", fmt="(e24.16)" ) phi(n,dfloc)%value(i,j,k)
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

                      write( file_unit, advance="no", fmt="(5e24.16)" ) x, y, z, beta3, time

                      if( enable_vel_df_vis )then
                         write( file_unit, advance="no", fmt="(e24.16)" ) phi(n,dfloc)%value(i,j,k)
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

  subroutine write_rot_levels( time, phi, molecule, vel_grid )

    use VelocityGrid
    use PhysicalGrid

    implicit none

    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    double precision, intent(in) :: time

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: modes, levels
    integer :: n, loc, i, j, k, l
    integer :: grid_ref, dfloc

    double precision :: energy_level, density_fraction

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

          do l = 1, levels
             energy_level = phi(n,dfloc)%rot_level(l)
             density_fraction = zero

             do k = k_min, k_max
                do j = j_min, j_max
                   do i = i_min, i_max
                      density_fraction = density_fraction + phi(n,dfloc)%rot(l,i,j,k)
                   end do
                end do
             end do
             
             write( file_unit, advance="yes", fmt="(e24.16,i3,3e24.16)" )&
                  time, (l-1), energy_level, density_fraction, density_fraction*energy_level
          end do

          close( file_unit )

       end do
    end do

    return
  end subroutine write_rot_levels

  subroutine write_vib_levels( time, phi, molecule, vel_grid )

    use VelocityGrid
    use PhysicalGrid

    implicit none

    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    double precision, intent(in) :: time

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: modes, levels
    integer :: n, loc, i, j, k, l
    integer :: grid_ref, dfloc

    double precision :: energy_level, density_fraction

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
                  time, (l-1), energy_level, density_fraction, density_fraction*energy_level
          end do

          close( file_unit )

       end do
    end do

    return
  end subroutine write_vib_levels

  subroutine write_rotvib_levels( time, phi, molecule, vel_grid )

    use VelocityGrid
    use PhysicalGrid

    implicit none

    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    double precision, intent(in) :: time

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: r_modes, v_modes, r_levels, v_levels
    integer :: n, loc, i, j, k, lr, lv
    integer :: grid_ref, dfloc

    double precision :: energy_level, density_fraction

    integer :: file_unit
    character(len=128) :: filename


    do n = 1, num_species

       r_modes = molecule(n)%rot_modes
       v_modes = molecule(n)%vib_modes

       do loc = 1, num_output_locations

          dfloc = output_locations(loc)

          file_unit = rotvib_levels_vis_unit + (n-1)*num_output_locations + (loc-1)
          filename = trim(adjustl(rotvib_levels_vis_filename))//"_"//&
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

          r_levels = phi(n,dfloc)%num_rot_levels
          v_levels = phi(n,dfloc)%num_vib_levels

          if( r_modes .eq. 0 .or. v_modes .eq. 0 ) cycle

          write( file_unit, advance="yes", fmt="(a7,i3,a4,i3)" )"zone i=",v_levels,", j=",r_levels

          do lr = 1, r_levels
             do lv = 1, v_levels
                energy_level = phi(n,dfloc)%int_energy_level(lr,lv)
                density_fraction = zero

                do k = k_min, k_max
                   do j = j_min, j_max
                      do i = i_min, i_max
                         density_fraction = density_fraction + phi(n,dfloc)%int_energy(lr,lv,i,j,k)
                      end do
                   end do
                end do

                write( file_unit, advance="yes", fmt="(e24.16,2i3,3e24.16)" )&
                     time, (lv-1), (lr-1), energy_level, density_fraction, density_fraction*energy_level
             end do
          end do

          close( file_unit )

       end do
    end do

    return
  end subroutine write_rotvib_levels

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

    integer :: nspace, np, ns, prop, i

    double precision :: delta_x

    integer :: file_unit
    character(len=128) :: filename

    call get_nspace( nspace )
    call get_delta_x( delta_x )

    call get_data_to_output( phi, molecule, vel_grid, properties )

    file_unit = properties_vis_unit
    filename = trim(adjustl(properties_vis_filename))//".dat"

    open( unit=file_unit, file=filename, status="old", position="append" )

    write( file_unit, advance="yes", fmt="(a7,i4)" ) "zone i=",nspace

    do np = 1, nspace
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

       write( file_unit, advance="yes", fmt="(a7,i4)" ) "zone i=",nspace

       do np = 1, nspace
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

    integer :: nspace, prop, i

    double precision :: delta_x

    integer :: file_unit
    character(len=128) :: filename

    call get_nspace( nspace )
    call get_delta_x( delta_x )

    call get_data_to_output( bkw, molecule, vel_grid, properties )

    file_unit = bkw_analytic_vis_unit
    filename = trim(adjustl(properties_vis_filename))//"_bkw.dat"

    open( unit=file_unit, file=filename, status="old", position="append" )
    
    write( file_unit, advance="yes", fmt="(a7,i4)" ) "zone i=",nspace

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
    integer :: np, nspace
    integer :: tot_collisions
    double precision :: delta_x

    integer :: file_unit
    character(len=128) :: filename

    call get_nspace( nspace )
    call get_delta_x( delta_x )

    file_unit = collisions_vis_unit
    filename = trim(adjustl(collisions_vis_filename))//".dat"

    open( unit=file_unit, file=filename, status="old", position="append" )

    write( file_unit, advance="yes", fmt="(a7,i4)" ) "zone i=", nspace

    do np = 1, nspace
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


  !===================================================================================================
  !  Calculate Data to Output
  !===================================================================================================
  subroutine get_data_to_output( phi, molecule, vel_grid, props )

    use VelocityGrid
    use PhysicalProperties
    use PhysicalGrid

    implicit none

    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    type(PropertiesType), dimension(:), intent(in) :: props

    integer :: np, ns, i, grid_ref
    integer :: nspace

    double precision :: prop_value

    call get_nspace( nspace )

    do np = 1, nspace

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
          do ns = 1, num_species
             species_data(ns,np,6) = zero
          end do
          total_data(np,6) = zero
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
             call compute_entropy( prop_value, props(np)%dens(ns), phi(ns,np), vel_grid(ns,grid_ref) )
             species_data(ns,np,18) = prop_value
          end do
          total_data(np,18) = zero
       end if

       ! Heat Flux X
       if( properties_to_output(19) )then
          do ns = 1, num_species
             call compute_heat_flux( prop_value, props(np)%x_vel(ns), props(np)%y_vel(ns), &
                  props(np)%z_vel(ns), x_dir, phi(ns,np), molecule(ns), vel_grid(ns,grid_ref) )
             species_data(ns,np,19) = prop_value
          end do
          total_data(np,19) = zero
       end if

       ! Heat Flux Y
       if( properties_to_output(20) )then
          do ns = 1, num_species
             call compute_heat_flux( prop_value, props(np)%x_vel(ns), props(np)%y_vel(ns), &
                  props(np)%z_vel(ns), y_dir, phi(ns,np), molecule(ns), vel_grid(ns,grid_ref) )
             species_data(ns,np,20) = prop_value
          end do
          total_data(np,20) = zero
       end if
       
       ! Heat Flux Z
       if( properties_to_output(21) )then
          do ns = 1, num_species
             call compute_heat_flux( prop_value, props(np)%x_vel(ns), props(np)%y_vel(ns), &
                  props(np)%z_vel(ns), z_dir, phi(ns,np), molecule(ns), vel_grid(ns,grid_ref) )
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

    end do

    return
  end subroutine get_data_to_output

end module Visualization
