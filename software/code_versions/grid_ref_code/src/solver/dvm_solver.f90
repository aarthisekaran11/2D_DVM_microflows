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
module DVMSolver

  use SpeciesAndReferenceData
  use DistFunc
  use VelocityGrid
  use PhysicalGrid
  use PhysicalProperties
  use TimeStepping
  use InitialConditions
  use ShockConditions
  use Visualization
  use ScreenOutput
  use RandomNumberGeneration
  use CollisionIntegral
  use CollisionUtilities
  use Convection
  use Constants
  use ErrorCheck
  use Restart
  
  implicit none

  private

  !=====================================================================
  ! Private data here
  !=====================================================================

  ! THE velocity distribution function. 
  ! Currently 0-D or 1-D in space. Defined for each species
  type(DistFuncType), allocatable, dimension(:,:) :: phi

  ! Ttemporary array used for updating phi
  type(DistFuncType), allocatable, dimension(:,:) :: delta_phi

  ! Molecular constants
  type(MoleculeType), allocatable, dimension(:) :: molecule

  ! Macroscopic properties at every spatial location
  type(PropertiesType), allocatable, dimension(:) :: properties

  ! The grid for velocity space
  type(VelocityGridType), allocatable, dimension(:,:) :: velocity_grid

  ! Shock conditions used for initialization
  type(NormalShock) :: shock_props

  ! Array for collision numbers
  integer, allocatable, dimension(:,:,:) :: colls_array

  ! Start time of simulation
  integer :: ntime_init

  !=====================================================================
  ! Declare public subroutines here
  !=====================================================================

  public :: initialize_solver
  public :: solve
  public :: destroy_solver

contains

  subroutine initialize_solver()

    use CommandLineProcessor
    use ReadInput
    
    implicit none

    integer :: status
    integer :: element, space, species, grid_ref
    integer :: num_grid_types, nx_space, ny_space
    integer :: num_time_steps
    logical :: restart

    ! Parse the command line
    call process_command_line_args()

    ! Read the input file
    call read_input()

    ! Screen output heading
    call title_message()
    call init_message()

    ! Allocate and set the species constants
    allocate( molecule( 1:num_species ), STAT=status )
    call allocate_error_check( status, "molecule" )

    call set_molecule_constants( molecule, num_species )
   
    ! Allocate and create the velocity grid
    call get_num_zones( num_grid_types )
    call get_nspace( nx_space, ny_space )
    
    allocate( velocity_grid( 1:num_species, 1:num_grid_types ), STAT=status )
    call allocate_error_check( status, "velocity_grid" )

    do element = 1, num_grid_types
       do species = 1, num_species
          call create_velocity_grid( velocity_grid( species, element ), species, element )
       end do
    end do

    ! Allocate the velocity distribution function
    allocate( phi( 1:num_species, 1:nx_space ), STAT=status )
    call allocate_error_check( status, "phi" )

    ! Allocate temporary velocity distribution function
    allocate( delta_phi( 1:num_species, 1:nx_space ), STAT=status )
    call allocate_error_check( status, "delta_phi" )

    ! Create velocity distribution functions
    do space = 1, nx_space
       call get_spatial_reference( space, grid_ref )

       do species = 1, num_species
          call create_dist_func( phi(species,space), velocity_grid(species,grid_ref), &
               molecule(species), species )
          call create_dist_func( delta_phi(species,space), velocity_grid(species,grid_ref), &
               molecule(species), species )
       end do

    end do

    ! Allocate and create macroscopic properties array
    allocate( properties( 1:nx_space ), STAT=status )
    call allocate_error_check( status, "properties" )

    do space = 1, nx_space
       call create_physical_properties( properties( space ) )
    end do

    ! Test restart and initialize velocity distribution functions
    call get_restart_flag( restart )
    
    ! Calculate normal shock properties - use domain of first species for test
    call initialize_normal_shock_props( shock_props, molecule )

    if( restart .eqv. .true. )then
       ! Read a restart file
       call read_restart( ntime_init, phi, delta_phi, shock_props, molecule )

    else
       ntime_init = 0

       ! Initialize the distribution function
       do species = 1, num_species
          call initialize_dist_func( phi( species, : ), velocity_grid( species, : ), &
               molecule(species), shock_props, species )
          call initialize_dist_func( delta_phi( species, : ), velocity_grid( species, : ), &
               molecule(species), shock_props, species )
       end do

    end if

    ! Find curve fit for vibrational and rotational temperature
    do species = 1, num_species
       call find_rotational_temperature_relation( molecule(species), phi(species, 1), species )
       call find_vibrational_temperature_relation( molecule(species), phi(species, 1), species )
    end do     

    ! Initialize visualization subroutines and data arrays
    call get_num_time_steps( num_time_steps )
    call initialize_visualization( velocity_grid, molecule, num_time_steps )

    ! Initialize collision integral subroutines and data arrays
    call create_collision_integral( velocity_grid, phi, molecule )

    ! Initialize the random number generator
    call initialize_rand()

    ! Initialize array for counting the number of collisions
    allocate( colls_array( 1:nx_space, 1:num_species, 1:num_species ), STAT=status )
    call allocate_error_check( status, "colls_array" )
    colls_array = 99999 ! preset value to number larger than cutoff number ( coll_limit )

    ! Destroy arrays used for initialization only
    call destroy_grid_input_arrays()

    return
  end subroutine initialize_solver

  subroutine solve()

    implicit none

    integer :: ntime, num_time_steps
    integer :: species, grid_ref

    integer :: nx_space, ny_space, nx

    double precision :: deltat, time

    logical :: convection_flag

    ! Get flags
    call get_convection_flag( convection_flag )

    ! Get time step data
    call get_num_time_steps( num_time_steps )
    call get_deltat( deltat )

    ntime = ntime_init
    if( ntime_init .ge. num_time_steps )then
       write(*,*) "Error: Start time is after/equal to end time. ntime, max_time = ", ntime, num_time_steps
       stop
    end if

    time = dble(ntime) * deltat

    ! Get spatial data
    call get_nspace( nx_space, ny_space )

    ! Pre-calculate macroscopic properties
    do nx = 1, nx_space
       call get_spatial_reference( nx, grid_ref )
       call calculate_physical_properties( properties(nx), phi(:,nx), molecule, velocity_grid(:,grid_ref) )
    end do

    ! Print initial macroscopic properties to screen for time zero
    call print_property_update( ntime, num_time_steps, properties, colls_array )

    ! Write intial vis data for time zero
    call write_vis_data( time, ntime, num_time_steps, phi, molecule, velocity_grid, properties, colls_array )

    ! Begin time loop
    do ntime = ntime_init+1, num_time_steps
       
       ! Find time
       time = dble(ntime) * deltat

       ! Collision integral
       do nx = 1, nx_space
          call get_spatial_reference( nx, grid_ref )

          ! Calculate changes in the distribution function due to collisions
          call compute_collision_integral( phi(:,nx), delta_phi(:,nx), molecule, &
               velocity_grid(:,grid_ref), properties(nx), nx, colls_array, ntime )
       end do

       ! Convection
       if( convection_flag .eqv. .true. )then
          do species = 1, num_species
             ! Calculate changes in the distribution function due to convection
             call compute_convection( phi(species,:), delta_phi(species,:), &
                  molecule(species), velocity_grid(species,:), species, deltat )
          end do
       end if

       ! Update the macroscopic properties
       do nx = 1, nx_space
          call get_spatial_reference( nx, grid_ref )
          call calculate_physical_properties( properties(nx), phi(:,nx), molecule, velocity_grid(:,grid_ref) )
       end do

       ! Output the macroscopic properties to the screen
       call print_property_update( ntime, num_time_steps, properties, colls_array )

       ! Write vis data
       call write_vis_data( time, ntime, num_time_steps, phi, molecule, velocity_grid, properties, colls_array )

       ! Write restart
       call write_restart( ntime, num_time_steps, phi, molecule )

    end do ! End of time step loop

    return
  end subroutine solve

  subroutine destroy_solver()

    implicit none

    integer :: space, species, grid_ref
    integer :: num_grid_types, nx_space, ny_space
    integer :: status

    call get_nspace( nx_space, ny_space )

    ! Destroy visualization
    call destroy_visualization()

    ! Destroy collision integral related data allocations
    call destroy_collision_integral( molecule )

    ! Destroy velocity distribution functions
    do species = 1, num_species
       do space = 1, nx_space
          call destroy_dist_func( phi( species, space ), molecule(species) )
          call destroy_dist_func( delta_phi( species, space ), molecule(species) )
       end do
    end do

    deallocate( phi, STAT=status )
    deallocate( delta_phi, STAT=status )

    ! Destroy velocity grid data structures
    call get_num_zones( num_grid_types )
    
    do grid_ref = 1, num_grid_types
       do species = 1, num_species
          call destroy_velocity_grid( velocity_grid( species, grid_ref ) )
       end do
    end do

    deallocate( velocity_grid, STAT=status )

    ! Destroy molecule data structure
    deallocate( molecule, STAT=status )

    ! Destroy random number generation related data allocations
    call destroy_rand()

    return
  end subroutine destroy_solver

end module DVMSolver
