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
  type(DistFuncType), allocatable, dimension(:,:,:) :: phi

  ! Temporary array used for updating phi
  type(DistFuncType), allocatable, dimension(:,:,:) :: delta_phi

  ! Molecular constants
  type(MoleculeType), allocatable, dimension(:) :: molecule

  ! Macroscopic properties at every spatial location
  type(PropertiesType), allocatable, dimension(:,:) :: properties

  ! The grid for velocity space
  type(VelocityGridType), allocatable, dimension(:) :: velocity_grid

  ! Shock conditions used for initialization
  type(NormalShock) :: shock_props

  ! Array for collision numbers
  integer, allocatable, dimension(:,:,:,:) :: colls_array

  ! Start time of simulation
  integer :: ntime_init

  ! Timing variables
  double precision :: tot_time, coll_time, conv_time, props_time
  double precision :: t0, t1, t2

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
    integer :: nx, ny, species
    integer :: nx_space, ny_space
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
    call get_nspace( nx_space, ny_space )
    
    allocate( velocity_grid( 1:num_species ), STAT=status )
    call allocate_error_check( status, "velocity_grid" )

    do species = 1, num_species
       call create_velocity_grid( velocity_grid( species ), species )
    end do

    ! Allocate the velocity distribution function
    allocate( phi( 1:num_species, 1:nx_space, 1:ny_space ), STAT=status )
    call allocate_error_check( status, "phi" )

    ! Allocate temporary velocity distribution function
    allocate( delta_phi( 1:num_species, 1:nx_space, 1:ny_space ), STAT=status )
    call allocate_error_check( status, "delta_phi" )

    ! Create velocity distribution functions
    do nx = 1, nx_space
       do ny = 1, ny_space
          do species = 1, num_species
             call create_dist_func( phi(species,nx,ny), velocity_grid(species), &
                  molecule(species), species )
             call create_dist_func( delta_phi(species,nx,ny), velocity_grid(species), &
                  molecule(species), species )
          end do
       end do
    end do

    ! Allocate and create macroscopic properties array
    allocate( properties( 1:nx_space, 1:ny_space ), STAT=status )
    call allocate_error_check( status, "properties" )

    do nx = 1, nx_space
       do ny = 1, ny_space
          call create_physical_properties( properties( nx, ny ) )
       end do
    end do

    ! Test restart and initialize velocity distribution functions
    call get_restart_flag( restart )
    
    ! Calculate normal shock properties - use domain of first species for test
    call initialize_normal_shock_props( shock_props, molecule )

    if( restart .eqv. .true. )then
       ! Read a restart file
!       call read_restart( ntime_init, phi, delta_phi, shock_props, molecule )

    else
       ntime_init = 0

       ! Initialize the distribution function
       do species = 1, num_species
          call initialize_dist_func( phi( species, :, : ), velocity_grid(species), &
               molecule(species), shock_props, species )
          call initialize_dist_func( delta_phi( species, :, : ), velocity_grid(species), &
               molecule(species), shock_props, species )
       end do

    end if

    ! Find curve fit for vibrational and rotational temperature
    do species = 1, num_species
       call find_rotational_temperature_relation( molecule(species), phi(species, 1, 1), species )
       call find_vibrational_temperature_relation( molecule(species), phi(species, 1, 1), species )
    end do

    ! Initialize visualization subroutines and data arrays
    call get_num_time_steps( num_time_steps )
    call initialize_visualization( velocity_grid, molecule, num_time_steps )

    ! Initialize collision integral subroutines and data arrays
    call create_collision_integral( velocity_grid, phi, molecule )

    ! Initialize the random number generator
    call initialize_rand()

    ! Initialize array for counting the number of collisions
    allocate( colls_array( 1:nx_space, 1:ny_space, 1:num_species, 1:num_species ), STAT=status )
    call allocate_error_check( status, "colls_array" )
    colls_array = 99999 ! preset value to number larger than cutoff number ( coll_limit )

    ! Destroy arrays used for initialization only
    call destroy_grid_input_arrays()

    ! TODO: temporary file creation
    if( nx_space .gt. 1 .and. ny_space .gt. 1 )then
       open( unit = 10, file = "left_wall.dat", status="unknown" )
       write(10,fmt="(a)")'title = "left wall fluxes"'
       write(10,fmt="(a)")'variables = "Time", "Location", "Heat Flux"'
       close( 10 )
       open( unit = 11, file = "right_wall.dat", status="unknown" )
       write(11,fmt="(a)")'title = "right wall fluxes"'
       write(11,fmt="(a)")'variables = "Time", "Location", "Heat Flux"'
       close( 11 )
       open( unit = 12, file = "bottom_wall.dat", status="unknown" )
       write(12,fmt="(a)")'title = "bottom wall fluxes"'
       write(12,fmt="(a)")'variables = "Time", "Location", "Heat Flux"'
       close( 12 )
       open( unit = 13, file = "top_wall.dat", status="unknown" )
       write(13,fmt="(a)")'title = "top wall fluxes"'
       write(13,fmt="(a)")'variables = "Time", "Location", "Heat Flux"'
       close( 13 )
    end if

    return
  end subroutine initialize_solver

  subroutine solve()

    implicit none

    integer :: ntime, num_time_steps
    integer :: species

    integer :: nx_space, ny_space, nx, ny

    double precision :: deltax, deltay
    double precision :: deltat, time
    double precision :: dens

    logical :: convection_flag

    ! Get flags
    call get_convection_flag( convection_flag )

    ! Get spatial data
    call get_delta_x( deltax, deltay )
    call get_nspace( nx_space, ny_space )

    ! Calculate convection time steps
    !TODO
    if( convection_flag )then
       do species = 1, num_species
          call calculate_t_conv( deltax, velocity_grid(species)%xv_min, velocity_grid(species)%xv_max )
       end do
    end if

    ! Get time step data
    call get_num_time_steps( num_time_steps )
    call get_deltat( deltat )

    ntime = ntime_init
    if( ntime_init .ge. num_time_steps )then
       write(*,*) "Error: Start time is after/equal to end time. ntime, max_time = ", ntime, num_time_steps
       stop
    end if

    time = dble(ntime) * deltat


    ! Pre-calculate macroscopic properties
    do nx = 1, nx_space
       do ny = 1, ny_space
          call calculate_physical_properties( properties(nx,ny), phi(:,nx,ny), molecule, velocity_grid )
       end do
    end do

    !dens = 0
    !do nx = 1, nx_space
    !   dens = dens + properties(nx)%mix_dens
    !end do
    !print *, "Density at time step zero is ", dens

    ! Print initial macroscopic properties to screen for time zero
    call print_property_update( ntime, num_time_steps, properties, colls_array )

    ! Write intial vis data for time zero
    call write_vis_data( time, ntime, num_time_steps, phi, molecule, velocity_grid, properties, colls_array )

    ! Start and initialize timers
    call cpu_time( t0 )
    coll_time = zero
    conv_time = zero

    ! Begin time loop
    do ntime = ntime_init+1, num_time_steps
       
       ! Find time
       time = dble(ntime) * deltat

       ! Collision integral
       do nx = 1, nx_space
          do ny = 1, ny_space
             ! Start collision timer
             call cpu_time( t1 )
             ! Calculate changes in the distribution function due to collisions
             call compute_collision_integral( phi(:,nx,ny), delta_phi(:,nx,ny), molecule, &
                  velocity_grid, properties(nx,ny), nx, ny, colls_array )
             call cpu_time( t2 )
             coll_time = coll_time + ( t2 - t1 )
          end do
       end do

       ! Convection
       if( convection_flag .eqv. .true. )then
          do species = 1, num_species
             ! Start convection timer
             call cpu_time( t1 )
             ! Calculate changes in the distribution function due to convection
             call compute_convection( phi(species,:,:), delta_phi(species,:,:), &
                  molecule(species), velocity_grid(species), species, deltat, &
                  properties, ntime )
             call cpu_time( t2 )
             conv_time = conv_time + ( t2 - t1 )
          end do
       end if

       ! Update the macroscopic properties
       ! Start calculate properties timer
       call cpu_time( t1 )
       do nx = 1, nx_space
          do ny = 1, ny_space
             call calculate_physical_properties( properties(nx,ny), phi(:,nx,ny), molecule, velocity_grid )
          end do
       end do
       call cpu_time( t2 )
       props_time = props_time + ( t2 - t1 )

       !dens = 0
       !do nx = 1, nx_space
       !   dens = dens + properties(nx)%mix_dens
       !end do
       !print *, "Density at time step ", ntime, " is ", dens

       ! Output the macroscopic properties to the screen
       call print_property_update( ntime, num_time_steps, properties, colls_array )

       ! Write vis data
       call write_vis_data( time, ntime, num_time_steps, phi, molecule, velocity_grid, properties, colls_array )

       ! Write restart
!       call write_restart( ntime, num_time_steps, phi, molecule )

    end do ! End of time step loop

    ! End timer
    call cpu_time( tot_time )

    tot_time = tot_time - t0
    print *, "Total time taken was: ",tot_time
    print *, "Time taken by collisions was: ",coll_time
    print *, "Time taken by convection was: ",conv_time
    print *, "Time taken by calculate properties was: ",props_time

    return
  end subroutine solve

  subroutine destroy_solver()

    implicit none

    integer :: nx, ny, species
    integer :: nx_space, ny_space
    integer :: status

    call get_nspace( nx_space, ny_space )

    ! Destroy visualization
    call destroy_visualization()

    ! Destroy collision integral related data allocations
    call destroy_collision_integral( molecule )

    ! Destroy velocity distribution functions
    do species = 1, num_species
       do nx = 1, nx_space
          do ny = 1, ny_space
             call destroy_dist_func( phi( species, nx, ny ), molecule(species) )
             call destroy_dist_func( delta_phi( species, nx, ny ), molecule(species) )
          end do
       end do
    end do

    deallocate( phi, STAT=status )
    deallocate( delta_phi, STAT=status )

    ! Destroy velocity grid data structures
    do species = 1, num_species
       call destroy_velocity_grid( velocity_grid( species ) )
    end do

    deallocate( velocity_grid, STAT=status )

    ! Destroy molecule data structure
    deallocate( molecule, STAT=status )

    ! Destroy random number generation related data allocations
    call destroy_rand()

    return
  end subroutine destroy_solver

end module DVMSolver
