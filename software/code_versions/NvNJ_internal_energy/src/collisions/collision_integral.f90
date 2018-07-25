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
module CollisionIntegral

  use ErrorCheck

  implicit none

  private

  integer :: collision_integral_solver = -1

  integer, parameter :: no_collisions      = 0
  integer, parameter :: n_squared          = 1
  integer, parameter :: full_monte_carlo   = 2
  integer, parameter :: variance_reduction = 3

  double precision, allocatable, dimension(:,:) :: coln_rms

  integer :: recalculate_equil_limit

  integer, parameter :: elastic   = 0
  integer, parameter :: rotation  = 0
  integer, parameter :: vibration = 0

  public :: set_collision_integral_solver
  public :: set_recalculate_equil_limit
  public :: initialize_coln_rms
  public :: create_collision_integral
  public :: destroy_collision_integral
  public :: compute_collision_integral

contains

  subroutine set_collision_integral_solver( coll_integral_solver_in )

    implicit none

    integer, intent(in) :: coll_integral_solver_in

    collision_integral_solver = coll_integral_solver_in

    return
  end subroutine set_collision_integral_solver

  subroutine set_recalculate_equil_limit( recalculate_equil_limit_in )
    
    implicit none

    integer, intent(in) :: recalculate_equil_limit_in

    recalculate_equil_limit = recalculate_equil_limit_in

    return
  end subroutine set_recalculate_equil_limit

  subroutine initialize_coln_rms( coln_rms_in )

    use DistFunc
    
    double precision, dimension(:), intent(in) :: coln_rms_in
    integer :: n, m, status

    allocate( coln_rms( 1:num_species, 1:num_species ), STAT=status )
    call allocate_error_check( status, "coln_rms" )

    do n = 1, num_species
       do m = 1, num_species
          coln_rms(n,m) = coln_rms_in( num_species*(n-1) + m )
       end do
    end do

    return
  end subroutine initialize_coln_rms

  subroutine create_collision_integral( vel_grid, phi, molecule )

    use VelocityGrid
    use ReplenishingCollisions
    use FullMonteCarlo
    use VarianceReduction
    use DistFunc
    use SpeciesAndReferenceData

    implicit none

    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(VelocityGridType), dimension(:,:), intent(in) :: vel_grid
    type(MoleculeType), dimension(:), intent(in) :: molecule

    select case( collision_integral_solver )
    case( no_collisions )
       ! Nothing to create for this case

    case( full_monte_carlo )
       call create_full_monte_carlo( vel_grid, molecule )

    case( variance_reduction )
       call create_variance_reduction( vel_grid, phi, molecule )

    case default
       write(*,*) "Error: Invalid value of collision_integral_solver: ", collision_integral_solver
       stop

    end select

    return
  end subroutine create_collision_integral

  subroutine destroy_collision_integral( molecule )

    use ReplenishingCollisions
    use FullMonteCarlo
    use VarianceReduction
    use DistFunc
    use SpeciesAndReferenceData

    implicit none

    type(MoleculeType), dimension(:), intent(in) :: molecule

    select case( collision_integral_solver )
    case( no_collisions )
       ! Nothing to destroy for this case

    case( full_monte_carlo )
       call destroy_full_monte_carlo( molecule )

    case( variance_reduction )
       call destroy_variance_reduction( molecule )

    case default
       write(*,*) "Error: Invalid value of collision_integral_solver: ", collision_integral_solver
       stop

    end select

    return
  end subroutine destroy_collision_integral

  subroutine compute_collision_integral( phi, delta_phi, molecule, vel_grid, properties, ns, colls_array )

    use DistFunc
    use VelocityGrid
    use PhysicalProperties
    use FullMonteCarlo
    use VarianceReduction
    use SpeciesAndReferenceData
    
    implicit none

    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(PropertiesType), intent(in) :: properties
    integer, intent(in) :: ns

    type(DistFuncType), dimension(:) :: phi, delta_phi
    integer, dimension(:,:,:) :: colls_array

    integer, parameter :: AA = 1
    integer, parameter :: AB = 2

    integer :: n, m
    integer :: colls
    integer :: energy_flag

    logical :: new_eq_flag

    energy_flag = vibration

    do n = 1, num_species
       call perform_dist_func_update( delta_phi(n), phi(n), molecule(n), vel_grid(n) )
    end do

    ! Construct necessary cumulative distribution functions
    select case( collision_integral_solver )
    case( no_collisions )
       ! Nothing to construct for this case

    case( full_monte_carlo )
       do n = 1, num_species
          call compute_psi_fmc( phi(n), molecule(n), vel_grid(n), n )
       end do

    case( variance_reduction )
       do n = 1, num_species

          new_eq_flag = .false.
          if( colls_array(ns,n,n) .gt. recalculate_equil_limit ) new_eq_flag = .true.

          ! Like species collisions
          call compute_psi_vr( phi(n), molecule(n), vel_grid(n), properties, new_eq_flag, &
               energy_flag, AA, n, ns )

          new_eq_flag = .false.
          do m = 1, num_species
             if( m .eq. n )cycle
             if( colls_array(ns,n,m) .gt. recalculate_equil_limit ) new_eq_flag = .true.
          end do

          ! Unlike species collisions
          call compute_psi_vr( phi(n), molecule(n), vel_grid(n), properties, new_eq_flag, &
               energy_flag, AB, n, ns )

       end do

    case default
       write(*,*) "Error: Invalid value of collision_integral_solver: ", collision_integral_solver
       stop

    end select

    ! Solve collisions integral
    do n = 1, num_species
       do m = n, num_species

          select case( collision_integral_solver )
          case( no_collisions )
             ! No collisions, just move on to convection

          case( full_monte_carlo )
             call full_monte_carlo_kernel( delta_phi, colls, n, m, coln_rms, molecule, vel_grid, properties )

          case( variance_reduction )
             call variance_reduction_kernel( delta_phi, colls, n, m, ns, coln_rms, &
                  molecule, vel_grid, properties )

          case default
             write(*,*) "Error: Invalid value of collision_integral_solver: ", collision_integral_solver
             stop

          end select

          colls_array(ns,n,m) = colls
          colls_array(ns,m,n) = colls

       end do
    end do

    do n = 1, num_species
       call perform_dist_func_update( phi(n), delta_phi(n), molecule(n), vel_grid(n) )
    end do

    return
  end subroutine compute_collision_integral

end module CollisionIntegral
