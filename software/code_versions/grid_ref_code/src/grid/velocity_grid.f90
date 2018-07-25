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
module VelocityGrid

  use ErrorCheck

  implicit none

  private

  type VelocityGridType
     integer :: num_points_x, num_points_y, num_points_z
     integer :: num_points

     integer :: i_min, i_max
     integer :: j_min, j_max
     integer :: k_min, k_max

     integer :: i_neg, i_pos
     integer :: j_neg, j_pos
     integer :: k_neg, k_pos

     integer :: i_zero, j_zero, k_zero

     double precision :: xv_min, xv_max
     double precision :: yv_min, yv_max
     double precision :: zv_min, zv_max

     double precision :: beta3_min, beta3_max, beta3_avg

     double precision, allocatable, dimension(:) :: x
     double precision, allocatable, dimension(:) :: y
     double precision, allocatable, dimension(:) :: z
     double precision, allocatable, dimension(:) :: beta_x
     double precision, allocatable, dimension(:) :: beta_y
     double precision, allocatable, dimension(:) :: beta_z

     ! Precalculated coefficients and values used in remapping scheme
     double precision :: xv_factor, yv_factor, zv_factor
     logical, dimension(3) :: uniform_grid
     double precision :: uniform_divisor, ud_xy, ud_xz, ud_yz, ud_x, ud_y, ud_z
     double precision :: ud_ix, ud_iy, ud_iz
     
  end type VelocityGridType

  integer, allocatable, dimension(:) :: x_vel_grid_type
  integer, allocatable, dimension(:) :: y_vel_grid_type
  integer, allocatable, dimension(:) :: z_vel_grid_type

  type(VelocityGridType), allocatable, dimension(:) :: inp_vel_grid

  integer, parameter :: file_read_grid = 0
  integer, parameter :: uniform_grid   = 1

  integer, parameter :: x_dir = 1
  integer, parameter :: y_dir = 2
  integer, parameter :: z_dir = 3

  integer :: num_species_grid

  ! Public variables
  public :: VelocityGridType

  ! Public subroutines
  public :: initialize_grid_input_arrays
  public :: initialize_file_read_arrays
  public :: destroy_grid_input_arrays
  public :: destroy_file_read_arrays
  public :: set_vel_grid_type
  public :: set_vel_grid_bounds
  public :: set_vel_grid_num_points
  public :: set_velocity_vectors
  public :: create_velocity_grid
  public :: destroy_velocity_grid
  public :: global2local_map
  public :: local2global_map

contains

  subroutine initialize_grid_input_arrays( num_elems )
    
    implicit none

    integer, intent(in) :: num_elems
    integer :: status

    num_species_grid = num_elems

    allocate( x_vel_grid_type( 1:num_elems ), STAT=status )
    call allocate_error_check( status, "x_vel_grid_type" )
    allocate( y_vel_grid_type( 1:num_elems ), STAT=status )
    call allocate_error_check( status, "y_vel_grid_type" )
    allocate( z_vel_grid_type( 1:num_elems ), STAT=status )
    call allocate_error_check( status, "z_vel_grid_type" )

    allocate( inp_vel_grid( 1:num_elems ), STAT=status )
    call allocate_error_check( status, "inp_vel_grid" )

    return
  end subroutine initialize_grid_input_arrays

  subroutine initialize_file_read_arrays( species )

    implicit none

    integer, intent(in) :: species
    integer :: status
    integer :: num_points_x, num_points_y, num_points_z

    num_points_x = inp_vel_grid(species)%num_points_x
    num_points_y = inp_vel_grid(species)%num_points_y
    num_points_z = inp_vel_grid(species)%num_points_z

    allocate( inp_vel_grid(species)%x(1:num_points_x), STAT=status )
    call allocate_error_check( status, "inp_vel_grid%x" )
    
    allocate( inp_vel_grid(species)%y(1:num_points_y), STAT=status )
    call allocate_error_check( status, "inp_vel_grid%y" )
    
    allocate( inp_vel_grid(species)%z(1:num_points_z), STAT=status )
    call allocate_error_check( status, "inp_vel_grid%z" )

    return
  end subroutine initialize_file_read_arrays

  subroutine destroy_grid_input_arrays()

    implicit none

    integer :: status

    deallocate( x_vel_grid_type, STAT=status )
    call deallocate_error_check( status, "x_vel_grid_type" )
    deallocate( y_vel_grid_type, STAT=status )
    call deallocate_error_check( status, "y_vel_grid_type" )
    deallocate( z_vel_grid_type, STAT=status )
    call deallocate_error_check( status, "z_vel_grid_type" )

    deallocate( inp_vel_grid, STAT=status )
    call deallocate_error_check( status, "inp_vel_grid" )

    return
  end subroutine destroy_grid_input_arrays

  subroutine destroy_file_read_arrays( num_elems )

    implicit none
    
    integer, intent(in) :: num_elems
    integer :: status, n

    do n = 1, num_elems
       deallocate( inp_vel_grid(n)%x, STAT=status )
       call deallocate_error_check( status, "inp_vel_grid%x" )

       deallocate( inp_vel_grid(n)%y, STAT=status )
       call deallocate_error_check( status, "inp_vel_grid%y" )

       deallocate( inp_vel_grid(n)%z, STAT=status )
       call deallocate_error_check( status, "inp_vel_grid%z" )
    end do

    return
  end subroutine destroy_file_read_arrays

  subroutine set_vel_grid_type( x_flag_in, y_flag_in, z_flag_in, element )

    use Conversion

    implicit none

    integer, intent(in) :: x_flag_in, y_flag_in, z_flag_in, element

    x_vel_grid_type( element ) = x_flag_in
    y_vel_grid_type( element ) = y_flag_in
    z_vel_grid_type( element ) = z_flag_in

    return
  end subroutine set_vel_grid_type

  subroutine set_vel_grid_bounds( xv_min_in, xv_max_in, yv_min_in, &
       yv_max_in, zv_min_in, zv_max_in, element )
    
    implicit none

    double precision, intent(in) :: xv_min_in, xv_max_in
    double precision, intent(in) :: yv_min_in, yv_max_in
    double precision, intent(in) :: zv_min_in, zv_max_in
    integer, intent(in) :: element

    inp_vel_grid(element)%xv_min = xv_min_in
    inp_vel_grid(element)%xv_max = xv_max_in
    inp_vel_grid(element)%yv_min = yv_min_in
    inp_vel_grid(element)%yv_max = yv_max_in
    inp_vel_grid(element)%zv_min = zv_min_in
    inp_vel_grid(element)%zv_max = zv_max_in

    return
  end subroutine set_vel_grid_bounds

  subroutine set_vel_grid_num_points( num_points_x_in, num_points_y_in, num_points_z_in, element )

    implicit none

    integer, intent(in) :: num_points_x_in, num_points_y_in, num_points_z_in
    integer, intent(in) :: element

    inp_vel_grid(element)%num_points_x = num_points_x_in
    inp_vel_grid(element)%num_points_y = num_points_y_in
    inp_vel_grid(element)%num_points_z = num_points_z_in

    return
  end subroutine set_vel_grid_num_points

  subroutine set_velocity_vectors( x_in, y_in, z_in, species )

    implicit none

    double precision, dimension(:), intent(in) :: x_in, y_in, z_in
    integer, intent(in) :: species

    integer :: max
    
    inp_vel_grid(species)%x = x_in
    inp_vel_grid(species)%y = y_in
    inp_vel_grid(species)%z = z_in

    inp_vel_grid(species)%xv_min = x_in(1)
    inp_vel_grid(species)%yv_min = y_in(1)
    inp_vel_grid(species)%zv_min = z_in(1)

    max = inp_vel_grid(species)%num_points_x
    inp_vel_grid(species)%xv_max = x_in(max)

    max = inp_vel_grid(species)%num_points_y
    inp_vel_grid(species)%yv_max = y_in(max)

    max = inp_vel_grid(species)%num_points_z
    inp_vel_grid(species)%zv_max = z_in(max)

    return
  end subroutine set_velocity_vectors

  subroutine create_velocity_grid( grid, species, grid_ref )

    implicit none

    type(VelocityGridType) :: grid

    integer, intent(in) :: species, grid_ref

    integer :: status

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: index

    ! location of grid dimensions in input array
    index = num_species_grid*( grid_ref -1 ) + species

    ! Number of grid points in each direction
    grid%num_points_x = inp_vel_grid(index)%num_points_x
    grid%num_points_y = inp_vel_grid(index)%num_points_y
    grid%num_points_z = inp_vel_grid(index)%num_points_z
    grid%num_points   = grid%num_points_x*grid%num_points_y*grid%num_points_z

    ! Velocity and index bounds
    grid%xv_min = inp_vel_grid(index)%xv_min
    grid%i_min  = 0

    grid%xv_max = inp_vel_grid(index)%xv_max
    grid%i_max  = grid%num_points_x - 1

    grid%xv_factor = dble( grid%num_points_x - 1 )/ ( grid%xv_max - grid%xv_min )

    grid%yv_min = inp_vel_grid(index)%yv_min
    grid%j_min  = 0

    grid%yv_max = inp_vel_grid(index)%yv_max
    grid%j_max  = grid%num_points_y - 1

    grid%yv_factor = dble( grid%num_points_y - 1 ) / ( grid%yv_max - grid%yv_min )

    grid%zv_min = inp_vel_grid(index)%zv_min
    grid%k_min  = 0

    grid%zv_max = inp_vel_grid(index)%zv_max
    grid%k_max  = grid%num_points_z - 1

    grid%zv_factor = dble( grid%num_points_z - 1 ) / ( grid%zv_max - grid%zv_min )

    i_min = grid%i_min
    i_max = grid%i_max
    j_min = grid%j_min
    j_max = grid%j_max
    k_min = grid%k_min
    k_max = grid%k_max

    allocate(grid%x(i_min:i_max), STAT=status)
    call allocate_error_check( status, "grid%x" )

    allocate(grid%y(j_min:j_max), STAT=status)
    call allocate_error_check( status, "grid%y" )

    allocate(grid%z(k_min:k_max), STAT=status)
    call allocate_error_check( status, "grid%z" )

    allocate(grid%beta_x(i_min:i_max), STAT=status)
    call allocate_error_check( status, "grid%beta_x" )

    allocate(grid%beta_y(j_min:j_max), STAT=status)
    call allocate_error_check( status, "grid%beta_y" )

    allocate(grid%beta_z(k_min:k_max), STAT=status)
    call allocate_error_check( status, "grid%beta_z" )

    ! Find the velocity and grid spacing values of the velocity grid
    call generate_grid( grid, index )

    return
  end subroutine create_velocity_grid

  subroutine destroy_velocity_grid( grid )

    use ErrorCheck

    implicit none

    type(VelocityGridType) :: grid

    integer :: status

    deallocate(grid%x, STAT=status)
    call deallocate_error_check(status, "grid%x")

    deallocate(grid%y, STAT=status)
    call deallocate_error_check(status, "grid%y")

    deallocate(grid%z, STAT=status)
    call deallocate_error_check(status, "grid%z")

    deallocate(grid%beta_x, STAT=status)
    call deallocate_error_check(status, "grid%beta_x")

    deallocate(grid%beta_y, STAT=status)
    call deallocate_error_check(status, "grid%beta_y")

    deallocate(grid%beta_z, STAT=status)
    call deallocate_error_check(status, "grid%beta_z")

    return
  end subroutine destroy_velocity_grid

  subroutine generate_grid( grid, index )

    use Constants

    implicit none

    integer, intent(in) :: index
    type(VelocityGridType) :: grid

    integer :: min, max, num_points
    
    double precision :: v_min, v_max

    ! X-velocity vector
    min = grid%i_min
    max = grid%i_max
    v_min = grid%xv_min
    v_max = grid%xv_max
    num_points = grid%num_points_x

    select case( x_vel_grid_type( index ) )
    case( file_read_grid )
       call generate_file_read_grid( grid%x, index, min, max, x_dir )

    case( uniform_grid )
       grid%uniform_grid(1) = .true.
       call generate_uniform_grid( grid%x, min, max, num_points, v_min, v_max )

    case default
       write(*,*) "Error: Invalid value of x_vel_grid_type. Given value is: ", x_vel_grid_type(index)
       stop

    end select

    call find_zero_point( grid%i_zero, grid%x, min, max )

    ! Y-velocity vector
    min = grid%j_min
    max = grid%j_max
    v_min = grid%yv_min
    v_max = grid%yv_max
    num_points = grid%num_points_y

    select case( y_vel_grid_type( index ) )
    case( file_read_grid )
       call generate_file_read_grid( grid%y, index, min, max, y_dir )

    case( uniform_grid )
       grid%uniform_grid(2) = .true.
       call generate_uniform_grid( grid%y, min, max, num_points, v_min, v_max )

    case default
       write(*,*) "Error: Invalid value of y_vel_grid_type. Given value is: ", y_vel_grid_type(index)
       stop

    end select

    call find_zero_point( grid%j_zero, grid%y, min, max )

    ! Z-velocity vector
    min = grid%k_min
    max = grid%k_max
    v_min = grid%zv_min
    v_max = grid%zv_max
    num_points = grid%num_points_z

    select case( z_vel_grid_type( index ) )
    case( file_read_grid )
       call generate_file_read_grid( grid%z, index, min, max, z_dir )

    case( uniform_grid )
       grid%uniform_grid(3) = .true.
       call generate_uniform_grid( grid%z, min, max, num_points, v_min, v_max )

    case default
       write(*,*) "Error: Invalid value of z_vel_grid_type. Given value is: ", z_vel_grid_type(index)
       stop

    end select

    call find_zero_point( grid%k_zero, grid%z, min, max )

    ! Calculate beta values
    call compute_beta( grid )

    call find_min_max_beta3( grid )

    if( all(grid%uniform_grid) )then
       grid%uniform_divisor = one / ( two * ( grid%beta_x(min) * grid%beta_x(min) + &
            grid%beta_y(min) * grid%beta_y(min) + grid%beta_z(min) * grid%beta_z(min) ) )

       grid%ud_xy = one / ( two * ( grid%beta_x(min) * grid%beta_x(min) + grid%beta_y(min) * grid%beta_y(min) ) )
       grid%ud_xz = one / ( two * ( grid%beta_x(min) * grid%beta_x(min) + grid%beta_z(min) * grid%beta_z(min) ) )
       grid%ud_yz = one / ( two * ( grid%beta_y(min) * grid%beta_y(min) + grid%beta_z(min) * grid%beta_z(min) ) )
       grid%ud_x  = one / ( two * ( grid%beta_x(min) * grid%beta_x(min) ) )
       grid%ud_y  = one / ( two * ( grid%beta_y(min) * grid%beta_y(min) ) )
       grid%ud_z  = one / ( two * ( grid%beta_z(min) * grid%beta_z(min) ) )
       grid%ud_ix = one / grid%beta_x(min)
       grid%ud_iy = one / grid%beta_y(min)
       grid%ud_iz = one / grid%beta_z(min)
    end if
    
    return
  end subroutine generate_grid

  subroutine compute_beta( grid )

    implicit none

    type(VelocityGridType) :: grid
    integer :: i,j,k
    integer :: i_min, i_max, j_min, j_max, k_min, k_max

    ! In each loop, compute_beta() is stepping through every grid point
    !  value and determining its respective beta. Each beta is determined
    !  as half the distance to the nearest point to the left added to
    !  half the distance to the nearest point to the right.
    !
    !        point 1              point 2    point 3
    !  |--------X----------[----------X----]----X-----------------|
    !                      |-----beta------|

    i_min = grid%i_min
    i_max = grid%i_max

    j_min = grid%j_min
    j_max = grid%j_max

    k_min = grid%k_min
    k_max = grid%k_max

    ! Calculate beta values on x-axis
    grid%beta_x(i_min) = grid%x(i_min+1) - grid%x(i_min)
    grid%beta_x(i_max) = grid%x(i_max) - grid%x(i_max-1)

    do i = (i_min + 1), (i_max - 1)
       grid%beta_x(i) = (grid%x(i+1) - grid%x(i))/2.0d0 &
            + (grid%x(i) - grid%x(i-1))/2.0d0
    end do

    ! Calculate beta values on y-axis
    grid%beta_y(j_min) = grid%y(j_min+1) - grid%y(j_min)
    grid%beta_y(j_max) = grid%y(j_max) - grid%y(j_max-1)

    do j = (j_min + 1), (j_max - 1)
       grid%beta_y(j) = (grid%y(j+1) - grid%y(j))/2.0d0 &
            + (grid%y(j)- grid%y(j-1))/2.0d0
    end do

    ! Calculate beta values on z-axis
    grid%beta_z(k_min) = grid%z(k_min+1) - grid%z(k_min)
    grid%beta_z(k_max) = grid%z(k_max) - grid%z(k_max-1)

    do k = (k_min + 1), (k_max - 1)
       grid%beta_z(k) = (grid%z(k+1) - grid%z(k))/2.0d0&
            + (grid%z(k) - grid%z(k-1))/2.0d0
    end do

    return
  end subroutine compute_beta

  subroutine find_min_max_beta3( grid )

    use Constants

    implicit none

    type(VelocityGridType) :: grid

    integer :: i,j,k
    integer :: i_min, i_max, j_min, j_max, k_min, k_max    

    double precision :: beta_x, beta_y, beta_z, beta3
    double precision :: beta3_min, beta3_max
    double precision :: num_points, beta3_avg

    i_min = grid%i_min
    i_max = grid%i_max
    j_min = grid%j_min
    j_max = grid%j_max
    k_min = grid%k_min
    k_max = grid%k_max

    beta_x = grid%beta_x(i_min)
    beta_y = grid%beta_y(j_min)
    beta_z = grid%beta_y(k_min)
    beta3  = beta_x * beta_y * beta_z

    beta3_min = beta3
    beta3_max = beta3
    beta3_avg = zero

    do k = k_min, k_max
       do j = j_min, j_max
          do i = i_min, i_max

             beta_x = grid%beta_x(i)
             beta_y = grid%beta_y(j)
             beta_z = grid%beta_y(k)
             beta3  = beta_x * beta_y * beta_z

             beta3_avg = beta3_avg + beta3

             if( beta3 .lt. beta3_min ) beta3_min = beta3
             if( beta3 .gt. beta3_max ) beta3_max = beta3

          end do
       end do
    end do

    num_points = dble( i_max - i_min ) * dble( j_max - j_min ) * dble( k_max - k_min )
    beta3_avg = beta3_avg / num_points

    grid%beta3_min = beta3_min
    grid%beta3_max = beta3_max
    grid%beta3_avg = beta3_avg ! The current average is the mean, but that is subject to reinterpretation

    return
  end subroutine find_min_max_beta3

  subroutine generate_file_read_grid( array, element, min, max, direction )

    implicit none

    integer, intent(in) :: element, min, max, direction

    double precision, dimension(min:max) :: array

    integer :: i

    select case( direction )
    case( x_dir )
       do i = min, max
          array(i) = inp_vel_grid(element)%x(i+1)
       end do

    case( y_dir )
       do i = min, max
          array(i) = inp_vel_grid(element)%y(i+1)
       end do

    case( z_dir )
       do i = min, max
          array(i) = inp_vel_grid(element)%z(i+1)
       end do

    case default
       write(*,*) "Error: Invalid value for direction: ", direction
       stop

    end select

    return
  end subroutine generate_file_read_grid

  subroutine generate_uniform_grid( array, min, max, num_points, v_min, v_max )

    implicit none

    integer, intent(in) :: min, max, num_points
    double precision, intent(in) :: v_min, v_max

    double precision, dimension(min:max) :: array

    integer :: i

    do i = min, max
       array(i) = v_min + dble( i - min ) * ( v_max - v_min ) / dble( num_points - 1 )
    end do

    return
  end subroutine generate_uniform_grid

  subroutine find_zero_point( zero_point, array, min, max )
    
    use Constants

    implicit none

    integer, intent(in) :: min, max
    double precision, dimension(min:max) :: array

    integer :: zero_point

    integer :: i
    
    if( array(min) .gt. zero )then
       zero_point = -1
       return
    end if

    do i = min+1, max
       if( array(i) .gt. zero )then
          zero_point = i-1
          return
       end if
    end do

    zero_point = max

    return
  end subroutine find_zero_point

  subroutine global2local_map( l, i_min, j_min, k_min, num_points_x, num_points_y, i, j, k )

    implicit none

    integer, intent(in) :: l, i_min, j_min, k_min, num_points_x, num_points_y
    integer, intent(out) :: i, j, k
    
    ! ( i_min, j_min, k_min ) is equivalent to l = 0 in this mapping scheme
    k = k_min + l / ( num_points_x * num_points_y )

    j = j_min + l / num_points_x - ( k - k_min ) * num_points_y

    i = i_min + l - ( j - j_min ) * num_points_x - ( k - k_min ) * num_points_x * num_points_y

    return
  end subroutine global2local_map

  subroutine local2global_map( i, j, k, i_min, j_min, k_min, num_points_x, num_points_y, l )

    implicit none

    integer, intent(in) :: i, j, k, i_min, j_min, k_min, num_points_x, num_points_y
    integer, intent(out) :: l

    ! The minimum value of l is 0 in this mapping scheme
    l = ( i - i_min ) + ( j - j_min ) * num_points_x + &
         ( k - k_min ) * num_points_x * num_points_y
 
    return
  end subroutine local2global_map

end module VelocityGrid
