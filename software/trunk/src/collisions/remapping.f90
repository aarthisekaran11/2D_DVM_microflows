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
module Remapping

  use MathUtilities
  use Constants
  use SortAndSearch

  implicit none

  private

  type MappingResultType
     integer, dimension(7) :: stencil
     double precision, dimension(8) :: mass
  end type MappingResultType

  integer, parameter :: x_only  = 1
  integer, parameter :: y_only  = 2
  integer, parameter :: z_only  = 4
  integer, parameter :: xy_only = 3
  integer, parameter :: xz_only = 5
  integer, parameter :: yz_only = 6
  integer, parameter :: xyz     = 7

  public :: MappingResultType

  public :: perform_remapping
  public :: compute_remapping_quantities

contains

  subroutine perform_remapping( mapping, u, v, w, interpolant, vel_grid )

    use VelocityGrid

    implicit none

    type(VelocityGridType), intent(in) :: vel_grid
    double precision, intent(in) :: u, v, w, interpolant
    
    type(MappingResultType) :: mapping

    integer :: o_i, o_j, o_k, a, b, c
    integer :: ext_map_flag
    double precision :: f_o, f_ix, f_iy, f_iz, f_ex, f_ey, f_ez

    call compute_remapping_quantities( u, v, w, interpolant, vel_grid, o_i, o_j, o_k, a, b, c, &
         f_o, f_ix, f_iy, f_iz, f_ex, f_ey, f_ez, ext_map_flag )

    mapping%stencil(1) = o_i
    mapping%stencil(2) = o_j
    mapping%stencil(3) = o_k
    mapping%stencil(4) = a
    mapping%stencil(5) = b
    mapping%stencil(6) = c
    mapping%stencil(7) = ext_map_flag

    mapping%mass(1) = f_o
    mapping%mass(2) = f_ix
    mapping%mass(3) = f_iy
    mapping%mass(4) = f_iz
    mapping%mass(5) = f_ex
    mapping%mass(6) = f_ey
    mapping%mass(7) = f_ez

    return
  end subroutine perform_remapping
    
  subroutine compute_remapping_quantities( u, v, w, f_off, vel_grid, o_i, o_j, o_k, sign_a, sign_b, sign_c, &
       f_o, f_ix, f_iy, f_iz, f_ex, f_ey, f_ez, ext_map_flag )

    use VelocityGrid

    implicit none

    type(VelocityGridType), intent(in) :: vel_grid
    double precision, intent(in) :: u, v, w, f_off

    integer :: o_i, o_j, o_k, sign_a, sign_b, sign_c
    integer :: ext_map_flag
    double precision :: f_o, f_ix, f_iy, f_iz, f_ex, f_ey, f_ez, f_e

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: num_points_x, num_points_y, num_points_z
    double precision :: xv_min, yv_min, zv_min
    double precision :: xv_factor, yv_factor, zv_factor
    integer :: corner

    integer :: sign_u, sign_v, sign_w
    integer :: x_ext_map_flag, y_ext_map_flag, z_ext_map_flag
    double precision :: a, b, c
    double precision :: o_x_vel, o_y_vel, o_z_vel
    double precision :: int_x_vel, int_y_vel, int_z_vel
    double precision :: ext_x_vel, ext_y_vel, ext_z_vel

    logical, dimension(3) :: uniform_grid

    uniform_grid = vel_grid%uniform_grid
    uniform_grid = .false.

    xv_min = vel_grid%xv_min
    yv_min = vel_grid%yv_min
    zv_min = vel_grid%zv_min

    xv_factor = vel_grid%xv_factor
    yv_factor = vel_grid%yv_factor
    zv_factor = vel_grid%zv_factor

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    num_points_x = vel_grid%num_points_x
    num_points_y = vel_grid%num_points_y
    num_points_z = vel_grid%num_points_z

    corner = 0

    ! TODO: speed up here by reorganizing if statement or figuring out how to do it differently
    ! Find nearest x-velocity point
    o_i = nint( xv_factor * ( u - xv_min ) )
    if( uniform_grid(1) .eqv. .false. )&
         call binary_search_with_first_guess( o_i, u, vel_grid%x, i_min, i_max )

    if( o_i .ge. i_max )then
       o_i = i_max
       corner = corner + 1
    else if( o_i .le. i_min )then
       o_i = i_min
       corner = corner + 1
    end if

    ! Find nearest y-velocity point
    o_j = nint( yv_factor * ( v - yv_min ) )
    if( uniform_grid(2) .eqv. .false. )&
         call binary_search_with_first_guess( o_j, v, vel_grid%y, j_min, j_max )

    if( o_j .ge. j_max )then
       o_j = j_max
       corner = corner + 1
    else if( o_j .le. j_min )then
       o_j = j_min
       corner = corner + 1
    end if

    ! Find nearest z-velocity point
    o_k = nint( zv_factor * ( w - zv_min ) )
    if( uniform_grid(3) .eqv. .false. )&
         call binary_search_with_first_guess( o_k, w, vel_grid%z, k_min, k_max )

    if( o_k .ge. k_max )then
       o_k = k_max
       corner = corner + 1
    else if( o_k .le. k_min )then
       o_k = k_min
       corner = corner + 1
    end if

    o_x_vel = vel_grid%x( o_i )
    o_y_vel = vel_grid%y( o_j )
    o_z_vel = vel_grid%z( o_k )

    ! o_i, o_j, o_k should be the coordinates for an origin on the grid.
    ! However, a corner of the cube can't be an origin, because there are no external points available.
    ! The following should correct the origin if it is on a corner.
    if( corner .eq. 3 )then
       sign_u = sign( 1, o_i - (i_min + 1) )
       sign_v = sign( 1, o_j - (j_min + 1) ) 
       sign_w = sign( 1, o_k - (k_min + 1) )
       o_i = o_i - sign_u
       o_j = o_j - sign_v
       o_k = o_k - sign_w
    end if

    o_x_vel = vel_grid%x( o_i )
    o_y_vel = vel_grid%y( o_j )
    o_z_vel = vel_grid%z( o_k )

    a = u - o_x_vel
    b = v - o_y_vel
    c = w - o_z_vel

    sign_a = isgn(a)
    sign_b = isgn(b)
    sign_c = isgn(c)

    ! Guarantee interior points are on the grid
    ! Interior points will be located at origin+(sign_a, sign_b, or sign_c)
    if( ( o_i + sign_a .gt. i_max ) .or. ( o_i + sign_a .lt. i_min ) )then
       sign_a = -1*sign_a
       call find_boundary_stencil( sign_a, sign_b, sign_c, a, o_i, o_j, o_k, vel_grid%x, vel_grid%y, vel_grid%z, &
            i_min, i_max, j_min, j_max, k_min, k_max, num_points_x, num_points_y, num_points_z )
    end if

    if( ( o_j + sign_b .gt. j_max ) .or. ( o_j + sign_b .lt. j_min ) )then
       sign_b = -1*sign_b
       call find_boundary_stencil( sign_b, sign_a, sign_c, b, o_j, o_i, o_k, vel_grid%y, vel_grid%x, vel_grid%z, &
            j_min, j_max, i_min, i_max, k_min, k_max, num_points_x, num_points_y, num_points_z )
    end if

    if( ( o_k + sign_c .gt. k_max ) .or. ( o_k + sign_c .lt. k_min ) )then
       sign_c = -1*sign_c
       call find_boundary_stencil( sign_c, sign_b, sign_a, c, o_k, o_j, o_i, vel_grid%z, vel_grid%y, vel_grid%x, &
            k_min, k_max, j_min, j_max, i_min, i_max, num_points_x, num_points_y, num_points_z )
    end if

    ! common quantities to all cases:
    int_x_vel = vel_grid%x( o_i + sign_a ) - o_x_vel
    int_y_vel = vel_grid%y( o_j + sign_b ) - o_y_vel
    int_z_vel = vel_grid%z( o_k + sign_c ) - o_z_vel

    ! Now, find possible exterior points.
    ext_map_flag   = 0
    x_ext_map_flag = 0
    y_ext_map_flag = 0
    z_ext_map_flag = 0

    ext_x_vel = zero
    ext_y_vel = zero
    ext_z_vel = zero

    if( ( o_i - sign_a .le. i_max ) .and. ( o_i - sign_a .ge. i_min ) )then
       x_ext_map_flag = 1
       ext_x_vel = vel_grid%x( o_i - sign_a ) - o_x_vel
    endif

    if( ( o_j - sign_b .le. j_max ) .and. ( o_j - sign_b .ge. j_min ) )then
       y_ext_map_flag = 2
       ext_y_vel = vel_grid%y( o_j - sign_b ) - o_y_vel
    endif

    if( ( o_k - sign_c .le. k_max ) .and. ( o_k - sign_c .ge. k_min ) )then
       z_ext_map_flag = 4
       ext_z_vel = vel_grid%z( o_k - sign_c ) - o_z_vel
    endif

    ext_map_flag = x_ext_map_flag + y_ext_map_flag + z_ext_map_flag

    select case( ext_map_flag )
    case( xyz )
       if( all(uniform_grid) )then
          f_e = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel ) * vel_grid%uniform_divisor
       else
          f_e = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel )/&
               ( ext_x_vel*( ext_x_vel - int_x_vel ) + &
               ext_y_vel*( ext_y_vel - int_y_vel ) + &
               ext_z_vel*( ext_z_vel - int_z_vel )    )
       end if
       f_ex = f_e
       f_ey = f_e
       f_ez = f_e

    case( x_only )
       if( all(uniform_grid) )then
          f_e = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel ) * vel_grid%ud_x
       else
          f_e = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel )/&
               ( ext_x_vel*( ext_x_vel - int_x_vel ) )
       end if
       f_ex = f_e
       f_ey = zero
       f_ez = zero

    case( y_only )
       if( all(uniform_grid) )then
          f_e = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel ) * vel_grid%ud_y
       else
          f_e = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel )/&
               ( ext_y_vel*( ext_y_vel - int_y_vel ) )
       end if
       f_ex = zero
       f_ey = f_e
       f_ez = zero

    case( z_only )
       if( all(uniform_grid) )then
          f_e = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel ) * vel_grid%ud_z
       else
          f_e = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel )/&
               ( ext_z_vel*( ext_z_vel - int_z_vel ) )
       end if
       f_ex = zero
       f_ey = zero
       f_ez = f_e

    case( xy_only )
       if( all(uniform_grid) )then
          f_e = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel ) * vel_grid%ud_xy
       else
          f_e = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel )/&
               ( ext_x_vel*( ext_x_vel - int_x_vel ) + &
               ext_y_vel*( ext_y_vel - int_y_vel ) )
       end if
       f_ex = f_e
       f_ey = f_e
       f_ez = zero

    case( xz_only )
       if( all(uniform_grid) )then
          f_e = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel ) * vel_grid%ud_xz
       else
          f_e = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel )/&
               ( ext_x_vel*( ext_x_vel - int_x_vel ) + &
               ext_z_vel*( ext_z_vel - int_z_vel )    )
       end if
       f_ex = f_e
       f_ey = zero
       f_ez = f_e

    case( yz_only )
       if( all(uniform_grid) )then
          f_e = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel ) * vel_grid%ud_yz
       else
          f_e = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel )/&
               ( ext_y_vel*( ext_y_vel - int_y_vel ) + &
               ext_z_vel*( ext_z_vel - int_z_vel )    )
       end if
       f_ex = zero
       f_ey = f_e
       f_ez = f_e

    case default
       write(*,*) "Error: Invalid value of ext_map_flag: ", ext_map_flag
       write(*,*)a, b, c, sign_a, sign_b, sign_c
       write(*,*)o_i, o_j, o_k
       write(*,*)i_min, i_max, j_min, j_max, k_min, k_max
       write(*,*)corner, sign_u, sign_v, sign_w
       stop

    end select

    if( all(uniform_grid) )then
       f_ix = ( f_off*a - f_ex*ext_x_vel ) * vel_grid%ud_ix * dble(sign_a)
       f_iy = ( f_off*b - f_ey*ext_y_vel ) * vel_grid%ud_iy * dble(sign_b)
       f_iz = ( f_off*c - f_ez*ext_z_vel ) * vel_grid%ud_iz * dble(sign_c)
    else
       f_ix = ( f_off*a - f_ex*ext_x_vel ) / int_x_vel
       f_iy = ( f_off*b - f_ey*ext_y_vel ) / int_y_vel
       f_iz = ( f_off*c - f_ez*ext_z_vel ) / int_z_vel
    end if

    f_o = f_off - f_ix - f_iy - f_iz - f_ex - f_ey - f_ez

!!$    if( ext_map_flag .eq. xy_only )then
!!$       write(*,*)f_off
!!$       write(*,*)f_ex, f_ey, f_ez
!!$       write(*,*)f_ix, f_iy, f_iz
!!$       write(*,*)f_o
!!$       write(*,*)a, b, c, sign_a, sign_b, sign_c
!!$       write(*,*)int_x_vel, int_y_vel, int_z_vel
!!$       write(*,*)ext_x_vel, ext_y_vel, ext_z_vel
!!$       write(*,*)"-------------------------------------------------------------------------"
!!$    end if

!!$    if( f_ix .lt. -one .or. f_o .gt. one )then
!!$       write(*,*)"================================================"
!!$       write(*,*)u, v, w, o_x_vel, o_y_vel, o_z_vel
!!$       write(*,*)o_i, o_j, o_k
!!$       write(*,*)a, b, c, sign_a, sign_b, sign_c
!!$       write(*,*)int_x_vel, int_y_vel, int_z_vel
!!$       write(*,*)ext_x_vel, ext_y_vel, ext_z_vel
!!$       write(*,*)f_o, f_ix, f_iy, f_iz, f_ex, f_ey, f_ez
!!$       write(*,*)ext_map_flag
!!$       write(*,*)"-------------------------------------------------"
!!$    end if

    return
  end subroutine compute_remapping_quantities

  subroutine find_boundary_stencil( sign_a, sign_b, sign_c, abc, o_i, o_j, o_k, x, y, z, &
       i_min, i_max, j_min, j_max, k_min, k_max, num_points_x, num_points_y, num_points_z )
    
    use VelocityGrid

    implicit none

    integer, intent(in) :: o_i, o_j, o_k
    integer, intent(in) :: i_min, i_max, j_min, j_max, k_min, k_max
    integer, intent(in) :: num_points_x, num_points_y, num_points_z
    double precision, intent(in) :: x(0:num_points_x-1), y(0:num_points_y-1), z(0:num_points_z-1)
    double precision, intent(in) :: abc

    integer :: sign_a, sign_b, sign_c
       
    do while( o_i+sign_a .gt. i_min .and. o_i+sign_a .lt. i_max )
       if( abs(abc) .gt. one_half * abs( x(o_i+sign_a) - x(o_i) ) )then
          sign_a = sign_a + sign(1,sign_a)
       else
          exit
       end if
    end do

    do while( o_j+sign_b .gt. j_min .and. o_j+sign_b .lt. j_max .and. &
         o_j-sign_b .gt. j_min .and. o_j-sign_b .lt. j_max )
       if( abs(abc) .gt. one_half * abs( y(o_j+sign_b) -  y(o_j) ) )then
          sign_b = sign_b + sign(1,sign_b)
       else
          exit
       end if
    end do

    do while( o_k+sign_c .gt. k_min .and. o_k+sign_c .lt. k_max .and. &
         o_k-sign_c .gt. k_min .and. o_k-sign_c .lt. k_max )
       if( abs(abc) .gt. one_half * abs( z(o_k+sign_c) - z(o_k) ) )then
          sign_c = sign_c + sign(1,sign_c)
       else
          exit
       end if
    end do

    return
  end subroutine find_boundary_stencil


end module Remapping
