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
     double precision, dimension(7) :: mass
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

    type(VelocityGridType), intent(in) :: vel_grid
    double precision, intent(in) :: u, v, w, f_off

    integer :: o_i, o_j, o_k, sign_a, sign_b, sign_c
    integer :: ext_map_flag
    double precision :: f_o, f_ix, f_iy, f_iz, f_ex, f_ey, f_ez

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: num_points_x, num_points_y, num_points_z
    double precision :: xv_min, xv_max, yv_min, yv_max, zv_min, zv_max

    integer :: sign_u, sign_v, sign_w
    integer :: x_ext_map_flag, y_ext_map_flag, z_ext_map_flag
    double precision :: a, b, c
    double precision :: o_x_vel, o_y_vel, o_z_vel
    double precision :: int_x_vel, int_y_vel, int_z_vel
    double precision :: ext_x_vel, ext_y_vel, ext_z_vel

    xv_min = vel_grid%xv_min
    xv_max = vel_grid%xv_max
    yv_min = vel_grid%yv_min
    yv_max = vel_grid%yv_max
    zv_min = vel_grid%zv_min
    zv_max = vel_grid%zv_max

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    num_points_x = vel_grid%num_points_x
    num_points_y = vel_grid%num_points_y
    num_points_z = vel_grid%num_points_z

    ! Find nearest x-velocity point
    if( u .gt. xv_max )then
       o_i = i_max

    else if( u .lt. xv_min )then
       o_i = i_min

    else
       o_i = nint( dble( num_points_x - 1 )*( u - xv_min )/( xv_max - xv_min ) )
       call binary_search_with_first_guess( o_i, u, vel_grid%x, i_min, i_max )

    end if

    ! Find nearest y-velocity point
    if( v .gt. yv_max )then
       o_j = j_max

    else if( v .lt. yv_min )then
       o_j = j_min

    else
       o_j = nint( dble( num_points_y - 1 )*( v - yv_min )/( yv_max - yv_min ) )
       call binary_search_with_first_guess( o_j, v, vel_grid%y, j_min, j_max )

    end if

    if( w .gt. zv_max )then
       o_k = k_max

    else if( w .lt. zv_min )then
       o_k = k_min

    else
       o_k = nint( dble( num_points_z - 1 )*( w - zv_min )/( zv_max - zv_min ) )
       call binary_search_with_first_guess( o_k, w, vel_grid%z, k_min, k_max )

    end if
    
    sign_u = isgn( u )
    sign_v = isgn( v ) 
    sign_w = isgn( w )

    o_x_vel = vel_grid%x( o_i )
    o_y_vel = vel_grid%y( o_j )
    o_z_vel = vel_grid%z( o_k )

    ! o_i, o_j, o_k should be the coordinates for an origin on the grid.
    ! However, a corner of the cube can't be an origin, because there are no external points available.
    ! The following should correct the origin if it is on a corner.
    if( (o_i .eq. i_max) .or. (o_i .eq. i_min) )then
       if( (o_j .eq. j_max) .or. (o_j .eq. j_min) )then
          if( (o_k .eq. k_max) .or. (o_k .eq. k_min) )then

             ! The origin is a corner, so make the origin the next closest point
             if( (abs(u- o_x_vel) .le. abs(v - o_y_vel)) .and. &
                  (abs(u - o_x_vel) .le. abs(w - o_z_vel))        )then
                o_i = o_i - sign_u

             else if( (abs(v - o_y_vel) .le. abs(u - o_x_vel)) .and. &
                  (abs(v - o_y_vel) .le. abs(w - o_z_vel))        )then
                o_j = o_j - sign_v

             else
                o_k = o_k - sign_w

             endif

          endif
       endif
    endif

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
    if( ( o_i + sign_a .gt. i_max ) .or. ( o_i + sign_a .lt. i_min ) ) sign_a = -1*sign_a

    if( ( o_j + sign_b .gt. j_max ) .or. ( o_j + sign_b .lt. j_min ) ) sign_b = -1*sign_b

    if( ( o_k + sign_c .gt. k_max ) .or. ( o_k + sign_c .lt. k_min ) ) sign_c = -1*sign_c

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
       f_ex = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel )/&
            ( ext_x_vel*( ext_x_vel - int_x_vel ) + &
            ext_y_vel*( ext_y_vel - int_y_vel ) + &
            ext_z_vel*( ext_z_vel - int_z_vel )    )
       f_ey = f_ex
       f_ez = f_ex

       f_ix = ( f_off*a - f_ex*ext_x_vel )/int_x_vel
       f_iy = ( f_off*b - f_ey*ext_y_vel )/int_y_vel
       f_iz = ( f_off*c - f_ez*ext_z_vel )/int_z_vel

       f_o = f_off - f_ix - f_iy - f_iz - f_ex - f_ey - f_ez

    case( x_only )
       f_ex = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel )/&
            ( ext_x_vel*( ext_x_vel - int_x_vel ) )
       f_ey = zero
       f_ez = zero

       f_ix = ( f_off*a - f_ex*ext_x_vel )/int_x_vel
       f_iy = ( f_off*b - f_ey*ext_y_vel )/int_y_vel
       f_iz = ( f_off*c - f_ez*ext_z_vel )/int_z_vel

       f_o = f_off - f_ix - f_iy - f_iz - f_ex - f_ey - f_ez


    case( y_only ) 
       f_ey = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel )/&
            ( ext_y_vel*( ext_y_vel - int_y_vel ) )
       f_ex = zero
       f_ez = zero

       f_ix = ( f_off*a - f_ex*ext_x_vel )/int_x_vel
       f_iy = ( f_off*b - f_ey*ext_y_vel )/int_y_vel
       f_iz = ( f_off*c - f_ez*ext_z_vel )/int_z_vel

       f_o = f_off - f_ix - f_iy - f_iz - f_ex - f_ey - f_ez

    case( z_only )
       f_ez = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel )/&
            ( ext_z_vel*( ext_z_vel - int_z_vel ) )
       f_ex = zero
       f_ey = zero

       f_ix = ( f_off*a - f_ex*ext_x_vel )/int_x_vel
       f_iy = ( f_off*b - f_ey*ext_y_vel )/int_y_vel
       f_iz = ( f_off*c - f_ez*ext_z_vel )/int_z_vel

       f_o = f_off - f_ix - f_iy - f_iz - f_ex - f_ey - f_ez

    case( xy_only )
       f_ex = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel )/&
            ( ext_x_vel*( ext_x_vel - int_x_vel ) + ext_y_vel*( ext_y_vel - int_y_vel ) )
       f_ey = f_ex
       f_ez = zero

       f_ix = ( f_off*a - f_ex*ext_x_vel )/int_x_vel
       f_iy = ( f_off*b - f_ey*ext_y_vel )/int_y_vel
       f_iz = ( f_off*c - f_ez*ext_z_vel )/int_z_vel

       f_o = f_off - f_ix - f_iy - f_iz - f_ex - f_ey - f_ez

    case( xz_only )
       f_ex = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel )/&
            ( ext_x_vel*( ext_x_vel - int_x_vel ) + ext_z_vel*( ext_z_vel - int_z_vel ) )
       f_ey = zero
       f_ez = f_ex

       f_ix = ( f_off*a - f_ex*ext_x_vel )/int_x_vel
       f_iy = ( f_off*b - f_ey*ext_y_vel )/int_y_vel
       f_iz = ( f_off*c - f_ez*ext_z_vel )/int_z_vel

       f_o = f_off - f_ix - f_iy - f_iz - f_ex - f_ey - f_ez

    case( yz_only )
       f_ex = zero
       f_ey = f_off*( a*a + b*b + c*c - a*int_x_vel - b*int_y_vel - c*int_z_vel )/&
            ( ext_y_vel*( ext_y_vel - int_y_vel ) + ext_z_vel*( ext_z_vel - int_z_vel ) )
       f_ez = f_ey

       f_ix = ( f_off*a - f_ex*ext_x_vel )/int_x_vel
       f_iy = ( f_off*b - f_ey*ext_y_vel )/int_y_vel
       f_iz = ( f_off*c - f_ez*ext_z_vel )/int_z_vel

       f_o = f_off - f_ix - f_iy - f_iz - f_ex - f_ey - f_ez

    case default
       write(*,*) "Error: Invalid value of ext_map_flag: ", ext_map_flag
       stop

    end select

    return
  end subroutine compute_remapping_quantities

end module Remapping
