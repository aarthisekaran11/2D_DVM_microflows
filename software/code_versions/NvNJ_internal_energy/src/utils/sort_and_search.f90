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
module SortAndSearch

  implicit none


  private

  public :: binary_search
  public :: binary_search_with_first_guess

contains

  subroutine binary_search( array, num_points, rand_val, index )

    use Constants
    
    implicit none

    integer, intent(in) :: num_points
    double precision, dimension(0:num_points), intent(in) :: array
    double precision, intent(in) :: rand_val

    integer :: index

    integer :: left_bound, right_bound

    left_bound  = 0
    right_bound = num_points

    index = left_bound + ( right_bound - left_bound )/2

    if( rand_val .lt. double_tol )then
       index = 0
       !NOTE: index cannot equal zero in this method! Therefore (l-1) always exists
       return
    end if

    do
       ! Loop termination condition
       if( ( rand_val .gt. array(index) ) .and. ( rand_val .le. array(index+1) ) ) exit

       if( rand_val .gt. array( index ) )then
          left_bound = index
       else
          right_bound = index
       endif

       index = left_bound + ( right_bound - left_bound )/2

    enddo

    return
  end subroutine binary_search

  subroutine binary_search_with_first_guess( index, val, array, min, max )

    use Constants

    implicit none

    integer, intent(in) :: min, max
    double precision, dimension(min:max), intent(in) :: array
    double precision, intent(in) :: val
    
    integer :: index

    integer :: left_bound, right_bound

    ! TODO: lots of if statements, maybe reduce the number of checks somehow?
    
    left_bound  = min
    right_bound = max

    if( abs( val - array(index) ) .lt. double_tol ) return

    if( index .eq. right_bound ) index = index - 1 

    do

       ! Loop termination condition
       ! If val is between index and index+1
       if( ( val .gt. array(index) ) .and. ( val .le. array(index+1) ) )then
          ! Calculate whether index or index+1 is closer to val
          if( abs(val - array(index)) .gt. abs(val - array(index+1)) ) index = index + 1
          exit
       end if

       if( val .gt. array( index ) )then
          left_bound = index
       else
          right_bound = index
       endif

       index = left_bound + ( right_bound - left_bound )/2

    enddo
    
    return
  end subroutine binary_search_with_first_guess

end module SortAndSearch
