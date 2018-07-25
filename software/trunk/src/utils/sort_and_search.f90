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

  public :: referenced_search
  public :: binary_search
  public :: binary_search_with_first_guess
  public :: ssort
  public :: ssort1

contains

  subroutine referenced_search( array, ref_array, num_points, rand_val, index )

    use Constants

    implicit none

    integer, intent(in) :: num_points
    double precision, dimension(0:num_points), intent(in) :: array
    integer, dimension(1:1000,2), intent(in) :: ref_array
    double precision, intent(in) :: rand_val

    integer :: i, index
    integer :: left_bound, right_bound

    double precision :: adj_rand

    if( rand_val .lt. double_tol )then
       index = 0
       return
    end if

    i = ceiling( rand_val * 1000.0d0 )
    
    left_bound  = ref_array( i, 1 )
    right_bound = ref_array( i, 2 )

    index = left_bound + ( right_bound - left_bound )/2

    if( left_bound .eq. right_bound ) return

    adj_rand = rand_val * array( num_points )

    do
       if( ( adj_rand .gt. array(index) ) .and. ( adj_rand .le. array(index+1) ) )  return

       if( adj_rand .gt. array(index) )then
          left_bound = index
       else
          right_bound = index
       end if

       index = left_bound + ( right_bound - left_bound )/2

    end do

    return
  end subroutine referenced_search

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

    if( array(min) .gt. val )then
       index = min
       return
    end if
    if( array(max) .lt. val )then
       index = max
       return
    end if
    
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

  subroutine ssort2( x, iy, n )
    implicit none

    integer :: n
    double precision :: x(*)
    integer :: iy(*)
    double precision :: temp
    integer :: i, j, jmax, itemp

    jmax=n-1
    do i=1, n-1
       temp=1.d38
       do j=1,jmax
          if(x(j) .gt. x(j+1)) cycle
          temp=x(j)
          x(j)=x(j+1)
          x(j+1)=temp
          itemp=iy(j)
          iy(j)=iy(j+1)
          iy(j+1)=itemp
       end do
       if(temp .eq. 1.d38) exit
       jmax=jmax-1
    end do
    return
  end subroutine ssort2

  subroutine ssort( s, x, y, z, n )

    implicit none

    integer, intent(in) :: n

    double precision, dimension(:) :: s
    integer, dimension(:) :: x, y, z

    double precision :: temp
    integer :: i, iswap(1), itemp, iswap1

    intrinsic maxloc

    do i = 1, n-1

       iswap = minloc( s(i:n) )
       iswap1 = iswap(1) + i - 1

       if( iswap1 .ne. i )then
          temp = s(i)
          s(i) = s(iswap1)
          s(iswap1) = temp

          itemp = x(i)
          x(i) = x(iswap1)
          x(iswap1) = itemp

          itemp = y(i)
          y(i) = y(iswap1)
          y(iswap1) = itemp

          itemp = z(i)
          z(i) = z(iswap1)
          z(iswap1) = itemp
       end if

    end do

    return
  end subroutine ssort

  subroutine ssort1( s, n )

    implicit none

    integer, intent(in) :: n

    double precision, dimension(:) :: s

    double precision :: temp
    integer :: i, iswap(1), iswap1

    intrinsic maxloc

    do i = 1, n-1

       iswap = minloc( s(i:n) )
       iswap1 = iswap(1) + i - 1

       if( iswap1 .ne. i )then
          temp = s(i)
          s(i) = s(iswap1)
          s(iswap1) = temp
       end if

    end do

    return
  end subroutine ssort1

end module SortAndSearch
