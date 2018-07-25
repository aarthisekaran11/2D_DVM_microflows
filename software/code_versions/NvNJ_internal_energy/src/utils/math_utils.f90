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
module MathUtilities

  implicit none

  private

  public :: isgn
  public :: dsgn
  public :: arithmetic_mean
  public :: gamma_calculator

contains

  function isgn(x) result(val)

    implicit none
    
    double precision, intent(in) :: x
    integer :: val

    if(x.lt.0.)then
       val = -1
    else
       val =  1
    endif
    
    return
  end function isgn

  function dsgn( x ) result(val)

    implicit none

    double precision, intent(in) :: x
    double precision :: val

    if(x.lt.0.)then
       val = -1.d0
    else
       val =  1.d0
    endif

    return
  end function dsgn

  function arithmetic_mean( x, elem ) result(val)

    implicit none
    
    double precision, dimension(:), intent(in) :: x
    integer, intent(in) :: elem ! Number of elements in x
    double precision :: val

    integer :: i

    val = 0.0d0
    
    do i = 1, elem
       val = val + x(i)
    end do

    val = val/dble(elem)

    return
  end function arithmetic_mean

  recursive function gamma_calculator( val_in ) result( val_out )
    ! gamma_calculator uses Lanczos approximation to calculate the solution to the
    ! gamma function

    use Constants

    implicit none

    double precision, intent(in) :: val_in
    double precision :: val_out

    double precision :: x, y, z
    integer :: i, iter

    ! the values for the following array were taken from the GNU Scientific Library
    double precision, dimension(0:8) :: array

    array(0) =  0.99999999999980993d0
    array(1) =  676.5203681218851d0
    array(2) = -1259.1392167224028d0
    array(3) =  771.32342877765313d0
    array(4) = -176.61502916214059d0
    array(5) =  12.507343278686905d0
    array(6) = -0.13857109526572012d0
    array(7) =  9.9843695780195716d-6
    array(8) =  1.5056327351493116d-7

    x = val_in

    iter = 7

    if ( x .lt. one_half ) then

       val_out = pi/( sin( pi*x ) * gamma_calculator( one - x ) )

    else

       x = x - one
       y = array(0)

       do i = 1, ( iter + 1 )

          y = y + array(i)/( x + dble( i ) )

       end do

       z = x + dble( iter ) + one_half

       val_out = sqrt( 2.0d0*pi )*z**( x + one_half )*exp( -z )*y

    end if

  end function gamma_calculator

end module MathUtilities
