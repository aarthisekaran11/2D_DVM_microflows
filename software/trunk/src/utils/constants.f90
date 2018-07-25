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
module Constants

  implicit none

  ! Constants
  double precision, parameter :: pi        = 3.1459265358979323d0 ! pi
  double precision, parameter :: kb        = 1.3805d-23           ! Boltzmann's constant [m^2-kg/s^2-K]
  double precision, parameter :: Navogadro = 6.023d23             ! Avogadro's constant [1/mol]
  double precision, parameter :: R_univ    = 8.3144621d0          ! Universal gas constant [J/K-mol]
  double precision, parameter :: kg2amu    = 1.66053892d-27       ! Kg to atomic mass units
  double precision, parameter :: Pa2atm    = 101325.0d0           ! Pa to atmospheres
  double precision, parameter :: angstrom  = 1.0d-10              ! meters in an angstrom
  double precision, parameter :: cm2J      = 5.03445d22           ! cm^-1 to Joules
  double precision, parameter :: planck    = 6.62606957d-34       ! Planck's constant [m^2-kg/s]

  ! Double precision numbers
  double precision, parameter :: zero        = 0.0d0
  double precision, parameter :: one         = 1.0d0
  double precision, parameter :: two         = 2.0d0
  double precision, parameter :: four        = 4.0d0
  double precision, parameter :: ten         = 10.0d0
  double precision, parameter :: one_sixth   = 1.0d0/6.0d0
  double precision, parameter :: one_half    = 0.5d0
  double precision, parameter :: one_third   = 1.0d0/3.0d0
  double precision, parameter :: one_fourth  = 0.25d0
  double precision, parameter :: two_thirds  = 2.0d0/3.0d0
  double precision, parameter :: four_thirds = 4.0d0/3.0d0
  double precision, parameter :: double_tol  = 1.0d-13
  double precision, parameter :: cnst99      = 0.99999999999999d0
  double precision, parameter :: bigdp       = 1.0d100

  ! Temperature used to calculate diameters
  double precision, parameter :: temp_diam = 273.0d0

end module Constants
