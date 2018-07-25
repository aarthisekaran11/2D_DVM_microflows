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
module Restart

  implicit none

  private

  character(len=128) :: read_restart_filename, write_restart_filename
  integer, parameter :: read_restart_file_unit  = 50
  integer, parameter :: write_restart_file_unit = 51

  logical :: read_restart_flag
  logical :: write_restart_flag

  integer :: write_restart_freq
  
  public :: set_restart_filename
  public :: set_restart_flags
  public :: get_restart_flag
  public :: write_restart
  public :: read_restart

  ! IMPORTANT:
  ! Restart data is written in unformatted form to reduce file size and speed up file write.
  ! This means only the computer it is written on can read it.

contains
  
  subroutine set_restart_filename( write_restart_filename_in, read_restart_filename_in )

    implicit none

    character(len=128), intent(in) :: write_restart_filename_in, read_restart_filename_in

    write_restart_filename = write_restart_filename_in
    read_restart_filename = read_restart_filename_in

    return
  end subroutine set_restart_filename

  subroutine set_restart_flags( write_restart_flag_in, write_restart_freq_in, read_restart_flag_in )

    use Conversion

    implicit none

    integer, intent(in) :: write_restart_flag_in, write_restart_freq_in, read_restart_flag_in

    write_restart_flag = int2logical( write_restart_flag_in )
    write_restart_freq = write_restart_freq_in

    read_restart_flag  = int2logical( read_restart_flag_in )

    return
  end subroutine set_restart_flags

  subroutine get_restart_flag( read_restart_flag_out )
    
    implicit none

    logical :: read_restart_flag_out

    read_restart_flag_out = read_restart_flag

    return
  end subroutine get_restart_flag

  subroutine write_restart( ntime, num_time_steps, phi, molecule )

    use DistFunc
    use PhysicalGrid
    use SpeciesAndReferenceData
    
    implicit none
    
    type(DistFuncType), dimension(:,:), intent(in) :: phi
    type(MoleculeType), dimension(:), intent(in) :: molecule
    integer, intent(in) :: ntime, num_time_steps

    integer :: r_modes, v_modes
    integer :: ns, species, nspace

    if( write_restart_flag .eqv. .false. ) return
    if( ntime .ne. num_time_steps )then
       if( write_restart_freq .eq. 0 ) return
       if( mod( ntime, write_restart_freq ) .ne. 0 ) return
    end if

    call get_nspace( nspace )

    open( unit=write_restart_file_unit, file=write_restart_filename, status='replace', form='unformatted' )

    write( write_restart_file_unit ) ntime

    do species = 1, num_species

       r_modes = molecule(species)%rot_modes
       v_modes = molecule(species)%vib_modes

       do ns = 1, nspace

          write( write_restart_file_unit ) phi(species,ns)%value

          if( r_modes .gt. 0 )then
             write( write_restart_file_unit ) phi(species,ns)%rot
             write( write_restart_file_unit ) phi(species,ns)%rot_level
          end if

          if( v_modes .gt. 0 )then
             write( write_restart_file_unit ) phi(species,ns)%vib
             write( write_restart_file_unit ) phi(species,ns)%vib_level
          end if

       end do
    end do

    close( write_restart_file_unit )

    return
  end subroutine write_restart

  subroutine read_restart( ntime, phi, molecule )
    
    use DistFunc
    use PhysicalGrid
    use SpeciesAndReferenceData

    implicit none

    type(MoleculeType), dimension(:), intent(in) :: molecule

    type(DistFuncType), dimension(:,:) :: phi
    integer :: ntime

    integer :: r_modes, v_modes
    integer :: ns, species, nspace
    
    if( read_restart_flag .eqv. .false. ) return

    call get_nspace( nspace )

    open( unit=read_restart_file_unit, file=read_restart_filename, status='old', form='unformatted' )

    read( read_restart_file_unit ) ntime

    do species = 1, num_species

       r_modes = molecule(species)%rot_modes
       v_modes = molecule(species)%vib_modes

       do ns = 1, nspace

          read( read_restart_file_unit ) phi(species,ns)%value

          if( r_modes .gt. 0 )then
             read( read_restart_file_unit ) phi(species,ns)%rot
             read( read_restart_file_unit ) phi(species,ns)%rot_level
          end if

          if( v_modes .gt. 0 )then
             read( read_restart_file_unit ) phi(species,ns)%vib
             read( read_restart_file_unit ) phi(species,ns)%vib_level
          end if

       end do
    end do

    close( read_restart_file_unit )

    return
  end subroutine read_restart

end module Restart
