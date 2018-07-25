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

  integer, parameter :: rankine_hugoniot = 4
  
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
    integer :: nx, species, nx_space, ny_space

    if( write_restart_flag .eqv. .false. ) return
    if( ntime .ne. num_time_steps )then
       if( write_restart_freq .eq. 0 ) return
       if( mod( ntime, write_restart_freq ) .ne. 0 ) return
    end if

    call get_nspace( nx_space, ny_space )

    open( unit=write_restart_file_unit, file=write_restart_filename, status='replace', form='unformatted' )

    write( write_restart_file_unit ) ntime

    do species = 1, num_species

       r_modes = molecule(species)%rot_modes
       v_modes = molecule(species)%vib_modes

       do nx = 1, nx_space

          write( write_restart_file_unit ) phi(species,nx)%value

          if( r_modes .gt. 0 )then
             write( write_restart_file_unit ) phi(species,nx)%rot
             write( write_restart_file_unit ) phi(species,nx)%rot_level
          end if

          if( v_modes .gt. 0 )then
             write( write_restart_file_unit ) phi(species,nx)%vib
             write( write_restart_file_unit ) phi(species,nx)%vib_level
          end if

       end do
    end do

    close( write_restart_file_unit )

    return
  end subroutine write_restart

  subroutine read_restart( ntime, phi, delta_phi, shock_props, molecule )
    
    use DistFunc
    use PhysicalGrid
    use SpeciesAndReferenceData
    use InitialConditions
    use BoundaryConditions
    use ShockConditions

    implicit none

    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(NormalShock), intent(in) :: shock_props

    type(DistFuncType), dimension(:,:) :: phi, delta_phi
    integer :: ntime

    integer :: r_modes, v_modes
    integer :: nx, species, nx_space, ny_space
    
    double precision :: dens_LW, u_LW, v_LW, w_LW, temp_LW, temp_rot_LW, temp_vib_LW
    double precision :: dens_RW, u_RW, v_RW, w_RW, temp_RW, temp_rot_RW, temp_vib_RW
    double precision :: dens_BW, u_BW, v_BW, w_BW, temp_BW, temp_rot_BW, temp_vib_BW
    double precision :: dens_TW, u_TW, v_TW, w_TW, temp_TW, temp_rot_TW, temp_vib_TW

    integer :: init_domain

    if( read_restart_flag .eqv. .false. ) return

    call get_nspace( nx_space, ny_space )

    open( unit=read_restart_file_unit, file=read_restart_filename, status='old', form='unformatted' )

    read( read_restart_file_unit ) ntime

    do species = 1, num_species

       r_modes = molecule(species)%rot_modes
       v_modes = molecule(species)%vib_modes

       ! Get boundary condition properties we need
       call get_left_wall_density( dens_LW, species )
       call get_right_wall_density( dens_RW, species )
       call get_left_wall_velocity( u_LW, v_LW, w_LW, species )
       call get_right_wall_velocity( u_RW, v_RW, w_RW, species )
       call get_left_wall_temp( temp_LW, species )
       call get_right_wall_temp( temp_RW, species )

       call get_spatial_domain_conditions( init_domain, species )
       ! Set Rankine-Hugoniot conditions at left wall if rankine hugoniot
       ! TODO: hard coded for left wall to be downstream
       if( init_domain .eq. rankine_hugoniot )then

          dens_LW    = shock_props%dens_ratio*dens_RW
          u_LW       = shock_props%u_down
          temp_LW    = shock_props%T_ratio*temp_RW

          ! Update boundary properties to shock properties
          call set_boundary_properties( dens_LW, u_LW, v_LW, w_LW, temp_LW, &
               dens_RW, u_RW, v_RW, w_RW, temp_RW, &
               dens_BW, u_BW, v_BW, w_BW, temp_BW, &
               dens_TW, u_TW, v_TW, w_TW, temp_TW,species )

          ! Internal energy boundaries
          if( r_modes .gt. 0 )then
             temp_rot_LW = temp_LW
             call set_temp_rot( temp_rot_LW, temp_rot_RW, &
                  temp_rot_BW, temp_rot_TW, species )
          end if

          if( v_modes .gt. 0 )then
             temp_vib_LW = temp_LW
             call set_temp_vib( temp_vib_LW, temp_vib_RW, &
                  temp_vib_LW, temp_vib_RW, species )
          end if

       end if

       do nx = 1, nx_space

          read( read_restart_file_unit ) phi(species,nx)%value
          delta_phi(species,nx)%value = phi(species,nx)%value

          if( r_modes .gt. 0 )then
             read( read_restart_file_unit ) phi(species,nx)%rot
             read( read_restart_file_unit ) phi(species,nx)%rot_level
             delta_phi(species,nx)%rot       = phi(species,nx)%rot
             delta_phi(species,nx)%rot_level = phi(species,nx)%rot_level
          end if

          if( v_modes .gt. 0 )then
             read( read_restart_file_unit ) phi(species,nx)%vib
             read( read_restart_file_unit ) phi(species,nx)%vib_level
             delta_phi(species,nx)%vib       = phi(species,nx)%vib
             delta_phi(species,nx)%vib_level = phi(species,nx)%vib_level
          end if

       end do
    end do

    close( read_restart_file_unit )

    return
  end subroutine read_restart

end module Restart
