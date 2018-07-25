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
module ScreenOutput

  implicit none
  
  private

  integer :: verbosity    = 0
  integer :: display_freq = 1

  integer, dimension(2) :: screen_out_loc

  integer, parameter :: pretty_messages_verbosity = 1
  integer, parameter :: properties_verbosity      = 2

  integer, parameter :: header_count = 10

  public :: set_verbosity
  public :: set_display_printing_freq
  public :: set_screen_output_location
  public :: title_message
  public :: init_message
  public :: solve_message
  public :: finalize_message
  public :: print_property_update

contains

  subroutine set_verbosity( verbosity_in )

    implicit none

    integer, intent(in) :: verbosity_in

    verbosity = verbosity_in

    return
  end subroutine set_verbosity

  subroutine set_display_printing_freq( display_freq_in )

    implicit none

    integer, intent(in) :: display_freq_in

    display_freq = display_freq_in

    return
  end subroutine set_display_printing_freq

  subroutine set_screen_output_location( screen_out_loc_in )
    
    use PhysicalGrid

    implicit none

    integer, dimension(2), intent(in) :: screen_out_loc_in

    integer :: nx_space, ny_space

    call get_nspace( nx_space, ny_space )

    screen_out_loc = screen_out_loc_in

    if( screen_out_loc(1) .gt. nx_space ) screen_out_loc(1) = nx_space
    if( screen_out_loc(2) .gt. ny_space ) screen_out_loc(2) = ny_space

    return
  end subroutine set_screen_output_location

  subroutine title_message( )

    implicit none

    if( verbosity .lt. pretty_messages_verbosity ) return

    write(*,"(a)")"==========================================================================================="
    write(*,"(a)")"               DVM: Discrete Velocity Method for the Boltzmann Equation"
    write(*,"(a)")"==========================================================================================="

    return
  end subroutine title_message

  subroutine init_message()

    implicit none

    if( verbosity .lt. pretty_messages_verbosity ) return

    write(*,"(a)")"==========================================================================================="
    write(*,"(a)")"                                  Initializing Solver"
    write(*,"(a)")"==========================================================================================="

    return
  end subroutine init_message

  subroutine solve_message()

    implicit none

    if( verbosity .lt. pretty_messages_verbosity ) return

    write(*,"(a)")"==========================================================================================="
    write(*,"(a)")"                                    Beginning Solve"
    write(*,"(a)")"==========================================================================================="

    return
  end subroutine solve_message

  subroutine finalize_message()

    implicit none

    if( verbosity .lt. pretty_messages_verbosity ) return

    write(*,"(a)")"==========================================================================================="
    write(*,"(a)")"                                   Finishing Solve"
    write(*,"(a)")"==========================================================================================="

    return
  end subroutine finalize_message

  subroutine print_property_update( ntime, num_time_steps, properties, colls_array )

    use PhysicalProperties
    use DistFunc

    implicit none

    type(PropertiesType), dimension(:,:), intent(in) :: properties
    integer, dimension(:,:,:,:) :: colls_array
    integer, intent(in) :: ntime, num_time_steps

    double precision :: dens, u, tr_temp, rot_temp, vib_temp, energy

    integer :: n, x_loc, y_loc

    x_loc = screen_out_loc(1)
    y_loc = screen_out_loc(2)

    if( verbosity .lt. properties_verbosity ) return

    if( mod(ntime,display_freq) .ne. 0 ) return

    if( ntime .eq. 0 .or. mod(ntime,header_count*display_freq) .eq. 0 )then
       write(*,advance="yes",fmt="(8x,8a20)") &
            "density  ", "u-velocity ", "trans. temp", &
            "rot. temp ", "vib. temp ", "energy ", "time step ", "collisions"
    end if

    dens     = properties(x_loc,y_loc)%mix_dens
    u        = properties(x_loc,y_loc)%mix_x_vel
    tr_temp  = properties(x_loc,y_loc)%mix_tr_temp
    rot_temp = properties(x_loc,y_loc)%mix_rot_temp
    vib_temp = properties(x_loc,y_loc)%mix_vib_temp
    energy   = properties(x_loc,y_loc)%mix_energy

    write(*,advance="no",fmt="(a11,6e20.10)") "total:", dens, u, tr_temp, rot_temp, vib_temp, energy
    write(*,advance="no",fmt="(i8,a6,i7)") ntime, " of ", num_time_steps
    write(*,advance="yes",fmt="(i14)") sum(colls_array(:,:,:,:))

    if( num_species .gt. 1 )then
       do n = 1, num_species

          dens     = properties(x_loc,y_loc)%dens(n)
          u        = properties(x_loc,y_loc)%x_vel(n)
          tr_temp  = properties(x_loc,y_loc)%tr_temp(n)
          rot_temp = properties(x_loc,y_loc)%rot_temp(n)
          vib_temp = properties(x_loc,y_loc)%vib_temp(n)
          energy   = properties(x_loc,y_loc)%total_energy(n)

          write(*,advance="no",fmt="(a8,i2,a1,6e20.10)") &
               "species ", n, ":", dens, u, tr_temp, rot_temp, vib_temp, energy

          write(*,advance="no",fmt="(a8,6x,a7)")"--N/A--","--N/A--"
          write(*,advance="yes",fmt="(i14)") colls_array(x_loc,y_loc,n,n)

       end do
    end if

    write(*,advance="yes",fmt="(1x)")
    
    return
  end subroutine print_property_update

end module ScreenOutput
