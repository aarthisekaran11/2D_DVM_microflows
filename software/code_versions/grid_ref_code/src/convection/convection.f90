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
module Convection

  use Constants
  use SpeciesAndReferenceData

  implicit none

  private

  logical :: enable_convection = .false.
  integer :: convection_update

  integer, parameter :: first_order_upwind = 1
  integer, parameter :: fourth_order_tan   = 4

  public :: set_convection_update
  public :: get_convection_flag
  public :: compute_convection

contains

  subroutine set_convection_update( convection_update_in )

    implicit none

    integer, intent(in) :: convection_update_in

    convection_update = convection_update_in

    if( convection_update .gt. 0 ) enable_convection = .true.

    return
  end subroutine set_convection_update

  subroutine get_convection_flag( enable_convection_out )

    implicit none

    logical :: enable_convection_out

    enable_convection_out = enable_convection

    return
  end subroutine get_convection_flag

  subroutine compute_convection( phi, phi_conv, molecule, vel_grid, species, deltat )

    use DistFunc
    use VelocityGrid
    use PhysicalGrid
    
    implicit none

    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: deltat
    integer, intent(in) :: species

    type(DistFuncType), dimension(:) :: phi, phi_conv

    integer :: nx_space, ny_space, nx, ny, grid_ref

    if( enable_convection .eqv. .false. ) return

    call get_nspace( nx_space, ny_space )

    select case( convection_update )
    case( first_order_upwind )
       call first_order_scheme( phi, phi_conv, molecule, vel_grid, species, deltat )

    case( fourth_order_tan )
       call fourth_order_scheme( phi, phi_conv, molecule, vel_grid, species, deltat )

    case default
       write(*,*) "Error: Invalid value of convection_update - :", convection_update
       stop
    end select

    do nx = 1, nx_space
       call get_spatial_reference( nx, grid_ref )
       call perform_dist_func_update( phi(nx), phi_conv(nx), molecule, vel_grid(grid_ref) )
    end do

    return
  end subroutine compute_convection

  subroutine first_order_scheme( phi, phi_conv, molecule, vel_grid, species, deltat )

    use DistFunc
    use VelocityGrid
    use PhysicalGrid
    use BoundaryConditions
    
    implicit none

    type(DistFuncType), dimension(:) :: phi_conv

    type(DistFuncType), dimension(:), intent(in) :: phi
    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: deltat
    integer, intent(in) :: species

    double precision :: deltax, deltay
    integer :: nx_space, ny_space

    double precision :: cfl, factor, x_vel
    double precision :: phi_0, phi_1, phi_m1

    integer :: r_modes, v_modes

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: i_zero, i_next
    integer :: r_levels, v_levels
    integer :: l, i, j, k, nx, gr, gr_1, gr_m1

    call get_delta_x( deltax, deltay )
    call get_nspace( nx_space, ny_space )

    factor = deltat/deltax

    r_modes = molecule%rot_modes
    v_modes = molecule%vib_modes

    do nx = 2, nx_space-1

       call get_spatial_reference( nx, gr )
       call get_spatial_reference( nx+1, gr_1 )
       call get_spatial_reference( nx-1, gr_m1 )

       r_levels = phi(nx)%num_rot_levels
       v_levels = phi(nx)%num_vib_levels

       i_min = vel_grid(gr)%i_min
       i_max = vel_grid(gr)%i_max
       j_min = vel_grid(gr)%j_min
       j_max = vel_grid(gr)%j_max
       k_min = vel_grid(gr)%k_min
       k_max = vel_grid(gr)%k_max

       i_zero = vel_grid(gr)%i_zero
       i_next = i_zero + 1

!!$       do i_min, i_max
!!$
!!$          x_vel = vel_grid(gr)%x(i)
!!$          cfl   = factor * x_vel
!!$
!!$          if( i .le. i_zero )then
!!$
!!$             phi_conv(nx)%value(i,:,:) = ( one + cfl ) * phi(nx)%value(i,:,:) - cfl * phi(nx+1)%value(i,:,:)
!!$
!!$             if( r_modes .gt. 0 )then
!!$                phi_conv(nx)%rot(:,i,:,:) = ( one + cfl ) * phi(nx)%rot(:,i,:,:) - cfl * phi(nx+1)%rot(:,i,:,:)
!!$             end if
!!$
!!$             if( v_modes .gt. 0 )then
!!$                phi_conv(nx)%vib(:,i,:,:) = ( one + cfl ) * phi(nx)%vib(:,i,:,:) - cfl * phi(nx+1)%vib(:,i,:,:)
!!$             end if
!!$
!!$          end if
!!$
!!$          if( i .gt. i_zero )then
!!$
!!$             phi_conv(nx)%value(i,:,:) = ( one - cfl ) * phi(nx)%value(i,:,:) + cfl * phi(nx-1)%value(i,:,:)
!!$
!!$             if( r_modes .gt. 0 )then
!!$                phi_conv(nx)%rot(:,i,:,:) = ( one - cfl ) * phi(nx)%rot(:,i,:,:) + cfl * phi(nx-1)%rot(:,i,:,:) )
!!$             end if
!!$
!!$             if( v_modes .gt. 0 )then
!!$                phi_conv(nx)%vib(:,i,:,:) = ( one - cfl ) * phi(nx)%vib(:,i,:,:) + cfl * phi(nx-1)%vib(:,i,:,:)
!!$             end if
!!$
!!$          end if
!!$
!!$       end do
!!$
!!$
       do k = k_min, k_max
          do j = j_min, j_max
             do i = i_min, i_max

                x_vel = vel_grid(gr)%x(i)
                cfl = factor*x_vel

                phi_0 = phi(nx)%value(i,j,k)

                if( x_vel .lt. zero )then

                   phi_1 = phi(nx+1)%value(i,j,k)

                   phi_conv(nx)%value(i,j,k) = phi_0 - cfl*( phi_1 - phi_0 )

                   if( r_modes .gt. 0 )then
                      phi_conv(nx)%rot(:,i,j,k) = phi(nx)%rot(:,i,j,k) - &
                           cfl*( phi(nx+1)%rot(:,i,j,k) - phi(nx)%rot(:,i,j,k) )
                   end if

                   if( v_modes .gt. 0 )then
                         phi_conv(nx)%vib(:,i,j,k) = phi(nx)%vib(:,i,j,k) - &
                              cfl*( phi(nx+1)%vib(:,i,j,k) - phi(nx)%vib(:,i,j,k) )
                   end if

                else
                   
                   phi_m1 = phi(nx-1)%value(i,j,k)

                   phi_conv(nx)%value(i,j,k) = phi_0 - cfl*( phi_0 - phi_m1 )

                   if( r_modes .gt. 0 )then
                         phi_conv(nx)%rot(:,i,j,k) = phi(nx)%rot(:,i,j,k) - &
                              cfl*( phi(nx)%rot(:,i,j,k) - phi(nx-1)%rot(:,i,j,k) )
                   end if

                   if( v_modes .gt. 0 )then
                         phi_conv(nx)%vib(:,i,j,k) = phi(nx)%vib(:,i,j,k) - &
                              cfl*( phi(nx)%vib(:,i,j,k) - phi(nx-1)%vib(:,i,j,k) )
                   end if

                end if

             end do
          end do
       end do

    end do

    call apply_boundary_conditions( phi, phi_conv, molecule, vel_grid, factor, species )

    return
  end subroutine first_order_scheme

  subroutine fourth_order_scheme( phi, phi_conv, molecule, vel_grid, species, deltat )

    use DistFunc
    use VelocityGrid
    use PhysicalGrid
    use BoundaryConditions

    implicit none

    type(DistFuncType), dimension(:) :: phi_conv

    type(DistFuncType), dimension(:), intent(in) :: phi
    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    type(MoleculeType), intent(in) :: molecule
    double precision, intent(in) :: deltat
    integer, intent(in) :: species

    double precision :: deltax, deltay
    integer :: nx_space, ny_space

    double precision :: c, factor, x_vel
    double precision :: a_m2, a_m1, a_0, a_1, a_2
    double precision :: phi_0, phi_1, phi_2, phi_m1, phi_m2

    integer :: r_modes, v_modes

    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    integer :: r_levels, v_levels
    integer :: l, i, j, k, nx, gr, gr_1, gr_2, gr_m1, gr_m2, n

    integer, dimension(2) :: ind

    call get_delta_x( deltax, deltay )
    call get_nspace( nx_space, ny_space )

    factor = deltat/deltax

    r_modes = molecule%rot_modes
    v_modes = molecule%vib_modes

    do nx = 3, nx_space-2

       call get_spatial_reference( nx-2, gr_m2 )
       call get_spatial_reference( nx-1, gr_m1 )
       call get_spatial_reference( nx, gr )
       call get_spatial_reference( nx+1, gr_1 )
       call get_spatial_reference( nx+2, gr_2 )

       r_levels = phi(nx)%num_rot_levels
       v_levels = phi(nx)%num_vib_levels

       i_min = vel_grid(gr)%i_min
       i_max = vel_grid(gr)%i_max
       j_min = vel_grid(gr)%j_min
       j_max = vel_grid(gr)%j_max
       k_min = vel_grid(gr)%k_min
       k_max = vel_grid(gr)%k_max

       do k = k_min, k_max
          do j = j_min, j_max
             do i = i_min, i_max

                x_vel = vel_grid(gr)%x(i)
                c = factor*x_vel

                phi_m2 = phi(nx-2)%value(i,j,k)
                phi_m1 = phi(nx-1)%value(i,j,k)
                phi_0  = phi(nx)%value(i,j,k)
                phi_1  = phi(nx+1)%value(i,j,k)
                phi_2  = phi(nx+2)%value(i,j,k)

                a_m2 = -c*( one - c*c )*( two + c )/24.0d0
                a_m1 = c*( one + c )*( four - c*c )/6.0d0
                a_0  = ( one - c*c )*( four - c*c )/4.0d0
                a_1  = -c*( one - c )*( four - c*c )/6.0d0
                a_2  = c*( one - c*c )*( two - c )/24.0d0

                phi_conv(nx)%value(i,j,k) = a_m2*phi_m2 + a_m1*phi_m1 + a_0*phi_0 + a_1*phi_1 + a_2*phi_2

                if( r_modes .gt. 0 )then
                   do l = 1, r_levels
                      phi_conv(nx)%rot(l,i,j,k) = &
                           a_m2*phi_conv(nx-2)%rot(l,i,j,k) + &
                           a_m1*phi_conv(nx-1)%rot(l,i,j,k) + &
                           a_0*phi_conv(nx)%rot(l,i,j,k) + &
                           a_1*phi_conv(nx+1)%rot(l,i,j,k) + &
                           a_2*phi_conv(nx+2)%rot(l,i,j,k)
                   end do
                end if

                if( v_modes .gt. 0 )then
                   do l = 1, v_levels
                      phi_conv(nx)%vib(l,i,j,k) = &
                           a_m2*phi_conv(nx-2)%vib(l,i,j,k) + &
                           a_m1*phi_conv(nx-1)%vib(l,i,j,k) + &
                           a_0*phi_conv(nx)%vib(l,i,j,k) + &
                           a_1*phi_conv(nx+1)%vib(l,i,j,k) + &
                           a_2*phi_conv(nx+2)%vib(l,i,j,k)
                   end do
                end if

             end do
          end do
       end do

    end do

    ind(1) = 2
    ind(2) = nx_space - 1

    do n = 1, 2
       nx = ind(n)

       call get_spatial_reference( nx, gr )
       call get_spatial_reference( nx+1, gr_1 )
       call get_spatial_reference( nx-1, gr_m1 )

       r_levels = phi(nx)%num_rot_levels
       v_levels = phi(nx)%num_vib_levels

       i_min = vel_grid(gr)%i_min
       i_max = vel_grid(gr)%i_max
       j_min = vel_grid(gr)%j_min
       j_max = vel_grid(gr)%j_max
       k_min = vel_grid(gr)%k_min
       k_max = vel_grid(gr)%k_max

       do k = k_min, k_max
          do j = j_min, j_max
             do i = i_min, i_max

                x_vel = vel_grid(gr)%x(i)
                c = factor*x_vel

                phi_0 = phi(nx)%value(i,j,k)

                if( x_vel .lt. zero )then

                   phi_1 = phi(nx+1)%value(i,j,k)

                   phi_conv(nx)%value(i,j,k) = phi_0 - c*( phi_1 - phi_0 )

                   if( r_modes .gt. 0 )then
                      do l = 1, r_levels
                         phi_conv(nx)%rot(l,i,j,k) = phi(nx)%rot(l,i,j,k) - &
                              c*( phi(nx+1)%rot(l,i,j,k) - phi(nx)%rot(l,i,j,k) )
                      end do
                   end if

                   if( v_modes .gt. 0 )then
                      do l = 1, v_levels
                         phi_conv(nx)%vib(l,i,j,k) = phi(nx)%vib(l,i,j,k) - &
                              c*( phi(nx+1)%vib(l,i,j,k) - phi(nx)%vib(l,i,j,k) )
                      end do
                   end if

                else

                   phi_m1 = phi(nx-1)%value(i,j,k)

                   phi_conv(nx)%value(i,j,k) = phi_0 - c*( phi_0 - phi_m1 )

                   if( r_modes .gt. 0 )then
                      do l = 1, r_levels
                         phi_conv(nx)%rot(l,i,j,k) = phi(nx)%rot(l,i,j,k) - &
                              c*( phi(nx)%rot(l,i,j,k) - phi(nx-1)%rot(l,i,j,k) )
                      end do
                   end if

                   if( v_modes .gt. 0 )then
                      do l = 1, v_levels
                         phi_conv(nx)%vib(l,i,j,k) = phi(nx)%vib(l,i,j,k) - &
                              c*( phi(nx)%vib(l,i,j,k) - phi(nx-1)%vib(l,i,j,k) )
                      end do
                   end if

                end if

             end do
          end do
       end do

    end do

    call apply_boundary_conditions( phi, phi_conv, molecule, vel_grid, factor, species )

    return
  end subroutine fourth_order_scheme

end module Convection
