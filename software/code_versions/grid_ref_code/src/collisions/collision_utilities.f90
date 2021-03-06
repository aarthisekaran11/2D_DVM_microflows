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
module CollisionUtilities

  use Constants

  implicit none

  private

  ! CumulativeDF is here because multiple collision solver methods use the
  ! data type and CollisionUtilities is accessible to all
  type CumulativeDF
     double precision :: neg_dens
     double precision :: neq_dens, rot_neq_dens, vib_neq_dens, eq_dens
     integer, dimension(1:1000,2) :: search_ref
     double precision, allocatable, dimension(:) :: cumul_df
     double precision, allocatable, dimension(:) :: kin_df, vib_df, rot_df
     double precision, allocatable, dimension(:) :: rot_eq, vib_eq
     double precision, allocatable, dimension(:) :: coeff, coeff2
     double precision, allocatable, dimension(:) :: sign

     ! Properties used to test for need to calculate a new Maxwellian
     double precision :: dens, u, v, w, temp
  end type CumulativeDF

  integer :: rot_method
  integer :: vib_method

  integer, parameter :: rigid_rotor_LB    = 1
  integer, parameter :: nonrigid_rotor_LB = 2

  integer, parameter :: SHO_LB = 1
  integer, parameter :: AHO_LB = 2

  integer, parameter :: file_read = 0

  integer, parameter :: homonuclear   = 1
  integer, parameter :: heteronuclear = 2

  integer :: Z_rot_method, Z_vib_method
  integer, parameter :: constant       = 0
  integer, parameter :: parker         = 1
  integer, parameter :: millikan_white = 1
  integer, parameter :: bird           = 2

  ! TODO: discrete translational collision probabilities test
  double precision, dimension(:), allocatable :: Et, PEt, PsiEt
  double precision :: PEt_max
  integer :: Et_length, PEt_imax

  public :: CumulativeDF

  ! Collision subroutines
  public :: compute_num_coll_pairs
  public :: pick_collision_partner
  public :: compute_relative_velocity
  public :: compute_center_of_mass_velocity
  public :: find_post_collision_velocities
  public :: pick_random_velocity_on_sphere  

  ! Internal energy subroutines
  public :: set_energy_methods
  public :: set_relaxation_rate_method
  public :: split_depletion
  public :: energy_fraction_array
  public :: pick_level
  public :: calculate_g_prime
  public :: compute_rot_level_change
  public :: compute_vib_level_change
  public :: compute_rot_vib_level_change

  ! Array manipulation
  public :: normalized_array
  public :: one_sign_normalized_array
  public :: absolute_normalized_array
  public :: pos_neg_split_array

  ! Et probability test
  public :: build_Et_probability_array

contains

  subroutine compute_num_coll_pairs( num_colls, dt, densA, densB, temp, mass, coln_rms, vel_grid, n, m )

    use VelocityGrid

    implicit none

    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    double precision, intent(in) :: dt, densA, densB
    double precision, dimension(2), intent(in) :: temp, mass
    double precision, dimension(:,:), intent(in) :: coln_rms
    integer, intent(in) :: n, m
    
    integer :: num_colls

    double precision :: rms_sq, beta3, beta3A, beta3B, coll_temp

    ! Find rms^2 value
    rms_sq = coln_rms(n,m) * coln_rms(n,m)

    ! Find adjustment for grid spacing
!!$    beta3 = vel_grid(n)%beta3_min
!!$    if( vel_grid(m)%beta3_min .lt. beta3 ) beta3 = vel_grid(m)%beta3_min

    beta3A = vel_grid(n)%beta3_avg * mass(1)**( 3.0d0/2.0d0 )
    beta3B = vel_grid(m)%beta3_avg * mass(2)**( 3.0d0/2.0d0 )
    beta3  = one_half * ( beta3A + beta3B )

    ! Find average temperature
    coll_temp = one_half * ( temp(1) + temp(2) )
    coll_temp = coll_temp**( 2.0d0 / 3.0d0 )

    num_colls = nint( ( one_half * dt * densA * densB )/( rms_sq * beta3 ) )

    return
  end subroutine compute_num_coll_pairs

  subroutine pick_collision_partner( i, j, k, global_index, psi, ref, vel_grid )

    use VelocityGrid
    use SortAndSearch
    use RandomNumberGeneration

    implicit none

    integer :: i, j, k, global_index

    type(VelocityGridType), intent(in) :: vel_grid
    double precision, dimension(:), intent(in) :: psi
    integer, dimension(:,:), intent(in) :: ref

    double precision :: Rf

    integer :: num_points_x, num_points_y, num_points, l
    integer :: i_min, j_min, k_min

    i_min = vel_grid%i_min
    j_min = vel_grid%j_min
    k_min = vel_grid%k_min

    num_points_x = vel_grid%num_points_x
    num_points_y = vel_grid%num_points_y
    num_points   = vel_grid%num_points

    ! Pick a scaled random number from psi
    Rf = get_rand()

    ! Do binary search for global index
!!$    call binary_search( psi, num_points, Rf*psi( num_points ), l )

    call referenced_search( psi, ref, num_points, Rf, l )

    ! Map back to local coordinates
    call global2local_map( l, i_min, j_min, k_min, num_points_x, num_points_y, i, j, k )

    global_index = l + 1

    return
  end subroutine pick_collision_partner

  subroutine compute_relative_velocity( i, j, k, g, vel_grid, n, m )

    use VelocityGrid
    
    implicit none

    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    integer, dimension(2), intent(in) :: i, j, k
    integer, intent(in) :: n, m

    double precision, dimension(2) :: x, y, z
    double precision :: g, gx, gy, gz, g_sq

    x(1) = vel_grid(n)%x( i(1) )
    y(1) = vel_grid(n)%y( j(1) )
    z(1) = vel_grid(n)%z( k(1) )

    x(2) = vel_grid(m)%x( i(2) )
    y(2) = vel_grid(m)%y( j(2) )
    z(2) = vel_grid(m)%z( k(2) )

    gx = x(2) - x(1)
    gy = y(2) - y(1)
    gz = z(2) - z(1)

    g_sq = gx * gx + gy * gy + gz * gz

    g = sqrt( g_sq )

    return
  end subroutine compute_relative_velocity

  subroutine compute_center_of_mass_velocity( x, y, z, mass, v_com )

    implicit none

    double precision, dimension(3) :: v_com

    double precision, dimension(2), intent(in) :: x, y, z, mass

    double precision, dimension(2) :: mass_frac

    mass_frac(1) = mass(1)/(mass(1) + mass(2))
    mass_frac(2) = mass(2)/(mass(1) + mass(2))

    v_com(1) = mass_frac(1)*x(1) + mass_frac(2)*x(2)
    v_com(2) = mass_frac(1)*y(1) + mass_frac(2)*y(2)
    v_com(3) = mass_frac(1)*z(1) + mass_frac(2)*z(2)

    return
  end subroutine compute_center_of_mass_velocity

  subroutine find_post_collision_velocities( xi, eta, x, y, z, v_com, g_prime, g_ratio, alpha )

    implicit none

    double precision, dimension(3) :: xi, eta

    double precision, dimension(3), intent(in) :: v_com
    double precision, dimension(2), intent(in) :: x, y, z
    double precision, intent(in) :: g_prime, g_ratio, alpha

    double precision, dimension(2) :: gx, gy, gz
    double precision, dimension(2) :: rgx, rgy, rgz

    gx(1) = abs( x(1) - v_com(1) )*g_ratio
    gy(1) = abs( y(1) - v_com(2) )*g_ratio
    gz(1) = abs( z(1) - v_com(3) )*g_ratio

    gx(2) = abs( x(2) - v_com(1) )*g_ratio
    gy(2) = abs( y(2) - v_com(2) )*g_ratio
    gz(2) = abs( z(2) - v_com(3) )*g_ratio

    call pick_random_velocity_on_sphere( g_prime, rgx, rgy, rgz, gx, gy, gz, alpha )

    xi(1) = v_com(1) + rgx(1)
    xi(2) = v_com(2) + rgy(1)
    xi(3) = v_com(3) + rgz(1)

    eta(1) = v_com(1) + rgx(2)
    eta(2) = v_com(2) + rgy(2)
    eta(3) = v_com(3) + rgz(2)

    return
  end subroutine find_post_collision_velocities

  subroutine pick_random_velocity_on_sphere( g, rg_x, rg_y, rg_z, gx, gy, gz, alpha )

    use RandomNumberGeneration

    implicit none

    double precision, intent(in) :: g, alpha
    double precision, dimension(2), intent(in) :: gx, gy, gz

    double precision, dimension(2) :: rg_x, rg_y, rg_z

    double precision, dimension(2) :: g_part
    double precision :: ur, vr, wr, r_x, r_y, r_z, g_ratio
    double precision :: phi, ctheta, stheta, cphi, sphi
    double precision :: denominator, ralpha

    ! set local, changeable variables
    ur = gx(1) + gx(2)
    vr = gy(1) + gy(2)
    wr = gz(1) + gz(2)

    g_part(1) = sqrt( gx(1)*gx(1) + gy(1)*gy(1) + gz(1)*gz(1) )
    g_part(2) = sqrt( gx(2)*gx(2) + gy(2)*gy(2) + gz(2)*gz(2) )

    phi = two*pi*get_rand()
    cphi = cos(phi)
    sphi = sin(phi)

    ralpha = one/alpha

    if( abs( ralpha - one ) .gt. double_tol )then

       ctheta = two*get_rand()**ralpha - one
       stheta = sqrt( one - ctheta*ctheta )
       denominator = sqrt( vr*vr + wr*wr )

       if( denominator .gt. double_tol*g )then
          r_x = ctheta*ur + stheta*sphi*denominator
          r_y = ctheta*vr + stheta*( g*wr*cphi - ur*vr*sphi )/denominator
          r_z = ctheta*wr - stheta*( g*vr*cphi + ur*wr*sphi )/denominator
       else
          r_x = ctheta*ur
          r_y = stheta*cphi*ur
          r_z = stheta*sphi*ur
       end if

    else

       ctheta = two*get_rand() - one
       stheta = sqrt( one - ctheta*ctheta )

       r_x = g*stheta*cphi
       r_y = g*stheta*sphi
       r_z = g*ctheta

    end if

    g_ratio = g_part(1)/g
    rg_x(1) = g_ratio*r_x
    rg_y(1) = g_ratio*r_y
    rg_z(1) = g_ratio*r_z

    g_ratio = g_part(2)/g
    rg_x(2) = -g_ratio*r_x
    rg_y(2) = -g_ratio*r_y
    rg_z(2) = -g_ratio*r_z

    return
  end subroutine pick_random_velocity_on_sphere

  subroutine set_energy_methods( rot_method_in, vib_method_in )

    implicit none

    integer, intent(in) :: rot_method_in, vib_method_in

    rot_method = rot_method_in
    vib_method = vib_method_in

    return
  end subroutine set_energy_methods

  subroutine set_relaxation_rate_method( Z_rot_method_in, Z_vib_method_in )

    implicit none

    integer, intent(in) :: Z_rot_method_in, Z_vib_method_in

    integer :: status

    Z_rot_method = Z_rot_method_in
    Z_vib_method = Z_vib_method_in

    return
  end subroutine set_relaxation_rate_method

  subroutine split_depletion( depl_frac, m_red, properties, molecule, n, m )

    use PhysicalProperties
    use SpeciesAndReferenceData

    implicit none

    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(PropertiesType), intent(in) :: properties
    integer, intent(in) :: n, m
    double precision, intent(in) :: m_red

    double precision, dimension(5) :: depl_frac

    double precision :: dens, temp, press, avg_speed
    double precision :: T_star, Zr_inf
    double precision :: Av, Bv, coll_freq
    double precision :: dof, omega_b, vib_cs, p_Zvscale
    double precision :: C1, C2

    double precision :: trn_dof, rot_dof, vib_dof, denom, omega

    double precision :: Zr1, Zr2, Zv1, Zv2
    double precision :: Pr1, Pr2, Pv1, Pv2
    double precision :: f_inelastic, f_rot, f_vib

    integer :: r_modes, v_modes

    temp = properties%mix_tr_temp

    omega = 0.5d0*( molecule(n)%omega + molecule(m)%omega )
    trn_dof = 5.0d0 - 2.0d0*omega

    r_modes = molecule(n)%rot_modes

    ! Adjustments to Zr1, Zr2, Zv1, Zv2 are made based on the work of ...

    if( r_modes .gt. 0 )then
       rot_dof = properties%rot_dof(n)
       select case( Z_rot_method )
       case( constant )
          Zr1 = molecule(n)%cZr

       case( parker )
          Zr_inf = molecule(n)%Zr_inf
          T_star = molecule(n)%T_star

          Zr1 = Zr_inf*((trn_dof+rot_dof)/trn_dof)/&
               ( one + sqrt( one_fourth*pi*pi*pi*T_star/temp ) + ( one_fourth*pi*pi + pi )*T_star/temp )

       case default
          write(*,*) "Error: Invalid value of Z_rot_method: ", Z_rot_method
          stop

       end select

       if( Zr1 .lt. 3.0d0 )then
          Zr1 = 3.0d0
       end if

       Pr1 = (one/Zr1)

    else
       Pr1 = zero

    end if
    
    r_modes = molecule(m)%rot_modes

    if( r_modes .gt. 0 )then
       rot_dof = properties%rot_dof(m)
       select case( Z_rot_method )
       case( constant )
          Zr2 = molecule(m)%cZr

       case( parker )
          Zr_inf = molecule(m)%Zr_inf
          T_star = molecule(m)%T_star

          Zr2 = Zr_inf*((trn_dof+rot_dof)/trn_dof)/&
               ( one + one_half * sqrt( pi*pi*pi*T_star/temp ) + ( one_fourth*pi*pi + pi )*T_star/temp )

       case default
          write(*,*) "Error: Invalid value of Z_rot_method: ", Z_rot_method
          stop

       end select

       if( Zr2 .lt. 3.0d0 )then
          Zr2 = 3.0d0
       end if

       Pr2 = (one/Zr2)

    else
       Pr2 = zero

    end if

    v_modes = molecule(n)%vib_modes

    if( v_modes .gt. 0 )then
       vib_dof = properties%vib_dof(n)*properties%vib_dof(n)*exp(molecule(n)%theta_v/temp)*0.5d0
       select case( Z_vib_method )
       case( constant )
          Zv1 = molecule(n)%cZv

       case( millikan_white )
          Av        = molecule(n)%Av*m_red**(0.5d0)
          Bv        = molecule(n)%Bv*m_red**(0.75d0) - 18.42d0
          vib_cs    = molecule(n)%vib_cs
          p_Zvscale = molecule(n)%p_Zvscale
          coll_freq = properties%coll_freq(n)/m_red**(0.5d0)
          dof       = properties%non_vib_dof(n)
          press     = properties%pressure(n)
          avg_speed = properties%mix_avg_speed
          dens      = properties%dens(n)

          Zv1 = coll_freq*( exp(Av*temp**(-one_third)+Bv)*p_Zvscale/press + one/(avg_speed*vib_cs*dens) )&
               *((trn_dof+vib_dof)/trn_dof)

       case( bird )
          C1 = molecule(n)%Zv_1
          C2 = molecule(n)%Zv_2
          omega_b = molecule(n)%omega

          Zv1 = C1*exp( C2*temp**(-one_third) )/temp**omega_b

       case default
          write(*,*) "Error: Invalid value of Z_vib_method: ", Z_vib_method
          stop

       end select

       Pv1 = (one/Zv1)

    else
       Pv1 = zero

    end if

    v_modes = molecule(m)%vib_modes

    if( v_modes .gt. 0 )then
       vib_dof = properties%vib_dof(m)*properties%vib_dof(m)*exp(molecule(m)%theta_v/temp)*0.5d0
       select case( Z_vib_method )
       case( constant )
          Zv2 = molecule(m)%cZv

       case( millikan_white )
          Av        = molecule(m)%Av*m_red**(0.5d0)
          Bv        = molecule(m)%Bv*m_red**(0.75d0) - 18.42d0
          vib_cs    = molecule(m)%vib_cs
          p_Zvscale = molecule(m)%p_Zvscale
          coll_freq = properties%coll_freq(m)/m_red**(0.5d0)
          dof       = properties%non_vib_dof(m)
          press     = properties%pressure(m)
          avg_speed = properties%mix_avg_speed
          dens      = properties%dens(m)

          Zv2 = coll_freq*( exp(Av*temp**(-one_third)+Bv)*p_Zvscale/press + one/(avg_speed*vib_cs*dens) )&
               *((trn_dof+vib_dof)/trn_dof)

       case( bird )
          C1 = molecule(m)%Zv_1
          C2 = molecule(m)%Zv_2
          omega_b = molecule(m)%omega

          Zv2 = C1*exp( C2*temp**(-one_third) )/temp**omega_b

       case default
          write(*,*) "Error: Invalid value of Z_vib_method: ", Z_vib_method
          stop

       end select

       Pv2 = (one/Zv2)

    else
       Pv2 = zero

    end if

    depl_frac(1) = one - Pr1 - Pr2 - Pv1 - Pv2
    depl_frac(2) = Pr1
    depl_frac(3) = Pr2
    depl_frac(4) = Pv1
    depl_frac(5) = Pv2

!!$    f_inelastic = one - one_half*( Pr1 + Pr2 )*one_half*( Pv1 + Pv2 )
!!$    f_rot = one - one_half*( Pr1 + Pr2 )
!!$    f_vib = one - one_half*( Pv1 + Pv2 )
!!$
!!$    depl_frac(1) = one - f_inelastic
!!$    depl_frac(2) = f_inelastic - f_vib
!!$    depl_frac(3) = f_inelastic - f_rot
!!$    depl_frac(4) = f_rot + f_vib - f_inelastic

    return
  end subroutine split_depletion

  subroutine calculate_g_prime( g_prime, g, d_E, m_red)

    implicit none

    double precision, intent(in) :: g, d_E, m_red
    
    double precision :: g_prime, squared

    squared = g*g - 2.0d0*d_E/m_red
    
    if( squared .lt. zero )then
       write(*,*)"Negative value of determinant in dE calculation"
       write(*,*)squared
       write(*,*)g*g, 2.0d0*d_E/m_red
       stop

    else
       g_prime = sqrt( squared ) ![sqrt(J/kg) = m/s]

    end if

    return
  end subroutine calculate_g_prime

  subroutine energy_fraction_array( f, array, levels, sgn )

    implicit none

    double precision, dimension(:), intent(in) :: array
    double precision, intent(in) :: sgn
    integer, intent(in) :: levels

    double precision, dimension(:) :: f
    double precision :: total

    integer :: n

    f = array
!!$    where( sgn*f .lt. zero ) f = zero
!!$    if( sum(f) .eq. zero )then
!!$       write(*,*)sum(f)
!!$       write(*,*)sgn
!!$       write(*,*)f
!!$       write(*,*)array
!!$       stop
!!$    end if
!!$    total = one/sum(f)
    total = one/sum(abs(f))

    f = f*total

    return
  end subroutine energy_fraction_array

  subroutine normalized_array( f, array, levels )

    implicit none

    double precision, dimension(:), intent(in) :: array
    integer, intent(in) :: levels

    double precision, dimension(:) :: f
    double precision :: total

    f = array

    total = sum(f)
    total = one/total

    f = f*total

    return
  end subroutine normalized_array

  subroutine absolute_normalized_array( f, array, levels )

    implicit none

    double precision, dimension(:), intent(in) :: array
    integer, intent(in) :: levels

    double precision, dimension(:) :: f
    double precision :: total

    f = array

    total = one/sum(abs(f))

    f = f*total

    return
  end subroutine absolute_normalized_array

  subroutine one_sign_normalized_array( f, array, levels, sgn )
    
    implicit none

    double precision, dimension(:), intent(in) :: array
    double precision, intent(in) :: sgn
    integer, intent(in) :: levels

    double precision, dimension(:) :: f
    double precision :: total

    f = array

    where( sgn*f .lt. zero ) f = zero

    total = sum(f)

    if( total .eq. zero )then
       total = zero
    else
       total = one/total
    end if

    f = f*total

    return
  end subroutine one_sign_normalized_array

  subroutine pos_neg_split_array( fpos, fneg, array, levels )
    
    implicit none

    double precision, dimension(:), intent(in) :: array
    integer, intent(in) :: levels

    double precision, dimension(:) :: fpos, fneg

    where( array .lt. zero )
       fpos = zero
       fneg = abs(array)
    elsewhere
       fpos = abs(array)
       fneg = zero
    end where

    return
  end subroutine pos_neg_split_array

  subroutine pick_level( l, array, phi, N )

    use RandomNumberGeneration
    use MathUtilities

    implicit none

    integer :: l

    double precision, dimension(:), intent(in) :: array
    double precision, intent(in) :: phi
    integer, intent(in) :: N

    double precision, allocatable, dimension(:) :: f
    double precision :: sgn

    double precision :: Rf, test, total

    total = sum(abs(array))    
    Rf = get_rand()*total

    test = zero

    do l = 1, N-1
       test = test + abs(array(l))
       if( test .ge. Rf ) return
    end do

    return
  end subroutine pick_level
  
  subroutine compute_rot_level_change( Jp, d_E, J, g, omega, m_red, mass, levels, &
       molecule, level_array )

    use SpeciesAndReferenceData

    implicit none

    type(MoleculeType), intent(in) :: molecule
    integer, intent(in) :: J, levels
    double precision, intent(in) :: m_red, g, omega, mass
    double precision, dimension(:), intent(in) :: level_array

    integer :: modes, type
    double precision :: theta, Y01, Y02

    integer :: Jp
    double precision :: d_E, E_c, E_t

    theta = molecule%theta_r
    modes = molecule%rot_modes
    type  = molecule%molecule_type
    Y01   = molecule%Y01
    Y02   = molecule%Y02

    select case( rot_method )
    case( file_read )


    case( rigid_rotor_LB )
       E_t = one_half*m_red*g*g
       E_c = E_t
       call rot_rigid_rotor_LB( Jp, E_c, J, omega, mass, theta, modes, levels, level_array, type )
       d_E = E_t - E_c

    case( nonrigid_rotor_LB )
       E_t = one_half*m_red*g*g
       E_c = E_t
       call rot_nonrigid_rotor_LB( Jp, E_c, J, omega, mass, theta, Y01, Y02, &
            modes, levels, level_array, type )
       d_E = E_t - E_c

    case default
       write(*,*) "Error: Invalid value of rot_method: ", rot_method
       stop

    end select

    return
  end subroutine compute_rot_level_change

  subroutine compute_vib_level_change( vp, d_E, v, g, omega, m_red, mass, levels, &
       molecule, level_array )

    use SpeciesAndReferenceData

    implicit none

    type(MoleculeType), intent(in) :: molecule
    integer, intent(in) :: v, levels
    double precision, intent(in) :: m_red, g, omega, mass
    double precision, dimension(:), intent(in) :: level_array

    integer :: vp
    double precision :: d_E

    integer :: modes
    double precision :: theta, Y10, Y20
    double precision :: E_c, E_t

    theta = molecule%theta_v
    modes = molecule%vib_modes
    Y10   = molecule%Y10
    Y20   = molecule%Y20
    
!!$    E_t = one_half*m_red*g*g
!!$    E_c = E_t
!!$    call vib_discrete_LB( vp, E_c, v, omega, mass, theta, modes, levels, level_array )
!!$    d_E = E_t - E_c
!!$
    select case( vib_method )
    case( file_read )


    case( SHO_LB )
       E_t = one_half*m_red*g*g
       E_c = E_t
       call vib_SHO_LB( vp, E_c, v, omega, mass, theta, modes, levels, level_array )
       d_E = E_t - E_c
       
    case( AHO_LB )
       E_t = one_half*m_red*g*g
       E_c = E_t
       call vib_AHO_LB( vp, E_c, v, omega, mass, Y10, Y20, modes, levels, level_array )
       d_E = E_t - E_c

    case default
       write(*,*) "Error: Invalid value of vib_method: ", vib_method
       stop

    end select

    return
  end subroutine compute_vib_level_change

  subroutine compute_rot_vib_level_change( Jp1, Jp2, vp1, vp2, d_E, J1, J2, v1, v2, &
       g, omega, m_red, mass, theta_r, theta_v, r_modes, r_levels, v_modes, v_levels, mol_type1, mol_type2 )

    integer, intent(in) :: J1, J2, v1, v2, mol_type1, mol_type2
    integer, dimension(2), intent(in) :: r_modes, v_modes, r_levels, v_levels
    double precision, intent(in) :: m_red, g, omega
    double precision, dimension(2), intent(in) :: mass, theta_r, theta_v

    integer :: Jp1, Jp2, vp1, vp2
    double precision :: d_E

    double precision :: E_c, E_t

!!$    if( rot_method .eq. boltzmann .and. vib_method .eq. boltzmann )then
!!$       E_t = one_half*m_red*g*g
!!$       E_c = E_t
!!$       call vib_boltz_exchange( vp1, E_c, v1, omega, mass(1), theta_v(1), v_modes(1), v_levels(1) )
!!$       call rot_boltz_exchange( Jp1, E_c, J1, omega, mass(1), theta_r(1), r_modes(1), r_levels(1), mol_type1 )
!!$       call vib_boltz_exchange( vp2, E_c, v2, omega, mass(2), theta_v(2), v_modes(2), v_levels(2) )
!!$       call rot_boltz_exchange( Jp2, E_c, J2, omega, mass(2), theta_r(2), r_modes(2), r_levels(2), mol_type2 )
!!$       d_E = E_t - E_c
!!$    end if

    return
  end subroutine compute_rot_vib_level_change

  subroutine rot_rigid_rotor_LB( Jp, E_c, J, omega, mass, theta, modes, levels, level_array, mol_type )

    use RandomNumberGeneration

    implicit none

    integer, intent(in) :: J, modes, levels, mol_type
    double precision, dimension(:), intent(in) :: level_array
    double precision, intent(in) :: omega, mass, theta

    integer :: Jp
    double precision :: E_c

    double precision :: l, factor
    double precision :: energy_level, max_level, star_level, new_level
    double precision :: star_E_r, E_r
    double precision :: R, P
    logical :: test = .true.

    if( modes .eq. 0 ) return

    select case( mol_type )
    case( homonuclear )

       if( mod((J-1),2) .eq. 0 )then
          ! Even numbered levels
          l = one_half*dble( J - 1 )
          energy_level = level_array(J)
          E_c = E_c + mass*energy_level

          max_level = floor( one_fourth*( -one + sqrt( one + 8.0d0*E_c/( mass*theta ) ) ) )
          if( (levels-1)/2 .lt. max_level ) max_level = floor( one_half*dble(levels-1) )

          star_level = ( one_fourth*( -one + &
               sqrt( (one + 8.0d0*E_c/(mass*theta) )/( 4.0d0 - two*omega ) ) ) )
          star_level = dble(nint(star_level))
          if( max_level .lt. star_level ) star_level = max_level

          energy_level = level_array( 2*nint(star_level) + 1 )
          star_E_r = mass*energy_level

          test = .true.
          do while( test )
             new_level = floor( get_rand()*( max_level + cnst99 ) )

             Jp = 2*nint(new_level) + 1
             energy_level = level_array( Jp )
             E_r = mass*energy_level

             R = get_rand()
             P = ( ( 4.0d0*new_level + one )/( 4.0d0*star_level + one ) )*&
                  ( ( E_c - E_r )/( E_c - star_E_r ) )**( 1.5d0 - omega )

             if( R .lt. P ) test = .false.

          end do

       else
          ! Odd numbered levels
          l = one_half*dble( J - 2 )
          energy_level = level_array(J)
          E_c = E_c + mass*energy_level

          max_level = floor( one_fourth*( -3.0d0 + sqrt( one + 8.0d0*E_c/( mass*theta ) ) ) )
          if( (levels-2)/2 .lt. max_level ) max_level = floor( one_half*dble(levels-2) )

          star_level = ( one_fourth*( -3.0d0 + &
               sqrt( ( one + 8.0d0*E_c/(mass*theta) )/( 4.0d0 - two*omega ) ) ) )
          star_level = dble(nint(star_level))
          if( max_level .lt. star_level ) star_level = max_level

          energy_level = level_array( 2*nint(star_level) + 2 )
          star_E_r = mass*energy_level

          test = .true.
          do while( test )
             new_level = floor( get_rand()*( max_level + cnst99 ) )

             Jp = 2*nint(new_level) + 2
             energy_level = level_array(Jp)
             E_r = mass*energy_level

             R = get_rand()
             P = ( ( 4.0d0*new_level + 3.0d0 )/( 4.0d0*star_level + 3.0d0 ) )*&
                  ( ( E_c - E_r )/( E_c - star_E_r ) )**( 1.5d0 - omega )

             if( R .lt. P ) test = .false.

          end do

       end if

    case( heteronuclear )

       l = dble( J - 1 )
       energy_level = level_array(J)
       E_c = E_c + mass*energy_level

       max_level = floor( one_half*( -one + sqrt( one + 8.0d0*E_c/( mass*theta ) ) ) )
       if( levels - 1 .lt. max_level ) max_level = dble( levels - 1 )

       star_level = ( one_half*( -one + &
            sqrt( ( one + 8.0d0*E_c/( mass*theta ) )/( 4.0d0 - 2.0d0*omega ) ) ) )
       star_level = dble(nint(star_level))
       if( max_level .lt. star_level ) star_level = max_level

       energy_level = level_array( nint(star_level) + 1 )
       star_E_r = mass*energy_level

       test = .true.
       do while( test )
          new_level = floor( get_rand()*( max_level + cnst99 ) )

          Jp = nint(new_level) + 1
          energy_level = level_array(Jp)
          E_r = mass*energy_level
          
          R = get_rand()
          P = ( ( two*new_level + one )/( two*star_level + one ) )*&
               ( ( E_c - E_r )/( E_c - star_E_r ) )**( 1.5d0 - omega )

          if( R .lt. P ) test = .false.

       end do

    case default
       write(*,*) "Error: Invalid value of mol_type. Value is: ", mol_type
       stop

    end select

    ! Take post-collision rotational energy away from available collision energy
    E_c = E_c - E_r

!!$!if( mode(J-1,2) .eq. 0 )then
!!$!    factor = even_spin
!!$!    l = dble( J - 1 )*one_half
!!$!else
!!$!    factor = odd_spin
!!$!    l = dble( J - 2 )*one_half 
!!$!end if
!!$ 
!!$    l = dble( J - 1 )!*factor
!!$    energy_level = one_half*mass*l*( l + 1 )*theta ! TODO: already exists
!!$    E_c = E_c + energy_level
!!$
!!$    max_level = floor( one_half*( -one + sqrt( one + 8.0d0*E_c/( mass*theta ) ) ) )
!!$    if( levels - 1 .lt. max_level ) max_level = dble( levels - 1 )
!!$ 
!!$    star_level = floor( one_half*( -one + &
!!$         sqrt( ( one_half + 3.0d0*omega + (16.0d0*E_c/( mass*theta ) ) )/( 6.5d0 - omega ) ) ) )
!!$
!!$    star_level = floor( one_half*( -one + &
!!$         sqrt( ( one + ( 16.0d0*E_c/( mass*theta ) ) )/( 7.0d0 - 4.0d0*omega ) ) ) )
!!$
!!$    star_level = floor( one_half*abs( -one + &
!!$         sqrt( ( one + 8.0d0*E_c/( mass*theta ) )/( 4.0d0 - 2.0d0*omega ) ) ) )
!!$
!!$    if( max_level .lt. star_level ) star_level = max_level
!!$    star_E_r = one_half*star_level*( star_level + one )*mass*theta
!!$
!!$    test = .true.
!!$    do while( test )
!!$       new_level = floor( get_rand()*( max_level + cnst99 ) )
!!$       E_r = one_half*new_level*( new_level + one )*mass*theta
!!$
!!$       R = get_rand()
!!$       P = ( ( two*new_level + one )/( two*star_level + one ) )*&
!!$            ( ( E_c - E_r )/( E_c - star_E_r ) )**( 1.5d0 - omega )
!!$
!!$       if( R .lt. P )then
!!$          Jp = nint( new_level ) + 1 ! levels run from 0 -> max-1, indices run from 1 -> max
!!$          test = .false.
!!$       end if
!!$    end do
!!$
!!$    E_c = E_c - one_half*mass*new_level*( new_level + one )*theta

    return
  end subroutine rot_rigid_rotor_LB

  subroutine rot_nonrigid_rotor_LB( Jp, E_c, J, omega, mass, theta, Y01, Y02, &
       modes, levels, level_array, mol_type )

    use RandomNumberGeneration

    implicit none

    integer, intent(in) :: J, modes, levels, mol_type
    double precision, dimension(:), intent(in) :: level_array
    double precision, intent(in) :: omega, mass, theta, Y01, Y02

    integer :: Jp
    double precision :: E_c

    double precision :: alpha, determinant
    double precision :: l, factor
    double precision :: energy_level, max_level, star_level, new_level
    double precision :: star_E_r, E_r
    double precision :: R, P
    logical :: test = .true.

    if( modes .eq. 0 ) return

    alpha = Y01/Y02
    
    energy_level = level_array(J)
    E_c = E_c + mass*energy_level

    determinant = alpha*alpha + 4.0d0*E_c/Y02

    select case( mol_type )
    case( homonuclear )

       if( mod((J-1),2) .eq. 0 )then
          ! Even numbered levels
          l = one_half*dble( J - 1 )

          if( determinant .lt. zero )then
             max_level = dble( levels - 1 )
          else
             max_level = one_fourth*(-one + sqrt(one + two*alpha + two*sqrt(determinant)))
          end if

          if( (levels-1)/2 .lt. max_level ) max_level = floor( one_half*dble(levels-1) )

          star_level = ( one_fourth*( -one + &
               sqrt( (one + 8.0d0*E_c/(mass*theta) )/( 4.0d0 - two*omega ) ) ) )
          star_level = dble(nint(star_level))
          if( max_level .lt. star_level ) star_level = max_level

          energy_level = level_array( 2*nint(star_level) + 1 )
          star_E_r = mass*energy_level

          test = .true.
          do while( test )
             new_level = floor( get_rand()*( max_level + cnst99 ) )

             Jp = 2*nint(new_level) + 1
             energy_level = level_array( Jp )
             E_r = mass*energy_level

             R = get_rand()
             P = ( ( 4.0d0*new_level + one )/( 4.0d0*star_level + one ) )*&
                  ( ( E_c - E_r )/( E_c - star_E_r ) )**( 1.5d0 - omega )

             if( R .lt. P ) test = .false.

          end do

       else
          ! Odd numbered levels
          l = one_half*dble( J - 2 )

          if( determinant .lt. zero )then
             max_level = dble( levels - 1 )
          else
             max_level = one_fourth*(-3.0d0 + sqrt(one + two*alpha + two*sqrt(determinant)))
          end if

          if( (levels-2)/2 .lt. max_level ) max_level = floor( one_half*dble(levels-2) )

          star_level = ( one_fourth*( -3.0d0 + &
               sqrt( ( one + 8.0d0*E_c/(mass*theta) )/( 4.0d0 - two*omega ) ) ) )
          star_level = dble(nint(star_level))
          if( max_level .lt. star_level ) star_level = max_level

          energy_level = level_array( 2*nint(star_level) + 2 )
          star_E_r = mass*energy_level

          test = .true.
          do while( test )
             new_level = floor( get_rand()*( max_level + cnst99 ) )

             Jp = 2*nint(new_level) + 2
             energy_level = level_array(Jp)
             E_r = mass*energy_level

             R = get_rand()
             P = ( ( 4.0d0*new_level + 3.0d0 )/( 4.0d0*star_level + 3.0d0 ) )*&
                  ( ( E_c - E_r )/( E_c - star_E_r ) )**( 1.5d0 - omega )

             if( R .lt. P ) test = .false.

          end do
          
       end if

    case( heteronuclear )
       l = dble( J - 1 )

       if( determinant .lt. zero )then
          max_level = dble( levels - 1 )
       else
          max_level = one_half*(-one + sqrt(one + two*alpha + two*sqrt(determinant)))
       end if

       if( levels - 1 .lt. max_level ) max_level = dble( levels - 1 )

       star_level = ( one_half*( -one + &
            sqrt( ( one + 8.0d0*E_c/( mass*theta ) )/( 4.0d0 - 2.0d0*omega ) ) ) )
       star_level = dble(nint(star_level))
       if( max_level .lt. star_level ) star_level = max_level

       energy_level = level_array( nint(star_level) + 1 )
       star_E_r = mass*energy_level

       test = .true.
       do while( test )
          new_level = floor( get_rand()*( max_level + cnst99 ) )

          Jp = nint(new_level) + 1
          energy_level = level_array(Jp)
          E_r = mass*energy_level
          
          R = get_rand()
          P = ( ( two*new_level + one )/( two*star_level + one ) )*&
               ( ( E_c - E_r )/( E_c - star_E_r ) )**( 1.5d0 - omega )

          if( R .lt. P ) test = .false.

       end do

    case default
       write(*,*) "Error: Invalid value of mol_type. Value is: ", mol_type
       stop

    end select

    ! Take post-collision rotational energy away from available collision energy
    E_c = E_c - E_r

    return
  end subroutine rot_nonrigid_rotor_LB

  subroutine vib_SHO_LB( vp, E_c, v, omega, mass, theta, modes, levels, level_array )

    use RandomNumberGeneration

    implicit none

    integer, intent(in) :: v, modes, levels
    double precision, intent(in) :: omega, mass, theta
    double precision, dimension(:), intent(in) :: level_array

    integer :: vp
    double precision :: E_c, E_v

    double precision :: l
    double precision :: energy_level, max_level, new_level
    double precision :: R, P
    logical :: test

    if( modes .eq. 0 ) return

    l = dble( v - 1 )
    energy_level = level_array(v)
    E_c = E_c + mass*energy_level

    max_level = floor( 2.0d0*E_c/( mass*theta ) )
    if( levels - 1 .lt. max_level ) max_level = dble( levels - 1 )

    test = .true.
    do while( test )
       new_level = floor( get_rand()*( max_level + cnst99 ) )

       vp = nint(new_level) + 1
       energy_level = level_array(vp)
       E_v = mass*energy_level

       R = get_rand()
       P = ( one - E_v/E_c )**( 1.5d0 - omega )

       if( R .lt. P ) test = .false.

    end do

    E_c = E_c - E_v

    return
  end subroutine vib_SHO_LB

  subroutine vib_AHO_LB( vp, E_c, v, omega, mass, Y10, Y20, modes, levels, level_array )

    use RandomNumberGeneration

    implicit none

    integer, intent(in) :: v, modes, levels
    double precision, intent(in) :: omega, mass, Y10, Y20
    double precision, dimension(:), intent(in) :: level_array

    integer :: vp
    double precision :: E_c

    double precision :: E_v, l
    double precision :: energy_level, max_level, new_level
    double precision :: determinant
    double precision :: R, P
    logical :: test

    if( modes .eq. 0 ) return

    l = dble( v - 1 )
    energy_level = level_array(v)
    E_c = E_c + mass*energy_level

    determinant = (Y20+Y10)*(Y20+Y10) + four*Y20*E_c
    if( determinant .lt. zero )then
       max_level = dble( levels - 1 )
    else
       max_level = floor(( -(Y20 + Y10) + sqrt(determinant) )/(two*Y20))
    end if

    if( levels - 1 .lt. max_level ) max_level = dble( levels - 1 )

    test = .true.
    do while( test )
       R = get_rand()
       new_level = floor( R*( max_level + cnst99 ) )

       vp = nint(new_level) + 1
       energy_level = level_array(vp)
       E_v = mass*energy_level

       R = get_rand()
       P = ( one - E_v/E_c )**( 1.5d0 - omega )

       if( R .lt. P ) test = .false.

    end do

    E_c = E_c - E_v

    return
  end subroutine vib_AHO_LB

  ! TODO: testing discrete translational collision probabilities method
  subroutine build_Et_probability_array( T, n, m, u, v, w, omega, vel_grid )

    ! This code is currently hardcoded to a uniform grid

    use DistFunc
    use VelocityGrid
    use SortAndSearch

    implicit none

    type(VelocityGridType), intent(in) :: vel_grid
    double precision, intent(in) :: T, n, m, u, v, w, omega
    
    double precision, dimension(:,:,:), allocatable :: phi
    double precision, dimension(:), allocatable :: Et_full

    integer :: num_points_x, num_points_y, num_points_z
    integer :: num_points
    integer :: i_min, i_max, j_min, j_max, k_min, k_max
    
    double precision :: beta_x, beta_y, beta_z, beta3
    double precision :: x_min, y_min, z_min, x, y, z
    double precision :: c_x, c_y, c_z, c2
    double precision :: g_x, g_y, g_z, E

    double precision :: coeff

    integer :: i, j, k
    integer :: ii, jj, kk
    integer :: l, ll, p

    integer :: count

    double precision :: sum

    write(*,*) "Building discrete E_t arrays..."

    ! Grid data
    num_points_x = vel_grid%num_points_x
    num_points_y = vel_grid%num_points_y
    num_points_z = vel_grid%num_points_z
    num_points   = num_points_x * num_points_y * num_points_z

    i_min = vel_grid%i_min
    i_max = vel_grid%i_max
    j_min = vel_grid%j_min
    j_max = vel_grid%j_max
    k_min = vel_grid%k_min
    k_max = vel_grid%k_max

    ! Allocate local arrays
    allocate( Et_full( 1:num_points ) )
    allocate( phi( i_min:i_max, j_min:j_max, k_min:k_max ) )

    ! Maxwell coefficient
    coeff = n * sqrt( m * m * m )/sqrt( pi * pi * pi * T * T * T )

    ! Minimum velocities for g calculation
    x_min = vel_grid%x(i_min)
    y_min = vel_grid%y(j_min)
    z_min = vel_grid%z(k_min)

    do l = 1, num_points

       ! Local index values
       call global2local_map( l-1, i_min, j_min, k_min, num_points_x, num_points_y, i, j, k )

       ! Velcoity and dV
       x = vel_grid%x(i)
       y = vel_grid%y(j)
       z = vel_grid%z(k)

       beta_x = vel_grid%beta_x(i)
       beta_y = vel_grid%beta_y(j)
       beta_z = vel_grid%beta_z(k)
       beta3  = beta_x * beta_y * beta_z

       ! Maxwellian
       c_x = x * x
       c_y = y * y
       c_z = z * z
       c2  = c_x + c_y + c_z

       phi(i,j,k) = coeff * exp( -c2 * m / T ) * beta3

       ! Relative velocity
       g_x = ( x - x_min ) * ( x - x_min )
       g_y = ( y - y_min ) * ( y - y_min )
       g_z = ( z - z_min ) * ( z - z_min )

       Et_full(l) = one_half * m * ( g_x + g_y + g_z )

    end do

    ! Sort relative velocity array
    call ssort1( Et_full, num_points )

    ! Count unique values
    count = 1
    do i = 2, num_points
       if( abs( Et_full(i) - Et_full(i-1) ) .lt. double_tol )then
          cycle
       else
          count = count + 1
       end if
    end do
    Et_length = count

    ! Allocate global arrays
    allocate( Et(1:Et_length), PEt(1:Et_length), PsiEt(0:Et_length) )

    ! Set Et
    count = 1
    Et(count) = Et_full(count)
    do i = 2, num_points
       if( abs( Et_full(i) - Et_full(i-1) ) .lt. double_tol )then
          cycle
       else
          count = count + 1
          Et(count) = Et_full(i)
       end if
    end do

    ! Build PEt
    PEt = zero
    sum = zero

    do l = 1, num_points

       call global2local_map( l-1, i_min, j_min, k_min, num_points_x, num_points_y, i, j, k )

       x = vel_grid%x(i)
       y = vel_grid%y(j)
       z = vel_grid%z(k)

       do ll = 1, num_points

          call global2local_map( ll-1, i_min, j_min, k_min, num_points_x, num_points_y, ii, jj, kk )

          g_x = ( x - vel_grid%x(ii) ) * ( x - vel_grid%x(ii) )
          g_y = ( y - vel_grid%y(jj) ) * ( y - vel_grid%y(jj) )
          g_z = ( z - vel_grid%z(kk) ) * ( z - vel_grid%z(kk) )
          E   = one_half * m * ( g_x + g_y + g_z )

          do p = 1, Et_length
             if( abs( E - Et(p) ) .lt. double_tol )then
                PEt(p) = PEt(p) + phi(i,j,k) * phi(ii,jj,kk) * ( g_x + g_y + g_z )**( one - omega )
                sum = sum + phi(i,j,k) * phi(ii,jj,kk) * ( g_x + g_y + g_z )**( one - omega )
                exit
             end if
          end do

       end do
    end do

    ! Find maximum and PsiEt
    !PEt(1) = zero ! remove probability of zero relative velocity
    PEt_max = zero
    PEt_imax = 1
    PsiEt(0) = zero

    do p = 1, Et_length
       PEt(p) = PEt(p)*exp( 2.0d0 * Et(p) / T )!/sum!/exp( -2.0d0 * Et(p) / T )
       PsiEt(p) = PsiEt(p-1) + PEt(p)/sum

       if( PEt(p) .gt. PEt_max )then
          PEt_max  = PEt(p)
          PEt_imax = p
       end if
    end do
!!$    PEt = PEt/PEt_max

    deallocate( phi, Et_full )

    write(*,*)"finished"

    return
  end subroutine build_Et_probability_array

  subroutine vib_discrete_LB( vp, E_c, v, omega, mass, theta, modes, levels, level_array )
    
    use RandomNumberGeneration
    use SortAndSearch

    implicit none

    integer, intent(in) :: v, modes, levels
    double precision, intent(in) :: omega, mass, theta
    double precision, dimension(:), intent(in) :: level_array

    integer :: vp, i, j
    double precision :: E_c, E_v, E_t

    double precision :: l, energy_level
    double precision :: max_level, new_level
    double precision :: P, R
    logical :: test
    double precision :: sum, max

    if( modes .eq. 0 )return

    l = dble( v - 1 )
    energy_level = level_array(v)
    E_c = E_c + mass * energy_level

    max_level = floor( 2.0d0 * E_c / ( mass * theta ) )
    if( levels - 1 .lt. max_level ) max_level = dble( levels - 1 )

    max = zero
    sum = zero
    j = 0
    do
       j = j + 1
       if( Et(j) .gt. E_c ) exit
       sum = sum + P
       if( PEt(j) .gt. max ) max = PEt(j)
    end do

    test = .true.
    do while( test )
       new_level = floor( get_rand() * ( max_level + cnst99 ) )
       vp = nint( new_level ) + 1
       
       energy_level = level_array( vp )
       E_v = mass * energy_level
       E_t = E_c - E_v

       i = (j-1)/2
       call binary_search_with_first_guess( i, E_t, Et, 1, j )
       
       P = PEt(i) + &
            ( PEt(i+1) - PEt(i) ) * &
            ( E_t - Et(i) ) / ( Et(i+1) - Et(i) )
       P = P / max
       R = get_rand()

       if( R .lt. P ) test = .false.

    end do

    E_c = E_c - E_v

    return
  end subroutine vib_discrete_LB

end module CollisionUtilities
