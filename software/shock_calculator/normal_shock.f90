program normal_shock

  implicit none

  double precision, parameter :: zero = 0.0d0
  double precision, parameter :: one = 1.0d0
  double precision, parameter :: two = 2.0d0
  double precision, parameter :: one_half = 0.5d0

  double precision, parameter :: boltz = 1.3806488e-23
  double precision, parameter :: avog = 6.022141e23

  integer :: file_unit = 1

  character(len=2) :: species, species_ref
  integer :: num_species
  double precision :: M_sh, n_up, T_up
  double precision :: mol_fraction

  integer :: s

  double precision :: mass, dof
  double precision :: theta_v, vib_dof
  double precision :: mol_mass, cp_mix, gamma

  double precision :: a_up, M_up

  integer :: stationary
  double precision :: u_shock, u_up

  double precision :: n_ratio, T_ratio, p_ratio
  double precision :: n_down, T_down, p_down
  double precision :: a_down, M_down, u_down

  double precision :: T_ref, m_ref, n_ref, eta_ref

  double precision :: x_min, x_max, y_min, y_max, z_min, z_max, test

  open( unit=file_unit, file="normal_shock_parameters.txt", status="unknown" )

  ! Reference properties
  write(*,*) "Input reference density [#/m^3]:"
  read(*,*) n_ref
  write(*,*) "Input reference temperature [K]:"
  read(*,*) T_ref
  write(*,*) "Input reference species:"
  read(*,*) species_ref
  call get_species_constants( species_ref, m_ref, dof, T_ref )
  write(*,*) ""

  eta_ref = sqrt( two * boltz * T_ref / m_ref )

  write(file_unit,fmt='(a8,e16.8,a8)') "n_ref = ", n_ref, " [#/m^3]"
  write(file_unit,fmt='(a8,e16.8,a4)') "T_ref = ", T_ref, " [K]"
  write(file_unit,fmt='(a14,a2)') "species_ref = ", species_ref
  write(file_unit,fmt='(a11,e16.8,a5)') "mass_ref = ", m_ref, " [kg]"
  write(file_unit,fmt='(a10,e16.8,a6)') "eta_ref = ", eta_ref, " [m/s]"
  write(file_unit,fmt='(x)')

  ! Upstream input properties
  write(*,*) "Upstream density [#/m^3]:"
  read(*,*) n_up
  write(*,*) "Upstream temperature [K]:"
  read(*,*) T_up
  write(*,*) "Upstream Mach number (shock reference frame):"
  read(*,*) M_sh
  write(*,*) ""

  ! Set up mixture
  mol_mass = zero
  cp_mix   = zero

  write(*,*) "Input number of species in mixture:"
  read(*,*) num_species
  write(*,*) ""

  write(*,*) "Input species name and mol fraction."
  write(*,*) "Species options include:"
  write(*,*) "He, Ne, Ar, Xe, H2, N2, O2, CO"
  do s = 1, num_species
     write(*,*) "species",s,":"
     read(*,*) species
     write(*,*) "mol fraction:"
     read(*,*) mol_fraction

     call get_species_constants( species, mass, dof, T_up )

     mol_mass = mol_mass + mol_fraction * mass
     cp_mix   = mol_fraction * n_up * boltz * ( dof * one_half + one )

     write(file_unit,fmt='(a2,a2,f8.6)') species, ": ", mol_fraction 
     write(file_unit,fmt='(a11,f13.10,a9,f7.5)') "    mass = ", mass/m_ref, " d.o.f = ", dof

  end do

  ! Currently using an infinite simple harmonic oscillator
  gamma  = cp_mix / ( cp_mix - n_up * boltz )

  write(*,*) ""
  write(*,*) "Molecular mass = ", mol_mass, "[kg]"
  write(*,*) "Specific heat ratio = ", gamma
  write(*,*) ""

  write(file_unit,fmt='(x)')
  write(file_unit,fmt='(a17,e16.8,a5)') "Molecular mass = ", mol_mass, " [kg]"
  write(file_unit,fmt='(a22,f7.5)') "Specific heat ratio = ", gamma
  write(file_unit,fmt='(x)')

  a_up = sqrt( gamma * T_up * boltz / mol_mass )

  write(*,*) "Upstream speed of sound = ", a_up, "[m/s]"
  write(*,*) ""

  ! Shock Ratios
  n_ratio = ( ( gamma + one ) * M_sh * M_sh ) / ( two + ( gamma - one ) * M_sh * M_sh )
  T_ratio = ( two + ( gamma - one ) * M_sh * M_sh ) * ( two * gamma * M_sh * M_sh - gamma + one ) &
       / ( ( gamma + one ) * ( gamma + one ) * M_sh * M_sh )
  p_ratio = one + ( two * gamma * ( M_sh * M_sh - one ) ) / ( gamma + one )

  ! Downstream speed quantities
  M_down = sqrt( ( one + one_half * ( gamma - one ) * M_sh * M_sh ) &
       / ( gamma * M_sh * M_sh - one_half * ( gamma - one ) ) )
  a_down = sqrt( gamma * T_ratio * T_up * boltz / mol_mass )
  u_down = M_down * a_down! - abs(u_shock)

  write(*,*) "Is the shock stationary? 1 = yes, 0 = no"
  read(*,*) stationary
  write(*,*) ""

  if( stationary == 0 )then
     u_shock = M_down * a_down
     !( gamma - ( gamma - one ) * M_sh * M_sh * a_up * a_up ) / ( ( gamma + one ) * M_sh * a_up )
     u_up    = M_sh * a_up - abs(u_shock)
  else
     u_shock = zero
     u_up    = M_sh * a_up
  end if

  M_up = u_up / a_up

  n_down = n_ratio * n_up
  T_down = T_ratio * T_up

  write(*,*) "Upstream Mach number (stationary reference frame) = ", M_up
  write(*,*) "Shock speed = ", u_shock, "[m/s]"
  write(*,*) "Upstream speed = ", u_up, "[m/s]"
  write(*,*) ""
  write(*,*) "Downstream density = ", n_down, "[#/m^3]"
  write(*,*) "Downstream temperature = ", T_down, "[K]"
  write(*,*) "Downstream Mach number (shock reference frame) = ", M_down
  write(*,*) "Downstream speed of sound = ", a_down, "[m/s]"
  write(*,*) "Downstream velocity = ", u_down, "[m/s]"
  write(*,*) ""  

  ! Scaled solutions
  write(*,*) "Reference speed = ", eta_ref, "[m/s]"
  write(*,*) ""
  write(*,*) "Scaled upstream density = ", n_up / n_ref
  write(*,*) "Scaled upstream temperature = ", T_up / T_ref
  write(*,*) "Scaled upstream speed = ", u_up / eta_ref
  write(*,*) ""
  write(*,*) "Scaled downstream density = ", n_down / n_ref
  write(*,*) "Scaled downstream temperature = ", T_down / T_ref
  write(*,*) "Scaled downstream speed = ", u_down / eta_ref

  if( gamma < 1.6 )then
     y_min = -3.5 * sqrt( T_down / T_ref )
     y_max = 3.5 * sqrt( T_down / T_ref )
     z_min = y_min
     z_max = y_max
     
     test  = -3.5 * sqrt( T_down / T_ref ) - u_down / eta_ref
     x_min = -3.5 * sqrt( T_up / T_ref ) - u_up / eta_ref
     if( test < x_min ) x_min = test

     test  = 3.5 * sqrt( T_down / T_ref ) - u_down / eta_ref
     x_max = 3.5 * sqrt( T_up / T_ref ) - u_up / eta_ref
     if( test > x_max ) x_max = test

  else
     y_min = -2.8 * sqrt( T_down / T_ref )
     y_max = 2.8 * sqrt( T_down / T_ref )
     z_min = y_min
     z_max = y_max
     
     test  = -2.8 * sqrt( T_down / T_ref ) - u_down / eta_ref
     x_min = -2.8 * sqrt( T_up / T_ref ) - u_up / eta_ref
     if( test < x_min ) x_min = test

     test  = 2.8 * sqrt( T_down / T_ref ) - u_down / eta_ref
     x_max = 2.8 * sqrt( T_up / T_ref ) - u_up / eta_ref
     if( test > x_max ) x_max = test

  end if

  ! Write to file
  write(file_unit,fmt='(a38,f8.6)') "Mach number (shock reference frame) = ", M_sh
  write(file_unit,fmt='(x)')
  write(file_unit,fmt='(a38,i1)') "Stationary shock ( yes = 1, no = 0 )? ", stationary 
  write(file_unit,fmt='(x)')
  write(file_unit,fmt='(a25,f13.10)') "Pressure ratio (P2/P1) = ", p_ratio
  write(file_unit,fmt='(a24,f13.10)') "Density ratio (n2/n1) = ", n_ratio
  write(file_unit,fmt='(a28,f13.10)') "Temperature ratio (T2/T1) = ", T_ratio
  write(file_unit,fmt='(x)')
  write(file_unit,fmt='(a7,f13.10)') "n_up = ", n_up / n_ref
  write(file_unit,fmt='(a7,f13.10)') "T_up = ", T_up / T_ref
  write(file_unit,fmt='(a7,f13.10)') "u_up = ", u_up / eta_ref
  write(file_unit,fmt='(a9,f13.10)') "n_down = ", n_down / n_ref
  write(file_unit,fmt='(a9,f13.10)') "T_down = ", T_down / T_ref
  write(file_unit,fmt='(a9,f13.10)') "u_down = ", u_down / eta_ref
  write(file_unit,fmt='(x)')
  write(file_unit,fmt='(a15,f14.10,a2,f14.10,a2)') "x-velocity = ( ", x_min, ", ", x_max, " )"
  write(file_unit,fmt='(a15,f14.10,a2,f14.10,a2)') "y-velocity = ( ", y_min, ", ", y_max, " )"
  write(file_unit,fmt='(a15,f14.10,a2,f14.10,a2)') "z-velocity = ( ", z_min, ", ", z_max, " )"
!!$  write(file_unit,fmt='(x)')
!!$  write(file_unit,fmt='()') 
!!$  write(file_unit,fmt='()') 
!!$  write(file_unit,fmt='()') 
!!$  write(file_unit,fmt='()') 

  close( file_unit )

end program normal_shock

subroutine get_species_constants( species, mass, dof, temp )

  implicit none

  double precision, parameter :: two = 2.0d0
  double precision, parameter :: one = 1.0d0

  character(len=2), intent(in) :: species
  double precision, intent(in) :: temp
  double precision :: mass, dof

  select case( species )
  case( 'He' )
     mass = 6.65e-27
     dof  = 3.0d0
  case( 'Ne' )
     mass = 33.5e-27
     dof  = 3.0d0
  case( 'Ar' )
     mass = 66.3e-27
     dof  = 3.0d0
  case( 'Xe' )
     mass = 218.e-27
     dof  = 3.0d0
  case( 'H2' )
     mass = 3.34e-27
     dof  = 5.0d0 + ( two * 6159.0d0 / temp ) / ( exp( 6159.0d0 / temp ) - one )
  case( 'N2' )
     mass = 46.5e-27
     dof  = 5.0d0 + ( two * 3371.0d0 / temp ) / ( exp( 3371.0d0 / temp ) - one )
  case( 'O2' )
     mass = 53.1e-27
     dof  = 5.0d0 + ( two * 2256.0d0 / temp ) / ( exp( 2256.0d0 / temp ) - one )
  case( 'CO' )
     mass = 46.5e-27
     dof  = 5.0d0 + ( two * 3103.0d0 / temp ) / ( exp( 3103.0d0 / temp ) - one )
  case default
     write(*,*) "Error: species not an option. species = ",species
     stop
  end select

  return
end subroutine get_species_constants
