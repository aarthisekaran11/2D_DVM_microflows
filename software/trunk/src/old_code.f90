
  subroutine variance_reduction_kernel2( phi, phi_old, tot_colls, n, m, ns, coln_rms, &
       molecule, vel_grid, properties )

    use DistFunc
    use VelocityGrid
    use PhysicalProperties
    use ReplenishingCollisions
    use TimeStepping
    use Scaling
    use RandomNumberGeneration

    implicit none

    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    type(DistFuncType), dimension(:), intent(in) :: phi_old
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(PropertiesType), intent(in) :: properties
    double precision, dimension(:,:), intent(in) :: coln_rms
    integer, intent(in) :: n, m, ns

    type(DistFuncType), dimension(:) :: phi
    integer :: tot_colls

    logical :: rot_flag, vib_flag

    double precision, dimension(:), allocatable :: frA, frB, fvA, fvB
    double precision, dimension(:), allocatable :: rotA, rotB, vibA, vibB

    double precision, dimension(:), allocatable :: rotA_pos, rotA_neg
    double precision, dimension(:), allocatable :: rotB_pos, rotB_neg
    double precision, dimension(:), allocatable :: vibA_pos, vibA_neg
    double precision, dimension(:), allocatable :: vibB_pos, vibB_neg

    double precision :: coeff
    double precision :: sumA, sumB
    double precision :: sum_rtA, sum_rtB, sum_vtA, sum_vtB
    double precision :: deplA, deplB, depl_adj

    double precision, dimension(5) :: depl_frac, depl_sign
    double precision, dimension(2) :: mass, diam, omega
    double precision, dimension(2) :: viscosity, temp_visc
    double precision, dimension(2) :: dens, eq_dens, neq_dens, temp, kin_neq_dens, kin_eq_dens
    double precision, dimension(2) :: depletion, kin_depl
    double precision :: m_red, omega_AB, vhs_exponent, sigma
    double precision :: dt, factor_coeff
    double precision :: factor!, depletion
    double precision :: g, g_sigma
    double precision :: Kn

    double precision :: densA, densB

    double precision :: phiA, phiB, phiA_eq, phiB_eq
    double precision :: phiA_kneq, phiB_kneq
    double precision :: phiA_rneq, phiB_rneq, phiA_vneq, phiB_vneq
    double precision :: signA, signB

    integer, dimension(4) :: rl, vl
    integer, dimension(2) :: r_modes, v_modes, r_levels, v_levels
    integer, dimension(2) :: i, j, k
    integer :: same
    integer :: coll, num_colls
    integer :: glA, glB
    integer :: status

    integer :: num_pointsA, num_pointsB

    num_pointsA = vel_grid(n)%num_points
    num_pointsB = vel_grid(m)%num_points

    ! Initialize flags for energy activation
    rot_flag = .false.
    vib_flag = .false.

    tot_colls = 0

    ! Species specific constants
    mass(1)  = molecule(n)%mass
    diam(1)  = molecule(n)%diam
    omega(1) = molecule(n)%omega
    viscosity(1) = molecule(n)%viscosity
    temp_visc(1) = molecule(n)%T_viscosity

    mass(2)  = molecule(m)%mass
    diam(2)  = molecule(m)%diam
    omega(2) = molecule(m)%omega
    viscosity(2) = molecule(m)%viscosity
    temp_visc(2) = molecule(m)%T_viscosity

    ! Species combination constants
    m_red = mass(1)*mass(2)/( mass(1) + mass(2) )

    omega_AB = one_half*( omega(1) + omega(2) )
    call calculate_cross_section( sigma, m_red, omega_AB, diam, viscosity, temp_visc )

    vhs_exponent = two - two*omega_AB

    ! Knudsen number and time step
    call get_knudsen_number( Kn )
    call get_deltat( dt )

    ! Factor coefficient accounts for double counting if both 
    ! collision partners are drawn from the same species
    if( n .eq. m )then
       factor_coeff = dt * one_half / Kn
       same = AA
    else
       factor_coeff = dt / Kn
       same = AB
    end if

    ! Equilibrium species properties
    dens(1) = properties%dens(n)
    temp(1) = properties%species_temp(n)
!!$
    dens(2) = properties%dens(m)
    temp(2) = properties%species_temp(m)

    ! Internal Structure
    r_modes(1) = molecule(n)%rot_modes
    v_modes(1) = molecule(n)%vib_modes

    r_modes(2) = molecule(m)%rot_modes
    v_modes(2) = molecule(m)%vib_modes

    ! Nonequilibrium density
    eq_dens(1) = cdf_eq(ns,n,same)%eq_dens
    eq_dens(2) = cdf_eq(ns,m,same)%eq_dens
    neq_dens(1) = cdf_neq(n,same)%neq_dens
    neq_dens(2) = cdf_neq(m,same)%neq_dens

    kin_neq_dens(1) = sum( abs( cdf_neq(n,same)%kin_df ) )
    kin_neq_dens(2) = sum( abs( cdf_neq(m,same)%kin_df ) )
    kin_eq_dens(1)  = properties%dens(n)
    kin_eq_dens(2)  = properties%dens(m)

!!$    write(*,*)"Non-equilibrium density = ", neq_dens
!!$
!!$    write(*,*)"total deviation: ",sum( cdf_neq(n,same)%kin_df(:) )

    ! Set levels and allocate energy arrays
    if( r_modes(1) .gt. 0 )then
       rot_flag = .true.
       r_levels(1) = phi(n)%num_rot_levels
       allocate( frA( 1:r_levels(1) ), rotA( 1:r_levels(1) ), STAT=status )
       call allocate_error_check( status, "frA, rotA (VR)" )
       allocate( rotA_pos( 1:r_levels(1) ), rotA_neg( 1:r_levels(1) ), STAT=status )
    end if

    if( r_modes(2) .gt. 0 )then
       rot_flag = .true.
       r_levels(2) = phi(m)%num_rot_levels
       allocate( frB( 1:r_levels(2) ), rotB( 1:r_levels(2) ), STAT=status )
       call allocate_error_check( status, "frB, rotB (VR)" )
       allocate( rotB_pos( 1:r_levels(2) ), rotB_neg( 1:r_levels(2) ), STAT=status )
    end if

    if( v_modes(1) .gt. 0 )then
       vib_flag = .true.
       v_levels(1) = phi(n)%num_vib_levels
       allocate( fvA( 1:v_levels(1) ), vibA( 1:v_levels(1) ), STAT=status )
       call allocate_error_check( status, "fvA, vibA (VR)" )
       allocate( vibA_pos( 1:v_levels(1) ), vibA_neg( 1:v_levels(1) ), STAT=status )
    end if

    if( v_modes(2) .gt. 0 )then
       vib_flag = .true.
       v_levels(2) = phi(m)%num_vib_levels
       allocate( fvB( 1:v_levels(2) ), vibB( 1:v_levels(2) ), STAT=status )
       call allocate_error_check( status, "fvB, vibB (VR)" )
       allocate( vibB_pos( 1:v_levels(2) ), vibB_neg( 1:v_levels(2) ), STAT=status )
    end if

    ! Separate depletion into elastic and inelastic parts (elastic, r-t A, r-t B, v-t A, v-t B)
    call split_depletion( depl_frac, m_red, properties, molecule, n, m )

    !======================================================================================================
    ! First collision integral: phi_eq*phi_noneq
    !                           A-B, not B-A
    !======================================================================================================
    ! Calculate the number of collisions
    ! TODO: the maximum number of collisions is bounded by 2(n^2)
    densA = eq_dens(1)/dens(1)
    densB = neq_dens(2)/dens(2)
    call compute_num_coll_pairs( num_colls, dt, Kn, densA, densB, &
         temp, mass, coln_rms, vel_grid, n, m )

    ! Number of collisions can't be zero
    if( num_colls .gt. 0 )then
       ! Calculate factor for depletion
       factor = factor_coeff * sigma * eq_dens(1) * neq_dens(2) / dble(num_colls)
    else
       factor = zero
    end if

    tot_colls = tot_colls + num_colls

    do coll = 1, num_colls
       ! Pick velocities to collide
       call pick_collision_partner( i(1), j(1), k(1), glA, cdf_eq(ns,n,same)%cumul_df, &
            cdf_eq(ns,n,same)%search_ref, vel_grid(n) )
       call pick_collision_partner( i(2), j(2), k(2), glB, cdf_neq(m,same)%cumul_df, &
            cdf_neq(m,same)%search_ref, vel_grid(m) )

       ! Calculate relative velocity between chosen collision partners
       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )

       ! Check for self collision
       if( g .lt. double_tol ) cycle

       ! Calculate g*sigma term
       g_sigma = g**vhs_exponent

       ! Distribution function values needed for calculations:
       phiA = phi_old(n)%value( i(1), j(1), k(1) ) 
       phiB = phi_old(m)%value( i(2), j(2), k(2) ) 

       phiB_eq   = cdf_eq(ns,m,same)%cumul_df( glB ) - cdf_eq(ns,m,same)%cumul_df( glB-1 )

       phiB_kneq = cdf_neq(m,same)%kin_df( glB )
       phiB_rneq = cdf_neq(m,same)%rot_df( glB )
       phiB_vneq = cdf_neq(m,same)%vib_df( glB )

       signB     = cdf_neq(m,same)%sign( glB )

       ! Equilibrium collision partner
       sumA = one

       if( r_modes(1) .gt. 0 )then
          frA = cdf_eq(ns,n,same)%rot_eq(:)

          sumA    = sum( frA )
          sum_vtA = sumA

          call pick_level( rl(1), frA, r_levels(1) )
          depl_sign(2) = one
       end if

       if( v_modes(1) .gt. 0 )then
          fvA = cdf_eq(ns,n,same)%vib_eq(:)

          sumA    = sum( fvA )
          sum_rtA = sumA

          call pick_level( vl(1), fvA, v_levels(1) )
          depl_sign(4) = one
       end if

       ! New Idea:
       ! Deviation collision partner:
       !    frB = P(rd,vE)*sum(fvE)*frd + P(rE,vd)*sum(fvd)*frE + P(rd,vd)*sum(fvd)*frd
       !    fvB = P(rd,vE)*sum(frd)*fvE + P(rE,vd)*sum(frE)*fvd + P(rd,vd)*sum(frd)*fvd
       !    rlevel - pick from |frB|, sign[frB(rlevel)]
       !    vlevel - pick from |fvB|, sign[fvB(vlevel)]
       !
       !    Elastic: dPhi*frB and dPhi*fvB, dPhi_B = dPhi*sum(frB) = dPhi*sum(fvB)
       !    R-T:     dPhi_B = dPhi*( P(rd,vE)*sum(fvE) + P(rE,vd)*sum(fvd) + P(rd,vd)*sum(fvd) )
       !    V-T:     dPhi_B = dPhi*( P(rd,vE)*sum(frd) + P(rE,vd)*sum(frE) + P(rd,vd)*sum(frd) )
       !
       ! Colliding partner:
       !    dPhi_A = sum(frE)*[dPhi*sum(frB)] = sum(fvE)*[dPhi*sum(fvB)] = 1*dPhi_B
       !


       coeff = cdf_neq(m,same)%coeff( glB )

       if( r_modes(2) .gt. 0 .and. v_modes(2) .gt. 0 )then
          frB = ( phi_old(m)%rot(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%rot_eq * phiB_eq ) / phiB_rneq
          fvB = ( phi_old(m)%vib(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%vib_eq * phiB_eq ) / phiB_vneq
          
          ! Elastic array - used for when particle does not exchange energy
          frB  = &
               ( phiB_eq * phiB_rneq * coeff )   * frB + &
               ( phiB_eq * phiB_vneq * coeff )   * ( phiB_kneq / phiB_vneq ) * cdf_eq(ns,m,same)%rot_eq + &
               ( phiB_vneq * phiB_rneq * coeff ) * abs( phiB_kneq / phiB_vneq ) * frB
          fvB  = &
               ( phiB_eq * phiB_rneq * coeff )   * ( phiB_kneq / phiB_rneq ) * cdf_eq(ns,m,same)%vib_eq + &
               ( phiB_eq * phiB_vneq * coeff )   * fvB + &
               ( phiB_vneq * phiB_rneq * coeff ) * abs( phiB_kneq / phiB_rneq ) * fvB

          call pick_level( rl(2), frB, r_levels(2) )
          call pick_level( vl(2), fvB, v_levels(2) )

          depl_sign(3) = dsgn( frB( rl(2) ) )
          depl_sign(5) = dsgn( fvB( vl(2) ) )

          ! Inelastic array - used for when particle undergoes r-t or v-t exchange
          rotB = &
               ( phiB_eq * phiB_vneq * coeff )       * depl_sign(5) * cdf_eq(ns,m,same)%rot_eq + &
               ( one - phiB_eq * phiB_vneq * coeff ) * frB
          vibB = &
               ( phiB_eq * phiB_rneq * coeff )       * depl_sign(3) * cdf_eq(ns,m,same)%vib_eq + &
               ( one - phiB_eq * phiB_rneq * coeff ) * fvB

          sumB    = sum( fvB )
          sum_vtB = sum( rotB )
          sum_rtB = sum( vibB )

       else if( r_modes(2) .gt. 0 )then
          frB = ( phi_old(m)%rot(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%rot_eq * phiB_eq ) / phiB_rneq

          sumB    = sum( frB )
          sum_vtB = sumB

          call pick_level( rl(2), frB, r_levels(2) )

          depl_sign(3) = dsgn( frB( rl(2) ) )

       else if( v_modes(2) .gt. 0 )then
          fvB = ( phi_old(m)%vib(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%vib_eq * phiB_eq ) / phiB_vneq

          sumB    = sum( fvB )
          sum_rtB = sumB

          call pick_level( vl(2), fvB, v_levels(2) )

          depl_sign(5) = dsgn( fvB( vl(2) ) )

       else
          sumB = signB

       end if

       !**ELASTIC*********************************************************************************************
       deplA = sumB
       deplB = sumA

       depletion(1) = deplA        * depl_frac(1) * factor * g_sigma
       depletion(2) = deplB        * depl_frac(1) * factor * g_sigma
       kin_depl(1)  = deplA * sumA * depl_frac(1) * factor * g_sigma
       kin_depl(2)  = deplB * sumB * depl_frac(1) * factor * g_sigma

       call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
            i(1), j(1), k(1), 0, frA, fvA, r_modes(1), v_modes(1) )
       call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
            i(2), j(2), k(2), 0, frB, fvB, r_modes(2), v_modes(2) )
       call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
            rl, i, j, k, g, n, m, vel_grid, m_red, elastic, A, molecule )

       !**ROTATIONAL - TRANSLATIONAL**************************************************************************
       if( rot_flag )then ! TODO: this call to rot_flag seems excessive....
          if( r_modes(1) .gt. 0 )then
             deplA = sumB
             deplB = sum_rtA

             depletion(1) = deplA           * depl_frac(2) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(2) * factor * g_sigma
             kin_depl(1)  = deplA * sum_rtA * depl_frac(2) * factor * g_sigma
             kin_depl(2)  = deplB * sumB    * depl_frac(2) * factor * g_sigma

             call vr_deplete( rot_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, vibA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, vibA, fvB, &
                  rl, i, j, k, g, n, m, vel_grid, m_red, rot_trans, A, molecule )

          end if

          if( r_modes(2) .gt. 0 )then
             deplA = sum_rtB
             deplB = sumA

             depletion(1) = deplA           * depl_frac(3) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(3) * factor * g_sigma
             kin_depl(1)  = deplA * sumA    * depl_frac(3) * factor * g_sigma
             kin_depl(2)  = deplB * sum_rtB * depl_frac(3) * factor * g_sigma

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( rot_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, vibB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, vibB, rl, i, j, k, &
                  g, n, m, vel_grid, m_red, rot_trans, B, molecule )

          end if
       end if

       !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
       if( vib_flag )then
          if( v_modes(1) .gt. 0 )then
             deplA = sumB
             deplB = sum_vtA

             depletion(1) = deplA           * depl_frac(4) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(4) * factor * g_sigma
             kin_depl(1)  = deplA * sum_vtA * depl_frac(4) * factor * g_sigma
             kin_depl(2)  = deplB * sumB    * depl_frac(4) * factor * g_sigma

             call vr_deplete( vib_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), rotA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, rotA, frB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, A, molecule )

          end if

          if( v_modes(2) .gt. 0 )then             
             deplA = sum_vtB
             deplB = sumA

             depletion(1) = deplA           * depl_frac(5) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(5) * factor * g_sigma
             kin_depl(1)  = deplA * sumA    * depl_frac(5) * factor * g_sigma
             kin_depl(2)  = deplB * sum_vtB * depl_frac(5) * factor * g_sigma

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( vib_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), rotB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, rotB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, B, molecule )

          end if
       end if

    end do


    !======================================================================================================
    ! Second collision integral: phi_noneq*phi_eq
    !                            B-A, not A-B
    !======================================================================================================
    ! Calculate the number of collisions
    densA = neq_dens(1)/dens(1)
    densB = eq_dens(2)/dens(2)
    call compute_num_coll_pairs( num_colls, dt, Kn, densA, densB, &
         temp, mass, coln_rms, vel_grid, n, m )
    
    ! Number of collisions can't be zero
    if( num_colls .gt. 0 )then
       ! Calculate factor for depletion
       factor = factor_coeff * sigma * neq_dens(1) * eq_dens(2) / dble(num_colls)
    else
       factor = zero
    end if

    tot_colls = tot_colls + num_colls

    do coll = 1, num_colls
       ! Pick velocities to collide
       call pick_collision_partner( i(1), j(1), k(1), glA, cdf_neq(n,same)%cumul_df, &
            cdf_neq(n,same)%search_ref, vel_grid(n) )
       call pick_collision_partner( i(2), j(2), k(2), glB, cdf_eq(ns,m,same)%cumul_df, &
            cdf_eq(ns,m,same)%search_ref, vel_grid(m) )

       ! Calculate relative velocity between chosen collision partners
       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )

       ! Check for self collision
       if( g .lt. double_tol ) cycle

       ! Calculate g*sigma term
       g_sigma = g**vhs_exponent

       ! Temporary arrays for phi
       phiA = abs( phi_old(n)%value( i(1), j(1), k(1) ) )
       phiB = abs( phi_old(m)%value( i(2), j(2), k(2) ) )

       phiA_eq  = cdf_eq(ns,n,same)%cumul_df( glA ) - cdf_eq(ns,n,same)%cumul_df( glA-1 )

       phiA_kneq  = cdf_neq(n,same)%kin_df( glA )
       phiA_rneq = cdf_neq(n,same)%rot_df( glA )
       phiA_vneq = cdf_neq(n,same)%vib_df( glA )

       signA = cdf_neq(n,same)%sign( glA )

       depl_sign = one
       sumA = one
       sumB = one
       sum_rtA = one
       sum_vtA = one
       sum_rtB = one
       sum_vtB = one


       ! Deviation collision partner
       coeff = cdf_neq(n,same)%coeff( glA )

       if( r_modes(1) .gt. 0 .and. v_modes(1) .gt. 0 )then
          frA = ( phi_old(m)%rot(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%rot_eq * phiA_eq ) / phiA_rneq
          fvA = ( phi_old(m)%vib(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%vib_eq * phiA_eq ) / phiA_vneq
          
          ! Elastic array - used for when particle does not exchange energy
          frA  = &
               ( phiA_eq * phiA_rneq * coeff )   * frA + &
               ( phiA_eq * phiA_vneq * coeff )   * ( phiA_kneq / phiA_vneq ) * cdf_eq(ns,n,same)%rot_eq + &
               ( phiA_vneq * phiA_rneq * coeff ) * abs( phiA_kneq / phiA_vneq ) * frA
          fvA  = &
               ( phiA_eq * phiA_rneq * coeff )   * ( phiA_kneq / phiA_rneq ) * cdf_eq(ns,n,same)%vib_eq + &
               ( phiA_eq * phiA_vneq * coeff )   * fvA + &
               ( phiA_vneq * phiA_rneq * coeff ) * abs( phiA_kneq / phiA_rneq ) * fvA

          call pick_level( rl(1), frA, r_levels(1) )
          call pick_level( vl(1), fvA, v_levels(1) )

          depl_sign(2) = dsgn( frA( rl(1) ) )
          depl_sign(4) = dsgn( fvA( vl(1) ) )

          ! Inelastic array - used for when particle undergoes r-t or v-t exchange
          rotA = &
               ( phiA_eq * phiA_vneq * coeff )       * depl_sign(5) * cdf_eq(ns,n,same)%rot_eq + &
               ( one - phiA_eq * phiA_vneq * coeff ) * frA
          vibA = &
               ( phiA_eq * phiA_rneq * coeff )       * depl_sign(3) * cdf_eq(ns,n,same)%vib_eq + &
               ( one - phiA_eq * phiA_rneq * coeff ) * fvA

          sumA    = sum( fvA )
          sum_vtA = sum( rotA )
          sum_rtA = sum( vibA )

       else if( r_modes(1) .gt. 0 )then
          frA = ( phi_old(m)%rot(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%rot_eq * phiA_eq ) / phiA_rneq

          sumA    = sum( frA )
          sum_vtA = sumA

          call pick_level( rl(1), frA, r_levels(1) )

          depl_sign(2) = dsgn( frA( rl(1) ) )

       else if( v_modes(1) .gt. 0 )then
          fvA = ( phi_old(m)%vib(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%vib_eq * phiA_eq ) / phiA_vneq

          sumA    = sum( fvA )
          sum_rtA = sumA

          call pick_level( vl(1), fvA, v_levels(1) )

          depl_sign(4) = dsgn( fvA( vl(1) ) )

       else
          sumA = signA

       end if

       ! Equilibrium collision partner
       sumB = one

       if( r_modes(2) .gt. 0 )then
          frB = cdf_eq(ns,m,same)%rot_eq(:)

          sumB    = sum( frB )
          sum_vtB = sumB

          call pick_level( rl(2), frB, r_levels(2) )
          depl_sign(3) = one
       end if

       if( v_modes(2) .gt. 0 )then
          fvB = cdf_eq(ns,m,same)%vib_eq(:)

          sumB    = sum( fvB )
          sum_rtB = sumB

          call pick_level( vl(2), fvB, v_levels(2) )
          depl_sign(5) = one
       end if


       !**ELASTIC*********************************************************************************************
       deplA = sumB
       deplB = sumA

       depletion(1) = deplA        * depl_frac(1) * factor * g_sigma
       depletion(2) = deplB        * depl_frac(1) * factor * g_sigma
       kin_depl(1)  = deplA * sumA * depl_frac(1) * factor * g_sigma
       kin_depl(2)  = deplB * sumB * depl_frac(1) * factor * g_sigma

       call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
            i(1), j(1), k(1), 0, frA, fvA, r_modes(1), v_modes(1) )
       call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
            i(2), j(2), k(2), 0, frB, fvB, r_modes(2), v_modes(2) )
       call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
            rl, i, j, k, g, n, m, vel_grid, m_red, elastic, A, molecule )

       !**ROTATIONAL - TRANSLATIONAL**************************************************************************
       if( rot_flag )then ! TODO: this call to rot_flag seems excessive....
          if( r_modes(1) .gt. 0 )then
             deplA = sumB
             deplB = sum_rtA

             depletion(1) = deplA           * depl_frac(2) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(2) * factor * g_sigma
             kin_depl(1)  = deplA * sum_rtA * depl_frac(2) * factor * g_sigma
             kin_depl(2)  = deplB * sumB    * depl_frac(2) * factor * g_sigma

             call vr_deplete( rot_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, vibA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, vibA, fvB, &
                  rl, i, j, k, g, n, m, vel_grid, m_red, rot_trans, A, molecule )

          end if


          if( r_modes(2) .gt. 0 )then
             deplA = sum_rtB
             deplB = sumA

             depletion(1) = deplA           * depl_frac(3) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(3) * factor * g_sigma
             kin_depl(1)  = deplA * sumA    * depl_frac(3) * factor * g_sigma
             kin_depl(2)  = deplB * sum_rtB * depl_frac(3) * factor * g_sigma

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( rot_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, vibB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, vibB, rl, i, j, k, &
                  g, n, m, vel_grid, m_red, rot_trans, B, molecule )

          end if
       end if

       !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
       if( vib_flag )then
          if( v_modes(1) .gt. 0 )then
             deplA = sumB
             deplB = sum_vtA

             depletion(1) = deplA           * depl_frac(4) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(4) * factor * g_sigma
             kin_depl(1)  = deplA * sum_vtA * depl_frac(4) * factor * g_sigma
             kin_depl(2)  = deplB * sumB    * depl_frac(4) * factor * g_sigma

             call vr_deplete( vib_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), rotA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, rotA, frB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, A, molecule )

          end if

          if( v_modes(2) .gt. 0 )then            
             deplA = sum_vtB
             deplB = sumA

             depletion(1) = deplA           * depl_frac(5) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(5) * factor * g_sigma
             kin_depl(1)  = deplA * sumA    * depl_frac(5) * factor * g_sigma
             kin_depl(2)  = deplB * sum_vtB * depl_frac(5) * factor * g_sigma

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( vib_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), rotB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, rotB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, B, molecule )

          end if
       end if

    end do


    !======================================================================================================
    ! Third collision integral: phi_noneq*phi_noneq
    !                           A-B and B-A
    !======================================================================================================
    ! Calculate the number of collisions
    densA = neq_dens(1)/dens(1)
    densB = neq_dens(2)/dens(2)
    call compute_num_coll_pairs( num_colls, dt, Kn, densA, densB, &
         temp, mass, coln_rms, vel_grid, n, m )

    ! Number of collisions can't be zero
    if( num_colls .gt. 0 )then
       ! Calculate factor for depletion
       factor = factor_coeff * sigma * neq_dens(1) * neq_dens(2) / dble(num_colls)
    else
       factor = zero
    end if
       
    tot_colls = tot_colls + num_colls

    do coll = 1, num_colls
       ! Pick velocities to collide
       call pick_collision_partner( i(1), j(1), k(1), glA, cdf_neq(n,same)%cumul_df, &
            cdf_neq(n,same)%search_ref, vel_grid(n) )
       call pick_collision_partner( i(2), j(2), k(2), glB, cdf_neq(m,same)%cumul_df, &
            cdf_neq(m,same)%search_ref, vel_grid(m) )

       ! Calculate relative velocity between chosen collision partners
       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )

       ! Check for self collision
       if( g .lt. double_tol ) cycle

       ! Calculate g*sigma term
       g_sigma = g**vhs_exponent

       ! Temporary arrays for phi
       phiA = abs( phi_old(n)%value( i(1), j(1), k(1) ) )
       phiB = abs( phi_old(m)%value( i(2), j(2), k(2) ) )

       phiA_eq  = cdf_eq(ns,n,same)%cumul_df( glA ) - cdf_eq(ns,n,same)%cumul_df( glA-1 )
       phiB_eq  = cdf_eq(ns,m,same)%cumul_df( glB ) - cdf_eq(ns,m,same)%cumul_df( glB-1 )

       phiA_kneq  = cdf_neq(n,same)%kin_df( glA )
       phiA_rneq = cdf_neq(n,same)%rot_df( glA )
       phiA_vneq = cdf_neq(n,same)%vib_df( glA )

       phiB_kneq = cdf_neq(m,same)%kin_df( glB )
       phiB_rneq = cdf_neq(m,same)%rot_df( glB )
       phiB_vneq = cdf_neq(m,same)%vib_df( glB )

       signA = cdf_neq(n,same)%sign( glA )
       signB = cdf_neq(m,same)%sign( glB )

       depl_sign = one
       sumA = one
       sumB = one
       sum_rtA = one
       sum_vtA = one
       sum_rtB = one
       sum_vtB = one
       
       ! Deviation collision partner - A
       coeff = cdf_neq(n,same)%coeff( glA )

       if( r_modes(1) .gt. 0 .and. v_modes(1) .gt. 0 )then
          frA = ( phi_old(m)%rot(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%rot_eq * phiA_eq ) / phiA_rneq
          fvA = ( phi_old(m)%vib(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%vib_eq * phiA_eq ) / phiA_vneq
          
          ! Elastic array - used for when particle does not exchange energy
          frA  = &
               ( phiA_eq * phiA_rneq * coeff )   * frA + &
               ( phiA_eq * phiA_vneq * coeff )   * ( phiA_kneq / phiA_vneq ) * cdf_eq(ns,n,same)%rot_eq + &
               ( phiA_vneq * phiA_rneq * coeff ) * abs( phiA_kneq / phiA_vneq ) * frA
          fvA  = &
               ( phiA_eq * phiA_rneq * coeff )   * ( phiA_kneq / phiA_rneq ) * cdf_eq(ns,n,same)%vib_eq + &
               ( phiA_eq * phiA_vneq * coeff )   * fvA + &
               ( phiA_vneq * phiA_rneq * coeff ) * abs( phiA_kneq / phiA_rneq ) * fvA

          call pick_level( rl(1), frA, r_levels(1) )
          call pick_level( vl(1), fvA, v_levels(1) )

          depl_sign(2) = dsgn( frA( rl(1) ) )
          depl_sign(4) = dsgn( fvA( vl(1) ) )

          ! Inelastic array - used for when particle undergoes r-t or v-t exchange
          rotA = &
               ( phiA_eq * phiA_vneq * coeff )       * depl_sign(5) * cdf_eq(ns,n,same)%rot_eq + &
               ( one - phiA_eq * phiA_vneq * coeff ) * frA
          vibA = &
               ( phiA_eq * phiA_rneq * coeff )       * depl_sign(3) * cdf_eq(ns,n,same)%vib_eq + &
               ( one - phiA_eq * phiA_rneq * coeff ) * fvA

          sumA    = sum( fvA )
          sum_vtA = sum( rotA )
          sum_rtA = sum( vibA )

       else if( r_modes(1) .gt. 0 )then
          frA = ( phi_old(m)%rot(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%rot_eq * phiA_eq ) / phiA_rneq

          sumA    = sum( frA )
          sum_vtA = sumA

          call pick_level( rl(1), frA, r_levels(1) )

          depl_sign(2) = dsgn( frA( rl(1) ) )

       else if( v_modes(1) .gt. 0 )then
          fvA = ( phi_old(m)%vib(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%vib_eq * phiA_eq ) / phiA_vneq

          sumA    = sum( fvA )
          sum_rtA = sumA

          call pick_level( vl(1), fvA, v_levels(1) )

          depl_sign(4) = dsgn( fvA( vl(1) ) )

       else
          sumA = signA

       end if

       ! Deviation collision partner - B
       coeff = cdf_neq(m,same)%coeff( glB )

       if( r_modes(2) .gt. 0 .and. v_modes(2) .gt. 0 )then
          frB = ( phi_old(m)%rot(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%rot_eq * phiB_eq ) / phiB_rneq
          fvB = ( phi_old(m)%vib(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%vib_eq * phiB_eq ) / phiB_vneq
          
          ! Elastic array - used for when particle does not exchange energy
          frB  = &
               ( phiB_eq * phiB_rneq * coeff )   * frB + &
               ( phiB_eq * phiB_vneq * coeff )   * ( phiB_kneq / phiB_vneq ) * cdf_eq(ns,m,same)%rot_eq + &
               ( phiB_vneq * phiB_rneq * coeff ) * abs( phiB_kneq / phiB_vneq ) * frB
          fvB  = &
               ( phiB_eq * phiB_rneq * coeff )   * ( phiB_kneq / phiB_rneq ) * cdf_eq(ns,m,same)%vib_eq + &
               ( phiB_eq * phiB_vneq * coeff )   * fvB + &
               ( phiB_vneq * phiB_rneq * coeff ) * abs( phiB_kneq / phiB_rneq ) * fvB

          call pick_level( rl(2), frB, r_levels(2) )
          call pick_level( vl(2), fvB, v_levels(2) )

          depl_sign(3) = dsgn( frB( rl(2) ) )
          depl_sign(5) = dsgn( fvB( vl(2) ) )

          ! Inelastic array - used for when particle undergoes r-t or v-t exchange
          rotB = &
               ( phiB_eq * phiB_vneq * coeff )       * depl_sign(5) * cdf_eq(ns,m,same)%rot_eq + &
               ( one - phiB_eq * phiB_vneq * coeff ) * frB
          vibB = &
               ( phiB_eq * phiB_rneq * coeff )       * depl_sign(3) * cdf_eq(ns,m,same)%vib_eq + &
               ( one - phiB_eq * phiB_rneq * coeff ) * fvB

          sumB    = sum( fvB )
          sum_vtB = sum( rotB )
          sum_rtB = sum( vibB )

       else if( r_modes(2) .gt. 0 )then
          frB = ( phi_old(m)%rot(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%rot_eq * phiB_eq ) / phiB_rneq

          sumB    = sum( frB )
          sum_vtB = sumB

          call pick_level( rl(2), frB, r_levels(2) )

          depl_sign(3) = dsgn( frB( rl(2) ) )

       else if( v_modes(2) .gt. 0 )then
          fvB = ( phi_old(m)%vib(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%vib_eq * phiB_eq ) / phiB_vneq

          sumB    = sum( fvB )
          sum_rtB = sumB

          call pick_level( vl(2), fvB, v_levels(2) )

          depl_sign(5) = dsgn( fvB( vl(2) ) )

       else
          sumB = signB

       end if

       depl_sign(1) = dsgn( sumA ) * dsgn( sumB )
       depl_adj     = sumA * sumB


       !**ELASTIC*********************************************************************************************
       deplA = sumB
       deplB = sumA

       depletion(1) = deplA        * depl_frac(1) * factor * g_sigma
       depletion(2) = deplB        * depl_frac(1) * factor * g_sigma
       kin_depl(1)  = deplA * sumA * depl_frac(1) * factor * g_sigma
       kin_depl(2)  = deplB * sumB * depl_frac(1) * factor * g_sigma

       call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
            i(1), j(1), k(1), 0, frA, fvA, r_modes(1), v_modes(1) )
       call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
            i(2), j(2), k(2), 0, frB, fvB, r_modes(2), v_modes(2) )
       call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
            rl, i, j, k, g, n, m, vel_grid, m_red, elastic, A, molecule )

       !**ROTATIONAL - TRANSLATIONAL**************************************************************************
       if( rot_flag )then ! TODO: this call to rot_flag seems excessive....
          if( r_modes(1) .gt. 0 )then
             deplA = sumB
             deplB = sum_rtA

             depletion(1) = deplA           * depl_frac(2) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(2) * factor * g_sigma
             kin_depl(1)  = deplA * sum_rtA * depl_frac(2) * factor * g_sigma
             kin_depl(2)  = deplB * sumB    * depl_frac(2) * factor * g_sigma

             call vr_deplete( rot_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, vibA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, vibA, fvB, &
                  rl, i, j, k, g, n, m, vel_grid, m_red, rot_trans, A, molecule )

          end if


          if( r_modes(2) .gt. 0 )then
             deplA = sum_rtB
             deplB = sumA

             depletion(1) = deplA           * depl_frac(3) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(3) * factor * g_sigma
             kin_depl(1)  = deplA * sumA    * depl_frac(3) * factor * g_sigma
             kin_depl(2)  = deplB * sum_rtB * depl_frac(3) * factor * g_sigma

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( rot_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, vibB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, vibB, rl, i, j, k, &
                  g, n, m, vel_grid, m_red, rot_trans, B, molecule )

          end if
       end if

       !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
       if( vib_flag )then
          if( v_modes(1) .gt. 0 )then
             deplA = sumB
             deplB = sum_vtA

             depletion(1) = deplA           * depl_frac(4) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(4) * factor * g_sigma
             kin_depl(1)  = deplA * sum_vtA * depl_frac(4) * factor * g_sigma
             kin_depl(2)  = deplB * sumB    * depl_frac(4) * factor * g_sigma

             call vr_deplete( vib_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), rotA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, rotA, frB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, A, molecule )

          end if

          if( v_modes(2) .gt. 0 )then            
             deplA = sum_vtB
             deplB = sumA

             depletion(1) = deplA           * depl_frac(5) * factor * g_sigma
             depletion(2) = deplB           * depl_frac(5) * factor * g_sigma
             kin_depl(1)  = deplA * sumA    * depl_frac(5) * factor * g_sigma
             kin_depl(2)  = deplB * sum_vtB * depl_frac(5) * factor * g_sigma

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( vib_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), rotB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, rotB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, B, molecule )

          end if
       end if

    end do


    ! deallocate energy arrays
    if( r_modes(1) .gt. 0 ) deallocate( frA, rotA )
    if( r_modes(2) .gt. 0 ) deallocate( frB, rotB )
    if( v_modes(1) .gt. 0 ) deallocate( fvA, vibA )
    if( v_modes(2) .gt. 0 ) deallocate( fvB, vibB )

    if( r_modes(1) .gt. 0 ) deallocate( rotA_pos, rotA_neg )
    if( r_modes(2) .gt. 0 ) deallocate( rotB_pos, rotB_neg )
    if( v_modes(1) .gt. 0 ) deallocate( vibA_pos, vibA_neg )
    if( v_modes(2) .gt. 0 ) deallocate( vibB_pos, vibB_neg )

    return
  end subroutine variance_reduction_kernel2


  subroutine vr_deplete( coll_type, phi, depletion, kin_depl, i, j, k, level, fr, fv, r_modes, v_modes )

    use DistFunc

    implicit none

    double precision, dimension(:), intent(in) :: fr, fv
    double precision, intent(in) :: depletion, kin_depl

    integer, intent(in) :: coll_type, level, i, j, k, r_modes, v_modes

    type(DistFuncType) :: phi

    select case( coll_type )
    case( vib_trans )
       if( v_modes .gt. 0 )then
          phi%vib( level, i, j, k ) = phi%vib( level, i, j, k ) - kin_depl
       end if
       if( r_modes .gt. 0 )then
          phi%rot( :, i, j, k ) = phi%rot( :, i, j, k ) - depletion*fr
       end if
       phi%value( i, j, k ) = phi%value( i, j, k ) - kin_depl

    case( rot_trans )
       if( r_modes .gt. 0 )then
          phi%rot( level, i, j, k ) = phi%rot( level, i, j, k ) - kin_depl
       end if
       if( v_modes .gt. 0 )then
          phi%vib( :, i, j, k ) = phi%vib( :, i, j, k ) - depletion*fv
       end if
       phi%value( i, j, k ) = phi%value( i, j, k ) - kin_depl

    case( elastic )
       if( v_modes .gt. 0 )then
          phi%vib( :, i, j, k ) = phi%vib( :, i, j, k ) - depletion*fv
       end if
       if( r_modes .gt. 0 )then
          phi%rot( :, i, j, k ) = phi%rot( :, i, j, k ) - depletion*fr
       end if
       phi%value( i, j, k ) = phi%value( i, j, k ) - kin_depl

    case default

    end select

    return
  end subroutine vr_deplete


subroutine variance_reduction_kernel( phi, phi_old, tot_colls, n, m, ns, coln_rms, &
       molecule, vel_grid, properties )

    use DistFunc
    use VelocityGrid
    use PhysicalProperties
    use ReplenishingCollisions
    use TimeStepping
    use Scaling
    use RandomNumberGeneration

    implicit none

    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    type(DistFuncType), dimension(:), intent(in) :: phi_old
    type(MoleculeType), dimension(:), intent(in) :: molecule
    type(PropertiesType), intent(in) :: properties
    double precision, dimension(:,:), intent(in) :: coln_rms
    integer, intent(in) :: n, m, ns

    type(DistFuncType), dimension(:) :: phi
    integer :: tot_colls

    logical :: rot_flag, vib_flag

    double precision, dimension(:), allocatable :: frA, frB, fvA, fvB
    double precision :: rotA, rotB, vibA, vibB

    double precision :: coeff
    double precision :: sumA, sumB

    double precision, dimension(5) :: depl_frac, depl_sign
    double precision, dimension(2) :: mass, diam, omega
    double precision, dimension(2) :: viscosity, temp_visc
    double precision, dimension(2) :: dens, eq_dens, neq_dens, temp, kin_neq_dens, kin_eq_dens, nc_dens_neq
    double precision, dimension(2) :: depletion, kin_depl
    double precision :: m_red, omega_AB, vhs_exponent, sigma
    double precision :: dt, factor_coeff
    double precision :: factor
    double precision :: g, g_sigma
    double precision :: Kn

    double precision :: densA, densB

    double precision :: phiA, phiB, phiA_eq, phiB_eq
    double precision :: phiA_kneq, phiB_kneq
    double precision :: phiA_rneq, phiB_rneq, phiA_vneq, phiB_vneq
    double precision :: signA, signB

    integer, dimension(4) :: rl, vl
    integer, dimension(2) :: r_modes, v_modes, r_levels, v_levels
    integer, dimension(2) :: i, j, k
    integer :: same
    integer :: coll, num_colls
    integer :: glA, glB
    integer :: status

    integer :: num_pointsA, num_pointsB

    num_pointsA = vel_grid(n)%num_points
    num_pointsB = vel_grid(m)%num_points

    ! Initialize flags for energy activation
    rot_flag = .false.
    vib_flag = .false.

    tot_colls = 0

    ! Species specific constants
    mass(1)  = molecule(n)%mass
    diam(1)  = molecule(n)%diam
    omega(1) = molecule(n)%omega
    viscosity(1) = molecule(n)%viscosity
    temp_visc(1) = molecule(n)%T_viscosity

    mass(2)  = molecule(m)%mass
    diam(2)  = molecule(m)%diam
    omega(2) = molecule(m)%omega
    viscosity(2) = molecule(m)%viscosity
    temp_visc(2) = molecule(m)%T_viscosity

    ! Species combination constants
    m_red = mass(1)*mass(2)/( mass(1) + mass(2) )

    omega_AB = one_half*( omega(1) + omega(2) )
    call calculate_cross_section( sigma, m_red, omega_AB, diam, viscosity, temp_visc )

    vhs_exponent = two - two*omega_AB

    ! Knudsen number and time step
    call get_knudsen_number( Kn )
    call get_deltat( dt )

    ! Factor coefficient accounts for double counting if both 
    ! collision partners are drawn from the same species
    if( n .eq. m )then
       factor_coeff = dt * one_half / Kn
       same = AA
    else
       factor_coeff = dt / Kn
       same = AB
    end if

    ! Equilibrium species properties
    dens(1) = properties%dens(n)
    temp(1) = properties%species_temp(n)

    dens(2) = properties%dens(m)
    temp(2) = properties%species_temp(m)

    ! Internal Structure
    r_modes(1) = molecule(n)%rot_modes
    v_modes(1) = molecule(n)%vib_modes

    r_modes(2) = molecule(m)%rot_modes
    v_modes(2) = molecule(m)%vib_modes

    ! Nonequilibrium density
    eq_dens(1) = cdf_eq(ns,n,same)%eq_dens
    eq_dens(2) = cdf_eq(ns,m,same)%eq_dens
    neq_dens(1) = cdf_neq(n,same)%neq_dens
    neq_dens(2) = cdf_neq(m,same)%neq_dens

    kin_neq_dens(1) = sum( abs( cdf_neq(n,same)%kin_df ) )
    kin_neq_dens(2) = sum( abs( cdf_neq(m,same)%kin_df ) )
    kin_eq_dens(1)  = properties%dens(n)
    kin_eq_dens(2)  = properties%dens(m)

    nc_dens_neq(1) = neq_dens(1)
    nc_dens_neq(2) = neq_dens(2)

    write(*,*)neq_dens
    write(*,*)kin_neq_dens

    ! Set levels and allocate energy arrays
    if( r_modes(1) .gt. 0 )then
       rot_flag = .true.
       r_levels(1) = phi(n)%num_rot_levels
       allocate( frA( 1:r_levels(1) ), STAT=status )
       call allocate_error_check( status, "frA(VR)" )
    end if

    if( r_modes(2) .gt. 0 )then
       rot_flag = .true.
       r_levels(2) = phi(m)%num_rot_levels
       allocate( frB( 1:r_levels(2) ), STAT=status )
       call allocate_error_check( status, "frB(VR)" )
    end if

    if( v_modes(1) .gt. 0 )then
       vib_flag = .true.
       v_levels(1) = phi(n)%num_vib_levels
       allocate( fvA( 1:v_levels(1) ), STAT=status )
       call allocate_error_check( status, "fvA(VR)" )
    end if

    if( v_modes(2) .gt. 0 )then
       vib_flag = .true.
       v_levels(2) = phi(m)%num_vib_levels
       allocate( fvB( 1:v_levels(2) ), STAT=status )
       call allocate_error_check( status, "fvB(VR)" )
    end if

    ! Separate depletion into elastic and inelastic parts (elastic, r-t A, r-t B, v-t A, v-t B)
    call split_depletion( depl_frac, m_red, properties, molecule, n, m )

    !======================================================================================================
    ! First collision integral: phi_eq*phi_noneq
    !                           A-B, not B-A
    !======================================================================================================
    ! Calculate the number of collisions
    ! TODO: the maximum number of collisions is bounded by 2(n^2)
    densA = eq_dens(1)/dens(1)
    densB = one_half*neq_dens(2)/dens(2)
    call compute_num_coll_pairs( num_colls, dt, Kn, densA, densB, &
         temp, mass, coln_rms, vel_grid, n, m )

    ! Number of collisions can't be zero
    if( num_colls .gt. 0 )then
       ! Calculate factor for depletion
       factor = factor_coeff * sigma * eq_dens(1) * neq_dens(2) / dble(num_colls)
    else
       factor = zero
    end if

    tot_colls = tot_colls + num_colls

    do coll = 1, num_colls
       ! Pick velocities to collide
       call pick_collision_partner( i(1), j(1), k(1), glA, cdf_eq(ns,n,same)%cumul_df, &
            cdf_eq(ns,n,same)%search_ref, vel_grid(n) )
       call pick_collision_partner( i(2), j(2), k(2), glB, cdf_neq(m,same)%cumul_df, &
            cdf_neq(m,same)%search_ref, vel_grid(m) )

       ! Calculate relative velocity between chosen collision partners
       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )

       ! Check for self collision
       if( g .lt. double_tol ) cycle

       ! Calculate g*sigma term
       g_sigma = g**vhs_exponent

       ! Distribution function values needed for calculations:
       phiA = phi_old(n)%value( i(1), j(1), k(1) ) 
       phiB = phi_old(m)%value( i(2), j(2), k(2) ) 

       phiB_eq   = cdf_eq(ns,m,same)%cumul_df( glB ) - cdf_eq(ns,m,same)%cumul_df( glB-1 )

       phiB_kneq = cdf_neq(m,same)%kin_df( glB )
       phiB_rneq = cdf_neq(m,same)%rot_df( glB )
       phiB_vneq = cdf_neq(m,same)%vib_df( glB )

       signB     = cdf_neq(m,same)%sign( glB )

       ! Equilibrium collision partner
       sumA = one

       if( r_modes(1) .gt. 0 )then
          frA = cdf_eq(ns,n,same)%rot_eq(:)

          sumA = one
          rotA = one

          call pick_level( rl(1), frA, r_levels(1) )
          depl_sign(2) = one
       end if

       if( v_modes(1) .gt. 0 )then
          fvA = cdf_eq(ns,n,same)%vib_eq(:)

          sumA = one
          vibA = one

          call pick_level( vl(1), fvA, v_levels(1) )
          depl_sign(4) = one
       end if

       ! New Idea:
       ! Deviation collision partner:
       !    frB = P(rd,vE)*sum(fvE)*frd + P(rE,vd)*sum(fvd)*frE + P(rd,vd)*sum(fvd)*frd
       !    fvB = P(rd,vE)*sum(frd)*fvE + P(rE,vd)*sum(frE)*fvd + P(rd,vd)*sum(frd)*fvd
       !    rlevel - pick from |frB|, sign[frB(rlevel)]
       !    vlevel - pick from |fvB|, sign[fvB(vlevel)]
       !
       !    Elastic: dPhi*frB and dPhi*fvB, dPhi_B = dPhi*sum(frB) = dPhi*sum(fvB)
       !    R-T:     dPhi_B = dPhi*( P(rd,vE)*sum(fvE) + P(rE,vd)*sum(fvd) + P(rd,vd)*sum(fvd) )
       !    V-T:     dPhi_B = dPhi*( P(rd,vE)*sum(frd) + P(rE,vd)*sum(frE) + P(rd,vd)*sum(frd) )
       !
       ! Colliding partner:
       !    dPhi_A = sum(frE)*[dPhi*sum(frB)] = sum(fvE)*[dPhi*sum(fvB)] = 1*dPhi_B
       !

       ! Deviation collision partner
       coeff = cdf_neq(m,same)%coeff( glB )

       if( r_modes(2) .gt. 0 .and. v_modes(2) .gt. 0 )then
          frB = ( phi_old(m)%rot(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%rot_eq * phiB_eq ) / phiB_rneq
          fvB = ( phi_old(m)%vib(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%vib_eq * phiB_eq ) / phiB_vneq
          
          ! Elastic array - used for when particle does not exchange energy
          frB  = &
               ( phiB * phiB_rneq * coeff )    * frB + &
               ( phiB_eq * phiB_kneq * coeff ) * cdf_eq(ns,m,same)%rot_eq
          fvB  = &
               ( phiB * phiB_vneq * coeff )    * fvB + &
               ( phiB_eq * phiB_kneq * coeff ) * cdf_eq(ns,m,same)%vib_eq

          call pick_level( rl(2), frB, r_levels(2) )
          call pick_level( vl(2), fvB, v_levels(2) )

          depl_sign(3) = dsgn( frB( rl(2) ) )
          depl_sign(5) = dsgn( fvB( vl(2) ) )

          ! Inelastic depletion - used for when particle undergoes r-t or v-t exchange
          rotB = ( phiB_eq * phiB_rneq ) * coeff + ( phiB_eq + phiB_rneq ) * phiB_kneq * coeff
          vibB = ( phiB_eq * phiB_vneq ) * coeff + ( phiB_eq + phiB_vneq ) * phiB_kneq * coeff

          sumB = sum( fvB )
          rotB = depl_sign(3) * rotB
          vibB = depl_sign(5) * vibB

       else if( r_modes(2) .gt. 0 )then
          frB = ( phi_old(m)%rot(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%rot_eq * phiB_eq ) / phiB_rneq

          call pick_level( rl(2), frB, r_levels(2) )

          depl_sign(3) = dsgn( frB( rl(2) ) )

          sumB = sum( frB )
          rotB = depl_sign(3)

       else if( v_modes(2) .gt. 0 )then
          fvB = ( phi_old(m)%vib(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%vib_eq * phiB_eq ) / phiB_vneq

          call pick_level( vl(2), fvB, v_levels(2) )

          depl_sign(5) = dsgn( fvB( vl(2) ) )

          sumB = sum( fvB )
          vibB = depl_sign(5)

       else
          sumB = signB

       end if

       ! TODO: repeated code
       !**ELASTIC*********************************************************************************************
       depletion(1) = sumB * depl_frac(1) * factor * g_sigma
       depletion(2) = sumA * depl_frac(1) * factor * g_sigma
       kin_depl(1)  = sumA * depletion(1)
       kin_depl(2)  = sumB * depletion(2)

       call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
            i(1), j(1), k(1), 0, frA, fvA, r_modes(1), v_modes(1) )
       call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
            i(2), j(2), k(2), 0, frB, fvB, r_modes(2), v_modes(2) )
       call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
            rl, i, j, k, g, n, m, vel_grid, m_red, elastic, A, molecule )

       !**ROTATIONAL - TRANSLATIONAL**************************************************************************
       if( rot_flag )then ! TODO: this call to rot_flag seems excessive....
          if( r_modes(1) .gt. 0 )then
             depletion(1) = sumB * depl_frac(2) * factor * g_sigma
             depletion(2) = rotA * depl_frac(2) * factor * g_sigma
             kin_depl(1)  = rotA * depletion(1)
             kin_depl(2)  = sumB * depletion(2)

             call vr_deplete( rot_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
                  rl, i, j, k, g, n, m, vel_grid, m_red, rot_trans, A, molecule )

          end if

          if( r_modes(2) .gt. 0 )then
             depletion(1) = rotB * depl_frac(3) * factor * g_sigma
             depletion(2) = sumA * depl_frac(3) * factor * g_sigma
             kin_depl(1)  = sumA * depletion(1)
             kin_depl(2)  = rotB * depletion(2)

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( rot_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, rl, i, j, k, &
                  g, n, m, vel_grid, m_red, rot_trans, B, molecule )

          end if
       end if

       !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
       if( vib_flag )then
          if( v_modes(1) .gt. 0 )then
             depletion(1) = sumB * depl_frac(4) * factor * g_sigma
             depletion(2) = vibA * depl_frac(4) * factor * g_sigma
             kin_depl(1)  = vibA * depletion(1)
             kin_depl(2)  = sumB * depletion(2)

             call vr_deplete( vib_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, A, molecule )

          end if

          if( v_modes(2) .gt. 0 )then
             depletion(1) = vibB * depl_frac(5) * factor * g_sigma
             depletion(2) = sumA * depl_frac(5) * factor * g_sigma
             kin_depl(1)  = sumA * depletion(1)
             kin_depl(2)  = vibB * depletion(2)

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( vib_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, B, molecule )

          end if
       end if

    end do


    !======================================================================================================
    ! Second collision integral: phi_noneq*phi_eq
    !                            B-A, not A-B
    !======================================================================================================
    ! Calculate the number of collisions
    densA = neq_dens(1)/dens(1)
    densB = one_half*eq_dens(2)/dens(2)
    call compute_num_coll_pairs( num_colls, dt, Kn, densA, densB, &
         temp, mass, coln_rms, vel_grid, n, m )
    
    ! Number of collisions can't be zero
    if( num_colls .gt. 0 )then
       ! Calculate factor for depletion
       factor = factor_coeff * sigma * neq_dens(1) * eq_dens(2) / dble(num_colls)
    else
       factor = zero
    end if

    tot_colls = tot_colls + num_colls

    do coll = 1, num_colls
       ! Pick velocities to collide
       call pick_collision_partner( i(1), j(1), k(1), glA, cdf_neq(n,same)%cumul_df, &
            cdf_neq(n,same)%search_ref, vel_grid(n) )
       call pick_collision_partner( i(2), j(2), k(2), glB, cdf_eq(ns,m,same)%cumul_df, &
            cdf_eq(ns,m,same)%search_ref, vel_grid(m) )

       ! Calculate relative velocity between chosen collision partners
       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )

       ! Check for self collision
       if( g .lt. double_tol ) cycle

       ! Calculate g*sigma term
       g_sigma = g**vhs_exponent

       ! Temporary arrays for phi
       phiA = abs( phi_old(n)%value( i(1), j(1), k(1) ) )
       phiB = abs( phi_old(m)%value( i(2), j(2), k(2) ) )

       phiA_eq  = cdf_eq(ns,n,same)%cumul_df( glA ) - cdf_eq(ns,n,same)%cumul_df( glA-1 )

       phiA_kneq = cdf_neq(n,same)%kin_df( glA )
       phiA_rneq = cdf_neq(n,same)%rot_df( glA )
       phiA_vneq = cdf_neq(n,same)%vib_df( glA )

       signA = cdf_neq(n,same)%sign( glA )

       depl_sign = one
       sumA = one
       sumB = one
       rotA = one
       rotB = one
       vibA = one
       vibB = one

       ! Deviation collision partner
       coeff = cdf_neq(n,same)%coeff( glA )

       if( r_modes(1) .gt. 0 .and. v_modes(1) .gt. 0 )then
          frA = ( phi_old(m)%rot(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%rot_eq * phiA_eq ) / phiA_rneq
          fvA = ( phi_old(m)%vib(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%vib_eq * phiA_eq ) / phiA_vneq
          
          ! Elastic array - used for when particle does not exchange energy
          frA  = &
               ( phiA * phiA_rneq * coeff )    * frA + &
               ( phiA_eq * phiA_kneq * coeff ) * cdf_eq(ns,n,same)%rot_eq
          fvA  = &
               ( phiA * phiA_vneq * coeff )    * fvA + &
               ( phiA_eq * phiA_kneq * coeff ) * cdf_eq(ns,n,same)%vib_eq

          call pick_level( rl(1), frA, r_levels(1) )
          call pick_level( vl(1), fvA, v_levels(1) )

          depl_sign(2) = dsgn( frA( rl(1) ) )
          depl_sign(4) = dsgn( fvA( vl(1) ) )

          ! Inelastic array - used for when particle undergoes r-t or v-t exchange
          rotA = ( phiA_eq * phiA_rneq ) * coeff + ( phiA_eq + phiA_rneq ) * phiA_kneq * coeff
          vibA = ( phiA_eq * phiA_vneq ) * coeff + ( phiA_eq + phiA_vneq ) * phiA_kneq * coeff

          sumA = sum( fvA )
          rotA = depl_sign(2) * rotA
          vibA = depl_sign(4) * vibA

       else if( r_modes(1) .gt. 0 )then
          frA = ( phi_old(m)%rot(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%rot_eq * phiA_eq ) / phiA_rneq

          call pick_level( rl(1), frA, r_levels(1) )

          depl_sign(2) = dsgn( frA( rl(1) ) )

          sumA = sum( frA )
          rotA = depl_sign(2)

       else if( v_modes(1) .gt. 0 )then
          fvA = ( phi_old(m)%vib(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%vib_eq * phiA_eq ) / phiA_vneq

          call pick_level( vl(1), fvA, v_levels(1) )

          depl_sign(4) = dsgn( fvA( vl(1) ) )

          sumA = sum( fvA )
          vibA = depl_sign(4)

       else
          sumA = signA

       end if

       ! Equilibrium collision partner
       sumB = one

       if( r_modes(2) .gt. 0 )then
          frB = cdf_eq(ns,m,same)%rot_eq(:)

          sumB = one
          rotB = one

          call pick_level( rl(2), frB, r_levels(2) )
          depl_sign(3) = one
       end if

       if( v_modes(2) .gt. 0 )then
          fvB = cdf_eq(ns,m,same)%vib_eq(:)

          sumB = one
          vibB = one

          call pick_level( vl(2), fvB, v_levels(2) )
          depl_sign(5) = one
       end if

       !**ELASTIC*********************************************************************************************
       depletion(1) = sumB * depl_frac(1) * factor * g_sigma
       depletion(2) = sumA * depl_frac(1) * factor * g_sigma
       kin_depl(1)  = sumA * depletion(1)
       kin_depl(2)  = sumB * depletion(2)

       call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
            i(1), j(1), k(1), 0, frA, fvA, r_modes(1), v_modes(1) )
       call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
            i(2), j(2), k(2), 0, frB, fvB, r_modes(2), v_modes(2) )
       call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
            rl, i, j, k, g, n, m, vel_grid, m_red, elastic, A, molecule )

       !**ROTATIONAL - TRANSLATIONAL**************************************************************************
       if( rot_flag )then ! TODO: this call to rot_flag seems excessive....
          if( r_modes(1) .gt. 0 )then
             depletion(1) = sumB * depl_frac(2) * factor * g_sigma
             depletion(2) = rotA * depl_frac(2) * factor * g_sigma
             kin_depl(1)  = rotA * depletion(1)
             kin_depl(2)  = sumB * depletion(2)

             call vr_deplete( rot_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
                  rl, i, j, k, g, n, m, vel_grid, m_red, rot_trans, A, molecule )

          end if

          if( r_modes(2) .gt. 0 )then
             depletion(1) = rotB * depl_frac(3) * factor * g_sigma
             depletion(2) = sumA * depl_frac(3) * factor * g_sigma
             kin_depl(1)  = sumA * depletion(1)
             kin_depl(2)  = rotB * depletion(2)

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( rot_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, rl, i, j, k, &
                  g, n, m, vel_grid, m_red, rot_trans, B, molecule )

          end if
       end if

       !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
       if( vib_flag )then
          if( v_modes(1) .gt. 0 )then
             depletion(1) = sumB * depl_frac(4) * factor * g_sigma
             depletion(2) = vibA * depl_frac(4) * factor * g_sigma
             kin_depl(1)  = vibA * depletion(1)
             kin_depl(2)  = sumB * depletion(2)

             call vr_deplete( vib_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, A, molecule )

          end if

          if( v_modes(2) .gt. 0 )then
             depletion(1) = vibB * depl_frac(5) * factor * g_sigma
             depletion(2) = sumA * depl_frac(5) * factor * g_sigma
             kin_depl(1)  = sumA * depletion(1)
             kin_depl(2)  = vibB * depletion(2)

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( vib_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, B, molecule )

          end if
       end if

    end do


    !======================================================================================================
    ! Third collision integral: phi_noneq*phi_noneq
    !                           A-B and B-A
    !======================================================================================================
    ! Calculate the number of collisions
    densA = one_half*neq_dens(1)/dens(1)
    densB = one_half*neq_dens(2)/dens(2)
    call compute_num_coll_pairs( num_colls, dt, Kn, densA, densB, &
         temp, mass, coln_rms, vel_grid, n, m )

    ! Number of collisions can't be zero
    if( num_colls .gt. 0 )then
       ! Calculate factor for depletion
       factor = factor_coeff * sigma * neq_dens(1) * neq_dens(2) / dble(num_colls)
    else
       factor = zero
    end if
       
    tot_colls = tot_colls + num_colls

    do coll = 1, num_colls
       ! Pick velocities to collide
       call pick_collision_partner( i(1), j(1), k(1), glA, cdf_neq(n,same)%cumul_df, &
            cdf_neq(n,same)%search_ref, vel_grid(n) )
       call pick_collision_partner( i(2), j(2), k(2), glB, cdf_neq(m,same)%cumul_df, &
            cdf_neq(m,same)%search_ref, vel_grid(m) )

       ! Calculate relative velocity between chosen collision partners
       call compute_relative_velocity( i, j, k, g, vel_grid, n, m )

       ! Check for self collision
       if( g .lt. double_tol ) cycle

       ! Calculate g*sigma term
       g_sigma = g**vhs_exponent

       ! Temporary arrays for phi
       phiA = abs( phi_old(n)%value( i(1), j(1), k(1) ) )
       phiB = abs( phi_old(m)%value( i(2), j(2), k(2) ) )

       phiA_eq  = cdf_eq(ns,n,same)%cumul_df( glA ) - cdf_eq(ns,n,same)%cumul_df( glA-1 )
       phiB_eq  = cdf_eq(ns,m,same)%cumul_df( glB ) - cdf_eq(ns,m,same)%cumul_df( glB-1 )

       phiA_kneq = cdf_neq(n,same)%kin_df( glA )
       phiA_rneq = cdf_neq(n,same)%rot_df( glA )
       phiA_vneq = cdf_neq(n,same)%vib_df( glA )

       phiB_kneq = cdf_neq(m,same)%kin_df( glB )
       phiB_rneq = cdf_neq(m,same)%rot_df( glB )
       phiB_vneq = cdf_neq(m,same)%vib_df( glB )

       signA = cdf_neq(n,same)%sign( glA )
       signB = cdf_neq(m,same)%sign( glB )
       
       depl_sign = one
       sumA = one
       sumB = one
       rotA = one
       rotB = one
       vibA = one
       vibB = one
       
        ! Deviation collision partner - A
       coeff = cdf_neq(n,same)%coeff( glA )

       if( r_modes(1) .gt. 0 .and. v_modes(1) .gt. 0 )then
          frA = ( phi_old(m)%rot(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%rot_eq * phiA_eq ) / phiA_rneq
          fvA = ( phi_old(m)%vib(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%vib_eq * phiA_eq ) / phiA_vneq
          
          ! Elastic array - used for when particle does not exchange energy
          frA  = &
               ( phiA * phiA_rneq * coeff )    * frA + &
               ( phiA_eq * phiA_kneq * coeff ) * cdf_eq(ns,n,same)%rot_eq
          fvA  = &
               ( phiA * phiA_vneq * coeff )    * fvA + &
               ( phiA_eq * phiA_kneq * coeff ) * cdf_eq(ns,n,same)%vib_eq

          call pick_level( rl(1), frA, r_levels(1) )
          call pick_level( vl(1), fvA, v_levels(1) )

          depl_sign(2) = dsgn( frA( rl(1) ) )
          depl_sign(4) = dsgn( fvA( vl(1) ) )

          ! Inelastic array - used for when particle undergoes r-t or v-t exchange
          rotA = ( phiA_eq * phiA_rneq ) * coeff + ( phiA_eq + phiA_rneq ) * phiA_kneq * coeff
          vibA = ( phiA_eq * phiA_vneq ) * coeff + ( phiA_eq + phiA_vneq ) * phiA_kneq * coeff

          sumA = sum( fvA )
          rotA = depl_sign(2) * rotA
          vibA = depl_sign(4) * vibA

       else if( r_modes(1) .gt. 0 )then
          frA = ( phi_old(m)%rot(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%rot_eq * phiA_eq ) / phiA_rneq

          call pick_level( rl(1), frA, r_levels(1) )

          depl_sign(2) = dsgn( frA( rl(1) ) )

          sumA = sum( frA )
          rotA = depl_sign(2)

       else if( v_modes(1) .gt. 0 )then
          fvA = ( phi_old(m)%vib(:,i(1),j(1),k(1)) - cdf_eq(ns,n,same)%vib_eq * phiA_eq ) / phiA_vneq

          call pick_level( vl(1), fvA, v_levels(1) )

          depl_sign(4) = dsgn( fvA( vl(1) ) )

          sumA = sum( fvA )
          vibA = depl_sign(4)

       else
          sumA = signA

       end if

       ! Deviation collision partner - B
       coeff = cdf_neq(m,same)%coeff( glB )

       if( r_modes(2) .gt. 0 .and. v_modes(2) .gt. 0 )then
          frB = ( phi_old(m)%rot(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%rot_eq * phiB_eq ) / phiB_rneq
          fvB = ( phi_old(m)%vib(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%vib_eq * phiB_eq ) / phiB_vneq
          
          ! Elastic array - used for when particle does not exchange energy
          frB  = &
               ( phiB * phiB_rneq * coeff )    * frB + &
               ( phiB_eq * phiB_kneq * coeff ) * cdf_eq(ns,m,same)%rot_eq
          fvB  = &
               ( phiB * phiB_vneq * coeff )    * fvB + &
               ( phiB_eq * phiB_kneq * coeff ) * cdf_eq(ns,m,same)%vib_eq

          call pick_level( rl(2), frB, r_levels(2) )
          call pick_level( vl(2), fvB, v_levels(2) )

          depl_sign(3) = dsgn( frB( rl(2) ) )
          depl_sign(5) = dsgn( fvB( vl(2) ) )

          ! Inelastic depletion - used for when particle undergoes r-t or v-t exchange
          rotB = ( phiB_eq * phiB_rneq ) * coeff + ( phiB_eq + phiB_rneq ) * phiB_kneq * coeff
          vibB = ( phiB_eq * phiB_vneq ) * coeff + ( phiB_eq + phiB_vneq ) * phiB_kneq * coeff

          sumB = sum( fvB )
          rotB = depl_sign(3) * rotB
          vibB = depl_sign(5) * vibB

       else if( r_modes(2) .gt. 0 )then
          frB = ( phi_old(m)%rot(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%rot_eq * phiB_eq ) / phiB_rneq

          call pick_level( rl(2), frB, r_levels(2) )

          depl_sign(3) = dsgn( frB( rl(2) ) )

          sumB = sum( frB )
          rotB = depl_sign(3)

       else if( v_modes(2) .gt. 0 )then
          fvB = ( phi_old(m)%vib(:,i(2),j(2),k(2)) - cdf_eq(ns,m,same)%vib_eq * phiB_eq ) / phiB_vneq

          call pick_level( vl(2), fvB, v_levels(2) )

          depl_sign(5) = dsgn( fvB( vl(2) ) )

          sumB = sum( fvB )
          vibB = depl_sign(5)

       else
          sumB = signB

       end if

       !**ELASTIC*********************************************************************************************
       depletion(1) = sumB * depl_frac(1) * factor * g_sigma
       depletion(2) = sumA * depl_frac(1) * factor * g_sigma
       kin_depl(1)  = sumA * depletion(1)
       kin_depl(2)  = sumB * depletion(2)

       call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
            i(1), j(1), k(1), 0, frA, fvA, r_modes(1), v_modes(1) )
       call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
            i(2), j(2), k(2), 0, frB, fvB, r_modes(2), v_modes(2) )
       call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
            rl, i, j, k, g, n, m, vel_grid, m_red, elastic, A, molecule )

       !**ROTATIONAL - TRANSLATIONAL**************************************************************************
       if( rot_flag )then ! TODO: this call to rot_flag seems excessive....
          if( r_modes(1) .gt. 0 )then
             depletion(1) = sumB * depl_frac(2) * factor * g_sigma
             depletion(2) = rotA * depl_frac(2) * factor * g_sigma
             kin_depl(1)  = rotA * depletion(1)
             kin_depl(2)  = sumB * depletion(2)

             call vr_deplete( rot_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
                  rl, i, j, k, g, n, m, vel_grid, m_red, rot_trans, A, molecule )

          end if

          if( r_modes(2) .gt. 0 )then
             depletion(1) = rotB * depl_frac(3) * factor * g_sigma
             depletion(2) = sumA * depl_frac(3) * factor * g_sigma
             kin_depl(1)  = sumA * depletion(1)
             kin_depl(2)  = rotB * depletion(2)

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), rl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( rot_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), rl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, rl, i, j, k, &
                  g, n, m, vel_grid, m_red, rot_trans, B, molecule )

          end if
       end if

       !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
       if( vib_flag )then
          if( v_modes(1) .gt. 0 )then
             depletion(1) = sumB * depl_frac(4) * factor * g_sigma
             depletion(2) = vibA * depl_frac(4) * factor * g_sigma
             kin_depl(1)  = vibA * depletion(1)
             kin_depl(2)  = sumB * depletion(2)

             call vr_deplete( vib_trans, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, A, molecule )

          end if

          if( v_modes(2) .gt. 0 )then
             depletion(1) = vibB * depl_frac(5) * factor * g_sigma
             depletion(2) = sumA * depl_frac(5) * factor * g_sigma
             kin_depl(1)  = sumA * depletion(1)
             kin_depl(2)  = vibB * depletion(2)

             call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
                  i(1), j(1), k(1), vl(1), frA, fvA, r_modes(1), v_modes(1) )
             call vr_deplete( vib_trans, phi(m), depletion(2), kin_depl(2), &
                  i(2), j(2), k(2), vl(2), frB, fvB, r_modes(2), v_modes(2) )
             call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
                  vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, B, molecule )

          end if
       end if

    end do


    ! deallocate energy arrays
    if( r_modes(1) .gt. 0 ) deallocate( frA )
    if( r_modes(2) .gt. 0 ) deallocate( frB )
    if( v_modes(1) .gt. 0 ) deallocate( fvA )
    if( v_modes(2) .gt. 0 ) deallocate( fvB )

    return
  end subroutine variance_reduction_kernel

  subroutine vr_depletion_routine( sumA, sumB, rotA, rotB, vibA, vibB, depl_frac, depl_sign, &
       phi, phi_old, frA, frB, fvA, fvB, rl, vl, i, j, k, g, n, m, vel_grid, m_red, molecule, &
       r_modes, v_modes, factor, g_sigma, rot_flag, vib_flag )

    use DistFunc
    use VelocityGrid
    use ReplenishingCollisions
    use PhysicalProperties

    implicit none

    type(VelocityGridType), dimension(:), intent(in) :: vel_grid
    type(DistFuncType), dimension(:), intent(in) :: phi_old
    type(MoleculeType), dimension(:), intent(in) :: molecule

    logical, intent(in) :: rot_flag, vib_flag

    integer, intent(in) :: n, m
    integer, dimension(2), intent(in) :: i, j, k, r_modes, v_modes
    integer, dimension(4), intent(in) :: rl, vl

    double precision, intent(in) :: sumA, sumB, rotA, rotB, vibA, vibB
    double precision, intent(in) :: g, m_red, factor, g_sigma
    double precision, dimension(5), intent(in) :: depl_frac, depl_sign
    double precision, dimension(:), intent(in) :: frA, frB, fvA, fvB

    type(DistFuncType), dimension(:) :: phi

    double precision, dimension(2) :: depletion, kin_depl

    !**ELASTIC*********************************************************************************************
    depletion(1) = sumB * depl_frac(1) * factor * g_sigma
    depletion(2) = sumA * depl_frac(1) * factor * g_sigma
    kin_depl(1)  = sumA * depletion(1)
    kin_depl(2)  = sumB * depletion(2)

    call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
         i(1), j(1), k(1), 0, frA, fvA, r_modes(1), v_modes(1) )
    call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
         i(2), j(2), k(2), 0, frB, fvB, r_modes(2), v_modes(2) )
    call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
         rl, i, j, k, g, n, m, vel_grid, m_red, elastic, A, molecule )

    !**ROTATIONAL - TRANSLATIONAL**************************************************************************
    if( rot_flag )then ! TODO: this call to rot_flag seems excessive....
       if( r_modes(1) .gt. 0 )then
          depletion(1) = sumB * depl_frac(2) * factor * g_sigma
          depletion(2) = rotA * depl_frac(2) * factor * g_sigma
          kin_depl(1)  = rotA * depletion(1)
          kin_depl(2)  = sumB * depletion(2)

          call vr_deplete( rot_trans, phi(n), depletion(1), kin_depl(1), &
               i(1), j(1), k(1), rl(1), frA, fvA, r_modes(1), v_modes(1) )
          call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
               i(2), j(2), k(2), rl(2), frB, fvB, r_modes(2), v_modes(2) )
          call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
               rl, i, j, k, g, n, m, vel_grid, m_red, rot_trans, A, molecule )

       end if

       if( r_modes(2) .gt. 0 )then
          depletion(1) = rotB * depl_frac(3) * factor * g_sigma
          depletion(2) = sumA * depl_frac(3) * factor * g_sigma
          kin_depl(1)  = sumA * depletion(1)
          kin_depl(2)  = rotB * depletion(2)

          call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
               i(1), j(1), k(1), rl(1), frA, fvA, r_modes(1), v_modes(1) )
          call vr_deplete( rot_trans, phi(m), depletion(2), kin_depl(2), &
               i(2), j(2), k(2), rl(2), frB, fvB, r_modes(2), v_modes(2) )
          call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, rl, i, j, k, &
               g, n, m, vel_grid, m_red, rot_trans, B, molecule )

       end if
    end if

    !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
    if( vib_flag )then
       if( v_modes(1) .gt. 0 )then
          depletion(1) = sumB * depl_frac(4) * factor * g_sigma
          depletion(2) = vibA * depl_frac(4) * factor * g_sigma
          kin_depl(1)  = vibA * depletion(1)
          kin_depl(2)  = sumB * depletion(2)

          call vr_deplete( vib_trans, phi(n), depletion(1), kin_depl(1), &
               i(1), j(1), k(1), vl(1), frA, fvA, r_modes(1), v_modes(1) )
          call vr_deplete( elastic, phi(m), depletion(2), kin_depl(2), &
               i(2), j(2), k(2), vl(2), frB, fvB, r_modes(2), v_modes(2) )
          call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
               vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, A, molecule )

       end if

       if( v_modes(2) .gt. 0 )then
          depletion(1) = vibB * depl_frac(5) * factor * g_sigma
          depletion(2) = sumA * depl_frac(5) * factor * g_sigma
          kin_depl(1)  = sumA * depletion(1)
          kin_depl(2)  = vibB * depletion(2)

          call vr_deplete( elastic, phi(n), depletion(1), kin_depl(1), &
               i(1), j(1), k(1), vl(1), frA, fvA, r_modes(1), v_modes(1) )
          call vr_deplete( vib_trans, phi(m), depletion(2), kin_depl(2), &
               i(2), j(2), k(2), vl(2), frB, fvB, r_modes(2), v_modes(2) )
          call replenish_collision( phi, depletion, kin_depl, frA, frB, fvA, fvB, &
               vl, i, j, k, g, n, m, vel_grid, m_red, vib_trans, B, molecule )

       end if
    end if

    return
  end subroutine vr_depletion_routine
