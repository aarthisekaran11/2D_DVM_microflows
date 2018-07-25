  do coll = 1, num_colls

     ! Pick velocities to collide
     call pick_collision_partner( i(1), j(1), k(1), glA, psi_eq(ns,n,same)%cumul_df, vel_grid(n) )
     call pick_collision_partner( i(2), j(2), k(2), glB, psi_neq(m,same)%cumul_rv_df, vel_grid(m) )

     ! Calculate relative velocity between chosen collision partners
     call compute_relative_velocity( i, j, k, g, vel_grid, n, m )

     ! Check for self collision
     if( g .lt. double_tol ) cycle

     ! Calculate g*sigma term
     g_sigma = sigma * g**vhs_exponent

     ! Temporary arrays for phi
     phiA = phi_old(n)%value( i(1), j(1), k(1) )
     phiB = phi_old(m)%value( i(2), j(2), k(2) )

     phiB_eq  = psi_eq(ns,m,same)%cumul_df( glB ) - psi_eq(ns,m,same)%cumul_df( glB-1 )
     phiB_neq = psi_neq(m,same)%cumul_rv_df( glB ) - psi_neq(m,same)%cumul_rv_df( glB-1 )

     phiB_rneq = psi_neq(m,same)%cumul_rot_df( glB )
     phiB_vneq = psi_neq(m,same)%cumul_vib_df( glB )

     signB = psi_neq(m,same)%sign( glB )
     phiB_neq = phiB_neq * signB

     ! Coefficients for inelastic collision partner
     B_rot = psi_neq(m,same)%c_rot_dev( glB )
     B_vib = psi_neq(m,same)%c_vib_dev( glB )

     ! Equilibrium collision partner
     if( r_modes(1) .gt. 0 )then
        rotA = psi_eq(ns,n,same)%rot(:)
        call pick_level( rl(1), rotA, one, r_levels(1) )

        ! rotA is strictly positive in the equilibrium case
        call normalized_array( frA, rotA, r_levels(1) )
     end if

     if( v_modes(1) .gt. 0 )then
        vibA = psi_eq(ns,n,same)%vib(:)
        call pick_level( vl(1), vibA, one, v_levels(1) )

        !vibA is strictly positive in the equilibrium case
        call normalized_array( fvA, vibA, v_levels(1) )
     end if


     ! Determine sign of non-equilibrium depletion part
     if( r_modes(2) .gt. 0 .or. v_modes(2) .gt. 0 )then
        if( r_modes(2) .gt. 0 .and. v_modes(2) .gt. 0 )then
           sign_test = one_half + phiB_neq / ( phiB_vneq + phiB_rneq )
        else if( r_modes(2) .gt. 0 )then
           sign_test = one_half * ( one + phiB_neq / phiB_rneq )
        else
           sign_test = one_half * ( one +  phiB_neq / phiB_vneq )
        end if

        array_sign = one
        if( get_rand() .gt. sign_test )then
           depl_sign(1) = -depl_sign(1)
           depl_sign(2) = -depl_sign(2)
           depl_sign(4) = -depl_sign(4)
           array_sign   = -one
        end if

     else
        depl_sign(1) = signB * depl_sign(1)

     end if

     ! Non-equilibrium collision partner
     if( r_modes(2) .gt. 0 )then
        rotB = phi_old(m)%rot(:,i(2),j(2),k(2)) - psi_eq(ns,m,same)%rot(:) * phiB_eq

        call pick_level( rl(2), rotB, one, r_levels(2) )
        depl_sign(3) = dsgn( rotB( rl(2) ) )

        call one_sign_normalized_array( frB, rotB, r_levels(2), array_sign )
        frB = B_vib * psi_eq(ns,m,same)%rot(:) + ( one - B_vib ) * frB
     end if

     if( v_modes(2) .gt. 0 )then
        vibB = phi_old(m)%vib(:,i(2),j(2),k(2)) - psi_eq(ns,m,same)%vib(:) * phiB_eq

        call pick_level( vl(2), vibB, one, v_levels(2) )
        depl_sign(5) = dsgn( vibB( vl(2) ) )

        call one_sign_normalized_array( fvB, vibB, v_levels(2), array_sign )
        fvB = B_rot * psi_eq(ns,m,same)%vib(:) + ( one - B_rot ) * fvB
     end if

     !**ELASTIC*********************************************************************************************
     depletion = depl_sign(1) * depl_frac(1) * factor * g_sigma

     call vr_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), 0, &
          frA, fvA, r_modes(1), v_modes(1) )
     call vr_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), 0, &
          frB, fvB, r_modes(2), v_modes(2) )
     call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
          g, n, m, vel_grid, m_red, elastic, A, molecule )

     !**ROTATIONAL - TRANSLATIONAL**************************************************************************
     if( rot_flag )then ! TODO: this call to rot_flag seems excessive....
        if( r_modes(1) .gt. 0 )then
           depletion = depl_sign(2) * depl_frac(2) * factor * g_sigma

           call vr_deplete( rot_trans, phi(n), depletion, i(1), j(1), k(1), rl(1), &
                frA, fvA, r_modes(1), v_modes(1) )
           call vr_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), rl(2), &
                frB, fvB, r_modes(2), v_modes(2) )
           call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
                g, n, m, vel_grid, m_red, rot_trans, A, molecule )

        end if


        if( r_modes(2) .gt. 0 )then 
           depletion = depl_sign(3) * depl_frac(3) * factor * g_sigma

           call vr_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), rl(1), &
                frA, fvA, r_modes(1), v_modes(1) )
           call vr_deplete( rot_trans, phi(m), depletion, i(2), j(2), k(2), rl(2), &
                frB, fvB, r_modes(2), v_modes(2) )
           call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
                g, n, m, vel_grid, m_red, rot_trans, B, molecule )

        end if
     end if

     !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
     if( vib_flag )then
        if( v_modes(1) .gt. 0 )then
           depletion = depl_sign(4) * depl_frac(4) * factor * g_sigma

           call vr_deplete( vib_trans, phi(n), depletion, i(1), j(1), k(1), vl(1), &
                frA, fvA, r_modes(1), v_modes(1) )
           call vr_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), vl(2), &
                frB, fvB, r_modes(2), v_modes(2) )
           call replenish_collision( phi, depletion, frA, frB, fvA, fvB, vl, i, j, k, &
                g, n, m, vel_grid, m_red, vib_trans, A, molecule )

        end if

        if( v_modes(2) .gt. 0 )then
           depletion = depl_sign(5) * depl_frac(5) * factor * g_sigma

           call vr_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), vl(1), &
                frA, fvA, r_modes(1), v_modes(1) )
           call vr_deplete( vib_trans, phi(m), depletion, i(2), j(2), k(2), vl(2), &
                frB, fvB, r_modes(2), v_modes(2) )
           call replenish_collision( phi, depletion, frA, frB, fvA, fvB, vl, i, j, k, &
                g, n, m, vel_grid, m_red, vib_trans, B, molecule )

        end if
     end if

  end do


! TWO==============================================================================================================
!==================================================================================================================
  do coll = 1, num_colls

     ! Pick velocities to collide
     call pick_collision_partner( i(1), j(1), k(1), glA, psi_neq(n,same)%cumul_rv_df, vel_grid(n) )
     call pick_collision_partner( i(2), j(2), k(2), glB, psi_eq(ns,m,same)%cumul_df, vel_grid(m) )

     ! Calculate relative velocity between chosen collision partners
     call compute_relative_velocity( i, j, k, g, vel_grid, n, m )

     ! Check for self collision
     if( g .lt. double_tol ) cycle

     ! Calculate g*sigma term
     g_sigma = sigma * g**vhs_exponent

     ! Temporary arrays for phi
     phiA = phi_old(n)%value( i(1), j(1), k(1) )
     phiB = phi_old(m)%value( i(2), j(2), k(2) )

     phiA_eq  = psi_eq(ns,n,same)%cumul_df( glA ) - psi_eq(ns,n,same)%cumul_df( glA-1 )
     phiA_neq = psi_neq(n,same)%cumul_rv_df( glA ) - psi_neq(n,same)%cumul_rv_df( glA-1 )

     phiA_rneq = psi_neq(n,same)%cumul_rot_df( glA )
     phiA_vneq = psi_neq(n,same)%cumul_vib_df( glA )

     signA = psi_neq(n,same)%sign( glA )
     phiA_neq = phiA_neq * signA

     ! Coefficients for inelastic collision partner
     A_rot = psi_neq(n,same)%c_rot_dev( glA )
     A_vib = psi_neq(n,same)%c_vib_dev( glA )

     ! Determine sign of non-equilibrium depletion part
     if( r_modes(1) .gt. 0 .or. v_modes(1) .gt. 0 )then
        if( r_modes(1) .gt. 0 .and. v_modes(1) .gt. 0 )then
           sign_test = one_half + phiA_neq / ( phiA_vneq + phiA_rneq )
        else if( r_modes(1) .gt. 0 )then
           sign_test = one_half * ( one + phiA_neq / phiA_rneq )
        else
           sign_test = one_half * ( one +  phiA_neq / phiA_vneq )
        end if

        array_sign = one
        if( get_rand() .gt. sign_test )then
           depl_sign(1) = -depl_sign(1)
           depl_sign(3) = -depl_sign(3)
           depl_sign(5) = -depl_sign(5)
           array_sign   = -one
        end if

     else
        depl_sign(1) = signA * depl_sign(1)

     end if

     ! Non-equilibrium collision partner
     if( r_modes(1) .gt. 0 )then
        rotA = phi_old(n)%rot(:,i(1),j(1),k(1)) - psi_eq(ns,n,same)%rot(:) * phiA_eq

        call pick_level( rl(1), rotA, one, r_levels(1) )
        depl_sign(2) = dsgn( rotA( rl(1) ) )

        call one_sign_normalized_array( frA, rotA, r_levels(1), array_sign )
        frA = A_vib * psi_eq(ns,n,same)%rot(:) + ( one - A_vib ) * frA
     end if

     if( v_modes(1) .gt. 0 )then
        vibA = phi_old(n)%vib(:,i(1),j(1),k(1)) - psi_eq(ns,n,same)%vib(:) * phiA_eq

        call pick_level( vl(1), vibA, one, v_levels(1) )
        depl_sign(4) = dsgn( vibA( vl(1) ) )

        call one_sign_normalized_array( fvA, vibA, v_levels(1), array_sign )
        fvA = A_rot * psi_eq(ns,n,same)%vib(:) + ( one - A_rot ) * fvA
     end if

     ! Determine sign of non-equilibrium depletion part
     if( r_modes(2) .gt. 0 .or. v_modes(2) .gt. 0 )then
        if( r_modes(2) .gt. 0 .and. v_modes(2) .gt. 0 )then
           sign_test = one_half + phiB_neq / ( phiB_vneq + phiB_rneq )
        else if( r_modes(2) .gt. 0 )then
           sign_test = one_half * ( one + phiB_neq / phiB_rneq )
        else
           sign_test = one_half * ( one +  phiB_neq / phiB_vneq )
        end if

        array_sign = one
        if( get_rand() .gt. sign_test )then
           depl_sign(1) = -depl_sign(1)
           depl_sign(2) = -depl_sign(2)
           depl_sign(4) = -depl_sign(4)
           array_sign   = -one
        end if

     else
        depl_sign(1) = signB * depl_sign(1)

     end if

     ! Non-equilibrium collision partner
     if( r_modes(2) .gt. 0 )then
        rotB = phi_old(m)%rot(:,i(2),j(2),k(2)) - psi_eq(ns,m,same)%rot(:) * phiB_eq

        call pick_level( rl(2), rotB, one, r_levels(2) )
        depl_sign(3) = depl_sign(3) * dsgn( rotB( rl(2) ) )

        call one_sign_normalized_array( frB, rotB, r_levels(2), array_sign )
        frB = B_vib * psi_eq(ns,m,same)%rot(:) + ( one - B_vib ) * frB
     end if

     if( v_modes(2) .gt. 0 )then
        vibB = phi_old(m)%vib(:,i(2),j(2),k(2)) - psi_eq(ns,m,same)%vib(:) * phiB_eq

        call pick_level( vl(2), vibB, one, v_levels(2) )
        depl_sign(5) = depl_sign(5) * dsgn( vibB( vl(2) ) )

        call one_sign_normalized_array( fvB, vibB, v_levels(2), array_sign )
        fvB = B_rot * psi_eq(ns,m,same)%vib(:) + ( one - B_rot ) * fvB
     end if

     !**ELASTIC*********************************************************************************************
     depletion = depl_sign(1) * depl_frac(1) * factor * g_sigma

     call vr_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), 0, &
          frA, fvA, r_modes(1), v_modes(1) )
     call vr_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), 0, &
          frB, fvB, r_modes(2), v_modes(2) )
     call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
          g, n, m, vel_grid, m_red, elastic, A, molecule )

     !**ROTATIONAL - TRANSLATIONAL**************************************************************************
     if( rot_flag )then ! TODO: this call to rot_flag seems excessive....
        if( r_modes(1) .gt. 0 )then
           depletion = depl_sign(2) * depl_frac(2) * factor * g_sigma

           call vr_deplete( rot_trans, phi(n), depletion, i(1), j(1), k(1), rl(1), &
                frA, fvA, r_modes(1), v_modes(1) )
           call vr_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), rl(2), &
                frB, fvB, r_modes(2), v_modes(2) )
           call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
                g, n, m, vel_grid, m_red, rot_trans, A, molecule )

        end if


        if( r_modes(2) .gt. 0 )then 
           depletion = depl_sign(3) * depl_frac(3) * factor * g_sigma

           call vr_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), rl(1), &
                frA, fvA, r_modes(1), v_modes(1) )
           call vr_deplete( rot_trans, phi(m), depletion, i(2), j(2), k(2), rl(2), &
                frB, fvB, r_modes(2), v_modes(2) )
           call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
                g, n, m, vel_grid, m_red, rot_trans, B, molecule )

        end if
     end if

     !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
     if( vib_flag )then
        if( v_modes(1) .gt. 0 )then
           depletion = depl_sign(4) * depl_frac(4) * factor * g_sigma

           call vr_deplete( vib_trans, phi(n), depletion, i(1), j(1), k(1), vl(1), &
                frA, fvA, r_modes(1), v_modes(1) )
           call vr_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), vl(2), &
                frB, fvB, r_modes(2), v_modes(2) )
           call replenish_collision( phi, depletion, frA, frB, fvA, fvB, vl, i, j, k, &
                g, n, m, vel_grid, m_red, vib_trans, A, molecule )

        end if

        if( v_modes(2) .gt. 0 )then
           depletion = depl_sign(5) * depl_frac(5) * factor * g_sigma

           call vr_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), vl(1), &
                frA, fvA, r_modes(1), v_modes(1) )
           call vr_deplete( vib_trans, phi(m), depletion, i(2), j(2), k(2), vl(2), &
                frB, fvB, r_modes(2), v_modes(2) )
           call replenish_collision( phi, depletion, frA, frB, fvA, fvB, vl, i, j, k, &
                g, n, m, vel_grid, m_red, vib_trans, B, molecule )

        end if
     end if

  end do

! THREE============================================================================================================
!==================================================================================================================
  do coll = 1, num_colls

     ! Pick velocities to collide
     call pick_collision_partner( i(1), j(1), k(1), glA, psi_neq(n,same)%cumul_rv_df, vel_grid(n) )
     call pick_collision_partner( i(2), j(2), k(2), glB, psi_neq(m,same)%cumul_rv_df, vel_grid(m) )

     ! Calculate relative velocity between chosen collision partners
     call compute_relative_velocity( i, j, k, g, vel_grid, n, m )

     ! Check for self collision
     if( g .lt. double_tol ) cycle

     ! Calculate g*sigma term
     g_sigma = sigma * g**vhs_exponent

     ! Temporary arrays for phi
     phiA = phi_old(n)%value( i(1), j(1), k(1) )
     phiB = phi_old(m)%value( i(2), j(2), k(2) )

     phiA_eq  = psi_eq(ns,n,same)%cumul_df( glA ) - psi_eq(ns,n,same)%cumul_df( glA-1 )
     phiA_neq = psi_neq(n,same)%cumul_rv_df( glA ) - psi_neq(n,same)%cumul_rv_df( glA-1 )
     phiB_eq  = psi_eq(ns,m,same)%cumul_df( glB ) - psi_eq(ns,m,same)%cumul_df( glB-1 )
     phiB_neq = psi_neq(m,same)%cumul_rv_df( glB ) - psi_neq(m,same)%cumul_rv_df( glB-1 )

     phiA_rneq = psi_neq(n,same)%cumul_rot_df( glA )
     phiA_vneq = psi_neq(n,same)%cumul_vib_df( glA )
     phiB_rneq = psi_neq(m,same)%cumul_rot_df( glB )
     phiB_vneq = psi_neq(m,same)%cumul_vib_df( glB )

     signA = psi_neq(n,same)%sign( glA )
     phiA_neq = phiA_neq * signA
     signB = psi_neq(m,same)%sign( glB )
     phiB_neq = phiB_neq * signB

     ! Coefficients for inelastic collision partner
     A_rot = psi_neq(n,same)%c_rot_dev( glA )
     A_vib = psi_neq(n,same)%c_vib_dev( glA )
     B_rot = psi_neq(m,same)%c_rot_dev( glB )
     B_vib = psi_neq(m,same)%c_vib_dev( glB )

     ! Determine sign of non-equilibrium depletion part
     if( r_modes(1) .gt. 0 .or. v_modes(1) .gt. 0 )then
        if( r_modes(1) .gt. 0 .and. v_modes(1) .gt. 0 )then
           sign_test = one_half + phiA_neq / ( phiA_vneq + phiA_rneq )
        else if( r_modes(1) .gt. 0 )then
           sign_test = one_half * ( one + phiA_neq / phiA_rneq )
        else
           sign_test = one_half * ( one +  phiA_neq / phiA_vneq )
        end if

        array_sign = one
        if( get_rand() .gt. sign_test )then
           depl_sign(1) = -depl_sign(1)
           depl_sign(3) = -depl_sign(3)
           depl_sign(5) = -depl_sign(5)
           array_sign   = -one
        end if

     else
        depl_sign(1) = signA * depl_sign(1)

     end if

     ! Non-equilibrium collision partner
     if( r_modes(1) .gt. 0 )then
        rotA = phi_old(n)%rot(:,i(1),j(1),k(1)) - psi_eq(ns,n,same)%rot(:) * phiA_eq

        call pick_level( rl(1), rotA, one, r_levels(1) )
        depl_sign(2) = dsgn( rotA( rl(1) ) )

        call one_sign_normalized_array( frA, rotA, r_levels(1), array_sign )
        frA = A_vib * psi_eq(ns,n,same)%rot(:) + ( one - A_vib ) * frA
     end if

     if( v_modes(1) .gt. 0 )then
        vibA = phi_old(n)%vib(:,i(1),j(1),k(1)) - psi_eq(ns,n,same)%vib(:) * phiA_eq

        call pick_level( vl(1), vibA, one, v_levels(1) )
        depl_sign(4) = dsgn( vibA( vl(1) ) )

        call one_sign_normalized_array( fvA, vibA, v_levels(1), array_sign )
        fvA = A_rot * psi_eq(ns,n,same)%vib(:) + ( one - A_rot ) * fvA
     end if


     ! Equilibrium collision partner
     if( r_modes(2) .gt. 0 )then
        rotB = psi_eq(ns,m,same)%rot(:)
        call pick_level( rl(2), rotB, one, r_levels(2) )

        ! rotB is strictly positive in the equilibrium case
        call normalized_array( frB, rotB, r_levels(2) )
     end if

     if( v_modes(2) .gt. 0 )then
        vibB = psi_eq(ns,m,same)%vib(:)
        call pick_level( vl(2), vibB, one, v_levels(2) )

        !vibB is strictly positive in the equilibrium case
        call normalized_array( fvB, vibB, v_levels(2) )
     end if

     !**ELASTIC*********************************************************************************************
     depletion = depl_sign(1) * depl_frac(1) * factor * g_sigma

     call vr_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), 0, &
          frA, fvA, r_modes(1), v_modes(1) )
     call vr_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), 0, &
          frB, fvB, r_modes(2), v_modes(2) )
     call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
          g, n, m, vel_grid, m_red, elastic, A, molecule )

     !**ROTATIONAL - TRANSLATIONAL**************************************************************************
     if( rot_flag )then ! TODO: this call to rot_flag seems excessive....
        if( r_modes(1) .gt. 0 )then
           depletion = depl_sign(2) * depl_frac(2) * factor * g_sigma

           call vr_deplete( rot_trans, phi(n), depletion, i(1), j(1), k(1), rl(1), &
                frA, fvA, r_modes(1), v_modes(1) )
           call vr_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), rl(2), &
                frB, fvB, r_modes(2), v_modes(2) )
           call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
                g, n, m, vel_grid, m_red, rot_trans, A, molecule )

        end if


        if( r_modes(2) .gt. 0 )then 
           depletion = depl_sign(3) * depl_frac(3) * factor * g_sigma

           call vr_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), rl(1), &
                frA, fvA, r_modes(1), v_modes(1) )
           call vr_deplete( rot_trans, phi(m), depletion, i(2), j(2), k(2), rl(2), &
                frB, fvB, r_modes(2), v_modes(2) )
           call replenish_collision( phi, depletion, frA, frB, fvA, fvB, rl, i, j, k, &
                g, n, m, vel_grid, m_red, rot_trans, B, molecule )

        end if
     end if

     !**VIBRATIONAL - TRANSLATIONAL*************************************************************************
     if( vib_flag )then
        if( v_modes(1) .gt. 0 )then
           depletion = depl_sign(4) * depl_frac(4) * factor * g_sigma

           call vr_deplete( vib_trans, phi(n), depletion, i(1), j(1), k(1), vl(1), &
                frA, fvA, r_modes(1), v_modes(1) )
           call vr_deplete( elastic, phi(m), depletion, i(2), j(2), k(2), vl(2), &
                frB, fvB, r_modes(2), v_modes(2) )
           call replenish_collision( phi, depletion, frA, frB, fvA, fvB, vl, i, j, k, &
                g, n, m, vel_grid, m_red, vib_trans, A, molecule )

        end if

        if( v_modes(2) .gt. 0 )then
           depletion = depl_sign(5) * depl_frac(5) * factor * g_sigma

           call vr_deplete( elastic, phi(n), depletion, i(1), j(1), k(1), vl(1), &
                frA, fvA, r_modes(1), v_modes(1) )
           call vr_deplete( vib_trans, phi(m), depletion, i(2), j(2), k(2), vl(2), &
                frB, fvB, r_modes(2), v_modes(2) )
           call replenish_collision( phi, depletion, frA, frB, fvA, fvB, vl, i, j, k, &
                g, n, m, vel_grid, m_red, vib_trans, B, molecule )

        end if
     end if

  end do
