#if defined(A04_3D)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Capture of one ice category by
!  another.
! Subroutine Interface:
      SUBROUTINE LSP_COLLECTION(                                        &
     &  points, timestep                                                &
                                          ! Number of points and tstep
     &, qcf1, qcf2, T                                                   &
                                          ! Water contents and temp
     &, area_mix, area_ice                                              &
     &, cficei                                                          &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
!    &, cf, cff                           ! Current cloud fractions for
!                                         ! updating
     &, rho, rhor, m0, tcg1, tcg1i, tcg2, tcg2i                         &
                                          ! Parametrization information
     &, corr, cx, constp, ai, bi, aic, bic, ice_type1, ice_type2        &
                                          ! Microphysical information
     &, l_update_cf, l_seq, l_psd_2 , l_use_area, l_no_t_check          &
                                          ! Code options
     &, ptransfer, iterations                                           &
                                          ! Mass transfer diagnostic
!    &, cftransfer, cfftransfer           ! Cloud transfer diagnostics
     &  )
!
      Implicit None
!
! Purpose:
!   Update cloud prognostics as a result of collisions between different
!   categories of ice particles
!
!  Method:
!   Solve the microphysical transfer equation for a specified
!   distribution of ice particles (species 1) sweeping out a different
!   specified distribution of ice particles (species 2) and adding
!   combining the mass into species 1.
!   Note that the current formulation is limited to the collection of
!   snow by graupel or the collection of ice by snow. If we want to do
!   graupel to ice then we need to look again at the prefactors
!   constp(16 to 19 +ice_offset1) because they depend on both the
!   species, not just species 1.
!   The commented out code referring to cloud fraction updates is to
!   highlight where in the subroutine changes would need to be made
!   if it was later decided to update cloud fractions in some way.
!   We note that we are only allowing use of the generic ice particle
!   size distribution for second species, intending its use for the
!   collection of aggregates (generic psd) by graupel and *not* the
!   collection of crystals by aggregates (this would go against the
!   idea that the generic psd represents both the crystals and the
!   aggregates). The logic in the call to the routine should prevent
!   any inconsistencies since l_mcr_qcf2 and l_psd should not both be
!   true.
!
! Current Owner of Code: Jonathan Wilkinson
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP 26.
!
!
! Subroutine Arguments
!
! Obtain the size for cx and constp and define ice_type_offset
#include "c_lspsiz.h"
#include "c_0_dg_c.h"
!
      Integer, Intent(In) ::                                            &
     &  points                                                          &
                          ! Number of points to calculate
     &, iterations                                                      &
                          ! Number of microphysics iterations
                          ! (for rate diagnostic)
     &, ice_type1                                                       &
     &, ice_type2
                          ! Type of ice (0 - crystals, 1 - aggregates
                          !              3 - graupel)
!
      Real, Intent(In) ::                                               &
     &  Timestep                                                        &
                          ! Timestep / s
     &, area_mix(points)                                                &
                          ! Fraction of gridbox with mixed phase cloud
     &, area_ice(points)                                                &
                          ! Fraction of gridbox with ice-only cloud
     &, cficei(points)                                                  &
                          ! 1/Fraction of gridbox with ice cloud
!    &, cf(points)        ! Current cloud fraction
!    &, cff(points)       ! Current ice cloud fraction
     &, rho(points)                                                     &
                          ! Air density / kg m-3
     &, rhor(points)                                                    &
                          ! 1/Air density / m3 kg-1
     &, m0                                                              &
                          ! Seed ice water content / kg kg-1
     &, tcg1(points)                                                    &
                          ! T dependent function in ice size dist'n
     &, tcg1i(points)                                                   &
                          ! 1/tcg (no units)
     &, tcg2(points)                                                    &
                          ! T dependent function in ice size dist'n
     &, tcg2i(points)                                                   &
                          ! 1/tcg (no units)
     &, corr(points)                                                    &
                          ! Fall velocity correction factor (no units)
     &, ai, bi, aic, bic
                          ! Mass size parametrization m(D) = ai D^bi
!
      Real, Intent(InOut) ::                                            &
     &  qcf1(points)                                                    &
                          ! Ice water content of 1st species / kg kg-1
     &, qcf2(points)                                                    &
                          ! Ice water content of 2nd species / kg kg-1
     &, T(points)                                                       &
                          ! Temperature / K
     &, ptransfer(points)
                          ! Mass rimed in this timestep / kg kg-1
!    &, cftransfer(points) ! Cumulative cloud fraction increment
!    &, cfftransfer(points)! Cumulative ice cloud fraction increment
!
      Logical, Intent(In) ::                                            &
     &  l_update_cf                                                     &
                          ! Update cloud fractions
     &, l_seq                                                           &
                          ! Carry out sequential updating
     &, l_psd_2                                                         &
                          ! Use generic ice particle size distribution
     &, l_use_area                                                      &
                          ! Use ice area amount in calculating
                          ! gridbox mean transfer rate
     &, l_no_t_check
                          ! Do not check that temperature is less
                          ! than 0 deg C before proceeding
!
! Local Variables
!
      Integer                                                           &
     &  i                                                               &
                          ! Loop counter for points
     &, ice_offset1                                                     &
                          ! Index offset for first ice species
     &, ice_offset2 
                          ! Index offset for second ice species
!
      Real                                                              &
     &  dqi(points)                                                     &
                          ! Transfer of mixing ratio  / kg kg-1
     &, vi1(points)                                                     &
                          ! Average fall speed of first ice
                          ! particle species  / m s-1
     &, vi2(points)                                                     &
                          ! Average fall speed of second ice
                          ! particle species  / m s-1
     &, fv1(points)                                                     &
                          ! Average fall speed difference between
                          ! ice species  / m s-1
     &, lami1(points)                                                   &
                          ! Reciprocal of slope parameter in
                          ! first ice particle size distribution  / m
     &, lami2(points)                                                   &
                          ! Reciprocal of slope parameter in
                          ! second ice particle size distribution  / m
     &, lamfac1(points)                                                 &
                          ! Combination of lamr1 and lamr2
     &, collision_eff(points)                                           &
                          ! Collision efficiency / no units
     &, m_0(points), m_1(points), m_2(points)                           &
                          ! zero, 1st and 2nd moments of the ice PSD
     &, m_bi_di(points)
                          ! bi+di moment of the generic ice size distn

      !-----------------------------------------------
      ! Select appropriate ice parametrization (see c_lspsiz)
      !-----------------------------------------------
      ice_offset1=ice_type1*ice_type_offset
      ice_offset2=ice_type2*ice_type_offset

      !-----------------------------------------------
      ! Calculate moments of size distribution if appropriate 
      !----------------------------------------------- 
      If (l_psd_2) then
        ! Calculate the 0th, 1st, 2nd and bi+di (cx(82)) moments of the
        ! ice particle size distribution 
! DEPENDS ON: lsp_moments 
        Call lsp_moments(points,rho,T,qcf2,cficei,                      & 
     &                   ai,bi,0.0,m_0)
! DEPENDS ON: lsp_moments 
        Call lsp_moments(points,rho,T,qcf2,cficei,                      & 
     &                   ai,bi,1.0,m_1)
! DEPENDS ON: lsp_moments 
        Call lsp_moments(points,rho,T,qcf2,cficei,                      & 
     &                   ai,bi,2.0,m_2)
! DEPENDS ON: lsp_moments  
        Call lsp_moments(points,rho,T,qcf2,cficei,                      &  
     &                   ai,bi,cx(82),m_bi_di) 
      End if  ! l_psd_2

      Do i = 1, points

        If (qcf1(i)  >   m0 .and. qcf2(i)  >   m0                       &
     &      .and. ( T(i) < zerodegc .or. l_no_t_check)                  &
     &      .and. ( (area_ice(i)+area_mix(i))  >   0.0                  &
     &             .or. (.not. l_use_area) ) ) Then

          !-----------------------------------------------
          ! Calculate first ice mass-weighted fallspeed
          !-----------------------------------------------
          ! Use size distribution based on intercepts
          vi1(i) = constp(4+ice_offset1) * corr(i) *                    &
     &           ( rho(i) * qcf1(i) * cficei(i)                         &
     &         * constp(5+ice_offset1) * tcg1i(i))**cx(3+ice_offset1)

          !-----------------------------------------------
          ! Calculate second ice mass-weighted fallspeed
          !-----------------------------------------------
          If (l_psd_2) then
            ! Use generic ice size distribution
            vi2(i) = constp(82) * corr(i) * m_bi_di(i)                  &
     &                 / (rho(i) * qcf1(i) * cficei(i))
          Else
            ! Use size distribution based on intercepts
            vi2(i) = constp(4+ice_offset2) * corr(i) *                  &
     &             ( rho(i) * qcf2(i) * cficei(i)                       &
     &           * constp(5+ice_offset2) * tcg2i(i))**cx(3+ice_offset2)
          End if  ! l_psd_2

          !-----------------------------------------------
          ! Estimate the mean absolute differences in velocities
          !-----------------------------------------------
          fv1(i) = max(abs(vi1(i)-vi2(i)),(vi1(i)+vi2(i))/8.0)

          !-----------------------------------------------
          ! Calculate reciprocal of lambda for first ice distribution
          !-----------------------------------------------
          lami1(i) = (rho(i) * qcf1(i) * cficei(i)                      &
     &      *constp(5+ice_offset1)*tcg1i(I))**(-cx(7+ice_offset1))

          !-----------------------------------------------
          ! Calculate reciprocal of lambda for second ice distribution
          !-----------------------------------------------
          If (.not. l_psd_2) then
            lami2(i) = (rho(i) * qcf2(i) * cficei(i)                    &
     &        *constp(5+ice_offset2)*tcg2i(I))**(-cx(7+ice_offset2))
          End if  ! l_psd_2

          !------------------------------------------------
          ! Calculate transfer
          !------------------------------------------------
          collision_eff(i) = 0.02*exp(0.08*(T(i)-zerodegc))

          If (l_psd_2) then
            ! Use the generic ice particle size distribution
            ! constp(80)=pi**2/24 x1g ag (80 = 20+ice_offset for graup)
            ! constp(91)=gamma(bi+3+x4i)
            ! constp(92)=2 gamma(bi+2+x4i)
            ! constp(93)=gamma(bi+1+x4i)
            dqi(i) = collision_eff(i) * fv1(i) *  timestep * rhor(i) *  &
     &               tcg1(i) * constp(20+ice_offset1) *                 &
     &               lami1(i)**(-cx(11+ice_offset1))  *                 &
     &            (constp(91) * m_0(i) * lami1(i)**CX(13+ice_offset1) + &
     &             constp(92) * m_1(i) * lami1(i)**CX(14+ice_offset1) + &
     &             constp(93) * m_2(i) * lami1(i)**CX(15+ice_offset1))
          Else
            ! Use the distribution defined by intercepts
            lamfac1(i) =                                                &
     &         constp(16+ice_offset1)*(lami2(i)**CX(8+ice_offset2)      &
     &                               * lami1(i)**CX(13+ice_offset1)) +  &
     &         constp(17+ice_offset1)*(lami2(i)**CX(9+ice_offset2)      &
     &                               * lami1(i)**CX(14+ice_offset1)) +  &
     &         constp(18+ice_offset1)*(lami2(i)**CX(10+ice_offset2)     &
     &                               * lami1(i)**CX(15+ice_offset1))

            dqi(i) = collision_eff(i) * tcg1(i) * tcg2(i)               &
     &               * constp(19+ice_offset1)                           &
     &               * lami1(i)**(-cx(11+ice_offset1))                  &
     &               * lami2(i)**(-cx(11+ice_offset2))                  &
     &               * fv1(i) * lamfac1(i) * timestep * rhor(i)

          End if  ! l_psd_2

          If (l_use_area) then
            ! The calculations above have used in-cloud quantities to
            ! obtain the transfer rates. Multiply the transfer by the
            ! ice cloud amount to get a gridbox mean.
            dqi(i) = dqi(i) * (area_mix(i) + area_ice(i))
          End if

          !------------------------------------------------
          ! Limit transfer to the mass of species 2 that is available
          !------------------------------------------------
          dqi(i) = min(dqi(i),qcf2(i))

          ptransfer(i) = ptransfer(i) + dqi(i) / (timestep*iterations)

          !------------------------------------------------
          ! Adjust ice species contents
          !------------------------------------------------
          If (l_seq) then
            qcf1(i)  = qcf1(i)   + dqi(i)
            qcf2(i)  = qcf2(i)   - dqi(i)
          End if

          !------------------------------------------------
          ! Update cloud fractions
          !------------------------------------------------
!         These are commented out since there is currently no
!         cloud fraction update associated with the collision terms.
!          cf_transfer_rate(i)  = 0.0 / (timestep*iterations)
!          cff_transfer_rate(i) = 0.0 / (timestep*iterations)
!
!          If (l_update_cf) then
!            cftransfer(i)  = cftransfer(i)  + cf_transfer_rate(i)
!            cfftransfer(i) = cfftransfer(i) + cff_transfer_rate(i)
!
!            If (l_seq) then
!              cf(i)  = cf(i)  +cf_transfer_rate(i) *timestep*iterations
!              cff(i) = cff(i) +cff_transfer_rate(i)*timestep*iterations
!            End if  ! l_seq
!
!          End if  ! l_update_cf

        End If ! qcf1(i) >  m0 etc.

      End Do  ! Points

      Return  ! End of the subroutine
      END SUBROUTINE LSP_COLLECTION
#endif
