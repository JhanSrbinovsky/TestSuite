#if defined(A09_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Area cloud parameterisation for use with PC2 Cloud Scheme.
!
      Subroutine pc2_arcld(                                             &
!      Pressure related fields
     & p_layer_centres,p_layer_boundaries, ccb, cumulus, rhcrit         &
!      Array dimensions
     &,levels, row_length,rows,rhc_row_length,rhc_rows                  &
     &,large_levels,levels_per_level                                    &
!      Prognostic Fields
     &,CF_area, T, CF, CFL, CFF, Q, qcl, qcf, RHTS                      &
     &,tlts,qtts,ptts                                                   &
!      Logical control
     &,l_mixing_ratio)
!
      IMPLICIT NONE
!
! Description:
!   Cusack-like vertical interpolation onto 3 sub-levels to calculate
!   cloud fraction and condensate in PC2 initiation. Also area cloud
!   parameterisation for use with radiation scheme.
!
! Method:
!   See the PC2 documentation.
!
! Current Code Owner: Cyril Morcrette
!
! History:
! Version   Date     Comment.
! -------   ----     --------
! 7.3       25/11/09 Original code.  Ian Boutle.
!
! Code Description:
!   Language: fortran 77 + common extensions
!   This code is written to UMDP3 v6 programming standards.
!
!  Global Variables:----------------------------------------------------
#include "c_lheat.h"
#include "c_r_cp.h"
!
!  Subroutine Arguments:------------------------------------------------
      INTEGER                                                           &
                        !, INTENT(IN)
     & levels                                                           &
!       No. of levels being processed. (wet_levels)
     &, row_length,rows                                                 &
!       Row length and number of rows being processed.
     &, rhc_row_length,rhc_rows                                         &
!       Dimensions of the rhcrit variable.
     &,ccb(row_length,rows)                                             &
!       convective cloud base
     &, large_levels                                                    &
!       Total no. of sub-levels being processed by cloud scheme.
!       Currently ((levels - 2)*levels_per_level) + 2
     &, levels_per_level
!       No. of sub-levels being processed by area cloud scheme.
!       Want an odd number of sublevels per level.
!       NB: levels_per_level = 3 is currently hardwired in the do loops
!
      LOGICAL                                                           &
                        !, INTENT(IN)
     & cumulus(row_length,rows)                                         &
!       Is this a boundary layer cumulus point
     &,l_mixing_ratio  ! Use mixing ratio formulation
!
      REAL                                                              &
                        !, INTENT(IN)
     &  p_layer_centres(row_length,rows,0:levels)                       &
!       pressure at all points, on theta levels (Pa).
     &, p_layer_boundaries(row_length,rows,0:levels)                    &
!       pressure at all points, on u,v levels (Pa).
     &,rhcrit(rhc_row_length,rhc_rows,levels)                           &
!       Critical relative humidity.  See the the paragraph incorporating
!       eqs P292.11 to P292.14; the values need to be tuned for the give
!       set of levels.
     &,CFF(row_length,rows,levels)                                      &
!       Ice cloud fraction (no units)
     &,qcf(row_length,rows,levels)
!       Cloud ice content at processed levels (kg water per kg air).
!
      REAL                                                              &
                        !, INTENT(INOUT)
     & T(row_length,rows,levels)                                        &
!       Temperature (K)
     &,CF(row_length,rows,levels)                                       &
!       Total cloud fraction (no units)
     &,CFL(row_length,rows,levels)                                      &
!       Liquid cloud fraction (no units)
     &,Q(row_length,rows,levels)                                        &
!       Vapour content (kg water per kg air)
     &,qcl(row_length,rows,levels)                                      &
!       Liquid content (kg water per kg air)
     &,RHTS(row_length,rows,levels)                                     &
!       Variable carrying initial RHT wrt TL from start of timestep
     &, tlts(row_length,rows,levels)                                    &
!       TL at start of timestep
     &, qtts(row_length,rows,levels)                                    &
!       qT at start of timestep
     &, ptts(row_length,rows,levels)                                    &
!       pressure at theta levels at start of timestep
     &,CF_AREA(row_length, rows, levels) !area cloud fraction
!
! --------------------------------------------------------------------
! Local variables
! ---------------------------------------------------------------------
      INTEGER i,j,k      ! Loop counters: k - vertical level index.
!                                       i,j - horizontal field index.
      INTEGER k_index    ! Extra loop counter for large arrays.
!
      REAL                                                              &
     &  inverse_level                                                   &
                       ! Set to (1. / levels_per_level)
     &, qt_norm_next                                                    &
                       ! Temporary space for qT_norm
     &, stretcher                                                       &
     &, delta_p        ! Layer pressure thickness * inverse_level
!
      REAL                                                              &
     &  qsl(row_length,rows)                                            &
!        Saturated specific humidity for temp TL or T.
     &, qcl_latest(row_length,rows,levels)                              &
!       Cloud liquid content at processed levels (kg water per kg air).
     &, TL(row_length,rows,levels)                                      &
!       liquid temperature (TL) (K).
     &, qt_norm(row_length,rows)                                        &
!        Total water content normalized to qSAT_WAT.
     &, rhcrit_large(rhc_row_length,rhc_rows,large_levels)              &
!        values of quantities on large_levels
     &, p_large(row_length,rows,large_levels)                           &
!
     &, T_large(row_length,rows,large_levels)                           &
!
     &, q_large(row_length,rows,large_levels)                           &
!
     &, qcl_large(row_length,rows,large_levels)                         &
!
     &, CF_large(row_length,rows,large_levels)                          &
!
     &, CFL_large(row_length,rows,large_levels)                         &
!
     &, CFF_large(row_length,rows,large_levels)                         &
!
     &, rhts_large(row_length,rows,large_levels)                        &
!
     &, tlts_large(row_length,rows,large_levels)                        &
!
     &, qtts_large(row_length,rows,large_levels)                        &
!
     &, ptts_large(row_length,rows,large_levels)
!
      LOGICAL                                                           &
     & linked(row_length,rows,levels)
!       True for sub-layers that have similar supersaturation properties
!
!  Local parameters and other physical constants------------------------
      REAL                                                              &
     & LCRCP
!       Latent heat of condensation divided by heat capacity of air.
!
      PARAMETER(                                                        &
     &          LCRCP=LC/CP                                             &
     &          )
      REAL drat_thresh            ! Test for continuity of sub-levels
      REAL tol_test               ! Tolerance for non-zero humidities
      PARAMETER (drat_thresh=3.0E-1, tol_test=1.0E-11)
!
! ---------------------------------------------------------------------
! Code starts
! ---------------------------------------------------------------------
      inverse_level = 1. / levels_per_level
!
! create new arrays for TL and current qcl
!
      Do k = 1, levels
        Do j = 1, rows
          Do i = 1, row_length
            TL(i,j,k) = T(i,j,k) - LCRCP*qcl(i,j,k)
            qcl_latest(i,j,k) = qcl(i,j,k)
          End Do
        End Do
      End Do
!
! Test for continuity between adjacent layers based on supersaturation
! (qt - qsl) / qsl : as we take differences the - qsl term drops out.
! DEPENDS ON: qsat_wat_mix
        CALL qsat_wat_mix( qsl, TL(1,1,1),                              &
     &       p_layer_centres(1,1,1), row_length*rows, l_mixing_ratio)
        DO j = 1, rows
          DO i = 1, row_length
            qt_norm(i,j) =(Q(i,j,1)+qcl_latest(i,j,1)                   &
     &                        +qcf(i,j,1))/qsl(i,j)
          END DO
        END DO
! DEPENDS ON: qsat_wat_mix
        CALL qsat_wat_mix( qsl, TL(1,1,2),                              &
     &       p_layer_centres(1,1,2), row_length*rows, l_mixing_ratio)
!
! Do nothing to top and bottom layers
        DO j = 1, rhc_rows
          DO i = 1, rhc_row_length
            rhcrit_large(i,j,1) = rhcrit(i,j,1)
            rhcrit_large(i,j,large_levels)= rhcrit(i,j,levels)
          END DO
        END DO
!
! Lsarc_do2:
        DO j = 1, rows
          DO i = 1, row_length
            p_large(i,j,1) = p_layer_centres(i,j,1)
            T_large(i,j,1) = T(i,j,1)
            q_large(i,j,1) = Q(i,j,1)
            qcl_large(i,j,1) = qcl_latest(i,j,1)
            CF_large(i,j,1) = CF(i,j,1)
            CFL_large(i,j,1) = CFL(i,j,1)
            CFF_large(i,j,1) = CFF(i,j,1)
            tlts_large(i,j,1) = tlts(i,j,1)
            qtts_large(i,j,1) = qtts(i,j,1)
            ptts_large(i,j,1) = ptts(i,j,1)
            rhts_large(i,j,1) = rhts(i,j,1)
!
       p_large(i,j,large_levels) = p_layer_centres(i,j,levels)
            T_large(i,j,large_levels) = T(i,j,levels)
            q_large(i,j,large_levels) = Q(i,j,levels)
            qcl_large(i,j,large_levels) = qcl_latest(i,j,levels)
            CF_large(i,j,large_levels) = CF(i,j,levels)
            CFL_large(i,j,large_levels) = CFL(i,j,levels)
            CFF_large(i,j,large_levels) = CFF(i,j,levels)
!           tlts_large(i,j,large_levels) = tlts(i,j,large_levels)
!           qtts_large(i,j,large_levels) = qtts(i,j,large_levels)
!           ptts_large(i,j,large_levels) = ptts(i,j,large_levels)
!           rhts_large(i,j,large_levels) = rhts(i,j,large_levels)
            tlts_large(i,j,large_levels) = tlts(i,j,levels)
            qtts_large(i,j,large_levels) = qtts(i,j,levels)
            ptts_large(i,j,large_levels) = ptts(i,j,levels)
            rhts_large(i,j,large_levels) = rhts(i,j,levels)
! Test for continuity (assumed if linked is .true.)
            qt_norm_next=(Q(i,j,2)+qcl_latest(i,j,2)                    &
     &                            +qcf(i,j,2))/qsl(i,j)
            linked(i,j,1) =                                             &
     &            (drat_thresh  >=  ABS(qt_norm(i,j) - qt_norm_next))
            qt_norm(i,j) = qt_norm_next
          END DO
        END DO ! Lsarc_do2
!
! Lsarc_do3:
        DO k = 2, (levels - 1)
          k_index = 3 + (levels_per_level * (k-2))
!
! DEPENDS ON: qsat_wat_mix
          CALL qsat_wat_mix( qsl, TL(1,1,(k+1)),                        &
     &   p_layer_centres(1,1,(k+1)), row_length*rows, l_mixing_ratio)
!
! Select associated rhcrit values
          DO j = 1, rhc_rows
            DO i = 1, rhc_row_length
              rhcrit_large(i,j,k_index-1) = rhcrit(i,j,k)
              rhcrit_large(i,j,k_index)   = rhcrit(i,j,k)
              rhcrit_large(i,j,k_index+1) = rhcrit(i,j,k)
            END DO
          END DO
!
! Lsarc_do3ij:
          DO j = 1, rows
            DO i = 1, row_length
! Test for continuity (assumed if linked = .true.)
              qt_norm_next=(Q(i,j,(k+1))+qcl_latest(i,j,(k+1))          &
     &                     +qcf(i,j,(k+1)))/ qsl(i,j)
              linked(i,j,k) =                                           &
     &            (drat_thresh  >=  ABS(qt_norm(i,j) - qt_norm_next))
              qt_norm(i,j) = qt_norm_next
!
! Select interpolated pressure levels
              delta_p = (p_layer_boundaries(i,j,(k-1)) -                &
     &                   p_layer_boundaries(i,j,k))    * inverse_level
              IF (p_layer_centres(i,j,k)  >=                            &
     &           (p_layer_boundaries(i,j,k) + delta_p)) THEN
                p_large(i,j,k_index) = p_layer_centres(i,j,k)
              ELSE
                p_large(i,j,k_index) = 0.5*(p_layer_boundaries(i,j,k) + &
     &                                 p_layer_boundaries(i,j,(k-1)))
              ENDIF
              p_large(i,j,(k_index-1)) = p_large(i,j,k_index)+delta_p
              p_large(i,j,(k_index+1)) = p_large(i,j,k_index)-delta_p
!
! Select variable values at layer centres
              T_large(i,j,k_index) = T(i,j,k)
              q_large(i,j,k_index) = Q(i,j,k)
              qcl_large(i,j,k_index) = qcl_latest(i,j,k)
              CF_large(i,j,k_index)=CF(i,j,k)
              CFL_large(i,j,k_index)=CFL(i,j,k)
              CFF_large(i,j,k_index)=CFF(i,j,k)
              tlts_large(i,j,k_index)=tlts(i,j,k)
              qtts_large(i,j,k_index)=qtts(i,j,k)
              ptts_large(i,j,k_index)=ptts(i,j,k)
              rhts_large(i,j,k_index)=rhts(i,j,k)
!
! Calculate increment in variable values, pressure interpolation
! NB: Using X_large(i,j,(k_index+1)) as store for X increments
! Lsarc_if2:
              IF ( linked(i,j,(k-1)) ) THEN
                IF ( linked(i,j,k) ) THEN
!               Interpolate from level k-1 to k+1
                  stretcher = delta_p /                                 &
     &           (p_layer_centres(i,j,k-1)-p_layer_centres(i,j,k+1))
                  T_large(i,j,(k_index+1)) = stretcher *                &
     &           (T(i,j,(k+1)) - T(i,j,(k-1)))
                  q_large(i,j,(k_index+1)) = stretcher *                &
     &           (Q(i,j,(k+1)) - Q(i,j,(k-1)))
                  qcl_large(i,j,(k_index+1)) = stretcher *              &
     &           (qcl_latest(i,j,(k+1)) - qcl_latest(i,j,(k-1)))
                  tlts_large(i,j,(k_index+1)) = stretcher *             &
     &           (tlts(i,j,(k+1)) - tlts(i,j,(k-1)))
                  qtts_large(i,j,(k_index+1)) = stretcher *             &
     &           (qtts(i,j,(k+1)) - qtts(i,j,(k-1)))
                  ptts_large(i,j,(k_index+1)) = stretcher *             &
     &           (ptts(i,j,(k+1)) - ptts(i,j,(k-1)))
                ELSE
!               Interpolate from level k-1 to k
                  stretcher = delta_p /                                 &
     &           (p_layer_centres(i,j,k-1) - p_large(i,j,k_index))
                  T_large(i,j,(k_index+1)) = stretcher *                &
     &           (T_large(i,j,k_index) - T(i,j,(k-1)))
                  q_large(i,j,(k_index+1)) = stretcher *                &
     &           (q_large(i,j,k_index) - Q(i,j,(k-1)))
                  qcl_large(i,j,(k_index+1)) = stretcher *              &
     &           (qcl_large(i,j,k_index) - qcl_latest(i,j,(k-1)))
                  tlts_large(i,j,(k_index+1)) = stretcher *             &
     &           (tlts_large(i,j,k_index) - tlts(i,j,(k-1)))
                  qtts_large(i,j,(k_index+1)) = stretcher *             &
     &           (qtts_large(i,j,k_index) - qtts(i,j,(k-1)))
                  ptts_large(i,j,(k_index+1)) = stretcher *             &
     &           (ptts_large(i,j,k_index) - ptts(i,j,(k-1)))
                ENDIF
!
              ELSE
                IF ( linked(i,j,k) ) THEN
!               Interpolate from level k to k+1
                  stretcher = delta_p /                                 &
     &            (p_large(i,j,k_index) - p_layer_centres(i,j,k+1))
                   T_large(i,j,(k_index+1)) = stretcher *               &
     &            (T(i,j,(k+1)) - T_large(i,j,k_index))
                   q_large(i,j,(k_index+1)) = stretcher *               &
     &            (Q(i,j,(k+1)) - q_large(i,j,k_index))
                   qcl_large(i,j,(k_index+1)) = stretcher *             &
     &            (qcl_latest(i,j,(k+1)) - qcl_large(i,j,k_index))
                   tlts_large(i,j,(k_index+1)) = stretcher *            &
     &            (tlts(i,j,(k+1)) - tlts_large(i,j,k_index))
                   qtts_large(i,j,(k_index+1)) = stretcher *            &
     &            (qtts(i,j,(k+1)) - qtts_large(i,j,k_index))
                   ptts_large(i,j,(k_index+1)) = stretcher *            &
     &            (ptts(i,j,(k+1)) - ptts_large(i,j,k_index))
                ELSE
!               No interpolation, freeze at level k
                  T_large(i,j,(k_index+1)) = 0.
                  q_large(i,j,(k_index+1)) = 0.
                  qcl_large(i,j,(k_index+1)) = 0.
                  tlts_large(i,j,(k_index+1)) = 0.
                  qtts_large(i,j,(k_index+1)) = 0.
                  ptts_large(i,j,(k_index+1)) = 0.
                ENDIF
!
              END IF ! Lsarc_if2
!
! Protect against q or qcl going negative (T would imply blow-up anyway)
              IF (q_large(i,j,k_index)  <                               &
     &                   (ABS(q_large(i,j,(k_index+1)))+tol_test))      &
     &                        q_large(i,j,(k_index+1)) = 0.
              IF (qcl_large(i,j,k_index)  <                             &
     &                     (ABS(qcl_large(i,j,(k_index+1)))+tol_test))  &
     &                          qcl_large(i,j,(k_index+1)) = 0.
              IF (qtts_large(i,j,k_index) <                             &
     &                    (ABS(qtts_large(i,j,(k_index+1)))+tol_test))  &
     &                         qtts_large(i,j,(k_index+1)) = 0.
!
! Select variable values at level below layer centre
              T_large(i,j,(k_index-1)) = T_large(i,j,k_index) -         &
     &                                   T_large(i,j,(k_index+1))
              q_large(i,j,(k_index-1)) = q_large(i,j,k_index) -         &
     &                                   q_large(i,j,(k_index+1))
              qcl_large(i,j,(k_index-1))=qcl_large(i,j,k_index)-        &
     &                                   qcl_large(i,j,(k_index+1))
              CF_large(i,j,(k_index-1))=CF(i,j,k)
              CFL_large(i,j,(k_index-1))=CFL(i,j,k)
              CFF_large(i,j,(k_index-1))=CFF(i,j,k)
              tlts_large(i,j,(k_index-1))=tlts_large(i,j,k_index) -     &
     &                                   tlts_large(i,j,(k_index+1))
              qtts_large(i,j,(k_index-1))=qtts_large(i,j,k_index) -     &
     &                                   qtts_large(i,j,(k_index+1))
              ptts_large(i,j,(k_index-1))=ptts_large(i,j,k_index) -     &
     &                                   ptts_large(i,j,(k_index+1))
!
! Select variable values at level above layer centre
! NB: CEASE using X_large(i,j,(k_index+1)) as store for X increments.
              T_large(i,j,(k_index+1)) = T_large(i,j,(k_index+1)) +     &
     &                                   T_large(i,j,k_index)
              q_large(i,j,(k_index+1)) = q_large(i,j,(k_index+1)) +     &
     &                                   q_large(i,j,k_index)
              qcl_large(i,j,(k_index+1))=qcl_large(i,j,(k_index+1)) +   &
     &                                   qcl_large(i,j,k_index)
              CF_large(i,j,(k_index+1))=CF(i,j,k)
              CFL_large(i,j,(k_index+1))=CFL(i,j,k)
              CFF_large(i,j,(k_index+1))=CFF(i,j,k)
              tlts_large(i,j,(k_index+1))=tlts_large(i,j,(k_index+1)) + &
     &                                    tlts_large(i,j,k_index)
              qtts_large(i,j,(k_index+1))=qtts_large(i,j,(k_index+1)) + &
     &                                    qtts_large(i,j,k_index)
              ptts_large(i,j,(k_index+1))=ptts_large(i,j,(k_index+1)) + &
     &                                    ptts_large(i,j,k_index)
!             
            END DO
          END DO ! Lsarc_do3ij
!
! calculate RH above and below layer centres
! DEPENDS ON: qsat_wat_mix
        call qsat_wat_mix(qsl,tlts_large(1,1,k_index-1),                &
     &                    ptts_large(1,1,k_index-1),                    &
     &                    row_length*rows,l_mixing_ratio)
!
          DO j = 1, rows
            DO i = 1, row_length
               rhts_large(i,j,k_index-1)=qtts_large(i,j,k_index-1)      &
     &                                     /qsl(i,j)
            END DO
          END DO
!
! DEPENDS ON: qsat_wat_mix
        call qsat_wat_mix(qsl,tlts_large(1,1,k_index+1),                &
     &                    ptts_large(1,1,k_index+1),                    &
     &                    row_length*rows,l_mixing_ratio)
!
          DO j = 1, rows
            DO i = 1, row_length
               rhts_large(i,j,k_index+1)=qtts_large(i,j,k_index+1)      &
     &                                     /qsl(i,j)
            END DO
          END DO
!
        END DO ! Lsarc_do3
!
! DEPENDS ON: pc2_initiate
          CALL PC2_INITIATE(p_large,ccb,cumulus,rhcrit_large,           &
     &      large_levels, row_length, rows, rhc_row_length,rhc_rows,    &
     &      t_large,cf_large,cfl_large,cff_large,q_large,qcl_large,     &
     &      rhts_large,l_mixing_ratio)
!
! Lsarc_do4:
        DO j = 1, rows
          DO i = 1, row_length
            T(i,j,1) = T_large(i,j,1)
            Q(i,j,1) = q_large(i,j,1)
!
            CF_area(i,j,1) = CF_large(i,j,1)
            CF(i,j,1) = CF_large(i,j,1)
!          
            qcl_latest(i,j,1) = qcl_large(i,j,1)
            CFL(i,j,1) = CFL_large(i,j,1)
!
            T(i,j,levels) = T_large(i,j,large_levels)
            Q(i,j,levels) = q_large(i,j,large_levels)
!
            CF_area(i,j,levels) = CF_large(i,j,large_levels)
            CF(i,j,levels) = CF_large(i,j,large_levels)
!
            qcl_latest(i,j,levels)=qcl_large(i,j,large_levels)
            CFL(i,j,levels) = CFL_large(i,j,large_levels)
          END DO
        END DO ! Lsarc_do4
!
! Output variables for remaining layers
! Lsarc_do5:
        DO k = 2, (levels - 1)
          k_index = 3 + (levels_per_level * (k-2))
! Lsarc_do5ij:
          DO j = 1, rows
            DO i = 1, row_length
! Area cloud fraction is maximum of sub-layer cloud fractions
              CF_area(i,j,k) =                                          &
     &        MAX( CF_large(i,j,k_index),                               &
     &             (MAX(CF_large(i,j,(k_index+1)),                      &
     &                  CF_large(i,j,(k_index-1)))) )
!
! Bulk cloud fraction is mean of sub-layer cloud fractions : strictly
! this is a pressure weighted mean being used to approximate a volume
! mean. Over the depth of a layer the difference should not be large.
              CF(i,j,k) = inverse_level *                               &
     &           ( CF_large(i,j,(k_index-1)) +                          &
     &             CF_large(i,j, k_index)    +                          &
     &             CF_large(i,j,(k_index+1)) )
!
! The pressure weighted mean of qcf is the input qcf: do not update.
!
! Qcl is the pressure weighted mean of qcl from each sub-layer.
              qcl_latest(i,j,k) = inverse_level *                       &
     &      ( qcl_large(i,j,(k_index-1)) +                              &
     &        qcl_large(i,j,k_index) + qcl_large(i,j,(k_index+1)) )
!
! Liq. cloud fraction is mean of sub-layer cloud fractions : strictly
! this is a pressure weighted mean being used to approximate a volume
! mean. Over the depth of a layer the difference should not be large.
              CFL(i,j,k) = inverse_level *                              &
     &      ( CFL_large(i,j,(k_index-1)) +                              &
     &        CFL_large(i,j,k_index)     +                              &
     &        CFL_large(i,j,(k_index+1)) )
!
! Update Q
! Update T
! Move qcl_latest into qcl.
              Q(i,j,k) = Q(i,j,k) +qcl(i,j,k) - qcl_latest(i,j,k)
              T(i,j,k) = T(i,j,k) -(qcl(i,j,k)*LCRCP) +                 &
     &                            (qcl_latest(i,j,k) * LCRCP)
              qcl(i,j,k) = qcl_latest(i,j,k)
!
            END DO
          END DO ! Lsarc_do5ij
!
        END DO ! Lsarc_do5
!
      RETURN
      END SUBROUTINE pc2_arcld
! ======================================================================
#endif
