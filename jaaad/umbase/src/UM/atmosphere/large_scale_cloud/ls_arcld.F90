#if defined(A09_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale Area (Vertical Gradient) Cloud Scheme.
! Subroutine Interface:
      SUBROUTINE LS_ARCLD(                                              &
!      Pressure related fields
     &  p_layer_centres, rhcrit, p_layer_boundaries                     &
!      Array dimensions
     &, model_levels, wet_model_levels, row_length, rows                &
     &, rhc_row_length, rhc_rows, bl_levels                             &
     &, cloud_fraction_method,overlap_ice_liquid                        &
     &, ice_fraction_method,ctt_weight,t_weight                         &
     &, qsat_fixed,sub_cld                                              &
     &, levels_per_level, large_levels                                  &
!      Switch on area cloud calculation and select which to use
     &, L_AREA_CLOUD, L_ACF_Cusack, L_ACF_Brooks, L_eacf                &
!      Needed for LS_ACF_Brooks
     &, halo_i, halo_j, off_x, off_y                                    &
     &, delta_lambda, delta_phi                                         &
     &, r_theta_levels, FV_cos_theta_latitude                           &
!      Convection diagnosis information (only used for A05_4A)
     &, ntml, cumulus, L_mixing_ratio                                   &
!      Prognostic Fields
     &, qcf_latest                                                      &
     &, T_latest, q_latest, qcl_latest                                  &
!      Various cloud fractions
     &, area_cloud_fraction, bulk_cloud_fraction                        &
     &, cloud_fraction_liquid, cloud_fraction_frozen                    &
     &, error_code, me)
!
      IMPLICIT NONE
!
! Purpose:
!   This subroutine calculates liquid and ice cloud fractional cover
!   for use with the enhanced precipitation microphysics scheme. It
!   also returns area and bulk cloud fractions for use in radiation.
!
! Method:
!   Statistical cloud scheme separates input moisture into specific
!   humidity and cloud liquid water. Temperature calculated from liquid
!   water temperature. Cloud fractions calculated from statistical
!   relation between cloud fraction and cloud liquid/ice water content.
!   Area cloud fraction calculated by subdividing layers, calculating
!   cloud on each sublayer and taking maximum (use mean for bulk cloud).
!
! Current Owner of Code: A. C. Bushell
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP No.29
!
!  Global Variables:----------------------------------------------------
#include "c_r_cp.h"
#include "c_lheat.h"
!
!  Subroutine Arguments:------------------------------------------------
!
! arguments with intent in. ie: input variables.
!
      INTEGER                                                           &
                        !, INTENT(IN)
     &  model_levels                                                    &
!       No. of levels being processed by model.
     &, wet_model_levels                                                &
!       No. of levels being processed by cloud scheme.
     &, bl_levels                                                       &
!       No. of boundary layer levels
     &, row_length, rows                                                &
!       Horizontal dimensions being processed by cloud scheme.
     &, rhc_row_length, rhc_rows                                        &
!       Horizontal dimensions of rh_crit diagnostic.
     &, levels_per_level                                                &
!       No. of sub-levels being processed by area cloud scheme.
!       Want an odd number of sublevels per level.
!       NB: levels_per_level = 3 is currently hardwired in the do loops
     &, large_levels                                                    &
!       Total no. of sub-levels being processed by cloud scheme.
!       Currently ((wet_model_levels - 2)*levels_per_level) + 2
     &, me                                                              &
!       My processor number
!
     &, ntml(row_length, rows)                                          &
                                  !  Height of diagnosed BL top
     &, halo_i                                                          &
                          ! Size of halo in i direction
     &, halo_j                                                          &
                          ! Size of halo in j direction
     &, off_x                                                           &
                          ! Size of small halo in i direction
     &, off_y
                          ! Size of small halo in j direction

      LOGICAL                                                           &
                        !, INTENT(IN)
     & L_AREA_CLOUD                                                     &
                           ! true if using area cloud fraction (ACF)
     &, L_ACF_Cusack                                                    &
                           ! ... and selected Cusack and PC2 off
     &, L_ACF_Brooks                                                    &
                           ! ... and selected Brooks
     &,L_eacf                                                           &
                       ! true if using empirically adjusted
                       ! cloud fraction
     &,L_mixing_ratio                                                    
                       ! true if using mixing rations

      LOGICAL                                                           &
                        !, INTENT(IN)
     &  cumulus(row_length, rows)   ! Logical indicator of convection

      REAL                                                              &
                        !, INTENT(IN)
     & qcf_latest(row_length,rows,wet_model_levels)                     &
!       Cloud ice content at processed levels (kg water per kg air).
     &,p_layer_centres(row_length,rows,0:model_levels)                  &
!       pressure at all points, on theta levels (Pa).
     &,rhcrit(rhc_row_length,rhc_rows,wet_model_levels)                 &
!       Critical relative humidity.  See the the paragraph incorporating
!       eqs P292.11 to P292.14; the values need to be tuned for a given
!       set of levels.
     &,p_layer_boundaries(row_length,rows,0:model_levels)               &
!       pressure at all points, on u,v levels (Pa).
     &, r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
                         ! height of theta levels (from centre of earth)
     &, FV_cos_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)                      &
                         ! Finite volume cos(lat)
     &, delta_lambda                                                    &
                         ! EW (x) grid spacing in radians
     &, delta_phi        ! NS (y) grid spacing in radians

!
      Integer                                                           &
                        !, Intent(IN)
     &  cloud_fraction_method                                           &
                               ! Method for calculating
                               ! total cloud fraction
     &, ice_fraction_method    ! Method for calculating ice cloud frac.

      Real                                                              &
                        !, Intent(IN)
     &  overlap_ice_liquid                                              &
                               ! Overlap between ice and liquid phases
     &, ctt_weight                                                      &
                               ! Weighting of cloud top temperature
     &, t_weight                                                        &
                               ! Weighting of local temperature
     &, qsat_fixed                                                      &
                               ! Fixed value of saturation humidity
     &, sub_cld                ! Scaling parameter
!
! arguments with intent in/out. ie: input variables changed on output.
!
      REAL                                                              &
                        !, INTENT(INOUT)
     & q_latest(row_length,rows,wet_model_levels)                       &
!       On input : Vapour + liquid water content (QW) (kg per kg air).
!       On output: Specific humidity at processed levels
!                   (kg water per kg air).
     &,T_latest(row_length,rows,model_levels)
!       On input : Liquid water temperature (TL) (K).
!       On output: Temperature at processed levels (K).
!
! arguments with intent out. ie: output variables.
!
!     Error Status:
      INTEGER error_code  !, INTENT(OUT)  0 if OK; 1 if bad arguments.
!
      REAL                                                              &
                          !, INTENT(OUT)
     & qcl_latest(row_length,rows,wet_model_levels)                     &
!       Cloud liquid content at processed levels (kg water per kg air).
     &,area_cloud_fraction(row_length,rows,wet_model_levels)            &
!       Area cloud fraction at processed levels (decimal fraction).
     &,bulk_cloud_fraction(row_length,rows,wet_model_levels)            &
!       Cloud fraction at processed levels (decimal fraction).
     &,cloud_fraction_liquid(row_length,rows,wet_model_levels)          &
!       Liquid cloud fraction at processed levels (decimal fraction).
     &,cloud_fraction_frozen(row_length,rows,wet_model_levels)
!       Frozen cloud fraction at processed levels (decimal fraction).
!
!  Local parameters and other physical constants------------------------
      REAL LCRCP                  ! Derived parameters.
      PARAMETER (LCRCP=LC/CP)     ! Lat ht of condensation/Cp.
      REAL drat_thresh            ! Test for continuity of sub-levels
      REAL tol_test               ! Tolerance for non-zero humidities
      PARAMETER (drat_thresh=3.0E-1, tol_test=1.0E-11)
!
!  Local scalars--------------------------------------------------------
!
!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
      REAL                                                              &
     &  qt_norm_next                                                    &
                       ! Temporary space for qT_norm
     &, stretcher                                                       &
                       !
     &, inverse_level                                                   &
                       ! Set to (1. / levels_per_level)
     &, delta_p        ! Layer pressure thickness * inverse_level
!
!  (b) Others.
      INTEGER i,j,k      ! Loop counters: k - vertical level index.
!                                       i,j - horizontal field index.
      INTEGER k_index    ! Extra loop counter for large arrays.
!
!  Local dynamic arrays-------------------------------------------------
!    11 blocks of real workspace are required.
      REAL                                                              &
     &  qsl(row_length,rows)                                            &
!        Saturated specific humidity for temp TL or T.
     &, qt_norm(row_length,rows)                                        &
!        Total water content normalized to qSAT_WAT.
     &, rhcrit_large(rhc_row_length,rhc_rows,large_levels)              &
!
     &, p_large(row_length,rows,large_levels)                           &
!
     &, T_large(row_length,rows,large_levels)                           &
!
     &, q_large(row_length,rows,large_levels)                           &
!
     &, qcf_large(row_length,rows,large_levels)                         &
!
     &, cloud_fraction_large(row_length,rows,large_levels)              &
!
     &, qcl_large(row_length,rows,large_levels)                         &
!
     &, cloud_fraction_liquid_large(row_length,rows,large_levels)       &
!
     &, cloud_fraction_frozen_large(row_length,rows,large_levels)
!
      INTEGER                                                           &
     &  ntml_large(row_length,rows)

      LOGICAL                                                           &
     & linked(row_length,rows,wet_model_levels)
!       True for sub-layers that have similar supersaturation properties
!
!  External subroutine calls: ------------------------------------------
      EXTERNAL QSAT_WAT,LS_CLD,LS_ACF_BROOKS
!- End of Header
! ----------------------------------------------------------------------
!  Check input arguments for potential over-writing problems.
! ----------------------------------------------------------------------
      error_code=0
!
! ==Main Block==--------------------------------------------------------
! Subroutine structure :
! Loop round levels to be processed.
! ----------------------------------------------------------------------
! Lsarc_if1:
      IF (.NOT. L_AREA_CLOUD) THEN
!     As before
!
! DEPENDS ON: ls_cld
        Call ls_cld( p_layer_centres(1,1,1), rhcrit,                    &
     &               wet_model_levels, bl_levels, row_length, rows,     &
     &               rhc_row_length, rhc_rows,                          &
     &               cloud_fraction_method,overlap_ice_liquid,          &
     &               ice_fraction_method,ctt_weight,t_weight,           &
     &               qsat_fixed,sub_cld,ntml,cumulus,L_eacf,            &
     &               L_mixing_ratio,T_latest, bulk_cloud_fraction,      &
     &               q_latest, qcf_latest, qcl_latest,                  &
     &               cloud_fraction_liquid,                             &
     &               cloud_fraction_frozen, error_code )
!
! Lsarc_do1:
        DO k = 1, wet_model_levels
          DO j = 1, rows
            DO i = 1, row_length
              area_cloud_fraction(i,j,k)=bulk_cloud_fraction(i,j,k)
            END DO
          END DO
        END DO ! Lsarc_do1
!
      ELSE ! L_area_cloud
        IF (L_ACF_CUSACK) THEN
!       Vertical gradient area cloud option
!
        inverse_level = 1. / levels_per_level
!
! Test for continuity between adjacent layers based on supersaturation
! (qt - qsl) / qsl : as we take differences the - qsl term drops out.
! DEPENDS ON: qsat_wat_mix
        CALL qsat_wat_mix( qsl, T_latest(1,1,1),                        &
     &       p_layer_centres(1,1,1), row_length*rows, l_mixing_ratio)
        DO j = 1, rows
          DO i = 1, row_length
            qt_norm(i,j) =(q_latest(i,j,1)+qcf_latest(i,j,1))/qsl(i,j)
          END DO
        END DO
! DEPENDS ON: qsat_wat_mix
        CALL qsat_wat_mix( qsl, T_latest(1,1,2),                        &
     &       p_layer_centres(1,1,2), row_length*rows, l_mixing_ratio)
!
! Do nothing to top and bottom layers
        DO j = 1, rhc_rows
          DO i = 1, rhc_row_length
            rhcrit_large(i,j,1) = rhcrit(i,j,1)
            rhcrit_large(i,j,large_levels)= rhcrit(i,j,wet_model_levels)
          END DO
        END DO
!
! Lsarc_do2:
        DO j = 1, rows
          DO i = 1, row_length
            p_large(i,j,1) = p_layer_centres(i,j,1)
            T_large(i,j,1) = T_latest(i,j,1)
            q_large(i,j,1) = q_latest(i,j,1)
            qcf_large(i,j,1) = qcf_latest(i,j,1)
!
       p_large(i,j,large_levels) = p_layer_centres(i,j,wet_model_levels)
            T_large(i,j,large_levels) = T_latest(i,j,wet_model_levels)
            q_large(i,j,large_levels) = q_latest(i,j,wet_model_levels)
          qcf_large(i,j,large_levels) = qcf_latest(i,j,wet_model_levels)
! Test for continuity (assumed if linked is .true.)
            qt_norm_next=(q_latest(i,j,2)+qcf_latest(i,j,2))/qsl(i,j)
            linked(i,j,1) =                                             &
     &            (drat_thresh  >=  ABS(qt_norm(i,j) - qt_norm_next))
            qt_norm(i,j) = qt_norm_next
          END DO
        END DO ! Lsarc_do2
!
! Lsarc_do3:
        DO k = 2, (wet_model_levels - 1)
          k_index = 3 + (levels_per_level * (k-2))
!
! DEPENDS ON: qsat_wat_mix
          CALL qsat_wat_mix( qsl, T_latest(1,1,(k+1)),                  &
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
              qt_norm_next=(q_latest(i,j,(k+1))+qcf_latest(i,j,(k+1)))  &
     &                     / qsl(i,j)
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
              T_large(i,j,k_index) = T_latest(i,j,k)
              q_large(i,j,k_index) = q_latest(i,j,k)
              qcf_large(i,j,k_index) = qcf_latest(i,j,k)
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
     &           (T_latest(i,j,(k+1)) - T_latest(i,j,(k-1)))
                  q_large(i,j,(k_index+1)) = stretcher *                &
     &           (q_latest(i,j,(k+1)) - q_latest(i,j,(k-1)))
                  qcf_large(i,j,(k_index+1)) = stretcher *              &
     &           (qcf_latest(i,j,(k+1)) - qcf_latest(i,j,(k-1)))
                ELSE
!               Interpolate from level k-1 to k
                  stretcher = delta_p /                                 &
     &           (p_layer_centres(i,j,k-1) - p_large(i,j,k_index))
                  T_large(i,j,(k_index+1)) = stretcher *                &
     &           (T_large(i,j,k_index) - T_latest(i,j,(k-1)))
                  q_large(i,j,(k_index+1)) = stretcher *                &
     &           (q_large(i,j,k_index) - q_latest(i,j,(k-1)))
                  qcf_large(i,j,(k_index+1)) = stretcher *              &
     &           (qcf_large(i,j,k_index) - qcf_latest(i,j,(k-1)))
                ENDIF
!
              ELSE
                IF ( linked(i,j,k) ) THEN
!               Interpolate from level k to k+1
                  stretcher = delta_p /                                 &
     &            (p_large(i,j,k_index) - p_layer_centres(i,j,k+1))
                   T_large(i,j,(k_index+1)) = stretcher *               &
     &            (T_latest(i,j,(k+1)) - T_large(i,j,k_index))
                   q_large(i,j,(k_index+1)) = stretcher *               &
     &            (q_latest(i,j,(k+1)) - q_large(i,j,k_index))
                   qcf_large(i,j,(k_index+1)) = stretcher *             &
     &            (qcf_latest(i,j,(k+1)) - qcf_large(i,j,k_index))
                ELSE
!               No interpolation, freeze at level k
                  T_large(i,j,(k_index+1)) = 0.
                  q_large(i,j,(k_index+1)) = 0.
                  qcf_large(i,j,(k_index+1)) = 0.
                ENDIF
!
              END IF ! Lsarc_if2
!
! Protect against q or qcf going negative (T would imply blow-up anyway)
              IF (q_large(i,j,k_index)  <                               &
     &                   (ABS(q_large(i,j,(k_index+1)))+tol_test))      &
     &                        q_large(i,j,(k_index+1)) = 0.
              IF (qcf_large(i,j,k_index)  <                             &
     &                     (ABS(qcf_large(i,j,(k_index+1)))+tol_test))  &
     &                          qcf_large(i,j,(k_index+1)) = 0.
!
! Select variable values at level below layer centre
              T_large(i,j,(k_index-1)) = T_large(i,j,k_index) -         &
     &                                   T_large(i,j,(k_index+1))
              q_large(i,j,(k_index-1)) = q_large(i,j,k_index) -         &
     &                                   q_large(i,j,(k_index+1))
              qcf_large(i,j,(k_index-1))=qcf_large(i,j,k_index)-        &
     &                                   qcf_large(i,j,(k_index+1))
!
! Select variable values at level above layer centre
! NB: CEASE using X_large(i,j,(k_index+1)) as store for X increments.
              T_large(i,j,(k_index+1)) = T_large(i,j,(k_index+1)) +     &
     &                                   T_large(i,j,k_index)
              q_large(i,j,(k_index+1)) = q_large(i,j,(k_index+1)) +     &
     &                                   q_large(i,j,k_index)
              qcf_large(i,j,(k_index+1))=qcf_large(i,j,(k_index+1)) +   &
     &                                   qcf_large(i,j,k_index)
            END DO
          END DO ! Lsarc_do3ij
!
        END DO ! Lsarc_do3
!
!
! Create an array of NTML values adjusted to the large levels

        DO  j = 1, rows
          DO i = 1, row_length
            ntml_large(i,j) = 3+(levels_per_level*(ntml(i,j)-1))
          END DO
        END DO

! DEPENDS ON: ls_cld
        CALL ls_cld( p_large, rhcrit_large,                             &
     &               large_levels, bl_levels, row_length, rows,         &
     &               rhc_row_length, rhc_rows,                          &
     &               cloud_fraction_method,overlap_ice_liquid,          &
     &               ice_fraction_method,ctt_weight,t_weight,           &
     &               qsat_fixed,sub_cld,                                &
     &               ntml_large, cumulus, L_eacf, l_mixing_ratio,       &
     &               T_large, cloud_fraction_large,                     &
     &               q_large, qcf_large, qcl_large,                     &
     &               cloud_fraction_liquid_large,                       &
     &               cloud_fraction_frozen_large, error_code )

! Lsarc_do4:
        DO j = 1, rows
          DO i = 1, row_length
            T_latest(i,j,1) = T_large(i,j,1)
            q_latest(i,j,1) = q_large(i,j,1)

            area_cloud_fraction(i,j,1) = cloud_fraction_large(i,j,1)
            bulk_cloud_fraction(i,j,1) = cloud_fraction_large(i,j,1)
          
            qcl_latest(i,j,1) = qcl_large(i,j,1)
            cloud_fraction_liquid(i,j,1) =                              &
     &        cloud_fraction_liquid_large(i,j,1)
            cloud_fraction_frozen(i,j,1) =                              &
     &          cloud_fraction_frozen_large(i,j,1)

            T_latest(i,j,wet_model_levels) = T_large(i,j,large_levels)
            q_latest(i,j,wet_model_levels) = q_large(i,j,large_levels)

            area_cloud_fraction(i,j,wet_model_levels) =                 &
     &          cloud_fraction_large(i,j,large_levels)
            bulk_cloud_fraction(i,j,wet_model_levels) =                 &
     &          cloud_fraction_large(i,j,large_levels)

            qcl_latest(i,j,wet_model_levels)=qcl_large(i,j,large_levels)
            cloud_fraction_liquid(i,j,wet_model_levels) =               &
     &          cloud_fraction_liquid_large(i,j,large_levels)
            cloud_fraction_frozen(i,j,wet_model_levels) =               &
     &          cloud_fraction_frozen_large(i,j,large_levels)
          END DO
        END DO ! Lsarc_do4
!
! Output variables for remaining layers
! Lsarc_do5:
        DO k = 2, (wet_model_levels - 1)
          k_index = 3 + (levels_per_level * (k-2))
! Lsarc_do5ij:
          DO j = 1, rows
            DO i = 1, row_length
! Area cloud fraction is maximum of sub-layer cloud fractions
              area_cloud_fraction(i,j,k) =                              &
     &        MAX( cloud_fraction_large(i,j,k_index),                   &
     &             (MAX(cloud_fraction_large(i,j,(k_index+1)),          &
     &                  cloud_fraction_large(i,j,(k_index-1)))) )
!
! Bulk cloud fraction is mean of sub-layer cloud fractions : strictly
! this is a pressure weighted mean being used to approximate a volume
! mean. Over the depth of a layer the difference should not be large.
              bulk_cloud_fraction(i,j,k) = inverse_level *              &
     &           ( cloud_fraction_large(i,j,(k_index-1)) +              &
     &             cloud_fraction_large(i,j, k_index)    +              &
     &             cloud_fraction_large(i,j,(k_index+1)) )
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
              cloud_fraction_liquid(i,j,k) = inverse_level *            &
     &      ( cloud_fraction_liquid_large(i,j,(k_index-1)) +            &
     &        cloud_fraction_liquid_large(i,j,k_index)     +            &
     &        cloud_fraction_liquid_large(i,j,(k_index+1)) )
!
! Froz cloud fraction is mean of sub-layer cloud fractions : strictly
! this is a pressure weighted mean being used to approximate a volume
! mean. Over the depth of a layer the difference should not be large.
              cloud_fraction_frozen(i,j,k) = inverse_level *            &
     &      ( cloud_fraction_frozen_large(i,j,(k_index-1)) +            &
     &        cloud_fraction_frozen_large(i,j,k_index)     +            &
     &        cloud_fraction_frozen_large(i,j,(k_index+1)) )
!
! Transform q_latest from qT(vapour + liquid) to specific humidity.
! Transform T_latest from TL(vapour + liquid) to temperature.
              q_latest(i,j,k) = q_latest(i,j,k) - qcl_latest(i,j,k)
              T_latest(i,j,k) = T_latest(i,j,k) +                       &
     &                            (qcl_latest(i,j,k) * LCRCP)
            END DO
          END DO ! Lsarc_do5ij
!
        END DO ! Lsarc_do5
!
        ELSE IF (L_ACF_Brooks) THEN  ! L_ACF_Cusack or Brooks
!
!       As before, update variables that would have been updated
!       without area cloud fraction on
!
! DEPENDS ON: ls_cld
          Call ls_cld( p_layer_centres(1,1,1), rhcrit,                  &
     &               wet_model_levels, bl_levels, row_length, rows,     &
     &               rhc_row_length, rhc_rows,                          &
     &               cloud_fraction_method,overlap_ice_liquid,          &
     &               ice_fraction_method,ctt_weight,t_weight,           &
     &               qsat_fixed,sub_cld,                                &
     &               ntml, cumulus, L_eacf, L_mixing_ratio,             &
     &               T_latest, bulk_cloud_fraction,                     &
     &               q_latest, qcf_latest, qcl_latest,                  &
     &               cloud_fraction_liquid,                             &
     &               cloud_fraction_frozen, error_code )
!
!      Calculate the area cloud fraction

! DEPENDS ON: ls_acf_brooks
          Call LS_ACF_Brooks (                                          &
     &             halo_i, halo_j, off_x, off_y                         &
     &,            row_length, rows, model_levels, wet_model_levels     &
     &,            r_theta_levels, delta_lambda, delta_phi              &
     &,            FV_cos_theta_latitude                                &
     &,            bulk_cloud_fraction, cloud_fraction_liquid           &
     &,            cloud_fraction_frozen, cumulus                       &
     &,            area_cloud_fraction )

!
        END IF ! L_ACF_Brooks
      END IF ! Lsarc_if1  L_area_cloud
!
      RETURN
      END SUBROUTINE LS_ARCLD
! ======================================================================
#endif
