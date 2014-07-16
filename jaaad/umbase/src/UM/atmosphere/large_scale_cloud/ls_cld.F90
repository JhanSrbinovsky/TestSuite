#if defined(A09_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale Cloud Scheme.
! Subroutine Interface:
      SUBROUTINE LS_CLD(                                                &
!      Pressure related fields
     & p_theta_levels, RHCRIT                                           &
!      Array dimensions
     &,LEVELS, BL_LEVELS, row_length, rows                              &
     &,rhc_row_length,rhc_rows                                          &
     &,cloud_fraction_method,overlap_ice_liquid                         &
     &,ice_fraction_method,ctt_weight,t_weight,qsat_fixed,sub_cld       &
!      From convection diagnosis (only used if A05_4A)
     &,ntml, cumulus, L_eacf, L_mixing_ratio                            &
!      Prognostic Fields
     &,T, CF, Q, QCF, QCL                                               &
!      Liquid and frozen ice cloud fractions
     &,CFL, CFF                                                         &
     &,ERROR)

      Use cv_run_mod, Only:                                             &
          l_conv4a

      IMPLICIT NONE
!
! Purpose:
!   This subroutine calculates liquid and ice cloud fractional cover
!   for use with the enhanced precipitation microphysics scheme.
!
! Method:
!   Statistical cloud scheme separates input moisture into specific
!   humidity and cloud liquid water. Temperature calculated from liquid
!   water temperature. Cloud fractions calculated from statistical
!   relation between cloud fraction and cloud liquid/ice water content.
!
! Current Owner of Code: A. C. Bushell
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP No. 29
!
!  Global Variables:----------------------------------------------------
#include "c_mdi.h"
#include "c_pi.h"
#include "c_cldsgs.h"
!
!  Subroutine Arguments:------------------------------------------------
      INTEGER                                                           &
                        !, INTENT(IN)
     & LEVELS                                                           &
!       No. of levels being processed.
     &,BL_LEVELS                                                        &
!       No. of boundary layer levels
     &, row_length,rows                                                 &
     &, rhc_row_length,rhc_rows
!
      REAL                                                              &
                        !, INTENT(IN)
     & QCF(row_length,rows,LEVELS)                                      &
!       Cloud ice content at processed levels (kg water per kg air).
     &,p_theta_levels(row_length,rows,levels)                           &
!       pressure at all points (Pa).
     &,RHCRIT(rhc_row_length,rhc_rows,LEVELS)
!       Critical relative humidity.  See the the paragraph incorporating
!       eqs P292.11 to P292.14; the values need to be tuned for the give
!       set of levels.
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
      INTEGER                                                           &
     & NTML(row_length,rows)     ! IN Height of diagnosed BL top

      LOGICAL                                                           &
     & L_eacf                                                           &
                       ! IN true if using empirically adjusted cloud
                       !    fraction parametrization
     &,L_mixing_ratio  ! IN true if using mixing ratios

      LOGICAL                                                           &
     & CUMULUS(row_length,rows)  ! IN Logical indicator of convection

      REAL                                                              &
                        !, INTENT(INOUT)
     & Q(row_length,rows,LEVELS)                                        &
!       On input : Total water content (QW) (kg per kg air).
!       On output: Specific humidity at processed levels
!                   (kg water per kg air).
     &,T(row_length,rows,LEVELS)
!       On input : Liquid/frozen water temperature (TL) (K).
!       On output: Temperature at processed levels (K).
!
      REAL                                                              &
                        !, INTENT(OUT)
     & CF(row_length,rows,LEVELS)                                       &
!       Cloud fraction at processed levels (decimal fraction).
     &,QCL(row_length,rows,LEVELS)                                      &
!       Cloud liquid water content at processed levels (kg per kg air).
     &,CFL(row_length,rows,LEVELS)                                      &
!       Liquid cloud fraction at processed levels (decimal fraction).
     &,CFF(row_length,rows,LEVELS)
!       Frozen cloud fraction at processed levels (decimal fraction).
!
!     Error Status:
      INTEGER ERROR     !, INTENT(OUT)  0 if OK; 1 if bad arguments.
!
!  Local parameters and other physical constants------------------------
      REAL ROOTWO       ! Sqrt(2.)
      REAL SUBGRID      ! Subgrid parameter in ice cloud calculation
!
!  Local scalars--------------------------------------------------------
!
!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
      REAL                                                              &
     & PHIQCF                                                           &
                        ! Arc-cosine term in Cloud ice fraction calc.
     &,COSQCF                                                           &
                        ! Cosine term in Cloud ice fraction calc.
     &,OVERLAP_MAX                                                      &
                        ! Maximum possible overlap
     &,OVERLAP_MIN                                                      &
                        ! Minimum possible overlap
     &,OVERLAP_RANDOM                                                   &
                        ! Random overlap
     &,TEMP0                                                            &
     &,TEMP1                                                            &
     &,TEMP2                                                            &
                        ! Temporaries for combinations of the
     &,QN_IMP                                                           &
     &,QN_ADJ
!                       ! overlap parameters
!
!  (b) Others.
      INTEGER K,I,J       ! Loop counters: K - vertical level index.
!                                        I,J - horizontal field indices.
!
      INTEGER QC_POINTS                                                 &
                        ! No. points with non-zero cloud
     &,       MULTRHC   ! Zero if (rhc_row_length*rhc_rows) le 1, else 1
!
!  Local dynamic arrays-------------------------------------------------
!    6 blocks of real workspace are required.
      REAL                                                              &
     & QCFRBS(row_length,rows)                                          &
                               ! qCF / bs
     &,QSL(row_length,rows)                                             &
!       Saturated specific humidity for temp TL or T.
     &,QSL_CTT(row_length,rows)                                         &
!       Saturated specific humidity wrt liquid at cloud top temperature
     &,QN(row_length,rows)                                              &
!       Cloud water normalised with BS.
     &,GRID_QC(row_length,rows,LEVELS)                                  &
!       Gridbox mean saturation excess at processed levels
!        (kg per kg air). Set to RMDI when cloud is absent.
     &,BS(row_length,rows,LEVELS)                                       &
!       Maximum moisture fluctuation /6*sigma at processed levels
!        (kg per kg air). Set to RMDI when cloud is absent.
     &,CTT(row_length,rows)
!       Ice cloud top temperature (K) - as coded it is really TL
      LOGICAL                                                           &
     & LQC(row_length,rows)      ! True for points with non-zero cloud
      INTEGER                                                           &
     & INDEX(row_length*rows,2)                                         &
                                 ! Index for points with non-zero cloud
     &,LLWIC(row_length,rows)
!       Last Level With Ice Cloud
      REAL RHCRITX              ! scalar copy of RHCRIT(I,J,K)
!
!  External subroutine calls: ------------------------------------------
      EXTERNAL QSAT,QSAT_WAT,LS_CLD_C
!- End of Header
! ----------------------------------------------------------------------
!  Check input arguments for potential over-writing problems.
! ----------------------------------------------------------------------
      ERROR=0
!
      IF ( (rhc_row_length * rhc_rows)  >   1) THEN
        MULTRHC = 1
      ELSE
        MULTRHC = 0
      END IF
!
! ==Main Block==--------------------------------------------------------
! Subroutine structure :
! Loop round levels to be processed.
! ----------------------------------------------------------------------
! Initialize cloud-top-temperature and last-level-with-ice-cloud arrays
!CDIR COLLAPSE
      DO j=1,rows
        DO i=1,row_length
          CTT(I,J)=0.0
          LLWIC(I,J)=0
        END DO
      END DO
! Levels_do1:
      DO K=LEVELS,1,-1
!
! ----------------------------------------------------------------------
! 1. Calculate QSAT at liquid/ice water temperature, TL, and initialize
!    cloud water, sub-grid distribution and fraction arrays.
!    This requires a preliminary calculation of the pressure.
!    NB: On entry to the subroutine 'T' is TL and 'Q' is QW.
! ----------------------------------------------------------------------
! Rows_do1:
!CDIR COLLAPSE
        DO j=1,rows
! Row_length_do1:
          DO I=1,row_length
          QCL(I,j,K) = 0.0
          CFL(I,j,K) = 0.0
          GRID_QC(I,j,K) = RMDI
          BS(I,j,K) = RMDI
          END DO ! Row_length_do1
        END DO ! Rows_do1
!
! DEPENDS ON: qsat_wat_mix
        CALL qsat_wat_mix(QSL,T(1,1,K),P_theta_levels(1,1,k),           &
     &                row_length*rows,l_mixing_ratio)
!
! Rows_do2:
        DO j=1,rows
! Row_length_do2:
          DO I=1,row_length
             IF (MULTRHC==1) THEN
                RHCRITX = RHCRIT(I,J,K)
             ELSE
                RHCRITX = RHCRIT(1,1,K)
             END IF

! Omit CUMULUS points below (and including) NTML+1

           IF ( .NOT. L_conv4a .OR. (L_conv4a .AND.                     &
     &        (.NOT. CUMULUS(I,j) .OR. ( CUMULUS(I,j)                   &
     &        .AND. (K  >   NTML(I,j)+1) ))) ) THEN

! Rhcrit_if:
            IF (RHCRITX  <   1.) THEN
! ----------------------------------------------------------------------
! 2. Calculate the quantity QN = QC/BS = (QW/QSL-1)/(1-RHcrit)
!    if RHcrit is less than 1
! ----------------------------------------------------------------------
!
              QN(I,j) = (Q(I,j,K) / QSL(I,j) - 1.) /                    &
     &                  (1. - RHCRITX)
!
! ----------------------------------------------------------------------
! 3. Set logical variable for cloud, LQC, for the case RHcrit < 1;
!    where QN > -1, i.e. qW/qSAT(TL,P) > RHcrit, there is cloud
! ----------------------------------------------------------------------
!
            LQC(I,j) = (QN(I,j)  >   -1.)
          ELSE
! ----------------------------------------------------------------------
! 2.a Calculate QN = QW - QSL if RHcrit equals 1
! ----------------------------------------------------------------------
!
            QN(I,j) = Q(I,j,K) - QSL(I,j)
!
! ----------------------------------------------------------------------
! 3.a Set logical variable for cloud, LQC, for the case RHcrit = 1;
!     where QN > 0, i.e. qW > qSAT(TL,P), there is cloud
! ----------------------------------------------------------------------
!
            LQC(I,j) = (QN(I,j)  >   0.)
          END IF ! Rhcrit_if
           ELSEIF (L_conv4a) THEN
            LQC(I,j) = .FALSE.
           END IF  ! Test on CUMULUS and NTML for A05_4A only
          END DO ! Row_length_do2
        END DO ! Rows_do2
!
! ----------------------------------------------------------------------
! 4. Form index of points where non-zero liquid cloud fraction
! ----------------------------------------------------------------------
!
        QC_POINTS=0
! Rows_do3:
        DO j=1,rows
! Row_length_do3:
         DO I=1,row_length
          IF (LQC(I,j)) THEN
            QC_POINTS = QC_POINTS + 1
            INDEX(QC_POINTS,1) = I
            INDEX(QC_POINTS,2) = j
          END IF
          END DO ! Row_length_do3
        END DO ! Rows_do3
!
! ----------------------------------------------------------------------
! 5. Call LS_CLD_C to calculate cloud water content, specific humidity,
!                  water cloud fraction and determine temperature.
! ----------------------------------------------------------------------
! Qc_points_if:
        IF (QC_POINTS  >   0) THEN
! DEPENDS ON: ls_cld_c
            CALL LS_CLD_C(P_theta_levels(1,1,K),RHCRIT(1,1,K),QSL,QN    &
     &     ,Q(1,1,K),T(1,1,K)                                           &
     &                 ,QCL(1,1,K),CFL(1,1,K),GRID_QC(1,1,K),BS(1,1,K)  &
     &     ,INDEX,QC_POINTS,row_length,rows,rhc_row_length,rhc_rows     &
     &     ,BL_LEVELS,K, L_eacf, l_mixing_ratio)
        END IF ! Qc_points_if
!
! ----------------------------------------------------------------------
! 6. Calculate cloud fractions for ice clouds.
!    THIS IS STILL HIGHLY EXPERIMENTAL.
!    Begin by calculating Qsat_wat(T,P), at Temp. T, for estimate of bs.
! ----------------------------------------------------------------------
        DO j=1,rows
          DO i=1,row_length
! Check for last level with cloud and update cloud top temperature
            IF (LLWIC(I,J)  /=  K+1) THEN
              CTT(I,J)=T(I,J,K)
            END IF
! Check for significant ice content and update last level with ice cloud
            IF (QCF(I,J,K)  >   QCFMIN) THEN
              LLWIC(I,J)=K
            ENDIF
          END DO
        END DO
! DEPENDS ON: qsat_wat_mix
        CALL qsat_wat_mix(QSL,T(1,1,K),P_theta_levels(1,1,k),           &
     &                row_length*rows,l_mixing_ratio)
        ROOTWO = SQRT(2.)
        IF (ICE_FRACTION_METHOD  ==  2) THEN
! Use cloud top temperature and a fixed qsat to give QCFRBS
! DEPENDS ON: qsat_wat_mix
          CALL qsat_wat_mix(QSL_CTT,CTT,p_theta_levels(1,1,k),          &
     &                  row_length*rows,l_mixing_ratio)
#if defined(VECTLIB)
          CALL POWR_V(row_length*rows,QSL_CTT,CTT_WEIGHT,QSL_CTT)
          CALL POWR_V(row_length*rows,QSL,T_WEIGHT,QSL)
#endif
          SUBGRID = SUB_CLD ** (1.0-T_WEIGHT)                           &
     &              / QSAT_FIXED ** (1.0-T_WEIGHT-CTT_WEIGHT)
        END IF ! ice_fraction_method eq 2
!
! Rows_do4:
        DO j=1,rows
! Row_length_do4:
         DO I=1,row_length
            IF ( MULTRHC== 1) THEN
               RHCRITX = RHCRIT(I,J,K)
            ELSE
               RHCRITX = RHCRIT(1,1,K)
            ENDIF
! ----------------------------------------------------------------------
! 6a Calculate qCF/bs.
! ----------------------------------------------------------------------
! Rhcrit_if2:
         IF (RHCRITX  <   1.) THEN
!
          IF (ICE_FRACTION_METHOD  ==  1) THEN
            QCFRBS(i,j)=  QCF(I,j,K) / ((1.-RHCRITX) * QSL(I,j))
          ELSEIF (ICE_FRACTION_METHOD  ==  2) THEN
#if defined(VECTLIB)
            QCFRBS(i,j) = SUBGRID * QCF(I,J,K) / ((1.-RHCRITX)          &
     &               * QSL_CTT(I,J)*QSL(I,J))
#else
            QCFRBS(i,j) = SUBGRID * QCF(I,J,K) / ((1.-RHCRITX)          &
     &               * QSL_CTT(I,J)**CTT_WEIGHT*QSL(I,J)**T_WEIGHT)
#endif
          ELSE
! No ice cloud fraction method defined
          END IF ! ice_fraction_method
!
! ----------------------------------------------------------------------
! 6b Calculate frozen cloud fraction from frozen cloud water content.
! ----------------------------------------------------------------------
          IF (QCFRBS(i,j)  <=  0.) THEN
            CFF(I,j,K) = 0.0
          ELSEIF (0.  <   QCFRBS(i,j)  .AND. (6. * QCFRBS(i,j))  <=  1.) THEN
            CFF(I,j,K) = 0.5 * ((6. * QCFRBS(i,j))**(2./3.))
          ELSEIF (1.  <   (6.*QCFRBS(i,j)) .AND. QCFRBS(i,j)  <   1.) THEN
            PHIQCF = ACOS(ROOTWO * 0.75 * (1. - QCFRBS(i,j)))
            COSQCF = COS((PHIQCF + (4. * PI)) / 3.)
            CFF(I,j,K) = 1. - (4. * COSQCF * COSQCF)
          ELSEIF (QCFRBS(i,j)  >=  1.) THEN
            CFF(I,j,K) = 1.
          END IF
          If (L_eacf) Then  ! Empirically adjusted cloud fraction
            ! Back out QN
            IF (0. <  QCFRBS(i,j)  .AND. (6. * QCFRBS(i,j))  <=  1.) THEN
              QN_IMP=SQRT(2.*CFF(I,j,K))-1.
            ELSEIF (1.  <   (6.*QCFRBS(i,j)) .AND. QCFRBS(i,j) <  1.) THEN
              QN_IMP=1.-SQRT((1.-CFF(I,j,K))*2.)
            ELSE
              QN_IMP = 1.
            ENDIF

            ! Modify QN with EACF relationship
            IF (K >  BL_LEVELS) THEN
              QN_ADJ=(QN_IMP+0.0955)/(1.-0.0955)
            ELSE
              QN_ADJ=(QN_IMP+0.184)/(1.-0.184)
            ENDIF

            ! Recalculate ice cloud fraction with modified QN
            IF (QCFRBS(i,j)  <=  0.) THEN
              CFF(I,j,K) = 0.0
            ELSEIF (QN_ADJ  <=  0.) THEN
              CFF(I,j,K) = 0.5 * (1. + QN_ADJ) * (1. + QN_ADJ)
            ELSEIF (QN_ADJ  <   1.) THEN
              CFF(I,j,K) = 1. - 0.5 * (1.-QN_ADJ) * (1.-QN_ADJ)
            ELSE
              CFF(I,j,K) = 1.
            ENDIF

          End If  ! L_eacf
!

         ELSE ! RHcrit = 1, set cloud fraction to 1 or 0
!
           IF (QCF(I,J,K)  >   0.0) THEN
             CFF(I,J,K) = 1.0
           ELSE
             CFF(I,J,K) = 0.0
           ENDIF
!
         ENDIF
          END DO ! Row_length_do4
        END DO ! Rows_do4

! ----------------------------------------------------------------------
! 6c Calculate combined cloud fraction.
! ----------------------------------------------------------------------
! Rows_do5:
        DO j=1,rows
! Row_length_do5:
         DO I=1,row_length
          IF (CLOUD_FRACTION_METHOD  ==  1) THEN
!           Use minimum overlap condition
            CF(I,j,K) = MIN(CFL(I,j,K)+CFF(I,j,K), 1.0)
          ELSEIF (CLOUD_FRACTION_METHOD  ==  2) THEN
! Calculate possible overlaps between ice and liquid in THIS layer
            OVERLAP_MAX=MIN(CFL(I,J,K),CFF(I,J,K))
            OVERLAP_MIN=MAX(CFL(I,J,K)+CFF(I,J,K)-1.0,0.0)
            OVERLAP_RANDOM=CFL(I,J,K)*CFF(I,J,K)
! Now use the specified degree of overlap to calculate the total
! cloud fraction (= cfice + cfliq - overlap)
            TEMP0=OVERLAP_RANDOM
            TEMP1=0.5*(OVERLAP_MAX-OVERLAP_MIN)
            TEMP2=0.5*(OVERLAP_MAX+OVERLAP_MIN)-OVERLAP_RANDOM
            CF(I,J,K)=CFL(I,J,K)+CFF(I,J,K)                             &
     &              -(TEMP0+TEMP1*OVERLAP_ICE_LIQUID                    &
     &              +TEMP2*OVERLAP_ICE_LIQUID*OVERLAP_ICE_LIQUID)
! Check that the overlap wasnt negative
            CF(I,J,K)=MIN(CF(I,J,K),CFL(I,J,K)+CFF(I,J,K))
          ELSE
! No total cloud fraction method defined
          END IF ! cloud_fraction_method
!
          END DO ! Row_length_do5
        END DO ! Rows_do5
!
      END DO ! Levels_do
!
 9999 CONTINUE ! Error exit
      RETURN
      END SUBROUTINE LS_CLD
#endif
