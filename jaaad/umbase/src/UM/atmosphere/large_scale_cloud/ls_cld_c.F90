#if defined(A09_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale Cloud Scheme.
! Subroutine Interface:
! ======================================================================
!
!+ Large-scale Cloud Scheme Compression routine (Cloud points only).
! Subroutine Interface:
      SUBROUTINE LS_CLD_C(                                              &
     & P_F,RHCRIT,QSL_F,QN_F,Q_F,T_F                                    &
     &,QCL_F,CF_F,GRID_QC_F,BS_F                                        &
     &,INDEX,POINTS,row_length,rows,rhc_row_length,rhc_rows             &
     &,BL_LEVELS,K, L_eacf, l_mixing_ratio)
!
      IMPLICIT NONE
!
! Purpose: Calculates liquid cloud water amounts and cloud amounts,
!          temperature and specific humidity from cloud-conserved and
!          other model variables. This is done for one model level.
!
! Current Owner of Code: A. C. Bushell
!
! History:
! Version   Date     Comment
!  4.4    14-11-96   Original Code (A. C. Bushell from Wilson/Ballard)
!  6.2    07-11-05   Include the EACF parametrization for cloud
!                    fraction.                            D. Wilson
!  6.4    08-08-06   Include mixing ratio formulation. D. Wilson
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP No.29
!
!  Global Variables:----------------------------------------------------
#include "c_r_cp.h"
#include "c_epslon.h"
#include "c_lheat.h"
!
!  Subroutine Arguments:------------------------------------------------
      INTEGER                                                           &
                        !, INTENT(IN)
     &  row_length,rows                                                 &
     &, rhc_row_length,rhc_rows                                         &
!       No. of gridpoints being processed.
     &, BL_LEVELS                                                       &
!       No. of boundary layer levels
     &, K                                                               &
!       Level no.
     &,POINTS                                                           &
!       No. of gridpoints with non-zero cloud
     &,INDEX(row_length*rows,2)
!       INDEX for  points with non-zero cloud from lowest model level.
!
      REAL                                                              &
                        !, INTENT(IN)
     & RHCRIT(rhc_row_length,rhc_rows)                                  &
!       Critical relative humidity.  See the paragraph incorporating
!       eqs P292.11 to P292.14.
     &,P_F(row_length,rows)                                             &
!       pressure (Pa).
     &,QSL_F(row_length,rows)                                           &
!       saturated humidity at temperature TL, and pressure P_F
     &,QN_F(row_length,rows)
!       Normalised super/subsaturation ( = QC/BS).
!
      LOGICAL                                                           &
                        !, INTENT(IN)
     & L_eacf                                                           &
                        !  Use empirically adjusted cloud fraction
     &,L_mixing_ratio   !  Use mixing ratio formulation
!
      REAL                                                              &
                        !, INTENT(INOUT)
     & Q_F(row_length,rows)                                             &
!       On input : Vapour + liquid water content (QW) (kg per kg air).
!       On output: Specific humidity at processed levels
!                   (kg water per kg air).
     &,T_F(row_length,rows)
!       On input : Liquid water temperature (TL) (K).
!       On output: Temperature at processed levels (K).
!
      REAL                                                              &
                        !, INTENT(OUT)
     & QCL_F(row_length,rows)                                           &
!       Cloud liquid water content at processed levels (kg per kg air).
     &,CF_F(row_length,rows)                                            &
!       Liquid cloud fraction at processed levels.
     &,GRID_QC_F(row_length,rows)                                       &
!       Super/subsaturation on processed levels. Input initially RMDI.
     &,BS_F(row_length,rows)
!       Value of bs at processed levels. Input initialized to RMDI.
!
!  Local parameters and other physical constants------------------------
      REAL ALPHL,LCRCP                  ! Derived parameters.
      PARAMETER (                                                       &
     & ALPHL=EPSILON*LC/R                                               &
                                        ! For liquid AlphaL calculation.
     &,LCRCP=LC/CP                                                      &
                                        ! Lat ht of condensation/Cp.
     &)
      REAL WTN                          ! Weighting for ALPHAL iteration
      INTEGER                                                           &
     & ITS                              ! Total number of iterations
      PARAMETER (ITS=5,WTN=0.75)
!
!  Local scalars--------------------------------------------------------
!
!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
      REAL                                                              &
     & AL                                                               &
                         ! LOCAL AL (see equation P292.6).
     &,ALPHAL                                                           &
                         ! LOCAL ALPHAL (see equation P292.5).
     &,QN_ADJ                                                           &
     &,RHCRITX          ! scalar copy of RHCRIT(I,J)
      INTEGER                                                           &
     & MULTRHC          ! Zero if (rhc_row_length*rhc_rows) le 1, else 1
!
!  (b) Others.
      INTEGER   I,II,IJ,N   ! Loop counters:I,II-horizontal field index.
!                                       : N - iteration counter.
!
!  Local dynamic arrays-------------------------------------------------
!    7 blocks of real workspace are required.
      REAL                                                              &
     & P(POINTS)                                                        &
!       Pressure  (Pa).
     &,QS(POINTS)                                                       &
!       Saturated spec humidity for temp T.
     &,QCN(POINTS)                                                      &
!       Cloud water normalised with BS.
     &,T(POINTS)                                                        &
!       temperature.
     &,Q(POINTS)                                                        &
!       specific humidity.
     &,BS(POINTS)                                                       &
!       Sigmas*sqrt(6): sigmas the parametric standard deviation of
!       local cloud water content fluctuations.
     &,ALPHAL_NM1(POINTS)
!       ALPHAL at previous iteration.
!
!  External subroutine calls: ------------------------------------------
      EXTERNAL QSAT_WAT
!
!- End of Header
!
! ==Main Block==--------------------------------------------------------
! Operate on INDEXed points with non-zero cloud fraction.
! ----------------------------------------------------------------------
      IF ( (rhc_row_length * rhc_rows)  >   1) THEN
        MULTRHC = 1
      ELSE
        MULTRHC = 0
      END IF
!
!        RHCRITX = RHCRIT(1,1)
! Points_do1:
!CDIR NODEP
      DO I=1, POINTS
        II = INDEX(I,1)
        IJ = INDEX(I,2)
        IF ( MULTRHC== 1) THEN
          RHCRITX = RHCRIT(II,IJ)
        ELSE
          RHCRITX = RHCRIT(1,1)
        ENDIF
        P(I)  = P_F(II,IJ)
        QCN(I)= QN_F(II,IJ)
! ----------------------------------------------------------------------
! 1. Calculate ALPHAL (eq P292.5) and AL (P292.6).
!    CAUTION: T_F acts as TL (input value) until update in final section
!    CAUTION: Q_F acts as QW (input value) until update in final section
! ----------------------------------------------------------------------
!
       ALPHAL = ALPHL * QSL_F(II,IJ) / (T_F(II,IJ) * T_F(II,IJ)) !P292.5
        AL = 1.0 / (1.0 + (LCRCP * ALPHAL))                    ! P292.6
        ALPHAL_NM1(I) = ALPHAL
!
! Rhcrit_if1:
        IF (RHCRITX  <   1.) THEN
! ----------------------------------------------------------------------
! 2. Calculate BS (ie. sigma*sqrt(6), where sigma is
!    as in P292.14) and normalised cloud water QCN=qc/BS, using eqs
!    P292.15 & 16 if RHcrit < 1.
! N.B. QN (input) is initially in QCN
! N.B. QN does not depend on AL and so CF and QCN can be calculated
!      outside the iteration (which is performed in LS_CLD_C).
!      QN is > -1 for all points processed so CF > 0.
! ----------------------------------------------------------------------
!
          BS(I) = (1.0-RHCRITX) * AL * QSL_F(II,IJ)  ! P292.14
          IF (QCN(I)  <=  0.) THEN
            CF_F(II,IJ) = 0.5 * (1. + QCN(I)) * (1. + QCN(I))
            QCN(I)= (1. + QCN(I)) * (1. + QCN(I)) * (1. + QCN(I)) / 6.
          ELSEIF (QCN(I)  <   1.) THEN
            CF_F(II,IJ) = 1. - 0.5 * (1. - QCN(I)) * (1. - QCN(I))
            QCN(I)=QCN(I) + (1.-QCN(I)) * (1.-QCN(I)) * (1.-QCN(I))/6.
          ELSE ! QN  >=  1
            CF_F(II,IJ) = 1.
          END IF ! QCN_if
!
! ----------------------------------------------------------------------
! 3.b If necessary, modify cloud fraction using empirically adjusted
!     cloud fraction parametrization, but keep liquid content the same.
! ----------------------------------------------------------------------
          If (L_eacf) Then
            ! Adjust QN according to EACF parametrization

            If (K <= BL_LEVELS) Then
              QN_ADJ=(QN_F(II,IJ)+0.184)/(1.-0.184)
            Else
              QN_ADJ=(QN_F(II,IJ)+0.0955)/(1.-0.0955)
            EndIf

!         Calculate cloud fraction using adjusted QN
            If (QN_ADJ  <=  0.) Then
              CF_F(II,IJ) = 0.5 * (1.0 + QN_ADJ) * (1.0 + QN_ADJ)
            ElseIf (QN_ADJ  <   1.) Then
              CF_F(II,IJ) = 1. - 0.5 * (1.0 - QN_ADJ) * (1.0 - QN_ADJ)
            Else ! QN_ADJ  >=  1
              CF_F(II,IJ) = 1.
            EndIf ! QN_ADJ_if

         EndIf  ! l_eacf
!
        ELSE ! i.e. if RHcrit = 1
! ----------------------------------------------------------------------
! 3.a If RHcrit = 1., all points processed have QN > 0 and CF = 1.
! ----------------------------------------------------------------------
          BS(I) = AL
          CF_F(II,IJ) = 1.
        END IF ! Rhcrit_if1
!
! ----------------------------------------------------------------------
! 3.1 Calculate 1st approx. to qc (store in QCL)
! ----------------------------------------------------------------------
!
        QCL_F(II,IJ) = QCN(I) * BS(I)
!
! ----------------------------------------------------------------------
! 3.2 Calculate 1st approx. specific humidity (total minus cloud water)
! ----------------------------------------------------------------------
!
        Q(I) = Q_F(II,IJ) - QCL_F(II,IJ)
!
! ----------------------------------------------------------------------
! 3.3 Calculate 1st approx. to temperature, adjusting for latent heating
! ----------------------------------------------------------------------
!
        T(I) = T_F(II,IJ) + LCRCP*QCL_F(II,IJ)
      END DO ! Points_do1
!
! ----------------------------------------------------------------------
! 4. Iteration to find better cloud water values.
! ----------------------------------------------------------------------
! Its_if:
      IF (ITS  >=  2) THEN
! Its_do:
        DO N=2, ITS
!
! DEPENDS ON: qsat_wat_mix
          CALL qsat_wat_mix(QS,T,P,POINTS,l_mixing_ratio)
! Points_do2:
        RHCRITX = RHCRIT(1,1)
          DO I=1, POINTS
            II = INDEX(I,1)
            IJ = INDEX(I,2)
            IF ( MULTRHC== 1) THEN
              RHCRITX = RHCRIT(II,IJ)
            ELSE
              RHCRITX = RHCRIT(1,1)
            ENDIF
! T_if:
            IF (T(I)  >   T_F(II,IJ)) THEN
!           NB. T > TL implies cloud fraction > 0.
              ALPHAL = (QS(I) - QSL_F(II,IJ)) / (T(I) - T_F(II,IJ))
              ALPHAL = WTN * ALPHAL + (1.0 - WTN) * ALPHAL_NM1(I)
              ALPHAL_NM1(I) = ALPHAL
              AL = 1.0 / (1.0 + (LCRCP * ALPHAL))
! Rhcrit_if2:
              IF (RHCRITX  <   1.) THEN
                BS(I) = (1.0-RHCRITX) * AL * QSL_F(II,IJ)
!                                                             P292.14
              ELSE
                BS(I) = AL
              END IF  ! Rhcrit_if2
!
! ----------------------------------------------------------------------
! 4.1 Calculate Nth approx. to qc (store in QCL).
! ----------------------------------------------------------------------
!
              QCL_F(II,IJ) = QCN(I) * BS(I)
!
! ----------------------------------------------------------------------
! 4.2 Calculate Nth approx. spec. humidity (total minus cloud water).
! ----------------------------------------------------------------------
!
              Q(I) = Q_F(II,IJ) - QCL_F(II,IJ)
!
! ----------------------------------------------------------------------
! 4.3 Calculate Nth approx. to temperature, adjusting for latent heating
! ----------------------------------------------------------------------
!
              T(I) = T_F(II,IJ) + LCRCP * QCL_F(II,IJ)
!
            END IF ! T_if
          END DO ! Points_do2
        END DO ! Its_do
      END IF ! Its_if
!
! ----------------------------------------------------------------------
! 5. Finally scatter back cloud point results to full field arrays.
!    CAUTION: T_F updated from TL (input) to T (output)
!    CAUTION: Q_F updated from QW (input) to Q (output)
! ----------------------------------------------------------------------
!
!DIR$ IVDEP
! Points_do3:
      DO I=1,POINTS
       II = INDEX(I,1)
       IJ = INDEX(I,2)
        Q_F(II,IJ) = Q(I)
        T_F(II,IJ) = T(I)
        GRID_QC_F(II,IJ) = BS(I) * QN_F(II,IJ)
        BS_F(II,IJ) = BS(I)
      END DO ! Points_do3
!
      RETURN
      END SUBROUTINE LS_CLD_C
! ======================================================================
#endif
