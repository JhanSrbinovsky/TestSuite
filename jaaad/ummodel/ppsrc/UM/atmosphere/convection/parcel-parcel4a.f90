
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE PARCEL-------------------------------------------------
!LL
!LL  PURPOSE : COMPLETES LIFTING OF THE PARCEL FROM LAYER K TO K+1
!LL
!LL            CALL SUBROUTINE DETRAIN, TERM_CON, CLOUD_W
!LL
!LL            AN INITIAL MASS FLUX IS CALCULATED
!LL
!LL            SUBROUTINE DETRAIN CARRIES OUT THE FORCED DETRAINMENT
!LL            CALCULATION
!LL
!LL            SUBROUTINE TERM_CON TESTS FOR ANY CONVECTION WHICH IS
!LL            TERMINATING IN LAYER K+1
!LL
!LL            SUBROUTINE CLOUD_W CARRIES OUT THE CLOUD MICROPHYSICS
!LL            CALCULATION
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
!LL                                   THE RELEVANT POINTS IN CONVECT
!LL
!LL  PROGRAMMING STANDARDS :
!LL
!LL  LOGICAL COMPONENTS COVERED: P27
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
      SUBROUTINE PARCEL (K,NPNTS,NLEV,PSTAR,THEKP1,THEK,QEKP1,QEK,      &
     &                   QSEK,QSEKP1,DQSK,DQSKP1,BLAND,BWKP1,           &
     &                   DELTAK,FLXK,THPK,QPK,THRK,QRK,                 &
     &                   BTERM,THPKP1,QPKP1,PREKP1,XPK,XPKP1,FLXKP1,    &
     &                   QCLPK,QCFPK,QCLPKP1,QCFPKP1,                   &
     &                   XSQKP1,THPI,QPI,EXPI,BGMK,BGMKP1,BLOWST,RBUOY, &
     &                   XSBMIN,CCA,ICCB,ICCT,TCW,DEPTH,                &
     &                   EKP14,EKP34,AMDETK,DELPKP1,PK,PKP1,            &
     &                   EXK,EXKP1,DELEXKP1,CCLWP,CCW,                  &
     &                   LCCA,LCBASE,LCTOP,LCCLWP,L_SHALLOW,            &
     &                   RBUOY_P_HERE,THE_HERE, THP_HERE,QE_HERE,       &
     &                   QP_HERE,RBUOY_P_OLD,AD_ON, NEW_TERMC,          &
     &                   L_Q_INTERACT,START_LEV,FLXKP12)

      Use cv_run_mod, Only:                                             &
          r_det

      IMPLICIT NONE
!
!----------------------------------------------------------------------
! MODEL CONSTANTS
!----------------------------------------------------------------------
!
      REAL C_DEEP,                                                      &
                      ! CONSTANTS USED TO DETERMINE INITIAL CONVECTIVE
     &     D_DEEP     ! MASS FLUX FROM PARCEL BUOYANCY FOR DEEP
                      ! CONVECTION
!
      REAL C_SHALLOW,                                                   &
                      ! CONSTANTS USED TO DETERMINE INITIAL CONVECTIVE
     &     D_SHALLOW  ! MASS FLUX FROM PARCEL BUOYANCY FOR SHALLOW
                      ! CONVECTION
!
      REAL C_MID,                                                       &
                      ! CONSTANTS USED TO DETERMINE INITIAL CONVECTIVE
     &     D_MID      ! MASS FLUX FROM PARCEL BUOYANCY FOR MID-LEVEL
                      ! CONVECTION
!
      PARAMETER (C_DEEP = 5.17E-4, D_DEEP = 0.0)
      PARAMETER (C_SHALLOW = 5.17E-4, D_SHALLOW = 0.0)
      PARAMETER (C_MID = 5.17E-4, D_MID = 0.0)
!
!----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!----------------------------------------------------------------------
!
      INTEGER NPNTS          ! IN VECTOR LENGTH
!
      INTEGER NLEV           ! IN NUMBER OF MODEL LEVELS
!
      INTEGER NDET           ! COMPRESSED VECTOR LENGTH FOR
                             ! FORCED DETRAINMENT CALCULATION
!
      INTEGER K              ! IN PRESENT MODEL LAYER
!
      INTEGER I              ! LOOP COUNTER
!
!
!----------------------------------------------------------------------
! VARAIBLES WHICH ARE INPUT
!----------------------------------------------------------------------
!
      REAL THEK(NPNTS)       ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K (K)
!
      REAL THEKP1(NPNTS)     ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (K)
!
      REAL QEK(NPNTS)        ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
!
      REAL QEKP1(NPNTS)      ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL QSEKP1(NPNTS)     ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL DQSKP1(NPNTS)     ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT IN LAYER K+1
                             !    (KG/KG/K)
!
      REAL PSTAR(NPNTS)      ! IN SURFACE PRESSURE (PA)
!
      REAL THPK(NPNTS)       ! IN PARCEL POTENTIAL TEMPERATURE
                             !    IN LAYER K (KG/KG)
!
      REAL QPK(NPNTS)        ! IN PARCEL MIXING RATIO IN LAYER K (KG/KG)
!
      REAL QCLPK(NPNTS)      ! IN PARCEL LIQUID CONDENSATE MIXING RATIO
!                                 IN LAYER K (KG/KG)
!
      REAL QCFPK(NPNTS)      ! IN PARCEL FROZEN CONDENSATE MIXING RATIO
!                                 IN LAYER K (KG/KG)
!
      REAL XSQKP1(NPNTS)     ! IN EXCESS PARCEL WATER AFER LIFTING FROM
                             !    LAYER K TO K+1 (KG/KG)
!
      REAL RBUOY(NPNTS)      ! IN PARCEL BUOYANCY IN LAYER K+1 (K)
!
      REAL QSEK(NPNTS)       ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
!
      REAL DQSK(NPNTS)       ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT OF LAYER K
                             !    (KG/KG/K)
!
      REAL THPI(NPNTS)       ! IN INITIAL PARCEL POTENTIAL TEMPERATURE
                             !    (K)
!
      REAL QPI(NPNTS)        ! IN INITIAL PARCEL MIXING RATIO (KG/KG)
!
      REAL EXPI(NPNTS)       ! IN INITIAL PARCEL EXNER PRESSURE
!
      REAL XPK(NPNTS)        ! IN PARCEL CLOUD WATER IN LAYER K (KG/KG)
!
      LOGICAL BWKP1(NPNTS)   ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
!
      LOGICAL BGMK(NPNTS)    ! IN MASK FOR PARCELS WHICH ARE
                             !    SATURATED IN LAYER K
!
      LOGICAL BLAND(NPNTS)   ! IN LAND/SEA MASK
!
      LOGICAL BLOWST(NPNTS)  ! IN MASK FOR THOSE POINTS AT WHICH
                             !    STABILITY IS LOW ENOUGH FOR
                             !    CONVECTION TO OCCUR
!
      LOGICAL L_SHALLOW(NPNTS) ! IN MASK FOR POINTS WHERE CONVECTION
                               !    IS EXPECTED TO BE SHALLOW
!
      LOGICAL L_Q_INTERACT   ! IN Switch allows overwriting of parcel
!                                 variables when calculating condensate
!                                 increments (will alter results).
!
      REAL EKP14(NPNTS)      ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+1/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
!
      REAL EKP34(NPNTS)      ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !   K+3/4 MULTIPLIED BY APPROPRIATE
                             !   LAYER THICKNESS
!
      REAL AMDETK(NPNTS)     ! IN MIXING DETRAINMENT COEFFICIENT
                             !    AT LEVEL K MULTIPLIED BY
                             !    APPROPORIATE LAYER THICKNESS
!
      REAL DELPKP1(NPNTS)    ! IN PRESSURE DIFFERENCE ACROSS
                             !    LAYER K+1 (PA)
!
      REAL PK(NPNTS)         ! IN PRESSURE AT LEVEL K (PA)
!
      REAL PKP1(NPNTS)       ! IN PRESSURE AT LEVEL K+1 (PA)
!
      REAL EXK(NPNTS)        ! IN EXNER FUNCTION AT LEVEL K
!
      REAL EXKP1(NPNTS)      ! IN EXNER FUNCTION AT LEVEL K+1
!
      REAL DELEXKP1(NPNTS)   ! IN DIFFERENCE IN EXNER FUNCTION ACROSS
                             !    LAYER K+1
!
      REAL XSBMIN(NPNTS)     ! IN THRESHOLD FOR FORCED DETRAINMENT
!                            !    Function of delta P
!
      INTEGER START_LEV(NPNTS)  ! IN LEVEL AT WHICH CONVECTION INITIATED
!
      REAL, INTENT(IN):: RBUOY_P_OLD(NPNTS)  ! RBUOY ON PREVIOUS LEVEL
                                             ! FOR ADAPTIVE CALC

      INTEGER, INTENT(IN):: AD_ON            ! FLAG FOR ADAPTIVE SWITCH
                                             !(1 FOR ON)
!
      INTEGER, INTENT(IN):: NEW_TERMC        ! FLAG FOR SIMPLIFIED
                                             !TERMINATION
                                             !OF CONVECTION (1 FOR ON)
!
!
!---------------------------------------------------------------------
! VARAIBLES WHICH ARE BOTH INPUT AND OUTPUT
!---------------------------------------------------------------------
!
      REAL THPKP1(NPNTS)     ! INOUT
                             ! IN  ESTIMATE OF PARCEL POTENTIAL
                             !     TEMPERATURE IN LAYER K+1 AFTER
                             !     ENTRAINMENT AND LATENT HEATING (K)
                             ! OUT FINAL PARCEL POTENTIAL TEMPERATURE
                             !     IN LAYER K+1 (AFTER FORCED
                             !     DETRAINEMNT) (K)
!
      REAL QPKP1(NPNTS)      ! INOUT
                             ! IN  ESTIMATE OF PARCEL MIXING RATIO
                             !     IN LAYER K+1 AFTER ENTRAINMENT AND
                             !     LATENT HEATING (KG/KG)
                             ! OUT FINAL PARCEL MIXING RATIO
                             !     IN LAYER K+1 (AFTER FORCED
                             !     DETRAINEMNT) (KG/KG)
!
      REAL QCLPKP1(NPNTS)    ! INOUT
!                              IN  ESTIMATED PARCEL LIQUID CONDENSATE
!                                  MIXING RATIO IN LAYER K+1 AFTER
!                                  ENTRAINMENT ONLY (KG/KG)
!                              OUT FINAL PARCEL LIQUID CONDENSATE MIXING
!                                  RATIO IN LAYER K+1 (KG/KG)
!
      REAL QCFPKP1(NPNTS)    ! INOUT
!                              IN  ESTIMATED PARCEL FROZEN CONDENSATE
!                                  MIXING RATIO IN LAYER K+1 AFTER
!                                  ENTRAINMENT ONLY (KG/KG)
!                              OUT FINAL PARCEL FROZEN CONDENSATE MIXING
!                                  RATIO IN LAYER K+1 (KG/KG)
!
      REAL FLXK(NPNTS)       ! INOUT
                             ! IN  PARCEL MASSFLUX IN LAYER K
                             !     (NON-ZERO IF CONVECTION IS NOT
                             !     INITIATED FROM LAYER K) (PA/S)
                             ! OUT PARCEL MASSFLUX IN LAYER K
                             !     (SET IF CONVECTION IS INITIATED
                             !     IN LAYER K) (PA/S)
!
      LOGICAL BGMKP1(NPNTS)  ! INOUT
                             ! IN  MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K+1
                             !     CALCULATED ON THE BASIS OF
                             !     INPUT PARCEL POTENTIAL TEMPERATURE
                             !     AND MIXING RATIO
                             ! OUT MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K+1 CALCULATED
                             !     FORM PARCEL TEMPERATURE AND
                             !     MIXING RATIO AFTER FORCED
                             !     DETARINMENT CALCULATION
!
      REAL TCW(NPNTS)        ! INOUT
                             ! IN  TOTAL CONDENSED WATER CONTENT
                             !     SUMMED UPTO LAYER K (KG/M**2/S)
                             ! OUT UPDATED TOTAL CONDENSED WATER
                             !     CONTENT SUMMED UPTO LAYER K+1
                             !     (KG/M**2/S)
!
      REAL DEPTH(NPNTS)      ! INOUT
                             ! IN  DEPTH OF CONVECTIVE CLOUD TO
                             !     LAYER K (M)
                             ! OUT UPDATED DEPTH OF CONVECTIVE
                             !     CLOUD TO LAYER K+1 (M)
!
      REAL CCLWP(NPNTS)      ! INOUT
                             ! IN  CONDENSED WATER PATH
                             !     SUMMED UPTO LAYER K (KG/M**2)
                             ! OUT UPDATED CONDENSED WATER PATH
                             !     SUMMED UPTO LAYER K+1 (KG/M**2)
!
!
!---------------------------------------------------------------------
! VARIABLES WHICH ARE OUTPUT
!---------------------------------------------------------------------
!
      LOGICAL BTERM(NPNTS)   ! OUT MASK FOR PARCELS WHICH TERMINATE IN
                             !     LAYER K+1
!
      REAL PREKP1(NPNTS)     ! OUT PRECIPITATION FROM PARCEL AS IT
                             !     RISES FROM LAYER K TO K+1 (KG/M**2/S)
!
      REAL THRK(NPNTS)       ! OUT PARCEL DETRAINMENT POTENTIAL
                             !     TEMPERATURE IN LAYER K (K)
!
      REAL QRK(NPNTS)        ! OUT PARCEL DETRAINMENT MIXING RATIO
                             !     IN LAYER K (KG/KG)
!
      REAL XPKP1(NPNTS)      ! OUT PARCEL CLOUD WATER IN LAYER K+1
                             !     (KG/KG)
!
      REAL FLXKP1(NPNTS)     ! OUT PARCEL MASSFLUX IN LAYER K+1 (PA/S)
!
      REAL DELTAK(NPNTS)     ! OUT PARCEL FORCED DETRAINMENT
                             !     COEFFICIENT IN LAYER K
                             !     MULTIPLIED BY APPROPRIATE
                             !     LAYER THICKNESS
!
      REAL CCA(NPNTS)        ! OUT CONVECTIVE CLOUD AMOUNT (%)
!
      INTEGER ICCB(NPNTS)    ! OUT CONVECTIVE CLOUD BASE LEVEL
!
      INTEGER ICCT(NPNTS)    ! OUT CONVECTIVE CLOUD TOP LEVEL
!
      REAL CCW(NPNTS)        ! OUT CONVECTIVE CLOUD LIQUID WATER
                             ! (G/KG) ON MODEL LEVELS
!
      REAL LCCA(NPNTS)       ! OUT LOWEST CONV.CLOUD AMOUNT (%)
!
      INTEGER LCBASE(NPNTS)  ! OUT LOWEST CONV.CLOUD BASE LEVEL
!
      INTEGER LCTOP(NPNTS)   ! OUT LOWEST CONV.CLOUD TOP LEVEL
!
      REAL LCCLWP(NPNTS)     ! OUT LOWEST CONV.CLOUD LIQ.WATER PATH
!
      REAL FLXKP12(NPNTS)    ! OUT HALF LEVEL MASS FLUX
!
      REAL, INTENT(OUT):: RBUOY_P_HERE(NPNTS) ! BUOYANCY FOR SCM DIAGS

      REAL, INTENT(OUT):: THE_HERE(NPNTS)   ! TH ENVIRON FOR SCM DIAGS

      REAL, INTENT(OUT):: THP_HERE(NPNTS)   ! TH PARCEL FOR SCM DIAGS

      REAL, INTENT(OUT):: QE_HERE(NPNTS)    ! Q ENVIRON FOR SCM DIAGS

      REAL, INTENT(OUT):: QP_HERE(NPNTS)    ! Q PARCEL FOR SCM DIAGS
!
!---------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!---------------------------------------------------------------------
!
      REAL XSBMIN_AD(NPNTS)  ! XSBMIN ADAPTIVE (NOTE, WILL BE DIFFERENT
                             ! AT DIFFERENT POINTS
!
      REAL THEK_C(NPNTS)     ! COMPRESSED POTENTIAL TEMPERATURE OF
                             ! CLOUD ENVIRONMENT IN LAYER K (K)
!
      REAL THEKP1_C(NPNTS)   ! COMPRESSED POTENTIAL TEMPERATURE OF
                             ! CLOUD ENVIRONMENT IN LAYER K+1 (K)
!
      REAL QEK_C(NPNTS)      ! COMPRESSED MIXING RATIO OF CLOUD
                             ! ENVIRONMENT IN LAYER K (KG/KG)
!
      REAL QEKP1_C(NPNTS)    ! COMPRESSED MIXING RATIO OF CLOUD
                             ! ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL QSEK_C(NPNTS)     ! COMPRESSED SATURATION MIXING RATIO OF
                             ! CLOUD ENVIRONMENT IN LAYER K (KG/KG)
!
      REAL DQSK_C(NPNTS)     ! COMPRESSED GRADIENT OF SATURATION MIXING
                             ! RATIO WITH POTENTIAL TEMPERATURE FOR THE
                             ! CLOUD ENVIRONMENT OF LAYER K (KG/KG/K)
!
      REAL QSEKP1_C(NPNTS)   ! COMPRESSED SATURATION MIXING RATIO OF
                             ! CLOUD ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL DQSKP1_C(NPNTS)   ! COMPRESSED GRADIENT OF SATURATION MIXING
                             ! RATIO WITH POTENTIAL TEMPERATURE FOR
                             ! THE CLOUD ENVIRONMENT IN LAYER K+1
                             ! (KG/KG/K)
!
      REAL THPK_C(NPNTS)     ! COMPRESSED PARCEL POTENTIAL
                             ! TEMPERATURE IN LAYER K (K)
!
      REAL QPK_C(NPNTS)      ! COMPRESSED PARCEL MIXING RATIO IN
                             ! LAYER K (KG/KG)
!
      REAL THPKP1_C(NPNTS)   ! COMPRESSED PARCEL POTENTIAL
                             ! TEMPERATURE IN LAYER K+1 (K)
!
      REAL QPKP1_C(NPNTS)    ! COMPRESSED PARCEL MIXING RATIO
                             ! IN LAYER K+1 (KG/KG)
!
      REAL XSQKP1_C(NPNTS)   ! EXCESS PARCEL WATER AFER LIFTING
                             ! FROM LAYER K TO K+1 (KG/KG)
!
      REAL THRK_C(NPNTS)     ! COMPRESSED PARCEL DETRAINMENT
                             ! POTENTIAL TEMPERATURE IN LAYER K (K)
!
      REAL QRK_C(NPNTS)      ! COMPRESSED PARCEL DETRAINMENT MIXING
                             ! RATIO IN LAYER K (KG/KG)
!
      REAL DELTAK_C(NPNTS)   ! COMPRESSED PARCEL FORCED DETRAINMENT
                             ! COEFFICIENT IN LAYER K
                             ! MULTIPLIED BY APPROPRIATE
                             ! LAYER THICKNESS
!
      REAL EKP14_C(NPNTS)    ! COMPRESSED IN ENTRAINMENT COEFFICIENT AT
                             ! LEVEL K+1/4 MULTIPLIED BY APPROPRIATE
                             ! LAYER THICKNESS
!
      REAL EKP34_C(NPNTS)    ! COMPRESSED ENTRAINMENT COEFFICIENT AT
                             ! LEVEL K+3/4 MULTIPLIED BY APPROPRIATE
                             ! LAYER THICKNESS
!
      REAL PK_C(NPNTS)       ! COMPRESSED PRESSURE AT LEVEL K (PA)
!
      REAL PKP1_C(NPNTS)     ! COMPRESSED PRESSURE AT LEVEL K+1 (PA)
!
      REAL XSBMIN_C(NPNTS)   ! COMPRESSED THRESHOLD FOR FORCED
!                            ! DETRAINMENT Function of delta P
!
      REAL EXK_C(NPNTS)      ! COMPRESSED EXNER FUNCTION AT LEVEL K
!
      REAL EXKP1_C(NPNTS)    ! COMPRESSED EXNER FUNCTION AT LEVEL K+1
!
      LOGICAL BWKP1_C(NPNTS) ! COMPRESSED MASK FOR WHETHER CONDENSATE
                             ! IS LIQUID IN LAYER K+1
!
      LOGICAL BGMK_C(NPNTS)  ! COMPRESSED MASK FOR PARCELS WHICH ARE
                             ! SATURATED IN LAYER K
!
      LOGICAL BGMKP1_C(NPNTS) ! COMPRESSED MASK FOR PARCELS
                              ! WHICH ARESATURATED IN LAYER K+1
!
      INTEGER INDEX1(NPNTS)  ! INDEX FOR COMPRESS AND EXPAND
!
      LOGICAL BDETK(NPNTS)   ! MASK FOR POINTS UNDERGOING
                             ! FORCED DETRAINMENT
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! EXTERNAL ROUTINES CALLED
!----------------------------------------------------------------------
!
      EXTERNAL DETRAIN,TERM_CON,CLOUD_W
!
!*--------------------------------------------------------------------
!
!

      DO 5 I=1,NPNTS
!L
!L---------------------------------------------------------------------
!L CALCULATE MASK FOR THOSE POINTS UNDERGOING FORCED DETRAINMENT
!L
!L UM DOCUMENTATION PAPER 27
!L SECTION (6), EQUATION (23)
!L---------------------------------------------------------------------
!L
       BDETK(I) = RBUOY(I)  <   XSBMIN(I)
!
!  Calculate XSBMIN_AD
!  Use adaptive detrainment if ADAPTIVE ON flag is set to 1
!  making sure that adaptive min. buoyancy doesn't fall below 0.0
       IF (AD_ON  /=  1) THEN
         XSBMIN_AD(I) = XSBMIN(I)
       ELSE
         XSBMIN_AD(I) = R_DET * RBUOY_P_OLD(I) + (1-R_DET)*RBUOY(I)
         XSBMIN_AD(I) = MAX(XSBMIN_AD(I),0.0)
       END IF
!
!  Set mask of points for forced detrainment
       BDETK(I) = RBUOY(I)  <   XSBMIN_AD(I)

      RBUOY_P_HERE(I)=RBUOY(I)
!
      THE_HERE(I)=THEK(I)
      QE_HERE(I)=QEK(I)
      THP_HERE(I)=THPK(I)
      QP_HERE(I)=QPK(I)

!
   5  CONTINUE
!L
!L
!L----------------------------------------------------------------------
!L  COMPRESS ALL INPUT ARRAYS FOR THE FORCED DETRAINMENT CALCULATIONS
!L----------------------------------------------------------------------
!L
      NDET = 0
      DO 10 I=1,NPNTS
       IF (BDETK(I))THEN
         NDET = NDET + 1
         INDEX1(NDET) = I
       END IF
  10  CONTINUE
!
      IF (NDET  /=  0) THEN
        DO I=1,NDET
          THEK_C(I)  = THEK(INDEX1(I))
          QEK_C(I)   = QEK(INDEX1(I))
        ENDDO
        DO I=1,NDET
          THPK_C(I)  = THPK(INDEX1(I))
          QPK_C(I)   = QPK(INDEX1(I))
        ENDDO
        DO I=1,NDET
          QSEK_C(I)  = QSEK(INDEX1(I))
          DQSK_C(I)  = DQSK(INDEX1(I))
        ENDDO
        DO I=1,NDET
          THEKP1_C(I)= THEKP1(INDEX1(I))
          QEKP1_C(I) = QEKP1(INDEX1(I))
        ENDDO
        DO I=1,NDET
          THPKP1_C(I)= THPKP1(INDEX1(I))
          QPKP1_C(I) = QPKP1(INDEX1(I))
        ENDDO
        DO I=1,NDET
          QSEKP1_C(I)= QSEKP1(INDEX1(I))
          DQSKP1_C(I)= DQSKP1(INDEX1(I))
        ENDDO
        DO I=1,NDET
          XSQKP1_C(I)= XSQKP1(INDEX1(I))
          EKP14_C(I) = EKP14(INDEX1(I))
        ENDDO
        DO I=1,NDET
          EKP34_C(I) = EKP34(INDEX1(I))
          PK_C(I)    = PK(INDEX1(I))
        ENDDO
        DO I=1,NDET
          PKP1_C(I)  = PKP1(INDEX1(I))
          EXK_C(I)   = EXK(INDEX1(I))
!           CHANGE XSBMIN_C TO ADAPTIVE VERSION
          XSBMIN_C(I)= XSBMIN_AD(INDEX1(I))
        ENDDO
        DO I=1,NDET
          EXKP1_C(I) = EXKP1(INDEX1(I))
!
          BGMK_C(I)  = BGMK(INDEX1(I))
        ENDDO
        DO I=1,NDET
          BGMKP1_C(I)= BGMKP1(INDEX1(I))
          BWKP1_C(I) = BWKP1(INDEX1(I))
        ENDDO
!L
!L-------------------------------------------------------------------
!L DETRAINMENT CALCULATION
!L
!L SUBROUTINE DETRAIN
!L
!L UM DOCUMENTATION PAPER 27
!L SECTION (6)
!L-------------------------------------------------------------------
!L
! DEPENDS ON: detrain
         CALL DETRAIN (NDET,THEK_C,QEK_C,THPK_C,QPK_C,                  &
     &                 QSEK_C,DQSK_C,BGMK_C,THEKP1_C,                   &
     &                 QEKP1_C,THPKP1_C,QPKP1_C,QSEKP1_C,               &
     &                 DQSKP1_C,BGMKP1_C,BWKP1_C,                       &
     &                 XSQKP1_C,DELTAK_C,                               &
     &                 THRK_C,QRK_C,EKP14_C,EKP34_C,                    &
     &                 PK_C,PKP1_C,XSBMIN_C,EXK_C,EXKP1_C)
!
!-----------------------------------------------------------------------
!   DECOMPRESS/EXPAND OUTPUT ARRAYS FROM THE DETRAINMENT CALCULATIONS
!-----------------------------------------------------------------------
!
!
        DO I=1,NDET
          THPKP1(INDEX1(I)) = THPKP1_C(I)
          QPKP1(INDEX1(I))  = QPKP1_C(I)
        ENDDO
        DO I=1,NDET
          XSQKP1(INDEX1(I)) = XSQKP1_C(I)
!
          BGMKP1(INDEX1(I)) = BGMKP1_C(I)
        ENDDO
      ENDIF
!
      DO 45 I=1,NPNTS
        DELTAK(I) = 0.0
        THRK(I) = 0.0
        QRK(I) = 0.0
  45  CONTINUE
!
      DO 50 I=1,NDET
        DELTAK(INDEX1(I)) = DELTAK_C(I)
        THRK(INDEX1(I))   = THRK_C(I)
        QRK(INDEX1(I))    = QRK_C(I)
  50  CONTINUE
!L
!L----------------------------------------------------------------------
!L  CALCULATE MASS FLUX AT LEVEL K+1.
!L
!L  UM DOCUMENTATION PAPER 27
!L  SECTION (2B), EQUATION (10A)
!L----------------------------------------------------------------------
!L
      DO 60 I=1,NPNTS
       FLXKP1(I) = FLXK(I)*(1.+EKP14(I))*(1.+EKP34(I))*(1.-DELTAK(I))*  &
     &                                           (1.-AMDETK(I))
       FLXKP12(I)= FLXK(I)*(1.0+EKP14(I))*(1.0-DELTAK(I))*              &
     &                    (1.0-AMDETK(I))
  60  CONTINUE
!L
!L---------------------------------------------------------------------
!L TEST FOR POINTS AT WHICH CONVECTION TERMINATES IN LAYER K+1
!L
!L SUBROUTINE TERM_CON
!L
!L UM DOCUMENTATION PAPER 27
!L SECTION (7)
!L---------------------------------------------------------------------
!L
! DEPENDS ON: term_con
      CALL TERM_CON (NPNTS,NLEV,K,BTERM,BWKP1,FLXKP1,THEKP1,QEKP1,THPI  &
                  ,  QPI,QSEKP1,DELTAK,EXPI,EKP14,EKP34,NEW_TERMC,PSTAR &
                  ,  PK,PKP1,XSBMIN_AD)
!L
!L----------------------------------------------------------------------
!L CLOUD MICROPHYSICS CALCULATION
!L
!L SUBROUTINE CLOUD_W
!L
!L UM DOCUMENTATION PAPER 27
!L SECTION (8), (9)
!L----------------------------------------------------------------------
!L
! DEPENDS ON: cloud_w
      CALL CLOUD_W (K,NPNTS,XPKP1,QCLPKP1,QCFPKP1,PREKP1,XSQKP1,BLOWST, &
     &              FLXKP1,XPK,QCLPK,QCFPK,THEKP1,QEKP1,BWKP1,BLAND,    &
     &     QSEKP1,BGMKP1,BTERM,CCA,ICCB,ICCT,TCW,DEPTH,EKP14,EKP34,     &
     &     DELEXKP1,CCLWP,DELPKP1,CCW,LCCA,LCBASE,LCTOP,LCCLWP,         &
     &     L_SHALLOW,L_Q_INTERACT,START_LEV)
!
      RETURN
      END SUBROUTINE PARCEL
