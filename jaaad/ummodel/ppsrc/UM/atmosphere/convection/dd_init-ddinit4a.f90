
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DD_INIT------------------------------------------------
!LL
!LL  PURPOSE : ROUTINE TO INITIALISE THE DOWNDRAUGHT
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
!LL VERSION  DATE
!LL   5.4   6/8/02    New deck created for version 4A of convection
!LL                   scheme. No CMT by downdraughts in this version.
!LL                                                       Gill Martin
!     6.2   03/02/05  Added section 5A R.A.Stratton
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL  VERSION NO. 4  DATED 5/2/92
!LL
!LL  LOGICAL COMPONENTS COVERED:
!LL
!LL  SYSTEM TASK : P27
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE DD_INIT(NPNTS,NP_FULL,TH_UD_K,Q_UD_K,THE_K,QE_K,PK,    &
     &                   EXK,THDD_K,QDD_K,DELTD,DELQD,BDD_START,K,BDDI, &
     &                   BDD_ON,                                        &
     &                   L_TRACER,NTRA,TRA_UD_K,                        &
     &                   TRAE_K,TRADD_K,DELTRAD)
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
! MODEL CONSTANTS
!-----------------------------------------------------------------------
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
!-----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!-----------------------------------------------------------------------
!
!
      INTEGER I,KTRA            ! LOOP COUNTERS
!
      INTEGER NPNTS             ! VECTOR LENGTH
!
      INTEGER NP_FULL           ! FULL VECTOR LENGTH
!
      INTEGER NTRA              ! NUMBER OF TRACER VARIABLES
!
      INTEGER K                 ! IN PRESENT MODEL LAYER
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!-----------------------------------------------------------------------
!
      REAL THE_K(NPNTS)         ! IN POTENTIAL TEMPERATURE OF
                                !    ENVIRONMENT IN LAYER K (K)
!
      REAL TH_UD_K(NPNTS)       ! IN PARCEL POTENTIAL TEMPERATURE OF
                                !    UPDRAUGHT, LAYER K (K)
!
      REAL QE_K(NPNTS)          ! IN MIXING RATIO OF ENVIRONMENT IN
                                !    LAYER K (KG/KG)
!
      REAL Q_UD_K(NPNTS)        ! IN PARCEL MIXING RATIO OF UPDRAUGHT,
                                !    LAYER K (KG/KG)
!
      REAL TRAE_K(NP_FULL,NTRA) ! IN TRACER CONTENT OF ENVIRONMENT
                                !    IN LAYER K (KG/KG)
!
      REAL TRA_UD_K(NP_FULL,                                            &
                                ! IN PARCEL TRACER CONTENT OF
     &              NTRA)       !    UPDRAUGHT IN LAYER K (KG/KG)
!
      REAL EXK(NPNTS)           ! IN EXNER RATIO OF LAYER K
!
      REAL PK(NPNTS)            ! IN PRESSURE OF LAYER K (PA)
!
      LOGICAL BDDI(NPNTS)       ! IN MASK FOR THOSE POINTS WHERE
                                !    DOWNDRAUGHT MAY INITIATE
!
      LOGICAL BDD_ON(NPNTS)     ! IN MASK FOR THOSE POINTS WHERE
                                !    DOWNDRAUGHT IS ON
!
      LOGICAL L_TRACER          ! IN SWITCH FOR INCLUSION OF TRACERS
!
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT AND OUTPUT
!-----------------------------------------------------------------------
!
      LOGICAL BDD_START(NPNTS)  ! INOUT
                                ! IN  MASK FOR THOSE POINT WHERE
                                !     DOWNDRAUGHT MAY START
                                ! OUT MASK FOR THOSE POINTS WHERE
                                !
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE OUTPUT
!-----------------------------------------------------------------------
!
      REAL THDD_K(NPNTS)        ! OUT DOWNDRAUGHT POTENTIAL TEMPERATURE
                                !     OF LAYER K
!
      REAL QDD_K(NPNTS)         ! OUT DOWNDRAUGHT MIXING RATIO OF
                                !     LAYER K
!
      REAL TRADD_K(NP_FULL,                                             &
                                ! OUT DOWNDRAUGHT TRACER CONTENT OF
     &             NTRA)        !     LAYER K
!
      REAL DELTD(NPNTS)         ! OUT COOLING NECESSARY TO ACHIEVE
                                !     SATURATION
!
      REAL DELQD(NPNTS)         ! OUT MOISTENING NECESSARY TO ACHIEVE
                                !     SATURATION
!
      REAL DELTRAD(NP_FULL,NTRA)! OUT DEPLETION OF ENVIRONMENT TRACER
                                !     DUE TO FORMATION OF DOWNDRAUGHT
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!-----------------------------------------------------------------------
!
!
      REAL TH_MEAN(NPNTS)       ! MEAN POTENTIAL TEMPERATURE USED IN
                                ! CALCULATION OF SATURATED DOWNDRAUGHT
                                ! POTENTIAL TEMPERATURE IN LAYER K
!
      REAL Q_MEAN(NPNTS)        ! MEAN MIXING RATIO USED IN CALCULATION
                                ! OF SATURATED DOWNDRAUGHT
                                ! MIXING RATIO FOR LAYER K
!
!
      REAL TRA_MEAN(NPNTS,NTRA) ! MEAN TRACER USED AS INITIAL TRACER
                                ! CONTENT OF DOWNDRAUGHT IN LAYER K
                                ! (KG/KG)
!
      REAL T_MEAN(NPNTS)        ! MEAN TEMPERATURE USED IN CALCULATION
                                ! OF SATURATED DOWNDRAUGHT POTENTIAL
                                ! TEMPERATURE OF LAYER K (K)
!
      REAL THDDS(NPNTS)         ! SATURATED DOWNDRAUGHT POTENTIAL
                                ! TEMPERATURE IN LAYER K (K)
!
      REAL QDDS(NPNTS)          ! SATURATED DOWNDRAUGHT MIXING RATIO
                                ! IN LAYER K (KG/KG)
!
      REAL BUOY(NPNTS)          ! BUOYANCY OF PARCEL IN LAYER K
!
!
      REAL THDD_V               ! VIRTUAL POTENTIAL TEMPERATURE OF
                                ! PARCEL IN LAYER K
!
      REAL THE_V                ! VIRTUAL POTENTIAL TEMPERATURE OF
                                ! ENVIRONMENT IN LAYER K
!
!-----------------------------------------------------------------------
! EXTERNAL ROUTINES CALLED
!-----------------------------------------------------------------------
!
      EXTERNAL SATCAL
!
!-----------------------------------------------------------------------
! CALCULATE MEAN TEMPERATURE, MIXING RATIO, U, V AND TRACER
!-----------------------------------------------------------------------
!
      DO I=1,NPNTS
       TH_MEAN(I) = (THE_K(I)+TH_UD_K(I))*0.5
       Q_MEAN(I) = (QE_K(I)+Q_UD_K(I))*0.5
       T_MEAN(I) = TH_MEAN(I)*EXK(I)
      END DO
!
      IF(L_TRACER)THEN
!
      DO KTRA=1,NTRA
        DO I=1,NPNTS
          TRA_MEAN(I,KTRA) = (TRAE_K(I,KTRA)+TRA_UD_K(I,KTRA))*0.5
        END DO
      END DO
!
      END IF
!
!
!-----------------------------------------------------------------------
! CALCULATE SATURATED DOWNDRAUGHT POTENTIAL TEMPERATURE FOR LAYER K
!-----------------------------------------------------------------------
!
! DEPENDS ON: satcal
      CALL SATCAL(NPNTS,T_MEAN,TH_MEAN,PK,QDDS,THDDS,K,EXK,Q_MEAN,      &
     &            THE_K)
!
!-----------------------------------------------------------------------
! IS SATURATED PARCEL NEGATIVELY BUOYANT COMPARED TO ENVIRONMENT
!-----------------------------------------------------------------------
!
      DO I=1,NPNTS
       IF (.NOT. BDD_ON(I) .AND. BDDI(I) ) THEN
          THDD_V = THDDS(I)*(1.0+C_VIRTUAL*QDDS(I))
          THE_V = THE_K(I)*(1.0+C_VIRTUAL*QE_K(I))
          BUOY(I) = THDD_V - THE_V
!
          IF (BUOY(I)  <   0.5 ) THEN
!
!-----------------------------------------------------------------------
! INITIATE DOWNDRAUGHT
!-----------------------------------------------------------------------
!
             THDD_K(I) = THDDS(I)
             QDD_K(I) = QDDS(I)
             BDD_START(I) = .TRUE.
!
!-----------------------------------------------------------------------
! CALCULATE COOLING AND MOISTENING TO ACHIEVE SATURATION
!-----------------------------------------------------------------------
!
             DELTD(I) = THDDS(I)-THE_K(I)
             DELQD(I) = QDDS(I)-QE_K(I)
          END IF
       END IF
      END DO
!
!
!
      IF(L_TRACER)THEN
!
        DO KTRA=1,NTRA
          DO I=1,NPNTS
            IF(.NOT.BDD_ON(I).AND.BDDI(I).AND.K >= 4)THEN
              IF(BUOY(I) <  0.5)THEN
              TRADD_K(I,KTRA) = TRA_MEAN(I,KTRA)
              DELTRAD(I,KTRA) = TRADD_K(I,KTRA)-TRAE_K(I,KTRA)
              END IF
            END IF
          END DO
        END DO
!
      END IF
      RETURN
      END SUBROUTINE DD_INIT
!
