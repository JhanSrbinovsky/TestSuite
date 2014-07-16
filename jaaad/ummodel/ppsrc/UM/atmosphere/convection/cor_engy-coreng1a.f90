
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE COR_ENGY-----------------------------------------------
!LL
!LL  PURPOSE : TO ADJUST THE POTENTIAL TEMPERATURE INCREMENTS
!LL            TO ENSURE THE CONSERVATION OF MOIST STATIC ENERGY
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  MODEL            MODIFICATION HISTORY:
!LL VERSION  DATE
!LL  4.5  22/7/98  Kill the IBM specific lines (JCThil)
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!    5.3 06/11/01  Correct declaration of 'P_LAYER_BOUNDARIES' S. Cusack
!     5.5  17/04/03 Removal of references to obsolete sections
!                   A05_2A,2C,3B. T.White
!    6.4  12/12/06 Remove def 3C. R A Stratton
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
!LL                  SECTION (12)
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
! cjj     SUBROUTINE COR_ENGY
! cjj    *          (NP_FIELD,NPNTS,NCORE,NLEV,DTHBYDT,DQBYDT,SNOW,
! cjj    *                   EXNER,PSTAR,DELAK,DELBK,AKH,BKH,INDEX4)

      SUBROUTINE COR_ENGY(                                              &
     &           NP_FIELD,NPNTS,NCORE,NLEV,DTHBYDT,DQBYDT,SNOW,         &
     &           EXNER_LAYER_CENTRES, P_LAYER_BOUNDARIES, INDEX4)
!
      IMPLICIT NONE
!
!
!----------------------------------------------------------------------
! MODEL CONSTANTS
!----------------------------------------------------------------------
!
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
!
!----------------------------------------------------------------------
! VECTOR LENGTH AND LOOP COUNTERS
!----------------------------------------------------------------------
!
      INTEGER NP_FIELD            ! LENGTH OF DATA (ALSO USED TO
                                  ! SPECIFY STARTING POINT OF
                                  ! DATA PASSED IN)
!
      INTEGER NCORE               ! IN VECTOR LENGTHS
!
      INTEGER NPNTS               ! IN FULL VECTOR LENGTH
!
      INTEGER NLEV                ! IN NUMBER OF MODEL LAYERS
!
      INTEGER I,K                 ! LOOP COUNTERS
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!----------------------------------------------------------------------
!
      INTEGER INDEX4(NPNTS)
      REAL DQBYDT(NP_FIELD,NLEV)  ! IN INCREMENT TO MODEL MIXING
                                  !    RATIO DUE TO CONVECTION
                                  !    (KG/KG/S)
!
      REAL SNOW(NP_FIELD)         ! IN SNOW AT SURFACE (KG/M**2/S)
!
! cjj      REAL EXNER(NP_FIELD,NLEV+1) ! IN EXNER RATIO
      REAL EXNER_LAYER_CENTRES(NP_FIELD,0:NLEV) ! IN EXNER RATIO
!
! cjj      REAL PSTAR(NP_FIELD)   ! IN SURFACE PRESSURE (PA)
!
! cjj      REAL DELAK(NLEV),      ! IN DIFFERENCE IN HYBRID CO-ORDINATE
! cjj     *     DELBK(NLEV)       !    COEFFICIENTS A AND B
                                  !    ACROSS LAYER K
!
! cjj      REAL AKH(NLEV+1)       ! IN Hybrid coordinate A at
                                  !    layer boundary
! cjj      REAL BKH(NLEV+1)       ! IN Hybrid coordinate B at

! cjj addition.
      REAL P_LAYER_BOUNDARIES(NP_FIELD,0:NLEV)
                                    !    layer boundary
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT AND OUTPUT
!----------------------------------------------------------------------
!
      REAL DTHBYDT(NP_FIELD,NLEV) ! INOUT
                                  ! IN  INCREMENT TO MODEL POTENTIAL
                                  !     TEMPERATURE DUE TO CONVECTION
                                  !     (K/S)
                                  ! OUT CORRECTED INCREMENT TO MODEL
                                  !     POTENTIAL TEMPERATURE DUE TO
                                  !     CONVECTION (K/S)
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE LOCALLY DEFINED
!
      REAL QSUM(NCORE)            ! SUMMATION OF INCREMENTS TO MODEL
                                  ! MIXING RATIO DUE TO CONVECTION
                                  ! IN THE VERTICAL, WEIGHTED
                                  ! ACCORDING TO THE MASS OF THE
                                  ! LAYER (KG/M**2/S)
!
      REAL TSPOS(NCORE)           ! SUMMATION OF POSITIVE INCREMENTS
                                  ! TO MODEL POTENTIAL TEMPERATURE
                                  ! DUE TO CONVECTION WITH HEIGHT,
                                  ! WEIGHTED ACCORDING TO THE MASS
                                  ! OF THE LAYER (K/M**2/S)
!
      REAL TSNEG(NCORE)           ! SUMMATION OF NEGATIVE INCREMENTS
                                  ! TO MODEL POTENTIAL TEMPERATURE
                                  ! DUE TO CONVECTION WITH HEIGHT,
                                  ! WEIGHTED ACCORDING TO THE MASS
                                  ! OF THE LAYER (K/M**2/S)
!
      REAL TERR(NCORE)            ! SUMMATION OF ALL INCREMENTS TO
                                  ! MODEL POTENTIAL TEMPERATURE
                                  ! DUE TO CONVECTION WITH HEIGHT,
                                  ! WEIGHTED ACCORDING TO THE MASS
                                  ! OF THE LAYER (K/M**2/S)
!
      LOGICAL BPOSER(NCORE)       ! MASK FOR POINTS IN LAYER K AT WHICH
                                  ! INCREMENTS TO MODEL POTENTIAL
                                  ! TEMPERATURE DUE TO CONVECTION ARE
                                  ! POSITIVE
!
      LOGICAL BCORR(NCORE)        ! MASK FOR POINTS AT WHICH ENTHALPY
                                  ! CORRECTION IS NECESSARY
!
      REAL DELPK                  ! DIFFERENCE IN PRESSURE ACROSS A
                                  ! LAYER (PA)
!
      REAL EXTEMPK                ! EXNER RATIO AT THE MID-POINT OF
                                  ! LAYER K
!

      REAL                                                              &
     &    PU,PL
! cjj *CALL P_EXNERC

!*----------------------------------------------------------------------
!L
!L----------------------------------------------------------------------
!L  SUM UP MIXING RATIO AND +VE AND -VE TEMPERATURE INCREMENTS
!L----------------------------------------------------------------------
!L
      DO 20 I=1,NCORE
       QSUM (I) = 0.0
       TSPOS(I) = 0.0
       TSNEG(I) = 0.0
   20  CONTINUE
!
      DO 40 K=1,NLEV
       DO 30 I=1,NCORE
!
! cjj        DELPK = -DELAK(K) - DELBK(K)*PSTAR(INDEX4(I))
        DELPK = - (P_LAYER_BOUNDARIES(INDEX4(I),K)  -                   &
     &          P_LAYER_BOUNDARIES(INDEX4(I),K-1))
!
! cjj        PU=PSTAR(INDEX4(I))*BKH(K+1) + AKH(K+1)
! cjj        PL=PSTAR(INDEX4(I))*BKH(K) + AKH(K)
! cjj        EXTEMPK  =
! cjj     &    P_EXNER_C(EXNER(INDEX4(I),K+1),
! cjj     &              EXNER(INDEX4(I),K),PU,PL,KAPPA)
        EXTEMPK  = EXNER_LAYER_CENTRES(INDEX4(I),K)
!
!
        QSUM(I) = QSUM(I) + DQBYDT(INDEX4(I),K)*DELPK
!
        IF (DTHBYDT(INDEX4(I),K)  >   0.0) THEN
           TSPOS(I) = TSPOS(I) +                                        &
     &                DTHBYDT(INDEX4(I),K)*(CP*DELPK*EXTEMPK)
        ELSE
           TSNEG(I) = TSNEG(I) +                                        &
     &                DTHBYDT(INDEX4(I),K)*(CP*DELPK*EXTEMPK)
        ENDIF
   30  CONTINUE
   40 CONTINUE
!L
!L----------------------------------------------------------------------
!L  CALCULATE THE ERROR AND APPLY THE NECESSARY CORRECTION
!L
!L  UM DOCUMENTATION PAPER 27
!L  SECTION (12), EQUATION (48), (49)
!L----------------------------------------------------------------------
!L
      DO 50 I=1,NCORE
!
       TERR(I) = LC*QSUM(I) - LF*G*SNOW(INDEX4(I)) +                    &
     &                                   TSPOS(I) + TSNEG(I)
!
       BPOSER(I) = TERR(I)  >   0.0
!
       IF (BPOSER(I) .AND. (TSPOS(I)  ==  0.0)) THEN
          BPOSER(I) = .FALSE.
       ELSE IF (.NOT.BPOSER(I) .AND. (TSNEG(I)  ==  0.0)) THEN
          BPOSER(I) = .TRUE.
       ENDIF
!
       BCORR(I) = (TSPOS(I)  /=  0.0) .OR. (TSNEG(I)  /=  0.0)
!
       IF (BPOSER(I) .AND. BCORR(I)) THEN
          TERR(I) = 1. - TERR(I)/TSPOS(I)
       ELSE IF (.NOT.BPOSER(I) .AND. BCORR(I)) THEN
          TERR(I) = 1. - TERR(I)/TSNEG(I)
       ENDIF
!
  50  CONTINUE
!
      DO K=1,NLEV
!DIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
       DO I=1,NCORE
        IF (BCORR(I) .AND. (( BPOSER(I) .AND.                           &
     &   (DTHBYDT(INDEX4(I),K)  >   0.0)) .OR. ( .NOT.BPOSER(I)         &
     &   .AND. (DTHBYDT(INDEX4(I),K)  <   0.0))))                       &
     &       DTHBYDT(INDEX4(I),K) = DTHBYDT(INDEX4(I),K)*TERR(I)
       ENDDO  ! NCORE
      ENDDO ! NLEV
!
      RETURN
      END SUBROUTINE COR_ENGY
