
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE LIFT_PAR-----------------------------------------------
!LL
!LL  PURPOSE : LIFTS THE PARCEL FROM LAYER K TO K+1
!LL            TAKING ENTRAINEMNT AND MOIST PROCESSES INTO ACOUNT
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
!LL  VERSION NO. 1
!LL
!LL  LOGICAL COMPONENTS COVERED: P27
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER NO. ##
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE LIFT_PAR (NPNTS,NP_FULL,                              &
                          THPKP1,QPKP1,XSQKP1,BGMKP1,BWKP1,BWK,        &
                          THPK,QPK,XPK,THEKP1,QEKP1,THEK,QEK,QSEKP1,   &
                          QCLPKP1, QCLPK, QCLEKP1, QCLEK, L_Q_INTERACT,&
                          QCFPKP1, QCFPK, QCFEKP1, QCFEK,              &
                          PK,PKP1,EXKP1,EKP14,EKP34,L_MOM_GK,UPKP1,    &
                          VPKP1,UPK,VPK,UEK,UEKP1,VEK,VEKP1,L_TRACER,  &
                          NTRA,TRAPKP1,TRAPK,TRAEKP1,TRAEK)
!
      IMPLICIT NONE
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
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
! RV gas constant for water vapour (J/kg/K)
      REAL,PARAMETER:: RV = 461.1
! RV end
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
!----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!----------------------------------------------------------------------
!
      INTEGER NPNTS          ! IN VECTOR LENGTH
!
      INTEGER NP_FULL        ! IN FULL VECTOR LENGTH
!
      INTEGER I,KTRA         ! LOOP COUNTERS
!
      INTEGER NTRA           ! IN NUMBER OF TRACER VARIABLES

!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
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
      REAL QCLEK(NPNTS)      ! IN LIQUID CONDENSATE MIXING RATIO OF
!                                 CLOUD ENVIRONMENT IN LAYER K (KG/KG)
!
      REAL QCLEKP1(NPNTS)    ! IN LIQUID CONDENSATE MIXING RATIO OF
!                                 CLOUD ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL QCFEK(NPNTS)      ! IN FROZEN CONDENSATE MIXING RATIO OF
!                                 CLOUD ENVIRONMENT IN LAYER K (KG/KG)
!
      REAL QCFEKP1(NPNTS)    ! IN FROZEN CONDENSATE MIXING RATIO OF
!                                 CLOUD ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL UEK(NPNTS)        ! IN U OF ENVIRONMENT IN LAYER K (M/S)
!
      REAL UEKP1(NPNTS)      ! IN U OF ENVIRONMENT IN LAYER K+1 (M/S)
!
      REAL VEK(NPNTS)        ! IN V OF ENVIRONMENT IN LAYER K (M/S)
!
      REAL VEKP1(NPNTS)      ! IN V OF ENVIRONMENT IN LAYER K+1 (M/S)
!
      REAL TRAEK(NP_FULL,                                               &
                             ! IN TRACER CONTENT OF CLOUD
     &           NTRA)       !    ENVIRONMENT IN LAYER K (KG/KG)
!
      REAL TRAEKP1(NP_FULL,                                             &
                             ! IN TRACER CONTENT OF CLOUD
     &             NTRA)     !    ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL QSEKP1(NPNTS)     ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL THPK(NPNTS)       ! IN PARCEL POTENTIAL TEMPERATURE IN
                             !    LAYER K (K)
!
      REAL QPK(NPNTS)        ! IN PARCEL MIXING RATIO IN LAYER K (KG/KG)
      REAL XPK(NPNTS)        ! IN PARCEL CONDENSATE IN LAYER K (KG/KG)
!
      REAL QCLPK(NPNTS)      ! IN PARCEL LIQUID CONDENSATE MIXING RATIO
!                                 IN LAYER K (KG/KG)
!
      REAL QCFPK(NPNTS)      ! IN PARCEL FROZEN CONDENSATE MIXING RATIO
!                                 IN LAYER K (KG/KG)
!
      REAL UPK(NPNTS)        ! IN PARCEL U IN LAYER K (M/S)
!
      REAL VPK(NPNTS)        ! IN PARCEL V IN LAYER K (M/S)
!
      REAL TRAPK(NP_FULL,                                               &
                             ! IN PARCEL TRACER CONTENT IN LAYER K
     &           NTRA)       !    (KG/KG)
!
      LOGICAL BWKP1(NPNTS)   ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
      LOGICAL BWK(NPNTS)     ! IN Ditto K
!
      REAL PKP1(NPNTS)       ! IN PRESSURE AT LEVEL K+1 (PA)
!
      REAL PK(NPNTS)         ! IN Pressure at centre of layer K
!
      REAL EXKP1(NPNTS)      ! IN EXNER RATIO AT MID-POINT OF LAYER K+1
!
      REAL EKP14(NPNTS)      ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+1/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
!
      REAL EKP34(NPNTS)      ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+3/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
!

      LOGICAL L_TRACER       ! IN LOGICAL SWITCH FOR INCLUSION OF
                             !    TRACERS

      LOGICAL L_MOM_GK       ! IN LOGICAL SWITCH FOR INCLUSION OF
                             !    MOMENTUM TRANSPORT (Gregory-Kershaw)
!
      LOGICAL L_Q_INTERACT   ! IN SWITCH FOR OVERWRITING OF PARCEL VARS
!                                 IF CALCULATING CONDENSATE INCREMENTS.
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT BUT WHICH ARE ALSO UPDATED IN THIS ROUTINE
!----------------------------------------------------------------------
!
      LOGICAL BGMKP1(NPNTS)  ! INOUT
                             ! IN  Is parcel saturated in layer K ?
                             ! OUT Is parcel saturated in layer K+1 ?
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE OUTPUT
!----------------------------------------------------------------------
!
      REAL THPKP1(NPNTS)     ! OUT PARCEL POTENTIAL TEMPERATURE IN
                             !     LAYER K+1 AFTER ENTRAINMENT AND
                             !     LATENT HEATING (K)
!
      REAL QPKP1(NPNTS)      ! OUT PARCEL MIXING RATIO IN LAYER K+1
                             !     AFTER ENTRAINMENT AND LATENT HEATING
                             !     (KG/KG)
!
      REAL QCLPKP1(NPNTS)    ! OUT PARCEL LIQUID CONDENSATE MIXING RATIO
!                                  IN LAYER K+1 AFTER ENTRAINMENT
!
      REAL QCFPKP1(NPNTS)    ! OUT PARCEL FROZEN CONDENSATE MIXING RATIO
!                                  IN LAYER K+1 AFTER ENTRAINMENT
!
      REAL UPKP1(NPNTS)      ! OUT PARCEL U IN LAYER K+1 AFTER
                             !     ENTRAINMENT (M/S)
!
      REAL VPKP1(NPNTS)      ! OUT PARCEL V IN LAYER K+1 AFTER
                             !     ENTRAINMENT (M/S)
!
      REAL TRAPKP1(NP_FULL,                                             &
                             ! OUT PARCEL TRACER CONTENT IN LAYER
     &             NTRA)     !     K+1 AFTER ENTRAINMENT. (KG/KG)
!
      REAL XSQKP1(NPNTS)     ! OUT EXCESS PARCEL WATER AFTER
                             !     LIFTING FROM LAYER K TO K+1
                             !     (KG/KG)
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!
! ON THE IBM ARRAYS ARE ALLOCATED USING A PARAMETER STATEMENT
!
! ON THE CRAY ARRAYS ARE DYNAMICALLY ALLOCATED
!----------------------------------------------------------------------
!
      REAL THPKP1T(NPNTS)    ! INITIAL ESTIMATE OF PARCEL TEMPERATURE
                             ! IN LAYER K+1 AFTER ENTRAINMENT (K)
!
      REAL TT(NPNTS)         ! TEMPORARY TEMPERATURE USED IN CALCULATION
                             ! OF SATURATION MIXING RATIO (K)
!
      REAL PT(NPNTS)         ! TEMPORARY PRESSURE USED IN CALCULATION
                             ! OF SATURATION MIXING RATIO (K)
!
      REAL QSPKP1(NPNTS)     ! SATURATION MIXING RATIO OF PARCEL
                             ! AFTER DRY ASCENT (KG/KG)
!
      REAL DQSDT,                                                       &
                             ! Rate of change of qsat with temperature
     &     EL                ! Latent heat of gas-to-whatever-condenses
                             ! PC2 defn
                             ! LATENT HEAT OF CONDENSATION OR
                             ! (CONDENSATION + FUSION) (J/KG)
!
      REAL REPSS             ! 1 / EPSS(K)
!
      REAL LFBYCP
      PARAMETER ( LFBYCP = LF / CP )
!----------------------------------------------------------------------
!
      DO I=1,NPNTS
!L
!L----------------------------------------------------------------------
!L  LIFT PARCEL MIXING RATIO, POTENTIAL TEMPERATURE, U, V AND TRACER
!L  TO THE NEXT LEVEL
!L----------------------------------------------------------------------
!L
!L----------------------------------------------------------------------
!L  INITIAL 'DRY' ASCENT
!L
!L  UM DOCUMENTATION PAPER P27
!L  SECTION (3), EQUATIONS (11B), (12B)
!L----------------------------------------------------------------------
!L
       REPSS = 1. / ((1.+EKP14(I))*(1.+EKP34(I)))
!
       THPKP1(I) = (  THPK(I)                                           &
     &             + EKP14(I)*THEK(I) + EKP34(I)*(1.+EKP14(I))*THEKP1(I)&
     &             ) / ((1.+EKP14(I))*(1.+EKP34(I)))
!
       QPKP1(I) = (  QPK(I)                                             &
     &             + EKP14(I)*QEK(I) + EKP34(I)*(1.+EKP14(I))*QEKP1(I)  &
     &             ) / ((1.+EKP14(I))*(1.+EKP34(I)))
!
      IF (L_Q_INTERACT) THEN
!
       QCLPKP1(I)= ( QCLPK(I)                                           &
     &           + EKP14(I)*QCLEK(I) + EKP34(I)*(1.+EKP14(I))*QCLEKP1(I)&
     &              ) * REPSS
!
       QCFPKP1(I)= ( QCFPK(I)                                           &
     &           + EKP14(I)*QCFEK(I) + EKP34(I)*(1.+EKP14(I))*QCFEKP1(I)&
     &              ) * REPSS
!
! ----------------------------------------------------------------------
! IMPORTANT NOTE: The condensate term has the same form as those for the
! vapour and temperature in the initial ascent. IF we had a clear idea
! on how to parametrize detrainment and mixing of condensate, it could
! be updated alongside the vapour process-by-process within convection.
! HOWEVER, we opt instead simply to STORE the humidity excess over Qsat
! in XSQKP1 as the vapour is updated and apply moisture conservation
! just before the increments to the environment variables are calculated
! (ie. XSQKP1 is added to QCLPKP1 or QCFPKP1 in CLOUDW1A).
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
!       CURRENTLY MIXED PHASE PARCEL IS FORBIDDEN. MELT OR FREEZE THE
!       ENTRAINED LAYER CLOUD AND ADJUST PARCEL TEMPERATURE ACCORDINGLY.
! ----------------------------------------------------------------------
!
          IF (BWKP1(I)) THEN
            QCLPKP1(I) = QCLPKP1(I) + QCFPKP1(I)
            THPKP1(I)  = ( THPKP1(I) - ( MAX( 0.0, QCFPKP1(I) ) * LF /  &
     &                                          (CP * EXKP1(I)) ) )
            QCFPKP1(I) = 0.0
          ELSE

            QCFPKP1(I) = QCLPKP1(I) + QCFPKP1(I)
            THPKP1(I)  = ( THPKP1(I) + ( MAX( 0.0, QCLPKP1(I) ) * LF /  &
     &                                          (CP * EXKP1(I)) ) )
            QCLPKP1(I) = 0.0
          END IF
!
      ELSE ! L_q_interact_if1
!
!
! Allow for release of latent heat of fusion if crossing the freezing
!   level.
        IF ( BWKP1(I) .AND. .NOT. BWK(I) ) THEN
          THPKP1(I) = THPKP1(I) - LFBYCP * XPK(I) / EXKP1(I)
        ELSE IF ( BWK(I) .AND. .NOT. BWKP1(I) ) THEN
          THPKP1(I) = THPKP1(I) + LFBYCP * XPK(I) / EXKP1(I)
        ENDIF
!
      ENDIF  ! L_q_interact_if1
!
       END DO
!
      IF (L_MOM_GK) THEN   ! l_mom_gk set according to type of CMT required
                           ! This code does Gregory-Kershaw CMT
!
         DO I=1,NPNTS
!
! CALCULATIONS FOR CUMULUS CONVECTION DONE AFTER TERMINATION
! 
           UPKP1(I) = (UPK(I)+                                            &
                      EKP14(I)*UEK(I) + EKP34(I)*(1.0+EKP14(I))*UEKP1(I)  &
                        ) / ((1.0+EKP14(I))*(1.0+EKP34(I)))

           VPKP1(I) = (VPK(I)+                                            &
                      EKP14(I)*VEK(I) + EKP34(I)*(1.0+EKP14(I))*VEKP1(I)  &
                        ) / ((1.0+EKP14(I))*(1.0+EKP34(I)))
           UPKP1(I) = UPKP1(I) - (0.7*(UEK(I)-UEKP1(I))/(1.0+EKP34(I)))
           VPKP1(I) = VPKP1(I) - (0.7*(VEK(I)-VEKP1(I))/(1.0+EKP34(I)))

         END DO

       END IF       ! l_mom_gk test
!
       IF(L_TRACER)THEN
!
       DO KTRA = 1,NTRA
       DO I = 1,NPNTS
!
       TRAPKP1(I,KTRA) = ( TRAPK(I,KTRA)                                &
     & + EKP14(I)*TRAEK(I,KTRA) + EKP34(I)*(1.+EKP14(I))*TRAEKP1(I,KTRA)&
     &             ) / ((1.+EKP14(I))*(1.+EKP34(I)))
!
       END DO
       END DO
!
       END IF
!
!-----------------------------------------------------------------------
!   CALCULATE WHERE THE PARCEL IS SUPERSATURATED (IE WHERE GAMMA(K+1)=1
!   SEE DCTN 29 PAGE 123)
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
! CONVERT POTENTIAL TEMPERATURE TO TEMPERATURE AND CALCULATE
! PRESSURE OF LAYER K FOR CALCULATION OF SATURATED
! MIXING RATIO
!-----------------------------------------------------------------------
!
      DO I=1,NPNTS
!
! For parcel saturated in layer K, use Exner pressure in K+1 with
! theta for layer K (gives partly compensating errors)
! For unsaturated parcel in layer K, use the K+1 to test whether
! saturation is occurring there
!
        IF ( BGMKP1(I) ) THEN
          TT(I) = THPK(I) * EXKP1(I)
        ELSE
          TT(I) = THPKP1(I) * EXKP1(I)
        ENDIF
      END DO
!
! DEPENDS ON: qsat
      CALL QSAT (QSPKP1,TT,PKP1,NPNTS)
!
      DO I=1, NPNTS

        BGMKP1(I) = BGMKP1(I) .OR. ( QSPKP1(I)  <   QPKP1(I) )

        IF ( BGMKP1(I) ) THEN
! ----------------------------------------------------------------------
!       ADJUST PARCEL TEMPERATURE TO ACCOUNT FOR LATENT HEATING
! ----------------------------------------------------------------------
          IF ( BWKP1(I) ) THEN
             EL = LC
           ELSE
             EL = LC + LF
          ENDIF
!
          DQSDT = EL * QSPKP1(I) / ( RV * TT(I)*TT(I) )
!
! Consistent with QSAT calculation above:
          THPKP1T(I) = THPKP1(I) + ( PKP1(I) - PK(I) ) * KAPPA * TT(I)  &
     &                * DQSDT  * ( CP * TT(I) - EPSILON * EL ) /        &
     &          ( EPSILON * exkp1(i)*( CP + EL * DQSDT ) * PKP1(I) )
!
! ----------------------------------------------------------------------
!       CALCULATE T FROM THETA
! ----------------------------------------------------------------------
          TT(I) = THPKP1T(I) * EXKP1(I)
        ENDIF
      ENDDO

! ----------------------------------------------------------------------
!     CALCULATE A MORE ACCURATE PARCEL SATURATED MIXING RATIO
! ----------------------------------------------------------------------
!
! DEPENDS ON: qsat
      CALL QSAT (QSPKP1, TT, PKP1, NPNTS)

      DO I=1,NPNTS
        IF ( BGMKP1(I) ) THEN
          IF ( BWKP1(I) ) THEN
             EL = LC
           ELSE
             EL = LC + LF
          ENDIF
          DQSDT = EL * QSPKP1(I) / ( RV * TT(I)*TT(I) )
          THPKP1T(I) = ( THPKP1(I) * CP   +                             &
     &    ((QPKP1(I) - QSPKP1(I))/exkp1(i) + THPKP1T(I) * DQSDT ) * EL) &
     &      /( CP  + DQSDT * EL )


          XSQKP1(I) = ( THPKP1T(I) - THPKP1(I) ) * CP * EXKP1(I) / EL
!
          BGMKP1(I) = ( XSQKP1(I)  >   0.0 ) .AND. (                    &
     &                ( QPKP1(I) - XSQKP1(I) )  >   0.85 * QSPKP1(I))
        ENDIF
      ENDDO
      DO I=1, NPNTS
        IF ( BGMKP1(I) ) THEN
           THPKP1(I) = THPKP1T(I)
           QPKP1(I)  = QPKP1(I) - XSQKP1(I)
         ELSE
           XSQKP1(I) = 0.
        ENDIF
      ENDDO
!
      RETURN
      END SUBROUTINE LIFT_PAR
