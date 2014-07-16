#if defined(A05_4A) || defined(A05_5A)
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
      SUBROUTINE EVAP_BCB_NODD(NP_FIELD,NPNTS,N_NODD,KCT,THE,QE,DTHBYDT,&
     &                         DQBYDT,EXNER_LAYER_CENTRES,              &
     &                         EXNER_LAYER_BOUNDARIES,RAIN,SNOW,        &
     &                         RAIN_3D, SNOW_3D, PRECIP,                &
     &                         PSTAR,P_LAYER_CENTRES,P_LAYER_BOUNDARIES,&
     &                         BWATER,CCA,ICCB,TIMESTEP,INDEX1)
!
      IMPLICIT NONE
!
!
!     Purpose: To calculate the convective precipitation reaching the
!              surface.
!
!     Method : The evaporation below cloud base follows that done in
!              the downdraught routine for the environmental part of the
!              column. The points which are gathered here are those
!              points which have an updraught, but no downdraught.
!
!
!     Compatible with version 3C of CONVEC3C
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   System component covered: P292
!
!   Documentation:
!
!
      INTEGER KCT                ! IN CONVECTIVE CLOUD TOP LAYER
!
      INTEGER N_NODD             ! IN COMPRESSED VECTOR LENGTH FOR
                                 ! CALCULATION
!
      INTEGER NP_FIELD           ! IN FULL VECTOR LENGTH
!
      INTEGER NPNTS              ! IN VECTOR LENGTH FOR SOME ARRAYS
!
!
      REAL P_LAYER_CENTRES(NP_FIELD,0:KCT+1)
                                 ! IN PRESSURE AT LAYER CENTRES
      REAL P_LAYER_BOUNDARIES(NP_FIELD,0:KCT+1)
                                 ! IN PRESSURE AT LAYER BOUNDARIES
      REAL EXNER_LAYER_CENTRES(NP_FIELD,0:KCT+1)
                                 ! IN EXNER FUNCTION AT LAYER CENTRES
                                 !    STARTING AT LEVEL K-1/2
!
      REAL EXNER_LAYER_BOUNDARIES(NP_FIELD,0:KCT+1)
                                 ! IN EXNER FUNCTION AT LAYER BOUNDARIES
                                 !    STARTING AT LEVEL K-1/2
!
      REAL THE(NP_FIELD,KCT+1)   ! IN MODEL ENVIRONMENTAL POTENTIAL
                                 !    TEMPERATURE (K)
!
      REAL QE(NP_FIELD,KCT+1)    ! IN ENVIRONMENT MIXING RATIO
                                 !    (KG/KG)
!
      REAL PSTAR(NP_FIELD)       ! IN SURFACE PRESSURE (PA)
!
      REAL PRECIP(NPNTS,KCT+1)   ! IN PRECIPITATION ADDED WHEN
                                 !    DESCENDING FROM LAYER K TO K-1
                                 !    (KG/M**2/S)
!
      INTEGER ICCB(NP_FIELD)     ! IN CLOUD BASE LEVEL
!
      REAL CCA(NP_FIELD)         ! IN CONVECTIVE CLOUD AMOUNT
!
      REAL TIMESTEP
!
      INTEGER INDEX1(NPNTS)      ! INDEX OF DOWNDRAUGHTS NOT POSSIBLE
!
      LOGICAL BWATER(NPNTS,2:KCT+1)!IN  MASK FOR THOSE POINTS AT WHICH
                                   !    CONDENSATE IS WATER IN LAYER K
!
      REAL RAIN(NP_FIELD)          ! OUT RAINFALL AT SURFACE (KG/M**2/S)
!
      REAL SNOW(NP_FIELD)          ! OUT SNOWFALL AT SURFACE (KG/M**2/S)

      REAL RAIN_3D(NP_FIELD,kct+1) ! OUT RAINFALL FLUX (KG/M**2/S)
!
      REAL SNOW_3D(NP_FIELD,kct+1) ! OUT SNOWFALL FLUX (KG/M**2/S)
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT AND OUTPUT
!-----------------------------------------------------------------------
!
      REAL DTHBYDT(NP_FIELD,KCT+1) ! INOUT
                                   ! IN  INCREMENT TO MODEL POTENTIAL
                                   !     TEMPERATURE (K/S)
                                   ! OUT UPDATED INCREMENT TO MODEL
                                   !     POTENTIAL TEMPERATURE (K/S)
!
      REAL DQBYDT(NP_FIELD,KCT+1)  ! INOUT
                                   ! IN  INCREMENT TO MODEL MIXING
                                   !     RATIO (KG/KG/S)
                                   ! OUT UPDATED INCREMENT TO MODEL
                                   !     MIXING RATIO (KG/KG/S)
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!-----------------------------------------------------------------------
!
      INTEGER I,K                ! LOOP COUNTERS
!
!
      REAL EXNER_KM12_C(N_NODD)  ! COMPRESSED EXNER FUNCTION AT
                                 ! LAYER K
!
      REAL EXNER_KP12_C(N_NODD)  ! COMPRESSED EXNER FUNCTION AT
                                 ! LAYER K+1
!
      REAL EXNER_KM32_C(N_NODD)  ! COMPRESSED EXNER FUNCTION AT
                                 ! LAYER K-1
!
      REAL PK(N_NODD)            ! PRESSURE OF LAYER K (PA)
!
      REAL PKM1_C(N_NODD)        ! PRESSURE OF LAYER K-1 (PA)
!
      REAL EXK_C(N_NODD)         ! EXNER RATIO FOR LAYER K
!
      REAL EXKM1_C(N_NODD)       ! EXNER RATIO FOR LAYER K-1
!
      REAL DELPKM1_C(N_NODD)     ! PRESSURE DIFFERENCE ACROSS
                                 ! LAYER K-1 (PA)
!
      REAL PRECIP_K_C(N_NODD)    ! COMPRESSED PRECIPITATION
                                 ! ADDED WHEN DESCENDING FROM
                                 ! LAYER K TO K-1 (KG/M**2/S)
!
      REAL QE_K_C(N_NODD)        ! COMPRESSED PARCEL MIXING RATIO
                                 ! OF LAYER K (KG/KG)
!
      REAL THE_K_C(N_NODD)       ! COMPRESSED PARCEL POTENTIAL
                                 ! TEMPERATURE OF LAYER K (K)
!
      REAL QE_KM1_C(N_NODD)      ! COMPRESSED PARCEL MIXING RATIO
                                 ! OF LAYER K (KG/KG)
!
      REAL THE_KM1_C(N_NODD)     ! COMPRESSED PARCEL POTENTIAL
                                 ! TEMPERATURE OF LAYER K (K)
!
      REAL PSTAR_C(N_NODD)       ! COMPRESSED SURFACE PRESSURE (PA)
!
      REAL DTHBYDT_KM1_C(N_NODD) ! COMPRESSED INCREMENT TO MODEL
                                 ! POTENTIAL TEMPERATURE OF LAYER K-1
                                 ! (K/S)
!
      REAL DQBYDT_KM1_C(N_NODD)  ! COMPRESSED INCREMENT TO MODEL
                                 ! MIXING RATIO OF LAYER K-1 (KG/KG/S)
!
      REAL RAIN_C(N_NODD)        ! AMOUNT OF RAINFALL PASSING THROUGH
                                 ! ENVIRONMENT (KG/M**2/S)
!
      REAL SNOW_C(N_NODD)        ! AMOUNT OF SNOWFALL PASSING THROUGH
                                 ! ENVIRONMENT (KG/M**2/S)
!
      INTEGER ICCB_C(N_NODD)     ! COMPRESSED CLOUD BASE LEVEL
!
      REAL CCA_C(N_NODD)         ! COMPRESSED CONVECTIVE CLOUD AMOUNT
!
      LOGICAL BWATER_K_C(N_NODD) ! COMPRESSED MASK FOR THOSE
                                 ! POINTS AT WHICH CONDENSATE
                                 ! IS WATER IN LAYER K
!
!
      REAL                                                              &
     &    PU,PL
!
!
!
!
#include "c_r_cp.h"
!
      DO K = KCT+1,2,-1
        DO I=1,N_NODD
          THE_K_C(I) = THE(INDEX1(I),K)
          THE_KM1_C(I) = THE(INDEX1(I),K-1)
        ENDDO
        DO I=1,N_NODD
          QE_K_C(I) = QE(INDEX1(I),K)
          QE_KM1_C(I) = QE(INDEX1(I),K-1)
        ENDDO
        DO I=1,N_NODD
          DTHBYDT_KM1_C(I) = DTHBYDT(INDEX1(I),K-1)
          DQBYDT_KM1_C(I) = DQBYDT(INDEX1(I),K-1)
        ENDDO
        DO I=1,N_NODD
          EXNER_KM12_C(I) = EXNER_LAYER_BOUNDARIES(INDEX1(I),K-1)
          EXNER_KP12_C(I) = EXNER_LAYER_BOUNDARIES(INDEX1(I),K)
        ENDDO
        DO I=1,N_NODD
          EXNER_KM32_C(I) = EXNER_LAYER_BOUNDARIES(INDEX1(I),K-2)
          PRECIP_K_C(I) = PRECIP(INDEX1(I),K)
        ENDDO
        IF (K == KCT+1) THEN
          DO I=1,N_NODD
            PSTAR_C(I) = PSTAR(INDEX1(I))
!            RAIN_C(I) = RAIN(INDEX1(I))
!            SNOW_C(I) = SNOW(INDEX1(I))
            RAIN_C(I) = 0.0
            SNOW_C(I) = 0.0
            ICCB_C(I) = ICCB(INDEX1(I))
            CCA_C(I) = CCA(INDEX1(I))
          ENDDO
        ENDIF
        IF (K == KCT+1) THEN
          DO I=1,N_NODD
            EXK_C(I) = EXNER_LAYER_CENTRES(INDEX1(I),K)
          ENDDO
        ELSE
          DO I=1,N_NODD
            EXK_C(I) = EXKM1_C(I)
          ENDDO
        ENDIF
        DO I=1,N_NODD
          PKM1_C(I) = P_LAYER_CENTRES(INDEX1(I),K-1)
          DELPKM1_C(I) = P_LAYER_BOUNDARIES(INDEX1(I),K-2) -            &
     &                               P_LAYER_BOUNDARIES(INDEX1(I),K-1)
          EXKM1_C(I) =EXNER_LAYER_CENTRES(INDEX1(I),K-1)
        ENDDO
!
        DO I=1,N_NODD
          IF (BWATER(INDEX1(I),K)) THEN
            RAIN_C(I) = RAIN_C(I) + PRECIP(INDEX1(I),K)
          ELSE
            SNOW_C(I) = SNOW_C(I) + PRECIP(INDEX1(I),K)
          ENDIF
        ENDDO
!
!----------------------------------------------------------------------
! CARRY OUT CHANGE OF PHASE CALCULATION FOR PRECIPITATION FALLING
! THROUGH ENVIRONMENT
!----------------------------------------------------------------------
!
! DEPENDS ON: chg_phse
        CALL CHG_PHSE (N_NODD,K,RAIN_C,SNOW_C,DTHBYDT_KM1_C,            &
                       EXK_C,EXKM1_C,DELPKM1_C,THE_K_C,THE_KM1_C,       &
                       TIMESTEP,CCA_C)
!
!----------------------------------------------------------------------
! RESET PRECIPITATION FALLING THROUGH ENVIRONMENT IF DOWNDRAUGHT
! TERMINATES
!----------------------------------------------------------------------
!
! DEPENDS ON: pevp_bcb
        CALL PEVP_BCB (N_NODD,K-1,ICCB_C,THE_KM1_C,PKM1_C,QE_KM1_C,     &
     &                 DELPKM1_C,RAIN_C,SNOW_C,DTHBYDT_KM1_C,           &
     &                 DQBYDT_KM1_C,EXKM1_C,TIMESTEP,CCA_C)
!
        DO I=1,N_NODD
          DTHBYDT(INDEX1(I),K-1) = DTHBYDT_KM1_C(I)
          DQBYDT(INDEX1(I),K-1) = DQBYDT_KM1_C(I)
! Zero precipitation, as is (slyly) done in downd3c
          PRECIP(INDEX1(I),K) = 0.
        ENDDO
        IF (K == 2) THEN
          DO I=1,N_NODD
            RAIN(INDEX1(I)) = RAIN(INDEX1(I)) + RAIN_C(I)
            SNOW(INDEX1(I)) = SNOW(INDEX1(I)) + SNOW_C(I)
          ENDDO
        ENDIF
 
        DO I=1,N_NODD
          RAIN_3D(INDEX1(I), k-1) = RAIN_3D(INDEX1(I), k-1) + RAIN_C(I)
          SNOW_3D(INDEX1(I), k-1) = SNOW_3D(INDEX1(I), k-1) + SNOW_C(I)
        ENDDO
!
!
      ENDDO      !  MAIN LOOP OVER LEVELS
!
      RETURN
      END SUBROUTINE EVAP_BCB_NODD
#endif
