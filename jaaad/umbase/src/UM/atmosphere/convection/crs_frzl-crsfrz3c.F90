#if defined(A05_4A) || defined(A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE CRS_FRZL------------------------------------------
!LL
!LL  PURPOSE : CHANGE OF PHASE ROUTINE WHERE PRECIPITATION
!LL            CROSSES A MELTING OR FREEZING LEVEL
!LL
!LL            ADD LATENT HEATING
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
!LL
!LL  MODEL            MODIFICATION HISTORY:
!LL VERSION  DATE
!LL   4.4   17/10/97  New version optimised for T3E.
!LL                   Single PE optimisations           D.Salmond
!LL   4.5   03/03/98  Insert missing brackets. Julie Gregory.
!     6.2   03/02/05  Added section 5A R.A.Stratton
!     6.4   12/12/06  Removed def 3C. R A Stratton
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL  VERSION NO. 4  DATED 5/2/92
!LL
!LL  LOGICAL COMPONENTS COVERED:
!LL
!LL  SYSTEM TASK : P27
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE CRS_FRZL (NPNTS,RAIN,SNOW,THDD_KM1,EXKM1,              &
     &                     FLX_DD_KM1,BDDWT_KM1)
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
! MODEL CONSTANTS
!----------------------------------------------------------------------
!
#include "c_lheat.h"
#include "c_r_cp.h"
#include "c_g.h"
#include "c_0_dg_c.h"
!
!----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!----------------------------------------------------------------------
!
      INTEGER NPNTS                ! IN VECTOR LENGTH
!
      INTEGER I                    ! LOOP COUNTER
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!----------------------------------------------------------------------
!
      LOGICAL BDDWT_KM1(NPNTS)     ! IN MASK FOR POINTS WHERE
                                   !    PRECIPITATION IS LIQUID
                                   !    IN LAYER K+1
!
      REAL EXKM1(NPNTS)            ! IN EXNER RATIO FOR LAYER K-1
!
      REAL FLX_DD_KM1(NPNTS)       ! IN MASS FLUX OF LAYER K-1 (PA/S)
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT AND OUTPUT
!----------------------------------------------------------------------
!
      REAL RAIN(NPNTS)             ! INOUT
                                   ! IN  AMOUNT OF FALLING RAIN
                                   !     DESCENDING FROM LAYER
                                   !     K-1 TO K-2 (KG/M**2/S)
                                   ! OUT UPDATED AMOUNT OF FALLING
                                   !     RAIN (KG/M**2/S)
!
      REAL SNOW(NPNTS)             ! INOUT
                                   ! IN  AMOUNT OF FALLING SNOW
                                   !     DESCENDING FROM LAYER
                                   !     K-1 TO K-2 (KG/M**2/S)
                                   ! OUT UPDATED AMOUNT OF FALLING
                                   !     SNOW (KG/M**2/S)
!
      REAL THDD_KM1(NPNTS)         ! INOUT
                                   ! IN  DOWNDRAUGHT POTENTIAL
                                   !     TEMPERATURE IN LAYER K-1 (K)
                                   ! OUT UPDATED DOWNDRAUGHT POTENTIAL
                                   !     TEMPERATURE IN LAYER K-1
                                   !     DUE TO CHANGE OF PHASE (K)
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!---------------------------------------------------------------------
!
      REAL FACTOR                  ! USED IN THE CALCULATION OF
                                   ! CHANGE OF PHASE OF FALLING
                                   ! PRECIPITATION
!
      REAL PRECIP_FRE              ! FREEZING PRECIPITATION
!
      REAL PRECIP_MELT             ! MELTING PRECIPITATION
!
!L
!L----------------------------------------------------------------------
!L  ADD LATENT HEATING WHERE PRECIP CROSSES A MELTING OR FREEZING LEVEL
!L
!L  UM DOCUMENTATION PAPER 27
!L  SECTION (11), EQUATION (42)
!L----------------------------------------------------------------------
!L
       DO I=1,NPNTS
!
        IF (.NOT.BDDWT_KM1(I).AND.RAIN(I) >  0.0.AND.THDD_KM1(I)        &
     &       *EXKM1(I) <  TM) THEN
! FREEZE
          FACTOR = (EXKM1(I)*CP*FLX_DD_KM1(I))/(LF*G)
          PRECIP_FRE = (TM/EXKM1(I)-THDD_KM1(I))* FACTOR
          PRECIP_FRE = MIN(RAIN(I),PRECIP_FRE)
          THDD_KM1(I) = THDD_KM1(I)+PRECIP_FRE/FACTOR
          RAIN(I) = RAIN(I)-PRECIP_FRE
          SNOW(I) = SNOW(I)+PRECIP_FRE
!
        ELSE IF (BDDWT_KM1(I).AND.SNOW(I) >  0.0) THEN
! MELT
          FACTOR = (EXKM1(I)*CP*FLX_DD_KM1(I))/(LF*G)
          PRECIP_MELT = (THDD_KM1(I)-TM/EXKM1(I))*FACTOR
          PRECIP_MELT = MIN(SNOW(I),PRECIP_MELT)
          THDD_KM1(I) = THDD_KM1(I)-PRECIP_MELT/FACTOR
          RAIN(I) = RAIN(I)+PRECIP_MELT
          SNOW(I) = SNOW(I)-PRECIP_MELT
        END IF
      END DO
!
      RETURN
      END SUBROUTINE CRS_FRZL
!
#endif
