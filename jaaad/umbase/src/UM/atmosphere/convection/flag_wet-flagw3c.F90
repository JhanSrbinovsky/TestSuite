#if defined(A05_4A) || defined(A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE FLAG_WET-----------------------------------------------
!LL
!LL  PURPOSE : CALCULATES A MASK FOR WHEN CONDENSATION IS LIQUID
!LL
!LL            IF 0.5 * (TK + TK+1) > TICE THEN ANY CONDENSATION
!LL                                        IN LAYER K+1 IS LIQUID
!LL
!LL            IF 0.5 * (TK + TK+1) < TICE THEN ANY CONDENSATION
!LL                                        IN LAYER K+1 IS ICE
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  PROGRAMMING STANDARDS :
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
!LL                  SECTION (2B)
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
! cjj      SUBROUTINE FLAG_WET (BWATER,TH,EXNER,PSTAR,AKH,BKH,
      SUBROUTINE FLAG_WET (BWATER,TH,EXNER_LAYER_CENTRES,               &
     &                     NP_FIELD,NPNTS,NLEV)

      Use cv_run_mod, Only:                                             &
          tice

!
!-----------------------------------------------------------------------
!   RETURNS 'BWATER' - A BIT VECTOR OF POINTS WHERE CONDENSATE IS WATER
!   RATHER THAN ICE.
!----------------------------------------------- AUTHOR: M FISHER 1987 -
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
! MODEL CONSTANTS
!----------------------------------------------------------------------
!
#include "c_r_cp.h"
!
!----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!----------------------------------------------------------------------
!
      INTEGER NP_FIELD           ! IN FULL VECTOR LENGTH
!
      INTEGER NPNTS              ! IN VECTOR LENGTH
!
      INTEGER NLEV               ! IN NUMBER OF MODEL LAYERS
!
      INTEGER I,K                ! LOOP COUNTERS
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!----------------------------------------------------------------------
!
      REAL TH(NP_FIELD,NLEV)        ! IN POTENTIAL TEMPERATURE (K)
!
! cjj      REAL EXNER(NP_FIELD,NLEV+1)   ! IN EXNER RATIO AT LAYER
      REAL EXNER_LAYER_CENTRES(NP_FIELD,0:NLEV)
                                    ! IN EXNER RATIO AT LAYER
                                    ! CENTRES (STARTING WITH THE
                                    ! SURFACE)
!
! cjj      REAL PSTAR(NPNTS)             ! IN Surface pressure
!
! cjj      REAL AKH(NLEV+1)              ! IN Hybrid coordinate A at
                                    !    layer boundary
! cjj      REAL BKH(NLEV+1)              ! IN Hybrid coordinate B at
                                    !    layer boundary
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE OUTPUT
!----------------------------------------------------------------------
!
      LOGICAL BWATER(NPNTS,2:NLEV)  ! OUT MASK FOR THOSE POINTS AT
                                    !     WHICH CONDENSATE IS LIQUID
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!----------------------------------------------------------------------
!
! cjj      REAL EXK                      ! EXNER RATIO FOR LEVEL K
! cjj      REAL EXKP1                    ! EXNER RATIO FOR LEVEL K+1
!

! cjj      REAL
! cjj     &    PU,PL,PU2
! cjj *CALL P_EXNERC

!*---------------------------------------------------------------------
!L
!L---------------------------------------------------------------------
!L  NO SIGNIFICANT STRUCTURE
!L---------------------------------------------------------------------
!L
      DO K=1,NLEV-1
       DO I=1,NPNTS
!
! cjj      PU2=PSTAR(I)*BKH(K+2) + AKH(K+2)
! cjj      PU=PSTAR(I)*BKH(K+1) + AKH(K+1)
! cjj      PL=PSTAR(I)*BKH(K) + AKH(K)
! cjj      EXK = P_EXNER_C(EXNER(I,K+1),EXNER(I,K),PU,PL,KAPPA)
! cjj      EXKP1 = P_EXNER_C(EXNER(I,K+2),EXNER(I,K+1),PU2,PU,KAPPA)
!
! cjj      BWATER(I,K+1) = 0.5*(TH(I,K)*EXK + TH(I,K+1)*EXKP1)  >   TICE
        BWATER(I,K+1) = 0.5*(TH(I,K)*EXNER_LAYER_CENTRES(I,K) +         &
     &  TH(I,K+1)*EXNER_LAYER_CENTRES(I,K+1))  >   TICE
       ENDDO ! NPNTS
      ENDDO ! NLEV-1

!
      RETURN
      END SUBROUTINE FLAG_WET
#endif
