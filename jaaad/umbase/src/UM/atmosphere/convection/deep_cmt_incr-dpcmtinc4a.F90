#if defined(A05_4A) || defined(A05_5A)
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DEEP_CMT_INCR------------------------------------
!LL
!LL  PURPOSE:  CALCULATES INCREMENTS TO U AND V DUE TO DEEP
!LL            CONVECTION
!LL
!LL  CALLED FROM: CONVECT
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL
!LL  MODEL            MODIFICATION HISTORY:
!LL VERSION  DATE
!LL
!LL   5.4  07/8/2002   New deck for revised convection scheme 4A
!LL                                     A.L.M. Grant
!     5.5  20/02/03    Replaced #ENDIF with #endif.      P.Dando
!     6.2  03/02/05   Added section 5A. R A Stratton
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL  DOCUMENTATION : http://hc0500/~hadag/cmt_param.ps.gz
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE DEEP_CMT_INCR(NP_FIELD,NPNTS,NCONV,NLEVS,NLCL,NTERM,   &
     &                    CU_TERM,CU_TEND,ZLCL,PHALF,P,RHO,             &
     &                    UW_BASE,VW_BASE,UW,VW,NTOP,                   &
                                     !OUTPUT
     &                    DUBYDT,DVBYDT)
!
      IMPLICIT NONE
!
! VARIABLES THAT ARE INPUT
!
      INTEGER NPNTS,                                                    &
                            ! TOTAL NUMBER OF POINTS IN SEGMENT
     &        NP_FIELD,                                                 &
                            ! TOTAL NUMBER OF MODEL POINTS FOR PE
     &        NCONV,                                                    &
                            ! NUMBER OF CUMULUS POINTS
     &        NLEVS,                                                    &
                            ! NUMBER OF MODEL LEVELS
     &        NLCL(NCONV),                                              &
                            ! LIFTING CONDENSATION LEVEL
     &        NTOP(NCONV),                                              &
     &        NTERM,                                                    &
                            ! NO OF TERMINATING POINTS,
     &        CU_TERM(NTERM),                                           &
                              ! INDEX OF TERMINATING POINTS
     &        CU_TEND(NTERM)  ! INDEX OF POINTS IN OUTPUT ARRAY
!
      REAL UW(NLEVS,NCONV+1),                                           &
                              ! U-COMP OF STRESS (PA M S-2)
     &     VW(NLEVS,NCONV+1),                                           &
                              ! V-COMP OF STRESS (PA M S-2)
     &     PHALF(NLEVS,NCONV+1),                                        &
                                 ! PRESSURE ON HALF LEVELS
                                 ! RELATIVE TO UV GRID (PA)
     &     P(NLEVS,NCONV+1),                                            &
                             ! PRESSURE ON LEVELS (UV-GRID) (PA)
     &     ZLCL(NPNTS)       ! HEIGHT OF LIFTING CONDENSATION LEVEL
                             ! (CLOUD BASE) (M)
!
      REAL RHO(NLEVS,NCONV+1),                                          &
                               ! AIR DENSITY (KG/M3)
     &     UW_BASE(NCONV+1),                                            &
                             ! U-COMP OF STRESS  AT CLOUD-BASE
     &     VW_BASE(NCONV+1)  ! V-COMP OF STRESS AT CLOUD-BASE
!
! VARIABLES THAT ARE OUTPUT
!
      REAL DUBYDT(NP_FIELD,NLEVS),                                      &
                                   ! U INCREMENT (MS-2)
     &     DVBYDT(NP_FIELD,NLEVS) ! V INCREMENT (MS-2)
!
! VARIABLES THAT ARE LOCAL
!
      INTEGER I,J,K,M,N         ! COUNTERS
!
      REAL DUDT_BL,                                                     &
                      ! U CMT TENDENCY IN SUBCLOUD-LAYER (MS-2)
     &     DVDT_BL,                                                     &
                      ! V CMT TENDENCY IN SUBCLOUD-LAYER (MS-2)
     &     DP         ! PRESSURE DIFFERENCE BETWEEN ADJACENT HALF LEVELS
                      ! (PA)
!
#include "c_g.h"
!
!
! CALCULATE U AND V WIND INCREMENTS BY DIFFERENTIATING STRESS PROFILE
!
      DO I=1,NTERM
        M=CU_TERM(I)
        N=CU_TEND(I)
        J=NLCL(M)
!
! CMT TENDENCIES IN THE SUBCLOUD LAYER ARE CONSTANT WITH HEIGHT
!
        DUDT_BL=-UW(J,M)/(G*ZLCL(N))
        DVDT_BL=-VW(J,M)/(G*ZLCL(N))
       DO K=1,J-1
        DUBYDT(N,K)=DUDT_BL
        DVBYDT(N,K)=DVDT_BL
       END DO
       DO K=J,NTOP(M)+1
        DP=-(PHALF(K+1,M)-PHALF(K,M))
        DUBYDT(N,K)=-(UW(K+1,M)-UW(K,M))/DP
        DVBYDT(N,K)=-(VW(K+1,M)-VW(K,M))/DP
       END DO
      END DO
      RETURN
      END SUBROUTINE DEEP_CMT_INCR
#endif
