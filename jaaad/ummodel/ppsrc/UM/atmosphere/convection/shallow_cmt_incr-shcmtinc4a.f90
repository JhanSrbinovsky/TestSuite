
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE SHALLOW_CMT_INCR------------------------------------
!LL
!LL  PURPOSE:  CALCULATES INCREMENTS TO U AND V DUE TO SHALLOW
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
!LL   5.4  7/8/2002   New deck for revised convection scheme 4A
!LL
!LL                                     A.L.M. Grant
!     5.5  20/02/03   Replaced #ENDIF with #endif.      P.Dando
!     6.2  02/09/05   Part of 5A version R.A.Stratton

!
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL  DOCUMENTATION : http://hc0500/~hadag/cmt_param.ps.gz
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE SHALLOW_CMT_INCR(NP_FIELD,NPNTS,N_CUMULUS,NLEVS,NTERM, &
     &                            CU_IND,CU_FULL,NLCL,NTOP,UW,VW,PHALF, &
     &                            RHO,ZLCL,                             &
                                        ! OUTPUT ARGUMENTS
     &                            DUBYDT,DVBYDT)
!
      IMPLICIT NONE
!
! VARIABLES THAT ARE INPUT
!
      INTEGER NP_FIELD,                                                 &
                            ! TOTAL NUMBER OF POINTS ON PE
     &        NPNTS,                                                    &
                            ! TOTAL NUMBER OF POINTS IN SEGMENT
     &        N_CUMULUS,                                                &
                            ! NUMBER OF POINTS DIAGNOSED CUMULUS
     &        NLEVS,                                                    &
                            ! NUMBER OF MODEL LEVELS
     &        NTERM,                                                    &
                            ! NUMBER OF TERMINATING POINTS
     &        CU_IND(NTERM),                                            &
                             ! INDEX ARRAY TO COMPRESSED ARRAYS
     &        CU_FULL(N_CUMULUS),                                       &
                                   ! INDEX ARRAY TO FULL ARRAYS
     &        NLCL(N_CUMULUS),                                          &
                                   ! MODEL HALF LEVEL OF LCL
     &        NTOP(N_CUMULUS)      ! MODEL HALF LEVEL OF TOP OF LAYER
!
      REAL UW(NLEVS,N_CUMULUS+1),                                       &
                                     ! U-COMPONENT OF STRESS PROFILE (M2
     &     VW(NLEVS,N_CUMULUS+1),                                       &
                                     ! V-COMPONENT OF STRESS PROFILE (M2
     &     PHALF(NLEVS,N_CUMULUS+1),                                    &
                                     ! PRESSURE OFMODEL HALF LEVELS (PA)
     &     RHO(NLEVS,N_CUMULUS+1),                                      &
                                     ! DENSITY ON MODEL LEVELS (KG M-3)
     &     ZLCL(NPNTS)               ! HEIGHT OF LCL (M)
!
! VARIABLES THAT ARE OUTPUT
!
      REAL DUBYDT(NP_FIELD,NLEVS),                                      &
                                        ! TENDENCY IN U (MS-2)
     &     DVBYDT(NP_FIELD,NLEVS)       ! TENDENCY IN V (MS-2)
!
! VARIABLES THAT ARE LOCAL
!
      INTEGER I,J,K,M            ! COUNTERS AND POINTERS
!
      REAL DUDT_BL,                                                     &
                    ! TENDENCY IN SUBCLOUD LAYER
     &     DVDT_BL,                                                     &
                    ! TENDENCY IN SUBCLOUD LAYER
     &     RHODZ    ! DENSITY TIMES HEIGHT DIFFERENCE BETWEEN HALF LEVEL
!
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!
      DO I=1,NTERM
       J=CU_IND(I)
       M=CU_FULL(I)
       DUDT_BL=-UW(NLCL(J),J)/ZLCL(M)
       DVDT_BL=-VW(NLCL(J),J)/ZLCL(M)
!
! MIXED LAYER INCREMENTS (ASSUMED CONSTANT IN MIXED LAYER)
!
       DO K=1,NLCL(J)-1
        DUBYDT(M,K)=DUDT_BL/RHO(K,J)
        DVBYDT(M,K)=DVDT_BL/RHO(K,J)
       END DO
!
! CLOUD LAYER INCREMENTS
!
       DO K=NLCL(J),NTOP(J)
        RHODZ=-(PHALF(K+1,J)-PHALF(K,J))/G
        DUBYDT(M,K)=-(UW(K+1,J)-UW(K,J))/RHODZ
        DVBYDT(M,K)=-(VW(K+1,J)-VW(K,J))/RHODZ
       END DO
      END DO
      RETURN
      END SUBROUTINE SHALLOW_CMT_INCR
