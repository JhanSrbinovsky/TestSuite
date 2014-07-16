
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DEEP_GRAD_STRESS------------------------------------
!LL
!LL  PURPOSE:  TO CALCULATE THE GRADIENT COMPONENT OF THE STRESS
!LL            PROFILE FOR DEEP CONVECTION. CALCULATION CAN BE
!LL            DONE EXPLICITLY OR IMPLICITLY.
!LL
!LL  CALLED FROM: CONVECT
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE DEEP_GRAD_STRESS(NP_FIELD,NPNTS,NCONV,NLEVS,NLCL,      &
     &                            NTERM,CU_TERM,CU_TEND,UE,VE,          &
     &                            VISC,PHALF,P,RHO,TIMESTEP,NTOP,       &
                                     !OUTPUT
     &                            UW,VW)
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
      REAL UE(NLEVS,NCONV+1),                                           &
                              ! ENVIRONMENT U-WIND COMPONENT (MS-1)
     &     VE(NLEVS,NCONV+1),                                           &
                              ! ENVIRONMENT V-WIND COMPONENT (MS-1)
     &     UP(NLEVS,NCONV+1),                                           &
                              ! IN CLOUD U-WIND COMPONENT (MS-1)
     &     VP(NLEVS,NCONV+1),                                           &
                              ! IN CLOUD V-WIND COMPONENT (MS-1)
     &     VISC(NLEVS,NCONV+1),                                         &
     &     PHALF(NLEVS,NCONV+1),                                        &
     &     P(NLEVS,NCONV+1),                                            &
     &     TIMESTEP             ! MODEL TIMESTEP (S)
!
      REAL RHO(NLEVS,NCONV+1)
!
! VARIABLES THAT ARE OUTPUT
!
      REAL UE_TP1(NLEVS,NCONV+1),                                       &
                                        ! IMPLICITLY UPDATED U WIND
     &     VE_TP1(NLEVS,NCONV+1),                                       &
                                        ! IMPLICITLY UPDATE V WIND
     &     UW(NLEVS,NCONV+1),                                           &
                                        ! U-COMP OF VISCOUS STRESS
     &     VW(NLEVS,NCONV+1)            ! V-COMP OF VISCOUS STRESS
!
! VARIABLES THAT ARE LOCAL
!
      INTEGER I,J,K,M,N         ! COUNTERS
!
! IMPLICIT SOLVER VARIABLES
!
      REAL A(NLEVS),B(NLEVS),C(NLEVS),                                  &
     &     U_T(NLEVS),U_TP1(NLEVS),                                     &
     &     V_T(NLEVS),V_TP1(NLEVS)
      INTEGER NLEV
!
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!
      EXTERNAL TRIDIAG

!
! USE IMPLICIT TIMESTEPPING
!

       DO I=1,NTERM
        J=CU_TERM(I)
        NLEV=0
!
! CALCULATE COMPONENTS OF TRIDIAGONAL MATRIX AND CONSTRUCT VECTOR OF
! CURRENT TIMESTEP WIND COMPONENTS
!
        DO K=NLCL(J),NTOP(J)+1
         NLEV=NLEV+1
         IF(K == NLCL(J)) THEN
           A(NLEV)=0.0
           C(NLEV)=-VISC(K+1,J)*TIMESTEP/                               &
     &              ((P(K+1,J)-P(K,J))*(PHALF(K+1,J)-PHALF(K,J)))
         ELSE IF(K <= NTOP(J)) THEN
           A(NLEV)=-VISC(K,J)*TIMESTEP/                                 &
     &              ((P(K,J)-P(K-1,J))*(PHALF(K+1,J)-PHALF(K,J)))
           C(NLEV)=-VISC(K+1,J)*TIMESTEP/                               &
     &              ((P(K+1,J)-P(K,J))*(PHALF(K+1,J)-PHALF(K,J)))
         ELSE IF(K == NTOP(J)+1) THEN
           A(NLEV)=-VISC(K,J)*TIMESTEP/                                 &
     &              ((P(K,J)-P(K-1,J))*(PHALF(K+1,J)-PHALF(K,J)))
           C(NLEV)=0.0
         ENDIF
         B(NLEV)=1.0-A(NLEV)-C(NLEV)
         U_T(NLEV)=UE(K,J)
         V_T(NLEV)=VE(K,J)
        END DO
!
! CALCULATE NEW TIMESTEP WIND COMPONENTS USING TRIDIAG
!
! DEPENDS ON: tridiag
        CALL TRIDIAG(A,B,C,U_T,U_TP1,NLEV)
! DEPENDS ON: tridiag
        CALL TRIDIAG(A,B,C,V_T,V_TP1,NLEV)
!
! STORE UPDATED WIND COMPONENTS FOR LATER
!
        NLEV=0
        DO K=NLCL(J),NTOP(J)+1
         NLEV=NLEV+1
         UE_TP1(K,J)=U_TP1(NLEV)
         VE_TP1(K,J)=V_TP1(NLEV)
        END DO
       END DO
!
! CALCULATE STRESS PROFILES
!
       DO I=1,NTERM
        M=CU_TERM(I)
        J=NLCL(M)
        UW(J,M)=0.0
        VW(J,M)=0.0
        DO K=J+1,NTOP(M)+1
         UW(K,M)=-VISC(K,M)*(UE_TP1(K,M)-UE_TP1(K-1,M))/                &
     &                      (P(K-1,M)-P(K,M))
         VW(K,M)=-VISC(K,M)*(VE_TP1(K,M)-VE_TP1(K-1,M))/                &
     &                      (P(K-1,M)-P(K,M))
        END DO
        UW(NTOP(M)+2,M)=0.0
        VW(NTOP(M)+2,M)=0.0
       END DO
      RETURN
      END SUBROUTINE DEEP_GRAD_STRESS
