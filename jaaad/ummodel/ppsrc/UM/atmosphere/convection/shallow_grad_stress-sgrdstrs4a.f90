
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE SHALLOW_GRAD_STRESS------------------------------------
!LL
!LL  PURPOSE:  CALCULATES THE GRADIENT COMPONENT OF THE STRESS
!LL            DUE TO SHALLOW CONVECTION
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
!     6.0   05/08/03   NEC optimisation - replace function F_W by MIN
!                      & rewrite UW,VW loop.  R Barnes.
!     6.2  02/09/05   Part of version 5A. R A Stratton
!
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL  DOCUMENTATION : http://hc0500/~hadag/cmt_param.ps.gz
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE SHALLOW_GRAD_STRESS(NPNTS,N_CUMULUS,NTERM,NLEVS,       &
     &                       CU_IND,NLCL,NTOP,MB,WSC,WSTR,ZCLD,PLCL,    &
     &                       PTOP,P,PHALF,RHO,UE,VE,TIMESTEP,           &
                                     ! OUTPUT ARGUMENTS
     &                       UW,VW)
!
      IMPLICIT NONE
!
! ARGUMENTS THAT ARE INPUT
!
      INTEGER NPNTS,                                                    &
                              ! TOTAL NUMBER OF POINTS IN SEGMENT
     &        N_CUMULUS,                                                &
                              ! NUMBER OF POINTS DIAGNOSED WITH CUMULUS
     &        NTERM,                                                    &
                              ! NUMBER OF TERMINATING POINTS
     &        NLEVS,                                                    &
                              ! NUMBER OF MODEL LEVELS
     &        CU_IND(NTERM),                                            &
                              ! INDEX OF TERMINATING POINTS
     &        NLCL(N_CUMULUS),                                          &
                               ! HALF LEVEL OF LCL
     &        NTOP(N_CUMULUS) ! TOP OF SHALLOW CUMULUS LAYER
!
      REAL MB(N_CUMULUS),                                               &
                              ! CLOUD BASE MASS FLUX FOR SHALLOW CU (MS-
     &     WSC(N_CUMULUS),                                              &
                              ! CLOUD-LAYER VELOCITY SCALE (MS-1)
     &     WSTR(N_CUMULUS),                                             &
                              ! MIXED LAYER VELCOITY SCALE (MS-1)
     &     ZCLD(N_CUMULUS),                                             &
                              ! CLOUD LAYER DEPTH (M)
     &     PLCL(N_CUMULUS),                                             &
                              ! PRESSURE AT LIFTING CONDENSATION LEVEL (
     &     PTOP(N_CUMULUS),                                             &
                              ! PRESSURE AT TOP OF CLOUD LAYER (PA)
     &     P(NLEVS,N_CUMULUS+1),                                        &
                                 ! MODEL LEVEL PRESSURES (PA)
     &     PHALF(NLEVS,N_CUMULUS+1),                                    &
                                     ! MODEL HALF LEVEL PRESSURES (PA
     &     RHO(NLEVS,N_CUMULUS+1),                                      &
                                   ! DENSITY ON MODEL LEVELS (KG M-3)
     &     UE(NLEVS,N_CUMULUS+1),                                       &
                                    ! U-COMPONENT OF MEAN WIND (MS-1)
     &     VE(NLEVS,N_CUMULUS+1),                                       &
                                     ! V-COMPONENT OF MEAN WIND (MS-1)
     &     TIMESTEP               ! MODEL TIMESTEP (S)
!
! VARIABLES THAT ARE OUTPUT
!
      REAL UW(NLEVS,N_CUMULUS+1),                                       &
                                    ! U-COMPONENT OF STRESS
     &     VW(NLEVS,N_CUMULUS+1)    ! V-COMPONENT OF STRESS
!
! VARIABLES THAT ARE LOCAL
!
      INTEGER I,J,K,M,NLEV            ! COUNTERS
!
      REAL W(NLEVS,N_CUMULUS),                                          &
                                 ! NON-DIMENSIONAL PLUME VERTICAL VELOCI
     &     MASS(NLEVS,N_CUMULUS),                                       &
                                  ! MASS FLUX PROFILE (MS-1)
     &     VISC(NLEVS,NTERM),                                           &
                                  ! VISCOSITY PROFILE (M2S-1)
     &     A(NLEVS),B(NLEVS),C(NLEVS),                                  &
                                       ! TRIDIAGONAL MATRIX ELEMENTS
     &     U_T(NLEVS),V_T(NLEVS),                                       &
                                  ! CURRENT VELOCITY VECTORS
     &     U_TP1(NLEVS),V_TP1(NLEVS),                                   &
                                      ! AFTER TIMESTEP VELOCITY VECTORS
     &     UE_TP1(NLEVS,N_CUMULUS+1),                                   &
                                      ! AFTER TIMESTEP VELOCITY VECTORS
     &     VE_TP1(NLEVS,N_CUMULUS+1),                                   &
                                      ! FOR SUBSEQUENT USE
     &     W02,P_DEPTH,ZETA,ENTR_SC,EXP_K,EXP_KP1,DZ,DZ12
!
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!
      REAL A_STRESS
      PARAMETER (A_STRESS=0.3)
      REAL A_W02
      PARAMETER (A_W02=10.24)
!
      REAL F_W
      EXTERNAL F_W
      EXTERNAL TRIDIAG
!
! CALCULATE VERTICAL VELOCITY PROFILE IN UPDRAUGHTS
!
      DO I=1,NTERM
       J=CU_IND(I)
       W02=A_W02*(MB(J)*WSTR(J)**2.0)**0.6667
       P_DEPTH=PTOP(J)-PLCL(J)
       DO K=NLCL(J),NTOP(J)
        ZETA=(PHALF(K,J)-PLCL(J))/P_DEPTH
! DEPENDS ON: f_w
        W(K,I)=SQRT(W02+WSC(J)**2*F_W(ZETA))/WSC(J)
       END DO
      END DO
!
! CALCULATE MASS FLUX PROFILE
! USES FRACTIONAL DETRAINMENT=1.3 FRACTIONAL ENTRAINMENT
!
      DO I=1,NTERM
       J=CU_IND(I)
       ENTR_SC=0.04*WSC(J)/(MB(J)*ZCLD(J))
       P_DEPTH=PTOP(J)-PLCL(J)
       MASS(NLCL(J),I)=MB(J)
       EXP_K=EXP(-(PHALF(NLCL(J),J)-PLCL(J))/P_DEPTH)
       DO K=NLCL(J),NTOP(J)-1
        EXP_KP1=EXP(-(PHALF(K+1,J)-PLCL(J))/P_DEPTH)
        ZETA=(1.0-1.3)*ENTR_SC*P_DEPTH*(EXP_KP1-EXP_K)/(G*RHO(K+1,J))
        MASS(K+1,I)=MASS(K,I)*EXP(ZETA)
        EXP_K=EXP_KP1
       END DO
      END DO
!
! CALCULATE THE EDDY VISCOSITY PROFILE
!
      DO I=1,NTERM
       J=CU_IND(I)
       DO K=NLCL(J)+1,NTOP(J)+1
        IF(K <  NTOP(J)) THEN
         VISC(K,I)=A_STRESS*MASS(K,I)*W(K,I)*ZCLD(J)
        ELSE IF(K == NTOP(J)) THEN
         VISC(K,I)=0.162*MB(J)*ZCLD(J)
        ELSE IF(K == NTOP(J)+1) THEN
         DZ=-(P(K,J)-P(K-1,J))/(G*(RHO(K,J)+RHO(K-1,J))/2.0)
         VISC(K,I)=0.09*MB(J)*DZ
        ENDIF
       END DO
      END DO
!
! CALCULATE GRADIENT COMPONENT OF STRESS
!
!
! USE IMPLICIT TIMESTEPPING
!

       DO I=1,NTERM
        J=CU_IND(I)
        NLEV=0
!
! CALCULATE COMPONENTS OF TRIDIAGONAL MATRIX AND CONSTRUCT VECTOR OF
! CURRENT TIMESTEP WIND COMPONENTS
!
        DO K=NLCL(J),NTOP(J)+1
         NLEV=NLEV+1
         IF(K == NLCL(J)) THEN
           DZ=-(PHALF(K+1,J)-PHALF(K,J))/(G*RHO(K,J))
           DZ12=-(P(K+1,J)-P(K,J))/(G*(RHO(K+1,J)+RHO(K,J))/2.0)
           A(NLEV)=0.0
           C(NLEV)=-VISC(K+1,I)*TIMESTEP/(DZ*DZ12)
         ELSE IF(K <= NTOP(J)) THEN
           DZ=-(PHALF(K+1,J)-PHALF(K,J))/(G*RHO(K,J))
           DZ12=-(P(K,J)-P(K-1,J))/(G*(RHO(K,J)+RHO(K-1,J))/2.0)
           A(NLEV)=-VISC(K,I)*TIMESTEP/(DZ*DZ12)
           DZ12=-(P(K+1,J)-P(K,J))/(G*(RHO(K+1,J)+RHO(K,J))/2.0)
           C(NLEV)=-VISC(K+1,I)*TIMESTEP/(DZ*DZ12)
         ELSE IF(K == NTOP(J)+1) THEN
           DZ=-(PHALF(K+1,J)-PHALF(K,J))/(G*RHO(K,J))
           DZ12=-(P(K,J)-P(K-1,J))/(G*(RHO(K,J)+RHO(K-1,J))/2.0)
           A(NLEV)=-VISC(K,I)*TIMESTEP/(DZ*DZ12)
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
        J=CU_IND(I)
        M=NLCL(J)
        UW(M,J)=0.0
        VW(M,J)=0.0
!CDIR NODEP
        DO K=M+1,NTOP(J)+1
         DZ=-(P(K,J)-P(K-1,J))/(G*(RHO(K,J)+RHO(K-1,J))/2.0)
         UW(K,J)=-VISC(K,I)*(UE_TP1(K,J)-UE_TP1(K-1,J))/DZ
         VW(K,J)=-VISC(K,I)*(VE_TP1(K,J)-VE_TP1(K-1,J))/DZ
        END DO
        UW(NTOP(J)+2,J)=0.0
        VW(NTOP(J)+2,J)=0.0
       END DO
      RETURN
      END SUBROUTINE SHALLOW_GRAD_STRESS
