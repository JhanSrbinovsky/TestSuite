#if defined(A05_4A) || defined(A05_5A)
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE SHALLOW_BASE_STRESS------------------------------------
!LL
!LL  PURPOSE:  TO CALCULATE CLOUD BASE STRESS FOR SHALLOW CU
!LL            (ALSO COMPLETES CALCULATION OF STRESS PROFILE)
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
!
!     6.2  09/09/04   Correction to density weighting. R A Stratton
!     6.1  17/06/04   Correct error in calculation of vw. R A Stratton
!     6.2  02/09/05   part of 5A version. R A Stratton
!     6.4  02/01/07   Extra ifs for SCM case of zero winds. R A Stratton
!
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL  DOCUMENTATION : http://hc0500/~hadag/cmt_param.ps.gz
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE SHALLOW_BASE_STRESS(NP_Field,NPNTS,N_CUMULUS,NLEVS,    &
     &                               NTERM,                             &
     &                               CU_IND,CU_FULL,NLCL,NTOP,MB,WSC,   &
     &                               ZLCL,ZCLD,UW0,VW0,PLCL,PTOP,UE,    &
     &                               VE,PHALF,P,RHO,TIMESTEP,           &
     &                               FLG_UW_SHALL,FLG_VW_SHALL,         &
                                           ! IN/OUT ARGUMENTS
     &                               UW,VW,                             &
                                           ! OUTPUT ARGUMENTS
     &                               UW_SHALL,VW_SHALL)
!
      IMPLICIT NONE
!
! VARIABLES THAT ARE INPUT
!
      INTEGER NPNTS,                                                    &
                           ! TOTAL NUMBER OF POINTS
     &        np_field,                                                 &
     &        N_CUMULUS,                                                &
                           ! NUMBER OF CUMULUS POINTS
     &        NLEVS,                                                    &
                           ! NUMBER OF MODEL LEVELS
     &        NTERM,                                                    &
                           ! NUMBER OF POINTS TERMINATING
     &        CU_IND(NTERM),                                            &
                            ! INDEX TO TERMINATING POINTS IN COMPRESSED
     &        CU_FULL(N_CUMULUS),                                       &
                                  ! INDEX TO CUMULUS POINTS IN FULL ARRA
     &        NLCL(N_CUMULUS),                                          &
                               ! LEVELS OF LCL
     &        NTOP(N_CUMULUS)  ! LEVELS OF TOP OF CLOUD LAYER
!
      LOGICAL FLG_UW_SHALL,                                             &
                               ! STASH FLAGS FOR SHALLOW
     &        FLG_VW_SHALL     ! CONVECTION STRESS DIAGNOSTIC
!
      REAL MB(N_CUMULUS),                                               &
                             ! CLOUD BASE MASS FLUX (MS-1)
     &     WSC(N_CUMULUS),                                              &
                             ! CLOUD-LAYER VELOCITY SCALE (MS-1)
     &     ZLCL(NPNTS),                                                 &
                             ! HEIGHT OF LCL (M)
     &     ZCLD(N_CUMULUS),                                             &
                             ! DEPTH OF CLOUD-LAYER TOP (M)
     &     PLCL(N_CUMULUS),                                             &
                             ! PRESSURE OF LCL(PA)
     &     PTOP(N_CUMULUS),                                             &
                             ! PRESSURE AT TOP OF CLOUD-LAYER (PA)
     &     UW0(NPNTS),                                                  &
                             ! U-COMPONENT OF SURFACE STRESS (M2S-2)
     &     VW0(NPNTS),                                                  &
                             ! V-COMPONENT OF SURFACE STRESS (M2S-2)
     &     PHALF(NLEVS,N_CUMULUS+1),                                    &
                                     ! PRESSURE ON MODEL HALF LEVELS (PA
     &     P(NLEVS,N_CUMULUS+1),                                        &
                                     ! PRESSURE ON MODEL LEVELS (PA)
     &     RHO(NLEVS,N_CUMULUS+1),                                      &
                                     ! DENSITY ON MODEL LEVELS (KG M-3)
     &     UE(NLEVS,N_CUMULUS+1),                                       &
                                     ! U-COMPONENT OF MEAN WIND (MS-1)
     &     VE(NLEVS,N_CUMULUS+1),                                       &
                                      ! V-COMPONENT OF MEAN WIND (MS-1)
     &     TIMESTEP                 ! MODEL TIMESTEP (S)
!
! VARIABLES THAT ARE INPUT/OUTPUT
!
      REAL UW(NLEVS,N_CUMULUS+1),                                       &
                                    ! U-COMPONENT OF STRESS PROFILE (M2S
     &     VW(NLEVS,N_CUMULUS+1)   ! V-COMPONENT OF STRESS PROFILE (M2S-
!
! VARIABLES THAT ARE OUTPUT
!
      REAL UW_SHALL(NP_FIELD,NLEVS),                                    &
                                      ! STASH DIAGNOSTIC FOR U-COMP STRE
     &     VW_SHALL(NP_FIELD,NLEVS)   ! STASH DIAGNOSTIC FOR V-COMP STRE
!
!
! VARIABLES THAT ARE LOCAL
!
      INTEGER I,J,K,n            ! COUNTERS
!
      REAL DELTA_Z(NTERM),                                              &
                                     ! LAYER THICKNESS ABOVE LCL (M)
     &     COEFF_1(NTERM),                                              &
                                     ! COEFFICIENT
     &     OMG2_JUMP(NTERM),                                            &
                                     ! JUMP IN Y COMPONENT OF VORTICITY
     &     OMG1_JUMP(NTERM),                                            &
                                     ! JUMP IN X-COMPONENT OF VORTICITY
     &     U_JUMP(NTERM),                                               &
     &     V_JUMP(NTERM),                                               &
     &     ZETA,DZ,P_DEPTH,A,B,C,T,DU,DV,DZ1,rho_h
!
#include "c_g.h"
!
      REAL BETA,DELTA,GAMMA
      PARAMETER (BETA=0.04,DELTA=2.3,GAMMA=1.63)
!
! CALCULATE JUMPS IN VORTICITY ACROSS CLOUD BASE (THIS IS DONE BY ASSUMI
! THAT DURING THE TIMESTEP DU AND DV VARY AS EXP(-T/TAU). NEEDS TO BE DO
! TO AVOID INSTABILITY AROUND CLOUD BASE).
!
      DO I=1,NTERM
       J=CU_IND(I)
       N=CU_FULL(I)
       DZ=-(P(NLCL(J),J)-PLCL(J))/(G*RHO(NLCL(J),J))
       DU=(UE(NLCL(J),J)-UE(NLCL(J)-1,J))
       DV=(VE(NLCL(J),J)-VE(NLCL(J)-1,J))
       DZ1=-(PHALF(NLCL(J)+1,J)-PHALF(NLCL(J),J))/(G*RHO(NLCL(J)+1,J))
       P_DEPTH=(PTOP(J)-PLCL(J))
       ZETA=BETA*WSC(J)*(PHALF(NLCL(J)+1,J)-PLCL(J))/(MB(J)*P_DEPTH)
       B=(1.0/ZLCL(N)-(EXP(-ZETA)-1.0)/DZ1)
       A=ZLCL(N)*MB(J)*B/(DELTA*DZ)
       C=(B*(1.0-GAMMA/DELTA)-1.0/ZLCL(N))*UW0(N)
#if defined(SCMA)
       IF (C == 0.0.AND.(DV == 0.0.OR.DU == 0.0)) THEN
         OMG2_JUMP=0.0
       ELSE      
         T=-LOG((C*(1.0-EXP(-A*TIMESTEP))/A+                            &
     &         DU*(EXP(-A*TIMESTEP)-1.0))/                              &
     &         ((C-A*DU)*TIMESTEP))/A
         OMG2_JUMP(I)=(C*(1.0-EXP(-A*T))/A+DU*EXP(-A*T))/DZ
       END IF 
#else
       T=-ALOG((C*(1.0-EXP(-A*TIMESTEP))/A+                             &
     &         DU*(EXP(-A*TIMESTEP)-1.0))/                              &
     &         ((C-A*DU)*TIMESTEP))/A
       OMG2_JUMP(I)=(C*(1.0-EXP(-A*T))/A+DU*EXP(-A*T))/DZ
#endif
       C=(B*(1.0-GAMMA/DELTA)-1.0/ZLCL(N))*VW0(N)
#if defined(SCMA)
       IF (C == 0.0.AND.(DV == 0.0.OR.DU == 0.0)) THEN
         OMG1_JUMP=0.0
       ELSE      
         T=-LOG((C*(1.0-EXP(-A*TIMESTEP))/A+                            &
     &         DV*(EXP(-A*TIMESTEP)-1.0))/                              &
     &         ((C-A*DV)*TIMESTEP))/A
         OMG1_JUMP(I)=-(C*(1.0-EXP(-A*T))/A+DV*EXP(-A*T))/DZ
       END IF 
#else
       T=-ALOG((C*(1.0-EXP(-A*TIMESTEP))/A+                             &
     &         DV*(EXP(-A*TIMESTEP)-1.0))/                              &
     &         ((C-A*DV)*TIMESTEP))/A
       OMG1_JUMP(I)=-(C*(1.0-EXP(-A*T))/A+DV*EXP(-A*T))/DZ
#endif
      END DO
!
! CALCULATE THE CLOUD-BASE STRESS COMPONENTS
!
      DO I=1,NTERM
       J=CU_IND(I)
       N=CU_FULL(I)
       UW(NLCL(J),J)=(ZLCL(N)*(-MB(J)*OMG2_JUMP(I)-                     &
     &                 GAMMA*UW0(N)/ZLCL(N))/DELTA+UW0(N))
       VW(NLCL(J),J)=(ZLCL(N)*(MB(J)*OMG1_JUMP(I)-                      &
     &                 GAMMA*VW0(N)/ZLCL(N))/DELTA+VW0(N))
      END DO
!
! CALCULATE NON-GRADENT STRESS PROFILE
!
      DO I=1,NTERM
       J=CU_IND(I)
       N=CU_FULL(I)
       P_DEPTH=(PTOP(J)-PLCL(J))
       DO K=NLCL(J)+1,NTOP(J)+1
        ZETA=BETA*WSC(J)*(PHALF(K,J)-PLCL(J))/(MB(J)*P_DEPTH)
        RHO_H=RHO(K-1,J)+(RHO(K,J)-RHO(K-1,J))/(P(K,J)-P(K-1,J))*       &
     &                   (PHALF(K,J)-P(K-1,J))
        IF(K <  NTOP(J)) THEN
         UW(K,J)=RHO_H*(UW(K,J)+UW(NLCL(J),J)*EXP(-ZETA))
         VW(K,J)=RHO_H*(VW(K,J)+VW(NLCL(J),J)*EXP(-ZETA))
        ELSE IF(K == NTOP(J)) THEN
         UW(K,J)=RHO_H*(UW(K,J)+UW(NLCL(J),J)*EXP(-BETA*WSC(J)/MB(J)))
         VW(K,J)=RHO_H*(VW(K,J)+VW(NLCL(J),J)*EXP(-BETA*WSC(J)/MB(J)))
        ELSE IF(K == NTOP(J)+1) THEN
         UW(K,J)=RHO_H*(UW(K,J)+UW(NLCL(J),J)*                          &
     &                          EXP(-1.75*BETA*WSC(J)/MB(J)))
         VW(K,J)=RHO_H*(VW(K,J)+VW(NLCL(J),J)*                          &
     &                          EXP(-1.75*BETA*WSC(J)/MB(J)))
        ENDIF
       END DO

! Weight cloud base stress by rho (omitted from above level loop)
! Needs to be done after level loop.

       k=nlcl(j)
       RHO_H=RHO(K-1,J)+(RHO(K,J)-RHO(K-1,J))/(P(K,J)-P(K-1,J))*        &
     &                   (PHALF(K,J)-P(K-1,J))
       UW(NLCL(J),J) = RHO_H*UW(NLCL(J),J)
       VW(NLCL(J),J) = RHO_H*VW(NLCL(J),J)


       IF(FLG_UW_SHALL) THEN
        DO K=NLCL(J),NTOP(J)+1
         UW_SHALL(N,K)=UW(K,J)
        END DO
       ENDIF
       IF(FLG_VW_SHALL) THEN
        DO K=NLCL(J),NTOP(J)+1
         VW_SHALL(N,K)=VW(K,J)
        END DO
       ENDIF
      END DO
      RETURN
      END SUBROUTINE SHALLOW_BASE_STRESS
#endif
