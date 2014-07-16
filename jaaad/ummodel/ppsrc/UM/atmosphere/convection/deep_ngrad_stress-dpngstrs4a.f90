
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DEEP_NGRAD_STRESS------------------------------------
!LL
!LL  PURPOSE:  TO CALCULATE CLOUD BASE STRESS FOR DEEP CONVECTION
!LL            AND COMPLETE CALCULATION OF STRESS PROFILE
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
!LL                                     A.L.M. Grant
!     5.5  20/02/03   Replaced #ENDIF with #endif.      P.Dando
!     6.2  03/02/05   Added section 5A. R A Stratton
!     6.4  02/01/07   Add SCM code for zero winds. R A Stratton.
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL  DOCUMENTATION : href="http://hc0500/~hadgi/cmods.html
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE DEEP_NGRAD_STRESS(NP_FIELD,NPNTS,NCONV,NTERM,NLEVS,    &
                                  NLCL,CU_TERM,CU_COMP,CU_TEND,         &
                                   pstar,                               &
                                  UW0,VW0,ZLCL,UE,VE,VISC,MASS,P,PHALF, &
                                  RHO,TIMESTEP,NTOP,FLG_UW_DP,          &
                                  FLG_VW_DP,                            &
                                       ! INPUT/OUTPUT
                                  UW,VW,                                &
                                        ! OUTPUT
                                  UW_BASE,VW_BASE,UW_DP,VW_DP)
!
      IMPLICIT NONE
!
! VARIABLES THAT ARE INPUT
!
      INTEGER NPNTS,                                                    &
                           ! TOTAL NUMBER OF POINTS
     &        NP_FIELD,                                                 &
     &        NCONV,                                                    &
                           ! NUMBER OF CONVECTING POINTS
     &        NTERM,                                                    &
                           ! NUMBER OF TERMINATING POINTS
     &        NLEVS,                                                    &
                           ! NUMBER OF MODEL LEVELS
     &        NLCL(NCONV),                                              &
                           ! LIFTING CONDENSATION LEVEL
     &        NTOP(NCONV),                                              &
     &        CU_COMP(NPNTS),                                           &
                             ! INDEX ARRAY FOR CONVECTING POINTS
     &        CU_TERM(NPNTS),                                           &
                              ! INDEX FOR TERMINATING POINTS
     &        CU_TEND(NPNTS)  ! INDEXES FULL ARRAYS
!
      LOGICAL FLG_UW_DP,                                                &
                             ! STASH FLAGS FOR STRESS
     &        FLG_VW_DP      ! DIAGNOSTICS
!
      REAL UW0(NPNTS),                                                  &
                           ! SURFACE SHEAR STRESS X-COMPONENT (M2S-2)
     &     VW0(NPNTS),                                                  &
                           ! SURFACE SHEAR STRESS Y-COMPONENT (M2S-2)
     &     ZLCL(NPNTS),                                                 &
                           ! HEIGHT OF LCL (M)
     &     pstar(NPNTS),                                                &
                           ! surface pressure
     &     UE(NLEVS,NCONV+1),                                           &
                              ! ENVIRONMENT WIND X-COMPONENT (MS-1)
     &     VE(NLEVS,NCONV+1),                                           &
                              ! ENVIRONMENT WIND Y-COMPONENT (MS-1)
     &     VISC(NLEVS,NCONV+1),                                         &
                                ! VISCOSITY
     &     MASS(NLEVS,NCONV+1),                                         &
                                ! UPDRAUGHT MASS FLUX
     &     P(NLEVS,NCONV+1),                                            &
                                ! PRESSURES ON MODEL LEVELS (PA)
     &     PHALF(NLEVS,NCONV+1),                                        &
                                ! PROSSURES ON HALF LEVELS (PA)
     &     RHO(NLEVS,NCONV+1),                                          &
                                 ! DENSITIES (KG M-3)
     &     TIMESTEP             ! MODEL TIMESTEP (S)
!
! VARIABLES THAT ARE INPUT AND OUTPUT
!
      REAL UW(NLEVS,NCONV+1),                                           &
                               ! U-COMP OF STRESS
     &     VW(NLEVS,NCONV+1)   ! V-COMPONENT OF STRESS

!
! VARIABLES THAT ARE OUTPUT
!
      REAL UW_BASE(NCONV+1),                                            &
     &     VW_BASE(NCONV+1),                                            &
     &     UW_DP(NP_FIELD,NLEVS),                                       &
                                       ! U STRESS COMPONENT FOR STASH
     &     VW_DP(NP_FIELD,NLEVS)       ! V STRESS COMPONENT FOR STASH
!
! VARIABLES THAT ARE LOCAL
!
      INTEGER I,                                                        &
                           ! LOCAL ARRAY INDEX
     &        K,                                                        &
                           ! LEVEL INDEX
     &        J,                                                        &
                           ! INDEXES POINTS IN COMPRESSED INPUT ARRAYS
     &        M,                                                        &
     &        N            ! INDEXES POINTS IN UNCOMPRESSED INPUT ARRAYS
!
      REAL A_0(NTERM),                                                  &
                           !
     &     A_U(NTERM),                                                  &
                           ! COEFFICIENTS NEEDED FOR EVALUATING
     &     A_V(NTERM),                                                  &
                           ! IN CLOUD WIND AT CLOUD BASE
     &     OMG2_JUMP(NTERM),                                            &
                             ! CLOUD BASE JUMP IN Y COMPONENT OF VORTICI
     &     OMG1_JUMP(NTERM),                                            &
                             ! X-COMPONENT OF VORTICITY JUMP
     &     MB(NTERM),                                                   &
                             ! CLOUD-BASE MASS FLUX (MS-1)
     &     DZ,DZP1,BETA,                                                &
                                   ! TEMPORARY CALCULATIONS
     &     DU,DV,DZ1,ZETA,A,B,C,T
!
! PARAMETERS
!
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!
      REAL DELTA,TOP_PRESS,GAMMA
      PARAMETER(GAMMA=1.63,DELTA=2.0,TOP_PRESS=15000.0)
!
! CONVERT CLOUD-BASE MASS FLUX FROM PA S-1 TO MS-1
!
      DO I=1,NTERM
       J=CU_TERM(I)
       N=CU_COMP(I)
       K=NLCL(J)
       MB(I)=MASS(K,J)/G
      END DO
!
! CALCULATE JUMPS IN VORTICITY COMPONENTS ACROSS CLOUD-BASE.
! 'IMPLICIT TECHNIQUE' ASSUMS DU,DV VARY AS EXP(-T/TAU) THROUGH TIMESTEP
! NEEDED BECAUSE EXPLICIT CALULATION CAN LEAD TO INSTABILITY IN DU,DV
! UNDER SOME CIRCUMSTANCES
!
      DO I=1,NTERM
       J=CU_TERM(I)
       N=CU_COMP(I)
       DZ=-(P(NLCL(J),J)-PHALF(NLCL(J),J))/(G*RHO(NLCL(J),J))
       DU=(UE(NLCL(J),J)-UE(NLCL(J)-1,J))
       DV=(VE(NLCL(J),J)-VE(NLCL(J)-1,J))
       DZ1=-(PHALF(NLCL(J)+1,J)-PHALF(NLCL(J),J))/(G*RHO(NLCL(J)+1,J))
       ZETA=-(PHALF(NLCL(J)+1,J)-PHALF(NLCL(J),J))/25000.0
       B=(1.0/ZLCL(N)-(EXP(-ZETA)-1.0)/DZ1)
       A=ZLCL(N)*MB(I)*B/(DELTA*DZ)
       C=(B*(1.0-GAMMA/DELTA)-1/ZLCL(N))*UW0(N)





       T=-ALOG((C*(1.0-EXP(-A*TIMESTEP))/A+                             &
     &         DU*(EXP(-A*TIMESTEP)-1.0))/                              &
     &         ((C-A*DU)*TIMESTEP))/A
       OMG2_JUMP(I)=(C*(1.0-EXP(-A*T))/A+DU*EXP(-A*T))/DZ

       C=(B*(1.0-GAMMA/DELTA)-1/ZLCL(N))*VW0(N)
       T=-ALOG((C*(1.0-EXP(-A*TIMESTEP))/A+                             &
     &         DV*(EXP(-A*TIMESTEP)-1.0))/                              &
     &         ((C-A*DV)*TIMESTEP))/A
       OMG1_JUMP(I)=-(C*(1.0-EXP(-A*T))/A+DV*EXP(-A*T))/DZ
      END DO
!
! CALCULATE CLOUD BASE STRESS COMPONENTS. NOTE FACTOR OF G TO CONVERT
! BACK TO PA S-1.
!
      DO I=1,NTERM
       J=CU_TERM(I)
       N=CU_COMP(I)
       UW_BASE(J)=G*(ZLCL(N)*(-MB(I)*OMG2_JUMP(I)-                      &
     &                 GAMMA*UW0(N)/ZLCL(N))/DELTA+UW0(N))
       VW_BASE(J)=G*(ZLCL(N)*(MB(I)*OMG1_JUMP(I)-                       &
     &                 GAMMA*VW0(N)/ZLCL(N))/DELTA+VW0(N))
      END DO
!
! CALCULATE TOTAL STRESS
!
!
! CALCULATE STRESS PROFILES THE FUNCTION BETA WAS AGAIN TUNED TO TOGA-CO
! CRM SIMULATION
!
       DO I=1,NTERM
        M=CU_TERM(I)
        N=CU_TEND(I)
        J=NLCL(M)
        UW(J,M)=UW_BASE(M)
        VW(J,M)=VW_BASE(M)

! below cloud base
        DO K=1,Nlcl(M)-1
         BETA=(PHALF(K,M)-pstar(m))/(PHALF(J,M)-pstar(m))
         UW(K,M)=UW(K,M)+BETA*UW_BASE(M)
         VW(K,M)=VW(K,M)+BETA*VW_BASE(M)
        END DO

        DO K=J+1,NTOP(M)+1
         BETA=EXP(((PHALF(K,M)-PHALF(J,M))/25000.))
         UW(K,M)=UW(K,M)+BETA*UW_BASE(M)
         VW(K,M)=VW(K,M)+BETA*VW_BASE(M)
        END DO
        UW(NTOP(M)+2,M)=0.0
        VW(NTOP(M)+2,M)=0.0
!
! STASH DIAGNOSTICS (NOTE DIAGNOSTICS ARE IN MS-1 FOR DIRECT COMPARISON
! WITH SHALLOW CONVECTION STRESSES)
!
        IF(FLG_UW_DP) THEN
         DO K=J,NTOP(M)+2
          UW_DP(N,K)=UW(K,M)/G
         END DO
        ENDIF
        IF(FLG_VW_DP) THEN
         DO K=J,NTOP(M)+2
          VW_DP(N,K)=VW(K,M)/G
        END DO
       ENDIF
      END DO
      RETURN
      END SUBROUTINE DEEP_NGRAD_STRESS
