#if defined(A05_4A) || defined (A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!LL  SUBROUTINE CMT_MASS------------------------------------
!LL
!LL  PURPOSE:  CALCULATES THE MASS FLUX PROFILE FOR DEEP CONVECTION
!LL            TO BE USED IN CMT CALCULATIONS. USES THE CLOUD-BASE
!LL            MASS FLUX FROM THE PLUME SCHEME, BUT PROFILE IS NOT
!LL            THE SAME AS USED FOR THE THERMODYNAMIC PART OF THE
!LL            CONVECTION SCHEME
!LL
!LL  CALLED FROM: CONVECT
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL  DOCUMENTATION : http://hc0500/~hadag/cmt_param.ps.gz
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE CMT_MASS(np_field,NCONV,NLEVS,NTERM,CU_TERM,MB,        &
                          P_0DEGC,KTERM,PLCL,PTOP,N_0DEGC,NLCL,NTOP,    &
                          PHALF,P,cu_tend,                              &
                                ! OUTPUT ARGUMENTS
                          MASS_UP,MASS_DWN,VISC)

      Use cv_run_mod, Only:                                             &
          deep_cmt_opt

      IMPLICIT NONE
!
! VARIABLES THAT ARE INPUT
!
      INTEGER, intent(in) :: &
              NCONV,                                                    &
                               ! NUMBER OF CONVECTING POINTS
     &        np_field,                                                 &
     &        NLEVS,                                                    &
                               ! NUMBER OF MODEL LEVELS
     &        NTERM,                                                    &
                               ! NUMBER OF POINTS TERMINATING
     &        KTERM(NP_FIELD),                                          &
                               ! TERMINATING LEVELS
     &        CU_TERM(NTERM),                                           &
                               ! INDICES FOR TERMINATING POINTS
     &        cu_tend(nterm),                                           &
     &        N_0DEGC(NCONV),                                           &
                               ! LEVEL CORRESPONDING TO MELTING LEVEL
     &        NLCL(NCONV),                                              &
                               ! LIFTING CONDENSATION LEVEL
     &        NTOP(NCONV)      ! TOP LEVEL OF CONVECTION
!
      REAL MB(NCONV),                                                   &
                               ! CLOUD BASE MASS FLUX
     &     P_0DEGC(NCONV),                                              &
                               ! PRESSURE OF MELTING LEVEL (HPA)
     &     PLCL(NCONV),                                                 &
                               ! PRESSURE OF LCL (HPA)
     &     PTOP(NCONV),                                                 &
                               ! PRESSURE AT TOP OF CONVECTION (HPA)
     &     PHALF(NLEVS,NCONV),                                          &
                               ! PRESSURE ON MODEL HALF LEVELS (HPA)
     &     P(NLEVS,NCONV)      ! PRESSURE ON MODEL LEVELS (HPA)
!
! VARIABLES THAT ARE OUTPUT
!
      REAL MASS_UP(NLEVS,NCONV),                                        &
                                    ! UPDRAUGHT MASS FLUX PROFILE
     &     MASS_DWN(NLEVS,NCONV),                                       &
                                    ! DOWNDRAUGHT MASS FLUX
     &     VISC(NLEVS,NCONV)        ! VISCOSITY
!
! VARIABLES THAT ARE LOCAL
!
      INTEGER I,J,K,N
!
      REAL BETA_L,                                                      &
     &     BETA_U,                                                      &
     &     A_M(NTERM),                                                  &
     &     ALPHA_UP,ALPHA_DWN,DELTA_P                                   &
     & ,   shape_factor
!
! PARAMETERS
!
      REAL A_MAX  ,                                                     &
                     ! RATIO MASS FLUX AT MELTING LEVEL TO CLOUD BASE FL
     &     A_TOP,                                                       &
                   ! RATIO MASS FLUX AT TOP TO CLOUD BASE FLUX
     &     P_REF,                                                       &
                     ! USED TO SCALE RATIOS AS P_ODEGC INCREASES
     &     P_MIN,                                                       &
     &     A_DWN,                                                       &
     &     TOP_PRESS,                                                   &
     &     ALPHA_VISC
!
      PARAMETER (A_MAX=1.5,A_TOP=0.2,A_DWN=0.0,P_REF=60000.0,           &
                P_MIN=10000.0,TOP_PRESS=15000.0,ALPHA_VISC=0.30)
!
! DETERMINE SHAPE OF MASS FLUX PROFILE. IF THE DEPTH OF CONVECTION BELOW
! THE ZERO DEGREE LEVEL IS LARGE ENOUGH, AND THE DEPTH ABOVE IS AS WELL,
! THE MASS FLUX PROFILE IS ASSUMED TO HAVE A MAXIMUM AT THE ZERO DEGREE
! LEVEL
!
      DO I=1,NTERM
       J=CU_TERM(I)
       IF((P_0DEGC(J)-PTOP(J)) >= P_MIN.AND.                            &
     &     (PLCL(J)-P_0DEGC(J)) >= P_MIN) THEN
!
! MASS FLUX PROFILE HAS AN ELEVATED MAXIMUM THE VALUE OF THE MAXIMUM
! MASS FLUX INCREASES WITH DECREASING PRESSURE OF THE ZERO DEGREE LEVEL
! WHILE IT IS BELOW P_REF. ABOVE P_REF THE MAXIMUM VALUE IS FIXED
!
        IF(P_0DEGC(J) >= P_REF) THEN
         A_M(I)=1.0+(A_MAX-1.0)*(P_0DEGC(J)-(PLCL(J)-P_MIN))/           &
     &                          (P_REF-(PLCL(J)-P_MIN))
        ELSE
         A_M(I)=A_MAX
        ENDIF
       ELSE
        A_M(I)=1.0
       ENDIF
      END DO
      DO I=1,NTERM
       J=CU_TERM(I)
       N=CU_TEND(I)
       MASS_UP(NLCL(J),J)=MB(J)
       MASS_DWN(NLCL(J),J)=A_DWN*MB(J)
       BETA_L=ALOG(A_M(I))
       BETA_U=ALOG(A_M(I)/A_TOP)
       IF((P_0DEGC(J)-PTOP(J)) >= P_MIN.AND.                            &
     &     (PLCL(J)-P_0DEGC(J)) >= P_MIN) THEN
        DO K=NLCL(J)+1,NTOP(J)+1
         IF(PHALF(K,J) >= P_0DEGC(J)) THEN
          MASS_UP(K,J)=A_M(I)*MB(J)*                                    &
     &                 EXP(-BETA_L*((PHALF(K,J)-P_0DEGC(J))/            &
     &                              (PLCL(J)-P_0DEGC(J)))**2)
          MASS_DWN(K,J)=A_DWN*MB(J)!*(P_0DEGC(J)-PHALF(K,J))/
!     +                              (P_0DEGC(J)-PLCL(J))
         ELSE
          MASS_UP(K,J)=A_M(I)*MB(J)*                                    &
     &                 EXP(-BETA_U*((PHALF(K,J)-P_0DEGC(J))/            &
     &                              (PTOP(J)-P_0DEGC(J)))**2)
          MASS_DWN(K,J)=0.0
         ENDIF
        END DO
       ELSE
        DO K=NLCL(J)+1,NTOP(J)+1
!
! CHOOSE EITHER THE DIAGNOSED TOP OF CONVECTION OR THE LEVEL AT
! WHICH SCHEME DETRAINS, WHICHEVER IS HIGHER (THIS PREVENTS ODD FAILURES
! IN THE CMT SCHEME AT MID-LATITUDES
!
         IF(PTOP(J) >= PHALF(KTERM(N)+1,J)) THEN
          MASS_UP(K,J)=A_M(I)*MB(J)*                                    &
     &                 EXP(-BETA_U*((PHALF(K,J)-PLCL(J))/               &
     &                              (PHALF(KTERM(N)+1,J)-PLCL(J)))**2)
         ELSE
          MASS_UP(K,J)=A_M(I)*MB(J)*                                    &
     &                 EXP(-BETA_U*((PHALF(K,J)-PLCL(J))/               &
     &                              (PTOP(J)-PLCL(J)))**2)
         ENDIF
         MASS_DWN(K,J)=0.0
        END DO
       ENDIF
      END DO
!
! CALCULATE EDDY VISCOSITY (NOTE VISCOSITY=0 AT NLCL AND KTERM+2
!
      DO I=1,NTERM
       J=CU_TERM(I)
       VISC(NLCL(J),J)=0.0
       VISC(NTOP(J)+2,J)=0.0
       DELTA_P=-(PTOP(J)-PLCL(J))
!
! ALPHA_UP IS A VERY CRUDE REPRESENTATION OF THE NON-DIMENSIONAL PROFILE
! OF UPDRAUGHT VERTICAL VELOCITY. TUNING OF THE VISCOSITY WAS DONE IN TH
! SCM, COMPARING MOMENTIUM FLUXES WITH THOSE DERIVED FROM CRM SIMULATION
! OF PERIODS DURING TOGA-COARE.
! ALLOWS FRO POSSIBILITY OF A DOWNDRAUGHT COMPONENT, BUT SET TO ZERO AT
! PRESENT.
!
        If (deep_cmt_opt == 0) then
 
          Do k=nlcl(j)+1,ntop(j)+1
            alpha_up=2.0*(phalf(k,j)-plcl(j))/(top_press-plcl(j)) 
            If(p_0degC(j) < Plcl(j) .and. phalf(k,j) > p_0degc(j)) then
              alpha_dwn =2.0*sqrt((phalf(k,j)-p_0degc(j))/(plcl(j)-p_0degc(j)))
            Else
              alpha_dwn =0.0
            End If 
            visc(k,j)=(alpha_up*mass_up(k,j)+alpha_dwn*mass_dwn(k,j))*  &
                        alpha_visc*delta_p
          End Do

        Else if (deep_cmt_opt == 1) then

! Added a shape factor to control gradient term
! Note as a_dwn=0.0 removed downward terms from calculation as waste of CPU

          Do k=nlcl(j)+1,ntop(j)+1
            alpha_up=2.0*(phalf(k,j)-plcl(j))/(top_press-plcl(j)) 
            shape_factor = 1.0-0.75*(phalf(k,j)-plcl(j))/(ptop(j)-plcl(j)) 
            visc(k,j)=alpha_up*mass_up(k,j)*alpha_visc*delta_p*shape_factor 
          End Do

        Else
          Write(6,*) 'CMT_MASS: Unacceptable value of deep_cmt_opt ',   &
                      deep_cmt_opt 
        End If

      END DO

      RETURN
      END SUBROUTINE CMT_MASS
#endif
