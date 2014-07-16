!*L------------------ COMDECK P_EXNERC ---------------------------------
!statement function to define exner pressure at model levels
!from exner pressures at model half-levels
      REAL P_EXNER_C
      REAL R_P_EXNER_C
      REAL P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM
#if defined(ABCALC4)
!arithmetic mean
      P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =          &
     & 0.5*(P_EXU_DUM + P_EXL_DUM)
      R_P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =        &
     & 2.0/(P_EXU_DUM + P_EXL_DUM)
#else
!consistent with geopotential see eqn 26 DOC PAPER 10
      P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =          &
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)/                           &
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )
      R_P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =        &
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )/                             &
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)
#endif

!*------------------- --------------------------------------------------
