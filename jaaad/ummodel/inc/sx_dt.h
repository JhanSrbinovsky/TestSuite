!========================== COMDECK SX_DT =======================
! Need to call sx_size or typsize BEFORE this deck
! Description:
!
!     This COMDECK contains data initialisation for variables in
!     SX.
!

      DATA LAND_FIELD / -99 /
      DATA TR_VARS / -99 /
#if !defined(MAKEBC)
      DATA TOT_LEVELS / -99 /
#endif
! ---------------------- End of comdeck SX_DT -------------------------
