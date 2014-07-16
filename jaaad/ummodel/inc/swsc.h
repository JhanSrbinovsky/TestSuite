!      REAL SC                           !  Solar constant
!      PARAMETER ( SC = 1365. )
!20110812: define SC in coupling control
#if defined(OASIS3)
!
! Now get SC from namelist of coupling control
       use auscom_cpl_data_mod, only: SC,VOLCTS_val
#endif
