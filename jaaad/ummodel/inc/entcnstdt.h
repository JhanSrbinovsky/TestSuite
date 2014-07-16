! ------------------------------------------------------------------
! Description
!   This comdeck includes the initialization of entcoef for STPH_RP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.2   18/07/06   New include for the
!                    Stochastic Physics Random Parameters scheme
!                    (STPH_RP). A. Arribas
!
! ------------------------------------------------------------------
!

#if defined(A05_4A) || defined(A05_5A) \
                    || defined(C70_1A)
      DATA ENTCOEF /3.0/
#endif
