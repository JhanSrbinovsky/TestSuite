#if defined(A05_4A) || defined(A05_5A) || defined(C70_1A)
! ENTCNST start

!  6.2   28/02/05  Modifications for STPH_RP
      ! coefficients used in calculation of entrainment rate
      REAL,PARAMETER:: AE1     = 1.0
      REAL,PARAMETER:: AE2     = 1.5
      REAL:: ENTCOEF
      COMMON/STPH_RP_1/entcoef
      REAL,PARAMETER:: SH_FAC  = 1.0

! ENTCNST end
#endif
