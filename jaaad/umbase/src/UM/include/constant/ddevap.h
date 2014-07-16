#if defined(A05_4A) || defined(A05_5A)
! DDEVAP start
      ! exponents used in calculation of evaporation of liquid
      REAL,PARAMETER :: P_LQ1 = 0.52
      REAL,PARAMETER :: P_LQ2 = 0.67

      ! exponents used in calculation of evaporation of ice
      REAL,PARAMETER :: P_ICE1 = 0.55
      REAL,PARAMETER :: P_ICE2 = 0.76

      ! exponents and constants associated with density term in
      ! evaporation of liquid
      REAL,PARAMETER :: RHO_LQP1 = 0.26
      REAL,PARAMETER :: RHO_LQP2 = 0.59
      REAL,PARAMETER :: RHO_LQA  = 108.80
      REAL,PARAMETER :: RHO_LQB  = 830.73

      ! exponents and constants associated with density term in
      ! evaporation of ice
      REAL,PARAMETER :: RHO_ICP1 = 0.28
      REAL,PARAMETER :: RHO_ICP2 = 0.63
      REAL,PARAMETER :: RHO_ICEA = 1569.52
      REAL,PARAMETER :: RHO_ICEB = 32069.02
! DDEVAP end
#endif
