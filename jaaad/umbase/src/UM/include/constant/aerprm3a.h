! AERPRM3A aerosol parametrizations for two-stream radiation code.

      INTEGER,PARAMETER:: IP_AEROSOL_PARAM_DRY=1
      INTEGER,PARAMETER:: IP_AEROSOL_PARAM_MOIST=2
      INTEGER,PARAMETER:: IP_AEROSOL_UNPARAMETRIZED=3 ! Observational
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      INTEGER,PARAMETER:: IP_AEROSOL_PARAM_PHF_DRY=4
      INTEGER,PARAMETER:: IP_AEROSOL_PARAM_PHF_MOIST=5
#endif

! AERPRM3A end
