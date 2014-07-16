! ROCPMA for ocean assimilation; documented in CODPMA
!
! History:
! Version  Date     Comment
! -------  ----     -------
!  5.5 27/02/03 Changes for upgrade to FOAM data assimilation. A. Hines
!  6.1 02/09/04 Implement Pressure Correction scheme. M. Martin.
!
      REAL :: O_COF_GA
      REAL :: O_NO_COR_SCL
      REAL :: O_COF_GEO
      REAL :: O_COF_GEO_DEP
      REAL :: O_COF_GEO_LAT
      REAL :: O_NON_DIV_COR_COF
      REAL :: O_TIM_WIN_STR
      REAL :: O_TIM_WIN_FIN
      REAL :: O_CS_AvVal_WE(NO_SPL_CV)
      REAL :: O_CS_AvVal_NS(NO_SPL_CV)
      REAL :: O_CS_LatRef_WE(NO_SPL_CV)
      REAL :: O_CS_LatRef_NS(NO_SPL_CV)
      REAL :: O_CS_LonRef_WE(NO_SPL_CV)
      REAL :: O_CS_LonRef_NS(NO_SPL_CV)
      REAL :: O_CS_DepRef_WE(NO_SPL_CV)
      REAL :: O_CS_DepRef_NS(NO_SPL_CV)
      REAL :: O_CS_PAmp_WE(NO_CSF_MAX,NO_SPL_CV)
      REAL :: O_CS_PAmp_NS(NO_CSF_MAX,NO_SPL_CV)
      REAL :: O_CS_PScl_WE(NO_CSF_MAX,NO_SPL_CV)
      REAL :: O_CS_PScl_NS(NO_CSF_MAX,NO_SPL_CV)
      REAL :: O_VRT_CS(NO_SPL_CV)
      REAL :: THRM_FACT
      REAL :: O_Nudge_Weight
      REAL :: WEIGHT_BIAS
      REAL :: MIN_LAT_BIAS
      REAL :: MAX_LAT_BIAS
      REAL :: BIAS_DECAY

      INTEGER :: JO_MOD_GET
      INTEGER :: JO_TYP_GA
      INTEGER :: JO_COL_GA
      INTEGER :: JO_ROW_GA
      INTEGER :: JO_LEV_GA
      INTEGER :: JO_ANA_GA
      INTEGER :: JO_CNV_FNC
      INTEGER :: JO_GEO_DEP
      INTEGER :: JO_GEO_LAT
      INTEGER :: IO_ADD_A
      INTEGER :: IO_ADD_B
      INTEGER :: JO_IT_NO
      INTEGER :: JO_COR_TIM_FNC
      INTEGER :: JO_WTS_TIM_FNC
      INTEGER :: JO_INC_TIM_FNC
      INTEGER :: JO_SST_ML_LEV
      INTEGER :: KO_LEV_LOW
      INTEGER :: JO_NPASS_RF
      INTEGER :: JO_CS_PFn_WE(NO_CSF_MAX,NO_SPL_CV)
      INTEGER :: JO_CS_PFn_NS(NO_CSF_MAX,NO_SPL_CV)
      INTEGER :: JO_NO_OUT_ITS

      LOGICAL :: LLO_STR_EDG
      LOGICAL :: LLO_FILTER
      LOGICAL :: LLO_THRMCLINE
      LOGICAL :: LLO_SAT_BIAS
! ROCPMA end
