! AOCPMA; for ocean assimilation; documented in CODPMA
!
! History:
! Version  Date     Comment
! -------  ----     -------
!  5.5 27/02/03 Changes for upgrade to FOAM data assimilation. A. Hines
!  6.1 02/09/04 Implement Pressure Correction scheme. M. Martin.
!
     &  LLO_STR_EDG, JO_MOD_GET, JO_TYP_GA, JO_COL_GA, JO_ROW_GA,       &
     &  JO_LEV_GA, JO_ANA_GA, JO_CNV_FNC, JO_GEO_DEP, JO_GEO_LAT,       &
     &  IO_ADD_A, IO_ADD_B, O_COF_GA, O_NO_COR_SCL, O_COF_GEO,          &
     &  O_COF_GEO_DEP, O_COF_GEO_LAT, O_NON_DIV_COR_COF, JO_IT_NO,      &
     &  JO_COR_TIM_FNC, JO_WTS_TIM_FNC, JO_INC_TIM_FNC, KO_LEV_LOW,     &
     &  JO_SST_ML_LEV,                                                  &
     &  O_TIM_WIN_STR, O_TIM_WIN_FIN, LLO_FILTER, JO_NPASS_RF,          &
     &  O_CS_AvVal_WE, O_CS_AvVal_NS,  O_CS_LatRef_WE,   O_CS_LatRef_NS,&
     &  O_CS_LonRef_WE,O_CS_LonRef_NS, O_CS_DepRef_WE,   O_CS_DepRef_NS,&
     &  O_CS_PAmp_WE,  O_CS_PAmp_NS,   O_CS_PScl_WE,     O_CS_PScl_NS,  &
     &  JO_CS_PFn_WE,  JO_CS_PFn_NS,                                    &
     & THRM_FACT, LLO_THRMCLINE, O_Nudge_Weight, LLO_SAT_BIAS,          &
     & JO_NO_OUT_ITS,O_VRT_CS,                                          &
     & WEIGHT_BIAS, MIN_LAT_BIAS, MAX_LAT_BIAS, BIAS_DECAY,             &
! AOCPMA end
