#if defined(ATMOS) || defined (A34_1A)
! ARGPTRA start
! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  05/10/94  Add murk and user ancillary pointers. RTHBarnes
!  4.1  04/12/95  Add pointers JSTHU and JSTHF. J.Smith
!  4.1  26/04/96  Add pointers for Sulphur Cycle variables (12)  MJW
!  4.3   18/3/97  And for HadCM2 sulphate loading patterns.  Will Ingram
!  4.4   05/8/97  And for Conv. cloud amt on model levs. Julie Gregory
!  4.5  04/03/98   Remove pointer SO2_HILEM (add to CARGPT_ATMOS)
!                  Add 1 pointers for NH3 in S Cycle
!                  Add 3 pointers for Soot              M. Woodage
!  4.5  08/05/98   Add 16 new pointers for User Anc.    D. Goddard
!  4.5  15/07/98   Add pointers for new 3D CO2 array.   C.D.Jones
!  4.5  17/08/98   Remove JSOIL_FLDS and JVEG_FLDS      D. Robinson
!  5.0    20/05/99 Extensive revision of pointers for new dynamics data
!                  variables. D.M. Goddard
!  5.1  28/02/00   Added ND LBC variables            P.Burton
!  5.2  31/08/00   Added orography LBC and changed order of LBC pointers
!                                                               P.Burton
!  5.2  25/09/00  Clear out RHcrit legacy code.    A.C.Bushell.
!
!  5.2    13/09/00 Remove JSOOT, JSO4, JH2SO4 pointers. P.Selwood.
!  5.3  19/06/01   Add pointer (JTPPSOZONE) for
!                  tropopause-based ozone        D. Tan
!  5.5  01/11/02   Add pointers JQCF2, JQRAIN, and JQGRAUP
!                  and associated LBCs for microphysics    R.M.Forbes
!  6.1  06/11/03   Add pointers for STOCHEM fields
!  6.0  30/07/03   Add pointers for cloud fractions and cloud
!                  fraction lbcs. Damian Wilson
! 6.1  13/08/04  Reduce number of continuation lines - prevent coupled
!                models hitting 99 line limit. R. Hill
!  6.2  01/03/06   Add pointer for direct component of total downward
!                  surface PAR flux.  M.G. Sanderson.
!  6.2  11/05/05   Put STOCHEM pointers in correct place.
!                      M.G. Sanderson
!  6.2  14/11/05   Add JTR_UKCA pointers for UKCA tracers. R Barnes
!  6.2  01/10/05   Add JMURK_LBC pointers.    R.M.Forbes
      ! Pointers for ATMOSPHERE model variables. Configuration dependent.

      ! Addresses in D1 array of primary variables held in primary and
      ! secondary space

      !  Array  variables (dimensions are resolution dependent)

      ! Data variables stored in primary space.
     &  JU, JV, JW, JRHO, JTHETA, JQ, JQCL, JQCF,                       &
     &  JEXNER_RHO_LEVELS, JU_ADV, JV_ADV, JW_ADV,                      &
      ! Data variables stored in secondary space.
     &  JP, JP_THETA_LEVELS,  JEXNER_THETA_LEVELS,                      &
      ! Cloud Fields
     &  JCCA, JCF_AREA, JCF_BULK, JCF_LIQUID, JCF_FROZEN,               &
      ! Soil Ancillary Fields
     &  J_DEEP_SOIL_TEMP,  JSMCL, JSTHU, JSTHF,                         &
      ! Radiation increments
     &  JSW_INCS, JLW_INCS,                                             &
      ! Ozone
     &  JOZONE,                                                         &
      ! Tracers and Aerosols
     &  JTRACER, JMURK, JMURK_SOURCE,                                   &
     &  JSO2, JDMS, JSO4_AITKEN, JSO4_ACCU, JSO4_DISS, JH2O2,           &
     &  JNH3, JSOOT_NEW, JSOOT_AGD, JSOOT_CLD, JSO2_NATEM,              &
     &  JOH, JHO2, JH2O2_LIMIT,JO3_CHEM, JHadCM2_SO4, JCO2,             &
      ! User Ancillary fields
     &  JUSER_MULT1, JUSER_MULT2, JUSER_MULT3, JUSER_MULT4,             &
     &  JUSER_MULT5, JUSER_MULT6, JUSER_MULT7, JUSER_MULT8,             &
     &  JUSER_MULT9, JUSER_MULT10, JUSER_MULT11, JUSER_MULT12,          &
     &  JUSER_MULT13, JUSER_MULT14, JUSER_MULT15, JUSER_MULT16,         &
     &  JUSER_MULT17, JUSER_MULT18, JUSER_MULT19, JUSER_MULT20,         &
      ! Lateral Boundary Conditions and tendencies
     &  JOROG_LBC,JU_LBC,JU_LBC_TEND,JV_LBC,JV_LBC_TEND,JW_LBC,         &
     &  JW_LBC_TEND,JRHO_LBC,JRHO_LBC_TEND,JTHETA_LBC,JTHETA_LBC_TEND,  &
     &  JQ_LBC,JQ_LBC_TEND,JQCL_LBC, JQCL_LBC_TEND,JQCF_LBC,            &
     &  JQCF_LBC_TEND,JEXNER_LBC,JEXNER_LBC_TEND,JU_ADV_LBC,            &
     &  JU_ADV_LBC_TEND,JV_ADV_LBC,JV_ADV_LBC_TEND,JW_ADV_LBC,          &
     &  JW_ADV_LBC_TEND,JTRACER_LBC,JTRACER_LBC_TEND,                   &
! Tropopause-based Ozone
     &  JTPPSOZONE,                                                     &
      ! Biomass aerosol
     &  JBMASS_NEW, JBMASS_AGD, JBMASS_CLD,                             &
      ! Additional microphysics fields and lbcs
     &  JQCF2,JQRAIN,JQGRAUP,JQCF2_LBC,JQCF2_LBC_TEND,JQRAIN_LBC,       &
     &  JQRAIN_LBC_TEND,JQGRAUP_LBC,JQGRAUP_LBC_TEND,JDUST_DIV1,        &
     &  JDUST_DIV2,JDUST_DIV3,                                          &
     &  JDUST_DIV4,JDUST_DIV5,JDUST_DIV6,                               &
     &  JCF_BULK_LBC,JCF_BULK_LBC_TEND,JCF_LIQUID_LBC,                  &
     &  JCF_LIQUID_LBC_TEND,JCF_FROZEN_LBC,JCF_FROZEN_LBC_TEND,         &
! Fields from STOCHEM
     &  JCH4_STOCH,JO3_STOCH,                                           &
! Pointer for direct PAR flux
     &  JDIRPAR,                                                        &
     &  JTR_UKCA, JMURK_LBC, JMURK_LBC_TEND,                            &
! Pointers for UKCA oxidant fields
     &  JOH_UKCA, JHO2_UKCA, JH2O2_UKCA, JO3_UKCA,                      &
! Convective Cloud Fields
     &  JLCBASE, JCCW_RAD,                                              &
! Ozone tracer and cariolle parameters
     &  JOZONE_TRACER,JO3_PROD_LOSS,JO3_P_L_VMR,JO3_VMR,JO3_P_L_TEMP,   &
     &  JO3_TEMP,JO3_P_L_COLO3,JO3_COLO3,                               &
! Pointers for Aerosol climatologies
     &  JARCLBIOG_BG, JARCLBIOM_FR, JARCLBIOM_AG, JARCLBIOM_IC,         &
     &  JARCLBLCK_FR, JARCLBLCK_AG, JARCLSSLT_FI, JARCLSSLT_JT,         &
     &  JARCLSULP_AC, JARCLSULP_AK, JARCLSULP_DI, JARCLDUST_B1,         &
     &  JARCLDUST_B2, JARCLDUST_B3, JARCLDUST_B4, JARCLDUST_B5,         & 
     &  JARCLDUST_B6, JARCLOCFF_FR, JARCLOCFF_AG, JARCLOCFF_IC,         &
     &  JARCLDLTA_DL,                                                   &
! Fossil-fuel organic carbon aerosol
     &  JOCFF_NEW, JOCFF_AGD, JOCFF_CLD,                                &
      ! CABLE
     &  JTSOIL_TILE, JSMCL_TILE, JSTHF_TILE, JSNOW_DEPTH3L,             &
     &  JSNOW_MASS3L, JSNOW_TMP3L, JSNOW_RHO3L, JSNOW_RHO1L, JSNOW_AGE, & 
     &  JSNOW_FLG3L,                                                    &
! Lestevens Sept 2012 - adding new progs
     &  JCPOOL_TILE,JNPOOL_TILE,JPPOOL_TILE,JSOIL_ORDER,                &
     &  JNIDEP,JNIFIX,JPWEA,JPDUST,JGLAI,JPHENPHASE,                    &
     &  JCABLE_LAI,                                                     &
! ARGPTRA end
#endif
