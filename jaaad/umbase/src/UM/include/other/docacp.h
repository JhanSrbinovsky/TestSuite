!L Description:
!L----------------------------------------------------------------------
!L  COMDECK DOCACP
!L  --------------
!L  DOCACP is list of variables in comdeck COMACP
!L
!L  COMACP contains parameters controlling the assimilation.
!L
!LHistory:    see DOCACP
!L  UMAC VN 2.8   LWBAL,WB_CUT_OFF_LEV,WB_PSTAR_METHOD added to
!L                control WINDBAL               Phil Andrews 15/9/92
!L                LAC_MES added to set Mes defaults
!L                and MRAMPFN for alternate timeramp S Bell 15/10/92
!L  UMAC VN 3.1   New variable MTHIN205,MVINT205     S Bell 18/11/92
!L                LWBAL,WB_CUT_OFF_LEV,WB_PSTAR_METHOD removed
!L                and replaced by WB_THETA_INC,WB_LAT_CC,LWBAL_SF,
!L                LWBAL_UA,WB_VERT_V,WB_LAND_SCALE,WB_LAND_FACTOR
!L                                      S Bell, Phil Andrews 8/1/93
!L  UMAC VN 3.2   New variables WB_THETA_SURF,MGLOSSFN S Bell 23/6/93
!L                MTHINxxx replaced by DEF_OBTHIN
!L  UMAC VN 3.3   New variables THRESH_** for MOPS B Macpherson 1/12/93
!L                New variables for limited area WINDBAL: WB_LonOffset,
!L                WB_LonPts, WB_LatOffset, WB_LatPts.
!L                WB_THETA_INC renamed WB_THETA_UA and
!L                WB_THETA_SURF renamed WB_THETA_SF Phil Andrews10/11/93
!L  UMAC VN 3.4   Remove L_REINITQC
!L                New variable LHYDROL for hydrology correction scheme
!L                New variables NRADARS,LRADAR,RADAR_RANGE,
!L                L_MOPS_EQUALS_RH,VERT_COR_AERO,
!L                L_LATLON_PRVER,NORTH/SOUTHLAT,WEST/EASTLON
!L                                                Bruce M  9/9/94
!L  UMAC VN 3.5   New variables LCHECK_GRID Greg B 21/2/95
!L  UMAC VN 4.0   Move variables RADAR_LAT, RADAR_LON from VANRAIN
!L                New variables RADAR_RANGE_MAX, F1_506, F2_506,
!L                F3_506, L_506_OBERR for 506 ob error calculation
!L                New variables EPSILON_LHN, LHN_RANGE, L_LHN,
!L                L_LHN_SCALE, L_LHN_SEARCH, LHN_DIAG, RELAX_CF_LHN.
!L                                            Chris Jones (3/5/95)
!L  UMAC VN 4.0   Add variable L_OBS_CHECK G.Bason
!L  UMAC VN 4.0   Change 'OBS_FORMAT' comment G.Bason
!L  UMAC VN 4.0   Add NON_DIV_COR_10M
!L                             Bruce M 21/3/95
!L  UMAC VN 4.1  : Introduce new LHN variables: ALPHA_LHN,
!L               :     LHN_LIMIT, FI_SCALE_LHN, NPASS_RF_LHN,
!L               :     L_LHN_LIMIT, L_LHN_FACT, L_LHN_FILT (Chris Jones)
!   UMAC VN 4.3  :Redefine CSCFACT_V & NSLABS_CSFACT (SB)
!   UMAC VN 6.2  : Introduce REMOVE_NEG_LH Mark Dixon 23/02/06
!L----------------------------------------------------------------------
!L
!L  AC_OBS_TYPES      - List of AC Obs Types to be processed.
!L                    - Order of processing is stored in LACT.
!L  AC_ORDER          - Order which AC Obs Types MUST be processed.
!L                    - (Coded as group*1000+Obstype)
!L  CSCFACT_H         - Horizontal Correlation Scale Factor.
!L                    - Varies with latitude.
!L  CSCFACT_V         - Vertical Correlation Scale Factor.
!L                    - Varies with level
!L                    - NSLABS_SCFACT(J-1)+1 TO NSLABS_SCFACT(J)
!L  DEF_AC_ORDER      - Default Order and groupings of Obs Types.
!L  DEF_CSCALE_START  - Default Correlation Scale at start of
!L  DEF_CSCALE_OBTIME - insertion period/observation time/end of
!L  DEF_CSCALE_END    - insertion period for each group. Use
!L                    - CSCALE_START/OBTIME/END to change defaults.
!L  DEF_MODE_HANAL    - Default Mode of Horizontal Analysis.
!L  DEF_FI_VAR_FACTOR - Default group dep scaling in FI
!L  DEF_NO_ANAL_LEVS  - Default No of Analysis Levels for each group.
!L  DEF_NO_WT_LEVS    - Default No of Weight Levels for each group.
!                    - Use N_ANAL_LEVS/N_WT_LEVS to change defaults.
!L  DEF_NO_ITERATIONS - Default No of Iterations for groups. Use
!L                    - NO_ITERATIONS in namelist to change defaults.
!L  DEF_INTERVAL_ITER - Default No of Iterations for groups. Use
!L                    - INTERVAL_ITER in namelist to change defaults.
!L  DEF_OBTHIN        - Default ob thinning (use OBTHIN in namelist)
!L                    - (values of 1 imply no thinning, N implies
!L                    - 1/N reports assimilated every Nth step)
#if defined(GLOBAL)
!L  DEF_NUDGE_NH      - Default Nudging Coeffs for NH for groups.
!L  DEF_NUDGE_TR      - Default Nudging Coeffs for TR for groups.
!L  DEF_NUDGE_SH      - Default Nudging Coeffs for SH for groups.
!                    - Use NUDGE_NH/TR/SH in namelist to change
!                    - defaults.
#else
!L  DEF_NUDGE_LAM     - Default Nudging Coeffs for LAM for groups.
!                    - Use NUDGE_LAM in namelist to change defaults.
#endif
!L  DEF_RADINF        - Default Max Normalised Influence Radius.
!L                    - Use RADINF in namelist to change defaults.
!L  DEF_TGETOBB )     - Default Time Window before/after obs time to
!L  DEF_TGETOBA )     - fetch obs for groups. Use TGETOBB/TGETOBA in
!L                    - namelist to change deafults.
!L  DEF_TIMEB )       - Default Insertion Period before/after obs time
!L  DEF_TIMEA )       - for groups. Use TIMEB/TIMEA in namelist
!L                    - to change defaults.
!L  DF_COEFF          - Coefficient for DIVFILT
!L  DF_SCALE          - DIVFILT scale (metres)
!L  DF_SCALE_LEV      - DIVFILT scale for each level
!L  DIAG_RDOBS        - Diagnostic Control for subroutine RDOBS.
!L  EPSILON_LHN       - Epsilon value for use in LHN
!L  F1_506           \                                                 .
!L  F2_506            } Parameters for 506 ob weight evaluation
!L  F3_506           /
!L  ALPHA_LHN         - Alpha value for use in LHN
!L  LHN_LIMIT         - Limit on + or - Theta incr rate in LHN (K/day)
!L  FI_SCALE_LHN      - Recursive filter scale in m
!L  NPASS_RF_LHN      - Number of passes through filter
!L  FI_SCALE          - FI (Filtered Increments) scale (metres)
!L  FI_SCALE_FACTOR   - FI Scale Factor
!L  GEOWT_H           - Latitude weights for Geostrophic Increments.
!L  GEOWT_V           - Vertical weights for Geostrophic Increments.
!L  GROUP_NO          - Group No of each obs type in LACT.
!L  GROUP_FIRST       - Position in LACT of first type of each group.
!L  GROUP_LAST        - Position in LACT of last  type of each group.
!L  GROUP_INDEX       - Corresponding group in DEF_AC_ORDER for
!L                    - groups in GROUP_NO.
!L  IOMITOBS          - List of Observations not to be assimilated.
!L                    - Use Model Observation Numbers to omit obs.
!L  IUNITNO           - Unit No of Cache file to store obs
!L                    - between timesteps.
!L  L_506_OBERR       - Logical switch to control 506 ob weight calc'n
!L  L_LHN             - Logical switch to perform latent heat nudging
!L  L_LHN_SCALE       - Logical switch to control scaling within LHN
!L  L_LHN_SEARCH      - Logical switch to control use of LHN_SEARCH
!L  L_VERIF_RANGE     - Logical switch to control verification range
!L  L_LHN_LIMIT       - Logical switch to control limiting of increments
!L  L_LHN_FACT        - Logical switch to control limiting by 1/alpha
!L  L_LHN_FILT        - Logical switch to control filtering of incrs
!L  LACT              - List of Obs Types to be processed in order
!L                    - of processing.
!L  LAC_UARS          - Logical switch for UARS assimilation.
!L  LAC_MES           - Logical switch for Mesoscale assimilation.
!L  LCHECK_GRID       - Logical switch to control CHECK_OBS call
!L  LENACT            - No of obs for each type in group fetched
!L                    - this timestep.
!L  LGEO  )           - Logical switches to calculate
!L  LHN_DIAG          - Logical switch for detailed LHN diagnostics
!L  LHN_RANGE         - Max search radius used in LHN_SEARCH. Do not set
!L                      to zero, set L_LHN_SEARCH=.FALSE. instead.
!L  LHYDR )           - Geostrophic/Hydrstatic Increments.
!L  LHYDROL           - Logical switch to calc hydrology incrs.
!L  LRADAR            - Logical array to determine which radars to use
!L  L_LATLON_PRVER    - Logical switch to verify precip in lat/lon area
!L    NORTHLAT        - Co-ords in
!L    SOUTHLAT        -            real lat/lon
!L    WESTLON         -                         for rain
!L    EASTLON         -                                  verif area.
!L  L_MOPS_EQUALS_RH  - If .TRUE. then MOPS cloud obs are
!L                    - rh values (%), else cloud fractions
!L  L_OBS_CHECK       - If .FALSE. then skip check to see if there
!L                    - are any obs to assimilate (non-oper run only)
!L  LSYN              - Logical switch for Synoptic Insertion.
!L  LWBAL_SF          - Controls use of WINDBAL routine for surface wind
!L  LWBAL_UA          - Controls use of WINDBAL routine for uair wind
!L  MASTER_AC_TYPES   - Master list of AC Observation Types
!L                    - known to AC Scheme. See DEF_TYPE.
!L  MACDIAG           - Diagnostics control for each timestep.
!L                    - See AC on use of this array.
!L  MDATADFN          - Mode for Data Density Formula.
!L  MGEOWT            - Mode of Latitude Weighting for
!L                    - Geostrophic Increments.
!L  MGLOSSFN          - GLOSS processing Function Option. (see VANLASS)
!L  MHCORFN           - Correlation Function Option.
!L  MODEACP           - No of timesteps allowed in dimensioning of
!L                    - MACDIAG. Code loops back to MACDIAG(1) if
!L                    - TIMESTEP_NO.GT.MODEACP
!L  MVINT205          - Options for vertical interp    (LASS/GLOSS  )
!L  MRAMPFN           - Mode for Time ramp in HORINF
!L  MWTFN             - Mode for Weights Formula.
!L  NACT              - No of obs types in LACT.
!L  N_GROUPS          - No of groups being processed.
!L  NO_OBS_FILES      - No of observation files to be used.
!L  NO_SCFACT         - List of obs types on which correlation
!L                    - scale factor is not to be applied.
!L  NON_DIV_COR       - Factor for Non-Divergent Correction.
!L  NON_DIV_COR_10M - As NON_DIV_COR but for 10m wind data
!L  NPASS_RF          - Number of passes in RFILT
!L  NPROG             - Number set by individual programmers
!L                    - doing test work. Numbers in use : 1001 (DR)
!L  NRADARS           - No of radars in network
!L  NSLABS_SCFACT     - Slab for each level
!L                    - (Slab is group of levels with same CSCFACT_V)
!L  OBS_FORMAT        - Format of AC Obs file (=2, only one format)
!L  OBS_UNITNO        - Unit No of first AC Obs file (=70)
!L  OBTIME_NOM        - Nominal Observation Time for Synoptic Insertion
!L                    - Mode. Relative time from start of assimilation.
!L  RADAR_LAT         - Coordinates of radars
!L  RADAR_LON         -     "        "    "
!L  RADAR_RANGE       - Max range (km) of reliable radar rain rates
!L  RADAR_RANGE_MAX   - Max. range of radar data used for LHN (km)
!L  RELAX_CF_LHN      - Relaxation coef used for theta incrs in LHN
!L  REMOVE_NEG_LH - Logical switch for removing -ve LH values
!L                                      Mark Dixon 23/02/06
!L  SPEED_LIMIT305    - Min speed of scatwinds for which observed
!L                    - direction is used. (below limit speed only
!L                    - is assimilated.
!L  THRESH_DL         - threshold (mm/hr) between dry/light rain
!L  THRESH_LM         - threshold (mm/hr) between light/moderate rain
!L  THRESH_DL         - threshold (mm/hr) between moderate/heavy rain
!L  THRESH_RMSF       - threshold (mm/hr) for calcn of rms factor score
!L  TIMEF_START       - Start    ) Values for
!L  TIMEF_OBTIME      - Obs time ) Time Factor
!L  TIMEF_END         - End      ) Ramp Function
!L  TROPLAT           - Latitude at which parameters start to change
!L                    - from their mid-latitude to tropical values.
!L  TROPINT           - Interval over which parameters start to change
!L                    - from their mid-latitude to tropical values.
!L  TYPE_INDEX        - Position of obs types in LACT in MASTER_AC_TYPES
!L  USE_CONV_IN_MOPS  - Logical switch for using convection in MOPS
!L  VERT_CUTOFF_SL    - No of scale hts over which single level obs used
!L  VERT_CUTOFF_BW    - as VERT_CUTOFF_SL but for bogus winds
!L  VERT_CUTOFF_BH    - as VERT_CUTOFF_SL but for bogus humidity
!L  VERT_COR_AERO     - Vertical correlation scale for aerosol incrs
!L  VERT_COR_SCALE    - Vertical Correlation Scale Coefficient.
!L                    - calculated from ACP namelist array CSCALE_VERT
!L                      (n,1)=extra tropical temps,(n,2)=tropical temps
!L                      (n,3)=extra tropical winds,(n,4)=tropical winds
!L  VERT_FILT         - Vertical Filtering of Increments
!L                    - from soundings.
!L  WB_THETA_UA       - If T WINDBAL will calc theta incs from upper air
!L                    - winds.
!L  WB_THETA_SF       - If T WINDBAL will calc theta incs from surface
!L                    - winds.
!L  WB_LAT_CC         - horizontal correlation coeff for WINDBAL
!L  WB_VERT_V         - vertical variation in WINDBAL correlations
!L  WB_LAND_FACTOR    - extra scaling factor of WINDBAL inc over land
!L  WB_LAND_SCALE     - apply WB_LAND_FACTOR scaling if true
#if !defined(GLOBAL)
!L  WB_LonOffset      -) These define a subset of a limited area model
!L  WB_LonPts         -) within which WINDBAL will work. They exist to
!L  WB_LatOffset      -) select a region on which a multigrid Poisson
!L  WB_LatPts         -) solver can be used efficiently. Offsets are
!L                    -) from start of full LAM grid (so offsets of
!L                    -) zero mean no offset). WB_LonPts & WB_LatPts
!L                    -) define the length of the subset in their
!L                    -) respective directions.
#endif
!L
!-----------------------------------------------------------------------
