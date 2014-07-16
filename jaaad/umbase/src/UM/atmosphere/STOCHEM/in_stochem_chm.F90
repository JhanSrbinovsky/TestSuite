#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      MODULE IN_STOCHEM_CHM
!
! Modification history from model version 6.0:
!  Model
! version  Date     Comment
!
!  6.0   17/09/03   Correct non-standard F90 continuation lines.
!                   Introduce standard UM modification history.
!                                                         P.Dando
!  6.1   20/08/04   Changes to allow vectorisation on SX6.
!                                                   M. Sanderson.
!  6.2   15/11/05  Extra variables added. M.G. Sanderson
!----------------------------------------------------------------------
! Pick up UM expt and job IDs.
#include "chsunits.h"
#include "cntlall.h"

! Directory names
!kk   Introduced IN_STOCHEM_CHM_GET_PATH to get Names from environment

      CHARACTER(LEN=80) :: datdir, emdir, emdir2

! Array length parameters
      INTEGER, PARAMETER :: nc = 79  ! No. of chemical species
      INTEGER, PARAMETER :: ndj = 17 ! No. Photolysis rates
      INTEGER, PARAMETER :: nr = 295 ! Max. no. reactions
      INTEGER, PARAMETER :: chunk = 1023
      INTEGER, PARAMETER :: chunkm1 = chunk - 1
      INTEGER, PARAMETER :: maxflux = 700 ! Max number of 3D fluxes
      INTEGER, PARAMETER :: n_spec_deposit = 21 ! No of species dry dep

      REAL, PARAMETER :: o2frac = 0.20948 ! o2 molar fraction in air
      REAL, PARAMETER :: n2frac = 0.78084 ! n2 molar fraction in air

      REAL, PARAMETER :: r_null  = 1.0e50 ! Null surface resistance
      REAL, PARAMETER :: r_nodep = 1.0e40
! Surface resistance must be less than
! this value for deposition to occur.

! Chemical names
      CHARACTER(LEN=8), DIMENSION(nc), PARAMETER :: cnames =            &
     &(/'OD      ','OP      ','OH      ','NO      ','NO2     ',         &
     &  'NO3     ','N2O5    ','CO      ','CH4     ','HCHO    ',         &
     &  'O3      ','H2      ','HNO3    ','H2O2    ','CH3O2   ',         &
     &  'HO2     ','C2H6    ','C2H5O2  ','CH3CHO  ','CH3COO2 ',         &
     &  'PAN     ','CH3OOH  ','NC4H10  ','SC4H9O2 ','CH3COE  ',         &
     &  'SO2     ','C2H4    ','C3H6    ','OXYL    ','SA      ',         &
     &  'C3H8    ','C3H7O2  ','C3H7OOH ','C2H5OOH ','C4H9OOH ',         &
     &  'CH3OH   ','ACETONE ','ACETO2  ','AMMSUL  ','CH3COX  ',         &
     &  'CH2O2C  ','MGLYOX  ','CH3CHX  ','GLYOX   ','OXYL1   ',         &
     &  'MEMALD  ','MEMALD1 ','C5H8    ','RO2IP1  ','MVK     ',         &
     &  'RO2IP2  ','ISOPOOH ','MVKOOH  ','TOLUEN  ','RNC2H4  ',         &
     &  'RNC3H6  ','RNC5H8  ','NAER    ','TOLP1   ','DMS     ',         &
     &  'DMSO2   ','DMSO    ','DMSP    ','CH3SO2  ','HO2NO2  ',         &
     &  'CH3SO   ','CH3SO3  ','MSA     ','ORGNIT  ','NH3     ',         &
     &  'O3_STRAT','CO_EMIT ','H2O     ','Be-7    ','Be-10   ',         &
     &  'Rn-222  ','Pb-210  ','STRAT   ','TROP    '/)

! Indicies for chemical species
      INTEGER, PARAMETER ::                                             &
     &  i_od      = 1, i_op      = 2, i_oh      = 3, i_no      = 4,     &
     &  i_no2     = 5, i_no3     = 6, i_n2o5    = 7, i_co      = 8,     &
     &  i_ch4     = 9, i_hcho    =10, i_o3      =11, i_h2      =12,     &
     &  i_hno3    =13, i_h2o2    =14, i_ch3o2   =15, i_ho2     =16,     &
     &  i_c2h6    =17, i_c2h5o2  =18, i_ch3cho  =19, i_ch3coo2 =20,     &
     &  i_pan     =21, i_ch3ooh  =22, i_nc4h10  =23, i_sc4h9o2 =24,     &
     &  i_ch3coe  =25, i_so2     =26, i_c2h4    =27, i_c3h6    =28,     &
     &  i_oxyl    =29, i_sa      =30, i_c3h8    =31, i_c3h7o2  =32,     &
     &  i_c3h7ooh =33, i_c2h5ooh =34, i_c4h9ooh =35, i_ch3oh   =36,     &
     &  i_acetone =37, i_aceto2  =38, i_ammsul  =39, i_ch3cox  =40,     &
     &  i_ch2o2c  =41, i_mglyox  =42, i_ch3chx  =43, i_glyox   =44,     &
     &  i_oxyl1   =45, i_memald  =46, i_memald1 =47, i_c5h8    =48,     &
     &  i_ro2ip1  =49, i_mvk     =50, i_ro2ip2  =51, i_isopooh =52,     &
     &  i_mvkooh  =53, i_toluen  =54, i_rnc2h4  =55, i_rnc3h6  =56,     &
     &  i_rnc5h8  =57, i_naer    =58, i_tolp1   =59, i_dms     =60,     &
     &  i_dmso2   =61, i_dmso    =62, i_dmsp    =63, i_ch3so2  =64,     &
     &  i_ho2no2  =65, i_ch3so   =66, i_ch3so3  =67, i_msa     =68,     &
     &  i_orgnit  =69, i_nh3     =70, i_o3_strat=71, i_co_emit =72,     &
     &  i_h2o     =73, i_be7     =74, i_be10    =75, i_rn222   =76,     &
     &  i_pb210   =77, i_ph      =-1, i_strat   =78, i_trop    =79,     &
     &  i_lnox    =80, i_acnox   =81

! Numbers of the depositing species
      INTEGER, DIMENSION(n_spec_deposit) :: dep_spec = (/               &
     &  i_no2, i_o3, i_pan, i_so2, i_nh3,                               &
     &  i_co, i_ch4, i_h2,                                              &
     &  i_hno3, i_h2o2, i_ch3ooh, i_c2h5ooh, i_c3h7ooh, i_c4h9ooh,      &
     &  i_isopooh, i_mvkooh,                                            &
     &  i_sa, i_ammsul, i_naer, i_msa, i_orgnit /)

! Field Code indexing
      INTEGER :: ifc
      INTEGER, DIMENSION(nc), PARAMETER :: fcodes =                     &
     &(/ (1700 + ifc, ifc = 1, nc) /)

! Select emission scenario. Choices are:
!  bu - IIASA "Current Legislation" scenario
!  mf - IIASA "Maximum Feasible Reductions" scenario
!  A2 - Special version of A2 scenario for IPCC 4AR runs (NB upper case)
!  **  Note that if using the 3 scenarios above, l_emiss_current must be
!  **  set to FALSE, and emiss_year must be either 2000 or 2030 only.
!  fi - SRES Fossil Fuel Intensive (A1FI) scenario
!  a2 - SRES A2 scenario
!  ab - SRES A1B scenario
!  b1 - SRES B1 scenario
!  b2 - SRES B2 scenario
!  pi - Pre-industrial scenario (nominally 1850).
!
      CHARACTER (LEN=2), PARAMETER :: scenario = 'b2'
!
! Set L_Emiss_Current true to use current time in scenarios,
! else use emiss_year
      LOGICAL, PARAMETER :: l_emiss_current = .FALSE.
      INTEGER, PARAMETER :: emiss_year = 2000
!
! Determines whether Aircraft and Lightning emissions are used. If using
! scenario 'pi' set aircrafton to FALSE
!
      LOGICAL, PARAMETER :: lightningon = .TRUE.
      LOGICAL, PARAMETER :: aircrafton = .TRUE.
      REAL,    PARAMETER :: lightmult = 7.0 / 1.5 ! Multiplier for LNOx
!
! Define arrays for annual emission totals:
! (the total in the CLASS array will determine the annual emission,
!  totals are read in from the monthly files to allow normalisation)
      REAL, SAVE :: annual_dms,                                         &
     &              annual_isop,                                        &
     &              annual_land_nvoc,                                   &
     &              annual_ocean_nvoc,                                  &
     &              annual_soil_nox
!
! masses in kg/mol
      REAL, PARAMETER ::                                                &
     &  mh2o    =0.018, mo3  =0.048, mnit   =0.014, mair   =0.02897,    &
     &  mno2    =0.046, msul =0.032, mco    =0.028, mch4   =0.016,      &
     &  mhcho   =0.030, mc2h6=0.030, mch3cho=0.044, mnc4h10=0.058,      &
     &  mc2h4   =0.028, mc3h6=0.042, mc3h8  =0.044, mch3oh =0.032,      &
     &  macetone=0.058, moxyl=0.106, mtoluen=0.092, mh2    =0.002,      &
     &  mc5h8   =0.068, mcarb=0.012, mbe7   =0.007, mbe10  =0.010,      &
     &  mrn222  =0.222

      CONTAINS

        SUBROUTINE IN_STOCHEM_CHM_GET_PATH

          IMPLICIT  NONE
          INTEGER :: getenv

!     Defult Directory names
          CHARACTER(LEN=*), PARAMETER ::                                &
     &      datdir_default = '/data/cr/cce/hadstoch/data/',             &
     &      emdir_default  = '/data/cr/cce/hadstoch/surf_emiss/',       &
     &      emdir2_default = '/data/cr/cce/hadstoch/other_emiss/'

          IF (getenv('STOCHEM_DATDIR',DATDIR) == 0) THEN
            datdir = datdir_default           ! In case is not found
          END IF
          IF (getenv('STOCHEM_EMDIR',EMDIR) == 0) THEN
            emdir  = emdir_default            ! In case is not found
          END IF
          IF (getenv('STOCHEM_EMDIR2',EMDIR2) == 0) THEN
            emdir2 = emdir2_default           ! In case is not found
          END IF

          WRITE(6,*) 'STOCHEM_DATDIR : ',TRIM(datdir)
          WRITE(6,*) 'STOCHEM_EMDIR  : ',TRIM(emdir)
          WRITE(6,*) 'STOCHEM_EMDIR2 : ',TRIM(emdir2)

          RETURN

        END SUBROUTINE IN_STOCHEM_CHM_GET_PATH

      END MODULE IN_STOCHEM_CHM
#endif
