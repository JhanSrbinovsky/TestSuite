#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ data module defining RECON namelist

Module Rcf_Recon_Mod

! Description:
!   Defines variables for the RECON namelist.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Use Rcf_Submodel_Mod, Only : &
    Submodel_Ident

Use Rcf_PrintStatus_Mod, Only : &
    L_IO_Timer

Implicit None

! Comdecks
#include "c_mdi.h"

Integer       :: DUMP_PACK
Integer       :: LAND_POINTS
Integer, Save :: NTILES
Integer, Save :: NICE
Integer       :: MODEL_LEVELS
Integer       :: WET_LEVELS
Integer, Save :: TR_VARS
Integer       :: TR_UKCA
Integer       :: TR_LEVELS
Integer       :: ST_LEVELS
Integer       :: SM_LEVELS
Integer       :: BL_LEVELS
Integer       :: OZONE_LEVELS
Integer       :: TPPS_OZONE_LEVELS
Integer       :: CONV_LEVELS
Integer       :: P_ROWS
Integer       :: ROW_LENGTH
Integer       :: RIVER_ROWS
Integer       :: RIVER_ROW_LENGTH
Integer       :: W_ZERO_START
Integer       :: W_ZERO_END
Integer, Save :: RIMWIDTHA
Integer, Save :: RIMWIDTHO

Integer, Save :: LEN_FIXHD     = IMDI
Integer, Save :: LEN_INTHD     = IMDI
Integer, Save :: LEN_REALHD    = IMDI
Integer, Save :: LEN2_LEVDEPC  = IMDI
Integer, Save :: LEN2_ROWDEPC  = IMDI
Integer, Save :: LEN2_COLDEPC  = IMDI
Integer, Save :: LEN2_FLDDEPC  = IMDI
Integer, Save :: LEN_EXTCNST   = IMDI
Integer, Save :: LEN_DUMPHIST  = IMDI

Real, Save    :: ANVIL_FACTOR
Real, Save    :: TOWER_FACTOR
Real, Save    :: Q_MIN
Real    :: PERTURBATION

Logical, Save :: GRIB
Logical, Save :: UARS
Logical, Save :: RESET
Logical, Save :: TRANS
Logical, Save :: LSPIRAL_S
Logical, Save :: LCAL360
Logical, Save :: LOZONE_ZONAL
Logical, Save :: L_TPPS_OZONE_ZONAL
Logical, Save :: LAMIPII
Logical, Save :: L_SSTANOM
Logical, Save :: L_Cloud_Deep
Logical, Save :: Var_Recon
Logical, Save :: USE_SMC_STRESS

! Namelist RECON.

Namelist /RECON/                                                  &
 SUBMODEL_IDENT, DUMP_PACK, LAND_POINTS, NTILES, ANVIL_FACTOR,    &
 TOWER_FACTOR, RIMWIDTHA, MODEL_LEVELS, WET_LEVELS, P_ROWS,       &
 ROW_LENGTH, TR_LEVELS, TR_VARS, ST_LEVELS,SM_LEVELS, BL_LEVELS,  &
 OZONE_LEVELS, CONV_LEVELS, LEN_FIXHD, LEN_INTHD, LEN_REALHD,     &
 LEN2_LEVDEPC, LEN2_ROWDEPC, LEN2_COLDEPC, LEN2_FLDDEPC,          &
 LEN_EXTCNST, LEN_DUMPHIST, GRIB, UARS, RESET, PERTURBATION,      &
 TRANS,  LSPIRAL_S, RIMWIDTHO, LCAL360, LOZONE_ZONAL, LAMIPII,    &
 L_SSTANOM, L_Cloud_Deep, Var_Recon, Q_Min, L_IO_Timer, NICE,     &
 W_Zero_Start, W_Zero_End, RIVER_ROWS, RIVER_ROW_LENGTH,          &
 TR_UKCA, Use_SMC_Stress

End Module Rcf_Recon_Mod
#endif
