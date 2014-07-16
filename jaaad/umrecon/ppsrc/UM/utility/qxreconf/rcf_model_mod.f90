
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Data module defining STSHCOMP namelist

Module Rcf_Model_Mod

! Description:
!   Data module to define the STSHCOMP namelist and related addressing
!   variables.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Use Rcf_Submodel_Mod, Only :      &
    N_Submodel_partition_max, &
    N_Internal_Model_Max

Use Rcf_Ppx_Info_Mod, Only :  &
    NsectP,                   &
    NitemP

Implicit None

! These 2 parameters used to be in rcf_version_mod, but that has been
!  scrapped
Integer, Parameter  ::   Outfile_S  = 20
Integer, Parameter  ::   Outfile_E  = 149

Integer, Parameter  ::   A_Max_UKCAVars=150 ! Max.no.of tracers allowed
Integer, Parameter  ::   A_Max_TrVars = 150 ! Max.no.of tracers allowed

Integer, Parameter  ::   AASSETS    = 9
Integer, Parameter  ::   MEAD_TYPES = 4

Real, Save          ::   H_A_EWSPACE ,H_A_NSSPACE
Real, Save          ::   H_A_FIRSTLAT,H_A_FIRSTLONG
Real, Save          ::   H_A_POLELAT ,H_A_POLELONG

Integer, Save       ::   H_A_GROUP
Integer, Save       ::   H_OROG_ROUGH
Integer, Save       ::   A_ASSMGRPS
Integer, Save       ::   NUM_PVPR

Logical, Save       ::   A_RECON
Logical, Save       ::   H_OROG_GRAD
Logical, Save       ::   ATMODS
Logical, Save       ::   CMODS
Logical, Save       ::   LMESO

Logical, Save       ::   TRACER_A    (0:A_Max_TrVars)
Logical, Save       ::   TR_UKCA_A    (0:A_Max_UKCAVars)
Logical, Save       ::   AASSET   (AASSETS)
Integer, Parameter  ::   MAX_AOBS  = 100
Integer, Save       ::   AOBINC   (MAX_AOBS)
Integer, Save       ::   AOBGRP   (MAX_AOBS)
Integer, Save       ::   AASPF    (AASSETS)
Integer, Save       ::   AASPL    (AASSETS)
Integer, Save       ::   RUN_TARGET_END( 6)

!Total data length for primary fields for each submodel data partition
Integer, Save       ::   LPRIM(N_SUBMODEL_PARTITION_MAX)
!Total data length for primary fields for each internal model
Integer, Save       ::   LPrimIM(N_INTERNAL_MODEL_MAX)
!Total data length for diagnostic flds for each submodel data partition
Integer, Save       ::   LDUMP(N_SUBMODEL_PARTITION_MAX)
!Total data length for diagnostic flds for each internal model
Integer, Save       ::   LDumpIM(N_INTERNAL_MODEL_MAX)
!Total data length for secondary flds for each submodel data partition
Integer, Save       ::   LSECD(N_SUBMODEL_PARTITION_MAX)
!Total data length for secondary flds for each internal model
Integer, Save       ::   LSecdIM(N_INTERNAL_MODEL_MAX)
!Total workspace length for each submodel data partition
Integer, Save       ::   LWORK(N_SUBMODEL_PARTITION_MAX)
!Total number of headers (i.e. levels) for each submodel data partition
Integer, Save       ::   NHeadSub(N_SUBMODEL_PARTITION_MAX)
!Total number of headers (i.e. levels) for each internal model
Integer, Save       ::   NHEAD(N_INTERNAL_MODEL_MAX)
!Total length of extra space for each submod. data part.
Integer, Save       ::   LEXTRA(N_SUBMODEL_PARTITION_MAX)
!Data length for dual-time level ocean fields
Integer, Save       ::   LPRIM_O2
Integer, Save       ::   ITEM_MAX_REQ
Integer, Save       ::   ITEM_MAX_ALL

Integer, Save       ::   NRECS_S
Integer, Save       ::   NTIMES_S
Integer, Save       ::   NSERBLK_S
Integer, Save       ::   NSERREC_S
Integer, Save       ::   NLEVL_S
Integer, Save       ::   NMAXLEV_S
Integer, Save       ::   NPSLISTS_S
Integer, Save       ::   NMAXPSL_S
Integer, Save       ::   NHEAD_FILE(OUTFILE_S:OUTFILE_E)
Logical, Save       ::   LSTUSER
Integer, Save       ::   ReconItems(NITEMP)

Character (Len = 1), Save  ::  H_ATMOS
Character (Len = 1), Save  ::  H_OCEAN
Character (Len = 1), Save  ::  H_SLAB
Character (Len = 1), Save  ::  H_WAVE
Character (Len = 1), Save  ::  H_FLOOR
Character (Len = 1), Save  ::  H_STRAT
Character (Len = 1), Save  ::  H_SLAB_CAL
Character (Len = 1), Save  ::  H_TOTEM
Character (Len = 1), Save  ::  H_GLOBAL(N_INTERNAL_MODEL_MAX)
Integer, Save              ::  H_VERS  &
                               (N_INTERNAL_MODEL_MAX,0:NSECTP)


! These are set in SETMODL:
Integer, Save       ::   MEAN_NUMBER(N_INTERNAL_MODEL_MAX)

Real, Save          ::   H_W_EWSPACE ,H_W_NSSPACE
Real, Save          ::   H_W_FIRSTLAT,H_W_FIRSTLONG

! Variables read in by namelist and used in SETMODL
Integer, Save       ::   OCAAA
Integer, Save       ::   NWTRAIN
Real, Save          ::   EWSPACEA, NSSPACEA
Real, Save          ::   FRSTLATA, FRSTLONA

Logical, Save       ::   ZonAvOzone
Logical, Save       ::   ZONAVTPPSOZONE !! control for tpps_ozone
Integer, Save       ::   IVDF
Real, Save          ::   LATS
Real, Save          ::   LONS
Integer, Save       ::   LWBND
Integer, Save       ::   LWINC
Integer, Save       ::   NECF(50)
Integer, Save       ::   OASFLDID(7)
Integer, Save       ::   OASLEV(6) ! dimensioned by max no of
                                   ! O-Assm groups
Integer, Save       ::   OBAS
Integer, Save       ::   OBS
Integer, Save       ::   OCALB
Integer, Save       ::   OCBOHaney
Integer, Save       ::   OICE
Integer, Save       ::   OIDYN
Integer, Save       ::   OMP(4)
Real, Save          ::   POLELATA
Real, Save          ::   POLELONA
Integer, Save       ::   PSA
Integer, Save       ::   StLevGWdrag
Integer, Save       ::   SWBND
Integer, Save       ::   SWINC
Integer, Save       ::   TCA(A_Max_TrVars)
Integer, Save       ::   TC_UKCA(A_Max_UKCAVars)
Integer, Save       ::   TCO(29)
Integer, Save       ::   BotVDiffLev
Integer, Save       ::   TopVDiffLev


Character (Len = 2), Save  ::  ATMOS_SR(0:NSECTP)
Character (Len = 2), Save  ::  OCEAN_SR(0:NSECTP)
Character (Len = 2), Save  ::  SLAB_SR (0:NSECTP)
Character (Len = 2), Save  ::  WAVE_SR (0:NSECTP)
Character (Len = 2), Save  ::  INDEP_SR(0:NSECTP)

Character (Len = 1), Save  ::  BSPMSL
Character (Len = 1), Save  ::  CCEW
Character (Len = 1), Save  ::  FLOOR
Character (Len = 1), Save  ::  IDO
Character (Len = 1), Save  ::  LOSSM
Character (Len = 1), Save  ::  MLMO
Character (Len = 1), Save  ::  OAFLD(7)
Character (Len = 1), Save  ::  OCARB
Character (Len = 1), Save  ::  OROGR
Character (Len = 1), Save  ::  OSFC
Character (Len = 1), Save  ::  SCAL
Character (Len = 1), Save  ::  SSTAnom
Character (Len = 1), Save  ::  SWMCR
Character (Len = 1), Save  ::  TOTAE
Character (Len = 1), Save  ::  TOTEM
Character (Len = 1), Save  ::  UPD175
Character (Len = 1), Save  ::  MESO


Namelist/STSHCOMP/                                               &
RUN_TARGET_END,                                                  &
INDEP_SR    ,ATMOS_SR    ,                                       &
OCAAA       ,EWSPACEA    ,POLELATA ,FRSTLATA,LATS   ,            &
             NSSPACEA    ,POLELONA ,FRSTLONA,LONS   ,            &
SWBND       ,LWBND       ,SWINC    ,LWINC   ,OROGR  ,            &
ZONAVOZONE  ,ZONAVTPPSOZONE, SWMCR       ,MESO     ,             &
StLevGWdrag ,BotVDiffLev, TopVDiffLev,                           &
OCALB       ,FLOOR       ,AOBINC   ,TOTAE   ,TOTEM  ,            &
TCA         ,TC_UKCA     ,                                       &
SSTAnom     ,SCAL        ,AOBGRP   ,                             &
NECF        ,BSPMSL      ,CCEW              ,UPD175 ,            &
OIDYN       ,OBAS        ,OCBOHaney,OBS     ,OICE   ,IVDF ,IDO,  &
OCARB       ,MLMO        ,PSA      ,OSFC    ,                    &
LOSSM       ,OASLEV      ,OAFLD    ,                             &
OMP         ,OASFLDID

End Module Rcf_Model_Mod
