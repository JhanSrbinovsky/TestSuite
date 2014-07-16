#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module to define setup and indices for gas and aerosol tracers.
!    Many of these are not used in UM but kept for consistency with
!    set-up in TOMCAT.
!    Contains public subroutines:
!      UKCA_ORGV1_SUSSBCOC_5MODE
!      UKCA_ORGV1_SUSSBCOCSO_5MODE
!      UKCA_SV1_SUSS_4MODE
!    which define particular setups.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Current Code Owner:       Graham Mann
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
! Subroutine Interface:
      MODULE UKCA_SETUP_INDICES
! ---------------------------------------------------------------------|
!+ Module to define setup and indices for gas and aerosol tracers.
!+ Many of these are not used in UM but kept for consistency with
!+ set-up in TOMCAT.
!
! History:
! Version    Date      Author         Description
! -----      --------  ---------      ------------
! vn6.1     05/VII/07  Graham Mann    Original code
!
! ---------------------------------------------------------------------|

      IMPLICIT NONE

      INTEGER :: NTRAER     ! # of aerosol advected tracers
      INTEGER :: NCHEMG     ! # of gas phase chemistry tracers
      INTEGER :: ICHEM      ! 1/0 = gas phase chemistry tracers on/off
      INTEGER :: NADVG      ! # of gas phase advected tracers
      INTEGER :: NOFFOX     ! # of offline oxidant species
      INTEGER :: NTRAG      ! total # of gas phase species
      INTEGER :: BUDGET     ! 1/0 = budget fields on/off
      INTEGER :: NBUDCHEM   ! # of gas chem budget fields
      INTEGER :: IDUSTDEP   ! 1/0 = size-resolved dust dep. fields on/off
      INTEGER :: NDUSTDEP   ! # of size-resolved dust deposition fields
      INTEGER :: NBUDAER    ! # of aerosol budget fields
      INTEGER :: NBUDAERTOT ! total # of aersol budget fields
      INTEGER :: NBUDGET    ! total # of budget terms
      INTEGER :: TRAQU      ! 1/0 = extra diagnostic terms on/off
      INTEGER :: NTRAQU     ! # of extra diagnostic terms

! .. (indices for advected gas phase tracers in array S0)
! .. (8 sulfur species, H2O2, MONOTER, SEC_ORG, Q3D, PT)
      INTEGER :: MSOTWO     ! for SO2
      INTEGER :: MMeSMe     ! for DMS
      INTEGER :: MH2SO4     ! for H2SO4
      INTEGER :: MDMSO      ! for DMSO
      INTEGER :: MMSA       ! for MSA
      INTEGER :: MCS2       ! for CS2
      INTEGER :: MH2S       ! for H2S
      INTEGER :: MCOS       ! for COS
      INTEGER :: MMONOTER   ! for lumped monoterpene species
      INTEGER :: MSEC_ORG   ! for involatile organic species
      INTEGER :: MH2O2      ! for H2O2
      INTEGER :: MQ3D       ! for water vapour
      INTEGER :: MPT        ! for potential temperature
!
! .. (indices for all gas phase tracers in array ST)
      INTEGER :: NSOTWO     ! for SO2
      INTEGER :: NMeSMe     ! for DMS
      INTEGER :: NH2SO4     ! for H2SO4
      INTEGER :: NDMSO      ! for DMSO
      INTEGER :: NMSA       ! for MSA
      INTEGER :: NCS2       ! for CS2
      INTEGER :: NH2S       ! for H2S
      INTEGER :: NCOS       ! for COS
      INTEGER :: NMONOTER   ! for lumped monoterpene species
      INTEGER :: NSEC_ORG   ! for involatile organic species
      INTEGER :: NH2O2      ! for H2O2
      INTEGER :: NO3        ! for O3
      INTEGER :: NOH        ! for OH
      INTEGER :: NNO3       ! for NO3
      INTEGER :: NQ3D       ! for water vapour
      INTEGER :: NPT        ! for potential temperature
!
! .. (indices for gas phase budget mass fluxes)
      INTEGER :: NDMSEMOC   ! for DMS emissions flux (oceanic sources)
      INTEGER :: NDMSTEND   ! for DMS ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NSO2EMAN   ! for SO2 emissions flux (anthropogenic sources)
      INTEGER :: NSO2EMBM   ! for SO2 emissions flux (biomass burning sources)
      INTEGER :: NSO2EMVL   ! for SO2 emissions flux (volcanic sources)
      INTEGER :: NSO2TEND   ! for SO2 ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NSO2DDEP   ! for SO2 dry deposition flux
      INTEGER :: NSO2WDEP   ! for SO2 wet deposition flux
      INTEGER :: NH2SO4TEND ! for H2SO4 ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NH2SO4DDEP ! for H2SO4 dry deposition flux
      INTEGER :: NCOSEMAN   ! for COS emissions flux (anthropogenic sources)
      INTEGER :: NCOSEMOC   ! for COS emissions flux (oceanic sources)
      INTEGER :: NCOSTEND   ! for COS ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NCS2EMAN   ! for CS2 emissions flux (anthropogenic sources)
      INTEGER :: NCS2EMOC   ! for CS2 emissions flux (oceanic sources)
      INTEGER :: NCS2TEND   ! for CS2 ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NDMSOTEND  ! for DMSO ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NDMSODDEP  ! for DMSO dry deposition flux
      INTEGER :: NMSATEND   ! for MSA ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NMSADDEP   ! for MSA dry deposition flux
      INTEGER :: NTERP_EM   ! for terpene emissions flux (biogenic sources)
      INTEGER :: NTERP_TEND ! for terpene ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NTERP_DDEP ! for terpene dry deposition flux
      INTEGER :: NSORG_TEND ! for sec_org ASAD tendency (combined chem,ddep,wdep)
      INTEGER :: NSORG_DDEP ! for sec_org dry deposition flux
      INTEGER :: NSORG_WDEP ! for sec_org wet deposition flux
!
! .. (indices for aerosol budget mass fluxes)
      INTEGER :: NMASPRIMSUAITSOL ! for primary SU ems (all srces) to Aitsol
      INTEGER :: NMASPRIMSUACCSOL ! for primary SU ems (all srces) to accsol
      INTEGER :: NMASPRIMSUCORSOL ! for primary SU ems (all srces) to corsol
      INTEGER :: NMASPRIMSSACCSOL ! for primary SS ems (ocean) to accsol
      INTEGER :: NMASPRIMSSCORSOL ! for primary SS ems (ocean) to corsol
      INTEGER :: NMASPRIMBCAITSOL ! for primary BC ems (all srces) to Aitsol
      INTEGER :: NMASPRIMBCAITINS ! for primary BC ems (all srces) to Aitins
      INTEGER :: NMASPRIMOCAITSOL ! for primary OC ems (all srces) to Aitsol
      INTEGER :: NMASPRIMOCAITINS ! for primary OC ems (all srces) to Aitins
      INTEGER :: NMASDDEPSUNUCSOL ! for drydep of SU from nucsol
      INTEGER :: NMASDDEPSUAITSOL ! for drydep of SU from Aitsol
      INTEGER :: NMASDDEPSUACCSOL ! for drydep of SU from accsol
      INTEGER :: NMASDDEPSUCORSOL ! for drydep of SU from corsol
      INTEGER :: NMASDDEPSSACCSOL ! for drydep of SS from accsol
      INTEGER :: NMASDDEPSSCORSOL ! for drydep of SS from corsol
      INTEGER :: NMASDDEPBCAITSOL ! for drydep of BC from Aitsol
      INTEGER :: NMASDDEPBCACCSOL ! for drydep of BC from accsol
      INTEGER :: NMASDDEPBCCORSOL ! for drydep of BC from corsol
      INTEGER :: NMASDDEPBCAITINS ! for drydep of BC from Aitins
      INTEGER :: NMASDDEPOCNUCSOL ! for drydep of OC from nucsol
      INTEGER :: NMASDDEPOCAITSOL ! for drydep of OC from Aitsol
      INTEGER :: NMASDDEPOCACCSOL ! for drydep of OC from accsol
      INTEGER :: NMASDDEPOCCORSOL ! for drydep of OC from corsol
      INTEGER :: NMASDDEPOCAITINS ! for drydep of OC from Aitins
      INTEGER :: NMASDDEPSONUCSOL ! for drydep of SO from nucsol
      INTEGER :: NMASDDEPSOAITSOL ! for drydep of SO from Aitsol
      INTEGER :: NMASDDEPSOACCSOL ! for drydep of SO from accsol
      INTEGER :: NMASDDEPSOCORSOL ! for drydep of SO from corsol
      INTEGER :: NMASNUSCSUNUCSOL ! for nuscav of SU from nucsol
      INTEGER :: NMASNUSCSUAITSOL ! for nuscav of SU from Aitsol
      INTEGER :: NMASNUSCSUACCSOL ! for nuscav of SU from accsol
      INTEGER :: NMASNUSCSUCORSOL ! for nuscav of SU from corsol
      INTEGER :: NMASNUSCSSACCSOL ! for nuscav of SS from accsol
      INTEGER :: NMASNUSCSSCORSOL ! for nuscav of SS from corsol
      INTEGER :: NMASNUSCBCAITSOL ! for nuscav of BC from Aitsol
      INTEGER :: NMASNUSCBCACCSOL ! for nuscav of BC from accsol
      INTEGER :: NMASNUSCBCCORSOL ! for nuscav of BC from corsol
      INTEGER :: NMASNUSCBCAITINS ! for nuscav of BC from Aitins
      INTEGER :: NMASNUSCOCNUCSOL ! for nuscav of OC from nucsol
      INTEGER :: NMASNUSCOCAITSOL ! for nuscav of OC from Aitsol
      INTEGER :: NMASNUSCOCACCSOL ! for nuscav of OC from accsol
      INTEGER :: NMASNUSCOCCORSOL ! for nuscav of OC from corsol
      INTEGER :: NMASNUSCOCAITINS ! for nuscav of OC from Aitins
      INTEGER :: NMASNUSCSONUCSOL ! for nuscav of SO from nucsol
      INTEGER :: NMASNUSCSOAITSOL ! for nuscav of SO from Aitsol
      INTEGER :: NMASNUSCSOACCSOL ! for nuscav of SO from accsol
      INTEGER :: NMASNUSCSOCORSOL ! for nuscav of SO from corsol
      INTEGER :: NMASIMSCSUNUCSOL ! for imscav of SU from nucsol
      INTEGER :: NMASIMSCSUAITSOL ! for imscav of SU from Aitsol
      INTEGER :: NMASIMSCSUACCSOL ! for imscav of SU from accsol
      INTEGER :: NMASIMSCSUCORSOL ! for imscav of SU from corsol
      INTEGER :: NMASIMSCSSACCSOL ! for imscav of SS from accsol
      INTEGER :: NMASIMSCSSCORSOL ! for imscav of SS from corsol
      INTEGER :: NMASIMSCBCAITSOL ! for imscav of BC from Aitsol
      INTEGER :: NMASIMSCBCACCSOL ! for imscav of BC from accsol
      INTEGER :: NMASIMSCBCCORSOL ! for imscav of BC from corsol
      INTEGER :: NMASIMSCBCAITINS ! for imscav of BC from Aitins
      INTEGER :: NMASIMSCOCNUCSOL ! for imscav of OC from nucsol
      INTEGER :: NMASIMSCOCAITSOL ! for imscav of OC from Aitsol
      INTEGER :: NMASIMSCOCACCSOL ! for imscav of OC from accsol
      INTEGER :: NMASIMSCOCCORSOL ! for imscav of OC from corsol
      INTEGER :: NMASIMSCOCAITINS ! for imscav of OC from Aitins
      INTEGER :: NMASIMSCSONUCSOL ! for imscav of SO from nucsol
      INTEGER :: NMASIMSCSOAITSOL ! for imscav of SO from Aitsol
      INTEGER :: NMASIMSCSOACCSOL ! for imscav of SO from accsol
      INTEGER :: NMASIMSCSOCORSOL ! for imscav of SO from corsol
      INTEGER :: NMASCLPRSUAITSOL1 ! for in-cl. ox. of SO2 by H2O2 to SU in Aitsol
      INTEGER :: NMASCLPRSUACCSOL1 ! for in-cl. ox. of SO2 by H2O2 to SU in accsol
      INTEGER :: NMASCLPRSUCORSOL1 ! for in-cl. ox. of SO2 by H2O2 to SU in corsol
      INTEGER :: NMASCLPRSUAITSOL2 ! for in-cl. ox. of SO2 by O3   to SU in Aitsol
      INTEGER :: NMASCLPRSUACCSOL2 ! for in-cl. ox. of SO2 by O3   to SU in accsol
      INTEGER :: NMASCLPRSUCORSOL2 ! for in-cl. ox. of SO2 by O3   to SU in corsol
      INTEGER :: NMASCONDSUNUCSOL ! for conden. of SU to nucsol
      INTEGER :: NMASCONDSUAITSOL ! for conden. of SU to Aitsol
      INTEGER :: NMASCONDSUACCSOL ! for conden. of SU to accsol
      INTEGER :: NMASCONDSUCORSOL ! for conden. of SU to corsol
      INTEGER :: NMASCONDSUAITINS ! for conden. of SU to Aitins
      INTEGER :: NMASNUCLSUNUCSOL ! for nucln. of SU to nucsol
      INTEGER :: NMASCONDOCNUCSOL ! for conden. of OC to nucsol
      INTEGER :: NMASCONDOCAITSOL ! for conden. of OC to Aitsol
      INTEGER :: NMASCONDOCACCSOL ! for conden. of OC to accsol
      INTEGER :: NMASCONDOCCORSOL ! for conden. of OC to corsol
      INTEGER :: NMASCONDOCAITINS ! for conden. of OC to Aitins
      INTEGER :: NMASCONDSONUCSOL ! for conden. of SO to nucsol
      INTEGER :: NMASCONDSOAITSOL ! for conden. of SO to Aitsol
      INTEGER :: NMASCONDSOACCSOL ! for conden. of SO to accsol
      INTEGER :: NMASCONDSOCORSOL ! for conden. of SO to corsol
      INTEGER :: NMASCONDSOAITINS ! for conden. of SO to Aitins
      INTEGER :: NMASCOAGSUINTR12 ! for inter-modal coag SU nucsol->Aitsol
      INTEGER :: NMASCOAGSUINTR13 ! for inter-modal coag SU nucsol->accsol
      INTEGER :: NMASCOAGSUINTR14 ! for inter-modal coag SU nucsol->corsol
      INTEGER :: NMASCOAGSUINTR15 ! for inter-modal coag SU nucsol->Aitins
      INTEGER :: NMASCOAGOCINTR12 ! for inter-modal coag OC nucsol->Aitsol
      INTEGER :: NMASCOAGOCINTR13 ! for inter-modal coag OC nucsol->accsol
      INTEGER :: NMASCOAGOCINTR14 ! for inter-modal coag OC nucsol->corsol
      INTEGER :: NMASCOAGOCINTR15 ! for inter-modal coag OC nucsol->Aitins
      INTEGER :: NMASCOAGSOINTR12 ! for inter-modal coag SO nucsol->Aitsol
      INTEGER :: NMASCOAGSOINTR13 ! for inter-modal coag SO nucsol->accsol
      INTEGER :: NMASCOAGSOINTR14 ! for inter-modal coag SO nucsol->corsol
      INTEGER :: NMASCOAGSOINTR15 ! for inter-modal coag SO nucsol->Aitins
      INTEGER :: NMASCOAGSUINTR23 ! for inter-modal coag SU Aitsol->accsol
      INTEGER :: NMASCOAGBCINTR23 ! for inter-modal coag BC Aitsol->accsol
      INTEGER :: NMASCOAGOCINTR23 ! for inter-modal coag OC Aitsol->accsol
      INTEGER :: NMASCOAGSOINTR23 ! for inter-modal coag SO Aitsol->accsol
      INTEGER :: NMASCOAGSUINTR24 ! for inter-modal coag SU Aitsol->corsol
      INTEGER :: NMASCOAGBCINTR24 ! for inter-modal coag BC Aitsol->corsol
      INTEGER :: NMASCOAGOCINTR24 ! for inter-modal coag OC Aitsol->corsol
      INTEGER :: NMASCOAGSOINTR24 ! for inter-modal coag SO Aitsol->corsol
      INTEGER :: NMASCOAGSUINTR34 ! for inter-modal coag SU accsol->corsol
      INTEGER :: NMASCOAGBCINTR34 ! for inter-modal coag BC accsol->corsol
      INTEGER :: NMASCOAGOCINTR34 ! for inter-modal coag OC accsol->corsol
      INTEGER :: NMASCOAGSSINTR34 ! for inter-modal coag SS accsol->corsol
      INTEGER :: NMASCOAGSOINTR34 ! for inter-modal coag SO accsol->corsol
      INTEGER :: NMASCOAGBCINTR53 ! for inter-modal coag BC Aitins->accsol
      INTEGER :: NMASCOAGOCINTR53 ! for inter-modal coag OC Aitins->accsol
      INTEGER :: NMASCOAGBCINTR54 ! for inter-modal coag BC Aitins->corsol
      INTEGER :: NMASCOAGOCINTR54 ! for inter-modal coag OC Aitins->corsol
      INTEGER :: NMASAGEDSUINTR52 ! for SU ageing flux Aitins->Aitsol
      INTEGER :: NMASAGEDBCINTR52 ! for BC ageing flux Aitins->Aitsol
      INTEGER :: NMASAGEDOCINTR52 ! for OC ageing flux Aitins->Aitsol
      INTEGER :: NMASAGEDSOINTR52 ! for SO ageing flux Aitins->Aitsol
      INTEGER :: NMASMERGSUINTR12 ! for SU mode-merging flux nucsol->Aitsol
      INTEGER :: NMASMERGOCINTR12 ! for OC mode-merging flux nucsol->Aitsol
      INTEGER :: NMASMERGSOINTR12 ! for SO mode-merging flux nucsol->Aitsol
      INTEGER :: NMASMERGSUINTR23 ! for SU mode-merging flux Aitsol->accsol
      INTEGER :: NMASMERGBCINTR23 ! for BC mode-merging flux Aitsol->accsol
      INTEGER :: NMASMERGOCINTR23 ! for OC mode-merging flux Aitsol->accsol
      INTEGER :: NMASMERGSOINTR23 ! for SO mode-merging flux Aitsol->accsol
      INTEGER :: NMASMERGSUINTR34 ! for SU mode-merging flux accsol->corsol
      INTEGER :: NMASMERGSSINTR34 ! for SS mode-merging flux accsol->corsol
      INTEGER :: NMASMERGBCINTR34 ! for BC mode-merging flux accsol->corsol
      INTEGER :: NMASMERGOCINTR34 ! for OC mode-merging flux accsol->corsol
      INTEGER :: NMASMERGSOINTR34 ! for SO mode-merging flux accsol->corsol
      INTEGER :: NMASPROCSUINTR23 ! for SU cloud-processing Aitsol->accsol
      INTEGER :: NMASPROCBCINTR23 ! for BC cloud-processing Aitsol->accsol
      INTEGER :: NMASPROCOCINTR23 ! for OC cloud-processing Aitsol->accsol
      INTEGER :: NMASPROCSOINTR23 ! for SO cloud-processing Aitsol->accsol
!
! .. (indices for additional diagnostics)
      INTEGER :: NWTCNT1
      INTEGER :: NWTCNT2
      INTEGER :: NWTCNT3
      INTEGER :: NWTCNT4
      INTEGER :: NLCFRAC
      INTEGER :: NLWC
      INTEGER :: NRAINLS
      INTEGER :: NRAINCS
      INTEGER :: NWINDSP
!
      INTEGER, PARAMETER :: NCHEMGMAX=11 ! max value for NCHEMG
!
! Indices of aerosol components into which condensable gases go to
      INTEGER :: CONDENSABLE_CHOICE(NCHEMGMAX)
!
! Gas phase species which are condensable (T/F)
      LOGICAL :: CONDENSABLE(NCHEMGMAX)
!
! Molecular masses of gas phase species (kg/mol)
      REAL :: MM_GAS(NCHEMGMAX)
!
! Molecular diameter of condensable gas phase species (others = 0)
      REAL :: DIMEN(NCHEMGMAX)
!
      CONTAINS

      SUBROUTINE UKCA_ORGV1_SUSSBCOC_5MODE

      IMPLICIT NONE

!-----------------------------------------------------------
!
! MAIN ARRAY LENGTHS AND SWITCHES
!
      NTRAER=20          ! # of aerosol advected tracers
      NCHEMG=11          ! # of gas phase chemistry tracers
      ICHEM=1            ! 1/0 = gas phase chemistry tracers on/off
      BUDGET=1           ! 1/0 = budget fields on/off
      TRAQU=1            ! 1/0 = extra diagnostic terms on/off
      IDUSTDEP=0         ! 1/0 = size-resolved dustdep fields on/off
      NOFFOX=3*ICHEM     ! # of offline oxidant species
      NDUSTDEP=8         ! # of size-resolved dustdep fields
      NBUDCHEM=26*ICHEM  ! # of gas chem budget fields
      NBUDAER=107        ! # of aerosol budget fields

      NADVG=2+NCHEMG*ICHEM ! # of gas phase advected tracers
      NTRAG=NADVG+NOFFOX   ! total # of gas phase species

      NBUDAERTOT=NBUDAER+NDUSTDEP*IDUSTDEP ! # of dustdep fields
      NBUDGET=NBUDCHEM+NBUDAERTOT ! total # of budget terms
      NTRAQU=TRAQU*9 ! # of extra diagnostic terms
!
! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
!                   NTRAER+NTRAG+NBUDGET+NTRAQU             = NVTOT
! .. ICHEM=1 SUSSBCOC  20    16    133     9 = 36 + 133 + 9 = 178 (orgv1)
!
!                   NTRAER+NADVG  NTRA
! .. ICHEM=1 SUSSBCOC 20    13  =  33
!
!
!-----------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
!
      IF(ICHEM == 0) THEN
! .. below are 2 fields if chemistry switched off
       MSOTWO  = 0 ! not included
       MMeSMe  = 0 ! not included
       MH2SO4  = 0 ! not included
       MDMSO   = 0 ! not included
       MMSA    = 0 ! not included
       MCS2    = 0 ! not included
       MH2S    = 0 ! not included
       MCOS    = 0 ! not included
       MMONOTER= 0 ! not included
       MSEC_ORG= 0 ! not included
       MH2O2   = 0 ! not included
       MQ3D    = 1
       MPT     = 2
      ENDIF
!
      IF(ICHEM == 1) THEN
! .. below are 13 fields for orgv1.d
       MSOTWO  = 1
       MMeSMe  = 2
       MH2SO4  = 3
       MDMSO   = 4
       MMSA    = 5
       MCS2    = 6
       MH2S    = 7
       MCOS    = 8
       MMONOTER= 9
       MSEC_ORG=10
       MH2O2   =11
       MQ3D    =12
       MPT     =13
      ENDIF
!
! .. molecular weights (kg/mol) for gases for orgv1
      MM_GAS=(/0.064,0.062,0.098,0.078,0.096,0.076,0.034,0.060,         &
               0.136,0.150,0.034/) ! 8 S species, terp, sec_org, H2O2
!
      CONDENSABLE_CHOICE=(/0,0,1,0,0,0,0,0,0,3,0/)
! .. H2SO4   to condense into 1st aerosol component (CP_SU)
! .. SEC_ORG to condense into 3rd aerosol component (CP_OC)

      CONDENSABLE=(CONDENSABLE_CHOICE > 0)
!
      DIMEN=(/0.0,0.0,4.5e-10,0.0,0.0,0.0,0.0,0.0,0.0,4.5e-10,0.0/)
! Molecular diameters of condensable gas phase species (m)
!
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY
!
      IF(ICHEM == 0) THEN
! .. below are 2 fields if chemistry switched off
       NSOTWO  = 0 ! not included
       NMeSMe  = 0 ! not included
       NH2SO4  = 0 ! not included
       NDMSO   = 0 ! not included
       NMSA    = 0 ! not included
       NCS2    = 0 ! not included
       NH2S    = 0 ! not included
       NCOS    = 0 ! not included
       NMONOTER= 0 ! not included
       NSEC_ORG= 0 ! not included
       NH2O2   = 0 ! not included
       NO3     = 0 ! not included
       NOH     = 0 ! not included
       NNO3    = 0 ! not included
       NQ3D    = 1
       NPT     = 2
      ENDIF
!
      IF(ICHEM == 1) THEN
! .. below are 16 fields for orgv1.d
       NSOTWO  = 1
       NMeSMe  = 2
       NH2SO4  = 3
       NDMSO   = 4
       NMSA    = 5
       NCS2    = 6
       NH2S    = 7
       NCOS    = 8
       NMONOTER= 9
       NSEC_ORG=10
       NH2O2   =11
       NO3     =12
       NOH     =13
       NNO3    =14
       NQ3D    =15
       NPT     =16
      ENDIF
!
!---------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 26 gas phase budget quantity indices for orgv1
      NDMSEMOC  = 1
      NDMSTEND  = 2
      NSO2EMAN  = 3
      NSO2EMBM  = 4
      NSO2EMVL  = 5
      NSO2TEND  = 6
      NSO2DDEP  = 7
      NSO2WDEP  = 8
      NH2SO4TEND= 9
      NH2SO4DDEP=10
      NCOSEMAN  =11
      NCOSEMOC  =12
      NCOSTEND  =13
      NCS2EMAN  =14
      NCS2EMOC  =15
      NCS2TEND  =16
      NDMSOTEND =17
      NDMSODDEP =18
      NMSATEND  =19
      NMSADDEP  =20
      NTERP_EM  =21
      NTERP_TEND=22
      NTERP_DDEP=23
      NSORG_TEND=24
      NSORG_DDEP=25
      NSORG_WDEP=26
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! .. below are 107 aerosol budget quantity indices for SUSSBCOC [SO in OC]
!
! .. redo these budgets to follow Dominicks approach as categorized under
!
! .. BUD_PRIM, BUD_DDEP, BUD_NUSC, BUD_IMSC which have NBOX,NMODES,NCP
! ..
! .. then also have:
! ..
! .. BUD_CLPR --- NBOX,NMODES,2 (production of sulfate by H2O2, O3)
! .. BUD_COND --- NBOX,NMODES,2 (condensation onto all modes by H2SO4,SEC_ORG)
! .. BUD_NUCL --- NBOX,2 (BHN and BLN)
! .. BUD_COAG --- NBOX,NMODES,NMODES,NCP
! .. BUD_AGED --- NBOX,3,NCP,2 (ageing of 3 insoluble modes by H2SO4, SEC_ORG)
! .. BUD_MERG --- NBOX,NMODES,NCP
! .. BUD_PROC --- NBOX,NCP (processing of Aitken mode to accumulation mode)
! ..
!
      NMASPRIMSUAITSOL= 1
      NMASPRIMSUACCSOL= 2
      NMASPRIMSUCORSOL= 3
      NMASPRIMSSACCSOL= 4
      NMASPRIMSSCORSOL= 5
      NMASPRIMBCAITSOL= 0 ! BC only emitted to insoluble
      NMASPRIMBCAITINS= 6
      NMASPRIMOCAITSOL= 7
      NMASPRIMOCAITINS= 8
!
      NMASDDEPSUNUCSOL= 9
      NMASDDEPSUAITSOL=10
      NMASDDEPSUACCSOL=11
      NMASDDEPSUCORSOL=12
      NMASDDEPSSACCSOL=13
      NMASDDEPSSCORSOL=14
      NMASDDEPBCAITSOL=15
      NMASDDEPBCACCSOL=16
      NMASDDEPBCCORSOL=17
      NMASDDEPBCAITINS=18
      NMASDDEPOCNUCSOL=19
      NMASDDEPOCAITSOL=20
      NMASDDEPOCACCSOL=21
      NMASDDEPOCCORSOL=22
      NMASDDEPOCAITINS=23
      NMASDDEPSONUCSOL= 0 ! SO stored in OC cpt
      NMASDDEPSOAITSOL= 0 ! SO stored in OC cpt
      NMASDDEPSOACCSOL= 0 ! SO stored in OC cpt
      NMASDDEPSOCORSOL= 0 ! SO stored in OC cpt
!
      NMASNUSCSUNUCSOL=24
      NMASNUSCSUAITSOL=25
      NMASNUSCSUACCSOL=26
      NMASNUSCSUCORSOL=27
      NMASNUSCSSACCSOL=28
      NMASNUSCSSCORSOL=29
      NMASNUSCBCAITSOL=30
      NMASNUSCBCACCSOL=31
      NMASNUSCBCCORSOL=32
      NMASNUSCBCAITINS=33
      NMASNUSCOCNUCSOL=34
      NMASNUSCOCAITSOL=35
      NMASNUSCOCACCSOL=36
      NMASNUSCOCCORSOL=37
      NMASNUSCOCAITINS=38
      NMASNUSCSONUCSOL= 0 ! SO stored in OC cpt
      NMASNUSCSOAITSOL= 0 ! SO stored in OC cpt
      NMASNUSCSOACCSOL= 0 ! SO stored in OC cpt
      NMASNUSCSOCORSOL= 0 ! SO stored in OC cpt
!
      NMASIMSCSUNUCSOL=39
      NMASIMSCSUAITSOL=40
      NMASIMSCSUACCSOL=41
      NMASIMSCSUCORSOL=42
      NMASIMSCSSACCSOL=43
      NMASIMSCSSCORSOL=44
      NMASIMSCBCAITSOL=45
      NMASIMSCBCACCSOL=46
      NMASIMSCBCCORSOL=47
      NMASIMSCBCAITINS=48
      NMASIMSCOCNUCSOL=49
      NMASIMSCOCAITSOL=50
      NMASIMSCOCACCSOL=51
      NMASIMSCOCCORSOL=52
      NMASIMSCOCAITINS=53
      NMASIMSCSONUCSOL= 0 ! SO stored in OC cpt
      NMASIMSCSOAITSOL= 0 ! SO stored in OC cpt
      NMASIMSCSOACCSOL= 0 ! SO stored in OC cpt
      NMASIMSCSOCORSOL= 0 ! SO stored in OC cpt
!
      NMASCLPRSUAITSOL1=54
      NMASCLPRSUACCSOL1=55
      NMASCLPRSUCORSOL1=56
      NMASCLPRSUAITSOL2=57
      NMASCLPRSUACCSOL2=58
      NMASCLPRSUCORSOL2=59
!
      NMASCONDSUNUCSOL=60
      NMASCONDSUAITSOL=61
      NMASCONDSUACCSOL=62
      NMASCONDSUCORSOL=63
      NMASCONDSUAITINS=64
      NMASNUCLSUNUCSOL=65
      NMASCONDOCNUCSOL=66
      NMASCONDOCAITSOL=67
      NMASCONDOCACCSOL=68
      NMASCONDOCCORSOL=69
      NMASCONDOCAITINS=70
      NMASCONDSONUCSOL= 0 ! SO stored in OC cpt
      NMASCONDSOAITSOL= 0 ! SO stored in OC cpt
      NMASCONDSOACCSOL= 0 ! SO stored in OC cpt
      NMASCONDSOCORSOL= 0 ! SO stored in OC cpt
      NMASCONDSOAITINS= 0 ! SO stored in OC cpt
!
      NMASCOAGSUINTR12=71
      NMASCOAGSUINTR13=72
      NMASCOAGSUINTR14=73
      NMASCOAGSUINTR15=74
      NMASCOAGOCINTR12=75
      NMASCOAGOCINTR13=76
      NMASCOAGOCINTR14=77
      NMASCOAGOCINTR15=78
      NMASCOAGSOINTR12= 0 ! stored in NMASCOAGOCINTR12
      NMASCOAGSOINTR13= 0 ! stored in NMASCOAGOCINTR13
      NMASCOAGSOINTR14= 0 ! stored in NMASCOAGOCINTR14
      NMASCOAGSOINTR15= 0 ! stored in NMASCOAGOCINTR15
      NMASCOAGSUINTR23=79
      NMASCOAGBCINTR23=80
      NMASCOAGOCINTR23=81
      NMASCOAGSOINTR23= 0 ! stored in NMASCOAGOCINTR23
      NMASCOAGSUINTR24=82
      NMASCOAGBCINTR24=83
      NMASCOAGOCINTR24=84
      NMASCOAGSOINTR24= 0 ! stored in NMASCOAGOCINTR24
      NMASCOAGSUINTR34=85
      NMASCOAGBCINTR34=86
      NMASCOAGOCINTR34=87
      NMASCOAGSSINTR34=88
      NMASCOAGSOINTR34= 0 ! stored in NMASCOAGOCINTR34
!
      NMASCOAGBCINTR53=89
      NMASCOAGOCINTR53=90
      NMASCOAGBCINTR54=91
      NMASCOAGOCINTR54=92
!
      NMASAGEDSUINTR52=93
      NMASAGEDBCINTR52=94
      NMASAGEDOCINTR52=95
      NMASAGEDSOINTR52= 0 ! stored in NMASAGEDOCINTR52
!
      NMASMERGSUINTR12=96
      NMASMERGOCINTR12=97
      NMASMERGSOINTR12= 0 ! stored in NMASMERGOCINTR12
      NMASMERGSUINTR23=98
      NMASMERGBCINTR23=99
      NMASMERGOCINTR23=100
      NMASMERGSOINTR23=  0 ! stored in NMASMERGOCINTR12
      NMASMERGSUINTR34=101
      NMASMERGSSINTR34=102
      NMASMERGBCINTR34=103
      NMASMERGOCINTR34=104
      NMASMERGSOINTR34=  0 ! stored in NMASMERGOCINTR34
      NMASPROCSUINTR23=105
      NMASPROCBCINTR23=106
      NMASPROCOCINTR23=107
      NMASPROCSOINTR23=  0 ! stored in NMASPROCOCINTR23
!
!----------------------------------------------------------------
!
! EXTRA DIAGNOSTIC INDICES
!
! .. below are 4 water content & 5 cloud, precip & windspd fields (TRAQU)
      NWTCNT1    = 1
      NWTCNT2    = 2
      NWTCNT3    = 3
      NWTCNT4    = 4
      NLCFRAC    = 5
      NLWC       = 6
      NRAINLS    = 7
      NRAINCS    = 8
      NWINDSP    = 9
!
!----------------------------------------------------------------

      END SUBROUTINE UKCA_ORGV1_SUSSBCOC_5MODE

      SUBROUTINE UKCA_ORGV1_SUSSBCOCSO_5MODE

      IMPLICIT NONE

!---------------------------------------------------------------
!
! MAIN ARRAY LENGTHS AND SWITCHES
!
      NTRAER=23          ! # of aerosol advected tracers
      NCHEMG=11          ! # of gas phase chemistry tracers
      ICHEM=1            ! 1/0 = gas phase chemistry tracers on/off
      BUDGET=1           ! 1/0 = budget fields on/off
      TRAQU=1            ! 1/0 = extra diagnostic terms on/off
      IDUSTDEP=0         ! 1/0 = size-resolved dustdep fields on/off
      NOFFOX=3*ICHEM     ! # of offline oxidant species
      NDUSTDEP=8         ! # of size-resolved dustdep fields
      NBUDCHEM=26*ICHEM  ! # of gas chem budget fields
      NBUDAER=123        ! # of aerosol budget fields

      NADVG=2+NCHEMG*ICHEM ! # of gas phase advected tracers
      NTRAG=NADVG+NOFFOX   ! total # of gas phase species

      NBUDAERTOT=NBUDAER+NDUSTDEP*IDUSTDEP ! # of dustdep fields
      NBUDGET=NBUDCHEM+NBUDAERTOT ! total # of budget terms
      NTRAQU=TRAQU*9 ! # of extra diagnostic terms

! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
!                   NTRAER+NTRAG+NBUDGET+NTRAQU             = NVTOT
! .. ICHEM=1 SUSSBCOC  23    16    149     9 = 39 + 149 + 9 = 197 (orgv1)
!
!                   NTRAER+NADVG  NTRA
! .. ICHEM=1 SUSSBCOC 23    13  =  36
!
!
!-----------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
!
      IF(ICHEM == 0) THEN
! .. below are 2 fields if chemistry switched off
       MSOTWO  = 0 ! not included
       MMeSMe  = 0 ! not included
       MH2SO4  = 0 ! not included
       MDMSO   = 0 ! not included
       MMSA    = 0 ! not included
       MCS2    = 0 ! not included
       MH2S    = 0 ! not included
       MCOS    = 0 ! not included
       MMONOTER= 0 ! not included
       MSEC_ORG= 0 ! not included
       MH2O2   = 0 ! not included
       MQ3D    = 1
       MPT     = 2
      ENDIF
!
      IF(ICHEM == 1) THEN
! .. below are 13 fields for orgv1.d
       MSOTWO  = 1
       MMeSMe  = 2
       MH2SO4  = 3
       MDMSO   = 4
       MMSA    = 5
       MCS2    = 6
       MH2S    = 7
       MCOS    = 8
       MMONOTER= 9
       MSEC_ORG=10
       MH2O2   =11
       MQ3D    =12
       MPT     =13
      ENDIF
!
! .. molecular weights (kg/mol) for gases for orgv1
      MM_GAS=(/0.064,0.062,0.098,0.078,0.096,0.076,0.034,0.060,         &
               0.136,0.150,0.034/) ! 8 S species, terp, sec_org, H2O2
!
      CONDENSABLE_CHOICE=(/0,0,1,0,0,0,0,0,0,6,0/)
! .. H2SO4   to condense into 1st aerosol component (CP_SU)
! .. SEC_ORG to condense into 6th aerosol component (CP_SO)

      CONDENSABLE=(CONDENSABLE_CHOICE > 0)
!
      DIMEN=(/0.0,0.0,4.5e-10,0.0,0.0,0.0,0.0,0.0,0.0,4.5e-10,0.0/)
! Molecular diameters of condensable gas phase species (m)
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY
!
      IF(ICHEM == 0) THEN
! .. below are 2 fields if chemistry switched off
       NSOTWO  = 0 ! not included
       NMeSMe  = 0 ! not included
       NH2SO4  = 0 ! not included
       NDMSO   = 0 ! not included
       NMSA    = 0 ! not included
       NCS2    = 0 ! not included
       NH2S    = 0 ! not included
       NCOS    = 0 ! not included
       NMONOTER= 0 ! not included
       NSEC_ORG= 0 ! not included
       NH2O2   = 0 ! not included
       NO3     = 0 ! not included
       NOH     = 0 ! not included
       NNO3    = 0 ! not included
       NQ3D    = 1
       NPT     = 2
      ENDIF
!
      IF(ICHEM == 1) THEN
! .. below are 16 fields for orgv1.d
       NSOTWO  = 1
       NMeSMe  = 2
       NH2SO4  = 3
       NDMSO   = 4
       NMSA    = 5
       NCS2    = 6
       NH2S    = 7
       NCOS    = 8
       NMONOTER= 9
       NSEC_ORG=10
       NH2O2   =11
       NO3     =12
       NOH     =13
       NNO3    =14
       NQ3D    =15
       NPT     =16
      ENDIF
!
!----------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 26 gas phase budget quantity indices for orgv1
      NDMSEMOC  = 1
      NDMSTEND  = 2
      NSO2EMAN  = 3
      NSO2EMBM  = 4
      NSO2EMVL  = 5
      NSO2TEND  = 6
      NSO2DDEP  = 7
      NSO2WDEP  = 8
      NH2SO4TEND= 9
      NH2SO4DDEP=10
      NCOSEMAN  =11
      NCOSEMOC  =12
      NCOSTEND  =13
      NCS2EMAN  =14
      NCS2EMOC  =15
      NCS2TEND  =16
      NDMSOTEND =17
      NDMSODDEP =18
      NMSATEND  =19
      NMSADDEP  =20
      NTERP_EM  =21
      NTERP_TEND=22
      NTERP_DDEP=23
      NSORG_TEND=24
      NSORG_DDEP=25
      NSORG_WDEP=26
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! .. below are 123 aerosol budget quantity indices for SUSSBCOCSO [SO in SO]
!
      NMASPRIMSUAITSOL= 1
      NMASPRIMSUACCSOL= 2
      NMASPRIMSUCORSOL= 3
      NMASPRIMSSACCSOL= 4
      NMASPRIMSSCORSOL= 5
      NMASPRIMBCAITSOL= 0 ! BC only emitted to insoluble
      NMASPRIMBCAITINS= 6
      NMASPRIMOCAITSOL= 7
      NMASPRIMOCAITINS= 8
!
      NMASDDEPSUNUCSOL= 9
      NMASDDEPSUAITSOL=10
      NMASDDEPSUACCSOL=11
      NMASDDEPSUCORSOL=12
      NMASDDEPSSACCSOL=13
      NMASDDEPSSCORSOL=14
      NMASDDEPBCAITSOL=15
      NMASDDEPBCACCSOL=16
      NMASDDEPBCCORSOL=17
      NMASDDEPBCAITINS=18
      NMASDDEPOCNUCSOL= 0  ! stored in NMASDDEPSONUCSOL
      NMASDDEPOCAITSOL=19
      NMASDDEPOCACCSOL=20
      NMASDDEPOCCORSOL=21
      NMASDDEPOCAITINS=22
      NMASDDEPSONUCSOL=23
      NMASDDEPSOAITSOL=24
      NMASDDEPSOACCSOL=25
      NMASDDEPSOCORSOL=26
!
      NMASNUSCSUNUCSOL=27
      NMASNUSCSUAITSOL=28
      NMASNUSCSUACCSOL=29
      NMASNUSCSUCORSOL=30
      NMASNUSCSSACCSOL=31
      NMASNUSCSSCORSOL=32
      NMASNUSCBCAITSOL=33
      NMASNUSCBCACCSOL=34
      NMASNUSCBCCORSOL=35
      NMASNUSCBCAITINS=36
      NMASNUSCOCNUCSOL= 0  ! stored in NMASNUSCSONUCSOL
      NMASNUSCOCAITSOL=37
      NMASNUSCOCACCSOL=38
      NMASNUSCOCCORSOL=39
      NMASNUSCOCAITINS=40
      NMASNUSCSONUCSOL=41
      NMASNUSCSOAITSOL=42
      NMASNUSCSOACCSOL=43
      NMASNUSCSOCORSOL=44
!
      NMASIMSCSUNUCSOL=45
      NMASIMSCSUAITSOL=46
      NMASIMSCSUACCSOL=47
      NMASIMSCSUCORSOL=48
      NMASIMSCSSACCSOL=49
      NMASIMSCSSCORSOL=50
      NMASIMSCBCAITSOL=51
      NMASIMSCBCACCSOL=52
      NMASIMSCBCCORSOL=53
      NMASIMSCBCAITINS=54
      NMASIMSCOCNUCSOL= 0 ! stored in NMASIMSCSONUCSOL
      NMASIMSCOCAITSOL=55
      NMASIMSCOCACCSOL=56
      NMASIMSCOCCORSOL=57
      NMASIMSCOCAITINS=58
      NMASIMSCSONUCSOL=59
      NMASIMSCSOAITSOL=60
      NMASIMSCSOACCSOL=61
      NMASIMSCSOCORSOL=62
!
      NMASCLPRSUAITSOL1=63
      NMASCLPRSUACCSOL1=64
      NMASCLPRSUCORSOL1=65
      NMASCLPRSUAITSOL2=66
      NMASCLPRSUACCSOL2=67
      NMASCLPRSUCORSOL2=68
!
      NMASCONDSUNUCSOL=69
      NMASCONDSUAITSOL=70
      NMASCONDSUACCSOL=71
      NMASCONDSUCORSOL=72
      NMASCONDSUAITINS=73
      NMASNUCLSUNUCSOL=74
      NMASCONDOCNUCSOL= 0 ! stored in NMASCONDSONUCSOL
      NMASCONDOCAITSOL= 0 ! stored in NMASCONDSOAITSOL
      NMASCONDOCACCSOL= 0 ! stored in NMASCONDSOACCSOL
      NMASCONDOCCORSOL= 0 ! stored in NMASCONDSOCORSOL
      NMASCONDOCAITINS= 0 ! stored in NMASCONDSOAITINS
      NMASCONDSONUCSOL=75
      NMASCONDSOAITSOL=76
      NMASCONDSOACCSOL=77
      NMASCONDSOCORSOL=78
      NMASCONDSOAITINS=79
!
      NMASCOAGSUINTR12=80
      NMASCOAGSUINTR13=81
      NMASCOAGSUINTR14=82
      NMASCOAGSUINTR15=83
      NMASCOAGOCINTR12= 0 ! stored in NMASCOAGSOINTR12
      NMASCOAGOCINTR13= 0 ! stored in NMASCOAGSOINTR13
      NMASCOAGOCINTR14= 0 ! stored in NMASCOAGSOINTR14
      NMASCOAGOCINTR15= 0 ! stored in NMASCOAGSOINTR15
      NMASCOAGSOINTR12=84
      NMASCOAGSOINTR13=85
      NMASCOAGSOINTR14=86
      NMASCOAGSOINTR15=87
      NMASCOAGSUINTR23=88
      NMASCOAGBCINTR23=89
      NMASCOAGOCINTR23=90
      NMASCOAGSOINTR23=91
      NMASCOAGSUINTR24=92
      NMASCOAGBCINTR24=93
      NMASCOAGOCINTR24=94
      NMASCOAGSOINTR24=95
      NMASCOAGSUINTR34=96
      NMASCOAGBCINTR34=97
      NMASCOAGOCINTR34=98
      NMASCOAGSSINTR34=99
      NMASCOAGSOINTR34=100
!
      NMASCOAGBCINTR53=101
      NMASCOAGOCINTR53=102
      NMASCOAGBCINTR54=103
      NMASCOAGOCINTR54=104
!
      NMASAGEDSUINTR52=105
      NMASAGEDBCINTR52=106
      NMASAGEDOCINTR52=107
      NMASAGEDSOINTR52=108
!
      NMASMERGSUINTR12=109
      NMASMERGOCINTR12=  0 ! separate SO component
      NMASMERGSOINTR12=110
      NMASMERGSUINTR23=111
      NMASMERGBCINTR23=112
      NMASMERGOCINTR23=113
      NMASMERGSOINTR23=114
      NMASMERGSUINTR34=115
      NMASMERGSSINTR34=116
      NMASMERGBCINTR34=117
      NMASMERGOCINTR34=118
      NMASMERGSOINTR34=119
      NMASPROCSUINTR23=120
      NMASPROCBCINTR23=121
      NMASPROCOCINTR23=122
      NMASPROCSOINTR23=123
!
!---------------------------------------------------------------
!
! EXTRA DIAGNOSTIC INDICES
!
! .. below are 4 water content & 5 cloud, precip & windspd fields (TRAQU)
      NWTCNT1    = 1
      NWTCNT2    = 2
      NWTCNT3    = 3
      NWTCNT4    = 4
      NLCFRAC    = 5
      NLWC       = 6
      NRAINLS    = 7
      NRAINCS    = 8
      NWINDSP    = 9
!
!----------------------------------------------------------------
      END SUBROUTINE UKCA_ORGV1_SUSSBCOCSO_5MODE
!
      SUBROUTINE UKCA_SV1_SUSS_4MODE
!
      IMPLICIT NONE
!
!-----------------------------------------------------------
!
! MAIN ARRAY LENGTHS AND SWITCHES
!
      NTRAER=10          ! # of aerosol advected tracers
      NCHEMG=9           ! # of gas phase chemistry tracers
      ICHEM=1            ! 1/0 = gas phase chemistry tracers on/off
      BUDGET=1           ! 1/0 = budget fields on/off
      TRAQU=1            ! 1/0 = extra diagnostic terms on/off
      IDUSTDEP=0         ! 1/0 = size-resolved dustdep fields on/off
      NOFFOX=2*ICHEM     ! # of offline oxidant species
      NDUSTDEP=8         ! # of size-resolved dustdep fields
      NBUDCHEM=20*ICHEM  ! # of gas chem budget fields
      NBUDAER=46         ! # of aerosol budget fields
!
      NADVG=2+NCHEMG*ICHEM ! # of gas phase advected tracers
      NTRAG=NADVG+NOFFOX   ! total # of gas phase species
!
      NBUDAERTOT=NBUDAER+NDUSTDEP*IDUSTDEP ! # of dustdep fields
      NBUDGET=NBUDCHEM+NBUDAERTOT ! total # of budget terms
      NTRAQU=TRAQU*9 ! # of extra diagnostic terms
!
! When used in TOMCAT need to set NVTOT and NTRA in main run script
!
!                   NTRAER+NTRAG+NBUDGET+NTRAQU             = NVTOT
! .. ICHEM=1 SUSS      10    13     66     9 = 23 +  66 + 9 =  98 (sv1)
!
!                   NTRAER+NADVG  NTRA
! .. ICHEM=1 SUSS      10   11  =  21
!
!
!-----------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR S0 ARRAY
!
!
      IF(ICHEM == 0) THEN
! .. below are 2 fields if chemistry switched off
       MSOTWO  = 0 ! not included
       MMeSMe  = 0 ! not included
       MH2SO4  = 0 ! not included
       MDMSO   = 0 ! not included
       MMSA    = 0 ! not included
       MCS2    = 0 ! not included
       MH2S    = 0 ! not included
       MCOS    = 0 ! not included
       MMONOTER= 0 ! not included
       MSEC_ORG= 0 ! not included
       MH2O2   = 0 ! not included
       MQ3D    = 1
       MPT     = 2
      ENDIF
!
      IF(ICHEM == 1) THEN
! .. below are 13 fields for sv1.d
       MSOTWO  = 1
       MMeSMe  = 2
       MH2SO4  = 3
       MDMSO   = 4
       MMSA    = 5
       MCS2    = 6
       MH2S    = 7
       MCOS    = 8
       MMONOTER= 0
       MSEC_ORG= 0
       MH2O2   = 9
       MQ3D    =10
       MPT     =11
      ENDIF
!
! .. molecular weights (kg/mol) for gases for sv1
      MM_GAS=(/0.064,0.062,0.098,0.078,0.096,0.076,0.034,0.060,         &
               0.034,0.000,0.000/) ! 8 S spec & H2O2 (+ 2 dummys -> 11)
!
      CONDENSABLE_CHOICE=(/0,0,1,0,0,0,0,0,0,0,0/)
! .. H2SO4   to condense into 1st aerosol component (CP_SU)
!
      CONDENSABLE=(CONDENSABLE_CHOICE > 0)
!
      DIMEN=(/0.0,0.0,4.5e-10,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
! Molecular diameters of condensable gas phase species (m)
!
!------------------------------------------------------------
!
! GAS PHASE TRACER INDICES FOR ST ARRAY

      IF(ICHEM == 0) THEN
! .. below are 2 fields if chemistry switched off
       NSOTWO  = 0 ! not included
       NMeSMe  = 0 ! not included
       NH2SO4  = 0 ! not included
       NDMSO   = 0 ! not included
       NMSA    = 0 ! not included
       NCS2    = 0 ! not included
       NH2S    = 0 ! not included
       NCOS    = 0 ! not included
       NMONOTER= 0 ! not included
       NSEC_ORG= 0 ! not included
       NH2O2   = 0 ! not included
       NO3     = 0 ! not included
       NOH     = 0 ! not included
       NNO3    = 0 ! not included
       NQ3D    = 1
       NPT     = 2
      ENDIF
!
      IF(ICHEM == 1) THEN
! .. below are 13 fields for sv1.d
       NSOTWO  = 1
       NMeSMe  = 2
       NH2SO4  = 3
       NDMSO   = 4
       NMSA    = 5
       NCS2    = 6
       NH2S    = 7
       NCOS    = 8
       NMONOTER= 0
       NSEC_ORG= 0
       NH2O2   = 9
       NO3     = 0
       NOH     =10
       NNO3    =11
       NQ3D    =12
       NPT     =13
      ENDIF
!
!---------------------------------------------------------------
!
! GAS PHASE BUDGET INDICES
!
! .. below are 20 gas phase budget quantity indices for sv1
      NDMSEMOC  = 1
      NDMSTEND  = 2
      NSO2EMAN  = 3
      NSO2EMBM  = 4
      NSO2EMVL  = 5
      NSO2TEND  = 6
      NSO2DDEP  = 7
      NSO2WDEP  = 8
      NH2SO4TEND= 9
      NH2SO4DDEP=10
      NCOSEMAN  =11
      NCOSEMOC  =12
      NCOSTEND  =13
      NCS2EMAN  =14
      NCS2EMOC  =15
      NCS2TEND  =16
      NDMSOTEND =17
      NDMSODDEP =18
      NMSATEND  =19
      NMSADDEP  =20
      NTERP_EM  = 0
      NTERP_TEND= 0
      NTERP_DDEP= 0
      NSORG_TEND= 0
      NSORG_DDEP= 0
      NSORG_WDEP= 0
!
!------------------------------------------------------------------
!
! AEROSOL PHASE BUDGET INDICES
!
! .. below are 46 aerosol budget variables for SUSS aerosol system
!
! .. below are 123 aerosol budget quantity indices for SUSSBCOCSO [SO in SO]
!
      NMASPRIMSUAITSOL= 1
      NMASPRIMSUACCSOL= 2
      NMASPRIMSUCORSOL= 3
      NMASPRIMSSACCSOL= 4
      NMASPRIMSSCORSOL= 5
      NMASPRIMBCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASPRIMBCAITINS= 0 ! no BC/OC/SO in this setup
      NMASPRIMOCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASPRIMOCAITINS= 0 ! no BC/OC/SO in this setup
!
      NMASDDEPSUNUCSOL= 6
      NMASDDEPSUAITSOL= 7
      NMASDDEPSUACCSOL= 8
      NMASDDEPSUCORSOL= 9
      NMASDDEPSSACCSOL=10
      NMASDDEPSSCORSOL=11
      NMASDDEPBCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPBCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPBCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPBCAITINS= 0 ! no BC/OC/SO in this setup
      NMASDDEPOCNUCSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPOCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPOCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPOCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPOCAITINS= 0 ! no BC/OC/SO in this setup
      NMASDDEPSONUCSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPSOAITSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPSOACCSOL= 0 ! no BC/OC/SO in this setup
      NMASDDEPSOCORSOL= 0 ! no BC/OC/SO in this setup
!
      NMASNUSCSUNUCSOL=12
      NMASNUSCSUAITSOL=13
      NMASNUSCSUACCSOL=14
      NMASNUSCSUCORSOL=15
      NMASNUSCSSACCSOL=16
      NMASNUSCSSCORSOL=17
      NMASNUSCBCAITSOL= 0
      NMASNUSCBCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCBCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCBCAITINS= 0 ! no BC/OC/SO in this setup
      NMASNUSCOCNUCSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCOCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCOCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCOCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCOCAITINS= 0 ! no BC/OC/SO in this setup
      NMASNUSCSONUCSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCSOAITSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCSOACCSOL= 0 ! no BC/OC/SO in this setup
      NMASNUSCSOCORSOL= 0 ! no BC/OC/SO in this setup
!
      NMASIMSCSUNUCSOL=18
      NMASIMSCSUAITSOL=19
      NMASIMSCSUACCSOL=20
      NMASIMSCSUCORSOL=21
      NMASIMSCSSACCSOL=22
      NMASIMSCSSCORSOL=23
      NMASIMSCBCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCBCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCBCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCBCAITINS= 0 ! no BC/OC/SO in this setup
      NMASIMSCOCNUCSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCOCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCOCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCOCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCOCAITINS= 0 ! no BC/OC/SO in this setup
      NMASIMSCSONUCSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCSOAITSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCSOACCSOL= 0 ! no BC/OC/SO in this setup
      NMASIMSCSOCORSOL= 0 ! no BC/OC/SO in this setup
!
      NMASCLPRSUAITSOL1=24
      NMASCLPRSUACCSOL1=25
      NMASCLPRSUCORSOL1=26
      NMASCLPRSUAITSOL2=27
      NMASCLPRSUACCSOL2=28
      NMASCLPRSUCORSOL2=29
!
      NMASCONDSUNUCSOL=30
      NMASCONDSUAITSOL=31
      NMASCONDSUACCSOL=32
      NMASCONDSUCORSOL=33
      NMASCONDSUAITINS= 0 ! only soluble modes in this setup
      NMASNUCLSUNUCSOL=34
      NMASCONDOCNUCSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDOCAITSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDOCACCSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDOCCORSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDOCAITINS= 0 ! no BC/OC/SO in this setup
      NMASCONDSONUCSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDSOAITSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDSOACCSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDSOCORSOL= 0 ! no BC/OC/SO in this setup
      NMASCONDSOAITINS= 0 ! no BC/OC/SO in this setup
!
      NMASCOAGSUINTR12=35
      NMASCOAGSUINTR13=36
      NMASCOAGSUINTR14=37
      NMASCOAGSUINTR15= 0 ! only soluble modes in this setup
      NMASCOAGOCINTR12= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR13= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR14= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR15= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR12= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR13= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR14= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR15= 0 ! no BC/OC/SO in this setup
      NMASCOAGSUINTR23=38
      NMASCOAGBCINTR23= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR23= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR23= 0 ! no BC/OC/SO in this setup
      NMASCOAGSUINTR24=39
      NMASCOAGBCINTR24= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR24= 0 ! no BC/OC/SO in this setup
      NMASCOAGSOINTR24= 0 ! no BC/OC/SO in this setup
      NMASCOAGSUINTR34=40
      NMASCOAGBCINTR34= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR34= 0 ! no BC/OC/SO in this setup
      NMASCOAGSSINTR34=41
      NMASCOAGSOINTR34= 0 ! no BC/OC/SO in this setup
!
      NMASCOAGBCINTR53= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR53= 0 ! no BC/OC/SO in this setup
      NMASCOAGBCINTR54= 0 ! no BC/OC/SO in this setup
      NMASCOAGOCINTR54= 0 ! no BC/OC/SO in this setup
!
      NMASAGEDSUINTR52= 0 ! only soluble modes in this setup
      NMASAGEDBCINTR52= 0 ! no BC/OC/SO in this setup
      NMASAGEDOCINTR52= 0 ! no BC/OC/SO in this setup
      NMASAGEDSOINTR52= 0 ! no BC/OC/SO in this setup
!
      NMASMERGSUINTR12=42
      NMASMERGOCINTR12= 0 ! no BC/OC/SO in this setup
      NMASMERGSOINTR12= 0 ! no BC/OC/SO in this setup
      NMASMERGSUINTR23=43
      NMASMERGBCINTR23= 0 ! no BC/OC/SO in this setup
      NMASMERGOCINTR23= 0 ! no BC/OC/SO in this setup
      NMASMERGSOINTR23= 0 ! no BC/OC/SO in this setup
      NMASMERGSUINTR34=44
      NMASMERGSSINTR34=45
      NMASMERGBCINTR34= 0 ! no BC/OC/SO in this setup
      NMASMERGOCINTR34= 0 ! no BC/OC/SO in this setup
      NMASMERGSOINTR34= 0 ! no BC/OC/SO in this setup
      NMASPROCSUINTR23=46
      NMASPROCBCINTR23= 0 ! no BC/OC/SO in this setup
      NMASPROCOCINTR23= 0 ! no BC/OC/SO in this setup
      NMASPROCSOINTR23= 0 ! no BC/OC/SO in this setup
!
!---------------------------------------------------------------
!
! EXTRA DIAGNOSTIC INDICES
!
! .. below are 4 water content & 5 cloud, precip & windspd fields (TRAQU)
      NWTCNT1    = 1
      NWTCNT2    = 2
      NWTCNT3    = 3
      NWTCNT4    = 4
      NLCFRAC    = 5
      NLWC       = 6
      NRAINLS    = 7
      NRAINCS    = 8
      NWINDSP    = 9
!
      END SUBROUTINE UKCA_SV1_SUSS_4MODE

      END MODULE UKCA_SETUP_INDICES
#endif
