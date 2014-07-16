#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Subroutine to define fields required from D1 by UKCA.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!  Set D1 section and item codes depending on the selected
!  chemistry.
!
! Current Code Owner:       Colin Johnson/Olaf Morgenstern
!                           Fiona O'Connor
!
!
!  Code Description:
!   Language:  FORTRAN 90 (formatted)
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
      SUBROUTINE UKCA_SETD1DEFS(row_length,rows,n_rows,model_levels,   &
                 bl_levels,tr_levels,wet_levels,land_pts,sm_levels,    &
                 ntiles,tr_ukca)

      USE UKCA_D1_DEFS
      IMPLICIT NONE

#include "c_mdi.h"
#include "version.h"
#include "csubmodl.h"
#include "cstash.h"
#include "cmaxsize.h"
#include "nstypes.h"
#include "model.h"
#include "cntlatm.h"
#include "cruntimc.h"

      INTEGER, INTENT(IN) :: row_length     ! length of row
      INTEGER, INTENT(IN) :: rows           ! number of rows
      INTEGER, INTENT(IN) :: n_rows         ! number of rows on v grid
      INTEGER, INTENT(IN) :: model_levels   ! number of levels
      INTEGER, INTENT(IN) :: bl_levels      ! number of levels in BL
      INTEGER, INTENT(IN) :: tr_levels      ! number of tracer levels
      INTEGER, INTENT(IN) :: wet_levels     ! number of wet levels
      INTEGER, INTENT(IN) :: land_pts       ! no of land points
      INTEGER, INTENT(IN) :: sm_levels      ! no of soil moisture levels
      INTEGER, INTENT(IN) :: ntiles         ! no of land tile types
      INTEGER, INTENT(IN) :: tr_ukca        ! no of activated tracers

      INTEGER :: I,J,idiag                  ! counters

      CHARACTER (LEN=70) :: cmessage        ! Error message

!     Set fluxdiags to zero initially
      n_BE_fluxdiags    = 0
      nmax_strat_fluxdiags = 0

      IF (L_UKCA_trop) THEN

! Standard tropospheric chemistry
! ===============================
        n_chem_emissions = 8
        n_3d_emissions = 1       ! aircraft NOX
        n_aero_tracers = 0
        IF (L_ukca_BEflux) n_BE_fluxdiags = 170
        ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
        IF (L_UKCA_family) THEN
          n_chem_tracers = 24
          n_chem_diags   = 20
          nr_therm       = 0
          nr_phot        = 0
          em_chem_spec =                                               &
        (/'NOx       ','CH4       ','CO        ','HCHO      ',         &
          'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',         &
          'NO_aircrft'/)
        ELSE
          n_chem_tracers =  26
          n_chem_diags   =  20
          nr_therm       = 102        ! thermal reactions
          nr_phot        = 27         ! photolytic ---"---
          em_chem_spec =                                               &
        (/'NO        ','CH4       ','CO        ','HCHO      ',         &
          'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',         &
          'NO_aircrft'/)
        ENDIF
      ELSE IF (L_UKCA_tropisop) THEN

! Std tropospheric chemistry with MIM isoprene scheme
! ===================================================
        n_chem_emissions = 9
        n_3d_emissions = 1       ! aircraft NOX
        n_aero_tracers = 0
        IF (L_ukca_BEflux) n_BE_fluxdiags = 170  ! Needs revision
        ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
        IF (L_UKCA_family) THEN
          n_chem_tracers = 36
          n_chem_diags   = 20
          nr_therm       = 0
          nr_phot        = 0
          em_chem_spec =                                               &
        (/'NOx       ','CH4       ','CO        ','HCHO      ',         &
          'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',         &
          'C5H8      ','NO_aircrft'/)
        ELSE
          n_chem_tracers = 38
          n_chem_diags   = 20
          nr_therm       = 137
          nr_phot        = 37
          em_chem_spec =                                               &
        (/'NO        ','CH4       ','CO        ','HCHO      ',         &
          'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',         &
          'C5H8      ','NO_aircrft'/)
        ENDIF
      ELSE IF (L_ukca_aerchem) THEN

! Std trop chem with SO2 and DMS
! ==============================
        IF (L_ukca_BEflux) n_BE_fluxdiags    = 170
        IF (L_UKCA_family) THEN
          cmessage='Aerosol chemistry not available with family'
          CALL EREPORT('UKCA_SETD1DEFS1',1,cmessage)
        ELSE
          n_chem_emissions = 14       ! Surface/ high-level emissions
          n_3d_emissions = 2          ! SO2_nat, aircraft NOX
          n_chem_tracers = 26         ! advected chemical tracers
          n_aero_tracers =  4         ! advected aerochem ---"---
          n_chem_diags   = 20         ! diagnosed species
          nr_therm       = 110        ! thermal reactions
          nr_phot        = 27         ! photolytic ---"---
          ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
          em_chem_spec =                                                &
        (/'NO        ','CH4       ','CO        ','HCHO      ',          &
          'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',          &
          'SO2_low   ','SO2_high  ','DMS       ','SO2_nat   ',          &
          'soot_high ','Bmass_low ','Bmass_high','NO_aircrft'/)
        ENDIF
      ELSE IF (L_UKCA_strat .OR. L_ukca_strattrop .OR. L_UKCA_stratcfc) &
        THEN

! Stratospheric chemistry
! =======================
        n_chem_emissions = 8
        n_3d_emissions = 1       ! aircraft NOX
        n_aero_tracers =  0
        ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
        IF (L_UKCA_family) THEN
          IF (L_ukca_strat) THEN
            n_chem_tracers = 25
          ELSE IF (L_ukca_strattrop) THEN
            n_chem_tracers = 39
          ELSE IF (L_ukca_stratcfc) THEN
            n_chem_tracers = 32
          END IF
          n_chem_diags   = 20
          nr_therm       = 0
          nr_phot        = 0
          em_chem_spec =                                                &
        (/'NOx       ','CH4       ','CO        ','HCHO      ',          &
          'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',          &
          'NO_aircrft'/)
        ELSE    ! non-family
          IF (L_ukca_strat) THEN
            n_chem_tracers = 36
          ELSE IF (L_ukca_strattrop) THEN
            n_chem_tracers = 56
          ELSE IF (L_ukca_stratcfc) THEN
            n_chem_tracers = 43
          END IF
          n_chem_diags   =  20
          nr_therm       = 102        ! NOT SURE ABOUT THIS (Olaf)!
          nr_phot        = 27         !
          em_chem_spec =                                                &
        (/'NO        ','CH4       ','CO        ','HCHO      ',          &
          'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',          &
          'NO_aircrft'/)
        ENDIF
      nmax_strat_fluxdiags = n_chem_tracers
      ELSE

! No chemistry
! ============
        n_chem_emissions = 0
        n_chem_tracers   = 0
        n_chem_diags     = 0
      ENDIF

      IF (l_UKCA_RnPb) THEN
        n_RnPb_emissions = 1
        n_RnPb_tracers   = 3
        ALLOCATE(em_RnPb_spec(n_RnPb_emissions))
        em_RnPb_spec(1:n_RnPb_emissions) = (/'Rn-222    '/)
      ELSE
        n_RnPb_emissions = 0
        n_RnPb_tracers   = 0
      ENDIF

      IF (L_UKCA_dust) THEN
        n_dust_emissions = 6
        n_dust_tracers   = 6
        ALLOCATE(em_dust_spec(n_dust_emissions))
        em_dust_spec(1:n_dust_emissions)=                              &
        (/'Dust_div_1','Dust_div_2','Dust_div_3','Dust_div_4',         &
          'Dust_div_5','Dust_div_6'/)
      ELSE
        n_dust_emissions = 0
        n_dust_tracers   = 0
      ENDIF

      IF (L_UKCA_MODE) THEN
        IF (.NOT. L_ukca_aerchem) THEN
          cmessage=' L_ukca_aerchem is required to run UKCA_MODE'
          CALL EREPORT('UKCA_SETD1DEFS',1,cmessage)
        ENDIF
        n_MODE_emissions = 0
        n_MODE_tracers   = 25
      ELSE
        n_MODE_emissions = 0
        n_MODE_tracers   = 0
      ENDIF

      n_use_tracers   = n_chem_tracers + n_mode_tracers +               &
                        n_aero_tracers + n_dust_tracers +               &
                        n_rnpb_tracers
      n_use_emissions = n_chem_emissions + n_mode_emissions +           &
                        n_3d_emissions   + n_dust_emissions +           &
                        n_RnPb_emissions

      n_in_progs     =  37     ! max no of progs reqd other than tracers/ems
      n_in_diags0    =   4     ! max no of diags (sect 0) reqd
      n_in_diags1    =   2     ! max no of diags (sect 1) reqd
      n_in_diags3    =  26     ! max no of diags (sect 3) reqd
      n_in_diags4    =   5     ! max no of diags (sect 4) reqd
      n_in_diags5    =   2     ! max no of diags (sect 5) reqd
      n_in_diags8    =   1     ! max no of diags (sect 8) reqd
      n_in_diags15   =   1     ! max no of diags (sect 15) reqd
      n_in_diags30   =   1     ! max no of diags (sect 30) reqd
      n_in_diags33   = n_chem_diags + n_BE_fluxdiags             &
        + nmax_strat_fluxdiags ! max no of diags (sect 33) reqd


      n_emiss_first  = 301     ! first stash no for UKCA emissions
      n_emiss_last   = 320     !  last stash no for UKCA emissions

!     Create species names array

      ALLOCATE(nm_spec(n_all_tracers))
      IF (L_UKCA_family) THEN
        nm_spec(1:n_all_tracers) = (/                                     &
       'Ox        ','NOx       ','XXX       ','XXX       ','N2O5      ',  &
       'HO2NO2    ','HNO3      ','H2O2      ','CH4       ','CO        ',  &
       'HCHO      ','MeOOH     ','HONO      ','C2H6      ','ETOOH     ',  &
       'MeCHO     ','PAN       ','C3H8      ','N-PrOOH   ','I-PrOOH   ',  &
       'EtCHO     ','Me2CO     ','MeCOCH2OOH','PPAN      ','MeONO2    ',  &
       'O3S       ','C5H8      ','ISOOH     ','ISON      ','MACR      ',  &
       'MACROOH   ','MPAN      ','HACET     ','MGLY      ','NALD      ',  &
       'HCOOH     ','MeCO3H    ','MeCO2H    ','XXX       ','XXX       ',  &
       'ClOx      ','XXX       ','XXX       ','OClO      ','BrOx      ',  &
       'XXX       ','BrCl      ','BrONO2    ','N2O       ','HCl       ',  &
       'HOCl      ','HBr       ','HOBr      ','ClONO2    ','CFCl3     ',  &
       'CF2Cl2    ','MeBr      ','XXX       ','XXX       ','XXX       ',  &
       'MeCl      ','CF2ClBr   ','CCl4      ','CF2ClCFCl2','CHF2Cl    ',  &
       'MeCCl3    ','CF3Br     ','H2OS      ','XXX       ','H2        ',  &
       'SO2       ','H2SO4     ','DMS       ','DMSO      ','MSA       ',  &
       'H2S       ','CS2       ','COS       ','NH3       ','XXX       ',  &
       'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  &
       'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  &
       'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  &
       'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  &
       'ND_Nuc_SOL','ND_Ait_SOL','ND_Acc_SOL','ND_Cor_SOL','ND_Ait_INS',  &
       'ND_Acc_INS','ND_Cor_INS','Nuc_SOL_Su','Ait_SOL_Su','Ait_SOL_BC',  &
       'Ait_SOL_OC','Acc_SOL_Su','Acc_SOL_BC','Acc_SOL_OC','Acc_SOL_SS',  &
       'Acc_SOL_Du','Cor_SOL_Su','Cor_SOL_BC','Cor_SOL_OC','Cor_SOL_SS',  &
       'Cor_SOL_Du','Ait_INS_BC','Ait_INS_OC','Acc_INS_Du','Cor_INS_Du',  &
       'NH4+      ','NO3-      ','XXX       ','XXX       ','XXX       ',  &
       'Dust_Div_1','Dust_Div_2','Dust_Div_3','Dust_Div_4','Dust_Div_5',  &
       'Dust_Div_6','XXX       ','XXX       ','Rn-222    ','Pb-210    ',  &
       'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  &
       'XXX       ','XXX       ','XXX       ','XXX       ','XXX       '   &
        /)
      ELSE  ! non family tracers
        nm_spec(1:n_all_tracers) = (/                                     &
       'O3        ','NO        ','NO2       ','NO3       ','N2O5      ',  &
       'HO2NO2    ','HNO3      ','H2O2      ','CH4       ','CO        ',  &
       'HCHO      ','MeOOH     ','HONO      ','C2H6      ','ETOOH     ',  &
       'MeCHO     ','PAN       ','C3H8      ','N-PrOOH   ','I-PrOOH   ',  &
       'EtCHO     ','Me2CO     ','MeCOCH2OOH','PPAN      ','MeONO2    ',  &
       'O3S       ','C5H8      ','ISOOH     ','ISON      ','MACR      ',  &
       'MACROOH   ','MPAN      ','HACET     ','MGLY      ','NALD      ',  &
       'HCOOH     ','MeCO3H    ','MeCO2H    ','XXX       ','XXX       ',  &
       'Cl        ','ClO       ','Cl2O2     ','OClO      ','Br        ',  &
       'BrO       ','BrCl      ','BrONO2    ','N2O       ','HCl       ',  &
       'HOCl      ','HBr       ','HOBr      ','ClONO2    ','CFCl3     ',  &
       'CF2Cl2    ','MeBr      ','N         ','O3P       ','XXX       ',  &
       'MeCl      ','CF2ClBr   ','CCl4      ','CF2ClCFCl2','CHF2Cl    ',  &
       'MeCCl3    ','CF3Br     ','H2OS      ','XXX       ','H2        ',  &
       'DMS       ','SO2       ','H2SO4     ','MSA       ','NH3       ',  &
       'H2S       ','CS2       ','COS       ','NH4NO3    ','H         ',  &
       'OH        ','HO2       ','MeOO      ','EtOO      ','MeCO3     ',  &
       'n-PrOO    ','i-PrOO    ','EtCO3     ','MeCOCH2OO ','XXX       ',  &
       'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  &
       'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  &
       'ND_Nuc_SOL','ND_Ait_SOL','ND_Acc_SOL','ND_Cor_SOL','ND_Ait_INS',  &
       'ND_Acc_INS','ND_Cor_INS','Nuc_SOL_Su','Ait_SOL_Su','Ait_SOL_BC',  &
       'Ait_SOL_OC','Acc_SOL_Su','Acc_SOL_BC','Acc_SOL_OC','Acc_SOL_SS',  &
       'Acc_SOL_Du','Cor_SOL_Su','Cor_SOL_BC','Cor_SOL_OC','Cor_SOL_SS',  &
       'Cor_SOL_Du','Ait_INS_BC','Ait_INS_OC','Acc_INS_Du','Cor_INS_Du',  &
       'NH4+      ','NO3-      ','XXX       ','XXX       ','XXX       ',  &
       'Dust_Div_1','Dust_Div_2','Dust_Div_3','Dust_Div_4','Dust_Div_5',  &
       'Dust_Div_6','XXX       ','XXX       ','Rn-222    ','Pb-210    ',  &
       'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  &
       'XXX       ','XXX       ','XXX       ','XXX       ','XXX       '   &
        /)
      ENDIF

!     Specify the section and item codes of prognostics and diagnostics
!     required from D1.  Set array dimensions, but ignore halos as these
!     are set in the call to UKCA_SET_ARRAY_BOUNDS from the halosize
!     array.  Set %prognostic and %required t/f.
!     Other components of UkcaD1Codes are read in from D1 address array
!     in UKCA_MAIN1.

      Nukca_D1items = n_use_tracers  + n_use_emissions +                &
                      n_in_progs     + n_in_diags0     +                &
                      n_in_diags1    + n_in_diags3     +                &
                      n_in_diags4    + n_in_diags5     +                &
                      n_in_diags8    + n_in_diags15    +                &
                      n_in_diags30   + n_in_diags33

      ALLOCATE(UkcaD1Codes(Nukca_D1items))

      UkcaD1Codes(:)%section=IMDI
      UkcaD1Codes(:)%item=IMDI
      UkcaD1Codes(:)%n_levels=IMDI
      UkcaD1Codes(:)%address=IMDI
      UkcaD1Codes(:)%length=IMDI
      UkcaD1Codes(:)%halo_type=IMDI
      UkcaD1Codes(:)%grid_type=IMDI
      UkcaD1Codes(:)%field_type=IMDI
      UkcaD1Codes(:)%len_dim1=IMDI
      UkcaD1Codes(:)%len_dim2=IMDI
      UkcaD1Codes(:)%len_dim3=IMDI
      UkcaD1Codes(:)%prognostic=.true.
      UkcaD1Codes(:)%required=.false.

!     Only retain names for active tracers and set required logical

      WRITE(6,*) ' UKCA: The following tracers were selected:'
      J=0
      DO I=1,n_all_tracers
        IF (TC_UKCA(I) /= 1) THEN
          nm_spec(I)='XXX       '
        ELSE
          J = J+1
          UkcaD1Codes(J)%section    = UKCA_sect
          UkcaD1Codes(J)%item       = I
          UkcaD1Codes(J)%len_dim1   = row_length
          UkcaD1Codes(J)%len_dim2   = rows
          UkcaD1Codes(J)%len_dim3   = tr_levels
          UkcaD1Codes(J)%prognostic = .true.
          UkcaD1Codes(J)%required   = .true.
          WRITE(6,*) nm_spec(I)
        ENDIF
      ENDDO

      IF (J /= tr_ukca) THEN
        cmessage='UKCA: wrong number of tracers'
        WRITE(6,*) 'Expected: ',tr_ukca,' Encountered: ',J
! DEPENDS ON: ereport
        CALL EREPORT('UKCA_SETD1DEFS',1,cmessage)
      ENDIF

!     Prognostics from section 0, emissions required set above

      J = n_use_tracers
      IF (n_chem_emissions+n_3d_emissions+n_mode_emissions > 0) THEN
        DO i=1,n_chem_emissions + n_3d_emissions
          UkcaD1Codes(J+i)%section    = 0
          UkcaD1Codes(J+i)%item       = 301+i-1          ! trop chemistry
          UkcaD1Codes(J+i)%len_dim1   = row_length       ! uses stash codes
          UkcaD1Codes(J+i)%len_dim2   = rows             ! 301-309 for
          UkcaD1Codes(J+i)%required   = .true.           ! surface emissions
          UkcaD1Codes(J+i)%prognostic = .true.           ! from Section 0
! Special cases, emissions already available in UM
          IF (em_chem_spec(i)(1:7) == 'SO2_low') THEN
            UkcaD1Codes(J+i)%item     = 58
            IF (.NOT. L_SO2_SURFEM .AND. L_ukca_aerchem) THEN
              cmessage='SO2 surface emissions from UM are not flagged'
              CALL EREPORT('UKCA_SETD1DEFS',58,cmessage)
            ENDIF
          ELSEIF (em_chem_spec(i)(1:7) == 'SO2_nat') THEN
            UkcaD1Codes(J+i)%item     = 121
            UkcaD1Codes(J+i)%len_dim3 = tr_levels
            IF (.NOT. L_SO2_NATEM .AND. L_ukca_aerchem) THEN
              cmessage='SO2 natural emissions from UM are not flagged'
              CALL EREPORT('UKCA_SETD1DEFS',121,cmessage)
            ENDIF
          ELSEIF (em_chem_spec(i)(1:8) == 'SO2_high') THEN
            UkcaD1Codes(J+i)%item     = 126
            IF (.NOT. L_SO2_HILEM .AND. L_ukca_aerchem) THEN
              cmessage='SO2 high-level emissions are not flagged'
              CALL EREPORT('UKCA_SETD1DEFS',126,cmessage)
            ENDIF
          ELSEIF (em_chem_spec(i)(1:3) == 'NH3') THEN
            UkcaD1Codes(J+i)%item     = 127
            IF (.NOT. L_NH3_EM .AND. L_ukca_aerchem) THEN
              cmessage='NH3 surface emissions from UM are not flagged'
              CALL EREPORT('UKCA_SETD1DEFS',127,cmessage)
            ENDIF
          ELSEIF (em_chem_spec(i)(1:8) == 'soot_low') THEN
            UkcaD1Codes(J+i)%item     = 128
            IF (.NOT. L_SOOT_SUREM .AND. L_ukca_aerchem) THEN
              cmessage='Soot surface emissions from UM are not flagged'
              CALL EREPORT('UKCA_SETD1DEFS',128,cmessage)
            ENDIF
          ELSEIF (em_chem_spec(i)(1:9) == 'soot_high') THEN
            UkcaD1Codes(J+i)%item     = 129
            IF (.NOT. L_SOOT_HILEM .AND. L_ukca_aerchem) THEN
              cmessage='Soot high-level emissions from UM not flagged'
              CALL EREPORT('UKCA_SETD1DEFS',129,cmessage)
            ENDIF
          ELSEIF (em_chem_spec(i)(1:9) == 'Bmass_low') THEN
            UkcaD1Codes(J+i)%item     = 130
            IF (.NOT. L_BMASS_SUREM .AND. L_ukca_aerchem) THEN
              cmessage='UM Biomass high-level emissions not flagged'
              CALL EREPORT('UKCA_SETD1DEFS',130,cmessage)
            ENDIF
          ELSEIF (em_chem_spec(i)(1:10) == 'Bmass_high') THEN
            UkcaD1Codes(J+i)%item     = 131
            IF (.NOT. L_BMASS_HILEM .AND. L_ukca_aerchem) THEN
              cmessage='UM Biomass high-level emissions not flagged'
              CALL EREPORT('UKCA_SETD1DEFS',131,cmessage)
            ENDIF
          ELSEIF (em_chem_spec(i)(1:3) == 'DMS') THEN
            UkcaD1Codes(J+i)%section  = 17
            UkcaD1Codes(J+i)%item     = 205
            UkcaD1Codes(J+i)%prognostic = .false.
            IF (.NOT. L_DMS_EM .AND. L_ukca_aerchem) THEN
              cmessage='DMS surface emissions from UM are not flagged'
              CALL EREPORT('UKCA_SETD1DEFS',17205,cmessage)
            ENDIF
          ELSEIF (em_chem_spec(i)(1:7) == 'NO_airc') THEN
            UkcaD1Codes(J+i)%item     = 340
            UkcaD1Codes(J+i)%len_dim3 = tr_levels
          ENDIF
        ENDDO
      ENDIF

      J = J + n_chem_emissions + n_3d_emissions
      IF (L_UKCA_RnPb) THEN
        DO i=1,n_RnPb_emissions
          UkcaD1Codes(J+i)%section    = 0
          UkcaD1Codes(J+i)%item       = 310+i-1          ! RnPb experiment
          UkcaD1Codes(J+i)%len_dim1   = row_length       ! uses stash code
          UkcaD1Codes(J+i)%len_dim2   = rows             ! 310 in Section 0
          UkcaD1Codes(J+i)%required   = .true.           ! for surface emissions
          UkcaD1Codes(J+i)%prognostic = .true.           ! of radon
        ENDDO
      ENDIF

      J = J + n_RnPb_emissions
      IF (L_UKCA_dust) THEN
        DO i=1,n_dust_emissions
          UkcaD1Codes(J+i)%section    = 0
          UkcaD1Codes(J+i)%item       = 311+i-1          ! Dust experiment
          UkcaD1Codes(J+i)%len_dim1   = row_length       ! uses stash code
          UkcaD1Codes(J+i)%len_dim2   = rows             ! 311-316 in Section 0
          UkcaD1Codes(J+i)%required   = .true.           ! surface emissions
          UkcaD1Codes(J+i)%prognostic = .true.           ! of dust
        ENDDO
      ENDIF

! Prognostic fields
      J = J + n_dust_emissions

      UkcaD1Codes(J+1:J+n_in_progs)%section    = 0
      UkcaD1Codes(J+1:J+n_in_progs)%prognostic = .true.
      UkcaD1Codes(J+1:J+n_in_progs)%required   = .true.
      UkcaD1Codes(J+1)%item=4              ! Potential Temperature
      UkcaD1Codes(J+1)%len_dim1=row_length
      UkcaD1Codes(J+1)%len_dim2=rows
      UkcaD1Codes(J+1)%len_dim3=model_levels
      UkcaD1Codes(J+2)%item=9              ! Soil Moisture
      UkcaD1Codes(J+2)%len_dim1=land_pts
      UkcaD1Codes(J+2)%len_dim2=sm_levels
      UkcaD1Codes(J+3)%item=10             ! Q
      UkcaD1Codes(J+3)%len_dim1=row_length
      UkcaD1Codes(J+3)%len_dim2=rows
      UkcaD1Codes(J+3)%len_dim3=wet_levels
      UkcaD1Codes(J+4)%item=12             ! QCF
      UkcaD1Codes(J+4)%len_dim1=row_length
      UkcaD1Codes(J+4)%len_dim2=rows
      UkcaD1Codes(J+4)%len_dim3=wet_levels
      UkcaD1Codes(J+5)%item=14             ! Conv cloud base level
      UkcaD1Codes(J+5)%len_dim1=row_length
      UkcaD1Codes(J+5)%len_dim2=rows
      UkcaD1Codes(J+6)%item=15             ! Conv cloud top level
      UkcaD1Codes(J+6)%len_dim1=row_length
      UkcaD1Codes(J+6)%len_dim2=rows
      UkcaD1Codes(J+7)%item=16             ! Conv cloud liquid water path
      UkcaD1Codes(J+7)%len_dim1=row_length
      UkcaD1Codes(J+7)%len_dim2=rows
      UkcaD1Codes(J+8)%item=24             ! Surface temperature
      UkcaD1Codes(J+8)%len_dim1=row_length
      UkcaD1Codes(J+8)%len_dim2=rows
      UkcaD1Codes(J+9)%item=25             ! Boundary layer height
      UkcaD1Codes(J+9)%len_dim1=row_length
      UkcaD1Codes(J+9)%len_dim2=rows
      UkcaD1Codes(J+10)%item=26             ! Roughness length
      UkcaD1Codes(J+10)%len_dim1=row_length
      UkcaD1Codes(J+10)%len_dim2=rows
      UkcaD1Codes(J+11)%item=30            ! Land sea mask
      UkcaD1Codes(J+11)%len_dim1=row_length
      UkcaD1Codes(J+11)%len_dim2=rows
      UkcaD1Codes(J+12)%item=31            ! Sea ice fraction
      UkcaD1Codes(J+12)%len_dim1=row_length
      UkcaD1Codes(J+12)%len_dim2=rows
      UkcaD1Codes(J+13)%item=60            ! Climatological ozone
      UkcaD1Codes(J+13)%len_dim1=1
      UkcaD1Codes(J+13)%len_dim2=rows
      UkcaD1Codes(J+13)%len_dim3=model_levels
      UkcaD1Codes(J+14)%item=103           ! SO4 Aitken Mode
      UkcaD1Codes(J+14)%len_dim1=row_length
      UkcaD1Codes(J+14)%len_dim2=rows
      UkcaD1Codes(J+14)%len_dim3=tr_levels
      UkcaD1Codes(J+15)%item=104           ! SO4 accumulation Mode
      UkcaD1Codes(J+15)%len_dim1=row_length
      UkcaD1Codes(J+15)%len_dim2=rows
      UkcaD1Codes(J+15)%len_dim3=tr_levels
      UkcaD1Codes(J+16)%item=211           ! Conv cloud amount
      UkcaD1Codes(J+16)%len_dim1=row_length
      UkcaD1Codes(J+16)%len_dim2=rows
      UkcaD1Codes(J+16)%len_dim3=model_levels
      UkcaD1Codes(J+17)%item=216           ! Fraction of surface types
      UkcaD1Codes(J+17)%len_dim1=land_pts
      UkcaD1Codes(J+17)%len_dim2=ntype
      UkcaD1Codes(J+18)%item=217           ! LAI of PFTs
      UkcaD1Codes(J+18)%len_dim1=land_pts
      UkcaD1Codes(J+18)%len_dim2=npft
      UkcaD1Codes(J+19)%item=218           ! Canopy heights of PFTs
      UkcaD1Codes(J+19)%len_dim1=land_pts
      UkcaD1Codes(J+19)%len_dim2=npft
      UkcaD1Codes(J+20)%item=229           ! Canopy water content on tiles
      UkcaD1Codes(J+20)%len_dim1=land_pts
      UkcaD1Codes(J+20)%len_dim2=ntype
      UkcaD1Codes(J+21)%item=233           ! Surface temperature on tiles
      UkcaD1Codes(J+21)%len_dim1=land_pts
      UkcaD1Codes(J+21)%len_dim2=ntype
      UkcaD1Codes(J+22)%item=234           ! Surface roughness lengths on t
      UkcaD1Codes(J+22)%len_dim1=land_pts
      UkcaD1Codes(J+22)%len_dim2=ntype
      UkcaD1Codes(J+23)%item=240           ! Snow depth on tiles
      UkcaD1Codes(J+23)%len_dim1=land_pts
      UkcaD1Codes(J+23)%len_dim2=ntiles
      UkcaD1Codes(J+24)%item=253           ! Rho_r2
      UkcaD1Codes(J+24)%len_dim1=row_length
      UkcaD1Codes(J+24)%len_dim2=rows
      UkcaD1Codes(J+24)%len_dim3=model_levels
      UkcaD1Codes(J+25)%item=254           ! QCL
      UkcaD1Codes(J+25)%len_dim1=row_length
      UkcaD1Codes(J+25)%len_dim2=rows
      UkcaD1Codes(J+25)%len_dim3=wet_levels
      UkcaD1Codes(J+26)%item=255           ! Exner pressure on rho levels
      UkcaD1Codes(J+26)%len_dim1=row_length
      UkcaD1Codes(J+26)%len_dim2=rows
      UkcaD1Codes(J+26)%len_dim3=model_levels+1
      UkcaD1Codes(J+27)%item=265           ! Area cloud fraction in each la
      UkcaD1Codes(J+27)%len_dim1=row_length
      UkcaD1Codes(J+27)%len_dim2=rows
      UkcaD1Codes(J+27)%len_dim3=wet_levels
      UkcaD1Codes(J+28)%item=266           ! Bulk Cloud fraction
      UkcaD1Codes(J+28)%len_dim1=row_length
      UkcaD1Codes(J+28)%len_dim2=rows
      UkcaD1Codes(J+28)%len_dim3=wet_levels
      IF (.NOT. L_ukca_aerchem) UkcaD1Codes(J+29)%required=.false.
      UkcaD1Codes(J+29)%item=418           ! Soil Clay
      UkcaD1Codes(J+29)%len_dim1=row_length
      UkcaD1Codes(J+29)%len_dim2=rows
      UkcaD1Codes(J+30)%item=421           ! Dust soil mass fraction for Di
      UkcaD1Codes(J+31)%item=422           ! Dust soil mass fraction for Di
      UkcaD1Codes(J+32)%item=423           ! Dust soil mass fraction for Di
      UkcaD1Codes(J+33)%item=424           ! Dust soil mass fraction for Di
      UkcaD1Codes(J+34)%item=425           ! Dust soil mass fraction for Di
      UkcaD1Codes(J+35)%item=426           ! Dust soil mass fraction for Di
      UkcaD1Codes(J+30:J+35)%len_dim1=row_length
      UkcaD1Codes(J+30:J+35)%len_dim2=rows
      IF (.NOT. L_ukca_dust) THEN   ! Required only when dust is ON
       UkcaD1Codes(J+29:J+35)%required = .false.
      ENDIF
      UkcaD1Codes(J+36)%item=505           ! Land fraction
      UkcaD1Codes(J+36)%len_dim1=land_pts
      UkcaD1Codes(J+37)%item=510           ! Land albedo
      UkcaD1Codes(J+37)%len_dim1=row_length
      UkcaD1Codes(J+37)%len_dim2=rows
      IF (.NOT. L_SW_Radiate) THEN   ! Required only on radiation TS
        UkcaD1Codes(J+37)%required = .false.
      ENDIF
      UkcaD1Codes(J+38)%item=338           ! (HadGAM1-ERA) temperature bias
      UkcaD1Codes(J+38)%len_dim1=row_length
      UkcaD1Codes(J+38)%len_dim2=rows
      UkcaD1Codes(J+38)%len_dim3=tr_levels
      UkcaD1Codes(J+39)%item=339           ! (HadGAM1-ERA) sp. humidity bias
      UkcaD1Codes(J+39)%len_dim1=row_length
      UkcaD1Codes(J+39)%len_dim2=rows
      UkcaD1Codes(J+39)%len_dim3=tr_levels
      IF (.NOT. L_ukca_Tbias) THEN
        UkcaD1Codes(J+38)%required = .false.
      ENDIF
      IF (.NOT. L_ukca_Qbias) THEN
        UkcaD1Codes(J+39)%required = .false.
      ENDIF

! Diagnostics from section zero
      J = J + n_in_progs

      UkcaD1Codes(J+1:J+n_in_diags0)%section    = 0
      UkcaD1Codes(J+1:J+n_in_diags0)%prognostic = .false.
      UkcaD1Codes(J+1:J+n_in_diags0)%required   = .true.
      UkcaD1Codes(J+1)%item=406           ! Exner Press on theta levels
      UkcaD1Codes(J+1)%len_dim1=row_length
      UkcaD1Codes(J+1)%len_dim2=rows
      UkcaD1Codes(J+1)%len_dim3=model_levels
      UkcaD1Codes(J+2)%item=407           ! P on Rho Levels
      UkcaD1Codes(J+2)%len_dim1=row_length
      UkcaD1Codes(J+2)%len_dim2=rows
      UkcaD1Codes(J+2)%len_dim3=model_levels
      UkcaD1Codes(J+3)%item=408           ! P on Theta Levels
      UkcaD1Codes(J+3)%len_dim1=row_length
      UkcaD1Codes(J+3)%len_dim2=rows
      UkcaD1Codes(J+3)%len_dim3=model_levels
      UkcaD1Codes(J+4)%item=409           ! Surface Pressure
      UkcaD1Codes(J+4)%len_dim1=row_length
      UkcaD1Codes(J+4)%len_dim2=rows

! Diagnostic items from section 1 (SW radiation)
      J = J + n_in_diags0

      UkcaD1Codes(J+1:J+n_in_diags1)%section    = 1
      UkcaD1Codes(J+1:J+n_in_diags1)%prognostic = .false.   ! always needed
      UkcaD1Codes(J+1:J+n_in_diags1)%required   = .true.
      UkcaD1Codes(J+1)%item=201           ! Net downward surface SW flux
      UkcaD1Codes(J+1)%len_dim1=row_length
      UkcaD1Codes(J+1)%len_dim2=rows
      UkcaD1Codes(J+2)%item=235           ! Total downward surface SW flux
      UkcaD1Codes(J+2)%len_dim1=row_length
      UkcaD1Codes(J+2)%len_dim2=rows

! Diagnostic variables in section 3 (Boundary Layer)
      J = J + n_in_diags1

      UkcaD1Codes(J+1:J+n_in_diags3)%section    = 3
      UkcaD1Codes(J+1:J+n_in_diags3)%prognostic = .false.
      UkcaD1Codes(J+1:J+n_in_diags3)%required   = .true.
      UkcaD1Codes(J+1)%item=217           ! Surface heat flux
      UkcaD1Codes(J+1)%len_dim1=row_length
      UkcaD1Codes(J+1)%len_dim2=rows
      UkcaD1Codes(J+2)%item=462            ! Stomatal conductance
      UkcaD1Codes(J+2)%len_dim1=row_length
      UkcaD1Codes(J+2)%len_dim2=rows
      Ukcad1codes(J+2)%len_dim3=npft
      UkcaD1Codes(J+3)%item=465            ! Surface friction velocity
      UkcaD1Codes(J+3)%len_dim1=row_length
      UkcaD1Codes(J+3)%len_dim2=rows
      UkcaD1Codes(J+4)%item=60            ! rhokh_mix
      UkcaD1Codes(J+4)%len_dim1=row_length
      UkcaD1Codes(J+4)%len_dim2=rows
      Ukcad1codes(J+4)%len_dim3=bl_levels
      UkcaD1Codes(J+5)%item=61            ! rho_aresist
      UkcaD1Codes(J+5)%len_dim1=row_length
      UkcaD1Codes(J+5)%len_dim2=rows
      UkcaD1Codes(J+6)%item=62            ! aresist
      UkcaD1Codes(J+6)%len_dim1=row_length
      UkcaD1Codes(J+6)%len_dim2=rows
      UkcaD1Codes(J+7)%item=63            ! resist_b
      UkcaD1Codes(J+7)%len_dim1=row_length
      UkcaD1Codes(J+7)%len_dim2=rows
      UkcaD1Codes(J+8)%item=64            ! dtrdz_charney_grid
      UkcaD1Codes(J+8)%len_dim1=row_length
      UkcaD1Codes(J+8)%len_dim2=rows
      Ukcad1codes(J+8)%len_dim3=bl_levels
      UkcaD1Codes(J+9)%item=65            ! kent
      UkcaD1Codes(J+9)%len_dim1=row_length
      UkcaD1Codes(J+9)%len_dim2=rows
      UkcaD1Codes(J+10)%item=66            ! we_lim
      UkcaD1Codes(J+10)%len_dim1=row_length
      UkcaD1Codes(J+10)%len_dim2=rows
      Ukcad1codes(J+10)%len_dim3=npft
      UkcaD1Codes(J+11)%item=67            ! t_frac
      UkcaD1Codes(J+11)%len_dim1=row_length
      UkcaD1Codes(J+11)%len_dim2=rows
      Ukcad1codes(J+11)%len_dim3=npft
      UkcaD1Codes(J+12)%item=68            ! zrzi
      UkcaD1Codes(J+12)%len_dim1=row_length
      UkcaD1Codes(J+12)%len_dim2=rows
      Ukcad1codes(J+12)%len_dim3=npft
      UkcaD1Codes(J+13)%item=69            ! kent_dsc
      UkcaD1Codes(J+13)%len_dim1=row_length
      UkcaD1Codes(J+13)%len_dim2=rows
      UkcaD1Codes(J+14)%item=70            ! we_lim_dsc
      UkcaD1Codes(J+14)%len_dim1=row_length
      UkcaD1Codes(J+14)%len_dim2=rows
      Ukcad1codes(J+14)%len_dim3=npft
      UkcaD1Codes(J+15)%item=71            ! t_frac_dsc
      UkcaD1Codes(J+15)%len_dim1=row_length
      UkcaD1Codes(J+15)%len_dim2=rows
      Ukcad1codes(J+15)%len_dim3=npft
      UkcaD1Codes(J+16)%item=72            ! zrzi_dsc
      UkcaD1Codes(J+16)%len_dim1=row_length
      UkcaD1Codes(J+16)%len_dim2=rows
      Ukcad1codes(J+16)%len_dim3=npft
      UkcaD1Codes(J+17)%item=73            ! zhsc
      UkcaD1Codes(J+17)%len_dim1=row_length
      UkcaD1Codes(J+17)%len_dim2=rows
      UkcaD1Codes(J+18)%item=74            ! r_b_dust div 1
      UkcaD1Codes(J+18)%len_dim1=row_length
      UkcaD1Codes(J+18)%len_dim2=rows
      UkcaD1Codes(J+19)%item=75            ! r_b_dust div 2
      UkcaD1Codes(J+19)%len_dim1=row_length
      UkcaD1Codes(J+19)%len_dim2=rows
      UkcaD1Codes(J+20)%item=76            ! r_b_dust div 3
      UkcaD1Codes(J+20)%len_dim1=row_length
      UkcaD1Codes(J+20)%len_dim2=rows
      UkcaD1Codes(J+21)%item=77            ! r_b_dust div 4
      UkcaD1Codes(J+21)%len_dim1=row_length
      UkcaD1Codes(J+21)%len_dim2=rows
      UkcaD1Codes(J+22)%item=78            ! r_b_dust div 5
      UkcaD1Codes(J+22)%len_dim1=row_length
      UkcaD1Codes(J+22)%len_dim2=rows
      UkcaD1Codes(J+23)%item=79            ! r_b_dust div 6
      UkcaD1Codes(J+23)%len_dim1=row_length
      UkcaD1Codes(J+23)%len_dim2=rows
      UkcaD1Codes(J+24)%item=25            ! ml_depth
      UkcaD1Codes(J+24)%len_dim1=row_length
      UkcaD1Codes(J+24)%len_dim2=rows
      UkcaD1Codes(J+25)%item=209           ! 10 m wind U component
      UkcaD1Codes(J+25)%len_dim1=row_length
      UkcaD1Codes(J+25)%len_dim2=rows
      IF (.NOT. L_UKCA_MODE) UkcaD1Codes(J+25)%required=.false.
      UkcaD1Codes(J+26)%item=210           ! 10 m wind V component
      UkcaD1Codes(J+26)%len_dim1=row_length
      UkcaD1Codes(J+26)%len_dim2=n_rows
      IF (.NOT. L_UKCA_MODE) UkcaD1Codes(J+26)%required=.false.

! Diagnostic variables in section 4 (LS Precipitation)
      J = J + n_in_diags3

      UkcaD1Codes(J+1:J+n_in_diags4)%section    = 4
      UkcaD1Codes(J+1:J+n_in_diags4)%prognostic = .false.
      UkcaD1Codes(J+1:J+n_in_diags4)%required   = .true.
      UkcaD1Codes(J+1)%item=205           ! Cloud Liquid Water after LS
      UkcaD1Codes(J+1)%len_dim1=row_length
      UkcaD1Codes(J+1)%len_dim2=rows
      UkcaD1Codes(J+1)%len_dim3=wet_levels
      UkcaD1Codes(J+2)%item=206           ! Cloud Ice Content after LS
      UkcaD1Codes(J+2)%len_dim1=row_length
      UkcaD1Codes(J+2)%len_dim2=rows
      UkcaD1Codes(J+2)%len_dim3=wet_levels
      UkcaD1Codes(J+3)%item=222           ! Rainfall rate out of model levs
      UkcaD1Codes(J+3)%len_dim1=row_length
      UkcaD1Codes(J+3)%len_dim2=rows
      UkcaD1Codes(J+3)%len_dim3=wet_levels
      UkcaD1Codes(J+4)%item=223           ! Snowfall rate out of model levs
      UkcaD1Codes(J+4)%len_dim1=row_length
      UkcaD1Codes(J+4)%len_dim2=rows
      UkcaD1Codes(J+4)%len_dim3=wet_levels
      UkcaD1Codes(J+5)%item=227           ! Rain fraction (3C ppn only)
      UkcaD1Codes(J+5)%len_dim1=row_length
      UkcaD1Codes(J+5)%len_dim2=rows
      UkcaD1Codes(J+5)%len_dim3=wet_levels

! Diagnostic variables in section 5 (Convection)
      J = J + n_in_diags4

      UkcaD1Codes(J+1:J+n_in_diags5)%section    = 5
      UkcaD1Codes(J+1:J+n_in_diags5)%prognostic = .false.
      UkcaD1Codes(J+1:J+n_in_diags5)%required   = .true.
      UkcaD1Codes(J+1)%item=227           ! 3D Convective rainfall rate
      UkcaD1Codes(J+1)%len_dim1=row_length
      UkcaD1Codes(J+1)%len_dim2=rows
      UkcaD1Codes(J+1)%len_dim3=model_levels
      UkcaD1Codes(J+2)%item=228           ! 3D Convective snowfall rate
      UkcaD1Codes(J+2)%len_dim1=row_length
      UkcaD1Codes(J+2)%len_dim2=rows
      UkcaD1Codes(J+2)%len_dim3=model_levels

! Diagnostic variables in section 8 (Hydrology)
      J = J + n_in_diags5

      UkcaD1Codes(J+1)%section    = 8
      UkcaD1Codes(J+1)%prognostic = .false.
      UkcaD1Codes(J+1)%required   = L_ukca_qch4inter
      UkcaD1Codes(J+1)%item       = 242        ! CH4 wetland flux
      UkcaD1Codes(J+1)%len_dim1   = row_length
      UkcaD1Codes(J+1)%len_dim2   = rows

! Diagnostic variables in section 15 (Processed Climate Diagnostics)
      J = J + n_in_diags8

      UkcaD1Codes(J+1)%section    = 15
      UkcaD1Codes(J+1)%prognostic = .false.
      UkcaD1Codes(J+1)%required   = .true.
      UkcaD1Codes(J+1)%item       = 218           ! PV on theta levels
      UkcaD1Codes(J+1)%len_dim1   = row_length
      UkcaD1Codes(J+1)%len_dim2   = rows
      UkcaD1Codes(J+1)%len_dim3   = model_levels

! Diagnostic variables in section 30 (Processed dynamics diagnostics)
      J = J + n_in_diags15

! Take tropopause height here. Needed only for volcanic SO2 emissions into
! the stratosphere. (Always required as in call to emission_ctl)
      UkcaD1Codes(J+1)%section = 30
      UkcaD1Codes(J+1)%prognostic = .false.
      UkcaD1Codes(J+1)%required   = .true.
      UkcaD1Codes(J+1)%item=453           ! Tropopause height
      UkcaD1Codes(J+1)%len_dim1=row_length
      UkcaD1Codes(J+1)%len_dim2=rows


! Diagnostic variables in section 34 (UKCA)
      J = J + n_in_diags30
      idiag_first = J+1
      idiag_last  = J+n_chem_diags

      DO I=1,n_chem_diags
        UkcaD1Codes(J+I)%section    = UKCA_sect
        UkcaD1Codes(J+I)%item       = 251+I-1
        UkcaD1Codes(J+I)%len_dim1   = row_length
        UkcaD1Codes(J+I)%len_dim2   = rows
        UkcaD1Codes(J+I)%len_dim3   = model_levels
        UkcaD1Codes(J+I)%required   = .true.
        UkcaD1Codes(J+I)%prognostic = .true.
      ENDDO

! Diagnostic variables for Backward-Euler fluxes in section 34
      J = J + n_chem_diags
      iflux_first = J+1
      iflux_last  = J+n_BE_fluxdiags

      IF (L_ukca_BEflux) THEN
        DO I=1,n_BE_fluxdiags
          UkcaD1Codes(J+I)%section    = UKCA_sect
          UkcaD1Codes(J+I)%item       = 301+I-1
          UkcaD1Codes(J+I)%len_dim1   = row_length
          UkcaD1Codes(J+I)%len_dim2   = rows
          UkcaD1Codes(J+I)%len_dim3   = model_levels
          UkcaD1Codes(J+I)%required   = .true.
          UkcaD1Codes(J+I)%prognostic = .false.
        ENDDO
      ENDIF

! Diagnostic variables for stratospheric fluxes in section 34
      J = J + n_BE_fluxdiags
      istrat_first = J+1
      istrat_last  = J+nmax_strat_fluxdiags

      n_strat_fluxdiags = 0
      DO I=1,nmax_strat_fluxdiags
        DO idiag=1,ndiag     ! search for stash requests
          IF (modl_b(idiag) == SUBMODEL_FOR_SM(atmos_im) .AND.         &
              isec_b(idiag) == UKCA_sect .AND.                         &
              item_b(idiag) == 471+I-1) THEN
            n_strat_fluxdiags = n_strat_fluxdiags + 1
            UkcaD1Codes(J+I)%section    = UKCA_sect
            UkcaD1Codes(J+I)%item       = item_b(idiag)
            UkcaD1Codes(J+I)%len_dim1   = row_length
            UkcaD1Codes(J+I)%len_dim2   = rows
            UkcaD1Codes(J+I)%len_dim3   = model_levels
            UkcaD1Codes(J+I)%required   = .true.
            UkcaD1Codes(J+I)%prognostic = .false.
           IF (.NOT.(L_ukca_stratflux)) L_ukca_stratflux = .true.
            EXIT
          ENDIF
        ENDDO
      ENDDO

      WRITE(6,*) 'n_use_tracers: ',n_use_tracers
      WRITE(6,*) 'n_use_emissions: ',n_use_emissions
      WRITE(6,*) 'Total no of items required = ',Nukca_D1items
      WRITE(6,*) ' UKCA: The following items are specified in D1:'
      WRITE(6,*) 'section, item, prognostic, required, dim1, dim2, dim3'
      DO I=1,Nukca_D1items
        WRITE(6,*) I,UkcaD1Codes(I)%section,                           &
                     UkcaD1Codes(I)%item,                              &
                     UkcaD1Codes(I)%prognostic,                        &
                     UkcaD1Codes(I)%required,                          &
                     UkcaD1Codes(I)%len_dim1,                          &
                     UkcaD1Codes(I)%len_dim2,                          &
                     UkcaD1Codes(I)%len_dim3
      ENDDO

      END SUBROUTINE UKCA_SETD1DEFS
#endif
