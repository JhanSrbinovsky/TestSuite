#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL--------------- SUBROUTINE INITPHYS --------------------------------
!LL
!LL Purpose : Calls the initialisation program for the data assimilation
!LL assimilation section (P3)
!LL  Control routine for CRAY YMP
!LL
!LL T.Johns     <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL 3.2    27/03/93 Change INITPHYS for dynamic allocation. Remove
!LL                 size parameters no longer needed.       R.Rawlins.
!   4.0    06/07/95 CNTLATM called to allow access to version numbers
!                   to pick up new routines for reading spectral
!                   namelists for version 3A of the SW or LW radiation.
!                                             J. M. Edwards
!   4.4    02/09/97 Error checking for different sections improved.
!                   Logical flags for aerosols passed down to allow
!                   earlier checking of the spectral data.
!                                             J. M. Edwards
!   4.5    18/05/98 Logicals for including gases are now passed down
!                   to lower levels to select gases from the
!                   spectral file.
!                                             J. M. Edwards
!  4.5  April 1998  Pass soot logical to R2_SW_SPECIN. Luke Robinson.
!LL 5.0    30/06/99 Add call to CRUNTIMC to get L_CLIMAT_AEROSOL, which
!LL                 used to be in cntlatm at 4.5. M.L.Gallani
!   5.2    15/11/00 Add flags for direct effect of sea-salt aerosol
!                   to calls to R2_SW_SPECIN and R2_LW_SPECIN.
!                                             A. Jones
!   5.3    04/04/01 Pass mesoscale aerosol logical L_MURK_RAD
!                   to lower levels.                      S. Cusack
!   5.5    05/02/03 Add logical for direct effect of biomass aerosol
!                   to calls to R2_SW_SPECIN and R2_LW_SPECIN.
!                                                P Davison
!   5.5    21/02/03 Add logical for direct effect of mineral dust
!                   to calls to R2_SW_SPECIN and R2_LW_SPECIN.
!                                                S Woodward
!   6.2    15/11/05 Add logical saying whether the aerosol optical
!                   depth block is needed to call to R2_LW_SPECIN.
!                                                N Bellouin
! 6.2  25/11/05 Functionality for improved time stepping, radiative
!               forcing and radiance code added for versions 3C
!               and 3Z of radiation code             (J.-C. Thelen)
!   6.2    21/03/06 Included nstypes.h J Ridley.
!LL
!LL Programming standard; U M Documentation Paper No. 3 version 1
!LL dated 15/01/90
!LL
!LL System components covered P0
!LL
!LL Documentation : UM documentation paper  no P0
!LL
!LL END


      SUBROUTINE INITPHYS(ICODE,CMESSAGE)

!
! If new version of radiation code is required then we
! need to use the modules for improved time-stepping
!
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)

      USE CONTROL_STRUC
      USE SW_CONTROL_STRUCT
      USE LW_CONTROL_STRUCT
      USE SPEC_SW_LW

#endif
!
      IMPLICIT NONE

      INTEGER    ICODE         ! Return code : 0 Normal exit
!                              !             : >0 Error
      CHARACTER*80 CMESSAGE    ! Error message if ICODE > 0

#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
!     Local variables
      INTEGER :: J
#endif
!*
!*L Subroutines called
      EXTERNAL R2_SW_SPECIN, R2_LW_SPECIN
!*
#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "cconsts.h"
#include "cntlatm.h"
#if defined(A01_3A) || defined(A02_3A)
#include "swopt3a.h"
#include "swcopt3a.h"
#include "lwopt3a.h"
#include "lwcopt3a.h"
#endif

!
#if defined(A01_3A) || defined(A02_3A)
!
!     ------------- Shortwave Radiation -------------------------
!
      IF (H_SECT(1) == '03A') THEN
!
! DEPENDS ON: r2_sw_specin
        CALL R2_SW_SPECIN(ICODE, CMESSAGE                               &
     &    , L_O2_SW                                                     &
     &    , L_CLIMAT_AEROSOL                                            &
     &    , L_USE_DUST, L_USE_ARCLDUST                                  &
     &    , L_USE_SULPC_DIRECT, L_USE_ARCLSULP                          &
     &    , L_USE_SOOT_DIRECT, L_USE_ARCLBLCK                           &
     &    , L_USE_BMASS_DIRECT, L_USE_ARCLBIOM                          &
     &    , L_USE_SEASALT_DIRECT, L_USE_ARCLSSLT                        &
     &    , L_USE_OCFF_DIRECT, L_USE_ARCLOCFF                           &
     &    , L_USE_BIOGENIC, L_USE_ARCLDLTA                              &
     &    , L_MURK_RAD                                                  &
     &   )
        IF (ICODE /= 0) RETURN
!
      ELSE
!
        ICODE=1
        CMESSAGE='Unknown version of SW radiation encountered.'
        RETURN
!
      ENDIF
!
!
!     ------------- Longwave Radiation -------------------------
!
      IF (H_SECT(2) == '03A') THEN
!
! DEPENDS ON: r2_lw_specin
        CALL R2_LW_SPECIN(ICODE, CMESSAGE                               &
     &    , L_CH4_LW, L_N2O_LW, L_CFC11_LW, L_CFC12_LW                  &
     &    , L_CFC113_LW, L_HCFC22_LW, L_HFC125_LW, L_HFC134A_LW         &
     &    , L_CLIMAT_AEROSOL                                            &
     &    , L_USE_DUST, L_USE_ARCLDUST                                  &
     &    , L_USE_SULPC_DIRECT, L_USE_ARCLSULP                          &
     &    , L_USE_SOOT_DIRECT, L_USE_ARCLBLCK                           &
     &    , L_USE_BMASS_DIRECT, L_USE_ARCLBIOM                          &
     &    , L_USE_SEASALT_DIRECT, L_USE_ARCLSSLT                        &
     &    , L_USE_OCFF_DIRECT, L_USE_ARCLOCFF                           &
     &    , L_USE_BIOGENIC, L_USE_ARCLDLTA                              &
     &    , L_MURK_RAD                                                  &
     &    , L_USE_AOD                                                   &
     &   )
!
      ELSE
!
        ICODE=1
        CMESSAGE='Unknown version of LW radiation encountered.'
        RETURN
!
      ENDIF
!
#endif
!
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
!
!------------- Shortwave Radiation -------------------------
!
      IF ((H_SECT(1) == '03C').OR.(H_SECT(1) == '03Z')) THEN
!
         DO J=1, N_SWCALL
! DEPENDS ON: r2_sw_specin
            CALL R2_SW_SPECIN(ICODE, CMESSAGE                           &
     &        , SW_CONTROL(J)%SPECTRAL_FILE                             &
     &        , SW_CONTROL(J)%L_O2                                      &
     &        , L_CLIMAT_AEROSOL                                        &
     &        , L_USE_DUST, L_USE_ARCLDUST                              &
     &        , L_USE_SULPC_DIRECT, L_USE_ARCLSULP                      &
     &        , L_USE_SOOT_DIRECT, L_USE_ARCLBLCK                       &
     &        , L_USE_BMASS_DIRECT, L_USE_ARCLBIOM                      &
     &        , L_USE_SEASALT_DIRECT, L_USE_ARCLSSLT                    &
     &        , L_USE_OCFF_DIRECT, L_USE_ARCLOCFF                       &
     &        , L_USE_BIOGENIC, L_USE_ARCLDLTA                          &
     &        , L_MURK_RAD                                              &
     &        , SW_CONTROL(J)%L_RAYLEIGH, SW_CONTROL(J)%L_GAS           &
     &        , SW_CONTROL(J)%L_CONTINUUM, SW_CONTROL(J)%L_DROP         &
     &        , SW_CONTROL(J)%L_AEROSOL, SW_CONTROL(J)%L_ICE            &
     &        , SW_SPECTRUM(J)                                          &
     &        )
         ENDDO

         IF (ICODE /= 0) RETURN
      ELSE
         ICODE=1
         CMESSAGE='Unknown version of SW radiation encountered.'
         RETURN
      ENDIF
!
!     ------------- Longwave Radiation -------------------------
!
      IF ((H_SECT(2) == '03C').OR.(H_SECT(2) == '03Z')) THEN

         DO J=1,N_LWCALL
! DEPENDS ON: r2_lw_specin
            CALL R2_LW_SPECIN(ICODE, CMESSAGE                           &
     &        , LW_CONTROL(J)%SPECTRAL_FILE                             &
     &        , LW_CONTROL(J)%L_CH4, LW_CONTROL(J)%L_N2O                &
     &        , LW_CONTROL(J)%L_CFC11, LW_CONTROL(J)%L_CFC12            &
     &        , LW_CONTROL(J)%L_CFC113                                  &
     &        , LW_CONTROL(J)%L_HCFC22, LW_CONTROL(J)%L_HFC125          &
     &        , LW_CONTROL(J)%L_HFC134A                                 &
     &        , L_CLIMAT_AEROSOL                                        &
     &        , L_USE_DUST, L_USE_ARCLDUST                              &
     &        , L_USE_SULPC_DIRECT, L_USE_ARCLSULP                      &
     &        , L_USE_SOOT_DIRECT, L_USE_ARCLBLCK                       &
     &        , L_USE_BMASS_DIRECT, L_USE_ARCLBIOM                      &
     &        , L_USE_SEASALT_DIRECT, L_USE_ARCLSSLT                    &
     &        , L_USE_OCFF_DIRECT, L_USE_ARCLOCFF                       &
     &        , L_USE_BIOGENIC, L_USE_ARCLDLTA                          &
     &        , L_MURK_RAD                                              &
     &        , L_USE_AOD                                               &
     &        , LW_CONTROL(J)%L_GAS, LW_CONTROL(J)%L_CONTINUUM          &
     &        , LW_CONTROL(J)%L_DROP, LW_CONTROL(J)%L_AEROSOL           &
     &        , LW_CONTROL(J)%L_ICE                                     &
     &        , LW_SPECTRUM(J)                                          &
     &        )
          ENDDO
          IF (ICODE /= 0) RETURN
      ELSE
         ICODE=1
         CMESSAGE='Unknown version of LW radiation encountered.'
         RETURN
      ENDIF
#endif
!
!
      RETURN
      END SUBROUTINE INITPHYS
#endif
