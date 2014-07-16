#if defined(A70_1A) || defined(A70_1B)
#if defined(A01_3A) || defined(A02_3A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to read a shortwave spectral namelist.
!
! Purpose:
!   To read a shortwave namelist into a spectral array.
!
! Method:
!   The spectrum is read into the dynamically allocated array
!   and then reduced to a more manageable size.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             14-05-96                Set lower limits
!                                               for reduced dimensions
!                                               to ensure that they
!                                               may never be 0.
!                                               (J. M. Edwards)
!       4.4             02-09-97                Aerosol flags passed
!                                               in to the code to
!                                               enable only those
!                                               required to be
!                                               selected. Spectral
!                                               data are now longer
!                                               compressed into a
!                                               single array.
!                                               Actual IOS code put
!                                               into CMESSAGE.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Coding to allow
!                                               selection of gases
!                                               from the spectral
!                                               file.
!                                               (J. M. Edwards)
!
!       4.5        April 1998   Allow soot spectral data to be read.
!                                                     Luke Robinson.
!       5.2             15-11-00                Coding to allow
!                                               sea-salt aerosol
!                                               to be read from
!                                               the spectral file
!                                               if required.
!                                               (A. Jones)
!       5.3             04-04-01                Allow aerosol data to be
!                                               read in if required by
!                                               mesoscale model.
!                                               (S. Cusack)
!       5.5             05-02-03                Allow biomass aerosol
!                                               spectral data to be
!                                               read in, if required.
!                                               (P. Davison)
!       5.5             21-02-03                Allow mineral dust
!                                               spectral data to be
!                                               read in, if required.
!                                               (S Woodward)
!       6.2             16-11-05                Argument list of
!                                               R2_COMPRESS_SPECTRUM
!                                               modified to match
!                                               changes in the LW.
!                                               (N Bellouin)
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SW_SPECIN(IERR, CMESSAGE                            &
     &   , L_O2                                                         &
     &   , L_CLIMAT_AEROSOL                                             &
     &   , L_USE_DUST, L_USE_ARCLDUST                                   &
     &   , L_USE_SULPC_DIRECT, L_USE_ARCLSULP                           &
     &   , L_USE_SOOT_DIRECT, L_USE_ARCLBLCK                            &
     &   , L_USE_BMASS_DIRECT, L_USE_ARCLBIOM                           &
     &   , L_USE_SEASALT_DIRECT, L_USE_ARCLSSLT                         &
     &   , L_USE_OCFF_DIRECT, L_USE_ARCLOCFF                            &
     &   , L_USE_BIOGENIC, L_USE_ARCLDLTA                               &
     &   , L_MURK_RAD                                                   &
     &   )
!
!
      IMPLICIT NONE
!
!
#include "mxsize3a.h"
#include "error3a.h"
#include "stdio3a.h"
!
!
!     DUMMY ARGUMENTS
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_O2                                                         &
!             ABSORPTION BY OXYGEN IS TO BE INCLUDED.
     &   , L_CLIMAT_AEROSOL                                             &
!             CLIMATOLOGICAL AEROSOLS ARE TO BE INCLUDED
     &   , L_USE_DUST                                                   &
                      ! direct effects of mineral dust aerosol
     &   , L_USE_ARCLDUST                                               &
                 ! direct effect of mineral dust aerosol from NWP clim
     &   , L_USE_SULPC_DIRECT                                           &
!             THE DIRECT EFFECTS OF SULPHATE AEROSOLS ARE
!             TO BE INCLUDED
     &   , L_USE_ARCLSULP                                               &
                 ! direct effect of sulphate aerosol from NWP clim
     &   , L_USE_SOOT_DIRECT                                            &
!             USE THE DIRECT RAD EFFECTS OF SOOT IN THE SW
     &   , L_USE_ARCLBLCK                                               &
                 ! direct effect of black-carbon aerosol from NWP clim
     &   , L_USE_BMASS_DIRECT                                           &
!             USE THE DIRECT RAD EFFECTS OF BIOMASS SMOKE IN THE SW
     &   , L_USE_ARCLBIOM                                               &
                 ! direct effect of biomass aerosol from NWP clim
     &   , L_USE_SEASALT_DIRECT                                         &
!             USE THE DIRECT RAD EFFECTS OF SEASALT IN THE SW
     &   , L_USE_ARCLSSLT                                               &
                 ! direct effect of sea-salt aerosol from NWP clim
     &   , L_USE_BIOGENIC                                               &
!             USE THE BIOGENIC AEROSOL DIRECT EFFECT IN THE SW
     &   , L_USE_OCFF_DIRECT                                            &
                 ! direct effect of fossil-fuel org.carb. in the SW
     &   , L_USE_ARCLOCFF                                               &
                 ! direct effect of fossil-fuel org.carb. from NWP clim
     &   , L_USE_ARCLDLTA                                               &
                 ! direct effect of delta aerosol from NWP climatology
     &   , L_MURK_RAD
!             MESOSCALE MODEL AEROSOLS ARE TO BE INCLUDED
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      CHARACTER*80                                                      &
                        !, INTENT(OUT)
     &     CMESSAGE
!
!
!
!     LOCAL VARIABLES.
!
!
!     RADIATIVE VARIABLES FOR REDUCING THE SPECTRUM
!
#include "aerprm3a.h"
#include "aercmp3a.h"
#include "gasid3a.h"
!
      CHARACTER*256                                                     &
     &     SW_SPECTRAL_FILE
!             NAME OF FILE CONTAINING THE SPECTRAL DATA
      INTEGER                                                           &
     &     IERR_GET_FILE                                                &
!             ERROR FLAG RETURNED BY GET_FILE (NOT NECESSARILY
!             CONSISTENT WITH THE FLAGS IN ERROR3A).
     &   , IOS
!             STATUS OF I/O
!
      LOGICAL                                                           &
     &     L_RETAIN_ABSORB(NPD_SPECIES)                                 &
!             FLAG SET TO .TRUE. IF THE ABSORBER IS TO BE RETAINED
     &   , L_GAS_INCLUDED(NPD_GASES)
!             LOGICAL TO TEST FOR ACTUAL GASES INCLUDED
      INTEGER                                                           &
     &     N_ABSORB_RETAIN                                              &
!             NUMBER OF ABSORBERS TO RETAIN
     &   , INDEX_ABSORB_RETAIN(NPD_SPECIES)                             &
!             INDICES OF ABSORBERS TO BE RETAINED
     &   , COMPRESSED_INDEX(NPD_SPECIES)                                &
!             MAPPING FROM ORIGINAL TO COMPRESSED INDICES OF ABSORBERS
     &   , N_AEROSOL_RETAIN                                             &
!             NUMBER OF AEROSOLS IN THE SPECTRAL FILE TO BE RETAINED
!             FOR THE RADIATIVE CALCULATION
     &   , INDEX_AEROSOL_RETAIN(NPD_AEROSOL_SPECIES)                    &
!             INDEXING NUMBERS OF THE RETAINED AEROSOLS
     &   , N_AEROSOL_FOUND
!             NUMBER OF AEROSOLS FOR THE CURRENT GROUP OF PROCESSES
!             FOUND IN THE SPECTRAL FILE
!
!
!
!     DECLARE THE ELEMENTS OF THE INITIAL SPECTRUM FOR DYNAMIC
!     ALLOCATION AND SET UP AN APPROPRIATE NAMELIST.
!
#include "spdec3a.h"
#include "swsp3a.h"
!
!
!     DECLARE THE REDUCED SW SPECTRAL FILE AND ITS HOLDING COMMON BLOCK.
!
#include "swspdl3a.h"
#include "swspcm3a.h"
!
!
!
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
!
      CHARACTER                                                         &
     &     CH_IOS*5
!             CHARACTER STRING FOR IOS ERROR
!
!
!     SUBROUTINES CALLED
      EXTERNAL                                                          &
     &     R2_COMPRESS_SPECTRUM
!
!
!     EACH BLOCK IS INITIALIZED AS MISSING:
      DATA L_PRESENT/.FALSE., NPD_TYPE*.FALSE./
!
!     INITIALIZE THE RANGE OF VALIDITY OF THE PARAMETRIZATIONS OF
!     DROPLETS AND ICE CRYSTALS. OLD SPECTRAL FILES WILL NOT CONTAIN
!     SUCH DATA, SO THE LIMITS FOR DROPLETS ARE INITIALIZED TO THOSE
!     FORMERLY SET IN THE MICROPHYSICAL SCHEME (MRF/UMIST
!     PARAMETRIZATION) TO ENSURE THAT THE RESULTS ARE BIT-REPRODUCIBLE.
!     VALUES FOR ICE COVER THE RANGE OF EFFECTIVE RADII USED IN
!     GENERATING THE DATA FOR THE ORIGINAL PARAMETRIZATION OF ICE
!     CRYSTALS.
!     AT SOME FUTURE RELEASE IT MAY BE DESIRABLE TO REMOVE DEFAULT
!     SETTINGS.
      DATA DROP_PARM_MIN_DIM/NPD_DROP_TYPE*3.5E-07/
      DATA DROP_PARM_MAX_DIM/NPD_DROP_TYPE*3.7E-05/
      DATA ICE_PARM_MIN_DIM/NPD_ICE_TYPE*3.75E-07/
      DATA ICE_PARM_MAX_DIM/NPD_ICE_TYPE*8.0E-05/
!
!
!
!     READ THE SHORTWAVE SPECTRUM AS A NAMELIST.
      CALL GET_FILE(57, SW_SPECTRAL_FILE, 256, IERR_GET_FILE)
      IF (IERR_GET_FILE /= 0) THEN
!        CONVERT THE ERROR FLAG FROM GET_FILE TO A FLAG RECOGNISED
!        BY THE RADIATION CODE.
         IERR=I_ERR_IO
         CMESSAGE='Error reading name of shortwave spectral file.'
         RETURN
      ENDIF
      OPEN(UNIT=57, FILE=SW_SPECTRAL_FILE, IOSTAT=IOS)
      IF (IOS /= 0) THEN
         IERR=I_ERR_IO
      WRITE(CH_IOS, '(I5)') IOS
         CMESSAGE='Error opening shortwave spectral file.'              &
     &      //' IOSTAT='//CH_IOS
         RETURN
      ENDIF
      READ(57, R2SWSP)
      CLOSE(57)
!
!     TEST FOR MINIMAL REQUISITE INFORMATION.
      IF ( .NOT.(L_PRESENT(0).AND.                                      &
     &           L_PRESENT(2) ) ) THEN
         CMESSAGE='Shortwave spectrum is deficient.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!
!
!     SET REDUCED DIMENSIONS, EITHER FROM THE SIZES OF THE FIXED ARRAYS
!     OR FROM THE ARRAYS READ IN.
!
      NPD_TYPE_SW=NPD_TYPE
      NPD_BAND_SW=MAX(N_BAND, 1)
      NPD_SPECIES_SW=MAX(N_ABSORB, 1)
      NPD_ALBEDO_PARM_SW=NPD_ALBEDO_PARM
      NPD_SCALE_FNC_SW=NPD_SCALE_FNC
      NPD_SCALE_VARIABLE_SW=NPD_SCALE_VARIABLE
      NPD_SURFACE_SW=NPD_SURFACE
      NPD_CONTINUUM_SW=NPD_CONTINUUM
      NPD_CLOUD_PARAMETER_SW=NPD_CLOUD_PARAMETER
      NPD_THERMAL_COEFF_SW=1
      NPD_AOD_WAVEL_SW=1
!
!
!     SEARCH THE SPECTRUM TO FIND MAXIMUM DIMENSIONS.
!
      NPD_EXCLUDE_SW=1
      IF (L_PRESENT(14)) THEN
         DO I=1, N_BAND
            NPD_EXCLUDE_SW=MAX(NPD_EXCLUDE_SW, N_BAND_EXCLUDE(I))
         ENDDO
      ENDIF
!
!     Search the spectrum to find those gases to be retained.
!     Water vapour, carbon dioxide and ozone are included
!     if present, but a warning is printed if they are
!     not included.
      DO I=1, NPD_GASES
         L_GAS_INCLUDED(I)=.FALSE.
      ENDDO
      N_ABSORB_RETAIN=0
!
      DO I=1, N_ABSORB
!
         L_RETAIN_ABSORB(I)=.FALSE.
         COMPRESSED_INDEX(I)=0
!
         IF ( (TYPE_ABSORB(I) == IP_H2O).OR.                            &
     &        (TYPE_ABSORB(I) == IP_CO2).OR.                            &
     &        (TYPE_ABSORB(I) == IP_O3).OR.                             &
     &        ( (TYPE_ABSORB(I) == IP_O2).AND.L_O2 ) ) THEN
            N_ABSORB_RETAIN=N_ABSORB_RETAIN+1
            INDEX_ABSORB_RETAIN(N_ABSORB_RETAIN)=I
            COMPRESSED_INDEX(I)=N_ABSORB_RETAIN
            L_RETAIN_ABSORB(I)=.TRUE.
            L_GAS_INCLUDED(TYPE_ABSORB(I))=.TRUE.
         ENDIF
!
      ENDDO
!
!
!     Print warning messages if those gases normally expected
!     are not present.
      IF (.NOT.L_GAS_INCLUDED(IP_H2O)) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** WARNING: Water vapour is not included in the '         &
     &      , 'shortwave spectral file.'
      ENDIF
!
      IF (.NOT.L_GAS_INCLUDED(IP_CO2)) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** WARNING: Carbon dioxide is not included in the '       &
     &      , 'shortwave spectral file.'
      ENDIF
!
      IF (.NOT.L_GAS_INCLUDED(IP_O3)) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** WARNING: Ozone is not included in the '                &
     &      , 'shortwave spectral file.'
      ENDIF
!
      IF ((.NOT.L_GAS_INCLUDED(IP_O2)).AND.L_O2) THEN
         WRITE(IU_ERR, '(/A, /A)')                                      &
     &      '*** ERROR: Oxygen is not included in the shortwave '       &
     &      , 'spectral file, but was requested in the run.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!     Set an appropriate reduced dimension.
      NPD_SPECIES_SW=MAX(N_ABSORB_RETAIN, 1)
!
!
      NPD_ESFT_TERM_SW=1
      IF (L_PRESENT(5)) THEN
         DO I=1, N_BAND
            DO J=1, N_BAND_ABSORB(I)
               IF (L_RETAIN_ABSORB(INDEX_ABSORB(J, I)))                 &
     &            NPD_ESFT_TERM_SW=MAX(NPD_ESFT_TERM_SW                 &
     &            , I_BAND_ESFT(I, INDEX_ABSORB(J, I)))
            ENDDO
         ENDDO
      ENDIF
!
      NPD_DROP_TYPE_SW=1
      IF (L_PRESENT(10)) THEN
         DO I=1, NPD_DROP_TYPE
            IF (L_DROP_TYPE(I)) THEN
               NPD_DROP_TYPE_SW=MAX(NPD_DROP_TYPE_SW, I)
            ENDIF
         ENDDO
      ENDIF
!
      NPD_ICE_TYPE_SW=1
      IF (L_PRESENT(12)) THEN
         DO I=1, NPD_ICE_TYPE
            IF (L_ICE_TYPE(I)) THEN
               NPD_ICE_TYPE_SW=MAX(NPD_ICE_TYPE_SW, I)
            ENDIF
         ENDDO
      ENDIF
!
!
!     Aerosols must be treated carefully to allow for various
!     different combinations without requiring the spectral file
!     to be too constrained. Only those required will be retained.
!
!     Basic initialization to safe values.
      NPD_HUMIDITIES_SW=1
      N_AEROSOL_RETAIN=0
!
!     Check the spectral file for climatological aerosols
      IF (L_CLIMAT_AEROSOL.OR.L_MURK_RAD) THEN
!
         IF (L_PRESENT(11)) THEN
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I) == IP_WATER_SOLUBLE).OR.           &
     &              (TYPE_AEROSOL(I) == IP_DUST_LIKE).OR.               &
     &              (TYPE_AEROSOL(I) == IP_OCEANIC).OR.                 &
     &              (TYPE_AEROSOL(I) == IP_SOOT).OR.                    &
     &              (TYPE_AEROSOL(I) == IP_SULPHURIC) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF

            ENDDO
!
            IF (N_AEROSOL_FOUND /= 5) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'climatological aerosols.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF

!     Check the spectral file for mineral dust classes.
!     (Only required for the direct effect).
!
      IF (L_USE_DUST .OR. L_USE_ARCLDUST) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I) == IP_DUST_1).OR.                  &
     &              (TYPE_AEROSOL(I) == IP_DUST_2).OR.                  &
     &              (TYPE_AEROSOL(I) == IP_DUST_3).OR.                  &
     &              (TYPE_AEROSOL(I) == IP_DUST_4).OR.                  &
     &              (TYPE_AEROSOL(I) == IP_DUST_5).OR.                  &
     &              (TYPE_AEROSOL(I) == IP_DUST_6) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 6) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'mineral dust aerosol.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!
!     Check the spectral file for sulphate aerosols. (These are
!     required only for the direct effect).
!
      IF (L_USE_SULPC_DIRECT .OR. L_USE_ARCLSULP) THEN
!
         IF (L_PRESENT(11)) THEN
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I) == IP_ACCUM_SULPHATE).OR.          &
     &              (TYPE_AEROSOL(I) == IP_AITKEN_SULPHATE) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF

            ENDDO
!
            IF (N_AEROSOL_FOUND /= 2) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'sulphate aerosols.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!
!     Check the spectral file for soot aerosol modes. (Also only
!     required for the direct effect).
!
      IF (L_USE_SOOT_DIRECT .OR. L_USE_ARCLBLCK) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I) == IP_FRESH_SOOT).OR.              &
     &              (TYPE_AEROSOL(I) == IP_AGED_SOOT) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 2) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'soot aerosol.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!     Check the spectral file for biomass aerosol modes.
!     (Only required for the direct effect).
!
      IF (L_USE_BMASS_DIRECT .OR. L_USE_ARCLBIOM) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I) == IP_BIOMASS_1).OR.               &
     &              (TYPE_AEROSOL(I) == IP_BIOMASS_2) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 2) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'biomass aerosol.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!
!     Check the spectral file for sea-salt aerosol modes. (Also only
!     required for the direct effect).
!
      IF (L_USE_SEASALT_DIRECT .OR. L_USE_ARCLSSLT) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosols required for this scheme.
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( (TYPE_AEROSOL(I) == IP_SEASALT_FILM).OR.            &
     &              (TYPE_AEROSOL(I) == IP_SEASALT_JET) ) THEN
                  N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                  INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                  N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 2) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'sea-salt aerosol.'
               RETURN
!
            ENDIF

         ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!     Check the spectral file for biogenic aerosol. (direct effect
!     only).
!
      IF (L_USE_BIOGENIC) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosol required for this scheme
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( TYPE_AEROSOL(I) == IP_BIOGENIC ) THEN
                 N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                 INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                 N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 1) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'biogenic aerosol.'
               RETURN
!
            ENDIF

          ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!     Check the spectral file for fossil-fuel organic carbon aerosol. 
!     (direct effect only).
!
      IF (L_USE_OCFF_DIRECT .OR. L_USE_ARCLOCFF) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosol required for this scheme
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF (  (TYPE_AEROSOL(I) == IP_OCFF_FRESH) .OR.            &
     &               (TYPE_AEROSOL(I) == IP_OCFF_AGED) ) THEN
                 N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                 INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                 N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 2) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'fossil fuel org.carb. aerosol.'
               RETURN
!
            ENDIF

          ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!     Check the spectral file for delta aerosol. (direct effect
!     only).
!
      IF (L_USE_ARCLDLTA) THEN
!
         IF (L_PRESENT(11)) THEN ! aerosol block present in spec file
!
!           Search for the aerosol required for this scheme
            N_AEROSOL_FOUND=0
            DO I=1, N_AEROSOL
!
               IF ( TYPE_AEROSOL(I) == IP_DELTA ) THEN
                 N_AEROSOL_RETAIN=N_AEROSOL_RETAIN+1
                 INDEX_AEROSOL_RETAIN(N_AEROSOL_RETAIN)=I
                 N_AEROSOL_FOUND=N_AEROSOL_FOUND+1
               ENDIF
!
            ENDDO
!
            IF (N_AEROSOL_FOUND /= 1) THEN
!
               IERR=I_ERR_FATAL
               CMESSAGE='The SW Spectral file lacks some '              &
     &            //'delta aerosol.'
               RETURN
!
            ENDIF

          ELSE
!
            IERR=I_ERR_FATAL
            CMESSAGE='SW Spectral file contains no aerosol data.'
            RETURN
!
         ENDIF
!
      ENDIF
!
!     Set an appropriate reduced dimension.
      NPD_AEROSOL_SPECIES_SW=MAX(N_AEROSOL_RETAIN, 1)
!
!     Set the allowed number of humidities from the number of
!     retained aerosols.
!
      IF (L_PRESENT(11)) THEN
         DO I=1, N_AEROSOL_RETAIN
            IF (I_AEROSOL_PARAMETRIZATION(INDEX_AEROSOL_RETAIN(I)) ==   &
     &         IP_AEROSOL_PARAM_MOIST) THEN
               NPD_HUMIDITIES_SW=MAX(NPD_HUMIDITIES_SW                  &
     &            , NHUMIDITY(INDEX_AEROSOL_RETAIN(I)))
            ENDIF
         ENDDO
      ENDIF
!
!
!
!
!     TRANSFER THE LARGE NAMELIST TO THE REDUCED SPECTRUM.
!
!
! DEPENDS ON: r2_compress_spectrum
      CALL R2_COMPRESS_SPECTRUM(                                        &
!                       Spectral Array in Namelist
     &     L_PRESENT                                                    &
     &   , N_BAND, WAVE_LENGTH_SHORT , WAVE_LENGTH_LONG                 &
     &   , N_BAND_EXCLUDE, INDEX_EXCLUDE                                &
     &   , SOLAR_FLUX_BAND, RAYLEIGH_COEFFICIENT                        &
     &   , N_ABSORB, N_BAND_ABSORB, INDEX_ABSORB, TYPE_ABSORB           &
     &   , L_RETAIN_ABSORB, N_ABSORB_RETAIN, INDEX_ABSORB_RETAIN        &
     &   , COMPRESSED_INDEX, I_BAND_ESFT, K_ESFT, W_ESFT, I_SCALE_ESFT  &
     &   , I_SCALE_FNC, SCALE_VECTOR, P_REFERENCE, T_REFERENCE          &
     &   , N_DEG_FIT, THERMAL_COEFFICIENT, T_REF_PLANCK                 &
     &   , I_SPEC_SURFACE, L_SURFACE, SURFACE_ALBEDO                    &
     &   , N_DIR_ALBEDO_FIT, DIRECT_ALBEDO_PARM, EMISSIVITY_GROUND      &
     &   , N_BAND_CONTINUUM, INDEX_CONTINUUM, INDEX_WATER               &
     &   , K_CONTINUUM, I_SCALE_FNC_CONT, SCALE_CONTINUUM               &
     &   , P_REF_CONTINUUM, T_REF_CONTINUUM                             &
     &   , L_DROP_TYPE, I_DROP_PARAMETRIZATION, DROP_PARAMETER_LIST     &
     &   , DROP_PARM_MIN_DIM, DROP_PARM_MAX_DIM                         &
     &   , L_ICE_TYPE, I_ICE_PARAMETRIZATION, ICE_PARAMETER_LIST        &
     &   , ICE_PARM_MIN_DIM, ICE_PARM_MAX_DIM                           &
     &   , N_AEROSOL, TYPE_AEROSOL                                      &
     &   , N_AEROSOL_RETAIN, INDEX_AEROSOL_RETAIN                       &
     &   , L_AEROSOL_SPECIES, AEROSOL_ABSORPTION                        &
     &   , AEROSOL_SCATTERING, AEROSOL_ASYMMETRY                        &
     &   , NHUMIDITY, HUMIDITIES, I_AEROSOL_PARAMETRIZATION             &
     &   , L_DOPPLER_PRESENT, DOPPLER_CORRECTION, .FALSE.               &
     &   , N_AOD_WAVEL, AOD_ABSORPTION, AOD_SCATTERING, I_AOD_TYPE      &
!                       Reduced Spectral Array
#include "swsarg3a.h"
     &   )
!
!
!
      RETURN
      END SUBROUTINE R2_SW_SPECIN
#endif
#endif
