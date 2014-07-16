#if defined(A70_1C) || defined(A70_1Z)
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------------
!
!+ Subroutine to calculate weights for ultraviolet fluxes.
!
! Purpose:
!   Weights to calculate the flux in a wavelength interval [a,b].
!   This is a straightforward extension of the subroutine
!   'R2_SET_690NM_WEIGHT'.
!
! Method:
!   Straightforward. The flux is assumed to be linearly distributed
!   across bands.
!
! Current Owner of Code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!       6.2             17-02-04                Original Code
!                                               (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!----------------------------------------------------------------------
      SUBROUTINE R2_SET_UV_WEIGHT(N_BAND                                &
     &   , L_PRESENT                                                    &
     &   , N_BAND_EXCLUDE, INDEX_EXCLUDE                                &
     &   , WAVE_LENGTH_SHORT, WAVE_LENGTH_LONG                          &
     &   , WEIGHT_UV                                                    &
     &   , NPD_BAND_SW, NPD_EXCLUDE_SW, NPD_TYPE_SW                     &
     &   )
!
!
      IMPLICIT NONE
!
#include "uvinterval.h"
!
!     DUMMY VARIABLES:
!
!     DIMENSIONS OF ARRAYS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_BAND_SW                                                  &
!             MAXIMUM NUMBER OF SPECTRAL BANDS
     &   , NPD_EXCLUDE_SW                                               &
!             MAXIMUM NUMBER OF EXCLUDED REGIONS
     &   , NPD_TYPE_SW
!             MAXIMUM NUMBER OF TYPES OF SPECTRAL DATA
!
!     ACTUAL SIZES USED:
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_BAND                                                       &
!             NUMBER OF SPECTRAL BANDS
     &   , N_BAND_EXCLUDE(NPD_BAND_SW)                                  &
!             NUMBER OF EXCLUDED REGIONS IN BANDS
     &   , INDEX_EXCLUDE(NPD_EXCLUDE_SW, NPD_BAND_SW)
!             INDICES OF EXCLUDED REGIONS IN BANDS
!
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_PRESENT(0: NPD_TYPE_SW)
!             FLAG FOR TYPES OF SPECTRAL DATA PRESENT
      REAL                                                              &
                !, INTENT(IN)
     &     WAVE_LENGTH_SHORT(NPD_BAND_SW)                               &
!             SHORT WAVELENGTH LIMITS OF BANDS
     &   , WAVE_LENGTH_LONG(NPD_BAND_SW)
!             LONG WAVELENGTH LIMITS OF BANDS
!
!
!     WEIGHTS SET.
      REAL                                                              &
                !, INTENT(OUT)
     &     WEIGHT_UV(NPD_BAND_SW)
!             WEIGHTS APPLYING TO EACH BAND
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE

      REAL                                                              &
     &     TOTAL_ENERGY_RANGE                                           &
!             TOTAL RANGE OF ENERGIES COVERED BY BAND
     &   , ENERGY_RANGE_UV
!             RANGE OF ENERGIES IN BAND IN THE WAVELENGHT INTERVAL
!
!
!
      WEIGHT_UV=0.0E+00

      DO I=1, N_BAND

         IF (WAVE_LENGTH_LONG(I) <  UV_INTERVAL_SHORT)                  &
     &                                   WEIGHT_UV(I)=0.0E+00

         IF (WAVE_LENGTH_SHORT(I) >  UV_INTERVAL_LONG)                  &
     &                                   WEIGHT_UV(I)=0.0E+00

         IF (WAVE_LENGTH_SHORT(I) <  UV_INTERVAL_SHORT) THEN

            TOTAL_ENERGY_RANGE=1.0E+00/WAVE_LENGTH_SHORT(I)             &
     &                        -1.0E+00/WAVE_LENGTH_LONG(I)

            IF (WAVE_LENGTH_LONG(I) >  UV_INTERVAL_LONG) THEN
               ENERGY_RANGE_UV=1.0/UV_INTERVAL_SHORT                    &
     &                        -1.0/UV_INTERVAL_LONG
               IF (L_PRESENT(14)) THEN
!                REMOVE CONTRIBUTIONS FROM EXCLUDED BANDS.
                 DO J=1, N_BAND_EXCLUDE(I)
                   IF (WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I)) <         &
     &                UV_INTERVAL_SHORT) THEN
                      ENERGY_RANGE_UV=ENERGY_RANGE_UV                   &
     &                 -1.0E+00/UV_INTERVAL_SHORT                       &
     &                 +1.0E+00/WAVE_LENGTH_LONG(INDEX_EXCLUDE(J, I))
                   ELSE IF (WAVE_LENGTH_LONG(INDEX_EXCLUDE(J, I)) >     &
     &                UV_INTERVAL_LONG) THEN
                      ENERGY_RANGE_UV=ENERGY_RANGE_UV                   &
     &                 -1.0E+00/WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I))  &
     &                 +1.0E+00/UV_INTERVAL_LONG
                   ELSE
                      ENERGY_RANGE_UV=ENERGY_RANGE_UV                   &
     &                 -1.0E+00/WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I))  &
     &                 +1.0E+00/WAVE_LENGTH_LONG(INDEX_EXCLUDE(J, I))
                   ENDIF
                   TOTAL_ENERGY_RANGE=TOTAL_ENERGY_RANGE                &
     &                -1.0E+00/WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I))   &
     &                +1.0E+00/WAVE_LENGTH_LONG(INDEX_EXCLUDE(J, I))
                 ENDDO
               ENDIF
              WEIGHT_UV(I)=ENERGY_RANGE_UV/TOTAL_ENERGY_RANGE
            ENDIF

            IF ((WAVE_LENGTH_LONG(I) <  UV_INTERVAL_LONG).AND.          &
     &          (WAVE_LENGTH_LONG(I) >  UV_INTERVAL_SHORT)) THEN
               ENERGY_RANGE_UV=1.0E+00/UV_INTERVAL_SHORT                &
     &                        -1.0E+00/WAVE_LENGTH_LONG(I)
               IF (L_PRESENT(14)) THEN
!                REMOVE CONTRIBUTIONS FROM EXCLUDED BANDS.
                 DO J=1, N_BAND_EXCLUDE(I)
                   IF (WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I)) >         &
     &                UV_INTERVAL_SHORT) THEN
                      ENERGY_RANGE_UV=ENERGY_RANGE_UV                   &
     &                 -1.0E+00/WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I))  &
     &                 +1.0E+00/WAVE_LENGTH_LONG(INDEX_EXCLUDE(J, I))
                   ELSE IF (WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I)) <    &
     &                UV_INTERVAL_SHORT) THEN
                      ENERGY_RANGE_UV=ENERGY_RANGE_UV                   &
     &                 -1.0E+00/UV_INTERVAL_SHORT                       &
     &                 +1.0E+00/WAVE_LENGTH_LONG(INDEX_EXCLUDE(J, I))
                   ENDIF
                   TOTAL_ENERGY_RANGE=TOTAL_ENERGY_RANGE                &
     &                -1.0E+00/WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I))   &
     &                +1.0E+00/WAVE_LENGTH_LONG(INDEX_EXCLUDE(J, I))
                 ENDDO
               ENDIF
               WEIGHT_UV(I)=ENERGY_RANGE_UV/TOTAL_ENERGY_RANGE
            ENDIF

         ELSE

            TOTAL_ENERGY_RANGE=1.0E+00/WAVE_LENGTH_SHORT(I)             &
     &                        -1.0E+00/WAVE_LENGTH_LONG(I)

            IF (WAVE_LENGTH_LONG(I) <= UV_INTERVAL_LONG)                &
     &                                     WEIGHT_UV(I)=1.0

            IF ((WAVE_LENGTH_SHORT(I) <  UV_INTERVAL_LONG).AND.         &
     &          (WAVE_LENGTH_LONG(I) >  UV_INTERVAL_LONG)) THEN
               ENERGY_RANGE_UV=1.0E+00/WAVE_LENGTH_SHORT(I)             &
     &                        -1.0E+00/UV_INTERVAL_LONG

               IF (L_PRESENT(14)) THEN
!                REMOVE CONTRIBUTIONS FROM EXCLUDED BANDS.
                 DO J=1, N_BAND_EXCLUDE(I)
                   IF (WAVE_LENGTH_LONG(INDEX_EXCLUDE(J, I)) <          &
     &                UV_INTERVAL_LONG) THEN
                      ENERGY_RANGE_UV=ENERGY_RANGE_UV                   &
     &                 -1.0E+00/WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I))  &
     &                 +1.0E+00/WAVE_LENGTH_LONG(INDEX_EXCLUDE(J, I))
                   ELSE IF (WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I)) <    &
     &                 UV_INTERVAL_LONG) THEN
                       ENERGY_RANGE_UV=ENERGY_RANGE_UV                  &
     &                  -1.0E+00/WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I)) &
     &                  +1.0E+00/UV_INTERVAL_LONG
                   ENDIF
                   TOTAL_ENERGY_RANGE=TOTAL_ENERGY_RANGE                &
     &                -1.0E+00/WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I))   &
     &                +1.0E+00/WAVE_LENGTH_LONG(INDEX_EXCLUDE(J, I))
                 ENDDO
               ENDIF
               WEIGHT_UV(I)=ENERGY_RANGE_UV/TOTAL_ENERGY_RANGE
            ENDIF

         ENDIF
!
      ENDDO
!
      RETURN
      END SUBROUTINE R2_SET_UV_WEIGHT
#endif
#endif
