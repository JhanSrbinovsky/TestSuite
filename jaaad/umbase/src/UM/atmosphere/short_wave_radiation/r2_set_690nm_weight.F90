#if defined(A01_3A) || defined(A01_3C) || defined(A01_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate weights for the flux below 690 nm.
!
! Purpose:
!   Weights to calculate the flux below 690 nm are set.
!
! Method:
!   Straightforward. The flux is assumed to be linearly distributed
!   across bands.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_690NM_WEIGHT(N_BAND                             &
     &   , L_PRESENT                                                    &
     &   , N_BAND_EXCLUDE, INDEX_EXCLUDE                                &
     &   , WAVE_LENGTH_SHORT, WAVE_LENGTH_LONG                          &
     &   , WEIGHT_690NM                                                 &
     &   , NPD_BAND_SW, NPD_EXCLUDE_SW, NPD_TYPE_SW                     &
     &   )
!
!
!
      IMPLICIT NONE
!
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
!
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
     &     WEIGHT_690NM(NPD_BAND_SW)
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
     &   , ENERGY_RANGE_BELOW_690NM
!             RANGE OF ENERGIES IN BAND BELOW 690 NM
!
!
!
      DO I=1, N_BAND
         IF (WAVE_LENGTH_LONG(I) <  6.9E-07) THEN
            WEIGHT_690NM(I)=1.0E+00
         ELSE IF (WAVE_LENGTH_SHORT(I) >  6.9E-07) THEN
            WEIGHT_690NM(I)=0.0E+00
         ELSE
!
            ENERGY_RANGE_BELOW_690NM=1.0E+00/WAVE_LENGTH_SHORT(I)       &
     &         -1.0E+00/6.9E-07
            TOTAL_ENERGY_RANGE=1.0E+00/WAVE_LENGTH_SHORT(I)             &
     &         -1.0E+00/WAVE_LENGTH_LONG(I)
            IF (L_PRESENT(14)) THEN
!              REMOVE CONTRIBUTIONS FROM EXCLUDED BANDS.
               DO J=1, N_BAND_EXCLUDE(I)
                  IF (WAVE_LENGTH_LONG(INDEX_EXCLUDE(J, I)) <           &
     &               6.9E-07) THEN
                     ENERGY_RANGE_BELOW_690NM=ENERGY_RANGE_BELOW_690NM  &
     &                  -1.0E+00/WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I)) &
     &                  +1.0E+00/WAVE_LENGTH_LONG(INDEX_EXCLUDE(J, I))
                  ELSE IF (WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I)) <     &
     &               6.9E-07) THEN
                     ENERGY_RANGE_BELOW_690NM=ENERGY_RANGE_BELOW_690NM  &
     &                  -1.0E+00/WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I)) &
     &                  +1.0E+00/6.9E-07
                  ENDIF
                  TOTAL_ENERGY_RANGE=TOTAL_ENERGY_RANGE                 &
     &               -1.0E+00/WAVE_LENGTH_SHORT(INDEX_EXCLUDE(J, I))    &
     &               +1.0E+00/WAVE_LENGTH_LONG(INDEX_EXCLUDE(J, I))
               ENDDO
            ENDIF
!
            WEIGHT_690NM(I)=ENERGY_RANGE_BELOW_690NM/TOTAL_ENERGY_RANGE
!
         ENDIF
!
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE R2_SET_690NM_WEIGHT
#endif
