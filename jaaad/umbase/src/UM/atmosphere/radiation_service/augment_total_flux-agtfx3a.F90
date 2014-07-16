#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to increment the total flux within a spectral band.
!
! Method:
!       The total flux is incremented by a multiple of the flux
!       flux within a spectral band. This routine is similar to
!       AUGMENT_FLUX, but here the Planckian flux must be
!       incremented in the IR.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE AUGMENT_TOTAL_FLUX(N_PROFILE, N_LAYER, N_AUGMENT       &
     &   , ISOLIR, L_CLEAR, L_NET                                       &
     &   , WEIGHT_BAND, PLANCK_SOURCE_BAND                              &
     &   , FLUX_DIRECT, FLUX_TOTAL                                      &
     &   , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND                            &
     &   , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                          &
     &   , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND                &
     &   , PLANCK_FLUX                                                  &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR                      &
     &   )
!
!
      IMPLICIT NONE
!
!
!     SIZES OF DUMMY ARRAYS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     INCLUDE COMDECKS
#include "spcrg3a.h"
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , N_AUGMENT                                                    &
!             LENGTH OF VECTOR TO AUGMENT
     &   , ISOLIR
!             SPECTRAL REGION
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CLEAR                                                      &
!             CLEAR FLUXES CALCULATED
     &   , L_NET
!             CALCULATE NET FLUXES
      REAL                                                              &
            !, INTENT(IN)
     &     WEIGHT_BAND                                                  &
!             WEIGHTING FACTOR FOR BAND
     &   , PLANCK_SOURCE_BAND(NPD_PROFILE, 0: NPD_LAYER)                &
!             PLANCK FUNCTION IN BAND
     &   , FLUX_DIRECT_BAND(NPD_PROFILE, 0: NPD_LAYER)                  &
!             DIRECT FLUX IN BAND
     &   , FLUX_TOTAL_BAND(NPD_PROFILE, 2*NPD_LAYER+2)                  &
!             TOTAL FLUX IN BAND
     &   , FLUX_DIRECT_CLEAR_BAND(NPD_PROFILE, 0: NPD_LAYER)            &
!             CLEAR DIRECT FLUX IN BAND
     &   , FLUX_TOTAL_CLEAR_BAND(NPD_PROFILE, 2*NPD_LAYER+2)
!             CLEAR TOTAL FLUX IN BAND
      REAL                                                              &
            !, INTENT(INOUT)
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)                       &
!             DIRECT FLUX
     &   , FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)                       &
!             TOTAL FLUX
     &   , FLUX_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)                 &
!             CLEAR DIRECT FLUX
     &   , FLUX_TOTAL_CLEAR(NPD_PROFILE, 2*NPD_LAYER+2)                 &
!             CLEAR TOTAL FLUX
     &   , PLANCK_FLUX(NPD_PROFILE, 0: NPD_LAYER)
!             PLANCKIAN FLUX AT EACH LAYER
!     VARIABLES SPECIFIC TO THE UM.
      REAL                                                              &
                !, INTENT(IN)
     &     ALBEDO_SURFACE_DIFF(NPD_PROFILE)                             &
!             DIFFUSE ALBEDO
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)
!             DIRECT ALBEDO
!
!     LOCAL ARGUMENTS.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
!
!     INCREMENT THE TOTAL FLUXES.
      IF (ISOLIR == IP_SOLAR) THEN
         DO I=0, N_LAYER
            DO L=1, N_PROFILE
               FLUX_DIRECT(L, I)=FLUX_DIRECT(L, I)                      &
     &            +WEIGHT_BAND*FLUX_DIRECT_BAND(L, I)
            ENDDO
         ENDDO
      ENDIF
      DO I=1, N_AUGMENT
         DO L=1, N_PROFILE
            FLUX_TOTAL(L, I)=FLUX_TOTAL(L, I)                           &
     &         +WEIGHT_BAND*FLUX_TOTAL_BAND(L, I)
         ENDDO
      ENDDO
!
      IF (L_CLEAR) THEN
         IF (ISOLIR == IP_SOLAR) THEN
            DO I=0, N_LAYER
               DO L=1, N_PROFILE
                  FLUX_DIRECT_CLEAR(L, I)=FLUX_DIRECT_CLEAR(L, I)       &
     &               +WEIGHT_BAND*FLUX_DIRECT_CLEAR_BAND(L, I)
               ENDDO
            ENDDO
         ENDIF
         DO I=1, N_AUGMENT
            DO L=1, N_PROFILE
               FLUX_TOTAL_CLEAR(L, I)=FLUX_TOTAL_CLEAR(L, I)            &
     &            +WEIGHT_BAND*FLUX_TOTAL_CLEAR_BAND(L, I)
            ENDDO
         ENDDO
!        BOTH UPWARD AND DOWNWARD FLUXES ARE NEEDED AT THE SURFACE
!        FOR DIAGNOSTICS. IF THE NET FLUX IS CALCULATED WE DETERMINE
!        THE DIFFUSE UPWARD FLUX AND PUT IT WHERE IT BELONGS IN THE
!        ARRAY. IF THE FULL FLUXES ARE CALCULATED NO ACTION IS NEEDED.
         IF (L_NET) THEN
            DO L=1, N_PROFILE
               FLUX_TOTAL_CLEAR(L, 2*N_LAYER+2)                         &
     &            =(ALBEDO_SURFACE_DIFF(L)                              &
     &            *FLUX_TOTAL_CLEAR_BAND(L, N_LAYER+1)                  &
     &            +ALBEDO_SURFACE_DIR(L)                                &
     &            *FLUX_DIRECT_CLEAR_BAND(L, N_LAYER))                  &
     &            /(1.0E+00-ALBEDO_SURFACE_DIFF(L))
            ENDDO
         ENDIF
      ENDIF
!
!     SUM THE PLANCKIAN FLUXES FOR LATER ADDITION TO THE DIFFERENTIAL
!     FLUXES. THIS IS NOT NECESSARY FOR A NET-FLUX SCHEME
      IF ( (ISOLIR == IP_INFRA_RED).AND.(.NOT.L_NET) ) THEN
         DO I=0, N_LAYER
            DO L=1, N_PROFILE
               PLANCK_FLUX(L, I)=PLANCK_FLUX(L, I)                      &
     &            +WEIGHT_BAND*PLANCK_SOURCE_BAND(L, I)
            ENDDO
         ENDDO
      ENDIF
!
!
      RETURN
      END SUBROUTINE AUGMENT_TOTAL_FLUX
#endif
#endif
