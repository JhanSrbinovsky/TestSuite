#if defined(A01_3A) || defined(A01_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the solar flux and source terms.
!
! Method:
!       Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             08-05-97                Formulation for
!                                               equivalent extinction
!                                               amended.
!       6.2             15-02-05                Direct flux at
!                                               ground corrected
!                                               for sloping terrain.
!                                               (James Manners)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLAR_SOURCE(N_PROFILE, N_LAYER                        &
     &   , FLUX_INC_DIRECT                                              &
     &   , TRANS_0, SOURCE_COEFF                                        &
     &   , L_SCALE_SOLAR, ADJUST_SOLAR_KE                               &
     &   , FLUX_DIRECT                                                  &
     &   , S_DOWN, S_UP                                                 &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
      USE solinc_data, ONLY: lg_orog_corr, L_orog
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
!     COMDECKS INCLUDED.
#include "dimfix3a.h"
#include "scfpt3a.h"
!
!     DUMMY VARIABLES.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
!
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_SCALE_SOLAR
!             SCALING APPLIED TO SOLAR BEAM
!
      REAL                                                              &
                !, INTENT(IN)
     &     FLUX_INC_DIRECT(NPD_PROFILE)                                 &
!             INCIDENT SOLAR FLUX
     &   , TRANS_0(NPD_PROFILE, NPD_LAYER)                              &
!             DIRECT TRANSMISSION COEFFICIENT
     &   , SOURCE_COEFF(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)       &
!             REFLECTION COEFFICIENT
     &   , ADJUST_SOLAR_KE(NPD_PROFILE, NPD_LAYER)
!             ADJUSTMENT TO SOLAR FLUX
!
!
      REAL                                                              &
                !, INTENT(OUT)
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)                       &
!             DIRECT FLUX
     &   , S_DOWN(NPD_PROFILE, NPD_LAYER)                               &
!             DOWNWARD SOURCE FUNCTION
     &   , S_UP(NPD_PROFILE, NPD_LAYER)
!             UPWARD SOURCE FUNCTION
!
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
!
!
      DO L=1, N_PROFILE
         FLUX_DIRECT(L, 0)=FLUX_INC_DIRECT(L)
      ENDDO
!
!     THE SOLAR FLUX MAY BE MULTIPLIED BY A SCALING FACTOR IF AN
!     EQUIVALENT EXTINCTION IS USED.
!
      IF (L_SCALE_SOLAR) THEN
!
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               FLUX_DIRECT(L, I)                                        &
     &            =FLUX_DIRECT(L, I-1)*TRANS_0(L, I)                    &
     &            *ADJUST_SOLAR_KE(L, I)
               S_UP(L, I)=SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP)           &
     &            *FLUX_DIRECT(L, I-1)
               S_DOWN(L, I)=(SOURCE_COEFF(L, I, IP_SCF_SOLAR_DOWN)      &
     &            -TRANS_0(L, I))*FLUX_DIRECT(L, I-1)                   &
     &            +FLUX_DIRECT(L, I)
            ENDDO
         ENDDO
!
      ELSE
!
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               FLUX_DIRECT(L, I)                                        &
     &            =FLUX_DIRECT(L, I-1)*TRANS_0(L, I)
               S_UP(L, I)=SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP)           &
     &            *FLUX_DIRECT(L, I-1)
               S_DOWN(L, I)=SOURCE_COEFF(L, I, IP_SCF_SOLAR_DOWN)       &
     &            *FLUX_DIRECT(L, I-1)
            ENDDO
         ENDDO
!
      ENDIF
!
!
!     CORRECT THE DIRECT FLUX AT THE GROUND FOR SLOPING TERRAIN

      IF (L_orog) THEN
         FLUX_DIRECT(1:N_PROFILE, N_LAYER) =                            &
     &      FLUX_DIRECT(1:N_PROFILE, N_LAYER) *                         &
     &      lg_orog_corr(1:N_PROFILE)

         S_DOWN(1:N_PROFILE, N_LAYER) =                                 &
     &         S_DOWN(1:N_PROFILE, N_LAYER) +                           &
     &         FLUX_DIRECT(1:N_PROFILE, N_LAYER) *                      &
     &         (lg_orog_corr(1:N_PROFILE) - 1.0)
      ENDIF

      RETURN
      END SUBROUTINE SOLAR_SOURCE
#endif
