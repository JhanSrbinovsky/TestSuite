


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to assign fluxes to the final arrays.
!
! Method:
!       The array FLUX_TOTAL holds the calculated total fluxes
!       (differential fluxes in the IR). These may be net or
!       upward and downward fluxes. Upward and downward fluxes
!       are passed to FLUX_UP and FLUX_DOWN and incremented by the
!       Planckian flux in the IR. Net downward fluxes are assigned
!       to FLUX_DOWN. Equivalent calculations may be performed for
!       the clear-sky fluxes.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE ASSIGN_FLUX(N_PROFILE, N_LAYER                         &
     &   , FLUX_TOTAL, FLUX_TOTAL_CLEAR                                 &
     &   , ISOLIR                                                       &
     &   , PLANCK_FLUX                                                  &
     &   , L_CLEAR, L_NET                                               &
     &   , FLUX_DOWN, FLUX_UP, FLUX_DOWN_CLEAR, FLUX_UP_CLEAR           &
     &   , NPD_PROFILE, NPD_LAYER                                       &
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
!
!     INCLUDE COMDECKS.
! SPCRG3A defines flags for different portions of the spectrum in
! two-stream radiation code.
      INTEGER,PARAMETER:: IP_SOLAR=1
      INTEGER,PARAMETER:: IP_INFRA_RED=2
! SPCRG3A end
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , ISOLIR
!             SPECTRAL REGION
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CLEAR                                                      &
!             CLEAR FLUX FLAG
     &   , L_NET
!             CALCULATE NET FLUXES
      REAL                                                              &
                !, INTENT(IN)
     &     FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)                       &
!             LONG VECTOR OF TOTAL FLUXES
     &   , FLUX_TOTAL_CLEAR(NPD_PROFILE, 2*NPD_LAYER+2)                 &
!             LONG VECTOR OF TOTAL CLEAR FLUXES
     &   , PLANCK_FLUX(NPD_PROFILE, 0: NPD_LAYER)
!             PLANCKIAN FLUXES AT BOUNDARIES
      REAL                                                              &
                !, INTENT(OUT)
     &     FLUX_DOWN(NPD_PROFILE, 0: NPD_LAYER)                         &
!             (NET) TOTAL DOWNWARD FLUXES
     &   , FLUX_UP(NPD_PROFILE, 0: NPD_LAYER)                           &
!             UPWARD FLUXES
     &   , FLUX_DOWN_CLEAR(NPD_PROFILE, 0: NPD_LAYER)                   &
!             CLEAR (NET) TOTAL DOWNWARD FLUXES
     &   , FLUX_UP_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR UPWARD FLUXES

!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
!
      IF (ISOLIR == IP_SOLAR) THEN
         IF (L_NET) THEN
            DO I=0, N_LAYER
               DO L=1, N_PROFILE
                  FLUX_DOWN(L, I)=FLUX_TOTAL(L, I+1)
               ENDDO
            ENDDO
         ELSE
            DO I=0, N_LAYER
               DO L=1, N_PROFILE
                  FLUX_UP(L, I)=FLUX_TOTAL(L, 2*I+1)
                  FLUX_DOWN(L, I)=FLUX_TOTAL(L, 2*I+2)
               ENDDO
            ENDDO
         ENDIF
      ELSE IF (ISOLIR == IP_INFRA_RED) THEN
         IF (L_NET) THEN
!           NO PLANCKIAN CORRECTION IS NECESSARY TO THE NET FLUX.
            DO I=0, N_LAYER
               DO L=1, N_PROFILE
                  FLUX_DOWN(L, I)=FLUX_TOTAL(L, I+1)
               ENDDO
            ENDDO
         ELSE
            DO I=0, N_LAYER
               DO L=1, N_PROFILE
                  FLUX_UP(L, I)=FLUX_TOTAL(L, 2*I+1)                    &
     &               +PLANCK_FLUX(L, I)
                  FLUX_DOWN(L, I)=FLUX_TOTAL(L, 2*I+2)                  &
     &               +PLANCK_FLUX(L, I)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
!
      IF (L_CLEAR) THEN
         IF (ISOLIR == IP_SOLAR) THEN
            IF (L_NET) THEN
               DO I=0, N_LAYER
                  DO L=1, N_PROFILE
                     FLUX_DOWN_CLEAR(L, I)=FLUX_TOTAL_CLEAR(L, I+1)
                  ENDDO
               ENDDO
            ELSE
               DO I=0, N_LAYER
                  DO L=1, N_PROFILE
                     FLUX_UP_CLEAR(L, I)=FLUX_TOTAL_CLEAR(L, 2*I+1)
                     FLUX_DOWN_CLEAR(L, I)=FLUX_TOTAL_CLEAR(L, 2*I+2)
                  ENDDO
               ENDDO
            ENDIF
!
         ELSEIF (ISOLIR == IP_INFRA_RED) THEN
            IF (L_NET) THEN
               DO I=0, N_LAYER
                  DO L=1, N_PROFILE
                     FLUX_DOWN_CLEAR(L, I)=FLUX_TOTAL_CLEAR(L, I+1)
                  ENDDO
               ENDDO
            ELSE
               DO I=0, N_LAYER
                  DO L=1, N_PROFILE
                     FLUX_UP_CLEAR(L, I)=FLUX_TOTAL_CLEAR(L, 2*I+1)     &
     &                  +PLANCK_FLUX(L, I)
                     FLUX_DOWN_CLEAR(L, I)=FLUX_TOTAL_CLEAR(L, 2*I+2)   &
     &                  +PLANCK_FLUX(L, I)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDIF
!
!     NET FLUXES ARE HABITUALLY USED IN THE UNIFIED MODEL, SO THE
!     FOLLOWING REDUCTION IS ALWAYS CARRIED OUT.
      IF (.NOT.L_NET) THEN
         DO I=0, N_LAYER
            DO L=1, N_PROFILE
               FLUX_DOWN(L, I)=FLUX_DOWN(L, I)-FLUX_UP(L, I)
            ENDDO
         ENDDO
         IF (L_CLEAR) THEN
            DO I=0, N_LAYER
               DO L=1, N_PROFILE
                  FLUX_DOWN_CLEAR(L, I)                                 &
     &               =FLUX_DOWN_CLEAR(L, I)-FLUX_UP_CLEAR(L, I)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
!
!
      RETURN
      END SUBROUTINE ASSIGN_FLUX
