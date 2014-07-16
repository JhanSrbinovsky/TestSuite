


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to initialize array of fluxes.
!
! Method:
!       A value is passed to the routine. This is used to initialize
!       arrays of total and direct, overall and clear fluxes as are
!       required for the calculation being performed.
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
      SUBROUTINE INITIALIZE_FLUX(N_PROFILE, N_LAYER, N_AUGMENT          &
     &   , ISOLIR                                                       &
     &   , FLUX_DIRECT, FLUX_TOTAL                                      &
     &   , L_CLEAR                                                      &
     &   , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                          &
     &   , VALUE                                                        &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   , L_NET                                                        &
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
!     INCLUDE COMDECKS.
! SPCRG3A defines flags for different portions of the spectrum in
! two-stream radiation code.
      INTEGER,PARAMETER:: IP_SOLAR=1
      INTEGER,PARAMETER:: IP_INFRA_RED=2
! SPCRG3A end
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , N_AUGMENT                                                    &
!             LENGTH OF LONG FLUX VECTOR
     &   , ISOLIR
!             SPECTRAL REGION
      LOGICAL                                                           &
                !,INTENT(IN)
     &     L_CLEAR
!             CLEAR-SKY FIELD FLAG
      REAL                                                              &
                !, INTENT(IN)
     &     VALUE
!             VALUE FOR INITIALIZATION
      REAL                                                              &
                !, INTENT(OUT)
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)                       &
!             DIRECT FLUX
     &   , FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)                       &
!             TOTAL FLUX
     &   , FLUX_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)                 &
!             CLEAR DIRECT FLUX
     &   , FLUX_TOTAL_CLEAR(NPD_PROFILE, 2*NPD_LAYER+2)
!             CLEAR TOTAL FLUX
!     VARIABLES SPECIFICALLY FOR THE UNIFIED MODEL.
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_NET
!             CALCULATE NET FLUXES
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
!
      IF (ISOLIR == IP_SOLAR) THEN
         DO I=0, N_LAYER
            DO L=1, N_PROFILE
               FLUX_DIRECT(L, I)=VALUE
            ENDDO
         ENDDO
      ENDIF
      DO I=1, N_AUGMENT
         DO L=1, N_PROFILE
            FLUX_TOTAL(L, I)=VALUE
         ENDDO
      ENDDO
!
      IF (L_CLEAR) THEN
         IF (ISOLIR == IP_SOLAR) THEN
            DO I=0, N_LAYER
               DO L=1, N_PROFILE
                  FLUX_DIRECT_CLEAR(L, I)=VALUE
               ENDDO
            ENDDO
         ENDIF
         DO I=1, N_AUGMENT
            DO L=1, N_PROFILE
               FLUX_TOTAL_CLEAR(L, I)=VALUE
            ENDDO
         ENDDO
!        INITIALIZE EXTRA AREAS IN THE CLEAR ARRAY FOR DIAGNOSTIC
!        CALCULATIONS:
         IF (L_NET) THEN
            DO L=1, N_PROFILE
               FLUX_TOTAL_CLEAR(L, 2*N_LAYER+2)=VALUE
            ENDDO
         ENDIF
      ENDIF
!
!
      RETURN
      END SUBROUTINE INITIALIZE_FLUX
