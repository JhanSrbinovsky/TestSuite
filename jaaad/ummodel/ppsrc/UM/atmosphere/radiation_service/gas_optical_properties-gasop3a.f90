


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the absorptive extinctions of gases.
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
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE GAS_OPTICAL_PROPERTIES(N_PROFILE, N_LAYER              &
     &   , N_GAS, I_GAS_POINTER, K_ESFT_MONO, GAS_MIX_RATIO             &
     &   , K_GAS_ABS                                                    &
     &   , NPD_PROFILE, NPD_LAYER, NPD_SPECIES                          &
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
     &   , NPD_LAYER                                                    &
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_SPECIES
!             MAXIMUM NUMBER OF GASEOUS SPECIES
!
!     DUMMY VARIABLES.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , N_GAS                                                        &
!             NUMBER OF GASES
     &   , I_GAS_POINTER(NPD_SPECIES)
!             POINTERS TO ACTIVE GASES
      REAL                                                              &
                !, INTENT(IN)
     &     K_ESFT_MONO(NPD_SPECIES)                                     &
!             ESFT EXPONENTS FOR EACH GAS
     &   , GAS_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER, NPD_SPECIES)
!             GAS MIXING RATIOS
      REAL                                                              &
                !, INTENT(OUT)
     &     K_GAS_ABS(NPD_PROFILE, NPD_LAYER)
!             CLEAR ABSORPTIVE EXTINCTION
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I_GAS                                                        &
!             TEMPORARY GAS INDEX
     &   , L                                                            &
!             LOOP VARIABLE
     &   , I                                                            &
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
!
!
!     CALCULATE THE ABSORPTION FOR THE FIRST GAS AND ADD ON THE REST.
      I_GAS=I_GAS_POINTER(1)
      DO J=1, N_LAYER
         DO L=1, N_PROFILE
            K_GAS_ABS(L, J)                                             &
     &         =K_ESFT_MONO(I_GAS)*GAS_MIX_RATIO(L, J, I_GAS)
         ENDDO
      ENDDO
      DO I=2, N_GAS
      I_GAS=I_GAS_POINTER(I)
         DO J=1, N_LAYER
            DO L=1, N_PROFILE
               K_GAS_ABS(L, J)=K_GAS_ABS(L, J)                          &
     &            +K_ESFT_MONO(I_GAS)*GAS_MIX_RATIO(L, J, I_GAS)
            ENDDO
         ENDDO
      ENDDO
!
!
      RETURN
      END SUBROUTINE GAS_OPTICAL_PROPERTIES
