#if defined(A70_1Z)
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
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE GAS_OPTICAL_PROPERTIES(N_PROFILE, N_LAYER              &
     &   , N_GAS, I_GAS_POINTER, K_ESFT_MONO, GAS_MIX_RATIO             &
     &   , K_GAS_ABS                                                    &
     &   , ND_PROFILE, ND_LAYER, ND_SPECIES                             &
     &   )
!
!
      IMPLICIT NONE
!
! Include header files
#include "c_kinds.h"
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &     ND_PROFILE                                                   &
!           Size allocated for atmospheric profiles
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ND_SPECIES
!           Size allocated for gaseous species
!
!     Dummy variables.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_GAS                                                         &
!           Number of gases
     &  , I_GAS_POINTER(ND_SPECIES)
!           Pointers to active gases
      REAL  (Real64), INTENT(IN) ::                                     &
     &    K_ESFT_MONO(ND_SPECIES)                                       &
!           ESFT exponents for each gas
     &  , GAS_MIX_RATIO(ND_PROFILE, ND_LAYER, ND_SPECIES)
!           Gas mixing ratios
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    K_GAS_ABS(ND_PROFILE, ND_LAYER)
!           Clear absorptive extinction
!
!     Local variables.
      INTEGER                                                           &
     &    I_GAS                                                         &
!           Temporary gas `index'
     &  , L                                                             &
!           Loop variable
     &  , I                                                             &
!           Loop variable
     &  , J
!           Loop variable
!
!
!     Calculate the absorption for the first gas and add on the rest.
      I_GAS=I_GAS_POINTER(1)
      DO J=1, N_LAYER
        DO L=1, N_PROFILE
          K_GAS_ABS(L, J)                                               &
     &      =K_ESFT_MONO(I_GAS)*GAS_MIX_RATIO(L, J, I_GAS)
        ENDDO
      ENDDO
      DO I=2, N_GAS
      I_GAS=I_GAS_POINTER(I)
        DO J=1, N_LAYER
          DO L=1, N_PROFILE
            K_GAS_ABS(L, J)=K_GAS_ABS(L, J)                             &
     &        +K_ESFT_MONO(I_GAS)*GAS_MIX_RATIO(L, J, I_GAS)
          ENDDO
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE GAS_OPTICAL_PROPERTIES
#endif
