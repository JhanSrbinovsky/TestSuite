#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to rescale optical depth and albedo.
!
! Method:
!       The standard rescaling formulae are applied.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE RESCALE_TAU_OMEGA(N_PROFILE                            &
     &  , I_LAYER_FIRST, I_LAYER_LAST                                   &
     &  , TAU, OMEGA, FORWARD_SCATTER                                   &
     &  , ND_PROFILE, ND_LAYER, ID_1                                    &
     &  )
!
!
      IMPLICIT NONE
!
! Include header files
#include "c_kinds.h"
!
!
!     Sizes of arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for profiles
     &  , ND_LAYER                                                      &
!           Size allocated for layers
     &  , ID_1
!           Topmost declared layer for optical properties
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , I_LAYER_FIRST                                                 &
!           First layer to rescale
     &  , I_LAYER_LAST
!           First layer to rescale
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FORWARD_SCATTER(ND_PROFILE, ID_1: ND_LAYER)
!           Forward scattering
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    TAU(ND_PROFILE, ID_1: ND_LAYER)                               &
!           Optical depth
     &  , OMEGA(ND_PROFILE, ID_1: ND_LAYER)
!           Albedo of single scattering
!
!     Local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L
!           Loop variable
!
!
!
      DO I=I_LAYER_FIRST, I_LAYER_LAST
        DO L=1, N_PROFILE
          TAU(L, I)=TAU(L, I)*(1.0E+00_Real64                           &
     &      -OMEGA(L, I)*FORWARD_SCATTER(L, I))
          OMEGA(L, I)=OMEGA(L, I)*(1.0E+00_Real64-FORWARD_SCATTER(L, I))&
     &      /(1.0E+00_Real64-OMEGA(L, I)*FORWARD_SCATTER(L, I))
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE RESCALE_TAU_OMEGA
#endif
