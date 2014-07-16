#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to increment radiances for radiation from the top.
!
! Purpose:
!   The contribution to the solution of radiances transmitted from
!   the top boundary is evaluated. In the IR where differential
!   radiances are used the radiance at the top will be the Planckian
!   radiance at that temperature. In idealized tests an incident
!   flux may be prescribed.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CALC_TOP_RAD(N_PROFILE, TAU                            &
     &  , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                  &
     &  , N_DIRECTION, MU_V                                             &
     &  , FLUX_INC_DOWN                                                 &
     &  , RADIANCE                                                      &
     &  , ND_PROFILE, ND_LAYER, ND_DIRECTION, ND_VIEWING_LEVEL          &
     &  )
!
!
      IMPLICIT NONE
!
!
!     Include header files
#include "c_kinds.h"
#include "c_pi.h"
!
!     Dummy arguments.
!     Sizes of arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ND_VIEWING_LEVEL                                              &
!           Size allocated for levels where radiances are calculated
     &  , ND_DIRECTION
!           Size allocated for viewing directions
!
!
!     The atmosphere:
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE
!           Number of atmospheric profiles
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU(ND_PROFILE, ND_LAYER)
!           Optical depths
!
!     Viewing geometry:
      INTEGER, INTENT(IN) ::                                            &
     &    N_DIRECTION                                                   &
!           Number of directions
     &  , N_VIEWING_LEVEL                                               &
!           Number of levels where the radiance is calculated
     &  , I_RAD_LAYER(ND_VIEWING_LEVEL)
!           Indices of layers containing viewing levels
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FRAC_RAD_LAYER(ND_VIEWING_LEVEL)                              &
!           Fraction optical depth into its layer of the
!           viewing level
     &  , MU_V(ND_PROFILE, ND_DIRECTION)
!           Cosines of polar viewing angles
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FLUX_INC_DOWN(ND_PROFILE)
!           Isotropic incident flux
!
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    RADIANCE(ND_PROFILE, ND_VIEWING_LEVEL, ND_DIRECTION)
!           Radiances (to be incremented by the contribution of
!           the particular integral)
!
!     Local variables
      INTEGER                                                           &
     &    L                                                             &
!           Loop variable (points)
     &  , IV                                                            &
!           Loop variable (viewing levels)
     &  , ID                                                            &
!           Loop variable (directions)
     &  , I                                                             &
!           Loop variable
     &  , LL
!           Loop variable
      INTEGER                                                           &
     &    N_LIST_DOWN                                                   &
!           Number of points where the viewing direction is upward
     &  , LIST_DOWN(ND_PROFILE)
!           List of points where the viewing direction is upward
      REAL  (Real64) ::                                                 &
     &    TAU_C(ND_PROFILE, ND_VIEWING_LEVEL)
!           Cumulative optical depths to the viewing level
!
!
!
!     Calculate the cumulative optical depth from the
!     top of the atmosphere to each viewing level.
      DO IV=1, N_VIEWING_LEVEL
!
        DO L=1, N_PROFILE
          TAU_C(L, IV)                                                  &
     &      =FRAC_RAD_LAYER(IV)*TAU(L, I_RAD_LAYER(IV))
        ENDDO
        DO I=I_RAD_LAYER(IV)-1, 1, -1
          DO L=1, N_PROFILE
            TAU_C(L, IV)=TAU_C(L, IV)+TAU(L, I)
          ENDDO
        ENDDO
!
      ENDDO
!
!
      DO ID=1, N_DIRECTION
!
!       Collect downward directions.
        N_LIST_DOWN=0
        DO L=1, N_PROFILE
          IF (MU_V(L, ID) <  0.0E+00_Real64) THEN
            N_LIST_DOWN=N_LIST_DOWN+1
            LIST_DOWN(N_LIST_DOWN)=L
          ENDIF
        ENDDO
!
        DO IV=1, N_VIEWING_LEVEL
          DO LL=1, N_LIST_DOWN
            L=LIST_DOWN(LL)
            RADIANCE(L, IV, ID)=RADIANCE(L, IV, ID)                     &
     &        +(FLUX_INC_DOWN(L)/PI)*EXP(TAU_C(L, IV)/MU_V(L, ID))
          ENDDO
        ENDDO
!
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE CALC_TOP_RAD
#endif
