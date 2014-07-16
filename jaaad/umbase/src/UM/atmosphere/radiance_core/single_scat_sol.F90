#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the singly scattered solar radiance.
!
! Purpose:
!   This subroutine is used to increment the radiances in the
!   required directions on the viewing levels with the singly
!   scattered solar radiance.
!
! Method:
!   Each direction is considered in turn. For each layer of the
!   atmosphere an angular factor involving the phase function
!   and for each viewing level a geometric factor involving the
!   optical depth between the layer in question and the viewing
!   level is calculated. The product of these with the solar
!   beam gives the contribution of that layer to the increment
!   to the radiance.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SINGLE_SCAT_SOL(N_PROFILE, N_LAYER                     &
     &  , N_DIRECTION, DIRECTION                                        &
     &  , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                  &
     &  , I_DIRECT, MU_0                                                &
     &  , TAU, OMEGA, PHASE_FNC_SOLAR                                   &
     &  , RADIANCE                                                      &
     &  , ND_PROFILE, ND_RADIANCE_PROFILE                               &
     &  , ND_LAYER, ND_DIRECTION, ND_VIEWING_LEVEL                      &
     &  )
!
!
!
      IMPLICIT NONE
!
!
!     Dummy arguments:
!     Sizes of arrays:
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_RADIANCE_PROFILE                                           &
!           Size allocated for atmospheric profiles where radiances
!           are calculated
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ND_VIEWING_LEVEL                                              &
!           Size allocated for viewing levels
     &  , ND_DIRECTION
!           Size allocated for viewing directions
!
!     Include header files:
#include "c_kinds.h"
#include "c_pi.h"
!
!     Atmospheric structure:
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of atmospheric profiles
     &  , N_LAYER
!           Number of atmospheric layers
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
     &    DIRECTION(ND_RADIANCE_PROFILE, ND_DIRECTION, 2)               &
!           Cosines of polar viewing angles
     &  , FRAC_RAD_LAYER(ND_VIEWING_LEVEL)
!           Fraction optical depth into its layer of the viewing level
!
!     Solar Radiances
      REAL  (Real64), INTENT(IN) ::                                     &
     &    I_DIRECT(ND_PROFILE, 0: ND_LAYER)                             &
!           Direct solar radiances
     &  , MU_0(ND_PROFILE)
!           Cosines of solar zenith angles
!
!     Optical properties of the atmosphere
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU(ND_PROFILE, ND_LAYER)                                     &
!           Optical depths
     &  , OMEGA(ND_PROFILE, ND_LAYER)                                   &
!           Albedos of single scattering
     &  , PHASE_FNC_SOLAR(ND_RADIANCE_PROFILE, ND_LAYER, ND_DIRECTION)
!           Phase function
!
!
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    RADIANCE(ND_RADIANCE_PROFILE, ND_VIEWING_LEVEL, ND_DIRECTION)
!           Radiances
!
!
!     Local variables
      INTEGER                                                           &
     &    L                                                             &
!           Loop variable (points)
     &  , LL                                                            &
!           Loop variable
     &  , I                                                             &
!           Loop variable
     &  , II                                                            &
!           Loop variable
     &  , IR                                                            &
!           Loop variable (radiative levels)
     &  , ID
!           Loop variable (directions)
      INTEGER                                                           &
     &    N_LIST_UP                                                     &
!           Numbers of points where the current viewing direction
!           points up
     &  , LIST_UP(ND_PROFILE)                                           &
!           List up points with where the current viewing direction
!           points up
     &  , N_LIST_DOWN                                                   &
!           Numbers of points where the current viewing direction
!           points down
     &  , LIST_DOWN(ND_PROFILE)
!           List up points with where the current viewing direction
!           points up
      REAL  (Real64) ::                                                 &
     &    GEOM(ND_PROFILE)                                              &
!           Geometrical factor
     &  , M_SLANT_DEPTH_NEAR(ND_PROFILE)                                &
!           Minus slantwise optical distance between the radiance
!           level and the nearer boundary of the current layer
     &  , M_SLANT_DEPTH_FAR(ND_PROFILE)                                 &
!           Minus slantwise optical distance between the radiance
!           level and the farther boundary of the current layer
     &  , TAU_I(ND_PROFILE)                                             &
!           Optical depth of the relevant part of the current layer
     &  , TRANS_D(ND_PROFILE)                                           &
!           Direct transmission from the layer containing the viewing
!           level to the viewung level
     &  , D_MU
!           Difference in cosines of directions
!
!     Variables related to the treatment of ill-conditioning
      REAL  (Real64) ::                                                 &
     &    EPS_R                                                         &
!           The smallest real number such that 1.0-EPS_R is not 1
!           to the computer's precision
     &  , SQ_EPS_R                                                      &
!           The square root of the above
     &  , ETA                                                           &
!           The conditioning weight
     &  , ETA_NM
!           The conditioning multiplier applied in the numerator
!
!
!
!     Set the tolerances used in avoiding ill-conditioning, testing
!     on any variable.
      EPS_R=EPSILON(MU_0(1))
      SQ_EPS_R=SQRT(EPS_R)
!
!     Consider each direction in turn. Collect points where the
!     viewing direction is upward and points where it is downward,
!     then calculate the angular and geometric factors.
      DO ID=1, N_DIRECTION
!
!       Collect points where the viewing direction is upward:
!       (Horizontal directions are not allowed).
        N_LIST_UP=0
        DO L=1, N_PROFILE
          IF (DIRECTION(L, ID, 1) >  0.0E+00_Real64) THEN
            N_LIST_UP=N_LIST_UP+1
            LIST_UP(N_LIST_UP)=L
          ENDIF
        ENDDO
!
!       Collect points where the viewing direction is downward:
        N_LIST_DOWN=0
        DO L=1, N_PROFILE
          IF (DIRECTION(L, ID, 1) <  0.0E+00_Real64) THEN
            N_LIST_DOWN=N_LIST_DOWN+1
            LIST_DOWN(N_LIST_DOWN)=L
          ENDIF
        ENDDO
!
!       Go through each atmospheric layer calculating the radiance
!       at each observing level.
        DO I=1, N_LAYER
!
!         Calculate the geometric factors:
          DO IR=1, N_VIEWING_LEVEL
!
!           Upward Radiances:
!           Contributions arise only from layers below the viewing
!           level.
            IF (I >= I_RAD_LAYER(IR)) THEN
!
!             Calculate minus the slantwise optical depths to the
!             boundaries of the layer. If the level where the radiance
!             is required lies in the current layer we perform the
!             calculation for a temporary layer reaching from the
!             viewing level to the bottom of the current layer.
              IF (I >  I_RAD_LAYER(IR)) THEN
!               Full layers are required.
                DO LL=1, N_LIST_UP
                  L=LIST_UP(LL)
                  M_SLANT_DEPTH_NEAR(L)                                 &
     &              =(1.0E+00_Real64-FRAC_RAD_LAYER(IR))                &
     &              *TAU(L, I_RAD_LAYER(IR))
                ENDDO
                DO II=I_RAD_LAYER(IR)+1, I-1
                  DO LL=1, N_LIST_UP
                    L=LIST_UP(LL)
                    M_SLANT_DEPTH_NEAR(L)                               &
     &               =M_SLANT_DEPTH_NEAR(L)+TAU(L, II)
                  ENDDO
                ENDDO
                DO LL=1, N_LIST_UP
                  L=LIST_UP(LL)
                  M_SLANT_DEPTH_NEAR(L)                                 &
     &              =-M_SLANT_DEPTH_NEAR(L)/DIRECTION(L, ID, 1)
                  M_SLANT_DEPTH_FAR(L)=M_SLANT_DEPTH_NEAR(L)            &
     &              -TAU(L, I)/DIRECTION(L, ID, 1)
!                 Collect the local optical depth to allow the use of
!                 generic code later.
                  TAU_I(L)=TAU(L, I)
                  TRANS_D(L)=1.0E+00_Real64
                ENDDO
              ELSE IF (I == I_RAD_LAYER(IR)) THEN
!               The viewing level lies in the current layer.
                DO LL=1, N_LIST_UP
                  L=LIST_UP(LL)
                  M_SLANT_DEPTH_NEAR(L)=0.0E+00_Real64
                  M_SLANT_DEPTH_FAR(L)                                  &
     &              =-(1.0E+00_Real64-FRAC_RAD_LAYER(IR))*TAU(L, I)     &
     &              /DIRECTION(L, ID, 1)
                  TAU_I(L)=(1.0E+00_Real64-FRAC_RAD_LAYER(IR))*TAU(L, I)
                  TRANS_D(L)                                            &
     &              =EXP(-FRAC_RAD_LAYER(IR)*TAU(L, I)/MU_0(L))
                ENDDO
              ENDIF
!
!
!             Set the geometrical term and increment the radiance.
              DO LL=1, N_LIST_UP
                L=LIST_UP(LL)
                GEOM(L)=(MU_0(L)/(MU_0(L)+DIRECTION(L, ID, 1)))         &
     &            *(EXP(M_SLANT_DEPTH_NEAR(L))                          &
     &            -EXP(M_SLANT_DEPTH_FAR(L)-TAU_I(L)/MU_0(L)))
                RADIANCE(L, IR, ID)=RADIANCE(L, IR, ID)                 &
     &            +I_DIRECT(L, I-1)*TRANS_D(L)*GEOM(L)*(OMEGA(L, I)     &
     &            /(4.0E+00_Real64*PI))*PHASE_FNC_SOLAR(L, I, ID)
              ENDDO
!
            ENDIF
!
!
!           Downward Radiances:
!           Contributions arise only from layers above the viewing
!           level.
            IF (I <= I_RAD_LAYER(IR)) THEN
!
!             Calculate the slantwise optical depths to the
!             boundaries of the layer. If the observing level lies
!             within the current layer we perform the calculation for
!             a layer reaching from the top of the current layer to
!             the observing level.
              IF (I <  I_RAD_LAYER(IR)) THEN
                DO LL=1, N_LIST_DOWN
                  L=LIST_DOWN(LL)
                  M_SLANT_DEPTH_NEAR(L)                                 &
     &            =FRAC_RAD_LAYER(IR)*TAU(L, I)
                ENDDO
                DO II=I_RAD_LAYER(IR)-1, I+1, -1
                  DO LL=1, N_LIST_DOWN
                    L=LIST_DOWN(LL)
                    M_SLANT_DEPTH_NEAR(L)                               &
     &                =M_SLANT_DEPTH_NEAR(L)+TAU(L, II)
                  ENDDO
                ENDDO
                DO LL=1, N_LIST_DOWN
                  L=LIST_DOWN(LL)
                  M_SLANT_DEPTH_NEAR(L)                                 &
     &              =M_SLANT_DEPTH_NEAR(L)/DIRECTION(L, ID, 1)
                  M_SLANT_DEPTH_FAR(L)=M_SLANT_DEPTH_NEAR(L)            &
     &              +TAU(L, I)/DIRECTION(L, ID, 1)
                  TAU_I(L)=TAU(L, I)
                ENDDO
              ELSE
!               The viewing level lies in the current layer.
                DO LL=1, N_LIST_DOWN
                  L=LIST_DOWN(LL)
                  TAU_I(L)=FRAC_RAD_LAYER(IR)*TAU(L, I)
                  M_SLANT_DEPTH_NEAR(L)=0.0E+00_Real64
                  M_SLANT_DEPTH_FAR(L)=TAU_I(L)/DIRECTION(L, ID, 1)
                ENDDO
              ENDIF
!
!
!             Set the geometrical terms for the solar integral.
              DO LL=1, N_LIST_DOWN
                L=LIST_DOWN(LL)
!               This may exhibit ill-conditioning, so it is perturbed
!               using L'Hopital's rule.
                D_MU=MU_0(L)+DIRECTION(L, ID, 1)
                ETA=EPS_R/(D_MU+SIGN(SQ_EPS_R, D_MU))
                ETA_NM=(1.0E+00_Real64-ETA*TAU_I(L)                     &
     &            /(MU_0(L)*DIRECTION(L, ID, 1)))
                GEOM(L)=(MU_0(L)/(D_MU+ETA))                            &
     &            *(EXP(M_SLANT_DEPTH_NEAR(L)-TAU_I(L)/MU_0(L))         &
     &            *ETA_NM                                               &
     &            -EXP(M_SLANT_DEPTH_FAR(L)))
                RADIANCE(L, IR, ID)=RADIANCE(L, IR, ID)                 &
     &            +I_DIRECT(L, I-1)*GEOM(L)                             &
     &            *(OMEGA(L, I)/(4.0E+00_Real64*PI))                    &
     &            *PHASE_FNC_SOLAR(L, I, ID)
              ENDDO
!
            ENDIF
!
          ENDDO
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SINGLE_SCAT_SOL
#endif
