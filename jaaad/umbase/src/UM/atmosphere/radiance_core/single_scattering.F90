#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to find the optical depth and single scattering albedo.
!
! Method:
!       Depending on the treatment of scattering, the optical and
!       and single scattering albedo are determined from the
!       extinctions supplied.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SINGLE_SCATTERING(I_SCATTER_METHOD_BAND                &
     &  , N_PROFILE, I_FIRST_LAYER, I_LAST_LAYER                        &
     &  , D_MASS                                                        &
     &  , K_GREY_TOT, K_EXT_SCAT, K_GAS_ABS                             &
     &  , TAU, OMEGA                                                    &
     &  , ND_PROFILE, ND_LAYER, ID_LT, ID_LB                            &
     &  )
!
!
!
      IMPLICIT NONE
!
!
!     Sizes of arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for profiles
     &  , ND_LAYER                                                      &
!           Size allocated for layers
     &  , ID_LT                                                         &
!           Topmost declared layer for optical properties
     &  , ID_LB
!           Bottom declared layer for optical properties
!
!     Inclusion of header files.
#include "c_kinds.h"
#include "scatter_method_pcf3z.h"
!
!     Dummy variables.
      INTEGER, INTENT(IN) ::                                            &
     &    I_SCATTER_METHOD_BAND
!           Treatment of scattering in this band
!
!                     Atmospheric properties
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , I_FIRST_LAYER                                                 &
!           First layer to consider
     &  , I_LAST_LAYER
!           Last layer to consider
      REAL  (Real64), INTENT(IN) ::                                     &
     &    D_MASS(ND_PROFILE, ND_LAYER)
!           Mass thickness of each layer
!
!                     Optical properties
      REAL  (Real64), INTENT(IN) ::                                     &
     &    K_GREY_TOT(ND_PROFILE, ID_LT: ID_LB)                          &
!           Absorptive extinction
     &  , K_EXT_SCAT(ND_PROFILE, ID_LT: ID_LB)                          &
!           Scattering extinction
     &  , K_GAS_ABS(ND_PROFILE, ND_LAYER)
!           Gaseous extinction
!
!                     Single scattering properties
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    TAU(ND_PROFILE, ID_LT: ID_LB)                                 &
!           Optical depth
     &  , OMEGA(ND_PROFILE, ID_LT: ID_LB)
!           Single scattering albedo
!
!
!
!     Local variables.
      INTEGER                                                           &
     &    L                                                             &
!           Loop variable
     &  , I
!           Loop variable
!
      REAL  (Real64) ::                                                 &
     &    K_TOTAL
!           Total extinction including gaseous contributions
!
!     Variables related to the treatment of ill-conditioning
      REAL  (Real64) ::                                                 &
     &    EPS_R
!           The smallest real number such that 1.0-EPS_R is not 1
!           to the computer's precision
!
!     Set the tolerances used in avoiding ill-conditioning, testing
!     on any variable.
      EPS_R=EPSILON(TAU(1, 1))
!
!     The machine tolerance is added to the denominator in the
!     expression for omega to prevent division by zero: this is
!     significant only if the total extinction is small, and thus
!     will not sensibly affect any physical results.
!
      IF (I_SCATTER_METHOD_BAND == IP_SCATTER_FULL) THEN
!
        DO I=I_FIRST_LAYER, I_LAST_LAYER
          DO L=1, N_PROFILE
            K_TOTAL=K_GREY_TOT(L, I)+K_GAS_ABS(L, I)
            TAU(L, I)=K_TOTAL*D_MASS(L, I)
            IF (K_TOTAL >  0.0E+00_Real64) THEN
              OMEGA(L, I)=K_EXT_SCAT(L, I)/K_TOTAL
            ELSE
              OMEGA(L, I)=0.0E+00_Real64
            ENDIF
            OMEGA(L, I)                                                 &
     &        =MIN(OMEGA(L, I), 1.0E+00_Real64-3.2E+01_Real64*EPS_R)
          ENDDO
        ENDDO
!
      ELSE IF (I_SCATTER_METHOD_BAND == IP_NO_SCATTER_ABS) THEN
!
!       The scattering extinction is ignored completely, so
!       only the absorptive contributions to the single
!       scattering properties are included. If full scattering
!       is not to be used in the IR this is normally the appropriate
!       approximation as scattering is still dominated by the
!       forward peak.
!
        DO I=I_FIRST_LAYER, I_LAST_LAYER
          DO L=1, N_PROFILE
            TAU(L, I)=(K_GREY_TOT(L, I)+K_GAS_ABS(L, I)                 &
     &        -K_EXT_SCAT(L, I))*D_MASS(L, I)
            OMEGA(L, I)=0.0
          ENDDO
        ENDDO
!
      ELSE IF (I_SCATTER_METHOD_BAND == IP_NO_SCATTER_EXT) THEN
!
!       The scattering extinction is added on to the absorption.
!       This option is usually a bad approximation to the effects
!       of scattering in the IR, but may occasionally be appropriate
!       if the asymmetry is low.
!
        DO I=I_FIRST_LAYER, I_LAST_LAYER
          DO L=1, N_PROFILE
            TAU(L, I)=(K_GREY_TOT(L, I)+K_GAS_ABS(L, I))                &
     &        *D_MASS(L, I)
            OMEGA(L, I)=0.0E+00_Real64
          ENDDO
        ENDDO
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE SINGLE_SCATTERING
#endif
