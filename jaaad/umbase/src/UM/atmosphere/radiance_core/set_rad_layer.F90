#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the layers in which radiances are required.
!
! Purpose:
!   This determines the layers of the atmosphere where the analytic
!   expression for the radiance must be intercepted to give values
!   on the correct levels.
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
      SUBROUTINE SET_RAD_LAYER(IERR                                     &
     &  , N_LAYER, N_VIEWING_LEVEL, VIEWING_LEVEL                       &
     &  , I_RAD_LAYER, FRAC_RAD_LAYER                                   &
     &  , ND_VIEWING_LEVEL                                              &
     &  )
!
!
      IMPLICIT NONE
!
!
      INTEGER, INTENT(IN) ::                                            &
     &    ND_VIEWING_LEVEL
!           Size allocated for levels where radiances are calculated
!
!     Header files
#include "c_kinds.h"
#include "error_pcf3z.h"
#include "def_std_io_icf3z.h"
!
!
!     Dummy arguments
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
!
      INTEGER, INTENT(IN) ::                                            &
     &    N_VIEWING_LEVEL                                               &
!           Number of levels on which to calculate the radiance
     &  , N_LAYER
!           Number of atmospheric layers
      REAL  (Real64), INTENT(IN) ::                                     &
     &    VIEWING_LEVEL(ND_VIEWING_LEVEL)
!           Levels where radiances are calculated
      INTEGER, INTENT(OUT) ::                                           &
     &    I_RAD_LAYER(ND_VIEWING_LEVEL)
!           Layers in which to intercept radiances
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FRAC_RAD_LAYER(ND_VIEWING_LEVEL)
!           Fractions below the tops of the layers where radiances
!           are calculated
!
!
!     Local Variables
      INTEGER                                                           &
     &    I
!           Loop variable
      REAL  (Real64) ::                                                 &
     &    TOL_BDY
!           The tolerance detecting the closeness of boundaries
!
!
!
!     Set the tolerance for detecting boundaries.
      TOL_BDY=1.6E+01_Real64*EPSILON(TOL_BDY)
!
      DO I=1, N_VIEWING_LEVEL
!
!       Check that a level is not above the top of the atmosphere.
        IF (VIEWING_LEVEL(I) <  0.0E+00_Real64) THEN
          WRITE(IU_ERR, '(/A)')                                         &
     &      '*** Error: A viewing level is above the TOA.'
          IERR=I_ERR_FATAL
          RETURN
        ENDIF
!
        I_RAD_LAYER(I)=INT(VIEWING_LEVEL(I))+1
        FRAC_RAD_LAYER(I)=1.0E+00_Real64+VIEWING_LEVEL(I)               &
     &    -REAL(I_RAD_LAYER(I), Real64)
!
!       At the bottom of the atmosphere this will give a value greater
!       than N_LAYER, so we reset, but check that an unreasonable
!       attempt to get radiances below the column has not been made:
!       this will give a fatal error.
        IF (I_RAD_LAYER(I) >  N_LAYER) THEN
          IF (FRAC_RAD_LAYER(I) <  TOL_BDY) THEN
            I_RAD_LAYER(I)=I_RAD_LAYER(I)-1
            FRAC_RAD_LAYER(I)=1.0E+00_Real64
          ELSE
            WRITE(IU_ERR, '(/A)')                                       &
     &        '*** Error: A viewing level is below the surface.'
            IERR=I_ERR_FATAL
            RETURN
          ENDIF
        ENDIF
!
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SET_RAD_LAYER
#endif
