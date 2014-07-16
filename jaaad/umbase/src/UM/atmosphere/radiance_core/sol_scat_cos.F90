#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate solar scattering angles.
!
! Purpose:
!   This routine returns the cosines of the angles of scattering
!   from the solar beam for each viewing direction.
!
! Method:
!   A scalar product of the solar and viewing directions is
!   evaluated. This routine is called only when radiances are to
!   be calculated, so ND_PROFILE can be used for all horizontal
!   dimensions.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOL_SCAT_COS(N_PROFILE, N_DIRECTION                    &
     &  , MU_0, DIRECTION, COS_SOL_VIEW                                 &
     &  , ND_PROFILE, ND_DIRECTION)
!
!
!
      IMPLICIT NONE
!
! Include header files
#include "c_kinds.h"
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_DIRECTION
!           Size allocated for viewing directions
!
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of atmospheric profiles
     &  , N_DIRECTION
!           Number of viewing directions
      REAL  (Real64), INTENT(IN) ::                                     &
     &    MU_0(ND_PROFILE)                                              &
!           Cosines of solar zenith angles
     &  , DIRECTION(ND_PROFILE, ND_DIRECTION, 2)
!           Viewing directions stored as the cosine of the polar
!           viewing angle and the azimuthal viewing angle itself
!           realative to the solar direction
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    COS_SOL_VIEW(ND_PROFILE, ND_DIRECTION)
!           Cosines of the angles between the solar beam and the
!           viewing directions
!
!
!     Local variables
      INTEGER                                                           &
     &    ID                                                            &
!           Loop variable
     &  , L
!           Loop variable
!
!
!
      DO ID=1, N_DIRECTION
        DO L=1, N_PROFILE
          COS_SOL_VIEW(L, ID)=-MU_0(L)*DIRECTION(L, ID, 1)              &
     &      +SQRT((1.0E+00_Real64-MU_0(L)*MU_0(L))                      &
     &      *(1.0E+00_Real64-DIRECTION(L, ID, 1)*DIRECTION(L, ID, 1)))  &
     &      *COS(DIRECTION(L, ID, 2))
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SOL_SCAT_COS
#endif
