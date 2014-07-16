#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set geometry of clouds.
!
! Method:
!       For use in multi-column mode arrays are set for each layer
!       pointing to profiles which have non-negligible clear or
!       cloudy fractions. The topmost cloudy layers are also
!       detected.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_CLOUD_GEOMETRY(N_PROFILE, N_LAYER                  &
     &  , L_GLOBAL_CLOUD_TOP, N_CLOUD_TOP_GLOBAL, W_CLOUD               &
     &  , N_CLOUD_TOP, N_CLOUD_PROFILE, I_CLOUD_PROFILE                 &
     &  , ND_PROFILE, ND_LAYER, ID_CT                                   &
     &  )
!
!
!
      IMPLICIT NONE
!
! Include header Files
#include "c_kinds.h"
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Maximum number of profiles
     &  , ND_LAYER                                                      &
!           Maximum number of layers
     &  , ID_CT
!           Topmost declared cloudy layer
!
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of atmospheric profiles
     &  , N_LAYER
!           Number of layers
!
      LOGICAL, INTENT(IN) ::                                            &
     &    L_GLOBAL_CLOUD_TOP
!           Flag to use a global value for the topmost cloudy layer
      INTEGER, INTENT(IN) ::                                            &
     &    N_CLOUD_TOP_GLOBAL
!           Global topmost cloudy layer
      REAL  (Real64), INTENT(IN) ::                                     &
     &    W_CLOUD(ND_PROFILE, ID_CT: ND_LAYER)
!           Amounts of cloud
!
      INTEGER, INTENT(OUT) ::                                           &
     &    N_CLOUD_TOP                                                   &
!           Topmost cloudy layer
     &  , N_CLOUD_PROFILE(ID_CT: ND_LAYER)                              &
!           Number of cloudy profiles
     &  , I_CLOUD_PROFILE(ND_PROFILE, ID_CT: ND_LAYER)
!           Profiles containing clouds
!
!
!     Local variables.
      INTEGER                                                           &
     &    L                                                             &
!           Loop variable
     &  , I
!           Loop variable
!
!
!
      DO I=ID_CT, N_LAYER
        N_CLOUD_PROFILE(I)=0
        DO L=1, N_PROFILE
          IF (W_CLOUD(L, I) >  0.0E+00_Real64) THEN
            N_CLOUD_PROFILE(I)=N_CLOUD_PROFILE(I)+1
            I_CLOUD_PROFILE(N_CLOUD_PROFILE(I), I)=L
          ENDIF
        ENDDO
      ENDDO
!
      IF (L_GLOBAL_CLOUD_TOP) THEN
        N_CLOUD_TOP=N_CLOUD_TOP_GLOBAL
      ELSE
        N_CLOUD_TOP=ID_CT
        DO WHILE ( (N_CLOUD_TOP <  N_LAYER).AND.                        &
     &             (N_CLOUD_PROFILE(N_CLOUD_TOP) == 0) )
          N_CLOUD_TOP=N_CLOUD_TOP+1
        ENDDO
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE SET_CLOUD_GEOMETRY
#endif
