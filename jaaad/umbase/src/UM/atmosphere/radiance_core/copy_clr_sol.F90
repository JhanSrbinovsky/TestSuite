#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set clear-sky solar phase function.
!
! Method:
!       The arrays of clear-sky forward scattering and the solar
!       phase function at the top of the column and of these
!       same properties from the total list lower down
!       are combined to give unified arrays of clear-sky
!       optical properties.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE COPY_CLR_SOL(N_PROFILE, N_LAYER, N_CLOUD_TOP           &
     &  , N_DIRECTION, L_RESCALE                                        &
     &  , FORWARD_SCATTER_CLR, PHASE_FNC_SOLAR_CLR                      &
     &  , FORWARD_SCATTER, PHASE_FNC_SOLAR                              &
     &  , FORWARD_SCATTER_CLR_F                                         &
     &  , PHASE_FNC_SOLAR_CLR_F                                         &
!                       Sizes of arrays
     &  , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_DIRECTION       &
     &  )
!
!
      IMPLICIT NONE
!
! Include header files
#include "c_kinds.h"
!
!     Sizes of arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for profiles
     &  , ND_LAYER                                                      &
!           Size allocated for layers
     &  , ND_LAYER_CLR                                                  &
!           Size allocated for totally clear layers
     &  , ID_CT                                                         &
!           Topmost declared layer for cloudy optical properties
     &  , ND_DIRECTION
!           SIze allocated for viewing directions
!
!                     Atmospheric properties
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of atmospheric layers
     &  , N_CLOUD_TOP                                                   &
!           Topmost cloudy layer
     &  , N_DIRECTION
!           Number of terms in the phase function
      LOGICAL, INTENT(IN) ::                                            &
     &    L_RESCALE
!           Flag for rescaling of the optical properties
!
!                     Optical properties
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FORWARD_SCATTER_CLR(ND_PROFILE, ND_LAYER_CLR)                 &
!           Forward scattering in the totally clear region
     &  , PHASE_FNC_SOLAR_CLR(ND_PROFILE, ND_LAYER_CLR, ND_DIRECTION)
!           Phase function in totally clear region
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FORWARD_SCATTER(ND_PROFILE, ID_CT: ND_LAYER)                  &
!           Forward scattering in the cloudy regions
     &  , PHASE_FNC_SOLAR(ND_PROFILE, ID_CT: ND_LAYER, ND_DIRECTION)
!           Phase function restricted to clear-sky regions
!
!                     Single scattering properties
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FORWARD_SCATTER_CLR_F(ND_PROFILE, ND_LAYER)                   &
!           Forward scattering expanded to the whole column
     &  , PHASE_FNC_SOLAR_CLR_F(ND_PROFILE, ND_LAYER, ND_DIRECTION)
!           Phase function expanded to the whole column
!
!
!
!     Local variables.
      INTEGER                                                           &
     &    L                                                             &
!           Loop variable
     &  , I                                                             &
!           Loop variable
     &  , K
!           Loop variable
!
!
!
!     Above cloud top.
      DO I=1, N_CLOUD_TOP-1
        IF (L_RESCALE) THEN
          DO L=1, N_PROFILE
            FORWARD_SCATTER_CLR_F(L, I)=FORWARD_SCATTER_CLR(L, I)
          ENDDO
        ENDIF
        DO K=1, N_DIRECTION
          DO L=1, N_PROFILE
            PHASE_FNC_SOLAR_CLR_F(L, I, K)=PHASE_FNC_SOLAR_CLR(L, I, K)
          ENDDO
        ENDDO
      ENDDO
!
!     Below cloud top.
      DO I=N_CLOUD_TOP, N_LAYER
        IF (L_RESCALE) THEN
          DO L=1, N_PROFILE
            FORWARD_SCATTER_CLR_F(L, I)=FORWARD_SCATTER(L, I)
          ENDDO
        ENDIF
        DO K=1, N_DIRECTION
          DO L=1, N_PROFILE
            PHASE_FNC_SOLAR_CLR_F(L, I, K)=PHASE_FNC_SOLAR(L, I, K)
          ENDDO
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE COPY_CLR_SOL
#endif
