#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the pentadiagonal matrix for the fluxes.
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
      SUBROUTINE SET_MATRIX_PENTADIAGONAL(N_PROFILE, N_LAYER            &
     &   , TRANS, REFLECT                                               &
     &   , S_DOWN, S_UP                                                 &
     &   , DIFFUSE_ALBEDO, DIRECT_ALBEDO                                &
     &   , FLUX_DIRECT_GROUND, FLUX_INC_DOWN                            &
     &   , D_PLANCK_FLUX_SURFACE                                        &
     &   , A5, B                                                        &
     &   , ND_PROFILE, ND_LAYER                                         &
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
     &    ND_PROFILE                                                    &
!           SIze allocated for atmospheric profiles
     &  , ND_LAYER
!           SIze allocated for atmospheric layers
!
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER
!           Number of layers
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TRANS(ND_PROFILE, ND_LAYER)                                   &
!           Transmission coefficient
     &  , REFLECT(ND_PROFILE, ND_LAYER)                                 &
!           Reflection coefficient
     &  , S_DOWN(ND_PROFILE, ND_LAYER)                                  &
!           Downward diffuse source
     &  , S_UP(ND_PROFILE, ND_LAYER)                                    &
!           Upward diffuse source
     &  , DIFFUSE_ALBEDO(ND_PROFILE)                                    &
!           Diffuse surface albedo
     &  , DIRECT_ALBEDO(ND_PROFILE)                                     &
!           Direct surface albedo
     &  , D_PLANCK_FLUX_SURFACE(ND_PROFILE)                             &
!           Difference in Planckian fluxes at the surface
     &  , FLUX_INC_DOWN(ND_PROFILE)                                     &
!           Incident total flux
     &  , FLUX_DIRECT_GROUND(ND_PROFILE)
!           Direct flux at ground level
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    A5(ND_PROFILE, 5, 2*ND_LAYER+2)                               &
!           Pentadiagonal matrix
     &  , B(ND_PROFILE, 2*ND_LAYER+2)
!           Source terms for equations
!
!     Declaration of local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L
!           Loop variable
!
!
!
!     The top boundary condition:
      DO L=1, N_PROFILE
        A5(L, 4, 2)=0.0E+00_Real64
        A5(L, 3, 2)=1.0E+00_Real64
        A5(L, 2, 2)=0.0E+00_Real64
        A5(L, 1, 2)=0.0E+00_Real64
        B(L, 2)=FLUX_INC_DOWN(L)
      ENDDO
!
!     Interior rows: odd and even rows correspond to different boundary
!     conditions.
      DO I=1, N_LAYER
        DO L=1, N_PROFILE
!
          A5(L, 5, 2*I-1)=0.0E+00_Real64
          A5(L, 4, 2*I-1)=0.0E+00_Real64
          A5(L, 3, 2*I-1)=-1.0E+00_Real64
          A5(L, 2, 2*I-1)=REFLECT(L, I)
          A5(L, 1, 2*I-1)=TRANS(L, I)
          B(L, 2*I-1)=-S_UP(L, I)
!
          A5(L, 5, 2*I+2)=TRANS(L, I)
          A5(L, 4, 2*I+2)=REFLECT(L, I)
          A5(L, 3, 2*I+2)=-1.0E+00_Real64
          A5(L, 2, 2*I+2)=0.0E+00_Real64
          A5(L, 1, 2*I+2)=0.0E+00_Real64
          B(L, 2*I+2)=-S_DOWN(L, I)
!
        ENDDO
      ENDDO
!
!     The surface boundary condition:
      DO L=1, N_PROFILE
        A5(L, 5, 2*N_LAYER+1)=0.0E+00_Real64
        A5(L, 4, 2*N_LAYER+1)=0.0E+00_Real64
        A5(L, 3, 2*N_LAYER+1)=1.0E+00_Real64
        A5(L, 2, 2*N_LAYER+1)=-DIFFUSE_ALBEDO(L)
        B(L, 2*N_LAYER+1)                                               &
     &    =(1.0E+00_Real64-DIFFUSE_ALBEDO(L))*D_PLANCK_FLUX_SURFACE(L)  &
     &    +(DIRECT_ALBEDO(L)-DIFFUSE_ALBEDO(L))                         &
     &    *FLUX_DIRECT_GROUND(L)
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SET_MATRIX_PENTADIAGONAL
#endif
