#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to solve for triple overlaps with approximate scattering.
!
! Method:
!       The flux is propagated downwards, ignoring reflection terms.
!       since the routine uses differential fluxes, this effectively
!       treats the upward flux as Planckian at this point. Upward
!       fluxes are calculated using the newly available approximate
!       downward fluxes in the reflected terms.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLVER_TRIPLE_APP_SCAT(N_PROFILE, N_LAYER, N_CLOUD_TOP &
     &   , T, R, S_DOWN, S_UP                                           &
     &   , T_STRAT, R_STRAT, S_DOWN_STRAT, S_UP_STRAT                   &
     &   , T_CONV, R_CONV, S_DOWN_CONV, S_UP_CONV                       &
     &   , V11, V12, V13, V21, V22, V23, V31, V32, V33                  &
     &   , U11, U12, U13, U21, U22, U23, U31, U32, U33                  &
     &   , FLUX_INC_DOWN                                                &
     &   , SOURCE_GROUND_FREE, SOURCE_GROUND_STRAT                      &
     &   , SOURCE_GROUND_CONV, ALBEDO_SURFACE_DIFF                      &
     &   , FLUX_TOTAL                                                   &
     &   , ND_PROFILE, ND_LAYER, ID_CT                                  &
     &   )
!
!
      IMPLICIT NONE
!
!
! Include Header files
#include "c_kinds.h"
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ID_CT
!           Topmost declared cloudy layer
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_CLOUD_TOP
!           Topmost cloudy layer
      REAL  (Real64), INTENT(IN) ::                                     &
     &    T(ND_PROFILE, ND_LAYER)                                       &
!           Clear-sky transmission
     &  , R(ND_PROFILE, ND_LAYER)                                       &
!           Clear-sky reflection
     &  , S_DOWN(ND_PROFILE, ND_LAYER)                                  &
!           Clear-sky downward source function
     &  , S_UP(ND_PROFILE, ND_LAYER)                                    &
!           Clear-sky upward source function
     &  , T_STRAT(ND_PROFILE, ND_LAYER)                                 &
!           Stratfiform transmission
     &  , R_STRAT(ND_PROFILE, ND_LAYER)                                 &
!           Stratfiform reflection
     &  , S_DOWN_STRAT(ND_PROFILE, ND_LAYER)                            &
!           Downward stratfiform source function
     &  , S_UP_STRAT(ND_PROFILE, ND_LAYER)                              &
!           Upward stratfiform source function
     &  , T_CONV(ND_PROFILE, ND_LAYER)                                  &
!           Convective transmission
     &  , R_CONV(ND_PROFILE, ND_LAYER)                                  &
!           Convective reflection
     &  , S_DOWN_CONV(ND_PROFILE, ND_LAYER)                             &
!           Downward convective source function
     &  , S_UP_CONV(ND_PROFILE, ND_LAYER)
!           Upward convective source function
      REAL  (Real64), INTENT(IN) ::                                     &
     &    V11(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for downward radiation
     &  , V12(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for downward radiation
     &  , V13(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for downward radiation
     &  , V21(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for downward radiation
     &  , V22(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for downward radiation
     &  , V23(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for downward radiation
     &  , V31(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for downward radiation
     &  , V32(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for downward radiation
     &  , V33(ND_PROFILE, ID_CT-1: ND_LAYER)
!           Energy transfer coefficient for downward radiation
      REAL  (Real64) ::                                                 &
     &    U11(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for upward radiation
     &  , U12(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for upward radiation
     &  , U13(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for upward radiation
     &  , U21(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for upward radiation
     &  , U22(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for upward radiation
     &  , U23(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for upward radiation
     &  , U31(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for upward radiation
     &  , U32(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient for upward radiation
     &  , U33(ND_PROFILE, ID_CT-1: ND_LAYER)
!           Energy transfer coefficient for upward radiation
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FLUX_INC_DOWN(ND_PROFILE)                                     &
!           Incident flux
     &  , SOURCE_GROUND_FREE(ND_PROFILE)                                &
!           Source from ground (clear sky)
     &  , SOURCE_GROUND_STRAT(ND_PROFILE)                               &
!           Source from ground (cloudy region)
     &  , SOURCE_GROUND_CONV(ND_PROFILE)                                &
!           Source from ground (cloudy region)
     &  , ALBEDO_SURFACE_DIFF(ND_PROFILE)
!           Diffuse albedo
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_TOTAL(ND_PROFILE, 2*ND_LAYER+2)
!           Total flux
!
!     Local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L
!           Loop variable
!
!
!     Temporary fluxes
      REAL  (Real64) ::                                                 &
     &    FLUX_DOWN_1(ND_PROFILE, 0: ND_LAYER)                          &
!           Downward fluxes outside clouds just below i'th level
     &  , FLUX_DOWN_2(ND_PROFILE, 0: ND_LAYER)                          &
!           Downward fluxes inside clouds just below i'th level
     &  , FLUX_DOWN_3(ND_PROFILE, 0: ND_LAYER)                          &
!           Downward fluxes inside clouds just below i'th level
     &  , FLUX_UP_1(ND_PROFILE, 0: ND_LAYER)                            &
!           Upward fluxes outside clouds just above i'th level
     &  , FLUX_UP_2(ND_PROFILE, 0: ND_LAYER)                            &
!           Upward fluxes inside clouds just above i'th level
     &  , FLUX_UP_3(ND_PROFILE, 0: ND_LAYER)                            &
!           Upward fluxes inside clouds just above i'th level
     &  , FLUX_PROPAG_1(ND_PROFILE)                                     &
!           Temporary fluxes for propagation across layers
     &  , FLUX_PROPAG_2(ND_PROFILE)                                     &
!           Temporary fluxes for propagation across layers
     &  , FLUX_PROPAG_3(ND_PROFILE)
!           Temporary fluxes for propagation across layers
!
!
!
!
!     The arrays flux_down and flux_up will eventually contain the total
!     fluxes, but initially they are used for the clear fluxes.
!     Note that downward fluxes refer to values just below the interface
!     and upward fluxes to values just above it.
!
!
!     Downward flux:
!
!     Region above clouds:
      DO L=1, N_PROFILE
        FLUX_TOTAL(L, 2)=FLUX_INC_DOWN(L)
      ENDDO
      DO I=1, N_CLOUD_TOP-1
        DO L=1, N_PROFILE
          FLUX_TOTAL(L, 2*I+2)=T(L, I)*FLUX_TOTAL(L, 2*I)               &
     &      +S_DOWN(L, I)
        ENDDO
      ENDDO
!
!     Pass into the cloudy region. here, downward fluxes hold values
!     just below the level and upward fluxes the values just above it.
!     Thus the fluxes impinging on the layer are held.
      I=N_CLOUD_TOP-1
      DO L=1, N_PROFILE
        FLUX_DOWN_1(L, I)=V11(L, I)*FLUX_TOTAL(L, 2*I+2)
        FLUX_DOWN_2(L, I)=V21(L, I)*FLUX_TOTAL(L, 2*I+2)
        FLUX_DOWN_3(L, I)=V31(L, I)*FLUX_TOTAL(L, 2*I+2)
      ENDDO
!
      DO I=N_CLOUD_TOP, N_LAYER-1
        DO L=1, N_PROFILE
!
!         Propagte the flux across the layer.
          FLUX_PROPAG_1(L)=T(L, I)*FLUX_DOWN_1(L, I-1)                  &
     &      +S_DOWN(L, I)
          FLUX_PROPAG_2(L)=T_STRAT(L, I)*FLUX_DOWN_2(L, I-1)            &
     &      +S_DOWN_STRAT(L, I)
          FLUX_PROPAG_3(L)=T_CONV(L, I)*FLUX_DOWN_3(L, I-1)             &
     &      +S_DOWN_CONV(L, I)
!
!         Transfer across the interface.
          FLUX_DOWN_1(L, I)=V11(L, I)*FLUX_PROPAG_1(L)                  &
     &      +V12(L, I)*FLUX_PROPAG_2(L)                                 &
     &      +V13(L, I)*FLUX_PROPAG_3(L)
          FLUX_DOWN_2(L, I)=V21(L, I)*FLUX_PROPAG_1(L)                  &
     &      +V22(L, I)*FLUX_PROPAG_2(L)                                 &
     &      +V23(L, I)*FLUX_PROPAG_3(L)
          FLUX_DOWN_3(L, I)=V31(L, I)*FLUX_PROPAG_1(L)                  &
     &      +V32(L, I)*FLUX_PROPAG_2(L)                                 &
     &      +V33(L, I)*FLUX_PROPAG_3(L)
!
        ENDDO
      ENDDO
!
!     Propagate across the bottom layer and form the reflected beam.
!     We do not transfer fluxes across the bottom interface, so as
!     to make the reflection consistent between regions.
      DO L=1, N_PROFILE
!
!       Propagte the flux through the layer.
        FLUX_DOWN_1(L, N_LAYER)                                         &
     &    =T(L, N_LAYER)*FLUX_DOWN_1(L, N_LAYER-1)                      &
     &    +S_DOWN(L, N_LAYER)
        FLUX_DOWN_2(L, N_LAYER)                                         &
     &    =T_STRAT(L, N_LAYER)*FLUX_DOWN_2(L, N_LAYER-1)                &
     &    +S_DOWN_STRAT(L, N_LAYER)
        FLUX_DOWN_3(L, N_LAYER)                                         &
     &    =T_CONV(L, N_LAYER)*FLUX_DOWN_3(L, N_LAYER-1)                 &
     &    +S_DOWN_CONV(L, N_LAYER)
!
!       Reflect from the surface.
        FLUX_UP_1(L, N_LAYER)                                           &
     &    =ALBEDO_SURFACE_DIFF(L)*FLUX_DOWN_1(L, N_LAYER)               &
     &    +SOURCE_GROUND_FREE(L)
        FLUX_UP_2(L, I)                                                 &
     &    =ALBEDO_SURFACE_DIFF(L)*FLUX_DOWN_2(L, N_LAYER)               &
     &    +SOURCE_GROUND_STRAT(L)
        FLUX_UP_3(L, I)                                                 &
     &    =ALBEDO_SURFACE_DIFF(L)*FLUX_DOWN_3(L, N_LAYER)               &
     &    +SOURCE_GROUND_CONV(L)
!
!       Propagate across the bottom layer.
        FLUX_PROPAG_1(L)                                                &
     &    =T(L, N_LAYER)*FLUX_UP_1(L, N_LAYER)+S_UP(L, N_LAYER)         &
     &    +R(L, N_LAYER)*FLUX_DOWN_1(L, N_LAYER-1)
        FLUX_PROPAG_2(L)                                                &
     &    =T_STRAT(L, N_LAYER)*FLUX_UP_2(L, N_LAYER)                    &
     &    +S_UP_STRAT(L, N_LAYER)                                       &
     &    +R_STRAT(L, N_LAYER)*FLUX_DOWN_2(L, N_LAYER-1)
        FLUX_PROPAG_3(L)                                                &
     &    =T_CONV(L, N_LAYER)*FLUX_UP_3(L, N_LAYER)                     &
     &    +S_UP_CONV(L, N_LAYER)                                        &
     &    +R_CONV(L, N_LAYER)*FLUX_DOWN_3(L, N_LAYER-1)
!
      ENDDO
!
!
!
!     Work back up through the column assigning the upward fluxes.
      DO I=N_LAYER-1, N_CLOUD_TOP, -1
        DO L=1, N_PROFILE
!
          FLUX_UP_1(L, I)=U11(L, I)*FLUX_PROPAG_1(L)                    &
     &      +U12(L, I)*FLUX_PROPAG_2(L)                                 &
     &      +U13(L, I)*FLUX_PROPAG_3(L)
          FLUX_UP_2(L, I)=U21(L, I)*FLUX_PROPAG_1(L)                    &
     &      +U22(L, I)*FLUX_PROPAG_2(L)                                 &
     &      +U23(L, I)*FLUX_PROPAG_3(L)
          FLUX_UP_3(L, I)=U31(L, I)*FLUX_PROPAG_1(L)                    &
     &      +U32(L, I)*FLUX_PROPAG_2(L)                                 &
     &      +U33(L, I)*FLUX_PROPAG_3(L)
!
          FLUX_PROPAG_1(L)=T(L, I)*FLUX_UP_1(L, I)+S_UP(L, I)           &
     &      +R(L, I)*FLUX_DOWN_1(L, I-1)
          FLUX_PROPAG_2(L)=T_STRAT(L, I)*FLUX_UP_2(L, I)                &
     &      +S_UP_STRAT(L, I)+R_STRAT(L, I)*FLUX_DOWN_2(L, I-1)
          FLUX_PROPAG_3(L)=T_CONV(L, I)*FLUX_UP_3(L, I)                 &
     &      +S_UP_CONV(L, I)+R_CONV(L, I)*FLUX_DOWN_3(L, I-1)
!
        ENDDO
      ENDDO
!
!     Propagate into the cloud-free region.
      I=N_CLOUD_TOP-1
      DO L=1, N_PROFILE
        FLUX_TOTAL(L, 2*I+1)=FLUX_PROPAG_1(L)+FLUX_PROPAG_2(L)          &
     &    +FLUX_PROPAG_3(L)
      ENDDO
!
!     Continue through the upper cloudy layers.
      DO I=N_CLOUD_TOP-1, 1, -1
        DO L=1, N_PROFILE
          FLUX_TOTAL(L, 2*I-1)=T(L, I)*FLUX_TOTAL(L, 2*I+1)             &
     &      +R(L, I)*FLUX_TOTAL(L, 2*I)+S_UP(L, I)
        ENDDO
      ENDDO
!
!     Assign the total fluxes on the intermediate cloudy layers.
      DO I=N_CLOUD_TOP, N_LAYER
        DO L=1, N_PROFILE
          FLUX_TOTAL(L, 2*I+1)=FLUX_UP_1(L, I)+FLUX_UP_2(L, I)          &
     &      +FLUX_UP_3(L, I)
          FLUX_TOTAL(L, 2*I+2)=FLUX_DOWN_1(L, I)+FLUX_DOWN_2(L, I)      &
     &      +FLUX_DOWN_3(L, I)
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SOLVER_TRIPLE_APP_SCAT
#endif
