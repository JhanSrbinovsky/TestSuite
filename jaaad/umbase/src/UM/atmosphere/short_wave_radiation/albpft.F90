#if defined(A01_3A) || defined(A01_3C) \
 || defined(A01_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine to calculate the spectral albedo of the land surface using
! the two stream approach of Sellers, 1995.
!
!**********************************************************************
      SUBROUTINE ALBPFT (P_FIELD,LAND_FIELD,                            &
     &                       LAND_INDEX,TILE_INDEX,TILE_PTS,ILAYERS,    &
     &                       ALBSOIL,COSZ,LAI,ALB_TYPE,                 &
     &                       FAPAR_DIR,FAPAR_DIF,CAN_RAD_MOD)


      IMPLICIT NONE

#include "nstypes.h"

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     & P_FIELD                                                          &
                                   ! Total number of grid points.
     &,LAND_FIELD                  ! Number of land points.

!   Array arguments with intent(in):
      INTEGER                                                           &
     & LAND_INDEX(LAND_FIELD)                                           &
                                   ! Index of land points.
     &,TILE_PTS(NTYPE)                                                  &
                                   ! Number of land points which
!                                  ! include the surface type.
     &,TILE_INDEX(LAND_FIELD,NTYPE)                                     &
                                   ! Indices of land points which
!                                  ! include the surface type.

     &,ILAYERS                                                          &
                                   ! IN Number of layers over which
                                   !    the PAR absorption profile is to
                                   !    be calculated.
     &,CAN_RAD_MOD                 ! IN Which canopy radiation model is 
                                   ! used (1/2)
      REAL                                                              &
     & ALBSOIL(LAND_FIELD)                                              &
                                   ! Soil albedo.
     &,COSZ(P_FIELD)                                                    &
                                   ! Cosine of the zenith angle.
     &,LAI(LAND_FIELD,NPFT)        ! Leaf area index.

!   Array arguments with intent(out):
      REAL                                                              &
     & ALB_TYPE(LAND_FIELD,NTYPE,4)! Albedos for surface types.
!                                  !   (*,*,1) - Direct beam visible
!                                  !   (*,*,2) - Diffuse visible
!                                  !   (*,*,3) - Direct beam near-IR
!                                  !   (*,*,4) - Diffuse near-IR

! Local arrays:
      REAL                                                              &
     & ALBUDIF(LAND_FIELD,2)                                            &
                                  ! Diffuse albedo of the underlying
!                                 ! surface.
     &,ALBUDIR(LAND_FIELD,2)                                            &
                                  ! Direct albedo of the underlying
!                                 ! surface.
     &,FAPAR_DIR(LAND_FIELD,NPFT,ILAYERS)                               &
!                                  ! OUT Profile of absorbed PAR
!                                  !     - Direct (fraction/LAI)
     &,FAPAR_DIF(LAND_FIELD,NPFT,ILAYERS)                               &
!                                  ! OUT Profile of absorbed PAR
!                                  !     - Diffuse (fraction/LAI)
     &,ALPL(2)                                                          &
                                  ! Leaf reflection coefficient.
     &,OM(2)                                                            &
                                  ! Leaf scattering coefficient.
     &,TAUL(2)                    ! Leaf transmission coefficient.

! Local scalars:
      REAL                                                              &
     & BETADIR                                                          &
                                  ! Upscatter parameter for direct beam.
     &,BETADIF                                                          &
                                  ! Upscatter parameter for diffuse beam
     &,COSZM                                                            &
                                  ! Mean value of COSZ.
     &,K                                                                &
                                  ! Optical depth per unit leaf area.
     &,G                                                                &
                                  ! Relative projected leaf area in
                                  ! direction cosz.
     &,SALB                                                             &
                                  ! Single scattering albedo.
     &,SQCOST                                                           &
                                  ! Cosine squared of the mean leaf
                                  ! angle to the horizontal.
     &,B,C,CA,D,F,H,U1                                                  &
                                  ! Miscellaneous variables from
     &,P1,P2,P3,P4,D1                                                   &
                                  ! Sellers (1985).
     &,H1,H2,H3,H7,H8                                                   &
                                  !
     &,S1,S2,SIG                  !

!-----------------------------------------------------------------------
! Additional work variables to calculate profiles of absorbed PAR
!-----------------------------------------------------------------------
      REAL                                                              &
     & DLAI                                                             &
                                      ! LAI increment.
     &,LA                                                               &
                                      ! Cumulative LAI through canopy.
     &,DRDIRD_DLAI,DRDIRU_DLAI                                          &
     &,DRDIFD_DLAI,DRDIFU_DLAI                                          &
     &,U2,U3,D2,H4,H5,H6,H9,H10
                                      ! Rate of change of fluxes with LAI
                                      ! (W/m2/LAI). 

      INTEGER                                                           &
     & BAND,I,J,L,N               ! Loop counters.

#include "trif.h"

      DO L=1,LAND_FIELD
        ALBUDIF(L,1) = ALBSOIL(L)
        ALBUDIF(L,2) = ALBSOIL(L)
        ALBUDIR(L,1) = ALBSOIL(L)
        ALBUDIR(L,2) = ALBSOIL(L)
      ENDDO

      DO N=1,NPFT

        OM(1) = OMEGA(N)
        OM(2) = OMNIR(N)
        ALPL(1) = ALPAR(N)
        ALPL(2) = ALNIR(N)

        DO BAND=1,2  ! Visible and near-IR bands
          TAUL(BAND) = OM(BAND) - ALPL(BAND)
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            I = LAND_INDEX(L)
            IF (ORIENT(N) == 0) THEN
              SQCOST = 1./3.
              G = 0.5
              COSZM = 1.0
              IF (CAN_RAD_MOD == 2) THEN
                K = G / 0.01
              ENDIF
              SALB = 0.5*OM(BAND)
              IF (COSZ(I) > 0.01) THEN
                IF (CAN_RAD_MOD == 2) THEN
                  K = G / COSZ(I)
              ENDIF
              SALB = 0.5*OM(BAND) *                                     &
     &                 ( 1. - COSZ(I)*LOG((COSZ(I)+1.)/COSZ(I)) )
          ENDIF
          ELSEIF (ORIENT(N) == 1) THEN
            SQCOST = 1.
            IF (CAN_RAD_MOD == 1) THEN
              G = COSZ(I)
            ENDIF
            COSZM = 1.
            IF (CAN_RAD_MOD == 2) THEN
!             K = G / COSZ(I) = 1.0 because G = COSZ(I)
!                                  for horizontal leaves
              K = 1.0
            ENDIF

            SALB = OM(BAND)/4.0
          ENDIF
          IF (CAN_RAD_MOD == 1) THEN
            K = G / 0.01
            IF (COSZ(I) >  0.01) K = G / COSZ(I)
          ENDIF

          BETADIR = (1. + COSZM*K)/(OM(BAND)*COSZM*K)*SALB
          C = 0.5*( ALPL(BAND) + TAUL(BAND) +                           &
     &               (ALPL(BAND) - TAUL(BAND))*SQCOST )
          BETADIF = C / OM(BAND)
          B = 1. - (1. - BETADIF)*OM(BAND)
          D = OM(BAND)*COSZM*K*BETADIR
          F = OM(BAND)*COSZM*K*(1. - BETADIR)
          H = SQRT(B*B - C*C) / COSZM
          SIG = (COSZM*K)**2 + C*C - B*B
          U1 = B - C/ALBUDIF(L,BAND)
          CA = C*ALBUDIR(L,BAND)/ALBUDIF(L,BAND)
          S1 = EXP(-H*LAI(L,N))
          S2 = EXP(-K*LAI(L,N))
          P1 = B + COSZM*H
          P2 = B - COSZM*H
          P3 = B + COSZM*K
          P4 = B - COSZM*K
          D1 = P1*(U1 - COSZM*H)/S1 - P2*(U1 + COSZM*H)*S1
          H1 = -D*P4 - C*F
          H2 = ( (D - P3*H1/SIG) * (U1 - COSZM*H) / S1 -                &
     &             (D - CA - (U1 + COSZM*K)*H1/SIG)*P2*S2 ) / D1
          H3 = - ( (D - P3*H1/SIG) * (U1 + COSZM*H)*S1 -                &
     &               (D - CA - (U1 + COSZM*K)*H1/SIG)*P1*S2 ) / D1
          H7 = (C/D1)*(U1 - COSZM*H) / S1
          H8 = - (C/D1)*(U1 + COSZM*H) * S1
          ALB_TYPE(L,N,2*BAND-1) = H1/SIG + H2 + H3   ! Direct beam
          ALB_TYPE(L,N,2*BAND) = H7 + H8              ! Diffuse
          IF (CAN_RAD_MOD == 2) THEN
!-----------------------------------------------------------------------
! If required calculate the profile of absorbed PAR through the canopy
! of direct and diffuse beams (BEWARE: assumes PAR is band 1)
!-----------------------------------------------------------------------
            IF (BAND==1) THEN
              U2 = B-C*ALBUDIF(L,BAND)
              U3 = F+C*ALBUDIF(L,BAND)
              D2 = (U2+COSZM*H)/S1-(U2-COSZM*H)*S1
              H4 = -F*P3-C*D
              H5 = -1./D2*(H4/SIG*(U2+COSZM*H)/S1                       &
     &           +(U3-H4/SIG*(U2-COSZM*K))*S2)
              H6 = 1./D2*(H4/SIG*(U2-COSZM*H)*S1                        &
     &          +(U3-H4/SIG*(U2-COSZM*K))*S2)
              H9 = 1./D2*(U2+COSZM*H)/S1
              H10 = -1./D2*(U2-COSZM*H)*S1
!-----------------------------------------------------------------------
! Two-stream equations for direct and diffuse upward and downward beams:
!             RDIRU(I)=H1/SIG*EXP(-K*LA)+H2*EXP(-H*LA)+H3*EXP(H*LA)
!             RDIRD(I)=(H4/SIG+1)*EXP(-K*LA)+H5*EXP(-H*LA)+H6*EXP(H*LA)
!             RDIFU(I)=H7*EXP(-H*LA)+H8*EXP(H*LA)
!             RDIFD(I)=H9*EXP(-H*LA)+H10*EXP(H*LA)
! Differentiate these equations to calculate PAR absorption per unit
! LAI down through the canopy. Centre derivatives in the centre of each
! LAI layer.
!-----------------------------------------------------------------------
              DLAI=LAI(L,N)/REAL(ILAYERS)
              LA=0.5*DLAI
              DO I=1,ILAYERS

                DRDIRU_DLAI=-K*H1/SIG*EXP(-K*LA)-H*H2*EXP(-H*LA)        &
     &                    +H*H3*EXP(H*LA)
                DRDIRD_DLAI=-K*(H4/SIG+1)*EXP(-K*LA)-H*H5*EXP(-H*LA)    &
     &                    +H*H6*EXP(H*LA)
                DRDIFU_DLAI=-H*H7*EXP(-H*LA)+H*H8*EXP(H*LA)
                DRDIFD_DLAI=-H*H9*EXP(-H*LA)+H*H10*EXP(H*LA)
  
                FAPAR_DIR(L,N,I)=-DRDIRD_DLAI+DRDIRU_DLAI
                FAPAR_DIF(L,N,I)=-DRDIFD_DLAI+DRDIFU_DLAI
                LA=LA+DLAI


              ENDDO  !layers
            ENDIF    !abs par
           ENDIF !Can_rad_mod=2 loop

          ENDDO
        ENDDO

      ENDDO

      RETURN
      END SUBROUTINE ALBPFT
#endif
