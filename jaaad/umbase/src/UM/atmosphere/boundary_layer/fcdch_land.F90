#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!   SUBROUTINES FCDCH_SEA AND FCDCH_LAND-----------------------------
!!!
!!!  Purpose: Calculate surface transfer coefficients at one or more
!!!           gridpoints.
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!  5.2   15/11/00   New Deck         M. Best
!    5.3  25/04/01  Add coastal tiling. Nic Gedney
!    5.5  12/02/03  Include code for mineral dust scheme. S Woodward
!  6.1  08/12/03  Add !CDIR NODEP to force vectorisation. R Barnes
!!!
!!!  Programming standard:
!!!
!!!  System component covered: Part of P243.
!!!
!!!  Documentation: UM Documentation Paper No 24, section P243.
!!!

!     SUBROUTINE FCDCH_SEA---------------------------------------------
!
!     Transfer coefficients for sea, sea-ice and leads
!
!     -----------------------------------------------------------------
!!  Arguments:---------------------------------------------------------


!!  Arguments:---------------------------------------------------------
      SUBROUTINE FCDCH_LAND(                                            &
     & ROW_LENGTH,ROWS,COR_MO_ITER,LAND_PTS,TILE_PTS,                   &
     & TILE_INDEX,LAND_INDEX,                                           &
     & DB,VSHR,Z0M,Z0H,ZH,Z1_UV,Z1_TQ,WIND_PROFILE_FACTOR,              &
     & CDV,CHV,CDV_STD,V_S,V_S_STD,RECIP_L_MO,U_S_STD,L_DUST,           &
     & LTIMER                                                           &
     &)
      IMPLICIT NONE

      INTEGER                                                           &
     & COR_MO_ITER          ! IN Switch for MO iteration correction

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                            ! IN Number of X points?
     &,ROWS                                                             &
                            ! IN Number of Y points?
     &,LAND_PTS                                                         &
                            ! IN Number of land points.
     &,TILE_PTS                                                         &
                            ! IN Number of tile points.
     &,TILE_INDEX(LAND_PTS)                                             &
                            ! IN Index of tile points.
     &,LAND_INDEX(LAND_PTS) ! IN Index of land points.

      LOGICAL                                                           &
     & LTIMER                                                           &
     &,L_DUST   !IN switch for mineral dust


      REAL                                                              &
     & DB(LAND_PTS)                                                     &
                     ! IN Buoyancy difference between surface and lowest
!                    !    temperature and humidity level in the
!                    !    atmosphere (m/s^2).
     &,VSHR(ROW_LENGTH,ROWS)                                            &
                             ! IN Wind speed difference between the surf
!                    !    the lowest wind level in the atmosphere (m/s).
     &,Z0M(LAND_PTS)                                                    &
                     ! IN Roughness length for momentum transport (m).
     &,Z0H(LAND_PTS)                                                    &
                     ! IN Roughness length for heat and moisture (m).
     &,ZH(ROW_LENGTH,ROWS)                                              &
                             ! IN Depth of boundary layer (m).
     &,Z1_UV(ROW_LENGTH,ROWS)                                           &
                             ! IN Height of lowest wind level (m).
     &,Z1_TQ(ROW_LENGTH,ROWS)                                           &
                             ! IN Height of lowest temperature and
!                    !    humidity level (m).
     &,WIND_PROFILE_FACTOR(LAND_PTS)
!                    ! IN for adjusting the surface transfer
!                    !    coefficients to remove form drag effects.

      REAL                                                              &
     & CDV(LAND_PTS)                                                    &
                     ! OUT Surface transfer coefficient for momentum
!                    !     including orographic form drag (m/s).
     &,CHV(LAND_PTS)                                                    &
                     ! OUT Surface transfer coefficient for
!                    !     heat, moisture & other scalars (m/s).
     &,CDV_STD(LAND_PTS)                                                &
!                    ! OUT Surface transfer coefficient for momentum
!                    !     excluding orographic form drag (m/s).
     &,V_S(LAND_PTS)                                                    &
                     ! OUT Surface layer scaling velocity
!                    !     including orographic form drag (m/s).
     &,V_S_STD(LAND_PTS)                                                &
!                    ! OUT Surface layer scaling velocity
!                    !     excluding orographic form drag (m/s).
     &,U_S_STD(LAND_PTS)                                                &
                     ! OUT Scaling velocity from middle of MO iteration
!                    !     - picked up in error by dust code!
     &,RECIP_L_MO(LAND_PTS)
!                    ! OUT Reciprocal of the Monin-Obukhov length
!                    !     (m^-1).


!*L  Workspace usage----------------------------------------------------
!
!     Local work arrays.
!
      REAL                                                              &
     & PHI_M(LAND_PTS)                                                  &
                       ! Monin-Obukhov stability function for momentum
!                      ! integrated to the model's lowest wind level.
     &,PHI_H(LAND_PTS) ! Monin-Obukhov stability function for scalars
!                      ! integrated to the model's lowest temperature
!                      ! and humidity level.
!
!*----------------------------------------------------------------------

      EXTERNAL PHI_M_H_LAND
      EXTERNAL TIMER

!*----------------------------------------------------------------------
!  Common and local constants.
#include "blopt8a.h"
#include "c_vkman.h"
      REAL BETA,THIRD
      PARAMETER (                                                       &
     & BETA=0.08,                                                       &
                    ! Tunable parameter in the surface layer scaling
!                   ! velocity formula (multiplying the turbulent
!                   ! convective scaling velocity).
     & THIRD=1./3.                                                      &
                    ! One third.
     &)
!
!  Define local variables
!
      INTEGER I,J,K,L ! Loop counter; horizontal field index.
      INTEGER IT    ! Iteration loop counter.
      INTEGER N_ITS ! Number of iterations for Monin-Obukhov length
!                   ! and stability functions.

      REAL                                                              &
     & B_FLUX                                                           &
                    ! Surface bouyancy flux over air density.
     &,U_S                                                              &
                    ! Iteration surface friction velocity
                    ! (effective value)
     &,U_S2                                                             &
                    ! Iteration surface friction velocity squared
                    ! (effective value)
     &,U_S_STD2                                                         &
                    ! Non-effective version of U_S2
     &,W_S          ! Surface turbulent convective scaling velocity.

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('FCDCH   ',3)
      ENDIF
!
!-----------------------------------------------------------------------
!! 1. Set initial values for the iteration.
!-----------------------------------------------------------------------
      IF (COR_MO_ITER == OFF) THEN
        N_ITS=5   ! original iteration count
!                 ! Found typically 2% from converged value
      ELSE 
        N_ITS=8   ! Found typically 0.2% from converged value
      END IF
!
!CDIR NODEP
      DO K=1,TILE_PTS
        L = TILE_INDEX(K)
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        IF (DB(L)  <   0.0 .AND. VSHR(I,J)  <   2.0) THEN
!-----------------------------------------------------------------------
!         Start the iteration from the convective limit.
!-----------------------------------------------------------------------
          RECIP_L_MO(L) = -VKMAN/(BETA*BETA*BETA*ZH(I,J))
        ELSE
!-----------------------------------------------------------------------
!         Start the iteration from neutral values.
!-----------------------------------------------------------------------
          RECIP_L_MO(L) = 0.0
        ENDIF
      ENDDO
! DEPENDS ON: phi_m_h_land
      CALL PHI_M_H_LAND (ROW_LENGTH,ROWS,LAND_PTS,TILE_PTS,             &
     &                   TILE_INDEX,LAND_INDEX,                         &
     &                   RECIP_L_MO,Z1_UV,Z1_TQ,Z0M,Z0H,                &
     &                   PHI_M,PHI_H,                                   &
     &                   LTIMER)
!
!CDIR NODEP
      DO K=1,TILE_PTS
        L = TILE_INDEX(K)
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        IF (DB(L)  <   0.0 .AND. VSHR(I,J)  <   2.0) THEN
!-----------------------------------------------------------------------
!         Start the iteration from the convective limit.
!-----------------------------------------------------------------------
          V_S_STD(L) = BETA *                                           &
     &        SQRT( BETA * ( VKMAN / PHI_H(L) ) * ZH(I,J) * (-DB(L)) )
          V_S(L) = V_S_STD(L)
        ELSE
!-----------------------------------------------------------------------
!         Start the iteration from neutral values.
!-----------------------------------------------------------------------
          V_S(L) = ( VKMAN / PHI_M(L) ) * VSHR(I,J)
          V_S_STD(L) = V_S(L) * WIND_PROFILE_FACTOR(L)
        ENDIF
        CHV(L) = ( VKMAN / PHI_H(L) ) * V_S_STD(L)
        CDV(L) = ( VKMAN / PHI_M(L) ) * V_S(L)
        CDV_STD(L) = CDV(L) * ( V_S_STD(L) / V_S(L) ) *                 &
     &                        WIND_PROFILE_FACTOR(L)
      ENDDO
!-----------------------------------------------------------------------
!! 2. Iterate to obtain sucessively better approximations for CD & CH.
!-----------------------------------------------------------------------
      DO IT = 1,N_ITS
!
       IF (COR_MO_ITER == OFF) THEN
!-----------------------------------------------------------------------
!        ! Original version with incorrect iteration of gustiness
!-----------------------------------------------------------------------
!CDIR NODEP
        DO K=1,TILE_PTS
          L = TILE_INDEX(K)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          B_FLUX = -CHV(L) * DB(L)
          U_S = SQRT( CDV(L) * VSHR(I,J) )
          U_S_STD(L) = SQRT( CDV_STD(L) * VSHR(I,J) )
          IF (DB(L)  <   0.0) THEN
            W_S = (ZH(I,J) * B_FLUX)**THIRD
            V_S(L) = SQRT(U_S*U_S + BETA*BETA*W_S*W_S)
              V_S_STD(L) =                                              &
     &             SQRT( U_S_STD(L)*U_S_STD(L) + BETA*BETA*W_S*W_S )
          ELSE
            V_S(L) = U_S
              V_S_STD(L) = U_S_STD(L)

          ENDIF
          RECIP_L_MO(L) = -VKMAN * B_FLUX /                             &
     &                     (V_S(L)*V_S(L)*V_S(L))
        ENDDO
       ELSE 
!-----------------------------------------------------------------------
!        ! Corrected version
!-----------------------------------------------------------------------
!CDIR NODEP
        DO K=1,TILE_PTS
          L = TILE_INDEX(K)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          B_FLUX = -CHV(L) * DB(L)
          U_S2   =  CDV(L) * VSHR(I,J)
          U_S_STD2 = CDV_STD(L) * VSHR(I,J)
          U_S_STD(L) = SQRT( U_S_STD2 ) 

          IF (DB(L)  <   0.0) THEN
            W_S = (ZH(I,J) * B_FLUX)**THIRD
!           ! Note that, during this iteration, CDV already includes
!           ! this gust enhancement and so U_S2 is not simply the 
!           ! friction velocity arising from the mean wind, hence:
            V_S(L) = SQRT( 0.5*( BETA*BETA*W_S*W_S                      &
     &                     + SQRT( (BETA*W_S)**4 + 4.0*U_S2*U_S2 )) )
            V_S_STD(L) = SQRT( 0.5*( BETA*BETA*W_S*W_S                  &
     &               + SQRT( (BETA*W_S)**4 + 4.0*U_S_STD2*U_S_STD2 )) )
          ELSE
            V_S(L) = SQRT( U_S2 )
            V_S_STD(L) = SQRT( U_S_STD2 )
          ENDIF

          RECIP_L_MO(L) = -VKMAN * B_FLUX /                             &
     &                     (V_S(L)*V_S(L)*V_S(L))
        ENDDO

       END IF  ! test on COR_MO_ITER

! DEPENDS ON: phi_m_h_land
        CALL PHI_M_H_LAND (ROW_LENGTH,ROWS,LAND_PTS,TILE_PTS,           &
     &                     TILE_INDEX,LAND_INDEX,                       &
     &                     RECIP_L_MO,Z1_UV,Z1_TQ,Z0M,Z0H,              &
     &                     PHI_M,PHI_H,                                 &
     &                     LTIMER)
!
!CDIR NODEP
        DO K=1,TILE_PTS
          L = TILE_INDEX(K)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          CHV(L) = ( VKMAN / PHI_H(L) ) * V_S_STD(L)
          CDV(L) = ( VKMAN / PHI_M(L) ) * V_S(L)
          CDV_STD(L) = CDV(L) * ( V_S_STD(L) / V_S(L) ) *               &
     &                          WIND_PROFILE_FACTOR(L)
        ENDDO
      ENDDO ! Iteration loop

!-----------------------------------------------------------------------
!! Set CD's and CH's to be dimensionless paremters
!-----------------------------------------------------------------------
!CDIR NODEP
      DO K=1,TILE_PTS
        L = TILE_INDEX(K)
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        CDV(L) = CDV(L) / VSHR(I,J)
        CDV_STD(L) = CDV_STD(L) / VSHR(I,J)
        CHV(L) = CHV(L) / VSHR(I,J)
      ENDDO


      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('FCDCH   ',4)
      ENDIF

      RETURN
      END SUBROUTINE FCDCH_LAND
#endif
