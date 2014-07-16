
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
      SUBROUTINE FCDCH_SEA(                                             &
     & ROW_LENGTH,ROWS,COR_MO_ITER,FLANDG,                              &
     & DB,VSHR,Z0M,Z0H,ZH,Z1_UV,Z1_TQ,                                  &
     & CDV,CHV,V_S,RECIP_L_MO,LTIMER                                    &
     &)
      IMPLICIT NONE

      INTEGER                                                           &
     & COR_MO_ITER          ! IN Switch for MO iteration correction

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                            ! IN Number of X points?
     &,ROWS                 ! IN Number of Y points?

      LOGICAL                                                           &
     & LTIMER

      REAL                                                              &
     & FLANDG(ROW_LENGTH,ROWS)                                          &
!                            ! IN Land fraction
     &,DB(ROW_LENGTH,ROWS)                                              &
                             ! IN Buoyancy difference between surface
!                            !    and lowest temperature and humidity
!                            !    level in the atmosphere (m/s^2).
     &,VSHR(ROW_LENGTH,ROWS)                                            &
                             ! IN Wind speed difference between the
!                            !    surface and the lowest wind level in
!                            ! the atmosphere (m/s).
     &,Z0M(ROW_LENGTH,ROWS)                                             &
                             ! IN Roughness length for momentum
!                            !    transport (m).
     &,Z0H(ROW_LENGTH,ROWS)                                             &
                             ! IN Roughness length for heat and
!                            !    moisture (m).
     &,ZH(ROW_LENGTH,ROWS)                                              &
                             ! IN Depth of boundary layer (m).
     &,Z1_UV(ROW_LENGTH,ROWS)                                           &
                             ! IN Height of lowest wind level (m).
     &,Z1_TQ(ROW_LENGTH,ROWS)! IN Height of lowest temperature and
!                            !    humidity level (m).

      REAL                                                              &
     & CDV(ROW_LENGTH,ROWS)                                             &
                             ! OUT Surface transfer coefficient for
!                    ! momentum including orographic form drag (m/s).
     &,CHV(ROW_LENGTH,ROWS)                                             &
                             ! OUT Surface transfer coefficient for
!                    !     heat, moisture & other scalars (m/s).
     &,V_S(ROW_LENGTH,ROWS)                                             &
                             ! OUT Surface layer scaling velocity
!                    !     including orographic form drag (m/s).
     &,RECIP_L_MO(ROW_LENGTH,ROWS)
!                    ! OUT Reciprocal of the Monin-Obukhov length
!                    !     (m^-1).

!*L  Workspace usage----------------------------------------------------
!
!     Local work arrays.
!
      REAL                                                              &
     & PHI_M(ROW_LENGTH,ROWS)                                           &
                              ! Monin-Obukhov stability function for
!                             ! momentum integrated to the model's
!                             ! lowest wind level.
     &,PHI_H(ROW_LENGTH,ROWS) ! Monin-Obukhov stability function for
!                             ! scalars integrated to the model's lowest
!                             ! temperature and humidity level.
!
!*----------------------------------------------------------------------

      EXTERNAL PHI_M_H_SEA
      EXTERNAL TIMER

!*----------------------------------------------------------------------
!  Common and local constants.
! Start blopt8a

! Description:
!   Permissible settings for BL options.
!
! Current Code Owner: J. M. Edwards
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 6.2      27/01/06 Original code.  J. M. Edwards
!
      INTEGER, PARAMETER :: Off = 0  ! Switch disabled
      INTEGER, PARAMETER :: On  = 1  ! Switch enabled
!
      INTEGER, PARAMETER :: BrownGrant97 = 1
      INTEGER, PARAMETER :: BrownGrant97_limited = 2
!     Options for non-gradient stress following
!     Brown and Grant (1997), version 2 including a limit on its size
!
!     Options for flux gradient formulation
      INTEGER, PARAMETER :: Locketal2000   = 0
!       Flux gradients as in Lock et al. (2000)
      INTEGER, PARAMETER :: HoltBov1993 = 1
!       Flux gradients as in Lock et al (2000) but using
!       coefficients from Holtslag and Boville (1993)
      INTEGER, PARAMETER :: LockWhelan2006 = 2
!       Flux gradients as in Lock and Whelan (2006)
!
!     Options for form drag
      INTEGER, PARAMETER :: No_drag         = 0
      INTEGER, PARAMETER :: Effective_z0    = 1
      INTEGER, PARAMETER :: Explicit_stress = 2
!
!     Options for marine boundary layers
      INTEGER, PARAMETER :: Fixed_Z0T = 0
!       Stanard flixed value of thermal roughness length over sea
      INTEGER, PARAMETER :: SurfDivZ0T = 1
!       Thermal roughness length over sea defined from surface
!       divergence theory
      INTEGER, PARAMETER :: DynDiag_ZL = 1
!       The ratio of the height of the inversion to the surface
!       Obukhov length is used a dynamic criterion in the
!       diagnosis of BL types
!
      INTEGER, PARAMETER :: Use_Correct_Ustar = 2
!       Option under the COR_MO_ITER switch for the dust scheme
!       to use the correct ustar
!
!     Options for stable boundary layers
      INTEGER, PARAMETER ::  Long_tails           = 0
      INTEGER, PARAMETER ::  Sharpest             = 1
      INTEGER, PARAMETER ::  Sharp_sea_long_land  = 2
      INTEGER, PARAMETER ::  Mes_tails            = 3
      INTEGER, PARAMETER ::  Louis_tails          = 4
      INTEGER, PARAMETER ::  Depth_based          = 5
      INTEGER, PARAMETER ::  Sharp_sea_mes_land   = 6
      INTEGER, PARAMETER ::  LEM_stability        = 7
      INTEGER, PARAMETER ::  Sharp_sea_Louis_land = 8

!     Options for Prandtl number (in local Ri scheme)
      INTEGER, PARAMETER ::  Constant_SBL = 0
      INTEGER, PARAMETER ::  LockMailhot2004 = 1

! End blopt8a
! C_VKMAN start
      REAL,PARAMETER:: VKMAN=0.4 ! Von Karman's constant
! C_VKMAN end
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
      INTEGER I,J   ! Loop counters; horizontal field index.
      INTEGER IT    ! Iteration loop counter.
      INTEGER N_ITS ! Number of iterations for Monin-Obukhov length
!                   ! and stability functions.

      REAL                                                              &
     & B_FLUX                                                           &
                    ! Surface bouyancy flux over air density.
     &,U_S                                                              &
                    ! Iteration surface friction velocity
     &,U_S2                                                             &
                    ! Iteration surface friction velocity squared
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

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        IF ( FLANDG(I,J) <  1.0 ) THEN
          IF (DB(I,J)  <   0.0 .AND. VSHR(I,J)  <   2.0) THEN
!-----------------------------------------------------------------------
!         Start the iteration from the convective limit.
!-----------------------------------------------------------------------
            RECIP_L_MO(I,J) = -VKMAN/(BETA*BETA*BETA*ZH(I,J))
          ELSE
!-----------------------------------------------------------------------
!         Start the iteration from neutral values.
!-----------------------------------------------------------------------
            RECIP_L_MO(I,J) = 0.0
          ENDIF
        ENDIF  ! SEA_MASK
       ENDDO
      ENDDO

! DEPENDS ON: phi_m_h_sea
      CALL PHI_M_H_SEA (ROW_LENGTH,ROWS,FLANDG,                         &
     &                  RECIP_L_MO,Z1_UV,Z1_TQ,Z0M,Z0H,                 &
     &                  PHI_M,PHI_H,                                    &
     &                  LTIMER)
!
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        IF ( FLANDG(I,J) <  1.0 ) THEN
          IF (DB(I,J)  <   0.0 .AND. VSHR(I,J)  <   2.0) THEN
!-----------------------------------------------------------------------
!         Start the iteration from the convective limit.
!-----------------------------------------------------------------------
            V_S(I,J) = BETA *                                           &
     &          SQRT( BETA * ( VKMAN / PHI_H(I,J) ) *                   &
     &                 ZH(I,J) * (-DB(I,J)) )
          ELSE
!-----------------------------------------------------------------------
!         Start the iteration from neutral values.
!-----------------------------------------------------------------------
            V_S(I,J) = ( VKMAN / PHI_M(I,J) ) * VSHR(I,J)
          ENDIF
          CHV(I,J) = ( VKMAN / PHI_H(I,J) ) * V_S(I,J)
          CDV(I,J) = ( VKMAN / PHI_M(I,J) ) * V_S(I,J)
        ENDIF  ! SEA_MASK
       ENDDO
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
         DO J=1,ROWS
         DO I=1,ROW_LENGTH
          IF ( FLANDG(I,J) <  1.0 ) THEN
            B_FLUX = -CHV(I,J) * DB(I,J)
            U_S = SQRT( CDV(I,J) * VSHR(I,J) )
            IF (DB(I,J)  <   0.0) THEN
              W_S = (ZH(I,J) * B_FLUX)**THIRD
              V_S(I,J) = SQRT(U_S*U_S + BETA*BETA*W_S*W_S)
            ELSE
              V_S(I,J) = U_S
            ENDIF
            RECIP_L_MO(I,J) = -VKMAN * B_FLUX /                         &
     &                       (V_S(I,J)*V_S(I,J)*V_S(I,J))
          ENDIF  ! SEA_MASK
         ENDDO
         ENDDO
       ELSE 
!-----------------------------------------------------------------------
!        ! Corrected version
!-----------------------------------------------------------------------
         DO J=1,ROWS
         DO I=1,ROW_LENGTH
          IF ( FLANDG(I,J) <  1.0 ) THEN
            B_FLUX = -CHV(I,J) * DB(I,J)
            U_S2   =  CDV(I,J) * VSHR(I,J)
            IF (DB(I,J)  <   0.0) THEN
              W_S = (ZH(I,J) * B_FLUX)**THIRD
!             ! Note that, during this iteration, CDV already includes
!             ! this gust enhancement and so U_S2 is not simply the 
!             ! friction velocity arising from the mean wind, hence:
              V_S(I,J) = SQRT( 0.5*( BETA*BETA*W_S*W_S                  &
     &                  + SQRT( (BETA*W_S)**4 + 4.0*U_S2*U_S2 )))
            ELSE
              V_S(I,J) = SQRT( U_S2 )
            ENDIF
            RECIP_L_MO(I,J) = -VKMAN * B_FLUX /                         &
     &                       (V_S(I,J)*V_S(I,J)*V_S(I,J))
          ENDIF  ! SEA_MASK
         ENDDO
         ENDDO

       END IF  ! test on COR_MO_ITER

! DEPENDS ON: phi_m_h_sea
        CALL PHI_M_H_SEA (ROW_LENGTH,ROWS,FLANDG,                       &
     &                    RECIP_L_MO,Z1_UV,Z1_TQ,Z0M,Z0H,               &
     &                    PHI_M,PHI_H,                                  &
     &                    LTIMER)
!

        DO J=1,ROWS
         DO I=1,ROW_LENGTH
           IF ( FLANDG(I,J) <  1.0 ) THEN
            CHV(I,J) = ( VKMAN / PHI_H(I,J) ) * V_S(I,J)
            CDV(I,J) = ( VKMAN / PHI_M(I,J) ) * V_S(I,J)
          ENDIF  ! SEA_MASK
        ENDDO
       ENDDO
      ENDDO ! Iteration loop

!-----------------------------------------------------------------------
!! Set CD's and CH's to be dimensionless paremters
!-----------------------------------------------------------------------
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        IF ( FLANDG(I,J) <  1.0 ) THEN
          CDV(I,J) = CDV(I,J) / VSHR(I,J)
          CHV(I,J) = CHV(I,J) / VSHR(I,J)
        ENDIF  ! SEA_MASK
       ENDDO
      ENDDO


      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('FCDCH   ',4)
      ENDIF

      RETURN
      END SUBROUTINE FCDCH_SEA


!!  Arguments:---------------------------------------------------------
