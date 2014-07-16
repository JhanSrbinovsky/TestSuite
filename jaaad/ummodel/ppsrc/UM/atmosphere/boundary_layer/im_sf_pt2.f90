
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE IM_SF_PT2 ---------------------------------------------
!!!
!!!  Purpose: Calculate implicit increments to surface variables
!!!           for the unconditionally stable and non-oscillatory 
!!!           BL numerical solver. 
!!!
!!!  Model           Modification history
!!! version  Date
!!!  6.4   10/01/07   New Deck         M. Diamantakis
!!!
!!!  Programming standard: UMDP4
!!!
!!!
!!!  Documentation: 
!!!          http://www-nwp/~frmd/DR/Reports/new_BLsolver_guide.ps
!!!
!!!---------------------------------------------------------------------
!!  Arguments :-
      SUBROUTINE IM_SF_PT2 (                                            &
     & OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,LAND_PTS                      &
     &,LAND_INDEX,NTILES,TILE_INDEX,TILE_PTS                            &
     &,FLANDG,TILE_FRAC,SNOW_TILE,ICE_FRACT                             &
     &,GAMMA_IN,GAMMA1_IN,GAMMA2_IN,ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE &
     &,RESFT,RHOKPM,RHOKPM_POT,RHOKPM_SICE                              &
     &,RHOKM_U_1,RHOKM_V_1,RHOKH_1,RHOKH1_SICE                          &
     &,CT_CTQ_1,CTCTQ1,DQW_1,DTL_1,DQW1_1,DTL1_1,CQ_CM_U_1              &
     &,CQ_CM_V_1,DU_1,DV_1,DU_STAR1,DV_STAR1,FLANDG_U,FLANDG_V          &
     &,FQW_GB,FTL_GB                                                    &
     &,TAUX_1,TAUX_LAND,TAUX_LAND_star,TAUX_SSI,TAUX_SSI_star,TAUY_1    &
     &,TAUY_LAND,TAUY_LAND_star,TAUY_SSI,TAUY_SSI_star                  & 
     &,FQW_TILE,EPOT_TILE,FTL_TILE,FQW_ICE,FTL_ICE,E_SEA,H_SEA          &
     &,L_CORRECT,L_FLUX_BC,LTIMER                                       &
     &)


      IMPLICIT NONE

      LOGICAL                                                           &
     & LTIMER                                                           &
     &,L_CORRECT                                                        &
                                   ! Flag for BL solver
     &,L_FLUX_BC                   ! Logical for prescribed
                                   ! surface fluxes (SCM)

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                                   ! IN Number of X points?
     &,ROWS                                                             &
                                   ! IN Number of Y points?
     &,N_ROWS                                                           &
                                   ! Local number of rows in a v field
     &,OFF_X                                                            &
                                   ! Size of small halo in i
     &,OFF_Y                                                            &
                                   ! Size of small halo in j.
     &,LAND_PTS                                                         &
                                   ! IN Total number of land points.
     &,LAND_INDEX(LAND_PTS)                                             &
                                   ! IN Index of land points.
     &,NTILES                                                           &
                                   ! IN Number of land surface tiles.
     &,TILE_INDEX(LAND_PTS,NTILES)                                      &
                                   ! IN Index of tile points.
     &,TILE_PTS(NTILES)            ! IN Number of tiles.


      REAL                                                              &
     & FLANDG(ROW_LENGTH,ROWS)                                          &
                                   ! IN Land fraction
     &,TILE_FRAC(LAND_PTS,NTILES)                                       &
                                   ! IN Tile fraction
     &,SNOW_TILE(LAND_PTS,NTILES)                                       &
                                   ! IN Lying snow on land tiles (kg/m2)
     &,ICE_FRACT(ROW_LENGTH,ROWS)                                       &
                                   ! IN Fraction of grid-box which is
!                                  !    sea-ice (decimal fraction).
     &,GAMMA_IN                                                         &
                                   ! IN Implicit weighting.
     &,GAMMA1_IN(ROW_LENGTH,ROWS)                                       &
     &,GAMMA2_IN(ROW_LENGTH,ROWS)                                       &
                                   ! IN Implicit weights for uncond. 
                                   !    stable non-oscillatory BL solver
     &,ALPHA1(LAND_PTS,NTILES)                                          &
                                   ! IN Gradient of saturated specific
!                                  !    humidity with respect to
!                                  !    temperature between the bottom
!                                  !    model layer and the surface.
     &,ALPHA1_SICE(ROW_LENGTH,ROWS)                                     &
                                   ! IN ALPHA1 for sea-ice
     &,ASHTF(ROW_LENGTH,ROWS)                                           &
                                   ! IN Coefficient to calculate surface
!                                  !    heat flux into soil or sea-ice
!                                  !    (W/m2/K).

     &,ASHTF_TILE(LAND_PTS,NTILES)                                      &
                                   ! IN Coefficient to calculate heat
!                                  !    flux into land tiles (W/m2/K).
     &,RESFT(LAND_PTS,NTILES)                                           &
                                   ! IN Total resistance factor
     &,RHOKPM(LAND_PTS,NTILES)                                          &
                                   ! IN Surface exchange coeff for tiles
     &,RHOKPM_POT(LAND_PTS,NTILES)                                      &
                                   ! IN Surface exchange coeff for
!                                       potential evaporation over tiles
!
     &,EPOT_TILE(land_pts,ntiles)                                       &
                                   ! IN surface tile potential
!                                  !    evaporation
     &,E_EPOT_TILE(land_pts,ntiles)                                     &
                                   ! Work ratio of explicit
!                                  !      EPOT/E
!                                  !    evaporation
!
     &,RHOKPM_SICE(ROW_LENGTH,ROWS)                                     &
                                   ! IN Sea-ice surface exchange coeff.
     &,RHOKM_U_1(ROW_LENGTH,ROWS)                                       &
                                   ! IN Level 1 exchange coefficient for
!                                  !    momentum
     &,RHOKM_V_1(ROW_LENGTH,N_ROWS)                                     &
                                   ! IN Level 1 exchange coefficient for
!                                  !    momentum
     &,RHOKH_1(LAND_PTS,NTILES)                                         &
                                   ! IN Surface exchange coeffs for FTL

     &,RHOKH1_SICE(ROW_LENGTH,ROWS)                                     &
                                   ! IN Sea and sea-ice surface exchange
     &,CT_CTQ_1(ROW_LENGTH,ROWS)                                        &
                                   ! IN Coefficient in T and q
!                                  !     tri-diagonal implicit matrix
     &,CTCTQ1(ROW_LENGTH,ROWS)                                          &
     &,CQ_CM_U_1(ROW_LENGTH,ROWS)                                       &
                                   ! IN Coefficient in U and V
!                                  !     tri-diagonal implicit matrix
     &,CQ_CM_V_1(ROW_LENGTH,N_ROWS)                                     &
                                   ! IN Coefficient in U and V
!                                  !     tri-diagonal implicit matrix
     &,DQW_1(ROW_LENGTH,ROWS)                                           &
                                   ! IN Level 1 increment to q field
     &,DTL_1(ROW_LENGTH,ROWS)                                           &
                                   ! IN Level 1 increment to T field
     &,DQW1_1(ROW_LENGTH,ROWS)                                          &
                                   ! IN Incr obtained by semi-implicit inte
     &,DTL1_1(ROW_LENGTH,ROWS)                                          &
                                   ! IN Incr obtained by semi-implicit inte
     &,DU_1(1-OFF_X:ROW_LENGTH+OFF_X,                                   &
     &      1-OFF_Y:ROWS+OFF_Y)                                         &
                                   ! IN Level 1 increment to u wind
!                                  !    field
     &,DV_1(1-OFF_X:ROW_LENGTH+OFF_X,                                   &
     &      1-OFF_Y:N_ROWS+OFF_Y)                                       &
                                   ! IN Level 1 increment to v wind
!                                  !    field
     &,DU_STAR1(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &,DV_STAR1(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y)         &
     &,FLANDG_U(ROW_LENGTH,ROWS)                                        &
                                   ! IN Land fraction on U grid.
     &,FLANDG_V(ROW_LENGTH,N_ROWS) ! IN Land fraction on V grid.


      REAL                                                              &
     & FQW_GB(ROW_LENGTH,ROWS)                                          &
                                   ! INOUT Grid-box value of QW flux at
!                                  !       Kg/sq m/s
     &,FTL_GB(ROW_LENGTH,ROWS)                                          &
                                   ! INOUT Grid-box value of TL flux at
!                                  !       i.e. H/Cp where H is sensible
!                                  !       in W per sq m).
     &,TAUX_1(ROW_LENGTH,ROWS)                                          &
                                   ! OUT   x-component of turbulent
!                                  !       stress at surface.
     &,TAUX_LAND(ROW_LENGTH,ROWS)                                       &
                                   ! INOUT x-component of turbulent
!                                  !       stress at land surface.
     &,TAUX_SSI(ROW_LENGTH,ROWS)                                        &
                                   ! INOUT x-component of turbulent
!                                  !       stress at sea surface.
     &,TAUY_1(ROW_LENGTH,N_ROWS)                                        &
                                   ! OUT   y-component of turbulent
!                                  !       stress at surface.
     &,TAUY_LAND(ROW_LENGTH,N_ROWS)                                     &
                                   ! INOUT y-component of turbulent
!                                  !       stress at land surface.
     &,TAUY_SSI(ROW_LENGTH,N_ROWS)                                      &
                                   ! INOUT y-component of turbulent
!                                  !       stress at sea surface.
     &,TAUX_LAND_star(ROW_LENGTH,ROWS)                                  &
     &,TAUX_SSI_star(ROW_LENGTH,ROWS)                                   &
     &,TAUY_LAND_star(ROW_LENGTH,N_ROWS)                                &
     &,TAUY_SSI_star(ROW_LENGTH,N_ROWS)                                 &
                                   ! INOUT as above, temporarily needed
                                   !       by predictor stage of the 
                                   !       new uncond stable BL solver
     &,FQW_TILE(LAND_PTS,NTILES)                                        &
                                   ! INOUT Tile flux of QW. Kg/sq m/s
     &,FTL_TILE(LAND_PTS,NTILES)                                        &
                                   ! INOUT Tile flux of TL
     &,E_SEA(ROW_LENGTH,ROWS)                                           &
                                   ! INOUT Evaporation from sea times
!                                  !       leads fraction (kg/m2/s).
!                                  !       Zero over land.
     &,H_SEA(ROW_LENGTH,ROWS)      ! INOUT Surface sensible heat flux ov
!                                  !       sea times leads fraction (W/m
!                                  !       Zero over land.


!  External references :-
      EXTERNAL TIMER


!  Local and other symbolic constants :-
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------


      REAL LS
      PARAMETER (                                                       &
     & LS=LC+LF                                                         &
                    ! Latent heat of sublimation (J per kg).
     &)

! Workspace :-
      REAL                                                              &
     & FQW_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! "Explicit" surface flux of QW for
!                                  !  sea-ice fraction of gridsquare.
     &,FTL_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! "Explicit" surface flux of TL for
!                                  !  sea-ice fraction of gridsquare.
     &,LAT_HT                                                           &
                ! Latent heat of evaporation for snow-free land
!               ! or sublimation for snow-covered land and ice.
     &,APART(ROW_LENGTH,ROWS,2)                                         &
                                   ! Temporary array
     &,BPART(ROW_LENGTH,ROWS,2)                                         &
                                   ! Temporary array
     &,RECIP(ROW_LENGTH,ROWS)                                           &
                                   ! Temporary array
     &,FTL_LAND(ROW_LENGTH,ROWS)                                        &
                                   ! Temporary array
     &,FQW_LAND(ROW_LENGTH,ROWS)
                                   ! Temporary array

!  Local scalars :-
      INTEGER                                                           &
     & I,J                                                              &
                ! Loop counter (horizontal field index).
     &,K                                                                &
                ! Loop counter (tile index).
     &,L                                                                &
                ! Loop counter (horizontal land index).
     &,N        ! Loop counter (tile counter).

      REAL                                                              &
     & FTL_OLD                                                          &
                ! Used to hold current value of FTL_GB before updating
     &,GAMMA_1                                                          &
     &,GAMMA1                                                           &
     &,GAMMA2


      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('IM_SF_PT2 ',3)
      ENDIF
!
!-------------------------------------------------------------------------
!
! Now compute surface stresses for the 2nd stage (predictor) of 
! the new scheme (BL solver) using its discretization.
!
!-------------------------------------------------------------------------
!
      IF ( .NOT. L_CORRECT ) THEN 
!-------------------------------------------------------------------------
!
! Compute surface stresses for the 1st stage (predictor) of 
! the new scheme (BL solver) using its discretization.
!
!-------------------------------------------------------------------------
        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            GAMMA1 = GAMMA1_IN(I,J)
            GAMMA2 = GAMMA2_IN(I,J)
            IF ( FLANDG_U(I,J) > 0.0 ) THEN
              TAUX_LAND_star(I,J) = ( GAMMA2*TAUX_LAND(I,J) +           &
     &                 GAMMA1*RHOKM_U_1(I,J)*DU_1(I,J) ) /              &
     &                ( 1.0 + GAMMA1*RHOKM_U_1(I,J)*CQ_CM_U_1(I,J) )
            ELSE
              TAUX_LAND_star(I,J) = 0.0
            ENDIF
            IF ( FLANDG_U(I,J) < 1.0 ) THEN
              TAUX_SSI_star(I,J) = ( GAMMA2*TAUX_SSI(I,J) +             &
     &                 GAMMA1*RHOKM_U_1(I,J)*DU_1(I,J) ) /              &
     &                ( 1.0 + GAMMA1*RHOKM_U_1(I,J)*CQ_CM_U_1(I,J) )
            ELSE
              TAUX_SSI_star(I,J) = 0.0
            ENDIF
            TAUX_1(I,J) = FLANDG_U(I,J)*TAUX_LAND_star(I,J)             &
     &                  + ( 1.0-FLANDG_U(I,J))*TAUX_SSI_star(I,J)
          ENDDO
        ENDDO

        DO J=1,N_ROWS
          DO I=1,ROW_LENGTH
            GAMMA1 = GAMMA1_IN(I,J)
            GAMMA2 = GAMMA2_IN(I,J)
            IF ( FLANDG_V(I,J) > 0.0 ) THEN
              TAUY_LAND_star(I,J) = ( GAMMA2*TAUY_LAND(I,J) +           &
     &                 GAMMA1*RHOKM_V_1(I,J)*DV_1(I,J) ) /              &
     &                ( 1.0 + GAMMA1*RHOKM_V_1(I,J)*CQ_CM_V_1(I,J) )
            ELSE
              TAUY_LAND_star(I,J) = 0.0
            ENDIF
            IF ( FLANDG_V(I,J) < 1.0 ) THEN
              TAUY_SSI_star(I,J) = ( GAMMA2*TAUY_SSI(I,J) +             &
     &                 GAMMA1*RHOKM_V_1(I,J)*DV_1(I,J) ) /              &
     &                ( 1.0 + GAMMA1*RHOKM_V_1(I,J)*CQ_CM_V_1(I,J) )
            ELSE
              TAUY_SSI_star(I,J) = 0.0
            ENDIF
            TAUY_1(I,J) = FLANDG_V(I,J)*TAUY_LAND_star(I,J)             &
     &                + ( 1.0-FLANDG_V(I,J))*TAUY_SSI_star(I,J)
          ENDDO
        ENDDO

        GAMMA_1 = GAMMA_IN
!
! time weights for specified scalar fluxes
!
        IF ( L_FLUX_BC ) GAMMA_1 = 0.0 
!
!-------------------------------------------------------------------------
!
! Compute scalar surface fluxes using the standard scheme described
! in MOSES 2.2 technical documentation (hadley centre tecnical note 30).
! This needs to be done only in the first stage of the scheme (predictor).
! The same fluxes will be used for the 2nd stage (corrector).
! NOTE: scalar surface fluxes could be computed using the new scheme
!       but with a penalty in code complexity.
! 
!------------------------------------------------------------------------
!
! Initialise APART and BPART to zero
        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            APART(I,J,1)=0.0
            APART(I,J,2)=0.0
            BPART(I,J,1)=0.0
            BPART(I,J,2)=0.0
            FTL_LAND(I,J)=0.0
            FQW_LAND(I,J)=0.0
          ENDDO
        ENDDO

! Land tiles
        DO N=1,NTILES
          DO K=1,TILE_PTS(N)
            L = TILE_INDEX(K,N)
            J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
            I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
            LAT_HT = LC
            IF (SNOW_TILE(L,N) >  0.) LAT_HT = LS

            APART(I,J,1)=APART(I,J,1) - TILE_FRAC(L,N) *                &
     &               GAMMA_1 * RHOKPM(L,N) *                            &
     &            ( LAT_HT*RESFT(L,N)*RHOKH_1(L,N)*ALPHA1(L,N) +        &
     &                         ASHTF_TILE(L,N) )
            APART(I,J,2)=APART(I,J,2) + TILE_FRAC(L,N) *                &
     &               GAMMA_1 * RHOKPM(L,N) *                            &
     &               LAT_HT*RESFT(L,N)*RHOKH_1(L,N)
            BPART(I,J,1)=BPART(I,J,1) + TILE_FRAC(L,N) *                &
     &               GAMMA_1 * RESFT(L,N)*RHOKPM(L,N) *                 &
     &               CP*RHOKH_1(L,N)*ALPHA1(L,N)
            BPART(I,J,2)=BPART(I,J,2) - TILE_FRAC(L,N) *                &
     &               GAMMA_1 * RESFT(L,N)*RHOKPM(L,N) *                 &
     &               ( CP*RHOKH_1(L,N) + ASHTF_TILE(L,N) )
          ENDDO
        ENDDO

! Sea points
        DO J=1,ROWS
          DO I=1,ROW_LENGTH

            IF ( FLANDG(I,J) < 1.0 .AND. ICE_FRACT(I,J) >  0.0 ) THEN
! Sea ice point
              APART(I,J,1)=FLANDG(I,J)*APART(I,J,1)                     &
     &       - GAMMA_1 * (1.0-FLANDG(I,J)) * ICE_FRACT(I,J)             &
     &      * RHOKPM_SICE(I,J) *                                        &
     &      ( LS*RHOKH1_SICE(I,J)*ALPHA1_SICE(I,J) + ASHTF(I,J) )       &
     &       - GAMMA_1 * (1.0-FLANDG(I,J)) * ( 1.0 - ICE_FRACT(I,J) )   &
     &      * RHOKH1_SICE(I,J)

              APART(I,J,2)=FLANDG(I,J)*APART(I,J,2)                     &
     &       + GAMMA_1 * (1.0-FLANDG(I,J)) *ICE_FRACT(I,J)              &
     &       * RHOKPM_SICE(I,J) * LS*RHOKH1_SICE(I,J)

              BPART(I,J,1)=FLANDG(I,J)*BPART(I,J,1)                     &
     &       + GAMMA_1 * ICE_FRACT(I,J) * ( 1.0 - FLANDG(I,J) )         &
     &       * RHOKPM_SICE(I,J) *CP*RHOKH1_SICE(I,J)*ALPHA1_SICE(I,J)

              BPART(I,J,2)=FLANDG(I,J)*BPART(I,J,2)                     &
     &       - GAMMA_1 * ICE_FRACT(I,J) * ( 1.0 - FLANDG(I,J) )         &
     &       * RHOKPM_SICE(I,J) * ( CP*RHOKH1_SICE(I,J) + ASHTF(I,J) )  &
     &       - GAMMA_1 * ( 1.0 - ICE_FRACT(I,J) )                       &
     &       * ( 1.0 - FLANDG(I,J) ) * RHOKH1_SICE(I,J)

            ELSEIF(FLANDG(I,J) <  1.0 .AND.                             &
     &         .NOT.ICE_FRACT(I,J) >  0.0) THEN
! Ordinary sea point
              APART(I,J,1)= FLANDG(I,J)*APART(I,J,1)                    &
     &       - GAMMA_1 * ( 1.0 - FLANDG(I,J) ) * RHOKH1_SICE(I,J)
              APART(I,J,2)= FLANDG(I,J)*APART(I,J,2)
              BPART(I,J,1)= FLANDG(I,J)*BPART(I,J,1)
              BPART(I,J,2)= FLANDG(I,J)*BPART(I,J,2)                    &
     &       - GAMMA_1*( 1.0 - FLANDG(I,J) )*RHOKH1_SICE(I,J)

            ENDIF
          ENDDO
        ENDDO

! Land tiles
      DO N=1,NTILES
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          E_EPOT_TILE(L,N)=1.0
          IF(EPOT_TILE(L,N) >  0.0.AND.FQW_TILE(L,N) >  0.0)            &
     &       E_EPOT_TILE(L,N)=EPOT_TILE(L,N)/FQW_TILE(L,N)
        ENDDO
      ENDDO

! Calculate grid-box fluxes of heat and moisture
        DO J=1,ROWS
          DO I=1,ROW_LENGTH

          RECIP(I,J)=( 1.0 + CTCTQ1(I,J)*APART(I,J,1) ) *               &
     &           ( 1.0 + CTCTQ1(I,J)*BPART(I,J,2) ) -                   &
     &           CTCTQ1(I,J)*APART(I,J,2)*CTCTQ1(I,J)*BPART(I,J,1)

          FTL_OLD=FTL_GB(I,J)

          FTL_GB(I,J) = ( ( 1.0 + CTCTQ1(I,J)*BPART(I,J,2) ) *          &
     &               ( FTL_OLD + APART(I,J,1)*DTL1_1(I,J) +             &
     &                 APART(I,J,2)*DQW1_1(I,J)) -                      &
     &                 CTCTQ1(I,J)*APART(I,J,2) * ( FQW_GB(I,J) +       &
     &                 BPART(I,J,1)*DTL1_1(I,J) +                       &
     &                 BPART(I,J,2)*DQW1_1(I,J)) ) / RECIP(I,J)

          FQW_GB(I,J) = ( ( 1.0 + CTCTQ1(I,J)*APART(I,J,1) ) *          &
     &                ( FQW_GB(I,J) + BPART(I,J,1)*DTL1_1(I,J) +        &
     &                  BPART(I,J,2)*DQW1_1(I,J)) -                     &
     &                  CTCTQ1(I,J)*BPART(I,J,1) * ( FTL_OLD +          &
     &                  APART(I,J,1)*DTL1_1(I,J) +                      &
     &                  APART(I,J,2)*DQW1_1(I,J)) ) / RECIP(I,J)

          ENDDO
        ENDDO

! Make implicit correction to tile fluxes

! Land tiles
        DO N=1,NTILES
          DO K=1,TILE_PTS(N)
            L = TILE_INDEX(K,N)
            J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
            I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
            LAT_HT = LC
            IF (SNOW_TILE(L,N) >  0.) LAT_HT = LS

            FTL_TILE(L,N)=FTL_TILE(L,N) -                               &
     &               GAMMA_1 * RHOKPM(L,N) *                            &
     &            ( LAT_HT*RESFT(L,N)*RHOKH_1(L,N)*ALPHA1(L,N) +        &
     &                       ASHTF_TILE(L,N) ) *                        &
     &         ( DTL1_1(I,J) - CTCTQ1(I,J)*FTL_GB(I,J) ) +              &
     &               GAMMA_1 * RHOKPM(L,N) *                            &
     &               LAT_HT*RESFT(L,N)*RHOKH_1(L,N) *                   &
     &         ( DQW1_1(I,J) - CTCTQ1(I,J)*FQW_GB(I,J) )

            FQW_TILE(L,N)=FQW_TILE(L,N) +                               &
     &               GAMMA_1 * RESFT(L,N)*RHOKPM(L,N) *                 &
     &               CP*RHOKH_1(L,N)*ALPHA1(L,N) *                      &
     &         ( DTL1_1(I,J) - CTCTQ1(I,J)*FTL_GB(I,J) ) -              &
     &               GAMMA_1 * RESFT(L,N)*RHOKPM(L,N) *                 &
     &               ( CP*RHOKH_1(L,N) + ASHTF_TILE(L,N) ) *            &
     &         ( DQW1_1(I,J) - CTCTQ1(I,J)*FQW_GB(I,J) )

            EPOT_TILE(L,N)=FQW_TILE(L,N)*E_EPOT_TILE(L,N)

            FQW_LAND(I,J)=FQW_LAND(I,J)+FQW_TILE(L,N)*TILE_FRAC(L,N)
            FTL_LAND(I,J)=FTL_LAND(I,J)+FTL_TILE(L,N)*TILE_FRAC(L,N)

          ENDDO
        ENDDO

! Sea points
        DO J=1,ROWS
          DO I=1,ROW_LENGTH

            IF(FLANDG(I,J) <  1.0 .AND. ICE_FRACT(I,J) >  0.0) THEN
! Sea ice point
              H_SEA(I,J)=H_SEA(I,J) - GAMMA_1*(1.0-ICE_FRACT(I,J)) *    &
     &                   CP * RHOKH1_SICE(I,J) *                        &
     &                     ( DTL1_1(I,J) - CTCTQ1(I,J) * FTL_GB(I,J) )
              E_SEA(I,J)=E_SEA(I,J) - GAMMA_1*(1.0-ICE_FRACT(I,J)) *    &
     &                   RHOKH1_SICE(I,J) *                             &
     &                     ( DQW1_1(I,J) - CTCTQ1(I,J) * FQW_GB(I,J) )
             FTL_ICE(I,J)=(FTL_GB(I,J)                                  &
     &         - FTL_LAND(I,J)*FLANDG(I,J))/(1.-FLANDG(I,J))            &
     &         - H_SEA(I,J)/CP
             FQW_ICE(I,J)=(FQW_GB(I,J)                                  &
     &        - FQW_LAND(I,J)*FLANDG(I,J))/(1.-FLANDG(I,J))             &
     &        - E_SEA(I,J)

            ELSEIF(FLANDG(I,J) <  1.0 .AND.                             &
     &             .NOT.ICE_FRACT(I,J) >  0.0) THEN
! Ordinary sea point
              H_SEA(I,J)=CP * (FTL_GB(I,J)                              &
     &             - FTL_LAND(I,J)*FLANDG(I,J))/(1.-FLANDG(I,J))
              E_SEA(I,J)=(FQW_GB(I,J)                                   &
     &             - FQW_LAND(I,J)*FLANDG(I,J))/(1.-FLANDG(I,J))
              FTL_ICE(I,J)=0.0
              FQW_ICE(I,J)=0.0

            ENDIF
          ENDDO
        ENDDO

      ELSE ! IF L_CORRECT==TRUE THEN: 

        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            GAMMA1 = GAMMA1_IN(I,J)
            GAMMA2 = GAMMA2_IN(I,J)
            IF ( FLANDG_U(I,J) >  0.0 ) THEN
              TAUX_LAND(I,J) = ( GAMMA2*(TAUX_LAND(I,J)+                &
     &                               RHOKM_U_1(I,J)*DU_STAR1(I,J)) +    &
     &                 GAMMA1*RHOKM_U_1(I,J)*DU_1(I,J) ) /              &
     &                ( 1.0 + GAMMA1*RHOKM_U_1(I,J)*CQ_CM_U_1(I,J) )
            ELSE
              TAUX_LAND(I,J) = 0.0
            ENDIF
            IF ( FLANDG_U(I,J) <  1.0 ) THEN
              TAUX_SSI(I,J) = ( GAMMA2*(TAUX_SSI(I,J)+                  &
     &                              RHOKM_U_1(I,J)*DU_STAR1(I,J)) +     &
     &                 GAMMA1*RHOKM_U_1(I,J)*DU_1(I,J) ) /              &
     &                ( 1.0 + GAMMA1*RHOKM_U_1(I,J)*CQ_CM_U_1(I,J) )
            ELSE
              TAUX_SSI(I,J) = 0.0
            ENDIF
            TAUX_1(I,J) = FLANDG_U(I,J)*TAUX_LAND(I,J)                  &
     &                  + ( 1.0-FLANDG_U(I,J))*TAUX_SSI(I,J)
          ENDDO
        ENDDO

        DO J=1,N_ROWS
          DO I=1,ROW_LENGTH
            GAMMA1 = GAMMA1_IN(I,J)
            GAMMA2 = GAMMA2_IN(I,J)
            IF ( FLANDG_V(I,J) >  0.0 ) THEN
              TAUY_LAND(I,J) = ( GAMMA2*(TAUY_LAND(I,J)+                &
     &                               RHOKM_V_1(I,J)*DV_STAR1(I,J)) +    &
     &                 GAMMA1*RHOKM_V_1(I,J)*DV_1(I,J) ) /              &
     &                ( 1.0 + GAMMA1*RHOKM_V_1(I,J)*CQ_CM_V_1(I,J) )
            ELSE
              TAUY_LAND(I,J) = 0.0
            ENDIF
            IF ( FLANDG_V(I,J) <  1.0 ) THEN
              TAUY_SSI(I,J) = ( GAMMA2*(TAUY_SSI(I,J)+                  &
     &                               RHOKM_V_1(I,J)*DV_STAR1(I,J)) +    &
     &                 GAMMA1*RHOKM_V_1(I,J)*DV_1(I,J) ) /              &
     &                ( 1.0 + GAMMA1*RHOKM_V_1(I,J)*CQ_CM_V_1(I,J) )
            ELSE
              TAUY_SSI(I,J) = 0.0
            ENDIF  
            TAUY_1(I,J) = FLANDG_V(I,J)*TAUY_LAND(I,J)                  &
     &                  + ( 1.0-FLANDG_V(I,J))*TAUY_SSI(I,J)
          ENDDO
        ENDDO

      ENDIF ! L_CORRECT

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('IM_SF_PT2 ',4)
      ENDIF

      RETURN
      END SUBROUTINE IM_SF_PT2
