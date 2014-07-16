
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE IM_SF_PT ----------------------------------------------
!!!
!!!  Purpose: Calculate implicit increments to surface variables
!!!
!!!
!!!  Model           Modification history
!!! version  Date
!!!  5.2   15/11/00   New Deck         M. Best
!    5.3  25/04/01  Add coastal tiling. Nic Gedney
!!!  6.1  01/09/04  Pass potential evaporation related variables.
!                                                          Nic Gedney
!  6.2  07/11/05  Calculate land potential evaporation.
!                                                    Nic Gedney
!!!  6.2  02/02/06  Passes L_flux_bc through argument list to allow
!!!                 settings for prescribed surface flux forcing
!!!                                                           R. Wong
!!!
!!!  Programming standard: UM Documentation Paper No 4, Version 2,
!!!                        dated 18/1/90
!!!
!!!  System component covered: P244
!!!
!!!  Project task: P24
!!!
!!!  Documentation: UM Documentation Paper No 24.
!!!
!!!---------------------------------------------------------------------
!!  Arguments :-
      SUBROUTINE IM_SF_PT (                                             &
     & OFF_X,OFF_Y,ROW_LENGTH,ROWS,N_ROWS,LAND_PTS                      &
     &,LAND_INDEX,NTILES,TILE_INDEX,TILE_PTS                            &
     &,FLANDG,TILE_FRAC,SNOW_TILE,ICE_FRACT                             &
     &,GAMMA_IN,ALPHA1,ALPHA1_SICE,ASHTF,ASHTF_TILE                     &
     &,RESFT,RHOKPM,RHOKPM_POT,RHOKPM_SICE                              &
     &,RHOKM_U_1,RHOKM_V_1,RHOKH_1,RHOKH1_SICE                          &
     &,CT_CTQ_1,DQW_1,DTL_1,CQ_CM_U_1,CQ_CM_V_1,DU_1,DV_1               &
     &,FLANDG_U,FLANDG_V                                                &
     &,FQW_GB,FTL_GB                                                    &
     &,TAUX_1,TAUX_LAND,TAUX_SSI,TAUY_1,TAUY_LAND,TAUY_SSI              &
     &,FQW_TILE,EPOT_TILE,FTL_TILE,FQW_ICE,FTL_ICE,E_SEA,H_SEA          &
     &,L_FLUX_BC,LTIMER                                                 &
! EAK
     &,l_cable                                                          &
     &)


      IMPLICIT NONE

      LOGICAL LTIMER                                                    &
     &,L_FLUX_BC                  &  ! Logical for prescribed
! EAK
     &, l_cable 
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
                                   ! INOUT surface tile potential
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
     &,DU_1(1-OFF_X:ROW_LENGTH+OFF_X,                                   &
     &      1-OFF_Y:ROWS+OFF_Y)                                         &
                                   ! IN Level 1 increment to u wind
!                                  !    field
     &,DV_1(1-OFF_X:ROW_LENGTH+OFF_X,                                   &
     &      1-OFF_Y:N_ROWS+OFF_Y)                                       &
                                   ! IN Level 1 increment to v wind
!                                  !    field
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
                                   ! Tempary array
     &,BPART(ROW_LENGTH,ROWS,2)                                         &
                                   ! Tempary array
     &,RECIP(ROW_LENGTH,ROWS)                                           &
                                   ! Tempary array
     &,FTL_LAND(ROW_LENGTH,ROWS)                                        &
                                      ! Tempary array
     &,FQW_LAND(ROW_LENGTH,ROWS)      ! Tempary array

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
! EAK
!     &,GAMMA_1     ! local implicit weight
     &,GAMMA_1(ROW_LENGTH,ROWS) ! local implicit weight

      integer ktau  
      save ktau
      LOGICAL, save :: LFIRST = .true.
    

      IF ( LFIRST ) THEN
       ktau = 0
       LFIRST = .false.
      ENDIF
      ktau = ktau  + 1

!-----------------------------------------------------------------------
!!  0.  Check that the scalars input to define the grid are consistent.
!       See comments to routine SF_EXCH for details.
!-----------------------------------------------------------------------

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('IM_SF_PT ',3)
      ENDIF
      GAMMA_1=GAMMA_IN
! EAK
!      IF (L_FLUX_BC) THEN
!       GAMMA_1 = 0.0       ! use GAMMA=0 for scalars (specified fluxes)
!      ENDIF
! EAK 
!     GAMMA_1 = 0.0  !*** for CABLE ***
!      IF( l_cable ) GAMMA_1 = 0.0 

! Initialise APART and BPART to zero
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
! EAK
        GAMMA_1(I,J)=GAMMA_IN
        IF(L_FLUX_BC) GAMMA_1(I,J)=0.0
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
!CDIR NODEP
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          LAT_HT = LC
          IF (SNOW_TILE(L,N) >  0.) LAT_HT = LS

! EAK     Gamma_1 set to zero for all land points (including land fraction <1)
          IF( l_cable ) GAMMA_1(I,J) = 0.0

          APART(I,J,1)=APART(I,J,1) - TILE_FRAC(L,N) *                  &
     &               GAMMA_1(I,J) * RHOKPM(L,N) *                       &
     &            ( LAT_HT*RESFT(L,N)*RHOKH_1(L,N)*ALPHA1(L,N) +        &
     &                         ASHTF_TILE(L,N) )
          APART(I,J,2)=APART(I,J,2) + TILE_FRAC(L,N) *                  &
     &               GAMMA_1(I,J) * RHOKPM(L,N) *                       &
     &               LAT_HT*RESFT(L,N)*RHOKH_1(L,N)
          BPART(I,J,1)=BPART(I,J,1) + TILE_FRAC(L,N) *                  &
     &               GAMMA_1(I,J) * RESFT(L,N)*RHOKPM(L,N) *            &
     &               CP*RHOKH_1(L,N)*ALPHA1(L,N)
          BPART(I,J,2)=BPART(I,J,2) - TILE_FRAC(L,N) *                  &
     &               GAMMA_1(I,J) * RESFT(L,N)*RHOKPM(L,N) *            &
     &               ( CP*RHOKH_1(L,N) + ASHTF_TILE(L,N) )

        ENDDO
      ENDDO

!      if( ktau .gt. 300) then
!      print 87,ftl_gb(4,2),fqw_gb(4,2), &
!      ftl_tile(100,9),fqw_tile(100,9),apart(4,2,1),apart(4,2,2), &
!      bpart(4,2,1),bpart(4,2,2)
!      print 88,ftl_gb(28,57),fqw_gb(28,57), &
!      ftl_tile(2132,1:5),ftl_tile(2132,8),fqw_tile(2132,1:5), &
!      fqw_tile(2132,8),apart(28,57,1),apart(28,57,2), &
!      bpart(28,57,1),bpart(28,57,2)
!      print 89,ftl_gb(35,56),fqw_gb(35,56), &
!      ftl_tile(2080,1:5),ftl_tile(2080,8),fqw_tile(2080,1:5), &
!      fqw_tile(2080,8),apart(35,56,1),apart(35,56,2), &
!      bpart(35,56,1),bpart(35,56,2)
! 87   format('imsfpt1a',4f9.6,x,4f7.2)
! 88   format('imsfpt2a',14f9.6,4f7.2)
! 89   format('imsfpt3a',14f9.6,4f7.2)
!      endif

! Sea points
      DO J=1,ROWS
       DO I=1,ROW_LENGTH

        IF(FLANDG(I,J) <  1.0 .AND. ICE_FRACT(I,J) >  0.0) THEN
! Sea ice point
          APART(I,J,1)=FLANDG(I,J)*APART(I,J,1)                         &
     &       - GAMMA_1(I,J) * (1.0-FLANDG(I,J)) * ICE_FRACT(I,J)        &
     &      * RHOKPM_SICE(I,J) *                                        &
     &      ( LS*RHOKH1_SICE(I,J)*ALPHA1_SICE(I,J) + ASHTF(I,J) )       &
     &       - GAMMA_1(I,J) * (1.0-FLANDG(I,J)) * (1.0 -ICE_FRACT(I,J)) &
     &      * RHOKH1_SICE(I,J)

          APART(I,J,2)=FLANDG(I,J)*APART(I,J,2)                         &
     &       + GAMMA_1(I,J) * (1.0-FLANDG(I,J)) *ICE_FRACT(I,J)         &
     &       * RHOKPM_SICE(I,J) * LS*RHOKH1_SICE(I,J)

          BPART(I,J,1)=FLANDG(I,J)*BPART(I,J,1)                         &
     &       + GAMMA_1(I,J) * ICE_FRACT(I,J) * ( 1.0 - FLANDG(I,J) )    &
     &       * RHOKPM_SICE(I,J) *CP*RHOKH1_SICE(I,J)*ALPHA1_SICE(I,J)

          BPART(I,J,2)=FLANDG(I,J)*BPART(I,J,2)                         &
     &       - GAMMA_1(I,J) * ICE_FRACT(I,J) * ( 1.0 - FLANDG(I,J) )    &
     &       * RHOKPM_SICE(I,J) * ( CP*RHOKH1_SICE(I,J) + ASHTF(I,J) )  &
     &       - GAMMA_1(I,J) * ( 1.0 - ICE_FRACT(I,J) )                  &
     &       * ( 1.0 - FLANDG(I,J) ) * RHOKH1_SICE(I,J)

        ELSEIF(FLANDG(I,J) <  1.0 .AND.                                 &
     &         .NOT.ICE_FRACT(I,J) >  0.0) THEN
! Ordinary sea point
          APART(I,J,1)= FLANDG(I,J)*APART(I,J,1)                        &
     &       - GAMMA_1(I,J) * ( 1.0 - FLANDG(I,J) ) * RHOKH1_SICE(I,J)
          APART(I,J,2)= FLANDG(I,J)*APART(I,J,2)

          BPART(I,J,1)= FLANDG(I,J)*BPART(I,J,1)
          BPART(I,J,2)= FLANDG(I,J)*BPART(I,J,2)                        &
     &       - GAMMA_1(I,J) * ( 1.0 - FLANDG(I,J) ) * RHOKH1_SICE(I,J)

        ENDIF
       ENDDO
      ENDDO

! Land tiles
      DO N=1,NTILES
!CDIR NODEP
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          E_EPOT_TILE(L,N)=1.0
          IF(EPOT_TILE(L,N) >  0.0.AND.FQW_TILE(L,N) >  0.0)            &
     &       E_EPOT_TILE(L,N)=EPOT_TILE(L,N)/FQW_TILE(L,N)
        ENDDO
      ENDDO


!      if( ktau .gt. 300) then
!
!      print 187,ftl_gb(4,2),fqw_gb(4,2), &
!      ftl_tile(100,9),fqw_tile(100,9),apart(4,2,1),apart(4,2,2), &
!      bpart(4,2,1),bpart(4,2,2)
!      print 188,ftl_gb(28,57),fqw_gb(28,57), &
!      ftl_tile(2132,1:5),ftl_tile(2132,8),fqw_tile(2132,1:5), &
!      fqw_tile(2132,8),apart(28,57,1),apart(28,57,2), &
!      bpart(28,57,1),bpart(28,57,2)
!      print 189,ftl_gb(35,56),fqw_gb(35,56), &
!      ftl_tile(2080,1:5),ftl_tile(2080,8),fqw_tile(2080,1:5), &
!      fqw_tile(2080,8),apart(35,56,1),apart(35,56,2), &
!      bpart(35,56,1),bpart(35,56,2)
!187   format('imsfpt1b',4f9.6,x,4f7.2)
!188   format('imsfpt2b',14f9.6,4f7.2)
!189   format('imsfpt3b',14f9.6,4f7.2)
!      endif


! Calculate grid-box fluxes of heat and moisture
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        RECIP(I,J)=( 1.0 + CT_CTQ_1(I,J)*APART(I,J,1) ) *               &
     &           ( 1.0 + CT_CTQ_1(I,J)*BPART(I,J,2) ) -                 &
     &           CT_CTQ_1(I,J)*APART(I,J,2)*CT_CTQ_1(I,J)*BPART(I,J,1)

        FTL_OLD=FTL_GB(I,J)

        FTL_GB(I,J) = ( ( 1.0 + CT_CTQ_1(I,J)*BPART(I,J,2) ) *          &
     &               ( FTL_OLD + APART(I,J,1)*DTL_1(I,J) +              &
     &                 APART(I,J,2)*DQW_1(I,J)) -                       &
     &                 CT_CTQ_1(I,J)*APART(I,J,2) * ( FQW_GB(I,J) +     &
     &                 BPART(I,J,1)*DTL_1(I,J) +                        &
     &                 BPART(I,J,2)*DQW_1(I,J)) ) / RECIP(I,J)

        FQW_GB(I,J) = ( ( 1.0 + CT_CTQ_1(I,J)*APART(I,J,1) ) *          &
     &                ( FQW_GB(I,J) + BPART(I,J,1)*DTL_1(I,J) +         &
     &                  BPART(I,J,2)*DQW_1(I,J)) -                      &
     &                  CT_CTQ_1(I,J)*BPART(I,J,1) * ( FTL_OLD +        &
     &                  APART(I,J,1)*DTL_1(I,J) +                       &
     &                  APART(I,J,2)*DQW_1(I,J)) ) / RECIP(I,J)

       ENDDO
      ENDDO
!      if( ktau .gt. 300) then
!      print 287,ftl_gb(4,2),fqw_gb(4,2), &
!      ftl_tile(100,9),fqw_tile(100,9),RECIP(4,2),ct_ctq_1(4,2), &
!      dtl_1(4,2),dqw_1(4,2) 
!      print 288,ftl_gb(28,57),fqw_gb(28,57), &
!      ftl_tile(2132,1:5),ftl_tile(2132,8),fqw_tile(2132,1:5), &
!      fqw_tile(2132,8),recip(28,57),ct_ctq_1(28,57), &
!      dtl_1(28,57),dqw_1(28,57)
!      print 289,ftl_gb(35,56),fqw_gb(35,56), &
!      ftl_tile(2080,1:5),ftl_tile(2080,8),fqw_tile(2080,1:5), &
!      fqw_tile(2080,8),recip(35,56),ct_ctq_1(35,56), &
!      dtl_1(35,56),dqw_1(35,56)
!287   format('imsfpt1c',4f9.6,x,4f9.4)
!288   format('imsfpt2c',14f9.6,4f9.4)
!289   format('imsfpt3c',14f9.6,4f9.4)
!      endif


! Make implicit correction to tile fluxes

! Land tiles
      DO N=1,NTILES
!CDIR NODEP
        DO K=1,TILE_PTS(N)
          L = TILE_INDEX(K,N)
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          LAT_HT = LC
          IF (SNOW_TILE(L,N) >  0.) LAT_HT = LS

          FTL_TILE(L,N)=FTL_TILE(L,N) -                                 &
     &               GAMMA_1(I,J) * RHOKPM(L,N) *                       &
     &            ( LAT_HT*RESFT(L,N)*RHOKH_1(L,N)*ALPHA1(L,N) +        &
     &                       ASHTF_TILE(L,N) ) *                        &
     &         ( DTL_1(I,J) - CT_CTQ_1(I,J)*FTL_GB(I,J) ) +             &
     &               GAMMA_1(I,J) * RHOKPM(L,N) *                       &
     &               LAT_HT*RESFT(L,N)*RHOKH_1(L,N) *                   &
     &         ( DQW_1(I,J) - CT_CTQ_1(I,J)*FQW_GB(I,J) )

          FQW_TILE(L,N)=FQW_TILE(L,N) +                                 &
     &               GAMMA_1(I,J) * RESFT(L,N)*RHOKPM(L,N) *            &
     &               CP*RHOKH_1(L,N)*ALPHA1(L,N) *                      &
     &         ( DTL_1(I,J) - CT_CTQ_1(I,J)*FTL_GB(I,J) ) -             &
     &               GAMMA_1(I,J) * RESFT(L,N)*RHOKPM(L,N) *            &
     &               ( CP*RHOKH_1(L,N) + ASHTF_TILE(L,N) ) *            &
     &         ( DQW_1(I,J) - CT_CTQ_1(I,J)*FQW_GB(I,J) )


          EPOT_TILE(L,N)=FQW_TILE(L,N)*E_EPOT_TILE(L,N)

          FQW_LAND(I,J)=FQW_LAND(I,J)+FQW_TILE(L,N)*TILE_FRAC(L,N)
          FTL_LAND(I,J)=FTL_LAND(I,J)+FTL_TILE(L,N)*TILE_FRAC(L,N)

        ENDDO
      ENDDO


!      if( ktau .gt. 300) then
!      print 387,ftl_gb(4,2),fqw_gb(4,2), &
!      ftl_tile(100,9),fqw_tile(100,9),apart(4,2,1),apart(4,2,2), &
!      bpart(4,2,1),bpart(4,2,2)
!      print 388,ftl_gb(28,57),fqw_gb(28,57), &
!      ftl_tile(2132,1:5),ftl_tile(2132,8),fqw_tile(2132,1:5), &
!      fqw_tile(2132,8),apart(28,57,1),apart(28,57,2), &
!      bpart(28,57,1),bpart(28,57,2)
!      print 389,ftl_gb(35,56),fqw_gb(35,56), &
!      ftl_tile(2080,1:5),ftl_tile(2080,8),fqw_tile(2080,1:5), &
!      fqw_tile(2080,8),apart(35,56,1),apart(35,56,2), &
!      bpart(35,56,1),bpart(35,56,2)
!387   format('imsfpt1c',4f9.6,x,4f7.2)
!388   format('imsfpt2c',14f9.6,4f7.2)
!389   format('imsfpt3c',14f9.6,4f7.2)
!      endif



! Sea points
      DO J=1,ROWS
       DO I=1,ROW_LENGTH

        IF(FLANDG(I,J) <  1.0 .AND. ICE_FRACT(I,J) >  0.0) THEN
! Sea ice point
          H_SEA(I,J)=H_SEA(I,J) - GAMMA_1(I,J) * (1.0-ICE_FRACT(I,J))*  &
     &                   CP * RHOKH1_SICE(I,J) *                        &
     &                     ( DTL_1(I,J) - CT_CTQ_1(I,J) * FTL_GB(I,J) )
          E_SEA(I,J)=E_SEA(I,J) - GAMMA_1(I,J) * (1.0-ICE_FRACT(I,J))*  &
     &                   RHOKH1_SICE(I,J) *                             &
     &                     ( DQW_1(I,J) - CT_CTQ_1(I,J) * FQW_GB(I,J) )
          FTL_ICE(I,J)=(FTL_GB(I,J)                                     &
     &       - FTL_LAND(I,J)*FLANDG(I,J))/(1.-FLANDG(I,J))              &
     &       - H_SEA(I,J)/CP
          FQW_ICE(I,J)=(FQW_GB(I,J)                                     &
     &       - FQW_LAND(I,J)*FLANDG(I,J))/(1.-FLANDG(I,J))              &
     &        - E_SEA(I,J)

        ELSEIF(FLANDG(I,J) <  1.0 .AND.                                 &
     &         .NOT.ICE_FRACT(I,J) >  0.0) THEN
! Ordinary sea point
          H_SEA(I,J)=CP * (FTL_GB(I,J)                                  &
     &       - FTL_LAND(I,J)*FLANDG(I,J))/(1.-FLANDG(I,J))
          E_SEA(I,J)=(FQW_GB(I,J)                                       &
     &       - FQW_LAND(I,J)*FLANDG(I,J))/(1.-FLANDG(I,J))
          FTL_ICE(I,J)=0.0
          FQW_ICE(I,J)=0.0

        ENDIF
       ENDDO
      ENDDO


      GAMMA_1=GAMMA_IN       ! use input GAMMA for winds

      DO J=1,ROWS
        DO I=1,ROW_LENGTH

        IF(FLANDG_U(I,J) >  0.0)THEN
          TAUX_LAND(I,J) = ( TAUX_LAND(I,J) +                           &
     &                 GAMMA_1(I,J)*RHOKM_U_1(I,J)*DU_1(I,J) ) /        &
     &               ( 1.0 + GAMMA_1(I,J)*RHOKM_U_1(I,J)*CQ_CM_U_1(I,J))
        ELSE
          TAUX_LAND(I,J) = 0.0
        ENDIF

        IF(FLANDG_U(I,J) <  1.0)THEN
          TAUX_SSI(I,J) = ( TAUX_SSI(I,J) +                             &
     &                 GAMMA_1(I,J)*RHOKM_U_1(I,J)*DU_1(I,J) ) /        &
     &               ( 1.0 + GAMMA_1(I,J)*RHOKM_U_1(I,J)*CQ_CM_U_1(I,J))
        ELSE
          TAUX_SSI(I,J) = 0.0
        ENDIF

        TAUX_1(I,J) = FLANDG_U(I,J)*TAUX_LAND(I,J)                      &
     &                + ( 1.0-FLANDG_U(I,J))*TAUX_SSI(I,J)

        ENDDO
      ENDDO


      DO J=1,N_ROWS
        DO I=1,ROW_LENGTH

        IF(FLANDG_V(I,J) >  0.0)THEN
          TAUY_LAND(I,J) = ( TAUY_LAND(I,J) +                           &
     &                 GAMMA_1(I,J)*RHOKM_V_1(I,J)*DV_1(I,J) ) /        &
     &              ( 1.0 + GAMMA_1(I,J)*RHOKM_V_1(I,J)*CQ_CM_V_1(I,J) )
        ELSE
          TAUY_LAND(I,J) = 0.0
        ENDIF

        IF(FLANDG_V(I,J) <  1.0)THEN
          TAUY_SSI(I,J) = ( TAUY_SSI(I,J) +                             &
     &                 GAMMA_1(I,J)*RHOKM_V_1(I,J)*DV_1(I,J) ) /        &
     &             ( 1.0 + GAMMA_1(I,J)*RHOKM_V_1(I,J)*CQ_CM_V_1(I,J) )
        ELSE
          TAUY_SSI(I,J) = 0.0
        ENDIF

        TAUY_1(I,J) = FLANDG_V(I,J)*TAUY_LAND(I,J)                      &
     &                + ( 1.0-FLANDG_V(I,J))*TAUY_SSI(I,J)

        ENDDO
      ENDDO

!      if( ktau .gt. 300) then
!      print 487,ftl_gb(4,2),fqw_gb(4,2), &
!      ftl_tile(100,9),fqw_tile(100,9),apart(4,2,1),apart(4,2,2),  &
!      bpart(4,2,1),bpart(4,2,2),taux_1(4,2),tauy_1(4,2)
!      print 488,ftl_gb(28,57),fqw_gb(28,57), &
!      ftl_tile(2132,1:5),ftl_tile(2132,8),fqw_tile(2132,1:5), &
!      fqw_tile(2132,8),apart(28,57,1),apart(28,57,2), &
!      bpart(28,57,1),bpart(28,57,2),taux_1(28,57),tauy_1(28,57)
!      print 489,ftl_gb(35,56),fqw_gb(35,56), &
!      ftl_tile(2080,1:5),ftl_tile(2080,8),fqw_tile(2080,1:5), &
!      fqw_tile(2080,8),apart(35,56,1),apart(35,56,2), &
!      bpart(35,56,1),bpart(35,56,2),taux_1(35,56),tauy_1(35,56)
!487   format('imsfpt1d',4f9.6,x,6f7.2)
!488   format('imsfpt2d',14f9.6,6f7.2)
!489   format('imsfpt3d',14f9.6,6f7.2)
!      endif




      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('IM_SF_PT ',4)
      ENDIF

      RETURN
      END SUBROUTINE IM_SF_PT
