
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!-----------------------------------------------------------------------
!
! Subroutines SF_FLUX_LAND and SF_FLUX_SEA to calculate explicit surface
! fluxes of heat and moisture
!
!
!    Model            Modification history
!   version  date
!    5.2   15/11/00   New Deck         M. Best
!    5.3  25/04/01  Add coastal tiling. Nic Gedney
!  6.1  01/09/04  Calculate potential evaporation related variables.
!                                                          Nic Gedney
!    6.2  07/11/05  Allow for the salinity of sea water in
!                     the evaporative flux.
!                     J. M. Edwards
!
!    Programming standard:
!
!-----------------------------------------------------------------------

!     SUBROUTINE SF_FLUX_LAND-------------------------------------------
!
!     Calculate explicit surface fluxes of heat and moisture over
!     land tiles
!
!     ------------------------------------------------------------------

!     SUBROUTINE SF_FLUX_SEA--------------------------------------------
!
!     Calculate explicit surface fluxes of heat and moisture over sea
!     and sea-ice
!
!     ------------------------------------------------------------------
      SUBROUTINE SF_FLUX_SEA (                                          &
     & ROW_LENGTH,ROWS,NSICE,SICE_INDEX,FLANDG,                         &
     & ICE_FRACT,QS1,QSTAR_ICE,QSTAR_SEA,QW_1,RADNET,RHOKH_1,TI,        &
     & TL_1,TSTAR_SICE,TSTAR_SEA,Z0H_ICE,Z0M_ICE,Z0H_SEA,Z0M_SEA,Z1_TQ, &
     & SeaSalinityFactor,                                               &
     & ALPHA1,ASHTF,E_SEA,FQW_ICE,FQW_1,FTL_ICE,FTL_1,H_SEA,RHOKPM,     &
     & LTIMER)

      IMPLICIT NONE

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                            ! IN Number of X points?
     &,ROWS                                                             &
                            ! IN Number of Y points?
     &,NSICE                                                            &
                            ! IN Number of sea-ice points.
     &,SICE_INDEX(ROW_LENGTH*ROWS,2)
!                           ! IN Index of sea-ice points

      LOGICAL                                                           &
     & LTIMER                      ! IN  Logical for TIMER

      REAL                                                              &
     & FLANDG(ROW_LENGTH,ROWS)                                          &
     &,ICE_FRACT(ROW_LENGTH,ROWS)                                       &
                                   ! IN Fraction of gridbox which is
!                                  !    sea-ice.
     &,QS1(ROW_LENGTH,ROWS)                                             &
                                   ! IN Sat. specific humidity
!                                  !    qsat(TL_1,PSTAR)
     &,QSTAR_ICE(ROW_LENGTH,ROWS)                                       &
                                   ! IN Surface qsat for sea-ice.
     &,QSTAR_SEA(ROW_LENGTH,ROWS)                                       &
                                   ! IN Surface qsat for sea or sea-ice
!                                  !    leads.
     &,QW_1(ROW_LENGTH,ROWS)                                            &
                                   ! IN Total water content of lowest
!                                  !    atmospheric layer
!                                  !    (kg per kg air).
     &,RADNET(ROW_LENGTH,ROWS)                                          &
                                   ! IN Net surface radiation (W/m2)
!                                  !    positive downwards
     &,RHOKH_1(ROW_LENGTH,ROWS)                                         &
                                   ! IN Surface exchange coefficient.
     &,TI(ROW_LENGTH,ROWS)                                              &
                                   ! IN Temperature of sea-ice surface
!                                  !    layer (K)
     &,TL_1(ROW_LENGTH,ROWS)                                            &
                                   ! IN Liquid/frozen water temperature
!                                  !    for lowest atmospheric layer (K)
     &,TSTAR_SICE(ROW_LENGTH,ROWS)                                      &
                                    ! IN Sea-ice surface temperature (K)
     &,TSTAR_SEA(ROW_LENGTH,ROWS)                                       &
                                   ! IN Sea surface temperature (K).
     &,Z0H_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! IN Sea-ice heat and moisture
!                                  !    roughness length (m).
     &,Z0M_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! IN Sea-ice momentum roughness
!                                  !    length (m).
     &,Z0H_SEA(ROW_LENGTH,ROWS)                                         &
                                   ! IN Sea and lead heat and moisture
!                                  !    roughness length (m).
     &,Z0M_SEA(ROW_LENGTH,ROWS)                                         &
                                   ! IN Sea and lead momentum roughness
!                                  !    length.
     &,Z1_TQ(ROW_LENGTH,ROWS)      ! IN Height of lowest atmospheric
!                                  !    level (m).

!
!
      REAL, Intent(IN) :: SeaSalinityFactor
!                                  ! Factor allowing for the
!                                  ! effect of the salinity of
!                                  ! sea water on the evaporative
!                                  ! flux.
      REAL                                                              &
     & ALPHA1(ROW_LENGTH,ROWS)                                          &
                                   ! OUT Gradient of saturated specific
!                                  !     humidity with respect to
!                                  !     temperature between the bottom
!                                  !     model layer and the surface.
     &,ASHTF(ROW_LENGTH,ROWS)                                           &
                                   ! OUT Coefficient to calculate
!                                  !     surface heat flux into sea-ice
!                                  !     (W/m2/K).
     &,E_SEA(ROW_LENGTH,ROWS)                                           &
                                   ! OUT Evaporation from sea times
!                                  !     leads fraction (kg/m2/s).
     &,FQW_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! OUT Surface flux of QW for sea-ice.
     &,FQW_1(ROW_LENGTH,ROWS)                                           &
                                   ! OUT GBM surface flux of QW
!                                  !     (kg/m2/s).
     &,FTL_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! OUT Surface flux of TL for sea-ice.
     &,FTL_1(ROW_LENGTH,ROWS)                                           &
                                   ! OUT GBM surface flux of TL.
     &,H_SEA(ROW_LENGTH,ROWS)                                           &
                                   ! OUT Surface sensible heat flux over
!                                  !     sea times leads fraction (W/m2)
     &,RHOKPM(ROW_LENGTH,ROWS)     ! OUT Modified surface exchange
!                                  !     coefficient.

!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
! C_KAPPAI start

! Thermal conductivity of sea-ice (W per m per K).
        REAL,PARAMETER:: KAPPAI=2.09

! Thermal conductivity of sea water (W per m per K).
        REAL,PARAMETER:: kappas=0.31

! Snow density (Kg per m**3)
        REAL,PARAMETER:: rhosnow=330.0

! Effective thickness of sea-ice surface layer (m).
        REAL,PARAMETER:: DE = 0.1

! C_KAPPAI end
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
! CSIGMA start
      ! Stefan-Boltzmann constant (W/m**2/K**4).
      REAL, PARAMETER ::  SBCON=5.67E-8
! CSIGMA end

! Derived local parameters.
      REAL GRCP,LS
      PARAMETER (                                                       &
     & GRCP=G/CP                                                        &
     &,LS=LF+LC                                                         &
                           ! Latent heat of sublimation.
     & )

! Scalars
      INTEGER                                                           &
     & I,J                                                              &
                           ! Horizontal field index.
     &,K                   ! Sea-ice field index.
      REAL                                                              &
     & DQ1                                                              &
                           ! (qsat(TL_1,PSTAR)-QW_1) + g/cp*alpha1*Z1
     &,D_T                                                              &
                           ! Temporary in calculation of alpha1.
     &,RAD_REDUC           ! Radiation term required for surface flux
!                          ! calcs.

      EXTERNAL TIMER

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SF_FLUX ',3)
      ENDIF

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        ALPHA1(I,J) = 0.
        E_SEA(I,J) = 0.
        H_SEA(I,J) = 0.
        FQW_ICE(I,J) = 0.
        FTL_ICE(I,J) = 0.
        RHOKPM(I,J) = 0.
       ENDDO
      ENDDO

!----------------------------------------------------------------------
!!  1 Calculate gradient of saturated specific humidity for use in
!!    calculation of surface fluxes - only required for sea-ice points
!----------------------------------------------------------------------
      DO K=1,NSICE
        I = SICE_INDEX(K,1)
        J = SICE_INDEX(K,2)
        D_T = TSTAR_SICE(I,J) - TL_1(I,J)
        IF (D_T  >   0.05 .OR. D_T  <   -0.05) THEN
          ALPHA1(I,J) = (QSTAR_ICE(I,J) - QS1(I,J)) / D_T
        ELSEIF (TL_1(I,J)  >   TM) THEN
          ALPHA1(I,J) = EPSILON*LC*QS1(I,J)*( 1.0+C_VIRTUAL*QS1(I,J) )  &
     &                                        /(R*TL_1(I,J)*TL_1(I,J))
        ELSE
          ALPHA1(I,J) = EPSILON*LS*QS1(I,J)*( 1.0+C_VIRTUAL*QS1(I,J) )  &
     &                                        /(R*TL_1(I,J)*TL_1(I,J))
        ENDIF
        ASHTF(I,J) = 2 * KAPPAI / DE + 4*SBCON*TI(I,J)**3
      ENDDO

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        IF (FLANDG(I,J) <  1.0 ) THEN

          E_SEA(I,J) = - (1. - ICE_FRACT(I,J)) *                        &
     &                        RHOKH_1(I,J)*(QW_1(I,J) -                 &
     &                        SeaSalinityFactor*QSTAR_SEA(I,J))
          H_SEA(I,J) = - (1. - ICE_FRACT(I,J))*CP*RHOKH_1(I,J) *        &
     &                 ( TL_1(I,J) - TSTAR_SEA(I,J)                     &
     &              + GRCP*(Z1_TQ(I,J) + Z0M_SEA(I,J) - Z0H_SEA(I,J)) )

          IF ( ICE_FRACT(I,J)  >   0. ) THEN
! Sea-ice
            RHOKPM(I,J) = RHOKH_1(I,J) / ( ASHTF(I,J) +                 &
     &                             RHOKH_1(I,J)*(LS*ALPHA1(I,J) + CP) )
            RAD_REDUC = RADNET(I,J) - ICE_FRACT(I,J) * ASHTF(I,J) *     &
     &                  ( TL_1(I,J) - TI(I,J) +                         &
     &                GRCP*(Z1_TQ(I,J) + Z0M_ICE(I,J) - Z0H_ICE(I,J)) )
            DQ1 = QS1(I,J) - QW_1(I,J) +                                &
     &                           GRCP*ALPHA1(I,J)*                      &
     &                       (Z1_TQ(I,J) + Z0M_ICE(I,J) - Z0H_ICE(I,J))
            FQW_ICE(I,J) = RHOKPM(I,J) * ( ALPHA1(I,J)*RAD_REDUC +      &
     &                        (CP*RHOKH_1(I,J) +                        &
     &                         ASHTF(I,J))*DQ1*ICE_FRACT(I,J) )
            FTL_ICE(I,J) = RHOKPM(I,J) * ( RAD_REDUC -                  &
     &                             ICE_FRACT(I,J)*LS*RHOKH_1(I,J)*DQ1 )

          ENDIF

          FTL_1(I,J) = (1.-FLANDG(I,J))*(FTL_ICE(I,J) + H_SEA(I,J) / CP)
          FQW_1(I,J) = (1.-FLANDG(I,J))*(FQW_ICE(I,J) + E_SEA(I,J))

        ENDIF
       ENDDO
      ENDDO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SF_FLUX ',4)
      ENDIF

      RETURN
      END SUBROUTINE SF_FLUX_SEA

