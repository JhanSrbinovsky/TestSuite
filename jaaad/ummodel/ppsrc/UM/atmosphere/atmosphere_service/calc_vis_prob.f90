
! SUBROUTINE CALC_VIS_PROB-------------------------------------------
!
! Purpose: Calculates visibility probability,
!          The visibility probability is similar to the cloud fraction
!          except it records the fraction of a grid box with RH
!          greater than that required for the critical visibility
!          (e.g 1 km), taking into account precipitation.
!
!          Suitable for single-column use.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.3   27/07/01   Original code.  Pete Clark.
!   5.5   17/04/03   Remove reference to obsolete section
!                    C90_1A. T.White
!   6.2   03/02/06   Move to section a71_1a. P.Selwood
!
!
! Programming standard:  Unified Model Documentation Paper No 3,
!                        Version 5, dated 08/12/92.
!
! Documentation:
!      Wright, B. J., 1997: Improvements to the Nimrod Visibility
!         Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.
!      Wright, B. J., 1997: A New Visibility Analysis/Forecast System
!         for Nimrod. Met. Office FR Tech Rep., No. 222.
!
!    NOTE: New UM Doc Paper to be produced soon (S.Cusack 5/11/01)
!
!END----------------------------------------------------------------
!
!
!Arguments:---------------------------------------------------------
      SUBROUTINE CALC_VIS_PROB(                                         &
     &  PSTAR,RHCRIT,LEVELS,POINTS,PFIELD                               &
                                                  !INPUT
     & ,T,AEROSOL,L_MURK,Q,QCL,QCF,VIS,NVIS                             &
                                                  !INPUT
     & ,LCA,CCA,PCT                                                     &
                                                  !INPUT
     & ,Beta_LS_Rain, Beta_LS_Snow                                      &
                                                  !INPUT
     & ,Beta_C_Rain, Beta_C_Snow                                        &
                                                  !INPUT
     & ,Prob_of_Vis                                                     &
                                                  !OUTPUT
     & ,ERROR                                                           &
                                                  !OUTPUT
     & )
      IMPLICIT NONE
      INTEGER                                                           &
     & LEVELS                                                           &
                           ! IN No. of levels being processed.
     &,POINTS                                                           &
                           ! IN No. of gridpoints being processed.
     &,PFIELD                                                           &
                           ! IN No. of points in global field (at one
!                          !    vertical level).
     &,NVIS                ! IN No. of visibility thresholds
      REAL                                                              &
     & PSTAR(PFIELD)                                                    &
                           ! IN Surface pressure (Pa).
     &,RHCRIT(LEVELS)                                                   &
                           ! IN Critical relative humidity.  See the
!                          !    the paragraph incorporating eqs P292.11
!                          !    to P292.14; the values need to be tuned
!                          !    for the given set of levels.
     &,Q(PFIELD,LEVELS)                                                 &
                           ! IN Specific Humidity
!                          !    (kg per kg air).
     &,QCL(PFIELD,LEVELS)                                               &
                           ! Cloud liquid water content at
!                          !     processed levels (kg per kg air).
     &,QCF(PFIELD,LEVELS)                                               &
                           ! Cloud ice content at processed levels
!                          !    (kg per kg air).
     &,T(PFIELD,LEVELS)                                                 &
                           ! IN Temperature (K).
     &,AEROSOL(PFIELD,LEVELS)                                           &
                              ! IN Aerosol mixing ratio(ug/kg)
     &,VIS(NVIS)                                                        &
                              ! Visibility thresholds
     &,LCA(PFIELD)                                                      &
                                       ! IN Total Layer Cloud.
     &,CCA(PFIELD)                                                      &
                                       ! IN Convective Cloud.
     &,Beta_LS_Rain(PFIELD,LEVELS)                                      &
                                       ! IN Scattering in LS Rain.
     &,Beta_LS_Snow(PFIELD,LEVELS)                                      &
                                       ! IN Scattering in LS Snow.
     &,Beta_C_Rain(PFIELD,LEVELS)                                       &
                                       ! IN Scattering in Conv Rain
     &,Beta_C_Snow(PFIELD,LEVELS)      ! IN Scattering in Conv Snow
      LOGICAL                                                           &
     &   L_MURK               ! IN : Aerosol present

      LOGICAL                                                           &
     & PCT                              ! IN T:Cloud amounts are in %

      REAL                                                              &
     & Prob_of_Vis(PFIELD,LEVELS,NVIS) ! OUT Vis prob at processed level
!                          !     (decimal fraction).
      INTEGER ERROR        ! OUT 0 if OK; 1 if bad arguments.
!
!---------------------------------------------------------------------
!  Workspace usage----------------------------------------------------
      REAL                                                              &
                           ! "Automatic" arrays on Cray.
     & Vis_Threshold(PFIELD,LEVELS,NVIS)                                &
     &,Prob_of_Vis_LS(PFIELD,LEVELS,NVIS)                               &
     &,Prob_of_Vis_C(PFIELD,LEVELS,NVIS)                                &
     &,P_LSP(PFIELD)                                                    &
     &,P_CP(PFIELD)                                                     &
     &,Temp ! Temporary
!  External subroutine called ----------------------------------------
      EXTERNAL FOG_FR
!
! Parameters
      REAL large_distance
!
      PARAMETER(                                                        &
     & large_distance = 1.0e6                                           &
     &         )
!
! Local, including SAVE'd, storage------------------------------------
!
      INTEGER K,I,J     ! Loop counters: K - vertical level index.
!                       !                I - horizontal field index.
                        !                J - Vis threshold index.
!--------------------------------------------------------------------
!  Local and other physical constants----------------------------------
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------
! C_VISBTY start
!LL Description:
!LL   This COMDECK contains declarations for constants used to diagnose
!LL visibility. Constants are set as PARAMTERs.
!LL
!LL
!LL  Model            Modification history:
!LL version  Date
!LL  3.2    29/04/93  CCN Parameters moved here from VISBTY so that
!LL                   they can also be used to compute fog fraction.
!LL                   Programmer: Pete Clark.
!LL  4.0 05/09/95  Variable AEROMAX used as upper limit to aerosol in
!LL                assimilation introduced. Programmer Pete Clark.
!LL  4.5 01/05/98  Completely re-written for NIMROD style diagnostic.
!    5.3 27/07/01  New parameters for the effect of precipitation
!                  visibility created.                 Pete Clark
!LL  5.3 17/10/01  Rename rho and rho_a. Adam Clayton
!LL
!LLEND----------------------------------------------------------------
      INTEGER, PARAMETER :: n_vis_thresh = 2

      ! Standard number density of the aerosol (/m3)
      REAL, PARAMETER :: N0 = 500.0E6

      ! Activation parameter
      REAL, PARAMETER :: B0= 0.5

      ! Radius of standard aerosol particle (m)
      REAL, PARAMETER :: radius0 = 0.16E-6

      REAL, PARAMETER :: FourThirds = 4.0/3.0     ! 4/3

      REAL :: vis_thresh(n_vis_thresh)
      DATA vis_thresh /1000.0,5000.0/
      ! Density of the the aerosol (Kg/m3)
      REAL, PARAMETER :: rho_aerosol = 1700.0

      ! Density of air (Kg/m3)
      REAL, PARAMETER :: rho_air = 1.0
      ! Standard aerosol mass mixing ratio (Kg/Kg)
      REAL, PARAMETER :: m0 = FourThirds * Pi *                         &
     &  radius0 * radius0 * radius0 *                                   &
     &                         (rho_aerosol/rho_air) * N0

      ! Aerosol particle radius/mass loading power
      REAL, PARAMETER :: power = 1.0/6.0

      ! Scattering coefficient normalisation
      REAL, PARAMETER :: Beta0  = 1.5 * Pi

      REAL, PARAMETER :: LiminalContrast   = 0.02

      ! Natural log of Liminal contrast
      REAL, PARAMETER :: LnLiminalContrast = -3.912023005

      ! Constant incorporating the scattering coefficient, normalisation
      ! transformation to visibility ( = ln(liminal contrast) / Beta0 )
      REAL, PARAMETER :: VisFactor = -LnLiminalContrast / Beta0

      ! Reciprocal of the clean air visibility
      REAL, PARAMETER :: RecipVisAir = 1.0E-5

      ! Constant involving surface energy of water
      REAL, PARAMETER :: A0 = 1.2E-9

      ! Visibility defining fog
      REAL, PARAMETER :: VISFOG = 1000.0

      ! Visibility defining mist
      REAL, PARAMETER :: VISMIST = 5000.0

      ! Minimum allowed aerosol
      REAL, PARAMETER :: AERO0 = 0.1

      ! maximum allowed aerosol
      REAL, PARAMETER :: AEROMAX = 200.0

      ! tunable parameter: the cumulative prob value at which vis is
      ! estimated
      REAL, PARAMETER :: calc_prob_of_vis = 0.4
!
!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
!-----------------------------------------------------------------------
!  Check input arguments for potential over-writing problems.
!-----------------------------------------------------------------------
      ERROR=0
      IF(POINTS >  PFIELD)THEN
        ERROR=1
        GOTO 9999
      ENDIF

      IF(PCT) THEN
        DO I=1,POINTS
          P_CP(I) =MAX(CCA(I)/100.0,0.0)
          P_LSP(I)=MAX((1.0-P_CP(I))*LCA(I)/100.0,0.0)
        ENDDO
      ELSE
        DO I=1,POINTS
          P_CP(I) =MAX(CCA(I),0.0)
          P_LSP(I)=MAX((1.0-P_CP(I))*LCA(I),0.0)
        ENDDO
      ENDIF


      DO K=1,LEVELS
        DO I=1,POINTS
          DO J=1,NVIS
            Vis_Threshold(I,K,J)=VIS(J)
          ENDDO
        ENDDO
      ENDDO

! DEPENDS ON: fog_fr
      CALL FOG_FR(                                                      &
     & PSTAR,RHCRIT,LEVELS,PFIELD,                                      &
     & T,AEROSOL,L_MURK,Q,QCL,QCF,Vis_Threshold,Prob_of_Vis,NVIS        &
     & )

      DO K=1,LEVELS
        DO I=1,POINTS
          DO J=1,NVIS
            IF (P_LSP(I)  >   0.0) THEN
              Temp = 1.0/VIS(J)+                                        &
     &         (Beta_LS_Rain(I,K)+Beta_LS_SNOW(I,K))/LnLiminalContrast
              IF (Temp  >   0.0) THEN
                Vis_Threshold(I,K,J)=MAX(1.0/Temp , VIS(J))
              ELSE
                Vis_Threshold(I,K,J)=large_distance
              ENDIF
            ELSE
              Vis_Threshold(I,K,J)=VIS(J)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

! DEPENDS ON: fog_fr
      CALL FOG_FR(                                                      &
     & PSTAR,RHCRIT,LEVELS,PFIELD,                                      &
     & T,AEROSOL,L_MURK,Q,QCL,QCF,Vis_Threshold,Prob_of_Vis_LS,NVIS     &
     & )

      DO K=1,LEVELS
        DO I=1,POINTS
          DO J=1,NVIS
            IF (P_CP(I)  >   0.0) THEN
              Temp = 1.0/VIS(J)+                                        &
     &         (Beta_C_Rain(I,K)+Beta_C_SNOW(I,K))/LnLiminalContrast
              IF (Temp  >   0.0) THEN
                Vis_Threshold(I,K,J)=MAX(1.0/Temp , VIS(J))
              ELSE
                Vis_Threshold(I,K,J)=large_distance
              ENDIF
            ELSE
              Vis_Threshold(I,K,J)=VIS(J)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

! DEPENDS ON: fog_fr
      CALL FOG_FR(                                                      &
     & PSTAR,RHCRIT,LEVELS,PFIELD,                                      &
     & T,AEROSOL,L_MURK,Q,QCL,QCF,Vis_Threshold,Prob_of_Vis_C,NVIS      &
     & )

      DO K=1,LEVELS
        DO I=1,POINTS
          DO J=1,NVIS
            Prob_of_Vis(I,K,J)=(1.0-P_CP(I)-P_LSP(I))*                  &
     &                            Prob_of_Vis(I,K,J) +                  &
     &                          P_LSP(I)*                               &
     &                            Prob_of_Vis_LS(I,K,J) +               &
     &                          P_CP(I)*                                &
     &                            Prob_of_Vis_C(I,K,J)
          ENDDO
        ENDDO
      ENDDO
!
 9999 Continue

      RETURN
      END SUBROUTINE CALC_VIS_PROB
