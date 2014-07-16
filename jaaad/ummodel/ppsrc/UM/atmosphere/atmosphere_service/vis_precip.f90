
!   5.5   17/04/03   Remove reference to obsolete section
!                    C90_1A. T.White
!  SUBROUTINE VIS_PRECIP ---------------------------------------------
!
!     PURPOSE:
! Process fields of precipitation intensity to give scattering coefft
! in 1/metres.
! Calculated at model level (eg bottom eta level 25m)
! or level within surface layer eg screen ht ( 1.5M )
!
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.3   27/07/01   Original code.  Pete Clark.
!   6.2   03/02/06   Moved to a71_1a. P.Selwood
!   6.2   22/08/05   Remove spaces from GOTOs. P.Selwood
!
!  Programming standard: U M Doc. Paper No. 4
!
!  Logical components covered :
!
!  Project task:
!
!  External documentation
!    Forecasting Research Scientific Paper NO.4
!    Diagnosis of visibility in the UK Met Office Mesoscale Model
!    and the use of a visibility analysis to constrain initial
!    conditions.  SP Ballard, BJ Wright, BW Golding    1992
!      NIMROD diagnostic:
!      Wright, B. J., 1997: Improvements to the Nimrod Visibility
!         Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.
!      Wright, B. J., 1997: A New Visibility Analysis/Forecast System
!         for Nimrod. Met. Office FR Tech Rep., No. 222.
!
!    NOTE: New UM Doc Paper to be produced soon (S.Cusack 5/11/01)
!
!END----------------------------------------------------------------
!
!  Arguments:-------------------------------------------------------
      SUBROUTINE VIS_PRECIP                                             &
     &           (Vis_No_Precip                                         &
                                                      !INPUT
     &           ,LCA,CCA,PCT                                           &
                                                      !INPUT
     &           ,Beta_LS_Rain, Beta_LS_Snow                            &
                                                      !INPUT
     &           ,Beta_C_Rain, Beta_C_Snow                              &
                                                      !INPUT
     &           ,P_FIELD,POINTS,K1STPT                                 &
                                                      !INPUT
     &           ,Vis_overall,Vis_LSP,Vis_CP                            &
                                                      !OUTPUT
     &           ,ERROR)                              !OUTPUT
      IMPLICIT NONE
!---------------------------------------------------------------------
! Workspace usage:----------------------------------------------------
! 3 real arrays of size P_FIELD
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! input variables-----------------------------------------------------
!---------------------------------------------------------------------
      INTEGER                                                           &
     & P_FIELD                                                          &
                                        ! IN NO. points in field.
     &,POINTS                                                           &
                ! IN Number of gridpoints being processed.
     &,K1STPT                                                           &
                ! IN First gridpoint processed within complete field.
     &,ERROR    ! OUT Error code
      REAL                                                              &
     & Vis_No_Precip(P_FIELD)                                           &
                                        ! IN Vis outside precip.
     &,LCA(P_FIELD)                                                     &
                                        ! IN Total Layer Cloud.
     &,CCA(P_FIELD)                                                     &
                                        ! IN Convective Cloud.
     &,Beta_LS_Rain(P_FIELD)                                            &
                                        ! IN Scattering in LS Rain.
     &,Beta_LS_Snow(P_FIELD)                                            &
                                        ! IN Scattering in LS Snow.
     &,Beta_C_Rain(P_FIELD)                                             &
                                        ! IN Scattering in Conv Rain
     &,Beta_C_Snow(P_FIELD)             ! IN Scattering in Conv Snow

      LOGICAL                                                           &
     & PCT                              ! IN T:Cloud amounts are in %
!---------------------------------------------------------------------
! output variables----------------------------------------------------
!---------------------------------------------------------------------
      REAL                                                              &
     & Vis_overall(P_FIELD)                                             &
                                       ! OUT Visibility overall
     &,Vis_LSP(P_FIELD)                                                 &
                                       ! OUT Visibility in LS Precip.
     &,Vis_CP(P_FIELD)                 ! OUT Visibility in Conv Precip.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Local varables:------------------------------------------------------
!  Define local variables ---------------------------------------------
      INTEGER I       ! Loop counters: I - horizontal field index;

      REAL                                                              &
     & Beta_No_Precip                                                   &
     &,P_LSP(P_FIELD)                                                   &
     &,P_CP(P_FIELD)
!---------------------------------------------------------------------
!  External subroutine called ----------------------------------------
!---------------------------------------------------------------------
! Local and other physical constants----------------------------------
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
      ERROR=0
      IF((K1STPT+POINTS-1) >  P_FIELD)THEN
        ERROR=1
        GOTO 9999
      ENDIF
      IF(PCT) THEN
        DO I=K1STPT,K1STPT+POINTS-1
          P_CP(I)=CCA(I)/100.0
          P_LSP(I)=(1.0-P_CP(I))*LCA(I)/100.0
        ENDDO
      ELSE
        DO I=K1STPT,K1STPT+POINTS-1
          P_CP(I)=CCA(I)
          P_LSP(I)=(1.0-P_CP(I))*LCA(I)
        ENDDO
      ENDIF


      DO I=K1STPT,K1STPT+POINTS-1

        Beta_No_Precip=-LnLiminalContrast/Vis_No_Precip(I)

        IF(P_LSP(I)  >   0.0) THEN
          Vis_LSP(I) = -LnLiminalContrast /                             &
     &      (Beta_No_Precip +                                           &
     &       Beta_LS_Rain(I) + Beta_LS_Snow(I))
        ELSE
          Vis_LSP(I)=Vis_No_Precip(I)
        ENDIF

        IF(P_CP(I)  >   0.0) THEN
          Vis_CP(I) = -LnLiminalContrast /                              &
     &      (Beta_No_Precip +                                           &
     &       Beta_C_Rain(I) + Beta_C_Snow(I))
        ELSE
          Vis_CP(I)=Vis_No_Precip(I)
        ENDIF

! Ensure no rounding problems lead to vis > vis in clear air
        Vis_overall(I) = MIN((1.0-P_CP(I)-P_LSP(I))*Vis_No_Precip(I) +  &
     &                   P_LSP(I)*Vis_LSP(I) +P_CP(I)*Vis_CP(I),        &
     &                   Vis_No_Precip(I))

      ENDDO

 9999 Continue

      RETURN
      END SUBROUTINE VIS_PRECIP
