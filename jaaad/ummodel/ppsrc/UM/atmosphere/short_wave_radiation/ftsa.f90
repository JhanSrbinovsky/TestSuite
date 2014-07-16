
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine FTSA -------------------------------------------------
!LL
!LL  It calculates (true) surface albedos for P234.
!LL    Release 2.8 of the UM allows for separate surface
!LL  albedos for direct and diffuse light over sea, calculating the
!LL  former from the formula of Briegleb and Ramanathan (1982, J. Appl.
!LL  Met., 21, 1160-1171) and passes the albedos out in different form
!LL  Land & ice albedos are still the same for direct and diffuse beams.
!LL                                            William Ingram 25/9/92
!LL  Suitable for single column model use.
!LL
!LL   Author: William Ingram
!LL
!LL
!LL  It conforms to programming standard A of UMDP 4, version 2.
!LL  It contains ! comments, but otherwise conforms to the FORTRAN 77
!LL  standard with no features deprecated by 8X.
!LL
!LL Logical components covered : P233
!LL   (ancillary calculations for the shortwave scheme)
!LL
!LL Project task : P23
!LL
!LL  Offline documentation is in UMDP 23, sections "True surface albedo
!LL  specification" and "Modifications to the radiation scheme to
!LL  accommodate the leads model"
!LLEND
!*L
      SUBROUTINE FTSA (                                                 &
        LAND, FLANDG, AICE, TSTAR, TSTAR_SICE, COSZ, S, S_SEA,          &
        ALPHAM,ALPHAC,ALPHAB,DTICE,L_MOSES_II,L_SSICE_ALBEDO,           &
        L_MOD_BARKER_ALBEDO,                                            &
        L_SICE_MELTPONDS,L_SICE_SCATTERING,L_SICE_HADGEM1A,             &
        DT_BARE,DALB_BARE_WET,PEN_RAD_FRAC,BETA,VERSION,                &
        L1, L2, SA_LAND, SA_SICE, SAOS)
!
      implicit none
!
      CHARACTER                                                         &
                !, INTENT(IN) ::
     &     VERSION*3                      ! Version of radiation scheme
      INTEGER                                                           &
              !, INTENT(IN) ::
     &     L1,                                                          &
                                          ! Full field dimension
     &     L2                             ! Number of points to treat
      LOGICAL                                                           &
              !, INTENT(IN) ::
     &     LAND(L1)                                                     &
                                          ! Land-sea mask (land .TRUE.)
     &    ,L_MOSES_II                                                   &
                                          ! .TRUE. if MOSES II land
!                                         ! surface is selected.
!         Switch on the effect of snow on sea-ice albedo
     &    ,L_SSICE_ALBEDO                                               &
     &    ,L_MOD_BARKER_ALBEDO                                          &
                                          ! Use modified Barker
!                                         ! albedo (open sea).
     &    ,L_SICE_MELTPONDS                                             &
                             ! switch on seaice albedo meltponds
     &    ,L_SICE_SCATTERING                                            &
                              ! switch on seaice albedo internal scatter
     &    ,L_SICE_HADGEM1A   ! switch for HadGEM1 bug correction
      REAL                                                              &
           !, INTENT(IN) ::
     &     FLANDG(L1),                                                  &
                                          ! Land fraction
     &     AICE(L1),                                                    &
                                          ! Sea-ice fraction
     &     TSTAR(L1),                                                   &
                                          ! Surface temperature
     &     TSTAR_SICE(L1),                                              &
                                          ! Seaice surface temperature
     &     COSZ(L1),                                                    &
                                          ! cos(solar zenith angle)
     &     S(L1)                                                        &
                                          ! Snow amount (mass/area)
     &    ,S_SEA(L1)                                                    &
                             ! Snow amount on sea ice (mass/area of ice)
!     Constants used to determine the albedo of sea-ice:
! Albedo of sea-ice at melting point (TM) if .not.l_ssice_albedo, or
! Albedo of snow on sea-ice at melting point (TM) if l_ssice_albedo
     &    ,ALPHAM                                                       &
                   ! "M" for "melting"
! Albedo of sea-ice at and below TM-DTICE if .not.l_ssice_albedo, or
! Albedo of snow on sea-ice at and below TM-DTICE if l_ssice_albedo
     &    ,ALPHAC                                                       &
                   ! "C" for "cold"
! Albedo of snow-free sea-ice if l_ssice_albedo
     &    ,ALPHAB                                                       &
                   ! "B" for "bare"
! Temperature range in which albedo of sea-ice, if .not.l_ssice_albedo,
! or of snow on sea-ice, if l_ssice_albedo, varies between its limits
     &    ,DTICE                                                        &
! Temperature range below TM over which meltponds form if
! l_sice_meltponds and l_ssice_albedo
     &    ,DT_BARE                                                      &
! Increment to albedo for each degree temperature rises above
! TM-DT_BARE. Only used if l_sice_meltponds and l_ssice_albedo
     &    ,DALB_BARE_WET                                                &
! Fraction of SW radiation that penetrates seaice and scatters
! back causing an addition to the albedo. Only active if l_ssice_albedo
! and l_sice_scattering
     &    ,PEN_RAD_FRAC                                                 &
! attenutation parameter for SW in seaice which controls the
! additional albedo due to internal scattering
     &   ,BETA
      REAL                                                              &
          !, INTENT(OUT)
     &     SA_LAND(L1),                                                 &
                                          ! Surface Albedos for Land.
!                                         ! (Not output for MOSESII).
     &     SA_SICE(L1),                                                 &
                                          ! Surface Albedos for seaice
     &     SAOS(L1,2)                     !  and Ice, and for Open Sea,
!     ! respectively, with zeroes for safety where no value applies
!
!     !  FTSA has no dynamically allocated workspace, no EXTERNAL calls
!     !  and no significant structure - just one loop and some IF blocks
!*
      INTEGER J                           ! Loops over points
      REAL DSA                           ! Deep-snow albedo (alphasubD)
      REAL BX                                                           &
                                         ! 1 - COSZ

! Temperature at which (snow on) sea-ice reaches its "cold" value
     &    ,TCICE                                                        &
! Slope and intercept of temperature dependence of the albedo of
! (snow on) sea-ice
     &    ,ICE1,ICE2
!     Local parameters
      REAL DTLAND, KLAND, TCLAND, ADIFC, FCATM,                         &
     &     SNOW_ALBEDO,                                                 &
                                          ! Snow albedo
     &     MASKD                          ! Masking depth (S in 3.6.1)
!     Note that the same masking depth is always used, both for land,
!     regardless of vegetation cover, and for sea-ice. This assumption
!     of constancy may be doubtful.
      PARAMETER ( MASKD = 0.2 )
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
!     ! Basic quantities for land CSSA calculations:
      PARAMETER ( DTLAND = 2.,                                          &
                                          ! delta(T) in 3.6.2
     &     FCATM = 0.3 )                  ! Fraction by which deep-snow
!     ! albedo changes (from "cold" value towards snow-free value) at TM
!     ! From these, 2 constants precalculated for efficiency in 3.6.2:
      PARAMETER ( KLAND = 0.3/DTLAND,                                   &
     &     TCLAND = TM-DTLAND )
!
      PARAMETER ( ADIFC = 0.06 )          ! Surface albedo of ice-free
!                                         !    sea for the diffuse beam
!     ! derive 3 constants from the basic quantities (supplied in the
!     ! namelist RUNCNST) for sea-ice CSSA calculations:
      TCICE = TM - DTICE
      ICE1 = (ALPHAM-ALPHAC) / DTICE
      ICE2 = ALPHAM - TM*ICE1

      DO J=1, L2
        SA_LAND(J)=0.0
        SA_SICE(J)=0.0
      ENDDO

      DO 100 J=1, L2
        IF (FLANDG(J) <  1.0) THEN
          IF ((VERSION /= '03C').AND.(VERSION /= '03Z').AND.            &
     &                               (VERSION /= '03A') ) THEN
             SAOS(J,1) = 0.05 / ( 1.1 * COSZ(J)**1.4 + 0.15 )
          ELSE IF (L_MOD_BARKER_ALBEDO) THEN
             bx=1.0-COSZ(J)
             SAOS(J, 1)= 0.0315 + 0.0421*bx**2                          &
     &             + 0.128*bx**3 - 0.04*bx**4                           &
     &             + (3.12/(5.68) + (0.074*bx/(1.0)))*bx**5
          ELSE
             SAOS(J, 1)=0.026/(COSZ(J)**1.7+0.065)                      &
     &          +0.15*(COSZ(J)-0.1)*(COSZ(J)-0.5)*(COSZ(J)-1.0)
          ENDIF
          SAOS(J,2) = ADIFC
!         ! Note that the following will add in ICE1*(TSTAR-TFS) to CSSA
!         ! if AICE#0 when it should be - even if only very small: for
!         ! large enough TSTAR this will give very large surface heating
!         ! and even negative atmospheric heating.  Check whether this
!         ! could occur.
          IF ( AICE(J)  ==  0. ) THEN
             SA_SICE(J) = 0.
          ELSE
             if (l_ssice_albedo) then
               if (s_sea(j) >  0.0) then   ! snow on sea ice
                                           ! Cox et al., Clim Dyn,1999
                 if (tstar_sice(j) >  tcice) then
                   snow_albedo=ice2+ice1*tstar_sice(j)
                 else
                   snow_albedo=alphac
                 endif
                 sa_sice(j)=alphab                                      &
     &           +(snow_albedo-alphab)*(1.0-exp(-maskd*s_sea(j)))
               else           ! no snow so bare ice only
                 if(l_sice_meltponds) then
                    ! bare ice, temperature dep. (Ebert &Curry,1993)
                    if (TSTAR_SICE(j)  >   (TM-DT_BARE))then
                      if(l_sice_hadgem1a) then
                        sa_sice(j)= alphab+(TSTAR_SICE(j)-TM+DT_BARE)   &
     &                                 *DALB_BARE_WET
                      else    ! Incorrect version used in HadGEM1
                        sa_sice(j)= alphab+(TM - TSTAR_SICE(j))         &
     &                                 *DALB_BARE_WET
                      endif   ! l_sice_hadgem1a
                    else      ! surface is dry
                      sa_sice(j)=alphab
                    endif     ! end melt ponds
                 else         ! just use bare ice albedo
                   sa_sice(j)=alphab
                 endif        ! l_sice_meltponds
                 if(l_sice_scattering) then
                     !Semtner modification dry albedo for internal
                     ! scattering (Semnter, 1976)
                     sa_sice(j)=sa_sice(j)+BETA*(1.0-sa_sice(j))        &
     &                              *PEN_RAD_FRAC
                 endif      ! l_sice_scattering
               endif        !  any snow on ice
             else           ! default to operational NWP scheme
!            !  3.5.1:
             IF ( TSTAR_SICE(J)  <   TCICE ) THEN
                SA_SICE(J) = ALPHAC
              ELSE
                SA_SICE(J) = ICE1 * TSTAR_SICE(J) + ICE2
             ENDIF
             endif
          ENDIF
       ENDIF

  100 CONTINUE

!
      RETURN
      END SUBROUTINE FTSA
