

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate differences in source functions.
!
! Method:
!       Using the polynomial fit to the Planck function, values
!       of this function at the boundaries of layers are found
!       and differences across layers are determined. If the
!       Planckian is being taken to vary quadratically across
!       the layer second differences are found.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             14-03-95                Explicit formulation
!                                               for sea ice introduced.
!
!       5.3             25-04-01         Allow for land and sea
!                                        to co-exist in same gridbox.
!                                        Replace fractional sea-ice
!                                        points with fraction open-sea
!                                        points so that fractional land
!                                        is catered for. Replace TFS
!                                        with general sea temperature.
!                                        (N. Gedney)
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE DIFF_PLANCK_SOURCE(N_PROFILE, N_LAYER                  &
     &   , N_DEG_FIT, THERMAL_COEFFICIENT                               &
     &   , T_REF_PLANCK, T_LEVEL, T_GROUND                              &
     &   , T_SOLID, T_SEA, L_CTILE                                      &
     &   , PLANCK_SOURCE, DIFF_PLANCK, PLANCK_GROUND                    &
     &   , L_IR_SOURCE_QUAD, T, DIFF_PLANCK_2                           &
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION             &
     &   , FLANDG, PLANCK_FREEZE_SEA, PLANCK_LEADS_SEA                  &
     &   , NPD_PROFILE, NPD_LAYER, NPD_THERMAL_COEFF                    &
     &   )
!

      USE auscom_cpl_data_mod,                                          &
     &    Only : access_tfs

!
      IMPLICIT NONE
!
!
!     DUMMY ARRAY SIZES
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER                                                    &
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_THERMAL_COEFF
!             NUMBER OF THERMAL COEFFICIENTS
!
!     COMDECKS INCLUDED.
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , N_DEG_FIT
!             DEGREE OF FITTING FUNCTION
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_IR_SOURCE_QUAD                                             &
!             IR-SOURCE QUADRATIC
     &   , L_CTILE
!             Coastal tiling switch
      REAL                                                              &
                !, INTENT(IN)
     &     THERMAL_COEFFICIENT(0: NPD_THERMAL_COEFF-1)                  &
!             COEFFICIENTS OF FIT TO PLANCK FNC
     &   , T_REF_PLANCK                                                 &
!             PLANCKIAN REFERENCE TEMPERATURE
     &   , T_LEVEL(NPD_PROFILE, 0: NPD_LAYER)                           &
!             TEMPERATURES ON LEVELS
     &   , T_GROUND(NPD_PROFILE)                                        &
!             TEMPERATURES AT GROUND
     &   , T_SOLID(NPD_PROFILE)                                         &
!             TEMPERATURES AT SOLID SURFACE
     &   , T_SEA(NPD_PROFILE)
!             SURFACE TEMPERATURE OVER OPEN SEA
      REAL                                                              &
                !, INTENT(OUT)
     &     PLANCK_SOURCE(NPD_PROFILE, 0: NPD_LAYER)                     &
!             PLANCK FUNCTION ON LEVELS
     &   , DIFF_PLANCK(NPD_PROFILE, NPD_LAYER)                          &
!             DIFFERENCES IN PLANCKIAN FNC
     &   , DIFF_PLANCK_2(NPD_PROFILE, NPD_LAYER)                        &
!             TWICE 2ND DIFFERENCES IN PLANCKIAN FUNCTION
     &   , T(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             TEMPERATURES AT CENTRES OF LAYERS
     &   , PLANCK_GROUND(NPD_PROFILE)                                   &
!             PLANCKIAN FUNCTION AT GROUND
     &   , FLANDG(NPD_PROFILE)
!             GATHERED LAND FRACTION
!
!     ARGUMENTS RELATING TO SEA ICE.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_FRAC_SOL_POINT                                             &
!             NUMBER OF POINTS WITH FRACTIONAL ICE COVER/LAND COVER
     &   , I_FRAC_SOL_POINT(NPD_PROFILE)
!             INDICES OF POINTS WITH FRACTIONAL ICE COVER/LAND COVER
      REAL                                                              &
                !, INTENT(IN)
     &     ICE_FRACTION(NPD_PROFILE)
!             FRACTION OF SEA-ICE IN SEA PORTION OF GRID BOX
!
      REAL                                                              &
                !, INTENT(OUT)
     &     PLANCK_FREEZE_SEA                                            &
!             PLANCK FUNCTION OVER FREEZING SEA
     &   , PLANCK_LEADS_SEA(NPD_PROFILE)
!             PLANCK FUNCTION OVER SEA LEADS
!
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
      REAL                                                              &
     &     T_RATIO(NPD_PROFILE)
!             TEMPERATURE RATIO
!
!     VARIABLES FOR FRACTIONAL ICE COVER.
      INTEGER                                                           &
     &     LG
!             GATHERED LOOP VARIABLE
      REAL                                                              &
     &     T_ICE_G(NPD_PROFILE)                                         &
!             TEMPERATURE OF ICE SURFACE (GATHERED OVER SEA-ICE POINTS)
     &   , PLANCK_GROUND_G(NPD_PROFILE)                                 &
!             PLANCKIAN OF SOLID SURFACE (GATHERED OVER SOLID POINTS)
     &   , SEAFRAC(NPD_PROFILE)
!             FRACTION OF OPEN SEA IN GRID-BOX
!
      REAL ltfs

      ltfs = access_tfs
!
!
!     CALCULATE THE CHANGE IN THE THERMAL SOURCE FUNCTION
!     ACROSS EACH LAYER FOR THE INFRA-RED PART OF THE SPECTRUM.
      DO L=1, N_PROFILE
         T_RATIO(L)=T_LEVEL(L, 0)/T_REF_PLANCK
         PLANCK_SOURCE(L, 0)                                            &
     &      =THERMAL_COEFFICIENT(N_DEG_FIT)
      ENDDO
      DO J=N_DEG_FIT-1, 0, -1
         DO L=1, N_PROFILE
            PLANCK_SOURCE(L, 0)                                         &
     &         =PLANCK_SOURCE(L, 0)                                     &
     &         *T_RATIO(L)+THERMAL_COEFFICIENT(J)
         ENDDO
      ENDDO
      DO I=1, N_LAYER
         DO L=1, N_PROFILE
            T_RATIO(L)=T_LEVEL(L, I)/T_REF_PLANCK
            PLANCK_SOURCE(L, I)                                         &
     &         =THERMAL_COEFFICIENT(N_DEG_FIT)
         ENDDO
         DO J=N_DEG_FIT-1, 0, -1
            DO L=1, N_PROFILE
               PLANCK_SOURCE(L, I)                                      &
     &            =PLANCK_SOURCE(L, I)                                  &
     &            *T_RATIO(L)+THERMAL_COEFFICIENT(J)
            ENDDO
         ENDDO
         DO L=1, N_PROFILE
            DIFF_PLANCK(L, I)=PLANCK_SOURCE(L, I)                       &
     &         -PLANCK_SOURCE(L, I-1)
         ENDDO
      ENDDO
!     CALCULATE THE SECOND DIFFERENCE IF REQUIRED.
      IF (L_IR_SOURCE_QUAD) THEN
         DO I=1, N_LAYER
!           USE THE SECOND DIFFERENCE FOR TEMPORARY STORAGE.
!           OF THE PLANCKIAN AT THE MIDDLE OF THE LAYER.
            DO L=1, N_PROFILE
               T_RATIO(L)=T(L, I)/T_REF_PLANCK
               DIFF_PLANCK_2(L, I)                                      &
     &            =THERMAL_COEFFICIENT(N_DEG_FIT)
            ENDDO
            DO J=N_DEG_FIT-1, 0, -1
               DO L=1, N_PROFILE
                  DIFF_PLANCK_2(L, I)                                   &
     &               =DIFF_PLANCK_2(L, I)                               &
     &               *T_RATIO(L)+THERMAL_COEFFICIENT(J)
               ENDDO
            ENDDO
            DO L=1, N_PROFILE
               DIFF_PLANCK_2(L, I)=2.0E+00*(PLANCK_SOURCE(L, I)         &
     &            +PLANCK_SOURCE(L, I-1)-2.0E+00*DIFF_PLANCK_2(L, I))
            ENDDO
         ENDDO
      ENDIF
!
      IF(L_CTILE)THEN
!     SOURCE AT THE OPEN SEA SURFACE.
!     CALCULATE OVER ALL POINTS EVEN THOUGH IT IS
!     OVERWRITTEN WHERE THERE IS LAND OR SEA-ICE.
      DO L=1, N_PROFILE
         T_RATIO(L)=T_SEA(L)/T_REF_PLANCK
         IF(FLANDG(L) >= 1.0.OR.ICE_FRACTION(L) >= 1.0)                 &
     &      T_RATIO(L)=T_SOLID(L)/T_REF_PLANCK

         PLANCK_GROUND(L)=THERMAL_COEFFICIENT(N_DEG_FIT)
      ENDDO
      DO J=N_DEG_FIT-1, 0, -1
         DO L=1, N_PROFILE
            PLANCK_GROUND(L)=PLANCK_GROUND(L)*T_RATIO(L)                &
     &         +THERMAL_COEFFICIENT(J)
         ENDDO
      ENDDO
!
!
! Initialise to zero:
      DO LG=1,NPD_PROFILE
        PLANCK_LEADS_SEA(LG)=0.0
      ENDDO
!     SET THE SOURCE FUNCTION OVER OPEN SEA LEADS.
!     DETERMINE THE TEMPERATURE OF THE NON-SEA FRACTION.
!     CALCULATE THE SOURCE FUNCTION AT POINTS WITH SOLID SURFACE
      DO L=1, N_FRAC_SOL_POINT
         LG=I_FRAC_SOL_POINT(L)
         PLANCK_LEADS_SEA(LG)=PLANCK_GROUND(LG)
         SEAFRAC(LG)=(1.-FLANDG(LG))*(1.0E+00-ICE_FRACTION(LG))
         T_RATIO(L)=T_SOLID(LG)/T_REF_PLANCK
         PLANCK_GROUND_G(L)=THERMAL_COEFFICIENT(N_DEG_FIT)
      ENDDO
      DO J=N_DEG_FIT-1, 0, -1
         DO L=1, N_FRAC_SOL_POINT
            PLANCK_GROUND_G(L)=PLANCK_GROUND_G(L)*T_RATIO(L)            &
     &         +THERMAL_COEFFICIENT(J)
         ENDDO
      ENDDO
!     DETERMINE THE OVERALL PLANCKIAN FUNCTION OF THE SURFACE.
      DO L=1, N_FRAC_SOL_POINT
         LG=I_FRAC_SOL_POINT(L)
         PLANCK_GROUND(LG)=(1.-SEAFRAC(LG))*PLANCK_GROUND_G(L)          &
     &      +PLANCK_LEADS_SEA(LG)*SEAFRAC(LG)
      ENDDO

      ELSE                      !End of L_CTILE

!     SOURCE AT THE SURFACE.
      DO L=1, N_PROFILE
         T_RATIO(L)=T_GROUND(L)/T_REF_PLANCK
         PLANCK_GROUND(L)=THERMAL_COEFFICIENT(N_DEG_FIT)
      ENDDO
      DO J=N_DEG_FIT-1, 0, -1
         DO L=1, N_PROFILE
            PLANCK_GROUND(L)=PLANCK_GROUND(L)*T_RATIO(L)                &
     &         +THERMAL_COEFFICIENT(J)
         ENDDO
      ENDDO
!
!     WHERE THERE IS FRACTIONAL SEA-ICE THE FORMULATION MUST BE
!     EXTENDED, BUT IT IS CONVENIENT TO CARRY OUT THE ABOVE OPERATION
!     AT ALL POINTS TO AVOID THE USE OF INDIRECT ADDRESSING.
!
!     CALCULATE THE SOURCE FUNCTION OVER OPEN SEA, ADOPTING THE MODEL'S
!     CONVENTION THAT THE TEMPAERTURE THERE IS FIXED.
      T_RATIO(1)=LTFS/T_REF_PLANCK
      PLANCK_FREEZE_SEA=THERMAL_COEFFICIENT(N_DEG_FIT)
      DO J=N_DEG_FIT-1, 0, -1
         PLANCK_FREEZE_SEA=PLANCK_FREEZE_SEA*T_RATIO(1)                 &
     &      +THERMAL_COEFFICIENT(J)
      ENDDO
!
!     DETERMINE THE TEMPERATURE OF THE ICE.
      DO L=1, N_FRAC_SOL_POINT
         LG=I_FRAC_SOL_POINT(L)
         T_ICE_G(L)=(T_GROUND(LG)                                       &
     &      -LTFS*(1.0E+00-ICE_FRACTION(LG)))/ICE_FRACTION(LG)
      ENDDO
!
!     CALCULATE THE SOURCE FUNCTION AT POINTS WITH FRACTIONAL ICE.
      DO L=1, N_FRAC_SOL_POINT
         T_RATIO(L)=T_ICE_G(L)/T_REF_PLANCK
         PLANCK_GROUND_G(L)=THERMAL_COEFFICIENT(N_DEG_FIT)
      ENDDO
      DO J=N_DEG_FIT-1, 0, -1
         DO L=1, N_FRAC_SOL_POINT
            PLANCK_GROUND_G(L)=PLANCK_GROUND_G(L)*T_RATIO(L)            &
     &         +THERMAL_COEFFICIENT(J)
         ENDDO
      ENDDO
!
!     DETERMINE THE OVERALL PLANCKIAN FUNCTION OF THE SURFACE.
      DO L=1, N_FRAC_SOL_POINT
         LG=I_FRAC_SOL_POINT(L)
         PLANCK_GROUND(LG)=ICE_FRACTION(LG)*PLANCK_GROUND_G(L)          &
     &      +PLANCK_FREEZE_SEA*(1.0E+00-ICE_FRACTION(LG))
      ENDDO
!
      ENDIF

!
!
      RETURN
      END SUBROUTINE DIFF_PLANCK_SOURCE
