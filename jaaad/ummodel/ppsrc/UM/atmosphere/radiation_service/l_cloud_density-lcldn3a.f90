


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to determine whether densities are required for clouds.
!
! Method:
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                L_CLOUD_DENSITY set
!                                               as .FALSE. initially
!                                               to provide a default.
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77 with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      FUNCTION L_CLOUD_DENSITY(N_CONDENSED, I_PHASE_CMP, L_CLOUD_CMP    &
     &   , I_CONDENSED_PARAM                                            &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     INCLUDE COMDECKS
! DIMFIX3A defines internal dimensions tied to algorithms for
! two-stream radiation code, mostly for clouds

      ! number of components of clouds
      INTEGER,PARAMETER:: NPD_CLOUD_COMPONENT=4

      ! number of permitted types of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_TYPE=4

      ! number of permitted representations of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_REPRESENTATION=4

      ! number of overlap coefficients for clouds
      INTEGER,PARAMETER:: NPD_OVERLAP_COEFF=18

      ! number of coefficients for two-stream sources
      INTEGER,PARAMETER:: NPD_SOURCE_COEFF=2

      ! number of regions in a layer
      INTEGER,PARAMETER:: NPD_REGION=3

! DIMFIX3A end
! ICLPRM3A defines numbers for ice cloud schemes in two-stream radiation
! code.
!
!   5.5   Feb 2003     Addition of the aggregate parametrization.
!                                                 John Edwards
!   6.2   Jan 2006     Various options for radiation code
!                      3Z added.   (J.-C. Thelen)
!

      ! number of cloud fitting schemes
      INTEGER,PARAMETER:: NPD_ICE_CLOUD_FIT=6

      ! parametrization of slingo and schrecker.
      INTEGER,PARAMETER:: IP_SLINGO_SCHRECKER_ICE=1

      ! unparametrized ice crystal data
       INTEGER,PARAMETER:: IP_ICE_UNPARAMETRIZED=3

      ! sun and shine's parametrization in the visible (version 2)
      INTEGER,PARAMETER:: IP_SUN_SHINE_VN2_VIS=4

      ! sun and shine's parametrization in the ir (version 2)
      INTEGER,PARAMETER:: IP_SUN_SHINE_VN2_IR=5

      ! scheme based on anomalous diffraction theory for ice crystals
      INTEGER,PARAMETER:: IP_ICE_ADT=6

      ! Provisional agregate parametrization.
      INTEGER,PARAMETER:: IP_ICE_AGG_DE=12

! ICLPRM3A end
! PHASE3A defines indices for phases in two-stream radiation code.
      INTEGER,PARAMETER:: IP_PHASE_WATER = 1
      INTEGER,PARAMETER:: IP_PHASE_ICE   = 2
! PHASE3A end
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_CONDENSED                                                  &
!             NUMBER OF TYPES OF CONDENSATE
     &   , I_PHASE_CMP(NPD_CLOUD_COMPONENT)                             &
!             PHASES OF COMPONENTS
     &   , I_CONDENSED_PARAM(NPD_CLOUD_COMPONENT)
!             PARAMETRIZATIONS OF COMPONENTS
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CLOUD_CMP(NPD_CLOUD_COMPONENT)
!             FLAGS FOR ENABLED COMPONENTS
      LOGICAL                                                           &
                !, INTENT(OUT)
     &     L_CLOUD_DENSITY
!             RETURNED FLAG FOR CALCULATING DENSITY
!
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     K
!             LOOP VARIABLE
!
!
      L_CLOUD_DENSITY=.FALSE.
!
!     DENSITIES MUST BE CALCULATED IF SUN & SHINE'S PARAMETRIZATIONS
!     ARE USED.
      DO K=1, N_CONDENSED
         L_CLOUD_DENSITY=L_CLOUD_DENSITY.OR.                            &
     &      (L_CLOUD_CMP(K).AND.(I_PHASE_CMP(K) == IP_PHASE_ICE).AND.   &
     &      ( (I_CONDENSED_PARAM(K) == IP_SUN_SHINE_VN2_VIS).OR.        &
     &        (I_CONDENSED_PARAM(K) == IP_SUN_SHINE_VN2_IR) ) )
      ENDDO
!
!
!
      RETURN
      END FUNCTION L_CLOUD_DENSITY
