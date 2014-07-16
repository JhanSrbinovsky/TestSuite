


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to set number of cloudy parameters.
!
! Method:
!       Straightforward
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.5             18-05-98                Code for new parametr-
!                                               ization of droplets
!                                               included.
!                                               (J. M. Edwards)
!       5.5             24-02-03                Code for aggregate
!                                               parametrization of
!                                               ice crystals included.
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77 with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      FUNCTION SET_N_CLOUD_PARAMETER(I_SCHEME, I_COMPONENT              &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     INCLUDE COMDECKS
! CLDCMP3A sets components of clouds for two-stream radiation code.

      ! stratiform water droplets
      INTEGER,PARAMETER:: IP_CLCMP_ST_WATER=1

      ! stratiform ice crystals
      INTEGER,PARAMETER:: IP_CLCMP_ST_ICE=2

      ! convective water droplets
      INTEGER,PARAMETER:: IP_CLCMP_CNV_WATER=3

      ! convective ice crystals
      INTEGER,PARAMETER:: IP_CLCMP_CNV_ICE=4

! CLDCMP3A end
! WCLPRM3A defines numbers for water cloud schemes in two-stream
! radiation code.
      ! number of cloud fitting schemes
      INTEGER,PARAMETER:: NPD_CLOUD_FIT=3

      ! parametrization of slingo-schrecker
      INTEGER,PARAMETER:: IP_SLINGO_SCHRECKER=1

      ! parametrization of ackerman & stephens
      INTEGER,PARAMETER:: IP_ACKERMAN_STEPHENS=2

      ! unparametrized droplet data
      INTEGER,PARAMETER:: IP_DROP_UNPARAMETRIZED=3

      ! pade approximation of the second order (third order for the
      ! extinction)
      INTEGER,PARAMETER:: IP_DROP_PADE_2=5
! WCLPRM3A end
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
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_SCHEME                                                     &
!             PARAMETRIZATION SCHEME
     &   , I_COMPONENT
!             COMPONENT IN CLOUD
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     SET_N_CLOUD_PARAMETER
!             RETURNED NUMBER OF COEFFICIENTS IN PARAMETRIZATION
!
!
!
      IF ( (I_COMPONENT == IP_CLCMP_ST_WATER).OR.                       &
     &     (I_COMPONENT == IP_CLCMP_CNV_WATER) ) THEN
!
         IF (I_SCHEME == IP_SLINGO_SCHRECKER) THEN
           SET_N_CLOUD_PARAMETER =6
         ELSE IF (I_SCHEME == IP_ACKERMAN_STEPHENS) THEN
            SET_N_CLOUD_PARAMETER=9
         ELSE IF (I_SCHEME == IP_DROP_PADE_2) THEN
            SET_N_CLOUD_PARAMETER=16
         ENDIF
!
      ELSE IF ( (I_COMPONENT == IP_CLCMP_ST_ICE).OR.                    &
     &          (I_COMPONENT == IP_CLCMP_CNV_ICE) ) THEN
!
         IF (I_SCHEME == IP_SLINGO_SCHRECKER_ICE) THEN
            SET_N_CLOUD_PARAMETER=6
         ELSE IF (I_SCHEME == IP_ICE_ADT) THEN
            SET_N_CLOUD_PARAMETER=30
         ELSE IF (I_SCHEME == IP_SUN_SHINE_VN2_VIS) THEN
            SET_N_CLOUD_PARAMETER=6
         ELSE IF (I_SCHEME == IP_SUN_SHINE_VN2_IR) THEN
            SET_N_CLOUD_PARAMETER=0
         ELSE IF (I_SCHEME == IP_ICE_AGG_DE) THEN
            SET_N_CLOUD_PARAMETER=14
         ENDIF
!
      ENDIF
!
!
!
      RETURN
      END FUNCTION SET_N_CLOUD_PARAMETER
