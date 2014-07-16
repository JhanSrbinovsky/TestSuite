


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate fluxes using Gaussian quadrature.
!
! Method:
!       Fluxes are calculated by using Gaussian quadrature for
!       the angular integration. This is not a full implementation
!       of Gaussian quadrature for multiple scattering, but is
!       intended only for non-scattering calculations in the
!       infra-red. In this case, the fluxes can be calculated as
!       a weighted sum of two-stream fluxes where the diffusivity
!       factors for the two-stream approximations are determined
!       from the Gaussian points.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE GAUSS_ANGLE(N_PROFILE, N_LAYER, L_NET, N_AUGMENT       &
     &   , N_ORDER_GAUSS                                                &
     &   , TAU                                                          &
     &   , FLUX_INC_DOWN                                                &
     &   , DIFF_PLANCK, SOURCE_GROUND, ALBEDO_SURFACE_DIFF              &
     &   , FLUX_DIFFUSE                                                 &
     &   , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                              &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
      IMPLICIT NONE
!
!
!     SIZES OF DUMMY ARRAYS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     INCLUDE COMDECKS.
! SPCRG3A defines flags for different portions of the spectrum in
! two-stream radiation code.
      INTEGER,PARAMETER:: IP_SOLAR=1
      INTEGER,PARAMETER:: IP_INFRA_RED=2
! SPCRG3A end
! GSSWTP3A defines arrays of gaussian points and weights for two-stream
! radiation code.

      ! max. order of gaussian quadrature
      INTEGER,PARAMETER:: NPD_GAUSS_ORD=10

      ! points of gaussian integration
      REAL :: GAUSS_POINT(NPD_GAUSS_ORD, NPD_GAUSS_ORD)

      ! weights of gaussian integration
      REAL :: GAUSS_WEIGHT(NPD_GAUSS_ORD, NPD_GAUSS_ORD)

      ! NB values defined in GSSWTD3A

! GSSWTP3A end
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
!
!     DUMMY VARIABLES.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , N_AUGMENT                                                    &
!             SIZE OF ARRAY TO INCREMENT
     &   , N_ORDER_GAUSS
!             ORDER OF GAUSSIAN INTEGRATION
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_NET                                                        &
!             NET FLUXES REQUIRED
     &   , L_IR_SOURCE_QUAD
!             USE QUADRATIC SOURCE TERM
      REAL                                                              &
                !, INTENT(IN)
     &     TAU(NPD_PROFILE, NPD_LAYER)                                  &
!             OPTICAL DEPTH
     &   , ALBEDO_SURFACE_DIFF(NPD_PROFILE)                             &
!             DIFFUSE ALBEDO
     &   , FLUX_INC_DOWN(NPD_PROFILE)                                   &
!             INCIDENT TOTAL FLUX
     &   , DIFF_PLANCK(NPD_PROFILE, NPD_LAYER)                          &
!             DIFFERENCE IN PI*PLANCK FUNCTION
     &   , SOURCE_GROUND(NPD_PROFILE)                                   &
!             GROUND SOURCE FUNCTION
     &   , DIFF_PLANCK_2(NPD_PROFILE, NPD_LAYER)
!             2x2ND DIFFERENCES OF PLANCKIAN
      REAL                                                              &
                !, INTENT(OUT)
     &     FLUX_DIFFUSE(NPD_PROFILE, 2*NPD_LAYER+2)
!             DIFFUSE FLUXES
!
!     LOCAL VARIABALES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , K
!             LOOP VARIABLE
      REAL                                                              &
     &     FLUX_STREAM(NPD_PROFILE, 2*NPD_LAYER+2)                      &
!             FLUX IN STREAM
     &   , FLUX_NULL(NPD_PROFILE, 2*NPD_LAYER+2)                        &
!             ARRAY OF NULL FLUXES
     &   , SECANT_RAY                                                   &
!             SECANT OF ANGLE WITH VERTICAL
     &   , DIFF_PLANCK_RAD(NPD_PROFILE, NPD_LAYER)                      &
!             DIFFERENCE IN PI*PLANCK FUNCTION
     &   , DIFF_PLANCK_RAD_2(NPD_PROFILE, NPD_LAYER)                    &
!             2x2ND DIFFERENCES OF PLANCKIAN
     &   , SOURCE_GROUND_RAD(NPD_PROFILE)                               &
!             GROUND SOURCE FUNCTION
     &   , RADIANCE_INC(NPD_PROFILE)                                    &
!             INCIDNET RADIANCE
     &   , WEIGHT_STREAM
!             WEIGHTING FOR STREAM
!
!     SET THE GAUSSIAN WEIGHTS FOR INTEGRATION.
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE SETTING GAUSSIAN POINTS AND WEIGHTS: ONLY EVEN
!     ORDERS ARE USED AT PRESENT.
!
!     FIRST ORDER:
      DATA GAUSS_POINT(1,1)/0.00000000000000E+00/
      DATA GAUSS_WEIGHT(1,1)/2.00000000000000E+00/
!
!     SECOND ORDER:
      DATA GAUSS_POINT(1,2)/-5.77350269189626E-01/
      DATA GAUSS_POINT(2,2)/5.77350269189626E-01/
      DATA GAUSS_WEIGHT(1,2)/1.00000000000000E+00/
      DATA GAUSS_WEIGHT(2,2)/1.00000000000000E+00/
!
!     THIRD ORDER:
      DATA GAUSS_POINT(1,3)/-7.74596669241484E-01/
      DATA GAUSS_POINT(2,3)/0.00000000000000E+00/
      DATA GAUSS_POINT(3,3)/7.74596669241484E-01/
      DATA GAUSS_WEIGHT(1,3)/5.55555555555556E-01/
      DATA GAUSS_WEIGHT(2,3)/8.88888888888889E-01/
      DATA GAUSS_WEIGHT(3,3)/5.55555555555556E-01/
!
!     FOURTH ORDER:
      DATA GAUSS_POINT(1,4)/-8.61136311594053E-01/
      DATA GAUSS_POINT(2,4)/-3.39981043584856E-01/
      DATA GAUSS_POINT(3,4)/3.39981043584856E-01/
      DATA GAUSS_POINT(4,4)/8.61136311594053E-01/
      DATA GAUSS_WEIGHT(1,4)/3.47854845137454E-01/
      DATA GAUSS_WEIGHT(2,4)/6.52145154862546E-01/
      DATA GAUSS_WEIGHT(3,4)/6.52145154862546E-01/
      DATA GAUSS_WEIGHT(4,4)/3.47854845137454E-01/
!
!     FIFTH ORDER:
      DATA GAUSS_POINT(1,5)/-9.06179845938664E-01/
      DATA GAUSS_POINT(2,5)/-5.38469310105683E-01/
      DATA GAUSS_POINT(3,5)/0.00000000000000E+00/
      DATA GAUSS_POINT(4,5)/5.38469310105683E-01/
      DATA GAUSS_POINT(5,5)/9.06179845938664E-01/
      DATA GAUSS_WEIGHT(1,5)/2.36926885056189E-01/
      DATA GAUSS_WEIGHT(2,5)/4.78628670499366E-01/
      DATA GAUSS_WEIGHT(3,5)/4.67913934572691E-01/
      DATA GAUSS_WEIGHT(4,5)/4.78628670499366E-01/
      DATA GAUSS_WEIGHT(5,5)/2.36926885056189E-01/
!
!     SIXTH ORDER:
      DATA GAUSS_POINT(1,6)/-9.32469514203152E-01/
      DATA GAUSS_POINT(2,6)/-6.61209386466265E-01/
      DATA GAUSS_POINT(3,6)/-2.38619186083197E-01/
      DATA GAUSS_POINT(4,6)/2.38619186083197E-01/
      DATA GAUSS_POINT(5,6)/6.61209386466265E-01/
      DATA GAUSS_POINT(6,6)/9.32469514203152E-01/
      DATA GAUSS_WEIGHT(1,6)/1.71324492379170E-01/
      DATA GAUSS_WEIGHT(2,6)/3.60761573048139E-01/
      DATA GAUSS_WEIGHT(3,6)/4.67913934572691E-01/
      DATA GAUSS_WEIGHT(4,6)/4.67913934572691E-01/
      DATA GAUSS_WEIGHT(5,6)/3.60761573048139E-01/
      DATA GAUSS_WEIGHT(6,6)/1.71324492379170E-01/
!
!     SEVENTH ORDER:
      DATA GAUSS_POINT(1,7)/-9.49107912342759E-01/
      DATA GAUSS_POINT(2,7)/-7.41531185599394E-01/
      DATA GAUSS_POINT(3,7)/-4.05845151377397E-01/
      DATA GAUSS_POINT(4,7)/0.00000000000000E+00/
      DATA GAUSS_POINT(5,7)/4.05845151377397E-01/
      DATA GAUSS_POINT(6,7)/7.41531185599394E-01/
      DATA GAUSS_POINT(7,7)/9.49107912342759E-01/
      DATA GAUSS_WEIGHT(1,7)/1.29484966168870E-01/
      DATA GAUSS_WEIGHT(2,7)/2.79705391489277E-01/
      DATA GAUSS_WEIGHT(3,7)/3.81830050505119E-01/
      DATA GAUSS_WEIGHT(4,7)/4.17959183673469E-01/
      DATA GAUSS_WEIGHT(5,7)/3.81830050505119E-01/
      DATA GAUSS_WEIGHT(6,7)/2.79705391489277E-01/
      DATA GAUSS_WEIGHT(7,7)/1.29484966168870E-01/
!
!     EIGHTH ORDER:
      DATA GAUSS_POINT(1,8)/-9.60289856497536E-01/
      DATA GAUSS_POINT(2,8)/-7.96666477413627E-01/
      DATA GAUSS_POINT(3,8)/-5.25532409916329E-01/
      DATA GAUSS_POINT(4,8)/-1.83434642495650E-01/
      DATA GAUSS_POINT(5,8)/1.83434642495650E-01/
      DATA GAUSS_POINT(6,8)/5.25532409916329E-01/
      DATA GAUSS_POINT(7,8)/7.96666477413627E-01/
      DATA GAUSS_POINT(8,8)/9.60289856497536E-01/
      DATA GAUSS_WEIGHT(1,8)/1.01228536290376E-01/
      DATA GAUSS_WEIGHT(2,8)/2.22381034453374E-01/
      DATA GAUSS_WEIGHT(3,8)/3.13706645877887E-01/
      DATA GAUSS_WEIGHT(4,8)/3.62683783378362E-01/
      DATA GAUSS_WEIGHT(5,8)/3.62683783378362E-01/
      DATA GAUSS_WEIGHT(6,8)/3.13706645877887E-01/
      DATA GAUSS_WEIGHT(7,8)/2.22381034453374E-01/
      DATA GAUSS_WEIGHT(8,8)/1.01228536290376E-01/
!
!     NINTH ORDER:
      DATA GAUSS_POINT(1,9)/-9.68160239507626E-01/
      DATA GAUSS_POINT(2,9)/-8.36031107326636E-01/
      DATA GAUSS_POINT(3,9)/-6.13371432700590E-01/
      DATA GAUSS_POINT(4,9)/-3.24253423403809E-01/
      DATA GAUSS_POINT(5,9)/0.00000000000000E+00/
      DATA GAUSS_POINT(6,9)/3.24253423403809E-01/
      DATA GAUSS_POINT(7,9)/6.13371432700590E-01/
      DATA GAUSS_POINT(8,9)/8.36031107326636E-01/
      DATA GAUSS_POINT(9,9)/9.68160239507626E-01/
      DATA GAUSS_WEIGHT(1,9)/8.1274388361574E-02/
      DATA GAUSS_WEIGHT(2,9)/1.80648160694857E-01/
      DATA GAUSS_WEIGHT(3,9)/2.60610696402935E-01/
      DATA GAUSS_WEIGHT(4,9)/3.12347077040003E-01/
      DATA GAUSS_WEIGHT(5,9)/3.30239355001260E-01/
      DATA GAUSS_WEIGHT(6,9)/3.12347077040003E-01/
      DATA GAUSS_WEIGHT(7,9)/2.60610696402935E-01/
      DATA GAUSS_WEIGHT(8,9)/1.80648160694857E-01/
      DATA GAUSS_WEIGHT(9,9)/8.1274388361574E-02/
!
!     TENTH ORDER:
      DATA GAUSS_POINT(1,10)/-9.73906528517172E-01/
      DATA GAUSS_POINT(2,10)/-8.65063366688985E-01/
      DATA GAUSS_POINT(3,10)/-6.79409568299024E-01/
      DATA GAUSS_POINT(4,10)/-4.33395394129427E-01/
      DATA GAUSS_POINT(5,10)/-1.48874338981631E-01/
      DATA GAUSS_POINT(6,10)/1.48874338981631E-01/
      DATA GAUSS_POINT(7,10)/4.33395394129427E-01/
      DATA GAUSS_POINT(8,10)/6.79409568299024E-01/
      DATA GAUSS_POINT(9,10)/8.65063366688985E-01/
      DATA GAUSS_POINT(10,10)/9.73906528517172E-01/
      DATA GAUSS_WEIGHT(1,10)/6.6671344308688E-02/
      DATA GAUSS_WEIGHT(2,10)/1.49451349150581E-01/
      DATA GAUSS_WEIGHT(3,10)/2.19086362515982E-01/
      DATA GAUSS_WEIGHT(4,10)/2.69266719309996E-01/
      DATA GAUSS_WEIGHT(5,10)/2.95524224714753E-01/
      DATA GAUSS_WEIGHT(6,10)/2.95524224714753E-01/
      DATA GAUSS_WEIGHT(7,10)/2.69266719309996E-01/
      DATA GAUSS_WEIGHT(8,10)/2.19086362515982E-01/
      DATA GAUSS_WEIGHT(9,10)/1.49451349150581E-01/
      DATA GAUSS_WEIGHT(10,10)/6.6671344308688E-02/
!
!     ------------------------------------------------------------------
!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     MONOCHROMATIC_IR_RADIANCE, AUGMENT_FLUX
!
!
!
!     SET THE SOURCE FUNCTION.
      DO L=1, N_PROFILE
         SOURCE_GROUND_RAD(L)=SOURCE_GROUND(L)/PI
         RADIANCE_INC(L)=FLUX_INC_DOWN(L)/PI
      ENDDO
      DO I=1, N_LAYER
         DO L=1, N_PROFILE
            DIFF_PLANCK_RAD(L, I)=DIFF_PLANCK(L, I)/PI
         ENDDO
      ENDDO
      DO I=1, 2*N_LAYER+2
         DO L=1, N_PROFILE
            FLUX_DIFFUSE(L, I)=0.0
         ENDDO
      ENDDO
      IF (L_IR_SOURCE_QUAD) THEN
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               DIFF_PLANCK_RAD_2(L, I)=DIFF_PLANCK_2(L, I)/PI
            ENDDO
         ENDDO
      ENDIF
!
!     CALCULATE THE FLUXES WITH A NUMBER OF DIFFUSIVITY FACTORS
!     AND SUM THE RESULTS.
      DO K=1, N_ORDER_GAUSS
         SECANT_RAY=2.0E+00/(GAUSS_POINT(K, N_ORDER_GAUSS)+1.0E+00)
!
!        CALCULATE THE RADIANCE AT THIS ANGLE.
! DEPENDS ON: monochromatic_ir_radiance
            CALL MONOCHROMATIC_IR_RADIANCE(N_PROFILE, N_LAYER           &
     &         , L_NET                                                  &
     &         , TAU                                                    &
     &         , RADIANCE_INC                                           &
     &         , DIFF_PLANCK_RAD, SOURCE_GROUND_RAD, ALBEDO_SURFACE_DIFF&
     &         , SECANT_RAY                                             &
     &         , FLUX_STREAM                                            &
     &         , NPD_PROFILE, NPD_LAYER                                 &
     &         )
!
!        AUGMENT THE FLUX BY THE AMOUNT IN THIS STREAM.
         WEIGHT_STREAM=5.0E-01*PI*GAUSS_WEIGHT(K, N_ORDER_GAUSS)        &
     &      *(GAUSS_POINT(K, N_ORDER_GAUSS)+1.0E+00)
! DEPENDS ON: augment_flux
         CALL AUGMENT_FLUX(N_PROFILE, N_LAYER, N_AUGMENT                &
     &      , IP_INFRA_RED, .FALSE.                                     &
     &      , WEIGHT_STREAM                                             &
     &      , FLUX_NULL, FLUX_DIFFUSE                                   &
     &      , FLUX_NULL, FLUX_STREAM                                    &
     &      , FLUX_NULL, FLUX_NULL                                      &
     &      , FLUX_NULL, FLUX_NULL                                      &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
!
      ENDDO
!
!
      RETURN
      END SUBROUTINE GAUSS_ANGLE
