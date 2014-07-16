#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate optical properties of water clouds.
!
! Method:
!       If the optical properties come from an observational
!       distribution a separate subroutine is called. Otherwise
!       appropriate mean quantities in the layer are calculated
!       as the parametrization requires and these values are
!       substituted into the parametrization to give the optical
!       properties.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             Oct. 96     T3E migration: HF functions
!                                   replaced by T3E vec_lib function
!                                   rtor_v      (S.J.Swarbrick)
!       4.5             18-05-98                New parametrization
!                                               of the optical
!                                               properties of cloud
!                                               droplets added.
!                                               (J. M. Edwards)
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!  6.0  21/08/03  NEC SX-6 optimisation - add vectorisation
!                 !CDIR NODEP directives.  R Barnes & J-C Rioual.
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE OPT_PROP_WATER_CLOUD(IERR                              &
     &   , N_PROFILE, N_LAYER, N_CLOUD_TOP                              &
     &   , N_CLOUD_PROFILE, I_CLOUD_PROFILE                             &
     &   , L_RESCALE, L_LAYER, L_CLOUD_LAYER                            &
     &   , I_PARAMETRIZATION_DROP, CLOUD_PARAMETER                      &
     &   , LIQ_WATER_MASS_FRAC, RADIUS_EFFECT                           &
     &   , K_EXT_TOT_CLOUD, K_EXT_SCAT_CLOUD                            &
     &   , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD                       &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   , NPD_CLOUD_PARAMETER                                          &
     &   )
!
!
      IMPLICIT NONE
!
!
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER                                                    &
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_CLOUD_PARAMETER
!             MAXIMUM NUMBER OF CLOUD PARAMETERS
!
!     INCLUDE COMDECKS
#include "stdio3a.h"
#include "wclprm3a.h"
#include "error3a.h"
!
!     DUMMY VARIABLES.
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , N_CLOUD_TOP                                                  &
!             TOPMOST CLOUDY LAYER
     &   , I_PARAMETRIZATION_DROP                                       &
!             TREATMENT OF DROPLETS
     &   , N_CLOUD_PROFILE(NPD_LAYER)                                   &
!             NUMBER OF CLOUDY PROFILES
     &   , I_CLOUD_PROFILE(NPD_PROFILE, NPD_LAYER)
!             PROFILES CONTAINING CLOUDS
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_LAYER                                                      &
!             VARIABLES GIVEN IN LAYERS
     &   , L_CLOUD_LAYER                                                &
!             CLOUD VARIABLES GIVEN IN LAYERS
     &   , L_RESCALE
!             FLAG FOR DELTA-RESCALING
      REAL                                                              &
                !, INTENT(IN)
     &     CLOUD_PARAMETER(NPD_CLOUD_PARAMETER)                         &
!             CLOUD PARAMETERS
     &   , LIQ_WATER_MASS_FRAC(NPD_PROFILE, 0: NPD_LAYER)               &
!             LIQUID WATER CONTENT
     &   , RADIUS_EFFECT(NPD_PROFILE, 0: NPD_LAYER)
!             EFFECTIVE RADIUS
      REAL                                                              &
                !, INTENT(OUT)
     &     K_EXT_SCAT_CLOUD(NPD_PROFILE, NPD_LAYER)                     &
!             SCATTERING EXTINCTION
     &   , K_EXT_TOT_CLOUD(NPD_PROFILE, NPD_LAYER)                      &
!             TOTAL EXTINCTION
     &   , ASYMMETRY_CLOUD(NPD_PROFILE, NPD_LAYER)                      &
!             CLOUDY ASYMMETRIES
     &   , FORWARD_SCATTER_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY FORWARD SCATTERING
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     L                                                            &
!             LOOP VARIABLE
     &   , LL                                                           &
!             LOOP VARIABLE
     &   , I                                                            &
!             LOOP VARIABLE
     &   , j !loop variable
      REAL                                                              &
     &     ASYMMETRY_PROCESS(NPD_PROFILE)                               &
!             ASYMMETRY FACTOR FOR CURRENT PROC.
     &   , RADIUS_AVE(3)                                                &
!             AVERAGE EFFECTIVE RADIUS IN LAYER
     &   , LIQ_MASS_FRAC_AVE                                            &
!             AVERAGE LIQUID WATER MASS FRACTION
     &   , cp(3),cpp(3)
!             workspace array
!     HALF-PRECISION FUNCTIONS FOR THE UNIFIED MODEL.
!
!
!
!
      IF (I_PARAMETRIZATION_DROP == IP_SLINGO_SCHRECKER) THEN
!
!        WE CALCULATE AVERAGE PROPERTIES FOR THE LAYER AND PUT THESE
!        INTO THE PARAMETRIZATION, RATHER THAN CALCULATING THE
!        PARAMETRIZATION AT EACH LEVEL: USUALLY THIS IS MORE ACCURATE.
!        IT ALSO FITS MORE NATURALLY WITH CASES WHERE DATA ARE GIVEN
!        IN LAYERS.
!
!
         DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
            DO LL=1, N_CLOUD_PROFILE(I)
               L=I_CLOUD_PROFILE(LL, I)
               LIQ_MASS_FRAC_AVE=LIQ_WATER_MASS_FRAC(L, I)
               RADIUS_AVE(1)=RADIUS_EFFECT(L, I)
               K_EXT_TOT_CLOUD(L, I)                                    &
     &            =LIQ_MASS_FRAC_AVE*(CLOUD_PARAMETER(1)                &
     &            +CLOUD_PARAMETER(2)/RADIUS_AVE(1))
               K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)             &
     &            *(1.0E+00-CLOUD_PARAMETER(3)                          &
     &            -CLOUD_PARAMETER(4)*RADIUS_AVE(1))
               ASYMMETRY_PROCESS(L)=                                    &
     &            CLOUD_PARAMETER(5)+CLOUD_PARAMETER(6)                 &
     &            *RADIUS_AVE(1)
               ASYMMETRY_CLOUD(L, I)=                                   &
     &            K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
            ENDDO
            IF (L_RESCALE) THEN
!CDIR NODEP
               DO LL=1, N_CLOUD_PROFILE(I)
                  L=I_CLOUD_PROFILE(LL, I)
                  FORWARD_SCATTER_CLOUD(L, I)                           &
     &               =K_EXT_SCAT_CLOUD(L, I)                            &
     &               *ASYMMETRY_PROCESS(L)**2
               ENDDO
            ENDIF
         ENDDO
!
      ELSE IF (I_PARAMETRIZATION_DROP == IP_ACKERMAN_STEPHENS) THEN
!
!     Set up CP array for use in rtor_v function
         CP(1)=CLOUD_PARAMETER(3)
         CP(2)=CLOUD_PARAMETER(6)
         CP(3)=CLOUD_PARAMETER(9)
!
         DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
            DO LL=1, N_CLOUD_PROFILE(I)
               L=I_CLOUD_PROFILE(LL, I)
               LIQ_MASS_FRAC_AVE=LIQ_WATER_MASS_FRAC(L, I)
               RADIUS_AVE(1)=RADIUS_EFFECT(L, I)
               RADIUS_AVE(2)=RADIUS_EFFECT(L, I)
               RADIUS_AVE(3)=RADIUS_EFFECT(L, I)
#if defined(VECTLIB)
               call rtor_v(3,radius_ave,cp,cpp)
#else
               do j=1,3
                 cpp(j)=radius_ave(j)**cp(j)
               end do
#endif
               K_EXT_TOT_CLOUD(L, I)=LIQ_MASS_FRAC_AVE                  &
     &            *( CLOUD_PARAMETER(1)+CLOUD_PARAMETER(2)              &
     &            *  CPP(1) )
               K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)             &
     &            *(1.0E+00-CLOUD_PARAMETER(4)-CLOUD_PARAMETER(5)       &
     &            *  CPP(2) )
               ASYMMETRY_PROCESS(L)=                                    &
     &            CLOUD_PARAMETER(7)+CLOUD_PARAMETER(8)                 &
     &            *  CPP(3)
               ASYMMETRY_CLOUD(L, I)=                                   &
     &            K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
            ENDDO
            IF (L_RESCALE) THEN
!CDIR NODEP
               DO LL=1, N_CLOUD_PROFILE(I)
                  L=I_CLOUD_PROFILE(LL, I)
                  FORWARD_SCATTER_CLOUD(L, I)                           &
     &               =K_EXT_SCAT_CLOUD(L, I)                            &
     &               *ASYMMETRY_PROCESS(L)**2
               ENDDO
            ENDIF
         ENDDO
!
      ELSE IF (I_PARAMETRIZATION_DROP == IP_DROP_PADE_2) THEN
!
         DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
            DO LL=1, N_CLOUD_PROFILE(I)
               L=I_CLOUD_PROFILE(LL, I)
               LIQ_MASS_FRAC_AVE=LIQ_WATER_MASS_FRAC(L, I)
               RADIUS_AVE(1)=RADIUS_EFFECT(L, I)
               K_EXT_TOT_CLOUD(L, I)=LIQ_MASS_FRAC_AVE                  &
     &            *(CLOUD_PARAMETER(1)+RADIUS_AVE(1)                    &
     &            *(CLOUD_PARAMETER(2)+RADIUS_AVE(1)                    &
     &            *CLOUD_PARAMETER(3)))                                 &
     &            /(1.0E+00+RADIUS_AVE(1)                               &
     &            *(CLOUD_PARAMETER(4)+RADIUS_AVE(1)                    &
     &            *(CLOUD_PARAMETER(5)+RADIUS_AVE(1)                    &
     &            *CLOUD_PARAMETER(6))))
               K_EXT_SCAT_CLOUD(L, I)=K_EXT_TOT_CLOUD(L, I)*(1.0E+00    &
     &            -(CLOUD_PARAMETER(7)+RADIUS_AVE(1)                    &
     &            *(CLOUD_PARAMETER(8)+RADIUS_AVE(1)                    &
     &            *CLOUD_PARAMETER(9)))                                 &
     &            /(1.0E+00+RADIUS_AVE(1)                               &
     &            *(CLOUD_PARAMETER(10)+RADIUS_AVE(1)                   &
     &            *CLOUD_PARAMETER(11))))
               ASYMMETRY_PROCESS(L)                                     &
     &            =(CLOUD_PARAMETER(12)+RADIUS_AVE(1)                   &
     &            *(CLOUD_PARAMETER(13)+RADIUS_AVE(1)                   &
     &            *CLOUD_PARAMETER(14)))                                &
     &            /(1.0E+00+RADIUS_AVE(1)                               &
     &            *(CLOUD_PARAMETER(15)+RADIUS_AVE(1)                   &
     &            *CLOUD_PARAMETER(16)))
               ASYMMETRY_CLOUD(L, I)=                                   &
     &            K_EXT_SCAT_CLOUD(L, I)*ASYMMETRY_PROCESS(L)
            ENDDO
            IF (L_RESCALE) THEN
!CDIR NODEP
               DO LL=1, N_CLOUD_PROFILE(I)
                  L=I_CLOUD_PROFILE(LL, I)
                  FORWARD_SCATTER_CLOUD(L, I)                           &
     &               =K_EXT_SCAT_CLOUD(L, I)                            &
     &               *ASYMMETRY_PROCESS(L)**2
               ENDDO
            ENDIF
         ENDDO
!
      ELSE
         WRITE(IU_ERR, '(/A)') '*** ERROR: AN INVALID PARAMETRIZATION ' &
     &      //'OF CLOUD DROPLETS HAS BEEN USED.'
         IERR=I_ERR_FATAL
         RETURN
      ENDIF
!
!
      RETURN
      END SUBROUTINE OPT_PROP_WATER_CLOUD
#endif
#endif
