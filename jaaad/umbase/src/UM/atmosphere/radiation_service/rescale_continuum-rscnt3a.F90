#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to apply a path-length scaling to the continuum.
!
! Method:
!       The scaling function is calculated. This is multpiled by a
!       suitable "amount" of continuum incorporating a broadening
!       density.
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
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE RESCALE_CONTINUUM(N_PROFILE, N_LAYER, I_CONTINUUM      &
     &   , P, T, L_LAYER, I_TOP                                         &
     &   , DENSITY, MOLAR_DENSITY_WATER, MOLAR_DENSITY_FRN              &
     &   , WATER_FRAC                                                   &
     &   , AMOUNT_CONTINUUM                                             &
     &   , I_FNC                                                        &
     &   , P_REFERENCE, T_REFERENCE, SCALE_PARAMETER                    &
     &   , NPD_PROFILE, NPD_LAYER, NPD_SCALE_FNC                        &
     &   , NPD_SCALE_VARIABLE                                           &
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
     &   , NPD_LAYER                                                    &
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_SCALE_FNC                                                &
!             NUMBER OF SCALING FUNCTIONS
     &   , NPD_SCALE_VARIABLE
!             MAX. NUMBER OF SCALING VARIABLES
!
!     INCLUDE COMDECKS
#include "phycn03a.h"
#include "cntuum3a.h"
#include "sclfnc3a.h"
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , I_CONTINUUM                                                  &
!             CONTINUUM TYPE
     &   , I_FNC                                                        &
!             SCALING FUNCTION
     &   , I_TOP
!             TOP INDEX OF ARRAYS
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_LAYER
!             DATA ARE SUPPLIED IN LAYERS
      REAL                                                              &
                !, INTENT(IN)
     &     WATER_FRAC(NPD_PROFILE, 0: NPD_LAYER)                        &
!             MASS FRACTION OF WATER
     &   , P(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             PRESSURE
     &   , T(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             TEMPERATURE
     &   , DENSITY(NPD_PROFILE, 0: NPD_LAYER)                           &
!             OVERALL DENSITY
     &   , MOLAR_DENSITY_WATER(NPD_PROFILE, 0: NPD_LAYER)               &
!             MOLAR DENSITY OF WATER VAPOUR
     &   , MOLAR_DENSITY_FRN(NPD_PROFILE, 0: NPD_LAYER)                 &
!             MOLAR DENSITY OF FOREIGN SPECIES
     &   , P_REFERENCE                                                  &
!             REFERENCE PRESSURE
     &   , T_REFERENCE                                                  &
!             REFERENCE PRESSURE
     &   , SCALE_PARAMETER(NPD_SCALE_VARIABLE)
!             SCALING PARAMTERS
      REAL                                                              &
                !, INTENT(OUT)
     &     AMOUNT_CONTINUUM(NPD_PROFILE, 0: NPD_LAYER)
!             AMOUNT OF CONTINUUM
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     L                                                            &
!             LOOP VARIABLE
     &   , I
!             LOOP VARIABLE
      REAL    PWK(N_PROFILE,N_LAYER-I_TOP+1)  ! Workspace
      REAL    TWK(N_PROFILE,N_LAYER-I_TOP+1)  ! Workspace
#if defined(VECTLIB)
      REAL    SP1(N_PROFILE,N_LAYER-I_TOP+1)  ! Workspace
      REAL    SP2(N_PROFILE,N_LAYER-I_TOP+1)  ! Workspace
      INTEGER n_input             ! No. of inputs for rtor_v function
#endif
!
!
#if defined(VECTLIB)
         do I=1, N_LAYER-I_TOP+1
            do L=1, N_PROFILE
              sp1(L,I)=SCALE_PARAMETER(1)
              sp2(L,I)=SCALE_PARAMETER(2)
            end do
         end do
         n_input=(N_LAYER-I_TOP+1)*N_PROFILE
#endif
         DO I=   1, N_LAYER-I_TOP+1
            DO L=1, N_PROFILE
              PWK(L,I)=P(L, I_TOP+I-1)/P_REFERENCE
            END DO
         END DO
#if defined(VECTLIB)
         call rtor_v(n_input,pwk,sp1,pwk)
#else
         DO I=   1, N_LAYER-I_TOP+1
            DO L=1, N_PROFILE
              PWK(L,I)=PWK(L,I)**SCALE_PARAMETER(1)
            ENDDO
         ENDDO
#endif
!
      IF (I_FNC == IP_SCALE_POWER_LAW) THEN
!
         DO I=   1, N_LAYER-I_TOP+1
            DO L=1, N_PROFILE
              TWK(L,I)=T(L, I_TOP+I-1)/T_REFERENCE
            END DO
         END DO
#if defined(VECTLIB)
         call rtor_v(n_input,twk,sp2,twk)
#else
         DO I=   1, N_LAYER-I_TOP+1
            DO L=1, N_PROFILE
              TWK(L,I)=TWK(L,I)**SCALE_PARAMETER(2)
            ENDDO
         ENDDO
#endif
!
         DO I=I_TOP, N_LAYER
            DO L=1, N_PROFILE
               AMOUNT_CONTINUUM(L, I)                                   &
     &                       =PWK(L,I-I_TOP+1)*TWK(L,I-I_TOP+1)
            ENDDO
         ENDDO
      ELSE IF(I_FNC == IP_SCALE_POWER_QUAD) THEN
         DO I=I_TOP, N_LAYER
            DO L=1, N_PROFILE
               AMOUNT_CONTINUUM(L, I)                                   &
     &            =PWK(L,I-I_TOP+1)                                     &
     &            *(1.0E+00+SCALE_PARAMETER(2)*(T(L, I)                 &
     &            /T_REFERENCE-1.0E+00)                                 &
     &            +SCALE_PARAMETER(3)*(T(L, I)                          &
     &            /T_REFERENCE-1.0E+00)**2)
            ENDDO
         ENDDO
      ENDIF
!
      IF (L_LAYER) THEN
         IF (I_CONTINUUM == IP_SELF_CONTINUUM) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  AMOUNT_CONTINUUM(L, I)=AMOUNT_CONTINUUM(L, I)         &
     &               *MOLAR_DENSITY_WATER(L, I)*WATER_FRAC(L, I)
               ENDDO
            ENDDO
         ELSE IF (I_CONTINUUM == IP_FRN_CONTINUUM) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  AMOUNT_CONTINUUM(L, I)=AMOUNT_CONTINUUM(L, I)         &
     &               *MOLAR_DENSITY_FRN(L, I)*WATER_FRAC(L, I)
               ENDDO
            ENDDO
         ELSE IF (I_CONTINUUM == IP_N2_CONTINUUM) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  AMOUNT_CONTINUUM(L, I)=AMOUNT_CONTINUUM(L, I)         &
     &               *N2_MASS_FRAC*DENSITY(L, I)
               ENDDO
            ENDDO
         ENDIF
      ELSE
!        IF VALUES ARE GIVEN ON LEVELS WE NOW INTERPOLATE TO AVERAGES
!        ACROSS THE LAYER.
         IF (I_CONTINUUM == IP_SELF_CONTINUUM) THEN
            DO I=N_LAYER, 1, -1
               DO L=1, N_PROFILE
                  AMOUNT_CONTINUUM(L, I)=0.5E+00                        &
     &               *(AMOUNT_CONTINUUM(L, I)                           &
     &               *MOLAR_DENSITY_WATER(L, I)*WATER_FRAC(L, I)        &
     &               +AMOUNT_CONTINUUM(L, I-1)                          &
     &               *MOLAR_DENSITY_WATER(L, I-1)*WATER_FRAC(L, I-1))
               ENDDO
            ENDDO
         ELSE IF (I_CONTINUUM == IP_FRN_CONTINUUM) THEN
            DO I=N_LAYER, 1, -1
               DO L=1, N_PROFILE
                  AMOUNT_CONTINUUM(L, I)=0.5E+00                        &
     &               *(AMOUNT_CONTINUUM(L, I)                           &
     &               *MOLAR_DENSITY_FRN(L, I)*WATER_FRAC(L, I)          &
     &               +AMOUNT_CONTINUUM(L, I-1)                          &
     &               *MOLAR_DENSITY_FRN(L, I-1)*WATER_FRAC(L, I-1))
               ENDDO
            ENDDO
         ELSE IF (I_CONTINUUM == IP_N2_CONTINUUM) THEN
            DO I=N_LAYER, 1, -1
               DO L=1, N_PROFILE
                  AMOUNT_CONTINUUM(L, I)=0.5E+00                        &
     &               *(AMOUNT_CONTINUUM(L, I)                           &
     &               *N2_MASS_FRAC*DENSITY(L, I)                        &
     &               +AMOUNT_CONTINUUM(L, I-1)                          &
     &               *N2_MASS_FRAC*DENSITY(L, I-1))
               ENDDO
            ENDDO
         ENDIF
      ENDIF
!
!
      RETURN
      END SUBROUTINE RESCALE_CONTINUUM
#endif
#endif
