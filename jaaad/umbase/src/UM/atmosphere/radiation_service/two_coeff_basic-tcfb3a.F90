#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate basic coefficients in two-stream equations.
!
! Method:
!       Depending on the two-stream equations employed, the
!       appropriate coefficients for the fluxes are calculated.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             08-08-96                PIFM85 restored to
!                                               original form. PIFM80
!                                               introduced.
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE TWO_COEFF_BASIC(IERR                                   &
     &   , N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST                       &
     &   , I_2STREAM                                                    &
     &   , ASYMMETRY, OMEGA                                             &
     &   , SUM, DIFF                                                    &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
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
#include "stdio3a.h"
#include "twostr3a.h"
#include "elsass3a.h"
#include "error3a.h"
!
!
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , I_LAYER_FIRST                                                &
!             FIRST LAYER TO CONSIDER
     &   , I_LAYER_LAST                                                 &
!             LAST LAYER TO CONSIDER
     &   , I_2STREAM
!             TWO STREAM SCHEME
!
!     OPTICAL PROPERTIES OF LAYER:
      REAL                                                              &
                !, INTENT(IN)
     &     ASYMMETRY(NPD_PROFILE, NPD_LAYER)                            &
!             ASYMMETRY FACTOR
     &   , OMEGA(NPD_PROFILE, NPD_LAYER)
!             ALBEDO OF SINGLE SCATTERING
!
!
!     COEFFICIENTS IN THE TWO-STREAM EQUATIONS:
      REAL                                                              &
                !, INTENT(OUT)
     &     SUM(NPD_PROFILE, NPD_LAYER)                                  &
!             SUM OF ALPHA_1 AND ALPHA_2
     &   , DIFF(NPD_PROFILE, NPD_LAYER)
!             DIFFERENCE OF ALPHA_1 AND ALPHA_2
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
      REAL                                                              &
     &     ROOT_3
!             SQUARE ROOT OF 3
!
!
      PARAMETER(                                                        &
     &     ROOT_3=1.7320508075688772E+00                                &
     &   )
!
!
!
!
      IF (I_2STREAM == IP_EDDINGTON) THEN
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               SUM(L, I)=1.5E+00*(1.0E+00                               &
     &            -OMEGA(L, I)*ASYMMETRY(L, I))
               DIFF(L, I)=2.0E+00*(1.0E+00-OMEGA(L, I))
            ENDDO
         ENDDO
!
      ELSE IF (I_2STREAM == IP_ELSASSER) THEN
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               SUM(L, I)=ELSASSER_FACTOR                                &
     &            -1.5E+00*OMEGA(L, I)*ASYMMETRY(L, I)
               DIFF(L, I)=ELSASSER_FACTOR*(1.0E+00-OMEGA(L, I))
            ENDDO
         ENDDO
!
      ELSE IF (I_2STREAM == IP_DISCRETE_ORD) THEN
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               SUM(L, I)=ROOT_3*(1.0E+00                                &
     &            -OMEGA(L, I)*ASYMMETRY(L, I))
               DIFF(L, I)=ROOT_3*(1.0E+00-OMEGA(L, I))
            ENDDO
         ENDDO
!
      ELSE IF (I_2STREAM == IP_PIFM85) THEN
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               SUM(L, I)=2.0E+00                                        &
     &            -1.5E+00*OMEGA(L, I)*ASYMMETRY(L, I)
               DIFF(L, I)=2.0E+00*(1.0E+00-OMEGA(L, I))
            ENDDO
         ENDDO
!
      ELSE IF (I_2STREAM == IP_2S_TEST) THEN
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               SUM(L, I)=1.5E+00                                        &
     &            -1.5E+00*OMEGA(L, I)*ASYMMETRY(L, I)
               DIFF(L, I)=1.5E+00*(1.0E+00-OMEGA(L, I))
            ENDDO
         ENDDO
!
      ELSE IF (I_2STREAM == IP_HEMI_MEAN) THEN
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               SUM(L, I)=2.0E+00                                        &
     &            *(1.0E+00-OMEGA(L, I)*ASYMMETRY(L, I))
               DIFF(L, I)=2.0E+00*(1.0E+00-OMEGA(L, I))
            ENDDO
         ENDDO
!
      ELSE IF (I_2STREAM == IP_PIFM80) THEN
         DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
               SUM(L, I)=2.0E+00                                        &
     &            -1.5E+00*OMEGA(L, I)*ASYMMETRY(L, I)                  &
     &            -0.5E+00*OMEGA(L, I)
               DIFF(L, I)=2.0E+00*(1.0E+00-OMEGA(L, I))
            ENDDO
         ENDDO
!
      ELSE IF (I_2STREAM == IP_IFM) THEN
         WRITE(IU_ERR, '(/A)')                                          &
     &      '*** ERROR: THE IMPROVED FLUX METHOD HAS '                  &
     &      //'NOT BEEN IMPLEMENTED.'
         IERR=I_ERR_FATAL
         RETURN
!
      ELSE IF (I_2STREAM == IP_ZDK_FLUX) THEN
         WRITE(IU_ERR, '(/A)')                                          &
     &      '*** ERROR: ZDUNKOWSKI''S FLUX METHOD HAS '                 &
     &      //'NOT BEEN IMPLEMENTED.'
         IERR=I_ERR_FATAL
         RETURN
!
      ELSE IF (I_2STREAM == IP_KRSCHG_FLUX) THEN
         WRITE(IU_ERR, '(/A)')                                          &
     &      '*** ERROR: KERSCHGEN''S FLUX METHOD HAS '                  &
     &      //'NOT BEEN IMPLEMENTED.'
         IERR=I_ERR_FATAL
         RETURN
!
      ELSE IF (I_2STREAM == IP_COAKLEY_CHYLEK_1) THEN
         WRITE(IU_ERR, '(/A)')                                          &
     &      '*** ERROR: COAKLEY-CHYLEK''S FIRST METHOD HAS '            &
     &      //'NOT BEEN IMPLEMENTED.'
         IERR=I_ERR_FATAL
         RETURN
!
      ELSE IF (I_2STREAM == IP_COAKLEY_CHYLEK_2) THEN
         WRITE(IU_ERR, '(/A)')                                          &
     &      '*** ERROR: COAKLEY-CHYLEK''S SECOND METHOD HAS '           &
     &      //'NOT BEEN IMPLEMENTED.'
         IERR=I_ERR_FATAL
         RETURN
!
      ELSE IF (I_2STREAM == IP_MEADOR_WEAVER) THEN
         WRITE(IU_ERR, '(/A)')                                          &
     &      '*** ERROR: MEADOR & WEAVER''S METHOD HAS '                 &
     &      //'NOT BEEN IMPLEMENTED.'
         IERR=I_ERR_FATAL
         RETURN
!
      ELSE
         WRITE(IU_ERR, '(/A)')                                          &
     &      '*** ERROR: AN ILLEGAL PARAMETER HAS BEEN SUPPLIED '        &
     &      //'TO DEFINE THE 2-STREAM SCHEME.'
         IERR=I_ERR_FATAL
         RETURN
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE TWO_COEFF_BASIC
#endif
#endif
