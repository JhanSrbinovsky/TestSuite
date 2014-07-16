#if defined(A70_1Z)
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
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE TWO_COEFF_BASIC(IERR                                   &
     &   , N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST                       &
     &   , I_2STREAM                                                    &
     &   , ASYMMETRY, OMEGA                                             &
     &   , SUM, DIFF                                                    &
     &   , ND_PROFILE, ID_LT, ID_LB                                     &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Maximum number of profiles
     &  , ID_LT                                                         &
!           Topmost declared layer
     &  , ID_LB
!           Bottom declared layer
!
!     Include header files.
#include "c_kinds.h"
#include "def_std_io_icf3z.h"
#include "diff_elsasser_ccf3z.h"
#include "two_stream_scheme_pcf3z.h"
#include "error_pcf3z.h"
!
!
!     Dummy arguments.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , I_LAYER_FIRST                                                 &
!           First layer to consider
     &  , I_LAYER_LAST                                                  &
!           Last layer to consider
     &  , I_2STREAM
!           Two stream scheme
!
!     Optical properties of layer:
      REAL  (Real64), INTENT(IN) ::                                     &
     &    ASYMMETRY(ND_PROFILE, ID_LT: ID_LB)                           &
!           Asymmetry factor
     &  , OMEGA(ND_PROFILE, ID_LT: ID_LB)
!           Albedo of single scattering
!
!
!     coefficients in the two-stream equations:
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    SUM(ND_PROFILE, ID_LT: ID_LB)                                 &
!           Sum of alpha_1 and alpha_2
     &  , DIFF(ND_PROFILE, ID_LT: ID_LB)
!           Difference of alpha_1 and alpha_2
!
!     Local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L
!           Loop variable
!
      REAL  (Real64) ::                                                 &
     &    ROOT_3
!           Square root of 3
!
!
      PARAMETER(                                                        &
     &    ROOT_3=1.7320508075688772E+00_Real64                          &
     &  )
!
!
!
!
      IF (I_2STREAM == IP_EDDINGTON) THEN
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1, N_PROFILE
            SUM(L, I)=1.5E+00_Real64*(1.0E+00_Real64                    &
     &        -OMEGA(L, I)*ASYMMETRY(L, I))
            DIFF(L, I)=2.0E+00_Real64*(1.0E+00_Real64-OMEGA(L, I))
          ENDDO
        ENDDO
!
      ELSE IF (I_2STREAM == IP_ELSASSER) THEN
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1, N_PROFILE
            SUM(L, I)=ELSASSER_FACTOR                                   &
     &        -1.5E+00_Real64*OMEGA(L, I)*ASYMMETRY(L, I)
            DIFF(L, I)=ELSASSER_FACTOR*(1.0E+00_Real64-OMEGA(L, I))
          ENDDO
        ENDDO
!
      ELSE IF (I_2STREAM == IP_DISCRETE_ORD) THEN
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1, N_PROFILE
            SUM(L, I)=ROOT_3*(1.0E+00_Real64                            &
     &        -OMEGA(L, I)*ASYMMETRY(L, I))
            DIFF(L, I)=ROOT_3*(1.0E+00_Real64-OMEGA(L, I))
          ENDDO
        ENDDO
!
      ELSE IF (I_2STREAM == IP_PIFM85) THEN
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1, N_PROFILE
            SUM(L, I)=2.0E+00_Real64                                    &
     &        -1.5E+00_Real64*OMEGA(L, I)*ASYMMETRY(L, I)
            DIFF(L, I)=2.0E+00_Real64*(1.0E+00_Real64-OMEGA(L, I))
          ENDDO
        ENDDO
!
      ELSE IF (I_2STREAM == IP_2S_TEST) THEN
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1, N_PROFILE
            SUM(L, I)=1.5E+00_Real64                                    &
     &        -1.5E+00_Real64*OMEGA(L, I)*ASYMMETRY(L, I)
            DIFF(L, I)=1.5E+00_Real64*(1.0E+00_Real64-OMEGA(L, I))
          ENDDO
        ENDDO
!
      ELSE IF (I_2STREAM == IP_HEMI_MEAN) THEN
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1, N_PROFILE
            SUM(L, I)=2.0E+00_Real64                                    &
     &        *(1.0E+00_Real64-OMEGA(L, I)*ASYMMETRY(L, I))
            DIFF(L, I)=2.0E+00_Real64*(1.0E+00_Real64-OMEGA(L, I))
          ENDDO
        ENDDO
!
      ELSE IF (I_2STREAM == IP_PIFM80) THEN
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1, N_PROFILE
            SUM(L, I)=2.0E+00_Real64                                    &
     &        -1.5E+00_Real64*OMEGA(L, I)*ASYMMETRY(L, I)               &
     &        -0.5E+00_Real64*OMEGA(L, I)
            DIFF(L, I)=2.0E+00_Real64*(1.0E+00_Real64-OMEGA(L, I))
          ENDDO
        ENDDO
!
      ELSE IF (I_2STREAM == IP_IFM) THEN
        WRITE(IU_ERR, '(/A)')                                           &
     &    '*** Error: The improved flux mathod has '                    &
     &    //'not been implemented.'
        IERR=I_ERR_FATAL
        RETURN
!
      ELSE IF (I_2STREAM == IP_ZDK_FLUX) THEN
        WRITE(IU_ERR, '(/A)')                                           &
     &    '*** Error: Zdunkowski''s flux method has '                   &
     &    //'not been implemented.'
        IERR=I_ERR_FATAL
        RETURN
!
      ELSE IF (I_2STREAM == IP_KRSCHG_FLUX) THEN
        WRITE(IU_ERR, '(/A)')                                           &
     &    '*** Error: Kerschgen''s flux method has '                    &
     &    //'not been implemented.'
        IERR=I_ERR_FATAL
        RETURN
!
      ELSE IF (I_2STREAM == IP_COAKLEY_CHYLEK_1) THEN
        WRITE(IU_ERR, '(/A)')                                           &
     &    '*** Error: Coakley and Chylek''s first method has '          &
     &    //'not been implemented.'
        IERR=I_ERR_FATAL
        RETURN
!
      ELSE IF (I_2STREAM == IP_COAKLEY_CHYLEK_2) THEN
        WRITE(IU_ERR, '(/A)')                                           &
     &    '*** Error: Coakley and Chylek''s second method has '         &
     &    //'not been implemented.'
        IERR=I_ERR_FATAL
        RETURN
!
      ELSE IF (I_2STREAM == IP_MEADOR_WEAVER) THEN
        WRITE(IU_ERR, '(/A)')                                           &
     &    '*** Error: Meador and Weaver''s method has '                 &
     &    //'not been implemented.'
        IERR=I_ERR_FATAL
        RETURN
!
      ELSE
        WRITE(IU_ERR, '(/A, /A)')                                       &
     &    '*** Error: An unrecognized value has been specified '        &
     &    //'to define the two-stream scheme.'
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
