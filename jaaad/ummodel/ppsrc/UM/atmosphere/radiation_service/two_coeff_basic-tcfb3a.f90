


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
! STDIO3A defines unit numbers for standard i/o in two-stream radiation
! code.
      INTEGER,PARAMETER:: IU_STDIN=5
      INTEGER,PARAMETER:: IU_STDOUT=6
      INTEGER,PARAMETER:: IU_ERR=6
! STDIO3A end
! TWOSTR3A defines the defining numbers for the two-stream schemes
      ! Eddington approximation
      INTEGER,PARAMETER::  IP_EDDINGTON=2

      ! discrete ordinate method
      INTEGER,PARAMETER:: IP_DISCRETE_ORD=4

      ! improved flux method
      INTEGER,PARAMETER:: IP_IFM=5

      ! practical improved flux method (version of Zdunkowski et al.
      ! 1985)
      INTEGER,PARAMETER:: IP_PIFM85=6

      ! Zdunkowski's flux method
      INTEGER,PARAMETER:: IP_ZDK_FLUX=7

      ! Lerschgen's flux method
      INTEGER,PARAMETER:: IP_KRSCHG_FLUX=8

      ! Coakley & Chylek's 1st method
      INTEGER,PARAMETER:: IP_COAKLEY_CHYLEK_1=9

      ! Coakley & Chylek's 2nd method
      INTEGER,PARAMETER:: IP_COAKLEY_CHYLEK_2=10

      ! Meador & Weaver's method
      INTEGER,PARAMETER:: IP_MEADOR_WEAVER=11

      ! Elsasser's diffusivity scheme
      INTEGER,PARAMETER:: IP_ELSASSER=12

      ! user's defined test approximation.
      INTEGER,PARAMETER:: IP_2S_TEST=14

      ! hemispheric mean approximation.
      INTEGER,PARAMETER:: IP_HEMI_MEAN=15

      ! practical improved flux method (version of Zdunkowski et al.
      ! 1980)
      INTEGER,PARAMETER:: IP_PIFM80=16
! TWOSTR3A end
! ELSASS3A defines Diffusivity for Elsasser's scheme in two-stream
! radiation code.

      REAL,PARAMETER:: ELSASSER_FACTOR=1.66E+00

! ELSASS3A end
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET ERROR FLAGS IN THE RADIATION CODE.
!
      INTEGER                                                           &
     &     I_NORMAL                                                     &
!             ERROR FREE CONDITION
     &   , I_ERR_FATAL                                                  &
!             FATAL ERROR: IMMEDIATE RETURN
     &   , I_ABORT_CALCULATION                                          &
!             CALCULATION ABORTED
     &   , I_MISSING_DATA                                               &
!             MISSING DATA ERROR: CONDITIONAL
     &   , I_ERR_IO                                                     &
!             I/O ERROR
     &   , I_ERR_RANGE                                                  &
!             INTERPOLATION RANGE ERROR
     &   , I_ERR_EXIST
!             EXISTENCE ERROR
!
      PARAMETER(                                                        &
     &     I_NORMAL=0                                                   &
     &   , I_ERR_FATAL=1                                                &
     &   , I_ABORT_CALCULATION=2                                        &
     &   , I_MISSING_DATA=3                                             &
     &   , I_ERR_IO=4                                                   &
     &   , I_ERR_RANGE=5                                                &
     &   , I_ERR_EXIST=6                                                &
     &   )
!
!     ------------------------------------------------------------------
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
