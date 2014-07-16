#if defined(A71_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Subroutine to calculate sea-salt aerosol number concentration.
!
! Purpose:
!   Sea-salt aerosol number concentration in 2 modes is calculated
!   from windspeed, height above surface and land fraction.
!
! Method:
!   10m windspeed is calculated assuming a neutral profile and a
!   roughness length appropriate for water. This is then used to
!   calculate sea-salt aerosol number concentration in 2 modes
!   following O'Dowd et al. (1999). An exponential decay with
!   height of these concentrations is then assumed, and they are
!   scaled depending on the fraction of open water in the gridbox.
!
! Current Owner of Code: A. Jones
!
! History:
!       Version         Date                    Comment
!       5.2             09-10-00                Original Code
!                                               (A. Jones)
!
!       5.4             25-04-02                Modified to take
!                                               account of coastal
!                                               tiling.
!                                               (A. Jones)
!
!       6.2             03-02-06                Moved to A71_1A
!                                               (P.Selwood)
!  6.2  21/11/05  make powers integers                   A.Malcolm
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_SEASALT(                                           &
     &     U_WIND                                                       &
     &   , V_WIND                                                       &
     &   , HEIGHT                                                       &
     &   , LAND_FRACT                                                   &
     &   , SEA_ICE                                                      &
     &   , ROW_LENGTH                                                   &
     &   , NUM_ROWS                                                     &
     &   , NUM_LEVELS                                                   &
     &   , N_LEVELS_BL                                                  &
     &   , SEA_SALT_FILM                                                &
     &   , SEA_SALT_JET                                                 &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!
      INTEGER                                                           &
     &     ROW_LENGTH                                                   &
!             Number of points in a row
     &   , NUM_ROWS                                                     &
!             Number of rows
     &   , NUM_LEVELS                                                   &
!             Number of levels
     &   , N_LEVELS_BL
!             Number of boundary-layer levels
!
      REAL                                                              &
     &     U_WIND(ROW_LENGTH, NUM_ROWS)                                 &
!             u-component of level 1 windspeed (ms-1)
     &   , V_WIND(ROW_LENGTH, NUM_ROWS)                                 &
!             v-component of level 1 windspeed (ms-1)
     &   , LAND_FRACT(ROW_LENGTH, NUM_ROWS)                             &
!             Land fraction
     &   , SEA_ICE(ROW_LENGTH, NUM_ROWS)                                &
!             Sea-ice fraction
     &   , HEIGHT(ROW_LENGTH, NUM_ROWS, NUM_LEVELS)
!             Height of layer centres above sea surface (m)
!
      REAL                                                              &
     &     SEA_SALT_FILM(ROW_LENGTH, NUM_ROWS, NUM_LEVELS)              &
!             Number concentration of sea-salt aerosol derived
!             from film-mode droplets (m-3)
     &   , SEA_SALT_JET(ROW_LENGTH, NUM_ROWS, NUM_LEVELS)
!             Number concentration of sea-salt aerosol derived
!             from jet-mode droplets (m-3)
!
!     Local variables:
!
      INTEGER                                                           &
     &     I                                                            &
              !
     &   , J                                                            &
              ! Loop counters
     &   , K  !
!
      REAL                                                              &
     &     Z0_SEA                                                       &
!             Roughness length over sea (m)
     &   , WINDSPEED_10M                                                &
!             10m windspeed (ms-1)
     &   , FILM_SCALE_HEIGHT                                            &
!             Scale height for film-mode aerosol (m)
     &   , JET_SCALE_HEIGHT
!             Scale height for jet-mode aerosol (m)
!
      PARAMETER(                                                        &
     &     Z0_SEA=2.5E-04                                               &
     &   , FILM_SCALE_HEIGHT=900.0                                      &
     &   , JET_SCALE_HEIGHT=900.0                                       &
     &   )
!
!
      DO K = 1, NUM_LEVELS
        DO J = 1, NUM_ROWS
          DO I = 1, ROW_LENGTH

            IF (LAND_FRACT(I, J) <  1.0 .AND. SEA_ICE(I, J) <  1.0      &
     &                                  .AND. K <= N_LEVELS_BL) THEN
               WINDSPEED_10M = (SQRT((U_WIND(I, J)**2)                  &
     &                                         +(V_WIND(I, J)**2)))     &
     &                          *(ALOG(10.0/Z0_SEA))                    &
     &                          /(ALOG(HEIGHT(I, J, 1)/Z0_SEA))
               SEA_SALT_FILM(I, J, K)                                   &
     &                          =10.0**((0.095*WINDSPEED_10M)+6.283)
               SEA_SALT_JET(I, J, K)                                    &
     &                          =10.0**((0.0422*WINDSPEED_10M)+5.7122)

               IF (WINDSPEED_10M  <   2.0) THEN
                 SEA_SALT_FILM(I, J, K)                                 &
     &                       =3.856E06*(1.0-EXP(-0.736*WINDSPEED_10M))
                 SEA_SALT_JET(I, J, K)                                  &
     &                       =0.671E06*(1.0-EXP(-1.351*WINDSPEED_10M))
               ENDIF

               IF (WINDSPEED_10M  >   17.5) THEN
                 SEA_SALT_FILM(I, J, K)=150.0E06*(1.0-(97.874           &
     &                                    *EXP(-0.313*WINDSPEED_10M)))
                 SEA_SALT_JET(I, J, K)=3.6E06*(1.0-(103.926             &
     &                                    *EXP(-0.353*WINDSPEED_10M)))
               ENDIF

               SEA_SALT_FILM(I, J, K)=SEA_SALT_FILM(I, J, K)            &
     &                  *EXP(10.0/FILM_SCALE_HEIGHT)                    &
     &                  *EXP(-HEIGHT(I, J, K)/FILM_SCALE_HEIGHT)        &
     &                  *(1.0-LAND_FRACT(I, J))*(1.0-SEA_ICE(I, J))

               SEA_SALT_JET(I, J, K)=SEA_SALT_JET(I, J, K)              &
     &                  *EXP(10.0/JET_SCALE_HEIGHT)                     &
     &                  *EXP(-HEIGHT(I, J, K)/JET_SCALE_HEIGHT)         &
     &                  *(1.0-LAND_FRACT(I, J))*(1.0-SEA_ICE(I, J))

            ELSE

               SEA_SALT_FILM(I, J, K)=0.0
               SEA_SALT_JET(I, J, K)=0.0

            ENDIF

          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE SET_SEASALT
#endif
