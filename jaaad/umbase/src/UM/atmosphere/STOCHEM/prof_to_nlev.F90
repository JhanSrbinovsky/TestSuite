#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE PROF_TO_NLEV(M_Storm,HT_IN_ETA,N_PROFILE)
! ----------------------------------------------------------------------
! Purpose:
!  Converts profiles of percentage mass NO production per km to
!  fractional mass NO production per vertical STOCHEM level
!
! Method:
!
! Original Programmer: Michael Sanderson
!
! Current code owner: Michael Sanderson
!
! History:
! Date        Version     Comment
! -------     -------     -------
!  1/8/02     1.0         Original    Michael Sanderson
!  5/8/02     1.1         ND version  Colin Johnson
! ----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
! ----------------------------------------------------------------------

! M_Storm   - Storm type (1-3); see below
! HT_IN_ETA - Profile heights in eta coordinates scaled to actual cloud
!             top height
! N_PROFILE - NO production profile (as a fraction of total) on STOCHEM
!             vertical levels

      INTEGER, INTENT(IN) :: M_Storm
      REAL, DIMENSION(0:NLKMS), INTENT(IN) :: HT_IN_ETA
      REAL, DIMENSION(NLEV), INTENT(OUT) :: N_PROFILE

      INTEGER :: K, L
      REAL :: ETALO, ETAHI, HTLO, HTHI, FRAC

! Define NO production profiles. Units are % of mass of N per level.
! 3 columns for Midlat continental, Tropical Marine and Tropical
! Continental storms. From Pickering et al. JGR 1998 Vol.103 No.D23,
! pp.31203-31216.  These are the standard profiles.  They have been
! scaled to particular cloud top heights in SUB LIGHTNOX, and the
! revised altitudes (in eta coords) are passed here in HT_IN_ETA().

!        Mid Con T.Mar. T.Con.
      REAL, DIMENSION(NLPROFS,NLKMS) :: N_PROD_PROF = RESHAPE((/        &
     &  (/ 20.1,  5.8,  8.2 /),                                         &
                                  !  0- 1 km
     &  (/  2.3,  2.9,  1.9 /),                                         &
                                  !  1- 2 km
     &  (/  0.8,  2.6,  2.1 /),                                         &
                                  !  2- 3 km
     &  (/  1.5,  2.4,  1.6 /),                                         &
                                  !  3- 4 km
     &  (/  3.4,  2.2,  1.1 /),                                         &
                                  !  4- 5 km
     &  (/  5.3,  2.1,  1.6 /),                                         &
                                  !  5- 6 km
     &  (/  3.6,  2.3,  3.0 /),                                         &
                                  !  6- 7 km
     &  (/  3.8,  6.1,  5.8 /),                                         &
                                  !  7- 8 km
     &  (/  5.4, 16.5,  7.6 /),                                         &
                                  !  8- 9 km
     &  (/  6.6, 14.1,  9.6 /),                                         &
                                  !  9-10 km
     &  (/  8.3, 13.7, 10.5 /),                                         &
                                  ! 10-11 km
     &  (/  9.6, 12.8, 12.3 /),                                         &
                                  ! 11-12 km
     &  (/ 12.8, 12.5, 11.8 /),                                         &
                                  ! 12-13 km
     &  (/ 10.0,  2.8, 12.5 /),                                         &
                                  ! 13-14 km
     &  (/  6.2,  0.9,  8.1 /),                                         &
                                  ! 14-15 km
     &  (/  0.3,  0.3,  2.3 /)                                          &
                                  ! 15-16 km
     &  /), SHAPE = (/ NLPROFS,NLKMS /))

! Convert the profiles from 16 eta levels to NLEV eta levels
      N_PROFILE = 0.0
      DO L = 0, NLEV-1
        ETALO = Eta_Stochem(L)     ! Bottom of STOCHEM vertical layer
        ETAHI = Eta_Stochem(L+1)   ! Top of STOCHEM vertical layer
        DO K = 1, NLKMS
          HTLO = HT_IN_ETA(K-1)  ! Bottom of scaled profile layer
          HTHI = HT_IN_ETA(K)    ! Top of scaled profile layer
          IF (HTHI < ETALO) CYCLE
          IF (HTLO > ETAHI) EXIT
          FRAC = 0.0
          IF (HTLO >= ETALO .AND. HTHI <= ETAHI) THEN
            FRAC = 1.0
          ELSE IF (HTLO <= ETALO.AND.HTHI > ETAHI) THEN
            FRAC = (ETALO-ETAHI) / (HTLO-HTHI)
          ELSE IF (HTLO < ETALO.AND.HTHI > ETALO.AND.HTHI < ETAHI) THEN
            FRAC = (ETALO-HTHI) / (HTLO-HTHI)
          ELSE IF (HTLO > ETALO.AND.HTLO < ETAHI.AND.HTHI > ETAHI) THEN
            FRAC = (HTLO-ETAHI) / (HTLO-HTHI)
          END IF
          N_PROFILE(L+1) = N_PROFILE(L+1)+FRAC*N_PROD_PROF(M_Storm,K)
        END DO
      END DO

! Convert profiles from percentages to fractions
      N_PROFILE = N_PROFILE / 100.0

      END SUBROUTINE PROF_TO_NLEV
#endif
