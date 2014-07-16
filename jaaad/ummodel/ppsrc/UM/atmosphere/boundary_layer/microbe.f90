
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!***********************************************************************
! Calculates the soil respiration based on a simplified version of the
! model of Raich et al. (1991).
!
!    Model            Modification history
!   version  date
!    5.2   15/11/00   New Deck         M. Best
!  6.2  01/03/06  adapt for use with RothC soil carbon model.
!                                                        C.D. Jones
!
!    Programming standard:
!
!***********************************************************************
      SUBROUTINE MICROBE (LAND_PTS,DIM_CS1,DIM_CS2,L_TRIFFID,L_Q10,CS,  &
     &                    STH_SOIL,V_SAT,V_WILT,TSOIL,RESP_S,VEG_FRAC)

      IMPLICIT NONE

      INTEGER                                                           &
     & LAND_PTS                                                         &
                                  ! IN Number of land points to be
!                                 !    processed.
     &,DIM_CS1, DIM_CS2           ! IN soil carbon dimensions
      LOGICAL                                                           &
     &        L_TRIFFID                                                 &
                                  ! TRUE if using TRIFFID
     &,       L_Q10               ! TRUE if using Q10 for soil resp

      REAL                                                              &
     & CS(LAND_PTS,DIM_CS1)                                             &
                                  ! IN Soil carbon (kg C/m2).
     &,VEG_FRAC(DIM_CS2)                                                &
                                  ! IN vegetated fraction.
     &,STH_SOIL(LAND_PTS)                                               &
                                  ! IN Top layer soil moisture as a
!                                 !    fraction of saturation (m3/m3).
     &,V_SAT(LAND_PTS)                                                  &
                                  ! IN Volumetric soil moisture
!                                 !    concentration at saturation
!                                 !    (m3 H2O/m3 soil).
     &,V_WILT(LAND_PTS)                                                 &
                                  ! IN Volumetric soil moisture
!                                 !    concentration below which
!                                 !    stomata close (m3 H2O/m3 soil).
                                  !    as a fraction of saturation.
     &,TSOIL(LAND_PTS)                                                  &
                                  ! IN Soil temperature (K).
     &,RESP_S(LAND_PTS,DIM_CS1)                                         &
                                  ! OUT Soil respiration (kg C/m2/s).
     &,FSTH(LAND_PTS)                                                   &
                                  ! WORK Factors describing the
     &,FPRF(LAND_PTS)                                                   &
                                  !      influence of soil moisture,
     &,FTEMP(LAND_PTS)                                                  &
                                  !      vegetation cover and
     &,STH_RESP_MIN                                                     &
                                  ! WORK soil moist giving min. resp
!                                 !      soil temperature respectively
!                                 !      on the soil respiration.
     &,STH_OPT                                                          &
                                  ! WORK Fractional soil moisture at
!                                 !      which respiration is maximum.
     &,STH_WILT                   ! WORK Wilting soil moisture as a
!                                 !      fraction of saturation.
      INTEGER                                                           &
     & L                          ! Loop counter

!-----------------------------------------------------------------------
! Local parameters
!-----------------------------------------------------------------------
      REAL                                                              &
     & KAPS(4)                                                          &
                                  ! Specific soil respiration rate
     &,MIN_FACTOR                                                       &
                                  ! FACTOR to scale WILT to get RESP_MIN
!                                 ! at 25 deg ! and optimum soil
!                                 ! moisture (/s).
     &,Q10                        ! Q10 factor for soil respiration.
      PARAMETER (Q10 = 2.0)

      MIN_FACTOR = 1.7
      KAPS(1) = 3.22E-7
      KAPS(2) = 9.65E-9
      KAPS(3) = 2.12E-8
      KAPS(4) = 6.43E-10

      IF (.NOT. L_TRIFFID) THEN
        L_Q10 = .TRUE.
        KAPS(1) = 0.5E-8
        MIN_FACTOR = 1.0
      ENDIF


! FSTH
      DO L=1, LAND_PTS
        FSTH(L)=0.0
      ENDDO
      DO L=1,LAND_PTS
        IF (V_SAT(L)  >   0.0) THEN

          STH_WILT = V_WILT(L) / V_SAT(L)
          STH_OPT = 0.5 * (1 + STH_WILT)
          STH_RESP_MIN = STH_WILT * MIN_FACTOR

          FSTH(L) = 0.2
          IF (STH_SOIL(L)  >   STH_RESP_MIN .AND.                       &
     &            STH_SOIL(L)  <=  STH_OPT) THEN
            FSTH(L) = 0.2 + 0.8 * ((STH_SOIL(L) - STH_RESP_MIN)         &
     &                          /  (STH_OPT - STH_RESP_MIN))
          ELSEIF (STH_SOIL(L)  >   STH_OPT) THEN
            FSTH(L) = 1 - 0.8 * (STH_SOIL(L) - STH_OPT)
          ENDIF
        ENDIF
      ENDDO

! FTEMP
      IF (L_Q10) THEN
! use original HadCM3LC Q10 formula
        DO L=1,LAND_PTS
          FTEMP(L) = Q10 ** (0.1 * (TSOIL(L) - 298.15))
        ENDDO
      ELSE
! use RothC temperature formula (using TSOIL for now...)
        DO L=1,LAND_PTS
          FTEMP(L)=0.0
          IF (TSOIL(L)  >   263.15)                                     &
     &        FTEMP(L) = 47.9 / (1.0+EXP(106.0/(TSOIL(L)-254.85)))
        ENDDO
      ENDIF

! FPRF - plant retainment factor
!          =0.6 for fully vegetated
!          =1.0 for bare soil
!        only set for RothC runs, (i.e. L_TRIFFID=TRUE)
!
      IF (L_TRIFFID) THEN
        DO L=1,LAND_PTS
          FPRF(L) = 0.6 + 0.4*(1-VEG_FRAC(L))
        ENDDO
      ENDIF

! set 1-D or 4-D soil resp depending on whether using RothC or not
      IF (L_TRIFFID) THEN
        DO L=1,LAND_PTS
          RESP_S(L,1) = KAPS(1)*CS(L,1) * FSTH(L)*FTEMP(L)*FPRF(L)
          RESP_S(L,2) = KAPS(2)*CS(L,2) * FSTH(L)*FTEMP(L)*FPRF(L)
          RESP_S(L,3) = KAPS(3)*CS(L,3) * FSTH(L)*FTEMP(L)*FPRF(L)
          RESP_S(L,4) = KAPS(4)*CS(L,4) * FSTH(L)*FTEMP(L)*FPRF(L)
        ENDDO
      ELSE
        DO L=1,LAND_PTS
          RESP_S(L,1) = KAPS(1) * CS(L,1) * FSTH(L) * FTEMP(L)
        ENDDO
      ENDIF

      RETURN

      END SUBROUTINE MICROBE
