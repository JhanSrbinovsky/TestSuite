
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE SICE_HTF-----------------------------------------------
!!!
!!!  Purpose: Updates sea-ice surface layer temperature.
!!!
!!!  Model            Modification history
!!! version  date
!!!
!!!  5.2   15/11/00   New Deck         M. Best
!!!  5.3   25/04/01  Add coastal tiling.   Nic Gedney
!!!  5.5   07/02/03  Added ice catagories. J. Ridley
!!!
!!!  Note: At present the formulation is so simple as to make this
!!!        routine fairly trivial; but in future the formulation may
!!!        be revised so as to make a separate routine more obviously
!!!        worthwhile.
!!!
!!!  Programming standard: Unified Model Documentation Paper No.4
!!!                        version no.2, dated 18/1/90.
!!!
!!!  System component covered: P241
!!!
!!!  Documentation: ??
!!!

! Arguments:---------------------------------------------------------
      SUBROUTINE SICE_HTF (                                             &
     & ROW_LENGTH,ROWS,FLANDG,SIMLT,NICE,                               &
     & DI_NCAT,ICE_FRACTION,ICE_FRACT_NCAT,SURF_HT_FLUX_SICE_NCAT,      &
     & TSTAR_SICE,TIMESTEP,                                             &
     & TI,SICE_MLT_HTF,SEA_ICE_HTF,L_SICE_HEATFLUX,                     &
     & LTIMER)


      USE auscom_cpl_data_mod,                                          &
     &    Only : auscom_salinity, ocn_sss, access_tfs



      IMPLICIT NONE

      LOGICAL LTIMER

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                            ! IN Number of X points?
     &,ROWS                                                             &
                            ! IN Number of Y points?
     &,NICE                 ! IN Number of ice catagories

      LOGICAL                                                           &
     & SIMLT                      ! IN STASH flag for sea-ice melting
!                                 !    ht flux.

      REAL                                                              &
     & FLANDG(ROW_LENGTH,ROWS)                                          &
                                    ! IN Land fraction.
     &,DI_NCAT(ROW_LENGTH,ROWS,NICE)                                    &
                                     ! IN "Equivalent thickness" of
!                                   !    sea-ice (m).
     &,ICE_FRACTION(ROW_LENGTH,ROWS)                                    &
                                    ! IN Fraction of gridbox covered by
!                                   !    sea-ice.
     &,ICE_FRACT_NCAT(ROW_LENGTH,ROWS,NICE)                             &
                                           ! IN Fraction of ice in
!                                 ! gridbox covered by each ice catagory
     &,SURF_HT_FLUX_SICE_NCAT(ROW_LENGTH,ROWS,NICE)                     &
!                                   ! IN Net downward heat flux at
!                                   !    sea-ice surface W/m2
     &,TSTAR_SICE(ROW_LENGTH,ROWS)                                      &
                                    ! INOUT Sea-ice surface
!                                   !    temperature (K).
     &,TIMESTEP                     ! IN Timestep (s).

      REAL                                                              &
     & TI(ROW_LENGTH,ROWS,NICE)                                         &
                                ! INOUT Sea-ice surface layer
!                                   !       temperature(K) Set to TSTAR
!                                   !       for unfrozen sea, missing
!                                   !       data for land.
     &,SICE_MLT_HTF(ROW_LENGTH,ROWS,NICE)                               &
                                         ! INOUT Heat flux due to melt
!                                   !       of sea-ice (W/m2).
     &,SEA_ICE_HTF(ROW_LENGTH,ROWS,NICE) ! OUT Heat flux through sea-ice
!                                   !   (W per sq m, positive downwards)

      LOGICAL                                                           &
     & L_SICE_HEATFLUX          ! IN T: semi-implicit update of TI

!  External routines called :-
      EXTERNAL TIMER

!  Common and local physical constants.

!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
! C_KAPPAI start

! Thermal conductivity of sea-ice (W per m per K).
        REAL,PARAMETER:: KAPPAI=2.09

! Thermal conductivity of sea water (W per m per K).
        REAL,PARAMETER:: kappas=0.31

! Snow density (Kg per m**3)
        REAL,PARAMETER:: rhosnow=330.0

! Effective thickness of sea-ice surface layer (m).
        REAL,PARAMETER:: DE = 0.1

! C_KAPPAI end
! C_SIECHC has constants for subroutine IMPL_CAL

      ! reciprocal effective areal heat capacity of sea-ice,
      !   ( 1 / (J per sq m per K)).
      REAL,PARAMETER:: AI  = 4.8E-6

! C_SIECHC end
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

      REAL                                                              &
     & ASURF,                                                           &
!                              ! Reciprocal areal heat capacity of
!                              ! sea-ice surface layer (Km2/J).
     & TSAVE                   ! Temporary temperature

      REAL ltfs                ! local copy of constant TFS

      INTEGER I,J,N            ! Loop Counter; horizontal field index.

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SICEHTF ',3)
      ENDIF

      ltfs = access_tfs

      DO J=1,ROWS
       DO I=1,ROW_LENGTH
        if (ocn_sss) THEN
            LTFS = ZeroDegC - 0.054 * auscom_salinity(I,J)
        END IF
        IF (FLANDG(I,J) == 1.0) THEN
          DO N=1,NICE
             SEA_ICE_HTF(I,J,N) = 0.0
             TI(I,J,N)     = RMDI
          ENDDO
        ELSE IF (ICE_FRACTION(I,J) <= 0.0) THEN
          DO N=1,NICE
             SEA_ICE_HTF(I,J,N) = 0.0
             TI(I,J,N)     = LTFS
          ENDDO
          TSTAR_SICE(I,J) = LTFS
        ELSE
          DO N=1,NICE
             IF(ICE_FRACT_NCAT(I,J,N)  >   0.0)THEN
                ASURF              = AI / ICE_FRACT_NCAT(I,J,N)
                IF (L_SICE_HEATFLUX) THEN
                  ! Semi-implicit update of TI
                  TSAVE=TI(I,J,N)
                  TI(I,J,N)=(  TI(I,J,N)+AI*TIMESTEP*(                  &
     &             (SURF_HT_FLUX_SICE_NCAT(I,J,N)/ICE_FRACT_NCAT(I,J,N))+&
     &             ((LTFS-TSAVE*0.5)*KAPPAI/DI_NCAT(I,J,N)) )  )  /     &
     &             ( 1.0+KAPPAI*AI*TIMESTEP*0.5/DI_NCAT(I,J,N) )
                  SEA_ICE_HTF(I,J,N) = ICE_FRACT_NCAT(I,J,N)*KAPPAI*    &
     &                      ((TI(I,J,N)+TSAVE)*0.5 - LTFS)/DI_NCAT(I,J,N)
                ELSE
                  ! Original explicit update of TI 
                  ! (unstable for small ice thicknesses)
                  SEA_ICE_HTF(I,J,N) = ICE_FRACT_NCAT(I,J,N)*KAPPAI*      &
     &                                 (TI(I,J,N) - LTFS)/DI_NCAT(I,J,N)
                  TI(I,J,N)          = TI(I,J,N) + ASURF*TIMESTEP*        &
     &               (SURF_HT_FLUX_SICE_NCAT(I,J,N) - SEA_ICE_HTF(I,J,N))
                END IF
                IF ( TI(I,J,N)  >   TM ) THEN
                  IF (SIMLT) SICE_MLT_HTF(I,J,N) = SICE_MLT_HTF(I,J,N)+ &
     &                        (TI(I,J,N) - TM)/(ASURF*TIMESTEP)
                  TI(I,J,N)        = TM
                ENDIF
             ELSE
                SEA_ICE_HTF(I,J,N) = 0.0
                TI(I,J,N)     = LTFS
          ENDIF
          ENDDO

        ENDIF
       ENDDO
      ENDDO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SICEHTF ',4)
      ENDIF

      RETURN
      END SUBROUTINE SICE_HTF
