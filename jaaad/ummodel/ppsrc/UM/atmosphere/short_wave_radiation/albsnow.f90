

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine to calculate spectral snow albedos for MOSES II, based on
! the Marshall (1989) parametrization of the Wiscombe and Warren (1980)
! model. Influence of contaminants in snow has not been included - see
! UM vn4.5 deck FTSA1A.
!
! History:
! Version  Date     Comment
! =======  ====     =======
! 5.2      3/00     New deck (Richard Essery)
! 5.4  28/08/02     Re-introduce snow soot.  M. Best
! 6.2  21/02/06     Introduced #if defined statements for
!                   version control of radiation code.
!                                  (J.-C.Thelen)
!
! *********************************************************************
      SUBROUTINE ALBSNOW (P_FIELD,LAND_FIELD,LAND_INDEX,                &
     &                    NTILES,TILE_INDEX,TILE_PTS,                   &
     &                    COSZ,RGRAIN,SNOW_TILE,SOOT,ALB_SNOW)

      IMPLICIT NONE

!------------------------ nstypes.h ----------------------------------
!jhan:further renovation of ths file may be necessary params are dependent on dataset
!jhan: ALSO nstypes_cable.h should be unecessary nsoil/soil is only difference
      !--- Number of non-vegetation surface types
      Integer, Parameter :: NNVG  = 4

      !--- Number of plant functional types.
      Integer, Parameter :: NPFT  = 13
      
      !--- Number of surface types.
      Integer, Parameter :: NTYPE =17 
      
      !--- Index of the surface type 'Soil'
      !Integer, Parameter :: SOIL  = 16 
      !dhb599, 20110615: change made as per Peter Vohralik, item 1:
      Integer, Parameter :: SOIL  = 14

!--- Land surface types :
!--- original veg. tiles 
!     1 - Broadleaf Tree
!     2 - Needleleaf Tree
!     3 - C3 Grass
!     4 - C4 Grass
!     5 - Shrub
!--- for testing these tiles are set = 1:5 
!     6 - Broadleaf Tree
!     7 - Needleleaf Tree
!     8 - C3 Grass
!     9 - C4 Grass
!    10 - Shrub
!--- for testing these tiles are set = 0
!    11 - 0 
!    11 - 0
!    11 - 0
!--- original non-veg tiles moved to these indices
!     14 - Urban
!     15 - Water
!     16 - Soil
!     17 - Ice



! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     & P_FIELD                                                          &
                                   ! Total number of grid points.
     &,LAND_FIELD                                                       &
                                   ! Number of land points.
     &,NTILES                      ! Number of surface tiles.

!   Array arguments with intent(in):
      INTEGER                                                           &
     & LAND_INDEX(LAND_FIELD)                                           &
                                   ! Index of land points.
     &,TILE_PTS(NTYPE)                                                  &
                                   ! Number tile points.
     &,TILE_INDEX(LAND_FIELD,NTYPE)! Index of tile points.

      REAL                                                              &
     & COSZ(P_FIELD)                                                    &
                                   ! Zenith cosine.
     &,RGRAIN(LAND_FIELD,NTILES)                                        &
                                   ! Snow grain size (microns).
     &,SNOW_TILE(LAND_FIELD,NTILES)                                     &
                                   ! Lying snow (kg/m2).
     &,SOOT(P_FIELD)               ! Snow soot content (kg/kg).

!   Array arguments with intent(out):
      REAL                                                              &
     & ALB_SNOW(LAND_FIELD,NTYPE,4)! Snow albedo.
!                                  !   (*,*,1) - Direct beam visible
!                                  !   (*,*,2) - Diffuse visible
!                                  !   (*,*,3) - Direct beam near-IR
!                                  !   (*,*,4) - Diffuse near-IR

! Local scalars:
      REAL                                                              &
     & REFF                                                             &
                                   ! Zenith effective grain size.
     &,R0                                                               &
                                   ! Grain size for fresh snow (microns)
     &,SIGMA                       ! Scaled soot content
      PARAMETER ( R0 = 50. )

      INTEGER                                                           &
     & BAND,I,J,L,N                ! Loop counters.

      REAL                                                              &
     & AMAX(2)                     ! Maximum albedo for fresh snow
!                 VIS     NIR
      DATA AMAX / 0.98,   0.7   /

      DO N=1,NTILES
        DO L=1,LAND_FIELD
          I = LAND_INDEX(L)
          IF (SNOW_TILE(L,N)  >   0.) THEN
            REFF = RGRAIN(L,N) * ( 1. + 0.77*(COSZ(I)-0.65) )**2
            ALB_SNOW(L,N,1) = AMAX(1) - 0.002*(SQRT(REFF) - SQRT(R0))
            ALB_SNOW(L,N,2) = AMAX(1) -                                 &
     &                        0.002*(SQRT(RGRAIN(L,N)) - SQRT(R0))
            ALB_SNOW(L,N,3) = AMAX(2) - 0.09*ALOG(REFF/R0)
            ALB_SNOW(L,N,4) = AMAX(2) - 0.09*ALOG(RGRAIN(L,N)/R0)
          ENDIF
        ENDDO
      ENDDO

! Adjust visible snow albedos for soot content
      DO N=1,NTILES
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          I = LAND_INDEX(L)
          SIGMA = SOOT(I) * RGRAIN(L,N) / 0.0017
          IF ( SIGMA  >   1. ) THEN
            ALB_SNOW(L,N,1) = 0.07 +                                    &
     &                        0.5*(ALB_SNOW(L,N,1) - 0.07)/(SIGMA**0.46)
            ALB_SNOW(L,N,2) = 0.07 +                                    &
     &                        0.5*(ALB_SNOW(L,N,2) - 0.07)/(SIGMA**0.46)
          ELSE
            ALB_SNOW(L,N,1) = ALB_SNOW(L,N,1) -                         &
     &                        0.5*(ALB_SNOW(L,N,1) - 0.07)*(SIGMA**0.6)
            ALB_SNOW(L,N,2) = ALB_SNOW(L,N,2) -                         &
     &                        0.5*(ALB_SNOW(L,N,2) - 0.07)*(SIGMA**0.6)
          ENDIF
        ENDDO
      ENDDO

! Adjust near-IR snow albedos for soot content
      DO N=1,NTILES
        DO J=1,TILE_PTS(N)
          L = TILE_INDEX(J,N)
          I = LAND_INDEX(L)
          SIGMA = SOOT(I) * RGRAIN(L,N) / 0.004
          IF ( SIGMA  >   1. ) THEN
            ALB_SNOW(L,N,3) = 0.06 +                                    &
     &                        0.5*(ALB_SNOW(L,N,3) - 0.06)/(SIGMA**0.6)
            ALB_SNOW(L,N,4) = 0.06 +                                    &
     &                        0.5*(ALB_SNOW(L,N,4) - 0.06)/(SIGMA**0.6)
          ELSE
            ALB_SNOW(L,N,3) = ALB_SNOW(L,N,3) -                         &
     &                        0.5*(ALB_SNOW(L,N,3) - 0.06)*(SIGMA**0.7)
            ALB_SNOW(L,N,4) = ALB_SNOW(L,N,4) -                         &
     &                        0.5*(ALB_SNOW(L,N,4) - 0.06)*(SIGMA**0.7)
          ENDIF
        ENDDO
      ENDDO

      IF (NTILES  ==  1) THEN
        DO BAND=1,4
          DO N=2,NTYPE
            DO J=1,TILE_PTS(N)
              L = TILE_INDEX(J,N)
              ALB_SNOW(L,N,BAND) = ALB_SNOW(L,1,BAND)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE ALBSNOW
