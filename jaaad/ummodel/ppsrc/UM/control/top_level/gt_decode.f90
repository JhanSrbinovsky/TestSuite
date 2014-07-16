
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Returns information about a given grid type

      SUBROUTINE GT_DECODE(                                             &
     &                      GRID_TYPE,                                  &
     &                      MODEL_TYPE,CONTENT,COVERAGE,DOMAIN,CYCLIC)






      IMPLICIT NONE

! Given an input GRIDTYPE, this routine will return a value
! for each of MODEL_TYPE, CONTENT, COVERAGE, DOMAIN and CYCLIC
! from the gt_* variables defined in the GRDTYPES comdeck


! Author: Paul Burton
! Current code owner: Paul Burton
!
! History
! Date       Version    Comment
! ----       -------    -------
! 12/10/99   5.0        New deck created.                P.Burton
! 19/09/00   5.2        Added descriptor for grid type 29 :
!                       Orography LBC                        P.Burton
! 06/01/03   5.5        River routing support. P.Selwood.

! Arguments:

      INTEGER                                                           &
     &  GRID_TYPE                                                       &
                   ! IN  : Grid type to investigate
     &, MODEL_TYPE                                                      &
                   ! OUT : What model type grid is for
     &, CONTENT                                                         &
                   ! OUT : What type of data the grid is for
     &, COVERAGE                                                        &
                   ! OUT : What type of points are on the grid
     &, DOMAIN                                                          &
                   ! OUT : What subset of points are on the grid
     &, CYCLIC     ! OUT : If the grid includes extra cyclic wrap
                   !       around points at the start and end of
                   !       each row

! Comdecks
! PARAMETERs defining the GRID_TYPE characteristics

! General
      INTEGER, PARAMETER :: gt_unset = -1
      ! Any value which is unset

! MODEL_TYPE
      INTEGER, PARAMETER :: gt_atmos = 1
      ! Atmosphere field

      INTEGER, PARAMETER :: gt_ocean = 2
      ! Ocean field

      INTEGER, PARAMETER :: gt_wave = 3
      ! Wave field

! CONTENT
      INTEGER, PARAMETER :: gt_thetamass = 1
      ! Contains theta or mass points

      INTEGER, PARAMETER :: gt_velocity = 2
      ! Contains velocity (B grid U,V) points

      INTEGER, PARAMETER :: gt_U_C = 3
      ! Contains U points on C grid

      INTEGER, PARAMETER :: gt_V_C = 4
      ! Contains V points on C grid

      INTEGER, PARAMETER :: gt_hybrid = 5
      ! Points on none of the above
      INTEGER, PARAMETER :: gt_river = 6
      ! River routing grid

! COVERAGE
      INTEGER, PARAMETER :: gt_allpts = 1
      ! All points

      INTEGER, PARAMETER :: gt_land = 2
      ! Land points

      INTEGER, PARAMETER :: gt_sea = 3
      ! Sea points

! DOMAIN
      INTEGER, PARAMETER :: gt_full = 1
      ! Full field

      INTEGER, PARAMETER :: gt_zonal = 2
      ! Zonal field

      INTEGER, PARAMETER :: gt_meridional = 3
      ! Meridional field

      INTEGER, PARAMETER :: gt_ozone = 4
      ! Ozone field

      INTEGER, PARAMETER :: gt_scalar = 5
      ! Single point

      INTEGER, PARAMETER :: gt_compressed = 6
      ! Compressed points

      INTEGER, PARAMETER :: gt_LBC = 7
      ! Lateral Boundary Condition Field

! CYCLIC

      INTEGER, PARAMETER :: gt_nocyclic = 1
      ! No cyclic columns

      INTEGER, PARAMETER :: gt_optcyclic = 2
      ! Optional cyclic columns

      INTEGER, PARAMETER :: gt_cyclic = 3
      ! Includes cyclic columns

! Local variables

      INTEGER max_grid_types
      PARAMETER (max_grid_types=65)

      INTEGER GRID_DATA(5,max_grid_types)
!             Array holding all the descriptions of the different
!             grid types.5 descriptors for each grid type.

      INTEGER                                                           &
     &  ICODE      ! Error code

      CHARACTER*80                                                      &
     &  CMESSAGE   ! Error message

! Definition of the field types
      DATA                                                              &
     &  GRID_DATA(1:5,1)                                                &
     &   /gt_atmos,gt_thetamass,gt_allpts,gt_full,gt_nocyclic/          &
     &, GRID_DATA(1:5,2)                                                &
     &   /gt_atmos,gt_thetamass,gt_land,gt_full,gt_nocyclic/            &
     &, GRID_DATA(1:5,3)                                                &
     &   /gt_atmos,gt_thetamass,gt_sea,gt_full,gt_nocyclic/             &
     &, GRID_DATA(1:5,4)                                                &
     &   /gt_atmos,gt_thetamass,gt_allpts,gt_zonal,gt_nocyclic/         &
     &, GRID_DATA(1:5,5)                                                &
     &   /gt_atmos,gt_thetamass,gt_allpts,gt_meridional,gt_nocyclic/    &
     &, GRID_DATA(1:5,6)  /5*gt_unset/                                  &
     &, GRID_DATA(1:5,7)  /5*gt_unset/                                  &
     &, GRID_DATA(1:5,8)  /5*gt_unset/                                  &
     &, GRID_DATA(1:5,9)  /5*gt_unset/                                  &
     &, GRID_DATA(1:5,10) /5*gt_unset/

      DATA                                                              &
     &  GRID_DATA(1:5,11)                                               &
     &   /gt_atmos,gt_velocity,gt_allpts,gt_full,gt_nocyclic/           &
     &, GRID_DATA(1:5,12)                                               &
     &   /gt_atmos,gt_velocity,gt_land,gt_full,gt_nocyclic/             &
     &, GRID_DATA(1:5,13)                                               &
     &   /gt_atmos,gt_velocity,gt_sea,gt_full,gt_nocyclic/              &
     &, GRID_DATA(1:5,14)                                               &
     &   /gt_atmos,gt_velocity,gt_allpts,gt_zonal,gt_nocyclic/          &
     &, GRID_DATA(1:5,15)                                               &
     &   /gt_atmos,gt_velocity,gt_allpts,gt_meridional,gt_nocyclic/     &
     &, GRID_DATA(1:5,16) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,17)                                               &
     &   /gt_atmos,gt_hybrid,gt_allpts,gt_scalar,gt_nocyclic/           &
     &, GRID_DATA(1:5,18)                                               &
     &   /gt_atmos,gt_U_C,gt_allpts,gt_full,gt_nocyclic/                &
     &, GRID_DATA(1:5,19)                                               &
     &   /gt_atmos,gt_V_C,gt_allpts,gt_full,gt_nocyclic/                &
     &, GRID_DATA(1:5,20) /5*gt_unset/

      DATA                                                              &
     &  GRID_DATA(1:5,21)                                               &
     &   /gt_atmos,gt_thetamass,gt_land,gt_compressed,gt_nocyclic/      &
     &, GRID_DATA(1:5,22)                                               &
     &   /gt_atmos,gt_thetamass,gt_allpts,gt_ozone,gt_nocyclic/         &
     &, GRID_DATA(1:5,23)                                               &
     &   /gt_atmos,gt_river,gt_allpts,gt_full,gt_nocyclic/              &
     &, GRID_DATA(1:5,24) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,25)                                               &
     &   /gt_atmos,gt_hybrid,gt_allpts,gt_LBC,gt_nocyclic/              &
     &, GRID_DATA(1:5,26)                                               &
     &   /gt_atmos,gt_thetamass,gt_allpts,gt_LBC,gt_nocyclic/           &
     &, GRID_DATA(1:5,27)                                               &
     &   /gt_atmos,gt_U_C,gt_allpts,gt_LBC,gt_nocyclic/                 &
     &, GRID_DATA(1:5,28)                                               &
     &   /gt_atmos,gt_V_C,gt_allpts,gt_LBC,gt_nocyclic/                 &
     &, GRID_DATA(1:5,29)                                               &
     &   /gt_atmos,gt_thetamass,gt_allpts,gt_LBC,gt_nocyclic/           &
     &, GRID_DATA(1:5,30) /5*gt_unset/

      DATA                                                              &
     &  GRID_DATA(1:5,31)                                               &
     &   /gt_ocean,gt_thetamass,gt_sea,gt_compressed,gt_nocyclic/       &
     &, GRID_DATA(1:5,32)                                               &
     &   /gt_ocean,gt_velocity,gt_sea,gt_compressed,gt_nocyclic/        &
     &, GRID_DATA(1:5,33) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,34) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,35) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,36)                                               &
     &   /gt_ocean,gt_thetamass,gt_allpts,gt_full,gt_optcyclic/         &
     &, GRID_DATA(1:5,37)                                               &
     &   /gt_ocean,gt_velocity,gt_allpts,gt_full,gt_optcyclic/          &
     &, GRID_DATA(1:5,38)                                               &
     &   /gt_ocean,gt_U_C,gt_allpts,gt_full,gt_optcyclic/               &
     &, GRID_DATA(1:5,39)                                               &
     &   /gt_ocean,gt_V_C,gt_allpts,gt_full,gt_optcyclic/               &
     &, GRID_DATA(1:5,40) /5*gt_unset/

      DATA                                                              &
     &  GRID_DATA(1:5,41)                                               &
     &   /gt_ocean,gt_thetamass,gt_allpts,gt_full,gt_cyclic/            &
     &, GRID_DATA(1:5,42)                                               &
     &   /gt_ocean,gt_velocity,gt_allpts,gt_full,gt_cyclic/             &
     &, GRID_DATA(1:5,43)                                               &
     &   /gt_ocean,gt_thetamass,gt_allpts,gt_zonal,gt_nocyclic/         &
     &, GRID_DATA(1:5,44)                                               &
     &   /gt_ocean,gt_velocity,gt_allpts,gt_zonal,gt_nocyclic/          &
     &, GRID_DATA(1:5,45)                                               &
     &   /gt_ocean,gt_thetamass,gt_allpts,gt_meridional,gt_nocyclic/    &
     &, GRID_DATA(1:5,46)                                               &
     &   /gt_ocean,gt_velocity,gt_allpts,gt_meridional,gt_nocyclic/     &
     &, GRID_DATA(1:5,47)                                               &
     &   /gt_ocean,gt_hybrid,gt_allpts,gt_scalar,gt_nocyclic/           &
     &, GRID_DATA(1:5,48) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,49) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,50) /5*gt_unset/

      DATA                                                              &
     &  GRID_DATA(1:5,51)                                               &
     &   /gt_ocean,gt_hybrid,gt_allpts,gt_LBC,gt_nocyclic/              &
     &, GRID_DATA(1:5,52) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,53) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,54) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,55) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,56) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,57) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,58) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,59) /5*gt_unset/

      DATA                                                              &
     &  GRID_DATA(1:5,60)                                               &
     &   /gt_wave,gt_thetamass,gt_allpts,gt_full,gt_nocyclic/           &
     &, GRID_DATA(1:5,61) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,62)                                               &
     &   /gt_wave,gt_thetamass,gt_sea,gt_compressed,gt_nocyclic/        &
     &, GRID_DATA(1:5,63) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,64) /5*gt_unset/                                  &
     &, GRID_DATA(1:5,65)                                               &
     &   /gt_wave,gt_hybrid,gt_allpts,gt_LBC,gt_nocyclic/

! And the code

      IF ((GRID_TYPE  <   1) .OR.                                       &
     &    (GRID_TYPE  >   max_grid_types)) THEN

        ICODE=1
        WRITE(CMESSAGE,*) 'Invalid GRID_TYPE ',GRID_TYPE,               &
     &                    ' passed to GT_DECODE'

        CALL EREPORT('GT_DECODE',ICODE,CMESSAGE)

      ENDIF

      MODEL_TYPE = GRID_DATA(1,GRID_TYPE)
      CONTENT    = GRID_DATA(2,GRID_TYPE)
      COVERAGE   = GRID_DATA(3,GRID_TYPE)
      DOMAIN     = GRID_DATA(4,GRID_TYPE)
      CYCLIC     = GRID_DATA(5,GRID_TYPE)

      IF (MODEL_TYPE  ==  gt_unset) THEN

        ICODE=2
        WRITE(CMESSAGE,*) 'GRID_TYPE ',GRID_TYPE,                       &
     &                    ' undefined in GT_DECODE'
        CALL EREPORT('GT_DECODE',ICODE,CMESSAGE)

      ENDIF

      RETURN

      END SUBROUTINE GT_DECODE
