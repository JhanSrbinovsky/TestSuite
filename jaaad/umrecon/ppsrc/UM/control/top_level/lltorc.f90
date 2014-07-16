
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Convert lat/long specification to row/column specification
! Subroutine Interface:

      SUBROUTINE LLTORC(GRID_CODE,                                      &
     &                  NORTH_LAT_IN,SOUTH_LAT_IN,                      &
     &                  WEST_LONG_IN,EAST_LONG_IN,                      &
     &                  START_ROW,END_ROW,START_COL,END_COL             &

     &                 ,RECON_GRID                                      &

     &                 )

      Use Rcf_Grid_Type_Mod, Only :                                     &
     &    Grid_Type

      Use Rcf_Model_Mod, Only :                                         &
     &    ZonAvOzone,           H_A_Polelat,                            &
     &    H_Global,                                                     &
     &    H_A_EWSpace,          H_A_FirstLong,                          &
     &    H_A_FirstLat,                                                 &
     &    H_A_PoleLong,         H_A_NSSpace

      Use Rcf_Submodel_Mod, Only :                                      &
     &    A_IM

      Use Rcf_Global_To_Local_Mod, Only :                               &
     &    Rcf_Get_Fld_Type

      Use Rcf_Parvars_Mod, Only :                                       &
     &    fld_type_p,                                                   &
     &    fld_type_u,                                                   &
     &    fld_type_v

      Use Ereport_Mod, Only :                                           &
     &    Ereport


      IMPLICIT NONE

! Description:
!   Uses the gridpoint code, the lat/long spec of the required area,
!   and the model area sizes to calculate the model row/column
!   numbers for the area.
!   Fix for the ocean model, since N is at the bottom.
!   Called by PRELIM, INPUTL, ADDRES.
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:

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

! Subroutine arguments:


      INTEGER                                                           &
     &  GRID_CODE                                                       &
                    ! IN : Grid type code from STASHmaster

     &, NORTH_LAT_IN                                                    &
                       ! IN : Latitude of Northern boundary
     &, SOUTH_LAT_IN                                                    &
                       ! IN : Latitude of Southern boundary
     &, WEST_LONG_IN                                                    &
                       ! IN : Longitude of Western boundary
     &, EAST_LONG_IN                                                    &
                       ! IN : Longitude of Eastern boundary

     &, START_ROW                                                       &
                    ! OUT : Row number of start of area
     &, END_ROW                                                         &
                    ! OUT : Row number of end of area
     &, START_COL                                                       &
                    ! OUT : Column number of start of area
     &, END_COL     ! OUT : Column number of end of area

      TYPE (GRID_TYPE) RECON_GRID

! Local variables
      INTEGER                                                           &
     &  MODEL_TYPE                                                      &
                    ! model type of grid
     &, CONTENT                                                         &
                    ! content type of grid
     &, COVERAGE                                                        &
                    ! coverage of grid
     &, DOMAIN                                                          &
                    ! domain of grid
     &, CYCLIC                                                          &
                    ! does grid contain cyclic wrap columns
     &, fld_type    ! P, U or V field type

      INTEGER                                                           &
     &  NORTH_LAT                                                       &
                   ! Modifiable copies
     &, SOUTH_LAT                                                       &
                   ! of the input arguments
     &, WEST_LONG                                                       &
                   ! so that the format can
     &, EAST_LONG                                                       &
                   ! be changed if necessary

     &, NROWS                                                           &
                   ! number of rows
     &, NCOLS      ! number of columns

      LOGICAL                                                           &
     &  LAM_MODEL  ! True if this is a LAM configuration

      REAL                                                              &
     &  POLE_LAT                                                        &
                   ! Latitude of rotated pole (LAM only)
     &, POLE_LONG  ! Longitude of rotated pole (LAM only)

      INTEGER                                                           &
     &  ROW_ORDERING                                                    &
                      ! ordering North->South or South->North
     &, North_to_South                                                  &
                       ! indicates North->South ordering
     &, South_to_North ! indicates South->North ordering

      PARAMETER                                                         &
     &  (North_to_South=1,South_to_North=2)

      REAL                                                              &
     &  START_LAT                                                       &
                   ! Starting latitude
     &, START_LONG                                                      &
                   ! Starting longitude
     &, DEL_LAT                                                         &
                   ! Latitude grid spacing
     &, DEL_LONG                                                        &
                   ! Longitude grid spacing
     &, R_START_ROW                                                     &
                    ! first row number
     &, R_END_ROW                                                       &
                    ! last row number
     &, R_START_COL                                                     &
                    ! first column number
     &, R_END_COL   ! last column number

! External and subroutine calls

      EXTERNAL LLTOLL,LLTOEQ

!-----------------------------------------------------------------------

! Copy input arguments
      NORTH_LAT=NORTH_LAT_IN
      SOUTH_LAT=SOUTH_LAT_IN
      WEST_LONG=WEST_LONG_IN
      EAST_LONG=EAST_LONG_IN

      IF (WEST_LONG  <   0) WEST_LONG=WEST_LONG+360
      IF (EAST_LONG  <=  0) EAST_LONG=EAST_LONG+360

! Get information about grid type

! DEPENDS ON: gt_decode
      CALL GT_DECODE(GRID_CODE,                                         &
     &               MODEL_TYPE,CONTENT,COVERAGE,DOMAIN,CYCLIC)

! Get information of field type for PARVARS variables

      fld_type=RCF_GET_FLD_TYPE(GRID_CODE)

! Find start latitude and longitude of full grid, grid spacing and
! size of grid

      IF (MODEL_TYPE  ==  gt_atmos) THEN

        START_LAT=H_A_FIRSTLAT
        START_LONG=H_A_FIRSTLONG
        DEL_LAT=H_A_NSSPACE
        DEL_LONG=H_A_EWSPACE
        IF (fld_type == fld_type_p) THEN
          NROWS=RECON_GRID % GLOB_P_ROWS
          NCOLS=RECON_GRID % GLOB_P_ROW_LENGTH
        ELSE IF (fld_type == fld_type_u) THEN
          NROWS=RECON_GRID % GLOB_U_ROWS
          NCOLS=RECON_GRID % GLOB_U_ROW_LENGTH
        ELSE IF (fld_type == fld_type_v) THEN
          NROWS=RECON_GRID % GLOB_V_ROWS
          NCOLS=RECON_GRID % GLOB_V_ROW_LENGTH
        END IF
        ROW_ORDERING=South_to_North

        IF (H_GLOBAL(A_IM) == 'N') THEN   ! LAM Configuration
          LAM_MODEL=.TRUE.
          POLE_LAT=H_A_POLELAT
          POLE_LONG=H_A_POLELONG
        ELSE
          LAM_MODEL=.FALSE.
        ENDIF

      ENDIF

! Ensure that DEL_LAT is always positive
! (ie. doesn't imply a direction)
      IF (DEL_LAT  <   0) DEL_LAT=-DEL_LAT

! This assumes the mass grid. The start latitude and longitude
! may be offset for velocity/U/V fields

      IF (CONTENT  ==  gt_velocity) THEN  ! B grid U/C
        START_LAT=START_LAT+DEL_LAT/2
        START_LONG=START_LONG+DEL_LONG/2
      ELSEIF (CONTENT  ==  gt_U_C) THEN   ! C grid U
        START_LONG=START_LONG+DEL_LONG/2
      ELSEIF( CONTENT  ==  gt_V_C) THEN   ! C grid V
        START_LAT=START_LAT-DEL_LAT/2
      ENDIF

! This assumes full domain. Now take account of zonal,meridional
! and scalar fields

      IF (DOMAIN  ==  gt_zonal) THEN
        NCOLS=1
      ELSEIF (DOMAIN  ==  gt_meridional) THEN
        NROWS=1
      ELSEIF (DOMAIN  ==  gt_scalar) THEN
        NCOLS=1
        NROWS=1
      ELSEIF (DOMAIN  ==  gt_ozone) THEN
        IF (ZonAvOzone) NCOLS=1
      ENDIF

! This is the global sizes - this may be all we need

      IF ((NORTH_LAT  ==  90) .AND. (SOUTH_LAT  ==  -90) .AND.          &
     &    (WEST_LONG  ==  0) .AND. (EAST_LONG  ==  360)) THEN

        START_ROW=1
        END_ROW=NROWS
        START_COL=1
        END_COL=NCOLS

      ELSE ! Not a simple case

! If this is a LAM configuration we need to transform the latitudes
! and longitudes relative to the rotated pole

        IF (LAM_MODEL) THEN
! DEPENDS ON: lltoll
          CALL LLTOLL(                                                  &
     &        NORTH_LAT,SOUTH_LAT,EAST_LONG,WEST_LONG,                  &
     &        POLE_LAT,POLE_LONG)
        ENDIF

! Make sure that DEL_LAT has a sign consistent with the
! grid ordering. Further up we have ensured that DEL_LAT
! is positive. Now we change it negative if necessary.
!               ( North->South => -ive
!                 South->North => +ive )

        IF (ROW_ORDERING  ==  North_to_South) DEL_LAT=-DEL_LAT

! Calculate the start and end rows

        IF (ROW_ORDERING  ==  North_to_South) THEN

          R_START_ROW=1.0   + (NORTH_LAT-START_LAT)/DEL_LAT
          R_END_ROW  =1.999 + (START_LAT-SOUTH_LAT)/DEL_LAT

        ELSE ! South->North

          R_START_ROW=1.0   + (SOUTH_LAT-START_LAT)/DEL_LAT
          R_END_ROW  =1.999 + (NORTH_LAT-START_LAT)/DEL_LAT

        ENDIF

        START_ROW=MAX(1.0,R_START_ROW)
        END_ROW=MIN(REAL(NROWS),R_END_ROW)

        IF (START_LONG  <   0) START_LONG=START_LONG+360.0

        IF ((REAL(START_LONG) + DEL_LONG*(NCOLS-1))  <=  360.0)         &
     &  THEN  ! If the total model domain doesn't cross
              ! the zero meridian

          R_START_COL=1.0 +   (WEST_LONG-START_LONG)/DEL_LONG
          R_END_COL=  1.999 + (EAST_LONG-START_LONG)/DEL_LONG

        ELSE ! model domain crosses the zero meridian

          IF (REAL(WEST_LONG)  <   START_LONG) THEN
            ! If the start is before the zero meridian

            R_START_COL=1.0  + (WEST_LONG + 360.0 - START_LONG)/        &
     &                         DEL_LONG

          ELSE ! the start lies after the meridian

            R_START_COL=1.0 + (REAL(WEST_LONG)-START_LONG)/             &
     &                        DEL_LONG

          ENDIF

          IF (REAL(EAST_LONG)  <   START_LONG) THEN
            ! If the end is after the zero meridian

            R_END_COL=1.0 + (EAST_LONG + 360.0 - START_LONG)/           &
     &                      DEL_LONG

          ELSE ! the end lies before the zero meridian

            R_END_COL=1.0 + (REAL(EAST_LONG) - START_LONG)/             &
     &                      DEL_LONG

          ENDIF

        ENDIF

        START_COL=MIN(REAL(NCOLS),MAX(1.0,R_START_COL))
        END_COL=MIN(REAL(NCOLS),MAX(1.0,R_END_COL))

      ENDIF

      RETURN
      END SUBROUTINE LLTORC

!- End of subroutine code ----------------------------------------


!+Use subroutine LLTOEQ to convert to equatorial lat/lon for LAM
! Subroutine Interface:


!- End of subroutine code ---------------------------------------------
