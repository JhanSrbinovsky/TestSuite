
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE MOSGRID ------------------------------------------------
!LL
!LL    PURPOSE:
!LL
!LL    THIS ROUTINE CREATES AN ARRAY WITH 1'S TO MARK THE MOS GRID
!LL    IN THE ARRAY DIMENSIONS OF THE INPUT GRID. THIS ARRAY IS USED
!LL    AS A MASK IN STASH TO EXTRACT THE MOS SUBGRID FIELDSFILE.
!LL
!LL    FOR THE MOS global GRID:
!LL
!LL    THE EUROPEAN 'BOX' FROM 69N TO 34.5N, 13.125W TO 26.25E IS
!LL    MARKED THEN 372 INDIVIDUAL 2 BY 2 GRID POINT BOXES AROUND
!LL    STATION POSITIONS. THE FIRST 40 STATIONS ARE PERIPHERAL TO THE
!LL    EUROPEAN MOS ARCHIVE, AND THE NEXT 332 ARE FOR THE WORLDWIDE
!LL    ARCHIVE. THE COMBINED GRID ARRAY IS EXTRACTED FROM THE NWP MODEL
!LL    FIELDS, BEFORE BEING SPLIT BY SUBSEQUENT JOBS IN MOS ARCHIVING
!LL    AND FORECASTING.
!LL
!LL    FOR THE MOS ELF GRID:
!LL
!LL    THE UK 'BOX' IS DEFINED AS A BOX BETWEEN THE ELF ROWS DEFINED BY
!LL    THE TRANSFORM OF POINTS (63.0 N 0.0 E), AND (48.0 N 0.0 E) AND
!LL    BETWEEN ELF COLUMNS DEFINED BY THE LATLON TO EQUATORIAL GRID
!LL    TRANSFORM OF (54.0 N 11.25 W) AND (54.0 N 7.5 E)
!LL    THIS, IN UNROTATED FIELDS WOULD BE THE BOX 64 TO 48 N, 11.25W TO
!LL    7.5 EAST.
!LL
!LL     TESTED UNDER COMPILER:    CFT77
!LL     TESTED UNDER OS VERSION:  UNICOS 5.0
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
!LL VERSION  DATE
!LL
!LL     AUTHOR:  GIL ROSS           DATE:          5 NOVEMBER 1990
!LL
!LL     CODE VERSION NO: 1.2      DATE 15/06/91
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!LL
!LL     PROGRAMMING STANDARD: UM DOC PAPER 3 VERSION 2 (10/8/90)
!LL
!LL     LOGICAL COMPONENTS COVERED: C4
!LL
!LL     PROJECT TASK: C4 (PART OF)
!LL
!LL     EXTERNAL DOCUMENTATION: ?? T.B.A.
!LL
!LL --------------------------------------------------------------------
!*L INTERFACE AND ARGUMENTS: -------------------------------------------
!
      SUBROUTINE MOSGRID(O_GLOBAL,O_ELF,POLE_LAT,POLE_LON,              &
     &                  DELTALAT,DELTALON,TOP_LAT,WEST_LON,NLATS,NLONS, &
     &                  MOSMAP,MAPSIZ,ICODE,CMESSAGE)
!
!*----------------------------------------------------------------------
      IMPLICIT NONE
!
!*L     SUBROUTINE CALLED:

      EXTERNAL LLTOEQ             ! LATLON TO EQUATORIAL GRID TRANSFORM
!*
!*L     ARGUMENT LIST:
!
      LOGICAL                                                           &
     & O_GLOBAL                                                         &
                    ! LOGICAL TRUE FOR global MODEL           INPUT
     &,O_ELF        ! LOGICAL TRUE ELF MODEL                  INPUT
      REAL                                                              &
     & POLE_LAT                                                         &
                    ! ACTUAL LATITUDE OF THE GRID POLE        INPUT
     &,POLE_LON                                                         &
                    ! ACTUAL LONGITUDE OF THE GRID POLE       INPUT
     &,DELTALAT                                                         &
                    ! LATITUDE STEP BETWEEN LATITUDE ROWS     INPUT
     &,DELTALON                                                         &
                    ! LONGITUDE STEP BETWEEN LONGITUDE COLS   INPUT
     &,TOP_LAT                                                          & 
                    ! MOST N'LY LAT, NEW GRID IF APPROPRIATE  INPUT
     &,WEST_LON     ! MOST W'LY LONG, IN DEG E NEW GRID       INPUT
!               ALL IN +VE NORTH AND +VE EAST DEGREES.
      INTEGER                                                           &
     & NLATS                                                            &
                    ! NUMBER OF LATITUDE ROWS                 INPUT
     &,NLONS        ! NUMBER OF LONGITUDE COLUMNS             INPUT
      INTEGER                                                           &
     & MOSMAP(NLONS,NLATS)                                              & 
                           ! MAP OF 1/0'S FOR MOS BIT MAPMASK OUTPUT
     &,MAPSIZ                                                           & 
                           ! NUMBER OF 1'S IN OUTPUT ARRAY    OUTPUT
     &,ICODE               ! ERROR RETURN CODE                OUTPUT
      CHARACTER*80                                                      &
     & CMESSAGE            ! ERROR RETURN MESSAGE             OUTPUT
!*
!L      LOCAL VARIABLES AND PARAMETERS:
!L
      INTEGER NSTATN,IONE
      REAL BOX_EU_N,BOX_EU_S,BOX_EU_W,BOX_EU_E
      REAL BOX_UK_N,BOX_UK_S,BOX_UK_W,BOX_UK_E
      REAL UK_LAM_LT,UK_LAM_LN
      PARAMETER (NSTATN = 372)             ! NO. OF global STATIONS
!L
      PARAMETER (BOX_EU_N = 69.0)          ! NORTHERN LATITUDE \       .
      PARAMETER (BOX_EU_S = 34.5)          ! SOUTHERN    "     | OF MOS
      PARAMETER (BOX_EU_W = 360.0-13.125)  ! WESTERN LONGITUDE | EU BOX
      PARAMETER (BOX_EU_E = 26.25)         ! EASTERN     "     /
!L
      PARAMETER (BOX_UK_N = 63.0)          ! NORTHERN LATITUDE \       .
      PARAMETER (BOX_UK_S = 48.0)          ! SOUTHERN    "     | OF MOS
      PARAMETER (BOX_UK_W = 360.0-11.25)   ! WESTERN LONGITUDE | UK BOX
      PARAMETER (BOX_UK_E = 7.5)           ! EASTERN     "     /
!L
      PARAMETER (UK_LAM_LT = 54.0)         ! LAT OF UK MOS GRID CENTRE
      PARAMETER (UK_LAM_LN = 0.0)          ! LONG OF UK MOS GRID CENTRE
      PARAMETER (IONE = 1)
!L                     THE ABOVE PARAMETERS ARE INITIALISING VALUES
!L                     HERE ARE THE WORKING VARIABLES:
      REAL BOX_NORTH,BOX_SOUTH,BOX_WEST,BOX_EAST ! ACTUAL BOX NAMES.
!L
      REAL OFF_RANGE                             ! OFFSET RANGE CHANGE
      REAL ZERO_LAT                              ! LOCAL PP ZERO LAT
      REAL ZERO_LON                              ! LOCAL PP ZERO LON
      REAL ELF_LAT,ELF_LON                       ! BOX CONVERSION STORES
      REAL PTLAT,PTLON                           ! POINT LAT LONS STORES
      REAL ALTLNG(2,372)                         ! ARRAY OF STATION PTS
      INTEGER IA,IB,IC,ID,IE,IF,IG,I,J           ! DO LOOP VARIABLES
      INTEGER LTPLUS,LTMNUS,LNPLUS,LNMNUS        ! 2X2 BOX INDICES
!L
!L     INDIVIDUAL STATION LAT/LONGS. THE FIRST 40 ARE FOR THE EXTRA
!L     STATIONS AROUND THE EUROPEAN BOX. THE OTHER 332 ARE THE
!L     WORLDWIDE STATION SET.
!L
      DATA ((ALTLNG(I,J),I=1,2),J=1,188)/                               &
     & 52.750,324.500, 57.000,340.000, 47.000,343.000, 65.167,336.433,  &
     & 63.967,337.400, 64.133,338.067, 65.683,341.917, 63.783,341.933,  &
     & 65.750,344.817, 65.267,346.417, 76.767,341.233, 70.417,338.033,  &
     & 65.600,322.366, 70.933,351.333, 74.500,019.017, 71.100,024.000,  &
     & 71.083,028.233, 69.717,029.883, 70.367,031.100, 67.367,026.650,  &
     & 62.667,029.633, 62.400,025.667, 68.967,033.050, 67.133,032.433,  &
     & 67.883,044.133, 64.983,034.783, 64.583,040.500, 61.717,030.717,  &
     & 61.817,034.267, 61.500,038.933, 59.967,030.300, 34.750,032.400,  &
     & 34.583,032.983, 33.817,035.483, 32.683,343.233, 39.450,328.866,  &
     & 38.517,331.283, 38.750,332.917, 37.750,334.283, 36.967,334.833,  &
     & 78.250,015.467, 67.000,309.200, 64.167,308.250, 61.183,314.583,  &
     & 16.733,337.050, 40.967,028.817, 40.117,032.983, 38.500,027.017,  &
     & 37.000,035.417, 56.217,043.817, 55.750,037.567, 56.800,060.633,  &
     & 53.250,050.450, 55.033,082.900, 52.267,104.350, 48.517,135.167,  &
     & 43.117,131.900, 50.400,030.450, 46.483,030.633, 44.500,034.167,  &
     & 49.933,036.283, 48.683,044.350, 47.250,039.817, 43.233,076.933,  &
     & 41.683,044.950, 40.133,044.467, 41.267,069.267, 33.417,036.517,  &
     & 32.000,034.900, 31.783,035.217, 29.550,034.950, 31.983,035.983,  &
     & 25.283,049.483, 24.717,046.717, 29.217,047.983, 34.550,069.217,  &
     & 21.667,039.150, 21.483,039.833, 25.250,051.567, 25.250,055.333,  &
     & 24.433,054.467, 23.583,058.283, 20.667,058.900, 12.833,045.033,  &
     & 34.017,071.583, 32.933,073.717, 31.550,074.333, 24.900,067.133,  &
     & 23.767,090.383, 34.083,074.833, 27.150,077.967, 26.750,080.883,  &
     & 23.067,072.633, 19.117,072.850, 18.533,073.850, 17.450,078.467,  &
     & 13.000,080.183, 12.967,077.583, 06.817,079.883, 47.933,106.983,  &
     & 27.700,085.367, 39.033,125.783, 37.567,126.967, 35.100,129.033,  &
     & 43.050,141.333, 39.717,140.100, 38.267,140.900, 33.583,130.383,  &
     & 31.567,130.550, 21.983,096.100, 16.900,096.183, 18.783,098.983,  &
     & 13.917,100.600, 03.117,101.550, 21.017,105.800, 10.817,106.667,  &
     & 45.750,126.767, 39.467,075.983, 36.050,103.883, 37.783,112.550,  &
     & 41.767,123.433, 39.100,117.167, 36.067,120.333, 29.667,091.133,  &
     & 30.667,104.017, 25.017,102.683, 34.300,108.933, 30.617,114.133,  &
     & 32.000,118.800, 31.167,121.433, 25.033,121.517, 23.133,113.317,  &
     & 28.617,342.250, 28.467,343.750, 28.050,343.433, 27.933,344.616,  &
     & 28.450,346.133, 28.950,346.400, 34.050,353.233, 33.567,352.333,  &
     & 31.617,351.967, 27.200,002.467, 13.483,002.167, 16.717,357.000,  &
     & 18.100,344.050, 14.733,342.500, 13.350,343.200,-15.933,354.333,  &
     &-20.883,055.517, 32.667,013.150, 32.083,020.267, 31.200,029.950,  &
     & 31.283,032.233, 27.050,031.017, 23.967,032.783, 19.583,037.217,  &
     & 15.600,032.550, 15.283,038.917, 11.550,043.150, 08.983,038.800,  &
     & -0.267,036.100, -3.417,037.067, -3.233,040.100, -4.033,039.617,  &
     & -6.217,039.217, -6.867,039.200, -4.667,055.517, -1.967,030.117,  &
     & -3.317,029.317, -4.817,011.900, -4.250,015.250, 00.450,009.417,  &
     & 04.400,018.517, 04.000,009.733, 06.350,002.383, 06.167,001.250,  &
     & 12.350,358.483, 07.733,354.933, 05.250,356.067, 06.250,349.650,  &
     & -8.850,013.233,-15.667,046.350,-18.800,047.483,-18.117,049.400,  &
     &-19.800,034.900,-25.917,032.567,-13.000,028.650,-13.783,033.767/
      DATA ((ALTLNG(I,J),I=1,2),J=189,372)/                             &
     &-15.317,028.450,-15.683,034.967,-17.917,031.133,-18.100,025.850,  &
     &-20.017,028.617,-22.567,017.100,-24.667,025.917,-26.133,028.233,  &
     &-28.800,024.767,-29.100,026.300,-29.967,030.950,-33.983,018.600,  &
     &-33.983,025.600,-33.033,027.833,-40.350,350.116, 71.300,203.217,  &
     & 64.817,212.133, 58.367,225.417, 82.500,297.667, 44.633,296.500,  &
     & 43.667,280.366, 45.467,286.250, 46.167,299.950, 46.800,288.616,  &
     & 46.617,279.200, 47.617,307.267, 53.317,299.583, 49.900,262.767,  &
     & 50.433,255.333, 53.967,258.900, 50.017,249.283, 51.100,245.983,  &
     & 49.183,236.833, 63.750,291.467, 58.750,265.933, 60.833,244.217,  &
     & 60.717,224.933, 25.817,279.717, 30.500,278.300, 32.900,279.967,  &
     & 27.967,277.467, 33.650,275.583, 32.300,273.600, 29.983,269.750,  &
     & 32.317,269.917, 29.300,265.200, 27.767,262.500, 29.533,261.533,  &
     & 32.900,262.967, 31.800,253.600, 32.117,249.067, 33.433,247.983,  &
     & 32.733,242.833, 36.900,283.800, 36.117,273.317, 34.833,267.750,  &
     & 35.400,262.400, 36.083,244.833, 39.183,283.333, 39.883,284.750,  &
     & 38.183,274.267, 39.900,275.800, 38.750,269.633, 39.733,273.733,  &
     & 39.317,265.283, 39.750,255.133, 37.617,237.617, 42.367,288.967,  &
     & 40.500,279.783, 41.417,278.133, 42.933,281.267, 41.983,272.100,  &
     & 42.233,276.667, 41.367,263.983, 40.783,248.033, 43.650,289.683,  &
     & 47.633,242.467, 47.450,237.700, 40.650,286.217, 22.217,262.150,  &
     & 20.667,256.616, 20.983,270.350, 19.150,263.883, 32.367,295.317,  &
     & 25.050,282.533, 23.167,277.650, 19.900,284.850, 18.500,282.083,  &
     & 18.567,287.700, 18.433,290.116, 18.433,294.000, 17.533,271.700,  &
     & 14.583,269.483, 13.700,270.883, 10.000,275.783, 08.917,280.400,  &
     & 16.267,298.483, 13.067,300.517, 10.617,298.650, 12.200,291.033,  &
     & 04.700,285.866, 03.550,283.616, 10.600,293.017, -3.133,299.983,  &
     & -8.067,325.150,-13.017,321.483,-19.833,316.067,-19.800,317.850,  &
     &-20.550,312.567,-23.000,312.866,-22.317,310.933,-22.917,316.833,  &
     &-23.383,308.817,-23.500,313.383,-30.000,308.817, -0.150,281.517,  &
     &-16.517,291.817,-33.383,289.217,-36.767,286.933,-25.267,302.366,  &
     &-34.833,304.000,-27.450,300.950,-75.500,333.350,-60.717,314.400,  &
     &-65.250,295.733,-17.550,210.383,-37.017,174.800,-41.283,174.767,  &
     &-43.483,172.550,-45.767,170.733, -9.433,147.217,-23.800,133.900,  &
     &-34.933,138.517,-33.950,151.183,-37.850,144.733,-35.300,149.183,  &
     &-42.833,147.500, 03.567,098.683, 05.933,116.050, -6.250,106.900,  &
     & -8.750,115.167, 14.517,121.000, 07.117,125.650, 06.900,122.067,  &
     & 62.083,129.750, 54.933,073.400, 22.300,114.167, 35.683,139.767,  &
     & 28.583,077.200, 22.650,088.450, 01.367,103.983, 43.783,087.617,  &
     & 39.933,116.283, 26.267,050.617, 16.900,042.583, 22.783,005.517,  &
     & -7.967,345.600, 12.533,352.050,-20.433,057.667, 30.133,031.400,  &
     & -1.317,036.917, -4.383,015.433,-49.350,070.250, 38.850,282.967,  &
     & 33.933,241.600, 44.883,266.783, 53.317,246.417, 61.167,209.983,  &
     & 19.433,260.917, 17.933,283.217, 21.350,202.067, 04.833,307.633,  &
     &-15.783,312.067,-34.833,301.467,-39.683,286.933,-12.000,282.883,  &
     &-51.817,301.550,-62.200,301.067,-77.850,166.667,-12.400,130.867,  &
     &-27.433,153.083,-31.933,115.950,-17.750,177.450,-27.167,250.567/
!L ---------------------------------------------------------------------
!L
!L
!L   SETTING LOCAL VARIABLES TO GIVE PP ZERO LATS AND LONS.
!L   ENSURING +VE DEGREES EAST.
!L
      ZERO_LAT = TOP_LAT - DELTALAT
      ZERO_LON = WEST_LON - DELTALON
      IF(ZERO_LON  <   0.0) THEN
         ZERO_LON = ZERO_LON + 360.0     ! ENSURING +VE DEG EAST
      ENDIF
!L     COMMENTED OUT DEBUGGING STATEMENTS ARE PREFACED BY CZ
!Z       WRITE(6,'(''1  SUBROUTINE MOSGRID'')')
!Z       WRITE(6,'(''0    GLOBAL MODEL ? ..........'',L5)') O_GLOBAL
!Z       WRITE(6,'(''     ELF MODEL ? .............'',L5)') O_ELF
!Z       WRITE(6,'(''0    LATITUDE OF GRID POLE....'',F10.4)') POLE_LAT
!Z       WRITE(6,'(''     LONGITUDE OF GRID POLE...'',F10.4)') POLE_LON
!Z       WRITE(6,'(''0    NORTHMOST LATITUDE.......'',F10.4)') TOP_LAT
!Z       WRITE(6,'(''     ZEROETH   LATITUDE (PP)..'',F10.4)') ZERO_LAT
!Z       WRITE(6,'(''     LATITUDE STEPS...........'',F10.4)') DELTALAT
!Z       WRITE(6,'(''     NUMBER OF LATITUDES......'',I6)')   NLATS
!Z       WRITE(6,'(''0    WESTMOST LONGITUDE.......'',F10.4)') WEST_LON
!Z       WRITE(6,'(''     ZEROETH  LONGITUDE (PP)..'',F10.4)') ZERO_LON
!Z       WRITE(6,'(''     LONGITUDE STEPS..........'',F10.4)') DELTALON
!Z       WRITE(6,'(''     NUMBER OF LONGITUDES.....'',I6)')   NLONS
!L ---------------------------------------------------------------------
!L
!L     CHOOSE BOX SIZES FOR ELF OR global. IF NEITHER,SET ERROR FLAG
!L     ZERO SIZES AND RETURN
!L
!L ---------------------------------------------------------------------
      ICODE = 0
      CMESSAGE = ' '
      MAPSIZ=0
!L
      IF(O_GLOBAL .AND. O_ELF) THEN
        WRITE(6,'(''1    -------  BOTH GLOBAL AND ELF PICKED -----'')')
         WRITE(6,'(''0    SUBROUTINE MOSGRID ERROR'')')
         MAPSIZ = 1
         ICODE = 1
         CMESSAGE = 'MOSGRID: BOTH GLOBAL AND ELF LOGICALS WERE SET '
         GOTO 99
      ELSE IF(O_GLOBAL) THEN
         IF(POLE_LAT  ==  90.0) THEN
!Z          WRITE(6,'(''0    GLOBAL MODEL HAS BEEN CHOSEN'')')
            BOX_NORTH = BOX_EU_N
             BOX_SOUTH = BOX_EU_S
              BOX_EAST = BOX_EU_E
               BOX_WEST = BOX_EU_W
         ELSE
           WRITE(6,'(''1    ----- GLOBAL MODEL POLE NOT 90.0 -----'')')
            WRITE(6,'(''0    SUBROUTINE MOSGRID ERROR'')')
            MAPSIZ = 0
            ICODE = 1
            CMESSAGE ='MOSGRID: GLOBAL MODEL WITH ROTATED GRID:- NO MOS'
            GOTO 99
         ENDIF
      ELSE IF (O_ELF) THEN
!Z       WRITE(6,'(''0    ELF MODEL HAS BEEN CHOSEN'')')
!L ---------------------------------------------------------------------
!L
!L    THE BOX IS DEFINED AS THE (ROTATED) ELF GRID BOX AROUND THE
!L    CENTRAL UK POINT. WE FIND THE THE BOX GRID_LATS AND GRID_LONS WITH
!L
!L    A)  (BOX_UK_N,UK_LAM_LN) AND TAKE THE ELF GRID_LAT AS BOX_NORTH
!L    B)  (BOX_UK_S,UK_LAM_LN) AND TAKE THE ELF GRID_LAT AS BOX_SOUTH
!L    C)  (UK_LAM_LT,BOX_UK_W) AND TAKE THE ELF GRID_LONG AS BOX_WEST
!L    D)  (UK_LAM_LT,BOX_UK_E) AND TAKE THE ELF GRID_LONG AS BOX_EAST
!L
!L
!L ---------------------------------------------------------------------
! DEPENDS ON: lltoeq
           CALL LLTOEQ(BOX_UK_N,UK_LAM_LN,                              &
     &                 ELF_LAT,ELF_LON,POLE_LAT,POLE_LON,IONE)
!Z         WRITE(6,'('' BOX N '',6F9.3)')
!Z   *                 BOX_UK_N,UK_LAM_LN,
!Z   *                 ELF_LAT,ELF_LON,POLE_LAT,POLE_LON
           BOX_NORTH = ELF_LAT
! DEPENDS ON: lltoeq
           CALL LLTOEQ(BOX_UK_S,UK_LAM_LN,                              &
     &                 ELF_LAT,ELF_LON,POLE_LAT,POLE_LON,IONE)
!Z         WRITE(6,'('' BOX S '',6F9.3)')
!Z   *                 BOX_UK_S,UK_LAM_LN,
!Z   *                 ELF_LAT,ELF_LON,POLE_LAT,POLE_LON
           BOX_SOUTH = ELF_LAT
! DEPENDS ON: lltoeq
           CALL LLTOEQ(UK_LAM_LT,BOX_UK_W,                              &
     &                 ELF_LAT,ELF_LON,POLE_LAT,POLE_LON,IONE)
!Z         WRITE(6,'('' BOX W '',6F9.3)')
!Z   *                 UK_LAM_LT,BOX_UK_W,
!Z   *                 ELF_LAT,ELF_LON,POLE_LAT,POLE_LON
           BOX_WEST = ELF_LON
! DEPENDS ON: lltoeq
           CALL LLTOEQ(UK_LAM_LT,BOX_UK_E,                              &
     &                 ELF_LAT,ELF_LON,POLE_LAT,POLE_LON,IONE)
!Z         WRITE(6,'('' BOX E '',6F9.3)')
!Z   *                 UK_LAM_LT,BOX_UK_E,
!Z   *                 ELF_LAT,ELF_LON,POLE_LAT,POLE_LON
           BOX_EAST = ELF_LON
      ELSE
        WRITE(6,'(''1    ----- NEITHER GLOBAL NOR ELF PICKED -----'')')
         WRITE(6,'(''0    SUBROUTINE MOSGRID ERROR'')')
         MAPSIZ = 0
         ICODE = 1
         CMESSAGE = 'MOSGRID: NEITHER GLOBAL NOR ELF LOGICALS WERE SET'
         GOTO 99
      ENDIF
!L ---------------------------------------------------------------------
!L
!L  SINCE THE BOX IS AROUND THE UK, THE PROBLEM OF ZERO LONGITUDE
!L  WILL ARISE. IF HOWEVER THE GRID IS ELF, THEN THE ACTUAL VALUE OF
!L  THE ROTATED GRID MAY WELL MAKE BOX_WEST LESS THAN BOX_EAST.
!L  TO RESOLVE THE INCLUSION PROBLEM, THEN IF BOX_WEST > BOX_EAST, THEN
!L  THE PRIME RANGE OF THE VALUES WILL BE 180 TO 540, RATHER THAN
!L  0 TO 360. ALL LONGITUDES LESS THAN 180 WILL HAVE AN ADDITION OF
!L  360 DEGREES AS A RANGE OFFSET.
!L
!L ---------------------------------------------------------------------
      IF(BOX_WEST  >   BOX_EAST) THEN
         BOX_EAST = BOX_EAST+360.0
         OFF_RANGE = 360.0
      ELSE
         OFF_RANGE = 0.0
      ENDIF
!L
!Z    WRITE(6,'(''0    GLOBAL MODEL ? ..........'',L5)')     O_GLOBAL
!Z    WRITE(6,'(''     ELF MODEL ? .............'',L5)')     O_ELF
!Z    WRITE(6,'(''0    ZEROETH   LATITUDE.......'',F10.4)')  ZERO_LAT
!Z    WRITE(6,'(''     LATITUDE STEPS...........'',F10.4)')  DELTALAT
!Z    WRITE(6,'(''     NUMBER OF LATITUDES......'',I6)')     NLATS
!Z    WRITE(6,'(''0    ZEROETH  LONGITUDE.......'',F10.4)')  ZERO_LON
!Z    WRITE(6,'(''     LONGITUDE STEPS..........'',F10.4)')  DELTALON
!Z    WRITE(6,'(''     NUMBER OF LONGITUDES.....'',I6)')     NLONS
!Z    WRITE(6,'(''0    NORTHMOST BOX LAT........'',F10.4)')  BOX_NORTH
!Z    WRITE(6,'(''     SOUTHMOST BOX LAT........'',F10.4)')  BOX_SOUTH
!Z    WRITE(6,'(''     WESTMOST  BOX LON........'',F10.4)')  BOX_WEST
!Z    WRITE(6,'(''     EASTMOST  BOX LON........'',F10.4)')  BOX_EAST
!Z    WRITE(6,'(''     OFFSET FOR RANGE CHANGE..'',F10.4)')  OFF_RANGE
!L ---------------------------------------------------------------------
!L
!L     INITIALISING MOSMAP TO ZERO.
!L
!L ---------------------------------------------------------------------
      DO IA = 1, NLATS
          DO IB = 1, NLONS
            MOSMAP(IB,IA) = 0
         ENDDO
      ENDDO
!L ---------------------------------------------------------------------
!L
!L     THIS LOOPS THE LATS. IF THE BOX CONTAINS THIS LAT, THE INNER DOES
!L     THE LONS.IF THE DEFINITION OF THE LONS GIVES THE POINT OUTSIDE
!L     THE RANGE 0 TO 360, IT IS RESET INTO THE RANGE, BEFORE ANY CHANGE
!L     OF RANGE FROM 180 TO 540
!L     CORRECTION 06/03/91 TO COPE WITH START OF FULL GRID OFFSET
!L        CHANGES  (IC-1) TO IC AND (ID-1) TO ID
!L
!L ---------------------------------------------------------------------
      DO IC = 1, NLATS                    ! TOP LATITUDE IS 0TH ROW
         PTLAT = ZERO_LAT +IC*DELTALAT ! LAT START + NUMBER OF STEPS
!Z       WRITE(6,'('' IC,LAT,BOXN,BOXS'',I8,3F8.3)')
!Z   *                IC,PTLAT,BOX_NORTH,BOX_SOUTH
!L
!L      PRECISION PROBLEMS. ALLOW *JUST* OUTSIDE BOX. WITH DELTALAT -VE
!L
         IF(((PTLAT - 0.01*DELTALAT)  >=  BOX_SOUTH) .AND.              &
     &      ((PTLAT + 0.01*DELTALAT)  <=  BOX_NORTH)) THEN
            DO ID = 1, NLONS                    ! WESTMOST LON- 0TH COL
              PTLON = ZERO_LON + ID*DELTALON ! START + NUMBER OF STEP S
               IF(PTLON  >   360.0) THEN
                  PTLON = PTLON - 360.0         ! WRAP AROUND
               ELSE IF(PTLON  <   0.0) THEN
                  PTLON = PTLON + 360.0
               ENDIF
               IF(PTLON  <=  180.0) THEN
                  PTLON = PTLON + OFF_RANGE     ! IF BOX CROSSES 0 E
               ENDIF
!Z             WRITE(6,'('' ID,LON,BOXW,BOXE'',I8,3F8.3)')
!Z   *                      ID,PTLON,BOX_WEST,BOX_EAST
!L
!L      PRECISION PROBLEMS AGAIN.
!L
               IF(((PTLON + 0.01*DELTALON)  >=  BOX_WEST) .AND.         &
     &            ((PTLON - 0.01*DELTALON)  <=  BOX_EAST)) THEN
!Z               WRITE(6,'(''   POINT SET IC,ID,LAT,LON'',2I5,2F10.4)')
!Z   *                                     IC,ID,PTLAT,PTLON
                  MOSMAP(ID,IC) = 1      ! SET VALUES
               ENDIF
            ENDDO
         ENDIF
      ENDDO
!Z    WRITE(6,'(''   BOX FINISHED '')')
!L ---------------------------------------------------------------------
!L
!L     THERE ARE STATIONS IN THE BITMAP WITH 2*2 GRID POINT BOXES
!L     *ONLY* FOR THE global GRID.
!L
!L ---------------------------------------------------------------------
      IF(O_GLOBAL) THEN
!L ---------------------------------------------------------------------
!L
!L     SO, FOR A FILLED global BITMAP.....
!L
!L     LOOP THE INDIVIDUAL STATIONS. THE 2 BY 2 GRID POINT BOXES
!L     AROUND THE STATION POSITION ARE SET IN MOSMAP.
!L
!L ---------------------------------------------------------------------
         DO IE = 1, NSTATN
          LTMNUS = (ZERO_LAT- ALTLNG(1,IE))/DELTALAT                    &
                                                     ! TRUNC NORTHWARDSS
     &                     -0.001                   ! ROUNDING PROTECT
           LTMNUS = -LTMNUS                         ! DELTALT NEGATIVES
            LTPLUS = LTMNUS + 1                     ! GIVES S_WARD ROW
          LNMNUS = (ALTLNG(2,IE) - ZERO_LON                             &
                                                    ! TRUNC WESTWARDS
     &                       +720.0)/DELTALON                           &
                                                    ! ENSURING +VE
     &                       +0.001                 ! ROUNDING PROTECT
            LNPLUS = LNMNUS + 1                     ! GIVES E_WARD COL
!L
!L            MODULO 360 (INDEXES IN THE RANGE 1 TO NLONS)
!L
            LNMNUS = MOD(LNMNUS,NLONS)
               IF(LNMNUS  ==  0) LNMNUS = NLONS
            LNPLUS = MOD(LNPLUS,NLONS)
               IF(LNPLUS  ==  0) LNPLUS = NLONS
!L
            MOSMAP(LNMNUS,LTMNUS) = 1         ! NORTHWEST CORNER
             MOSMAP(LNPLUS,LTMNUS) = 1        ! NORTHEAST CORNER
              MOSMAP(LNMNUS,LTPLUS) = 1       ! SOUTHWEST CORNER
               MOSMAP(LNPLUS,LTPLUS) = 1      ! SOUTHEAST CORNER
!Z              WRITE(6,'(''  POINT '',I5,'' BOX GRID '',4I5)')
!Z   *                IE,LTMNUS,LTPLUS,LNMNUS,LNPLUS
         ENDDO
      ENDIF
!L ---------------------------------------------------------------------
!L
!L     COUNTING THE NUMBER OF POINTS SET.
!L
!L ---------------------------------------------------------------------
      MAPSIZ=0
      DO IF = 1, NLATS
         DO IG = 1, NLONS
          IF(MOSMAP(IG,IF)  /=  0) MAPSIZ=MAPSIZ+1
         ENDDO
      ENDDO
!Z    WRITE(6,'(''  NUMBER OF POINTS SET '',I8)') MAPSIZ
   99 RETURN
!L ---------------------------------------------------------------------
      END SUBROUTINE MOSGRID
