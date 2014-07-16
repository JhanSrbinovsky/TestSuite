
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      SUBROUTINE RIV_ROUT(SURF_RUNOFFIN, SUB_RUNOFFIN,                  &
     &          ICOLS,JROWS,FRAC,RMDI,                                  &
     &          INVERT_ATMOS,XTA,YTA,XUA,YUA,XVA,YVA,                   &
     &          GRIDBOX_AREAS,                                          &
     &          REGRID, CYCLIC_TRIP, GLOBAL_TRIP, INVERT_TRIP,          &
     &          IMTR, JMTR, RU, RATMED, DT,                             &
     &          R_RIVDIR, R_RIVSEQ, TWATSTOR,                           &
     &          RIVEROUT_TRIP, RUNOFF_TRIP_KGS,                         &
!         put new inland basin variables in river1a arguments
     &              RIVEROUT_INV,INLANDOUT_TRIP,INLANDOUT_ATM)




!
! Purpose:
!
! Perform the routing of surface and sub-surface runoff.
!
! Method:
! This routine regrids the total runoff to the river routing grid
! and passes it to the TRIP routines to be routed, then maps the
! outflow to seapoints on the atmos grid and passes back the
! updated water storage. The subroutine DO_MAP_MAX contained in
! this deck uses the mapping information obatained by PREAV1.dk
! to map the river outflow (Kg/s) produced on the river routing
! grid onto sea points in the atmos grid. This is then converted
! to the usual Kg/m2/s using the atmos gridbox areas.
!
!  Author  C.B.Bunton 28.02.03
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!-----------------------------------------------------------------

      IMPLICIT NONE

      INTEGER                                                           &
     & ICOLS                                                            &
                             ! IN NO. OF COLUMNS IN ATMOSPHERE
     &,JROWS                                                            &
                            ! IN NO. OF ROWS IN ATMOSPHERE (TP GRID)
     &,IMTR                                                             &
                             ! IN NO. OF COLUMNS IN RIVER GRID
     &,JMTR                  ! NO. OF ROWS IN RIVER GRID

      REAL                                                              &
     & RMDI                                                             &
                            ! IN real missing data indicator
     &,RU                                                               &
                            ! IN river velocity (m/s)
     &,RATMED                                                           &
                            ! IN meandering coefficient
     &,DT                   ! IN river routing timestep(s)
      REAL                                                              &
     & SURF_RUNOFFIN(ICOLS,JROWS)                                       &
                                  ! IN SURFACE RUNOFF (KG/M2/S=mm/s)
     &,SUB_RUNOFFIN(ICOLS,JROWS)                                        &
                                 ! IN SUB_SURFACE RUNOFF (KG/M2/S=mm/s)
     &,R_RIVDIR(IMTR,JMTR)                                              &
                            ! IN Real values of river direction
     &,R_RIVSEQ(IMTR,JMTR)                                              &
                            ! IN Real values of river sequence
     &,TWATSTOR(IMTR,JMTR)                                              &
                            ! IN Initial Water Storage file
     &,XUA(0:ICOLS)                                                     &
                            ! IN Atmosphere UV longitude coordinates
     &,YUA(JROWS)                                                       &
                          ! IN Atmosphere UV latitude coordinates
     &,XTA(ICOLS+1)                                                     &
                            ! IN Atmosphere TP longitude coordinates
     &,YTA(JROWS)                                                       &
                            ! IN Atmosphere TP latitude coordinates
     &,XVA(ICOLS+1)                                                     &
                           ! Atmosphere V longitude coordinates
     &,YVA(0:JROWS)                                                     &
                           ! Atmosphere V latitude coordinates
     &,FRAC(ICOLS,JROWS)                                                &
                           ! Fractional land/sea mask
     &,GRIDBOX_AREAS(ICOLS,JROWS) !IN atmos gridbox areas

      LOGICAL                                                           &
     & INVERT_ATMOS                                                     &
                           ! IN state of atmos fields needed for
!                        ! regridding runoff from atmos.
     &, INVERT_TRIP                                                     &
                           ! TRUE WHEN ROW INVERSION IS REQUIRED
     &, REGRID                                                          &
                           ! TRUE if TRIP grid different to atmos
     &, CYCLIC_TRIP                                                     &
                           ! TRUE WHEN THE TRIP MODEL HAS CYCLIC
     &, GLOBAL_TRIP        ! TRUE WHEN TRIP GRID SURFACE IS SPHER

!      REAL RMDI_TRIP
      REAL                                                              &
     & RIVEROUT_INV(ICOLS,JROWS)                                        &
                                ! OUT river flow out from each gridbox
!                           ! (KG/m2/S)
     &,RIVEROUT_TRIP(IMTR, JMTR)                                        &
                                 ! OUT gridbox outflow on river grid
!                                ! (Kg/s)
     &,RUNOFF_TRIP_KGS(IMTR, JMTR)                                      &
                                  ! OUT gridbox runoff on river grid
!                                ! (Kg/s)
! add variables for inland basins
     &,INLANDOUT_TRIP(IMTR,JMTR)                                        &
                                           !OUT TRIP OUTFLOW
!                            FROM INLAND BASINS ON TRIP GRID kg/s
     &,INLANDOUT_ATM(ICOLS,JROWS)        !OUT TRIP OUTFLOW
!                   FROM INLAND BASINS ON atmos GRID kg/m2/s


      EXTERNAL ROUTEDBL

!      LOCAL

      INTEGER                                                           &
     & IERR,ICODE
!
      CHARACTER*256                                                     &
     & CMessage         ! Error message if return code >0
      Character(*) RoutineName
      Parameter ( RoutineName='RIV_INTCTL')
!
      INTEGER NDEV                ! No. TRIP timesteps/day

       REAL LATBDY                ! SOUTHERN BOUNDARY OF ROUTING
       PARAMETER(LATBDY=-90.0)
      INTEGER                                                           &
     & I,J,K,L,                                                         &
                             ! LOOP COUNTERS
     & IRT,IRU               ! NUMBER OF COLUMNS OF DISTINCT VALUES
                             ! on ocean TS and UV grids
      LOGICAL                                                           &
     & AMASK(ICOLS,JROWS)                                               &
                          ! IN atmos MODEL LAND-SEA MASK FOR atmos grid
                          ! TRUE if FRAC >0.0
     &,TRMASK_ALL(IMTR,JMTR)                                            &
                              ! Set TRIP mask to all sea for interp.
     &,AMASK_ALL(ICOLS,JROWS) ! Set atmos mask to all sea for interp.

      INTEGER                                                           &
     & AMINT(ICOLS*JROWS)  ! INTEGER LAND-SEA MASK ON ATMOSPHERE GRID
!
      LOGICAL TRMASK(IMTR,JMTR) ! TRUE IF POINT IS LAND IN TRIP MODEL
!                               ! (Includes River mouths)
      LOGICAL AMASK_SEA(ICOLS,JROWS)  ! atmos L/sea mask, F for any sea

      INTEGER                                                           &
     & TRIVDIR(IMTR,JMTR)                                               &
                           ! TRIP River Direction file
     &,TRIVSEQ(IMTR,JMTR)  ! TRIP River Sequence file

      REAL                                                              &
     & TRIP_OUTFLOW(IMTR,JMTR)                                          &
! add local variables for inland basins
     &,TRIP_INLANDOUT(IMTR,JMTR)                                        &
                                                     !TRIP OUTFLOW
!                  FROM  INLAND BASINS ON NON_INVERTED TRIP
!                                       GRID KG/S

     &,RIVEROUT_TRIP_MOUTH(IMTR, JMTR)                                  &
                                       ! gridbox outflow on river grid
!                                ! (Kg/s)
     &,ADJUSTMENTT(IMTR,JMTR)                                           &
                              ! RATIO/DELTA stored for 2nd do_areaver
                             ! call. Assumes atm size:ICOLS & JROWS
     &,DUMMYA(ICOLS,JROWS)                                              &
     &,RUNOFFIN_INV(ICOLS,JROWS)                                        &
                                  ! Inverted runoff
     &,YVA_INV(0:JROWS)                                                 &
                                ! INVERTED Atmosphere UV lat coords
     &,RIVEROUT(ICOLS,JROWS)                                            &
                            ! river flow out from each gridbox
!                           ! (KG/S)
! add local variable for inland basins
     &,INLANDOUT(ICOLS,JROWS)                                           &
                                         !OUT TRIP OUTFLOW
!                              FROM INLAND BASINS ON atmos GRID kg/s


     &,RUNOFFIN(ICOLS,JROWS) ! TOTAL RATE OF RUNOFF (KG/M2/S=mm/s)
      REAL                                                              &
     & runofftrip(imtr,jmtr)                                            &
                                  ! regridded runoffin
     & ,runofftrip_inv(imtr,jmtr)                                       &
                                   ! regridded runoff (N-S inverted)
     &,worklat                   ! work for

      INTEGER                                                           &
     & jmax                      ! max. no. of j values for TRIP

      REAL                                                              &
     & XTT(IMTR)                                                        &
                          ! TRIP TS longitude coordinates
     &,YTT(JMTR)                                                        &
                          ! TRIP TS latitude coordinates
     &,XUT(0:IMTR)                                                      &
                          ! TRIP UV longitude coordinates
     &,YUT(0:JMTR)        ! TRIP UV latitude coordinates

! Variables used for regridding runoff
      REAL                                                              &
     & ATLAMBDA(ICOLS*JROWS),                                           &
                                ! LONGITUDE COORDS OF COLUMNS IN
                                !  ATMOSPHERE GRID, IN DEGREES.
     & ATPHI(ICOLS*JROWS),                                              &
                                ! LATITUDE COORDS OF ROWS IN ATMOSC
                                !   GRID, IN DEGREES.
     & TRLAMBDA(IMTR*JMTR),                                             &
                                ! LONG OF SEA POINTS ON OCEAN GRID
     & TRPHI(IMTR*JMTR),                                                &
                                ! LAT OF SEA POINTS ON OCEAN GRID
     & WEIGHTTR(ICOLS,JROWS),                                           & 
                                ! WEIGHTS OF 'TOP RIGHT' CORNERS
     & WEIGHTTL(ICOLS,JROWS),                                           & 
                                ! WEIGHTS OF 'TOP LEFT' CORNERS
     & WEIGHTBR(ICOLS,JROWS),                                           & 
                                ! WEIGHTS OF 'BOTTOM RIGHT' CORNERS
     & WEIGHTBL(ICOLS,JROWS)    ! WEIGHTS OF 'BOTTOM LEFT' CORNERS
!
      INTEGER                                                           &
     & INDEXA(IMTR*JMTR),                                               &
                                ! INDEX OF CORRESPONDING atmos POINTS.
     & INDEXTR(IMTR*JMTR),                                              &
                                ! INDEX OF CORRESPONDING atmos POINTS
     & NCOASTAL                                                         &
                                ! NUMBER OF COASTAL POINTS
     &,TRPOINTS                                                         &
                                ! NUMBER OF POINTS IN TRIP GRID
     &,TRLANDPOINTS                                                     &
                                ! No. land pts on TRIP grid
     &,TRPOINT(IMTR*JMTR)                                               &
                                ! LIST OF LAND POINTS IN TRIP GRID
     &,ATPOINTS                                                         &
                                ! NUMBER OF POINTS IN atmos GRID
     &,ATPOINT(ICOLS*JROWS)                                             &
                                ! LIST OF LAND POINTS IN atmos GRID
     &,SEAPOINTS                !seapoints on the atmos grid
!
!     THE NEXT FOUR VARIABLES ARE NECESSARY TO SATISFY THE ARGUMENT
!     LIST OF COAST_AJ. THEY ARE NOT USED FOR ANYTHING IN THIS ROUTINE.
!
      INTEGER                                                           &
     & IDUMMY1(ICOLS,JROWS),                                            &
                                 ! DUMMY ARGUMENT FOR CALL TO COAST_AJ.
     & IDUMMY2(ICOLS,JROWS),                                            &
                                 ! DUMMY ARGUMENT FOR CALL TO COAST_AJ.
     & N1,                                                              &
                                 ! DUMMY ARGUMENT FOR CALL TO COAST_AJ.
     & N2                        ! DUMMY ARGUMENT FOR CALL TO COAST_AJ.
!
      INTEGER                                                           &
     & INDEXBL(ICOLS,JROWS),                                            &
                                 ! GATHER INDICES FOR INTERPOLATION.
     & INDEXBR(ICOLS,JROWS)      ! GATHER INDICES FOR INTERPOLATION.

!
!     Local arrays required for area-averaging
!
      integer                                                           &
     & lenl                      ! length of lists on the target grid

!      parameter(maxl=(IMTR+ICOLS)*(JMTR+JROWS))
      integer                                                           &
     & index_arav((IMTR+ICOLS)*(JMTR+JROWS))                            &
                                               ! list of source boxes
     &,index_arav1((IMTR+ICOLS)*(JMTR+JROWS))                           &
     &,index_back((IMTR+ICOLS)*(JMTR+JROWS))

      real                                                              &
     & weight((IMTR+ICOLS)*(JMTR+JROWS))                                &
                                               ! weights for source b
     &,weight1((IMTR+ICOLS)*(JMTR+JROWS))                               &
     &,weight2((IMTR+ICOLS)*(JMTR+JROWS))                               &
              ! weights used in o2a coastal flux regridding
     &,backweight((IMTR+ICOLS)*(JMTR+JROWS))                            &
     &,yut_inv(jmtr+1)                                                  &
                                   ! TRIP latitudes in decr. order
!                                  ! if invert_TRIP=.TRUE.
     &,ytt_inv(jmtr)               ! TRIP latitudes in decr. order
!                                  ! if invert_TRIP=.TRUE.
      logical                                                           &
     & trmask_inv(imtr,jmtr)       ! TRIP mask with lats in decr. order
!                                  ! if invert_TRIP=.TRUE.
      integer                                                           &
     & count_tr(imtr*jmtr)                                              &
                                   ! number of A boxes per TRIP box
! add variable for inland basins regridding
     &,COUNT_TARGET(IMTR,JMTR)                                          &
                                                     !COUNT_TR
!                                             DIMENSIONED ON  IMTR,JMTR

     &,base_tr(imtr*jmtr)          ! first index in A box list
      integer                                                           &
     & count_a(icols*jrows)                                             &
                                   ! number of TRIP boxes per A box
     &,base_a(icols*jrows)         ! first index in TRIP box list
!

       REAL VALMAX,VALMIN             ! Field A==>T test variables
       REAL WRKMAX,WRKMIN
       PARAMETER(WRKMAX=-1.0e10,WRKMIN=1.0e10)

! initialise inland basin  variables
       INLANDOUT_ATM=0.0
       inlandout_trip=0.0
       trip_inlandout=0.0
       inlandout=0.0


! Sum the surface and subsurface runoffs
! Sum the two types of runoff
      DO J = 1, JROWS
        DO I = 1, ICOLS
          IF(SURF_RUNOFFIN(i,j) /= RMDI                                 &
     &       .and.SUB_RUNOFFIN(i,j) /= rmdi)THEN
            RUNOFFIN(i,j) = (SURF_RUNOFFIN(i,j)                         &
     &                          + SUB_RUNOFFIN(i,j))*FRAC(i,j)
          ELSE
            RUNOFFIN(i,j) = RMDI
          ENDIF
        ENDDO
      ENDDO
!******************************************************************
! Change river direction and sequence arrays to integer as expected by
! river routing scheme
       DO J=1, JMTR
         DO I=1,IMTR
           TRIVDIR(I,J) = INT(R_RIVDIR(I,J))
           TRIVSEQ(I,J) = INT(R_RIVSEQ(I,J))
         ENDDO
       ENDDO

!******************************************************************
!    SECTION 1: PREPARATIONS FOR atmos to TRIP GRIDS
!
! Set up 2 logical land/sea masks.

      DO J=1,JROWS
        DO I=1,ICOLS
! .TRUE. for gridboxes with any land
          AMASK(I,J) = (FRAC(I,J) >  0.0)
! .FALSE. for gridboxes with any sea
          AMASK_SEA(I,J) = (FRAC(I,J) == 1.0)
        ENDDO
      ENDDO


!     USE THE RIVER DIRECTION FILE TO IDENTIFY THE LAND PTS ON THE
!    set up the TRIP grid from the river direction file
! (0.0 for seapts)

         DO J=1,JMTR
             DO I=1,IMTR
               TRMASK(I,J) = TRIVDIR(I,J) /= RMDI
             ENDDO
         ENDDO

! Set the number of river routing timestes per day to 1
      ndev = 1

      write (6, *) 'number of routing timesteps per day = ', ndev
      write (6, *) 'routing timestep =  ', DT, ' secs.'
      write (6, *) 'Meandering Ratio  = ', ratmed
      write (6, *) 'Effective Velocity ru =', ru, ' (m/s)'

! Southern boundary latitude for TRIP routines
      worklat = int ((latbdy + 90.0) * real(jmtr) / 180.0 + 0.5)

      IF(worklat == 0)THEN
      jmax=jmtr
!      write (6, *) '90N -- ',latbdy,' jmax set to ', jmax
      ELSE
      jmax = jmtr - worklat + 1
!      write (6, *) '90N -- ',latbdy,' jmax=', jmax
      ENDIF


!-------------------------------------------------------------------
! Regrid the runoff if required
       IF(REGRID.AND.(IMTR == ICOLS.AND.JMTR == JROWS))THEN
         CMESSAGE =' RIVER1A : river row length and river rows wrong'
         ICODE = 1
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ICODE,Cmessage)
       ENDIF

       IF(REGRID)THEN

           AMASK_ALL = .FALSE.
           TRMASK_ALL = .FALSE.

! Set up the TRIP coords
       DO J = 0,JMTR
         YUT(J)= 90. - J
         IF(J /= 0)YTT(J)= 90. + 0.5 - J
       ENDDO

       DO J = 1,JMTR+1
         IF(INVERT_ATMOS)THEN
           YUT_INV(J)=YUT(JMTR+1-J)
           IF(J /= JMTR+1)YTT_INV(J)=YTT(JMTR+1-J)
         ELSE
           YUT_INV(J)=YUT(J-1)
           IF(J /= JMTR+1)YTT_INV(J)=YTT(J)
         ENDIF
       ENDDO
       DO I = 0,IMTR
         XUT(I)= 0.0 + I
         IF(I /= 0)XTT(I)= 0.0 - 0.5 + I
       ENDDO
!
! Check for any negative runoff values
       DO J=1,JROWS
         DO I=1,ICOLS
           IF(RUNOFFIN(I,J) <  0.0)THEN
! Check if any values are negative (but not RMDI)
             IF(RUNOFFIN(I,J) /= RMDI)THEN
               WRITE(6,*)' In RIVER1A Negative RUNOFFIN(',I,',',J,')= ' &
     &          ,RUNOFFIN(I,J)
             ENDIF
           ENDIF
         ENDDO
       ENDDO
!

       DO J=1,JMTR
         DO I=1,IMTR
           RUNOFFTRIP_INV(I,J)=RMDI
         ENDDO
       ENDDO


! Regrid atmos runoff to River routing grid by Area averaging as it is
! sufficient to correctly regrid to ALL TRIP gridpoints (any runoff
! 'appearing' in the TRIP sea is added in as outflow to conserve runoff


      lenl=(IMTR+ICOLS)*(JMTR+JROWS)
!
! DEPENDS ON: pre_areaver
      call pre_areaver(icols,xua,jrows,yva,.true.,icols                 &
     &,.false.,amask_all,imtr,xut,jmtr,yut_inv,.true.,.true.            &
     &,lenl,count_tr,base_tr,index_arav,weight,weight2,icode,cmessage)

! DEPENDS ON: do_areaver
      call do_areaver(icols,jrows,icols,.false.,runoffin                &
     &,imtr,jmtr                                                        &
     &,count_tr,base_tr,imtr,.false.,trmask_all,index_arav,weight,0     &
     &,runofftrip_INV,adjustmentt,icode,cmessage)


!-------------------------------------------------------------------
      ELSE                        ! TRIP same resolution as atmos
        DO J=1,JMTR
          DO I=1,IMTR
            RUNOFFTRIP_INV(I,J)=RUNOFFIN(I,J)
          ENDDO
        ENDDO
      ENDIF                       ! REGRID



! Call TRIP routines to route runoff, first convert runoff from mm/s to
! mm/day .

      DO J=1,JMTR
        DO I=1,IMTR
! Change the runoff from Kg/m2/sec (mm/sec) to mm/day for TRIP routines
! Allow Runoff to be at TRIP seapoints for conserve.

          IF(RUNOFFTRIP_INV(I,J) >  RMDI)THEN
             RUNOFFTRIP_INV(I,J)=RUNOFFTRIP_INV(I,J)*86400
          ELSE
             RUNOFFTRIP_INV(I,J)=0.0
          ENDIF
        ENDDO
      ENDDO
! Then invert the runoff values to the 'normal' TRIP grid
      IF(INVERT_ATMOS)THEN
        do j=1,jmtr
          do i=1,imtr
            runofftrip(i,j)=runofftrip_INV(i,jmtr+1-j)
          enddo
        enddo
      ELSE
        do j=1,jmtr
          do i=1,imtr
            runofftrip(i,j)=runofftrip_INV(i,j)
          enddo
        enddo
      ENDIF
!*****************************************************************

! DEPENDS ON: routedbl
      CALL ROUTEDBL(RUNOFFTRIP, RU, RATMED, NDEV, DT                    &
     &,      IMTR, JMTR, TRIVDIR, TRIVSEQ, TWATSTOR, JMAX, RMDI         &
     &,      RIVEROUT_TRIP,RUNOFF_TRIP_KGS)


! Invert the TRIP OUTFLOW for 'mapping' and remove TRIP_RMDI values

!AJW fix bug where water is only regridded for ocean pour points thus excluding inland
!basins. TRIVDIR(I,J) now >= 9 rather than ==9.

! Use only outflow at river mouths and sea for the mapping onto atmos
! grid
       DO J=1,JMTR
         DO I=1,IMTR
           IF(TRIVDIR(I,J) == RMDI.or.                                  &
     &         TRIVDIR(I,J) >= 9.or.TRIVDIR(I,J) == 0)THEN
             RIVEROUT_TRIP_MOUTH(I,J) = RIVEROUT_TRIP(I,J)
           ELSE
             RIVEROUT_TRIP_MOUTH(I,J) = 0.0
           ENDIF
         ENDDO
       ENDDO

! Invert the TRIP OUTFLOW for 'mapping' and remove TRIP_RMDI values

       DO J=1,JMTR
         DO I=1,IMTR
! For DO_MAP_MAX, do not allow RMDITRIP values
           IF(RIVEROUT_TRIP_MOUTH(I,J) /= RMDI)THEN
             IF(INVERT_ATMOS)THEN
               TRIP_OUTFLOW(I,JMTR+1-J)=RIVEROUT_TRIP_MOUTH(I,J)
             ELSE
               TRIP_OUTFLOW(I,J)=RIVEROUT_TRIP_MOUTH(I,J)
             ENDIF
           ELSE
             IF(INVERT_ATMOS)THEN
               TRIP_OUTFLOW(I,JMTR+1-J)=0.0
             ELSE
               TRIP_OUTFLOW(I,J)=0.0
             ENDIF
           ENDIF
         ENDDO
       ENDDO

! Map River outflow directly to atmos using do_max_map
!
!*******************************************************

       IF (REGRID) THEN

! Prepare for mapping and count river mouths on TRIP grid as sea.
! Instead of regridding, use count_tr, base_tr and weight to 'map' the
! TRIP sea points onto the atmos grid where the major part of the TRIP
! box lies using do_map_max

       IF(INVERT_ATMOS)THEN
        DO J=1,JMTR
         DO I=1,IMTR
           TRMASK_INV(I,J)=(RIVEROUT_TRIP(I,JMTR+1-J) == RMDI)
         ENDDO
        ENDDO
       ELSE
        DO J=1,JMTR
         DO I=1,IMTR
           TRMASK_INV(I,J)=(RIVEROUT_TRIP(I,J) == RMDI)
         ENDDO
        ENDDO
       ENDIF
! Initialise Riverout
       RIVEROUT = 0.0
! First Regrid by AVER_TAO
        lenl=(IMTR+ICOLS)*(JMTR+JROWS)
! DEPENDS ON: pre_areaver
        call pre_areaver(icols,xua,jrows,yva,.true.,icols               &
     &,.false.,amask_sea,imtr,xut,jmtr,yut_inv,.true.,.true.            &
     &,lenl,count_tr,base_tr,index_arav,weight,weight2,icode,cmessage)

! DEPENDS ON: do_map_max
        call DO_MAP_MAX(icols,jrows,icols,.false.,RIVEROUT              &
     &,imtr,jmtr                                                        &
     &,count_tr,base_tr,imtr,.false.,trmask_inv,index_arav,weight,0     &
!  add variables to do_map_max call for inland basin identification
     &,TRIP_OUTFLOW,ADJUSTMENTT,ICODE,CMESSAGE                          &
     &,COUNT_TARGET)



!********************************************************************
! Inland basin now 10
!  calculate  inland basin outflow using grid points where
!  count_target=0, trivdir=10, rmdi or 0, and count_target=0

         DO J=1,JMTR
         DO I=1,IMTR

         IF(INVERT_ATMOS) THEN

            IF((TRIVDIR(I,JMTR+1-J) >= 9.OR.TRIVDIR(I,JMTR+1-J)        &
     &       == RMDI.OR.TRIVDIR(I,JMTR+1-J) == 0)                      &
     &       .AND.COUNT_TARGET(I,J) == 0) THEN

                TRIP_INLANDOUT(I,J)=TRIP_OUTFLOW(I,J)
                INLANDOUT_TRIP(I,JMTR+1-J)=TRIP_OUTFLOW(I,J)

            ELSE

                TRIP_INLANDOUT(I,J)=0.0
                INLANDOUT_TRIP(I,JMTR+1-J)=0.0

            ENDIF

         ELSE


            IF((TRIVDIR(I,J) >= 9.OR.TRIVDIR(I,J) == RMDI.OR.           &
     &      TRIVDIR(I,J) == 0).AND.COUNT_TARGET(I,J) == 0) THEN

                TRIP_INLANDOUT(I,J)=TRIP_OUTFLOW(I,J)
                INLANDOUT_TRIP(I,J)=TRIP_OUTFLOW(I,J)

             ELSE

                TRIP_INLANDOUT(I,J)=0.0
                INLANDOUT_TRIP(I,J)=0.0

            ENDIF

         ENDIF

         ENDDO
         ENDDO


! call pre_areaver to regrid inland basin outflow
        lenl=(IMTR+ICOLS)*(JMTR+JROWS)
! DEPENDS ON: pre_areaver
        call PRE_AREAVER(icols,xua,jrows,yva,.true.,icols               &
     &,.false.,amask_all,imtr,xut,jmtr,yut_inv,.true.,.true.            &
     &,lenl,count_tr,base_tr,index_arav,weight,weight2,icode,cmessage)


! send inland basin outputs to inlandout/inlandout_atm
! set trmask_inv check to false

! DEPENDS ON: do_map_max
        call DO_MAP_MAX(icols,jrows,icols,.true.,INLANDOUT              &
     &,imtr,jmtr                                                        &
     &,count_tr,base_tr,imtr,.false.,trmask_all,index_arav,weight,0     &
     &,INLANDOUT_TRIP,adjustmentt,icode,cmessage,count_target)


! output inland basin outflow for case where grids are the same


      ELSE                        ! TRIP same resolution as atmos
        DO J=1,JROWS
          DO I=1,ICOLS
! Save only seapoints,
            IF(AMASK_sea(I,J).eqv..FALSE.)THEN
              RIVEROUT(I,J)=TRIP_OUTFLOW(I,J)
! set rmdi values to zero since hydrol7a cannot use them
      if(trip_outflow(i,j) == rmdi)then
         inlandout_atm(i,j)=0.0
      else
         INLANDOUT_ATM(I,J)=TRIP_OUTFLOW(I,J)
      endif

            ELSE
              RIVEROUT(I,J)=0.0
      INLANDOUT_ATM(I,J)=0.0

            ENDIF
          ENDDO
        ENDDO
      ENDIF                       ! REGRID

!******************************************************************

! Change river outflow to mm/s (as expected by the ocean later)

      DO J=1,JROWS
        DO I=1,ICOLS
          IF(RIVEROUT(I,J) /= RMDI)THEN
              RIVEROUT_INV(I,J) =                                       &
     &          RIVEROUT(I,J)/GRIDBOX_AREAS(I,J)
          ELSE
              RIVEROUT_INV(I,J) = RMDI
          ENDIF
        ENDDO
      ENDDO


! convert inland basin outflow from kg/s to kg/m2/s
! as expected by hydrol7a

       DO J=1,JROWS
         DO I=1,ICOLS

          IF(INLANDOUT(I,J) /= RMDI)THEN
             INLANDOUT_ATM(I,J)=                                        &
     &       INLANDOUT(I,J)/GRIDBOX_AREAS(I,J)
          ELSE

!      set to zero, not rmdi since hydrol7a cannot
!      use rmdi values

             INLANDOUT_ATM(I,J)=0.0
          ENDIF

         ENDDO
       ENDDO

! Set land points values of inlandout_trip to rmdi post regridding

       DO J=1,JMTR
         DO I=1,IMTR

            IF(INVERT_ATMOS) THEN

             IF(RIVEROUT_TRIP(I,JMTR+1-J) == RMDI)                      &
     &        INLANDOUT_TRIP(I,JMTR+1-J)=RMDI
            ELSE

             IF(RIVEROUT_TRIP(I,J) == RMDI)                             &
     &        INLANDOUT_TRIP(I,J)=RMDI
            ENDIF

         ENDDO
       ENDDO


! Set RMDI values in river outflow and runoff on river grid to 0.0 as
! points with RMDI are not constant as runoff greater than 0.0 is
! 'added in' at some TRIP seapoints due to regridding from
! non-coincident grids.
! Otherwise postprocessing of fields can give negative values.

      DO J=1,JMTR
       DO I=1,IMTR
        IF( RIVEROUT_TRIP(I,J) == RMDI)RIVEROUT_TRIP(I,J)=0.0
       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE RIV_ROUT

!*********************************************************************

