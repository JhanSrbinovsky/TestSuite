#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE BOX_BND
!LL
!LL  Routine sets up arrays of indexes and boundaries of target grid
!LL  boxes relative to a source grid for use  by BOX_SUM so that
!LL  area weighted means can be calculated.
!LL
!LL  NOT SUITABLE FOR SINGLE COLUMN USE
!LL
!LL  SUITABLE FOR ROTATED GRIDS
!LL
!LL  ORIGINAL VERSION FOR CRAY Y-MP/IBM
!LL  WRITTEN 06/09/91 BY C. WILSON
!LL
!LL  CODE REVIEWED BY R.SMITH ??/??/??
!LL
!LL  VERSION NO. 1 DATED 06/09/91
!LL         COSMOS DSN MS15.CWUM.JOBS(BOXBND1)
!LL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
!LL  VERSION 1, DATED 12/09/89
! History:
! Version   Date     Comment
! -------   ----     -------
! 4.0      12/04/95  Imported into Unified model. D.M. Goddard
! 4.1      12/06/96  Corrections to longitude indexes at boundary.
!                    D.M. Goddard
! 4.5      12/10/98  Stops LONG_L and COLAT_T becoming negative in
!                    top row, preventing out of bound in index array.
!                    Author D.M Goddard
!LL  5.1  10/04/00   New Reconfiguration & S->N ordering support.
!LL                  P.Selwood.
!LL
!LL  SYSTEM TASK:  S1 (part,extension for area mean interpolation)
!LL
!LL  PURPOSE: To set up for grid-boxes on a target grid the longitude,
!LL           colatitude,and indexes of the overlapping source grid-
!LL           boxes at left hand side and top of target grid-boxes.
!LL           Both grids are regular lat-long with the same pole
!LL           and orientation. Either may be a 'p-grid' (ie with
!LL           half-size boxes at the poles) or a 'u-grid' (ie with
!LL           regular size boxes everywhere.)
!LL
!LL           NB The units used are "source grid-box lengths"
!LL           The area of the target grid_boxes in squared "source
!LL           grid-box units" is also returned.
!LL
!LL  DOCUMENTATION:  UNIFIED MODEL DOCUMENTATION S1
!LL                  BY A.DICKINSON/C WILSON VERSION ??DATED ??/??/91
!LL
!LL
!LLEND-------------------------------------------------------------

!
!*L  ARGUMENTS:---------------------------------------------------
      SUBROUTINE BOX_BND                                                &
     &  (I_L,LONG_L,J_T,COLAT_T,AREA_BOX,                               &
     &   ROW_LENGTH,ROWS,ROW_LENGTH_SRCE,ROWS_SRCE,                     &
     &   DELTA_LONG,DELTA_LAT,START_LONG,START_LAT,                     &
     &   DELTA_LONG_SRCE,DELTA_LAT_SRCE,START_LONG_SRCE,START_LAT_SRCE, &
     &   IGRID,IGRID_SRCE,GLOBAL)


#if defined(RECON)
      Use Rcf_PrintStatus_Mod, Only :                                   &
     &    PrintStatus,                                                  &
     &    PrStatus_Diag
#endif

      IMPLICIT NONE

      INTEGER                                                           &
     &  ROW_LENGTH                                                      &
                           !IN    Number of points per row target area
     &, ROWS                                                            &
                           !IN    Number of rows of target area
     &, ROW_LENGTH_SRCE                                                 &
                           !IN    Number of points per row source area
     &, ROWS_SRCE          !IN    Number of rows of source area
!
      INTEGER                                                           &
     & I_L(ROW_LENGTH+1)                                                &
                        !OUT Index of source box overlapping lhs of
                        ! target grid-box
     &,J_T(ROWS+1)      !OUT Index of source box overlapping top of
                        ! target grid-box
!
      REAL                                                              &
     & LONG_L(ROW_LENGTH +1)                                            &
                             !OUT Left longitude of target grid-box (in
                             ! units of DELTA_LONG_SRCE)
     &,COLAT_T(ROWS+1)                                                  &
                             !OUT Colatitude of top of target grid-box
                             ! (in units of DELTA_LAT_SRCE)
     &,  AREA_BOX !OUT area of grid box in sq units of source grid

      REAL                                                              &
     & DELTA_LONG                                                       &
                       !IN   Longitude increment of target grid (deg)
     &,DELTA_LAT                                                        &
                       !IN   Latitude increment of target grid (deg)
     &,START_LONG                                                       &
                       !IN   start longitude of centre of first grid-
                       !     box in target area
     &,START_LAT                                                        &
                       !IN   start latitude of centre of first grid-
                       !     box in target area
     &,DELTA_LAT_SRCE                                                   &
                       !IN   Latitude increment of source grid
     &,DELTA_LONG_SRCE                                                  &
                       !IN   Longitude increment of source grid
     &,START_LONG_SRCE                                                  &
                       !IN   start longitude of centre of first grid-
                       !     box in source area
     &,START_LAT_SRCE  !IN   start latitude of centre of first grid
                       !     box in source area

      INTEGER                                                           &
     & IGRID                                                            &
                       !IN   Grid indicator 1=p-grid,2=u-grid
     &,IGRID_SRCE      !IN   Grid indicator 1=p-grid,2=u-grid

      LOGICAL GLOBAL   !IN    true if global area required
!*---------------------------------------------------------------------
#if !defined(RECON)
#include "cprintst.h"
#endif

!*L  WORKSPACE USAGE:-------------------------------------------------
!   NONE
!*---------------------------------------------------------------------
!
!*L EXTERNAL SUBROUTINES CALLED---------------------------------------
!     EXTERNAL NONE
!*------------------------------------------------------------------
!----------------------------------------------------------------------
!    DEFINE LOCAL VARIABLES
      REAL EW_BOX                                                       &
                  ! length of grid box in units of source
                             ! (DELTA_LONG_SRCE)
     &,    NS_BOX                                                       &
                  ! height of grid box in units of source
                             ! (DELTA_LAT_SRCE)
     &,START_LONG_BOX                                                   &
                      ! start longitude of first grid box left edge
     &,START_COLAT_BOX                                                  &
                      ! start colatitude of first grid box top edge
     &,START_LONG_BOX_SRCE                                              &
                           ! start longitude of first grid box left
                           ! edge ( source)
     &,START_COLAT_BOX_SRCE                                             &
                           ! start colatitude of first grid box top
                           ! edge ( source)
     &,LONG_OFFSET                                                      &
                           ! start longitude difference
     &,COLAT_OFFSET        ! start colatitude difference

      INTEGER I,J ! loop counters
      REAL P1,P2
      LOGICAL LNER
      LNER(P1,P2) = ((ABS(P1-P2))  >   (1.E-5*ABS(P1+P2)))

!L *********************************************************************
!L 1.0 Set target gridbox length,height and area in units of source grid
!L     Set start 'longitude' and 'colatitude' of first grid box
!L     also in source units.
!L *********************************************************************

      EW_BOX= DELTA_LONG/DELTA_LONG_SRCE
      NS_BOX= DELTA_LAT/DELTA_LAT_SRCE
      AREA_BOX= EW_BOX*NS_BOX

!L *********************************************************************
!L 1.1 Set start colatitude of top and start longitude of LHS of first
!L     boxes on both target and source grids.
!L *********************************************************************
      If ( PrintStatus >= PrStatus_Diag ) Then
        Write(6,*) DELTA_LAT_SRCE,DELTA_LAT
        Write(6,*) DELTA_LONG_SRCE,DELTA_LONG
        Write(6,*) START_LAT_SRCE,START_LAT
        Write(6,*) START_LONG_SRCE,START_LONG
      End If
      START_LONG_BOX = START_LONG - 0.5*DELTA_LONG
      START_COLAT_BOX = (-90. - START_LAT) - 0.5*DELTA_LAT
      START_LONG_BOX_SRCE = START_LONG_SRCE - 0.5*DELTA_LONG_SRCE
      START_COLAT_BOX_SRCE = (-90. - START_LAT_SRCE) -0.5*DELTA_LAT_SRCE

!##   START_LONG_BOX = START_LONG/DELTA_LONG_SRCE - 0.5*EW_BOX
!##   START_COLAT_BOX = (90. - START_LAT)/DELTA_LAT_SRCE - 0.5*NS_BOX
!##   START_LONG_BOX_SRCE = START_LONG_SRCE/DELTA_LONG_SRCE - 0.5
!##   START_COLAT_BOX_SRCE = (90. - START_LAT_SRCE)/DELTA_LAT_SRCE - 0.5

      IF(GLOBAL) THEN
       IF(IGRID_SRCE == 1.AND.LNER(START_LAT_SRCE,-90.)) THEN
         WRITE(6,*)' BOX_BND: source grid not global'
         STOP
       ENDIF
       IF(IGRID_SRCE == 2.AND.                                          &
     &   LNER(START_LAT_SRCE,(-90.+DELTA_LAT_SRCE*0.5))) THEN
         WRITE(6,*)' BOX_BND: source grid not global'
       ENDIF
      ELSE
       WRITE (6,*) 'BOX_BND'
       WRITE (6,*) 'start_colat_box      ',START_COLAT_BOX
       WRITE (6,*) 'start_colat_box_srce ',START_COLAT_BOX_SRCE
       IF(LNER(START_COLAT_BOX,START_COLAT_BOX_SRCE)) THEN
         WRITE(6,*)' BOX_BND: target area larger than source area'
         STOP
       ENDIF
      ENDIF

      LONG_OFFSET  =(START_LONG_BOX-START_LONG_BOX_SRCE)/               &
     &                            DELTA_LONG_SRCE
      COLAT_OFFSET = (START_COLAT_BOX - START_COLAT_BOX_SRCE)/          &
     &               DELTA_LAT_SRCE

      IF (.NOT.GLOBAL) THEN
      IF (LONG_OFFSET <  0.0) THEN
        WRITE (6,*) ' LONG_OFFSET = ',LONG_OFFSET,' ; Reset to 0.0'
        LONG_OFFSET = 0.0
      ENDIF
      IF (COLAT_OFFSET <  0.0) THEN
        WRITE (6,*) ' COLAT_OFFSET = ',COLAT_OFFSET,' ; Reset to 0.0'
        COLAT_OFFSET = 0.0
      ENDIF
      ENDIF
!L *********************************************************************
!L 2.0 Set grid box left longitudes, top colatitudes and indices
!L *********************************************************************

      DO 220 I=1,ROW_LENGTH + 1
        LONG_L(I) = LONG_OFFSET + (I-1)*EW_BOX
        IF(GLOBAL.AND.LONG_L(I) <  0.0)THEN
          LONG_L(I)=LONG_L(I)+REAL(ROW_LENGTH_SRCE)
        ELSE IF(GLOBAL.AND.LONG_L(I) >= ROW_LENGTH_SRCE)THEN
          LONG_L(I)=LONG_L(I)-REAL(ROW_LENGTH_SRCE)
        END IF
        I_L(I) = LONG_L(I) +1
 220  CONTINUE
      IF(LONG_L(1) <  0.0)LONG_L(1)=0.0

      If ( PrintStatus >= PrStatus_Diag ) Then
      WRITE(6,*) ' I_L'
      WRITE(6,*) I_L
      WRITE(6,*) ' LONG_L'
      WRITE(6,*) LONG_L
      End If

!##   COLAT_T(1) = START_COLAT_BOX - START_COLAT_BOX_SRCE
      COLAT_T(1) = COLAT_OFFSET
      IF(GLOBAL.AND.IGRID == 1)THEN
!##     COLAT_T(1) = 0.0  - START_COLAT_BOX_SRCE
        COLAT_T(1) = (0.0  - START_COLAT_BOX_SRCE)/DELTA_LAT_SRCE
      ENDIF
        IF(COLAT_T(1) <  0.0)COLAT_T(1)=0.0
      J_T(1) = COLAT_T(1) + 1
      DO 230 J=2,ROWS+1
        COLAT_T(J) = COLAT_OFFSET  + (J-1)*NS_BOX
        J_T(J) = COLAT_T(J) + 1
 230  CONTINUE
!    ROWS+1 ie bottom boundary
      IF(GLOBAL) THEN
        IF(IGRID == 1)COLAT_T(ROWS+1) = COLAT_T(ROWS+1)-0.5*NS_BOX
        J_T(ROWS+1) = ROWS_SRCE
      ELSE
        IF(J_T(ROWS+1) >  ROWS_SRCE) THEN
          IF(COLAT_T(ROWS+1) >  REAL(ROWS_SRCE)) THEN
            WRITE(6,*)' BOX_BND: target area larger than source area'
            STOP
          ELSE
            J_T(ROWS+1) = ROWS_SRCE
          ENDIF
        ENDIF
      ENDIF

      If ( PrintStatus >= PrStatus_Diag ) Then
      WRITE(6,*) ' J_T'
      WRITE(6,*) J_T
      WRITE(6,*) ' COLAT_T'
      WRITE(6,*) COLAT_T
      End If

      RETURN
      END SUBROUTINE BOX_BND
#endif
