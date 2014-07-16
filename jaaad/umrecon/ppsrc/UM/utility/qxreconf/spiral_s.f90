
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL SUBROUTINE SPIRAL_S-----------------------------------------------
!LL
!LL Written by C.P. Jones
!LL
!LL Reviewed by ??
!LL
!LL Programming standard:
!LL    Unified Model Documentation Paper No 3
!LL
!LL System component: S121
!LL
!LL System task: S1
!LL
!LL Purpose:
!LL   Attempts to set a value at points which are unresolved when
!LL   interpolating between one grid and another.  A value is set
!LL   by finding the mean of surrounding points which do have data
!LL   set within a search radius determined by NSEARCH.
!LL
!LL  Modification History:
!LL
!LL  Model
!LL  Version  Date
!LL  4.0      02/11/95  The Fortran for squaring an array needed to be
!LL     redefined for the DecAlpha because of lexcon. (N.Farnon)
!LL  4.4      23/09/97 amendment to calculation of row number and
!LL                    to check for missing data (M. J. Bell)
!LL Documentation:
!LL   UMDP S1
!LL
!LL -------------------------------------------------------------
      SUBROUTINE SPIRAL_S(LAND_SEA_MASK,INDEX_UNRES,NO_POINT_UNRES,     &
     &           POINTS_PHI,POINTS_LAMBDA,DATA_FIELD,NSEARCH,SEA_LAND,  &
     &           CYCLIC)

      IMPLICIT NONE

!*L ARGUMENTS:---------------------------------------------------

      INTEGER                                                           &
     & POINTS_PHI                                                       &
                        !IN number of rows in grid
     &,POINTS_LAMBDA                                                    &
                        !IN number of columns in grid
     &,NSEARCH                                                          &
                        !IN number of points in each direction to search
     &,NO_POINT_UNRES                                                   &
                        !INOUT number of unresolved points
     &,LAND_SEA_MASK(POINTS_LAMBDA*POINTS_PHI)                          &
                        !IN land sea mask
     &,INDEX_UNRES(POINTS_LAMBDA*POINTS_PHI)                            &
                        !INOUT index to unresolved pts
     &,SEA_LAND         !IN =0 for sea field  =1/-1 for land field

      REAL                                                              &
     & DATA_FIELD(POINTS_LAMBDA*POINTS_PHI) !INOUT field

      LOGICAL                                                           &
     & CYCLIC           ! IN =T if data covers complete latitude circle

!*L Parameters
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

!*L LOCAL VARIABLES

      INTEGER                                                           &
     & I,J,JJ,JJJ,K,JK                                                  &
                        ! indices
     &,IPT,IROW,ICOL                                                    &
                      ! coordinate of unresolved point
     &,IPOINT,IUNRES                                                    &
                      ! do loop variable
     &,NPOINTS                                                          &
                      ! number of points in serach box
     &,IR((1+2*NSEARCH)*(1+2*NSEARCH))                                  &
                         ! row numbers of points to serach
     &,IC((1+2*NSEARCH)*(1+2*NSEARCH))                                  &
                         ! col numbers of points to search
     &,IND_SEARCH((1+2*NSEARCH)*(1+2*NSEARCH))                          &
                         ! index to points to search
!LL
     &,NOT_YET_SET                                                      &
                                    ! number of points still to set
     &,IND_YET_SET(POINTS_LAMBDA*POINTS_PHI)                            &
                                             ! index of points
                     ! still unresolved after calling this subroutine
     &,ISUM_MASK     ! number of surrounding points which have data

      REAL                                                              &
     & SUM_DATA                                                         &
                     ! sum of data surrounding unresolved points
     &, RMDI_TOL ! values within this tolerance counted as missing

!*L   EXTERNAL ROUTINES
!     None
!---------------------------------------------------------------------
      RMDI_TOL  = ABS (RMDI) * 0.0001
!

! Calculate number of points in search box
      NPOINTS=(1+2*NSEARCH)**2 ! number of grid points in search box

! Loop around unresolved points
      NOT_YET_SET=0
      DO 10 IUNRES=1,NO_POINT_UNRES

! find unresolved point coordinate in terms of rows and cols
        IPT=INDEX_UNRES(IUNRES)
        IROW= (IPT - 1)/POINTS_LAMBDA   +   1
        ICOL=IPT-(IROW-1)*POINTS_LAMBDA

! calculate surrounding points' coords in terms of rows and cols
        JJJ=1
        DO 20 J=-NSEARCH,NSEARCH
          DO 30 JJ=JJJ,JJJ+2*NSEARCH
            IR(JJ)=IROW+J
 30       CONTINUE
        JJJ=JJJ+1+2*NSEARCH
 20     CONTINUE

        JJJ=1+2*NSEARCH
        JK=1
        DO 40 J=-NSEARCH,NSEARCH
          DO 50 JJ=0,2*NSEARCH
            IC(JK+JJ*JJJ)=ICOL+J
 50       CONTINUE
        JK=JK+1
 40     CONTINUE

! Check that col and rows are in range of grid
        DO 70 IPOINT=1,NPOINTS
          IF(IC(IPOINT) >  POINTS_LAMBDA) THEN
            IF(CYCLIC) THEN
              IC(IPOINT)=IC(IPOINT)-POINTS_LAMBDA
            ELSE
              IC(IPOINT)=POINTS_LAMBDA
            ENDIF
          ENDIF
          IF(IC(IPOINT) <  1) THEN
            IF(CYCLIC) THEN
              IC(IPOINT)=IC(IPOINT)+POINTS_LAMBDA
            ELSE
              IC(IPOINT)=1
            ENDIF
          ENDIF
          IF(IR(IPOINT) <  1) IR(IPOINT)=1
          IF(IR(IPOINT) >  POINTS_PHI) IR(IPOINT)=POINTS_PHI
 70     CONTINUE

! Form index search array
        DO 80 IPOINT=1,NPOINTS
          IND_SEARCH(IPOINT)=(IR(IPOINT)-1)*POINTS_LAMBDA+IC(IPOINT)
 80     CONTINUE

! search for data around this point. If no data is found the point
! remains unresolved

        ISUM_MASK=0   ! number of points with data found
        SUM_DATA=0.0  ! sum of data of surrounding grid points

        DO 90 IPOINT=1,NPOINTS
          IF(IABS(LAND_SEA_MASK(IND_SEARCH(IPOINT))) == IABS(SEA_LAND)  &
     &    .AND.DATA_FIELD(IND_SEARCH(IPOINT))  >   RMDI+RMDI_TOL) THEN
            SUM_DATA=SUM_DATA+DATA_FIELD(IND_SEARCH(IPOINT))
            ISUM_MASK=ISUM_MASK+1
          ENDIF
 90     CONTINUE

       IF(ISUM_MASK  >   0) THEN
! data found - take mean
          DATA_FIELD(IPT)=SUM_DATA/REAL(ISUM_MASK)
        ELSE
! data not found - point remains unresolved
          NOT_YET_SET=NOT_YET_SET+1
          IND_YET_SET(NOT_YET_SET)=IPT
        ENDIF

 10   CONTINUE


! amend output array with points remaining unresolved
      IF(NOT_YET_SET >  0) THEN
        DO 100 IPOINT=1,NOT_YET_SET
          INDEX_UNRES(IPOINT)=IND_YET_SET(IPOINT)
 100    CONTINUE
        NO_POINT_UNRES=NOT_YET_SET
      ELSE
        NO_POINT_UNRES=0
      ENDIF
      RETURN
      END SUBROUTINE SPIRAL_S
