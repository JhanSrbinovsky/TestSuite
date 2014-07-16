#if defined(C80_1A) || defined(UTILIO) || defined(FLDOP)               \
 || defined(VAROPSVER)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE PR_LFLD----------------------------------------
!LL
!LL  Purpose: Prints out selected values from logical data
!LL           using information from associated PP header.
!LL
!LL  Written by D. Robinson
!LL
!LL  Model            Modification history:
!LL version  date
!LL   3.3  22/11/93  New routine (adapted from deck PRIFLD1A)
!LL   5.1  31/03/00  Added comma missing from write format. D.P.Matthews
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!LL
!LL  System component: R30/W30
!LL
!LL  System task: F3
!LL
!LL  Programming standard:
!LL           Unified Model Documentation Paper No 3
!LL           Version No 1 15/1/90
!LL
!LL  Documentation:
!LL           Unified Model Documentation Paper No F3
!LL           Version No 5 9/2/90
!LL
!LL------------------------------------------------------------
!*L Arguments:-------------------------------------------------
      SUBROUTINE PR_LFLD(LOOKUP,RLOOKUP,LEN1_LOOKUP,LD1,K)

      IMPLICIT NONE

      INTEGER                                                           &
     & K                                                                &
                     !IN Field number ie position in 2nd dim
                     !   of LOOKUP
     &,LEN1_LOOKUP                                                      &
                     !IN First dimension of LOOKUP table
     &,LOOKUP(LEN1_LOOKUP,*)  !IN Integer equivalence of PP LOOKUP

      REAL                                                              &
     & RLOOKUP(LEN1_LOOKUP,*) !IN Real equivalence of PP LOOKUP

      LOGICAL                                                           &
     & LD1(*)        !IN Kth field in data array
! -------------------------------------------------------------
!*L External subroutines called:-------------------------------
! None
!--------------------------------------------------------------
!*L Local control constants:-----------------------------------
      INTEGER                                                           &
     & NS_PTS                                                           &
                     !PARAM No of points down to print
     &,EW_PTS        !PARAM No of points across to print
      PARAMETER(NS_PTS=6,EW_PTS=5)
! -------------------------------------------------------------
! Workspace usage:---------------------------------------------
      REAL LON(EW_PTS)     ! Longitudes printed out
      INTEGER I(EW_PTS)    ! Index of values printed out
      CHARACTER*12 DASH(EW_PTS)  !Stores dashed lines
!*-------------------------------------------------------------
! Local variables:---------------------------------------------
      INTEGER                                                           &
     & N_ROWS                                                           &
                   ! No of rows in field
     &,N_COLS                                                           &
                   ! No of colums in field
     &,ROW                                                              &
                   ! Row number
     &,R_INC,F_INC                                                      &
                   ! No of rows/points between printed lines
     &,J,L                                                              &
                   ! Loop counts
     &,EW_PRINT                                                         &
                   ! No of E-W values printed out
     &,POS_MIN                                                          &
                   ! Position of Minimum value of field
     &,POS_MAX                                                          &
                   ! Position of Maximum value of field
     &,F_MIN                                                            &
                   ! Minimum value of field
     &,F_MAX       ! Maximum value of field

      REAL                                                              &
     & LAT         ! Latitude
!--------------------------------------------------------------

#include "clookadd.h"
#include "c_mdi.h"

!L Internal structure: None

! Initialise string used to create table boundaries
      DO J=1,EW_PTS
        DASH(J)='------------'
      ENDDO

      IF(LOOKUP(LBCODE,K) == IMDI) THEN
!       IF LBCODE IS MISSING DATA, ASSUME THAT THE FIELD IN DUMP
!       HAS NOT BEEN WRITTEN TO BY STASH.
!       THIS SHOULD ONLY OCCUR TO DIAGNOSTIC PARTS OF THE DUMP BEFORE
!       FIRST WRITE BY STASH TO THAT AREA/HEADER.
        WRITE(6,*) 'MESSAGE FROM PR_LFLD'
        WRITE(6,*) 'LBCODE NOT SET; ASSUME DATA NOT SET. NO PRINT'
        RETURN
      END IF

! No of rows and columns in field
      N_ROWS=LOOKUP(LBROW,K)
      N_COLS=LOOKUP(LBNPT,K)


      IF(N_COLS /= 0.AND.N_COLS /= IMDI)THEN

! No of E-W values to be printed
      EW_PRINT=MIN(N_COLS,EW_PTS)

! Calculate longitudes and addresses of values to be printed from 1st ro
      I(1)=1
      LON(1)=RLOOKUP(BZX,K)+RLOOKUP(BDX,K)
      DO 100 J=1,EW_PTS-2
      I(J+1)=I(J)+N_COLS/(EW_PTS-1)
      LON(J+1)=LON(J)+RLOOKUP(BDX,K)*(N_COLS/(EW_PTS-1))
100   CONTINUE
      I(EW_PTS)=N_COLS
      LON(EW_PTS)=RLOOKUP(BZX,K)+RLOOKUP(BDX,K)*N_COLS

! Initialise row and field pointers
      ROW=1
      LAT=RLOOKUP(BZY,K)+RLOOKUP(BDY,K)
      R_INC=N_ROWS/(NS_PTS-1)
      F_INC=R_INC*N_COLS

! Print 1st row
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)
      WRITE(6,'('' FIELD NO'',I4,'':'',9(F10.3,2X))')                   &
     &K,(LON(J),J=1,EW_PRINT)
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)

! Print remaining rows except last
      DO 200 L=1,NS_PTS-1
      WRITE(6,'(1X,I3,'':'',F8.3,'':'',3X,9(L9,3X))')ROW,LAT,           &
     &(LD1(I(J)),J=1,EW_PRINT)
      DO 300 J=1,EW_PTS
      I(J)=I(J)+F_INC
300   CONTINUE
      ROW=ROW+R_INC
      LAT=LAT+R_INC*RLOOKUP(BDY,K)
200   CONTINUE

! Calculate addresses used to print values for last row
      I(1)=1+(N_ROWS-1)*N_COLS
      DO 400 J=1,EW_PTS-2
      I(J+1)=I(J)+N_COLS/(EW_PTS-1)
400   CONTINUE
      I(EW_PTS)=N_ROWS*N_COLS

! Set row pointers to last row
      LAT=RLOOKUP(BZY,K)+RLOOKUP(BDY,K)*N_ROWS
      ROW=N_ROWS

! Print last row
      WRITE(6,'(1X,I3,'':'',F8.3,'':'',3X,9(L9,3X))')ROW,LAT,           &
     &(LD1(I(J)),J=1,EW_PRINT)
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)
      ELSE

! Print out summary of non standard fields

      EW_PRINT=MIN(EW_PTS,LOOKUP(LBLREC,K))
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)
      WRITE(6,'('' FIELD NO'',I4,'':  DATA NOT ON MODEL GRID''          &
     &,'' SO FIRST FEW VALUES PRINTED'')')K
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)
      WRITE(6,'(1X,3X,'':'',8X,'':'',3X,9(L9,3X))')                     &
     &(LD1(J),J=1,EW_PRINT)
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)

      ENDIF

      WRITE(6,'('' '')')

      RETURN
      END SUBROUTINE PR_LFLD

#endif
