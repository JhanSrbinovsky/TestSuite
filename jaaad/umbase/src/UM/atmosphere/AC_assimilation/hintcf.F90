#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE HINTCF -------------------------------------------------
!LL
!LL  Purpose :
!LL
!LL THIS CALCULATES INTERPOLATION COEFFICIENTS FOR BI-LINEAR
!LL INTERPOLATION OF VALUES FROM MODEL GRID POINTS TO OBSERVATION
!LL POINTS.IT STORES THE COEFFICIENTS IN ANALYSIS WORK ARRAY ANWORK.
!LL IN CALCULATING THE COEFFICIENTS,DISTANCES ARE MEASURED RELATIVE
!LL TO THE NEAREST MODEL POINT TO THE OBSERVATION.INCREMENT VECTORS
!LL (OVER OBSERVATIONS) ARE CALCULATED TO DEFINE THE POSITION
!LL RELATIVE TO THE NEAREST POINT,OF THE FOUR POINTS SURROUNDING
!LL THE OBSERVATION.
!LL
!LL  For use on Cray Y-MP
!LL  For Cray - Global  ; Enable defs GLOBAL
!LL
!LL A.Lorenc    <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  4.1  4/09/96:  Port to CRAY T3E  Deborah Salmond
!LL  4.3  31/1/97:  Bugfix for LAM T3E. Stuart Bell
!LL  5.2  21/10/99: Fix for problem of obs on a boundary. Adam Maycock
!LL  5.2  12/12/00:  change atbase,atright to elements of at_extremity
!LL                  -amend lasize for extra dimensions
!LL                  -remove refs to daco-type routines   B Macpherson
!    5.3  12/07/01:  amend final NP#PT pointers for S->N grid
!                    order in ND.  Correct treatment of
!                    pointers outside grid   B Macpherson
!    6.2  26/05/06:  Remove redundent defs. P.Selwood
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LL WORKS FOR ARAKAWA 'B' GRID ONLY
!LL
!LLEND------------------------------------------------------------------
      SUBROUTINE HINTCF (LWIND,LENOB,OBS_LAT,OBS_LONG,OBS_NO,           &
     &                   ROW_LENGTH,P_ROWS,                             &
     &                   CF1PT,CF2PT,CF3PT,CF4PT,                       &
     &                   NP1PT,NP2PT,NP3PT,NP4PT,                       &
     &                   ICODE,CMESSAGE)
!L    --------------------------------------------------
!
!
      IMPLICIT NONE
#include "acparm.h"
#include "commg.h"
#include "comacdg.h"
#include "comag.h"
#include "comacp.h"
#if defined(MPP)
#include "parvars.h"
#include "mppac.h"
#endif
      EXTERNAL TIMER
!     -----------------------------------------------------------
      LOGICAL LWIND              !IN switch to identify wind grid
      INTEGER LENOB              !IN number of obs
      INTEGER ROW_LENGTH,P_ROWS  !IN size of model grid
!     -----------------------------------------------------------
      REAL     OBS_LAT(LENOB)    !IN ob co-latitudes between 0 & PI
      REAL     OBS_LONG(LENOB)   !IN ob longitudes between 0 & 2*PI
!     NP#PT are pointers to the 4 gridpoints surrounding each ob.
      INTEGER  NP1PT(LENOB)      !OUT pointer to nearest model point
      INTEGER  NP2PT(LENOB)      !OUT pointer to same row, +or-1 pt.
      INTEGER  NP3PT(LENOB)      !OUT pointer to +or-1 row, same pt.
      INTEGER  NP4PT(LENOB)      !OUT pointer to +or-1 row, +or-1 pt.
!     CF#PT are the weights given to each pt pointed to by NP#PT.
      REAL     CF1PT(LENOB)      !OUT interpolation coeffs
      REAL     CF2PT(LENOB)      !OUT interpolation coeffs
      REAL     CF3PT(LENOB)      !OUT interpolation coeffs
      REAL     CF4PT(LENOB)      !OUT interpolation coeffs
      INTEGER  OBS_NO(LENOB)     !IN pointers to obs
!     -----------------------------------------------------------
      INTEGER ICODE              !OUT error code and message
      CHARACTER*256 CMESSAGE
!     -----------------------------------------------------------
!**   DYNAMIC ALLOCATION WITH LENOB
#if !defined(GLOBAL)
      REAL           WKLON(LENOB)
#endif
#if defined(MPP)
#if defined(GLOBAL)
      REAL           WKLON(LENOB)
#endif
#endif
      REAL           WORK5(LENOB)
!     CF#PT and NP1PT are used for workspace during the calculation
      INTEGER N_ROW(LENOB) ! nearest row to ob     minus 1
      INTEGER N_PNT(LENOB) ! nearest point to ob   minus 1
      INTEGER I_ROW(LENOB) ! increment to row other side of ob   (+or-1)
      INTEGER I_PNT(LENOB) ! increment to point other side of ob (+or-1)
!     -----------------------------------------------------------
      REAL                                                              &
     &               ZDLAT,ZDLONG,ZLATN,ZLONGW,R_ZDLAT,R_ZDLONG
#if !defined(GLOBAL)
      REAL           ZLONGE
#endif
      INTEGER JOB
!     -----------------------------------------------------------
#include "c_pi.h"
      REAL PI2P
      PARAMETER (PI2P = 2.0*PI)
!     -----------------------------------------------------------
      REAL Tiny
      PARAMETER (Tiny = 1.6E-7)
!     -----------------------------------------------------------
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('HINTCF  ',3)
!
! PART 1
! ------
! SET UP LATITUDE OF NORTHERN ROW AND LONGITUDE OF WESTERN POINT
! FOR EACH MODEL GRID
!
! PARAMETERS FOR MODEL GRID
!
!     GET MODEL GRID INFORMATION FROM COMMG AND CONVERT TO RADIANS
!     CONVERT XLATN TO CO-LATITUDE FIRST
!
      ZDLAT  = DLAT         * PI_OVER_180
      R_ZDLAT = 1.0/ZDLAT
      ZDLONG = DLONG        * PI_OVER_180
      R_ZDLONG = 1.0/ZDLONG
#if !defined(MPP)
      ZLATN  = (90.0-XLATN) * PI_OVER_180
#else
      ZLATN  = (90.0-LAT_N) * PI_OVER_180
#endif
#if !defined(MPP)
      ZLONGW = XLONGW       * PI_OVER_180
#else
      ZLONGW = LONG_W_MODEL * PI_OVER_180
#endif
#if !defined(GLOBAL)
#if !defined(MPP)
      ZLONGE = (XLONGW +(ROW_LENGTH-1)*DLONG ) * PI_OVER_180
#else
      ZLONGE = LONG_E_MODEL * PI_OVER_180
#endif
#endif
!     OFFSET WINDS
      IF (LWIND) THEN
        ZLATN = ZLATN +0.5*ZDLAT
        ZLONGW = ZLONGW+0.5*ZDLONG
#if !defined(GLOBAL)
        ZLONGE = ZLONGE+0.5*ZDLONG
#endif
      END IF
#if !defined(GLOBAL)
!     MAKE SURE ZLONGW,ZLONGE ARE IN RANGE 0-->2*PI
      IF (ZLONGW <  0.0)  THEN
        ZLONGW=ZLONGW+PI2P
      ELSEIF (ZLONGW >  PI2P) THEN
        ZLONGW=ZLONGW-PI2P
      ENDIF
      IF (ZLONGE <  0.0)  THEN
        ZLONGE=ZLONGE+PI2P
      ELSEIF (ZLONGE >  PI2P) THEN
        ZLONGE=ZLONGE-PI2P
      ENDIF
#endif
!
!
! PART 2
! ------
! CALCULATE (NEAREST-1) MODEL GRID PT TO OBSERVATION
!
!L    INPUT  VECTOR IN OBS_LAT  - COLATITUDE (RADIANS) SET UP IN AC
!L                  IN OBS_LONG - LONGITUDE  (RADIANS) SET UP IN AC
!
!L    OUTPUT VECTORS
!
!L    OUTPUT VECTORS OF INTERPOLATION COEFFICIENTS ARE IN
!L    CF1PT CF2PT CF3PT CF4PT NP1PT NP2PT NP3PT NP4PT
!
! 2.1
! ---
!
      DO JOB=1,LENOB
!       FIND NEAREST MODEL ROW minus 1.  in N_ROW
        N_ROW(JOB) = NINT((OBS_LAT(JOB)-ZLATN)*R_ZDLAT)
!       Wind obs at poles might just, because of rounding, have their
!       nearest row outside the grid.  Avoid this:
#if !defined(MPP)
        IF(LWIND)THEN
#else
        IF(LWIND.and.at_extremity(PSouth))THEN
#endif
          N_ROW(JOB) = MAX(N_ROW(JOB),0)
          N_ROW(JOB) = MIN(N_ROW(JOB),P_ROWS-2)
        ENDIF
!       N.B. checks that other obs are within the grid are made earlier
!
! 2.2
! ---
! FIND NEAREST POINT minus 1 ON ROW.   in N_PNT
!
#if !defined(GLOBAL)
#endif
!
#if defined(GLOBAL)
#if !defined(MPP)
        N_PNT(JOB) = NINT((OBS_LONG(JOB)           -ZLONGW)*R_ZDLONG)
#else
        WKLON(JOB)=0.0
        IF (OBS_LONG(JOB) <  ZLONGW.AND.at_extremity(PEast))THEN
        WKLON(JOB)=PI2P
        ENDIF
        N_PNT(JOB) = NINT((OBS_LONG(JOB)+WKLON(JOB)-ZLONGW)*R_ZDLONG)
#endif
#else
!       SPECIAL TREATMENT WHEN LIMITED AREA SPANS MERIDIAN.
!       OBS EAST OF THE MERIDIAN HAVE LONGITUDES < WESTERN BOUNDARY.
!       CALCULATE THEIR LONGITUDE DISPLACEMENT WITH EXTRA 2*PI SHIFT.
!
!       INITIALISE WORK ARRAY OF LONGITUDE SHIFTS
        WKLON(JOB)=0.0
        IF (ZLONGW  >   ZLONGE) THEN
!         ASSUME AREA STRADDLES MERIDIAN
!         SET UP LONGITUDE SHIFT FOR OBS BETWEEN MERIDIAN
!                                      AND EASTERN BOUNDARY.
          IF (OBS_LONG(JOB) >= 0.0.AND.OBS_LONG(JOB) <= ZLONGE)         &
     &                                     WKLON(JOB)= PI2P
          ENDIF
        N_PNT(JOB) = NINT((OBS_LONG(JOB)+WKLON(JOB)-ZLONGW)*R_ZDLONG)
#endif
!
! PART 3
! ------
! CALCULATE INTERPOLATION COEFFS. In lat in CF1PT, in long in CF2PT.
!
        CF1PT(JOB) = OBS_LAT(JOB) -(ZLATN +N_ROW(JOB)*ZDLAT )
        CF2PT(JOB) = OBS_LONG(JOB)-(ZLONGW+N_PNT(JOB)*ZDLONG)           &
     &               + WKLON(JOB)
!
! GET INCREMENT VECTORS TO OBTAIN FOUR SURROUNDING GRID POINTS
! INCLUDE POSSIBILITY OF OBS BEING ON SOUTHERN OR EASTERN HALO
! BOUNDARIES. (Cannot be on northern or western because of the way
! obs are distributed in RDOBS)
!
!   N.B Value of Tiny is chosen so that observations within
!   approx 1 metre of grid edge are affected.

!
        IF (CF1PT(JOB) >= 0.0) THEN
          I_ROW(JOB)=1
!         Set I_ROW to 0 if CF1PT is very small.
          IF (CF1PT(JOB)  <   Tiny) THEN
            I_ROW(JOB) = 0
          END IF
        ELSE
          I_ROW(JOB)=-1
        ENDIF
!
        IF (CF2PT(JOB) >= 0.0) THEN
          I_PNT(JOB) = 1
!         Set I_PNT to 0 if CF2PT is very small.
          IF (CF2PT(JOB)  <   Tiny) THEN
            I_PNT(JOB) = 0
          END IF
        ELSE
          I_PNT(JOB) = -1
        ENDIF
!
        CF1PT(JOB) = ABS(CF1PT(JOB))
        CF2PT(JOB) = ABS(CF2PT(JOB))
!
!       Make sure that all pointers are within grid, and different
        NP1PT(JOB)=N_ROW(JOB)+I_ROW(JOB) ! use NP1PT as workspace
        IF(NP1PT(JOB) <  0)THEN          ! row-1 is outside grid
          I_ROW(JOB)=1                   ! point instead to row+1
          CF1PT(JOB)=0.                  ! but do not use row+1.
        ENDIF
#if !defined(MPP)
        IF(LWIND)NP1PT(JOB)=NP1PT(JOB)+1 ! allow for 1 fewer rows
#else
        IF(LWIND.and.at_extremity(PSouth))NP1PT(JOB)=NP1PT(JOB)+1
!                                        ! allow for 1 fewer rows
#endif
        IF(NP1PT(JOB) >= P_ROWS)THEN     ! row+1 is outside grid
          I_ROW(JOB)=-1                  ! point instead to row-1
          CF1PT(JOB)=0.                  ! but do not use row-1.
        ENDIF
        IF ((N_PNT(JOB) + I_PNT(JOB)) <  1) THEN
          I_PNT(JOB) = 1
          CF2PT(JOB) = 0.
        ENDIF
        IF ((N_PNT(JOB) + I_PNT(JOB)) >  ROW_LENGTH) THEN
          I_PNT(JOB) = -1
          CF2PT(JOB) = 0.
        ENDIF
!
#if defined(GLOBAL)
! WRAP AROUND CHECK
!
        IF (N_PNT(JOB) >= ROW_LENGTH)THEN
             N_PNT(JOB) = N_PNT(JOB)-ROW_LENGTH
        ELSE IF (N_PNT(JOB) <  0)THEN
             N_PNT(JOB) = N_PNT(JOB)+ROW_LENGTH
        ENDIF
        NP1PT(JOB)=N_PNT(JOB)+I_PNT(JOB) ! use NP1PT as workspace
        IF (NP1PT(JOB) >= ROW_LENGTH)THEN
             I_PNT(JOB) = 1-ROW_LENGTH
        ELSE IF (NP1PT(JOB) <  0)THEN
          I_PNT(JOB) = ROW_LENGTH-1
        ENDIF
!
#endif
!       Convert to 2-D coeffs for 4 surrounding gridpoints
        CF4PT(JOB)  = CF1PT(JOB)*R_ZDLAT   ! lat coeff for row +or-1
        WORK5 (JOB) = CF2PT(JOB)*R_ZDLONG  ! long coeff for pt +or-1
        CF2PT(JOB)  = 1.0-CF4PT(JOB)       ! lat  coeff for row
        CF3PT(JOB)  = 1.0-WORK5(JOB)       ! long coeff for pt
!
        CF1PT(JOB) = CF2PT(JOB)*CF3PT(JOB) ! row, pt
        CF2PT(JOB) = CF2PT(JOB)*WORK5(JOB) ! row, pt +or-1
        CF3PT(JOB) = CF4PT(JOB)*CF3PT(JOB) ! row +or-1, pt
        CF4PT(JOB) = CF4PT(JOB)*WORK5(JOB) ! row +or-1, pt +or-1
!
        N_ROW(JOB)=MAX(N_ROW(JOB),0)
        N_PNT(JOB)=MAX(N_PNT(JOB),0)
!  up to now have assumed north most row is first
!  but ND order is reversed, so adjust use of N_ROW
      NP1PT(JOB) =N_PNT(JOB) + 1 + (P_ROWS -1 - N_ROW(JOB) )*ROW_LENGTH
        NP2PT(JOB) =NP1PT(JOB) + I_PNT(JOB)
!  respecting S->N order of ND fields,
!  use -I_ROW to point to other row adjacent to ob
      NP3PT(JOB) =NP1PT(JOB) - I_ROW(JOB)*ROW_LENGTH
        NP4PT(JOB) =NP3PT(JOB) + I_PNT(JOB)
#if defined(MPP)
      if(NP1PT(JOB) >  ROW_LENGTH*P_ROWS)then
          CMESSAGE = 'NP1PT(JOB) too big!!'
        write(0,*)'NP1PT(JOB) too big!!',NP1PT(JOB)
          ICODE=1
        endif
        if(NP1PT(JOB) <  1)then
          CMESSAGE = 'NP1PT(JOB) too small!!'
        write(0,*)'NP1PT(JOB) too small!!',NP1PT(JOB)
          ICODE=1
        endif

      if(NP2PT(JOB) >  ROW_LENGTH*P_ROWS)then
          CMESSAGE = 'NP2PT(JOB) too big!!'
          write(0,*)'NP2PT(JOB) too big!!',NP2PT(JOB)
          ICODE=1
        endif
        if(NP2PT(JOB) <  1)then
          CMESSAGE = 'NP2PT(JOB) too small!!'
          write(0,*)'NP2PT(JOB) too small!!',NP2PT(JOB)
          ICODE=1
        endif

      if(NP3PT(JOB) >  ROW_LENGTH*P_ROWS)then
          CMESSAGE = 'NP3PT(JOB) too big!!'
          write(0,*)'NP3PT(JOB) too big!!',NP3PT(JOB)
          ICODE=1
        endif
        if(NP3PT(JOB) <  1)then
          CMESSAGE = 'NP3PT(JOB) too small!!'
          write(0,*)'NP3PT(JOB) too small!!',NP3PT(JOB)
          ICODE=1
        endif

      if(NP4PT(JOB) >  ROW_LENGTH*P_ROWS)then
          CMESSAGE = 'NP4PT(JOB) too big!!'
          write(0,*)'NP4PT(JOB) too big!!',NP4PT(JOB)
          ICODE=1
        endif
        if(NP4PT(JOB) <  1)then
          CMESSAGE = 'NP4PT(JOB) too small!!'
          write(0,*)'NP4PT(JOB) too small!!',NP4PT(JOB)
          ICODE=1
        endif
#endif
      ENDDO ! JOB
!
!DIAG
!
!
!
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('HINTCF  ',4)
      RETURN
      END SUBROUTINE HINTCF
#endif
