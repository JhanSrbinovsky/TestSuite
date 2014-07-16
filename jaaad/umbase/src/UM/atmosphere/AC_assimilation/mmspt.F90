#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE MMSPT   MMSPTW   ---------------------------------------
!LL
!LL   2 Subroutines in deck : MMSPT and MMSPTW
!LL   MMSPTW is same as MMSPT for Limited Area Winds.
!LL
!LL  Purpose : Provide weighted Mean and S.D. plus extremes
!LL            on Model grids. Used for Weights and Increments
!LL
!LL  For Global  ; Enable defs GLOBAL
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL  3.1   15/1/93    Amended to be similar to PP routine MSXN  S Bell
!LL  4.1   31/05/96     The number of v points to be processed on a
!LL                     C grid differs from u by row_length. u,v
!LL                     dimensioned separately in call to WLLTOEQ.
!LL                     Requirement for VAR.
!LL                     Author I.Edmond       Reviewer D. Goddard
!    4.2 25/11/96: T3E mods Stuart Bell
!LL
!LL  4.1  4/09/96:  Port to CRAY T3E  Deborah Salmond
!LL  4.3      17/04/97     Tidy DEFS and code so that blank source is no
!LL                        produced (A. Brady)
!LL  4.3  12/2/97:  MMSPTW bugfix to rms-v,change format,
!LL                 + Add safety barriers Stuart Bell
!LL  4.4  23/7/97:  Comdeck revision to suit IAU diags. S Bell
!LL  5.2  13/12/00:  remove mmsptw, amend glsize  B Macpherson
!    5.3  08/10/01:  fix diagnostics S Cusack
!    5.3  05/12/01: Remove reference to the shmcomm and iovarsac include
!                   files, use a local dynamic array rather than common
!                   block for array 'PVALS_G'.                S. Cusack
!    6.0  10/10/03: Replace SHMEM with GCOM for SX6. Clive Jones
!    6.2  15/02/05: Optimisation for the SX-6 and removal of
!                   MPP def. P.Selwood
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
      SUBROUTINE MMSPT (PVALS,KLEV,KGRID,PNTLAB,                        &
     &                  ROW_LENGTH,P_ROWS)
!
!L    CALCULATE AREA-WEIGHTED MEAN & MEAN-SQUARE
!L    OF A FIELD ON THE MODEL GRID
!L
!L    KLEV : LEVEL OF FIELD (USED IN PRINT OUT ONLY)
!L
!L    IF KGRID=0 DATA ON P*/THETA MODEL GRID
!L            =1 DATA ON WIND MODEL GRID
!L
!L    16 CHARACTER PNTLAB IS PRINTED TO IDENTIFY MEAN & MEAN SQ.
!L
      IMPLICIT NONE
!
      EXTERNAL TIMER
!
#include "acparm.h"
! only for LTIMER_AC
#include "comacp.h"
#include "commg.h"
#include "c_pi.h"
!-----------------------------------------------------------------------
      INTEGER KLEV,KGRID                !IN LEVEL AND GRID IDENTIFER
      INTEGER ROW_LENGTH,P_ROWS         !IN MODEL DIMENSIONS
      REAL    PVALS(ROW_LENGTH,*)       !IN MODEL FIELD
      CHARACTER*16 PNTLAB               !IN CHARACTER IDENTIFER OF FIELD

#include "parvars.h"
#include "mppac.h"

!     LOCAL VARIABLES
      REAL PVALS_G(Max2DFieldSize)
      INTEGER i, j
      INTEGER row_length_global
      INTEGER p_rows_global
      REAL WT,SUMWT,SUM,ZROW1,ZLAT,ZDLAT,ZM,ZMS,ZMAX,ZMIN,ZRMS
      INTEGER JPTF,JPTL,JROWF,JROWL,JROW,JPT,NPTS
      INTEGER IMAXPT,IMAXRO,IMINPT,IMINRO

      INTEGER ICODE
      CHARACTER (LEN=80) :: CMESSAGE
!-----------------------------------------------------------------------
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('MMSPT   ',3)
      row_length_global=glsize(1,fld_type_p)
      p_rows_global=glsize(2,fld_type_p)


!     JPTF  = FIRST POINT IN ROW
!     JPTL  = LAST  POINT IN ROW
!     JROWF = FIRST ROW
!     JROWL = LAST  ROW

#if defined(GLOBAL)
!       FULL GLOBAL GRID USED
        JPTF  = 1
        JPTL  = row_length_global
        JROWF = 1
        JROWL = p_rows_global
#else
!       OUTSIDE TWO BOUNDARY POINTS OF ELF GRID NOT USED
        JPTF  = 3
        JPTL  = row_length_global-2
        JROWF = 3
        JROWL = p_rows_global-2
#endif

! Gather P_VALS onto a global field PVALS_G
! DEPENDS ON: gather_field
      Call Gather_Field( PVALS, PVALS_G, row_length, p_rows,            &
     &                   row_length_global, p_rows_global,              &
     &                   fld_type_p, halo_type_no_halo,                 &
     &                   0, gc_all_proc_group,                          &
     &                   icode, cmessage)


      if(mype == 0)then

      IF (KGRID == 0) THEN
!       P* GRID
        ZROW1 = XLATN - DLAT*(JROWL-JROWF+1)
      ELSEIF (KGRID == 1) THEN
!       WIND GRID - ONE LESS ROW
        ZROW1 = XLATN - DLAT*(JROWL-JROWF+0.5)
        JROWL = JROWL-1
      ENDIF

      ZLAT   = ZROW1*PI_OVER_180
      ZDLAT  = DLAT *PI_OVER_180

!L SET ACCUMULATORS
      ZMAX=pvals_g(JPTF+(JROWF-1)*row_length_global)
      ZMIN=pvals_g(JPTF+(JROWF-1)*row_length_global)

      ZM    = 0.0
      ZMS   = 0.0
      ZRMS  = 0.0
      SUMWT = 0.0
      WT    = 0.0
      IMAXPT= JPTF
      IMAXRO= JROWF
      IMINPT= JPTF
      IMINRO= JROWF

      DO 4 JROW = JROWF,JROWL
      WT = COS(ZLAT)
      SUMWT = SUMWT+WT

!     CALCULATE MEAN

      SUM = 0.0
      DO JPT=JPTF,JPTL
        SUM = SUM + pvals_g(JPT+(JROW-1)*row_length_global)
      ENDDO
      ZM = ZM + SUM*WT

!     CALCULATE MEAN SQUARE

      SUM = 0.0
      DO JPT=JPTF,JPTL
        SUM = SUM + pvals_g(JPT+(JROW-1)*row_length_global)*            &
     &              pvals_g(JPT+(JROW-1)*row_length_global)
      ENDDO
      ZMS = ZMS + SUM*WT

!     CALCULATE MAX

      DO JPT=JPTF,JPTL
        IF(pvals_g(JPT+(JROW-1)*row_length_global) >  ZMAX)THEN
          ZMAX=pvals_g(JPT+(JROW-1)*row_length_global)
          IMAXPT=JPT
          IMAXRO=JROW
        ENDIF
      ENDDO

!     CALCULATE MIN

      DO JPT=JPTF,JPTL
        IF(pvals_g(JPT+(JROW-1)*row_length_global) <  ZMIN)THEN
          ZMIN=pvals_g(JPT+(JROW-1)*row_length_global)
        IMINPT=JPT
        IMINRO=JROW
        ENDIF
      ENDDO

      ZLAT = ZLAT + ZDLAT
 4    CONTINUE

!     EVALUATE STATS AND WRITE OUT

      NPTS = JPTL-JPTF+1
      IF(SUMWT /= 0.) WT   = 1.0/(SUMWT*NPTS)
      ZM   = ZM  * WT
      ZMS  = ZMS * WT
      IF(ZMS >  0.)   ZRMS = SQRT(ZMS)

      PRINT 62,                                                         &
     &  PNTLAB,KLEV,ZM,ZRMS,ZMAX,IMAXRO,IMAXPT,ZMIN,IMINRO,IMINPT

62    FORMAT(1X,A16,2X,I4,' MEAN=',G12.5,' RMS=',G12.5,                 &
     & ' MAX=',G12.5,' AT (',I4,',',I4,')',                             &
     & ' MIN=',G12.5,' AT (',I4,',',I4,')')

      endif      ! mype == 0
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('MMSPT   ',4)
      RETURN
      END SUBROUTINE MMSPT
#endif
