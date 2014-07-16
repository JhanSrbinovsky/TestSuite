#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Construct preliminary STASH list of user requests
!
! Subroutine Interface:

      SUBROUTINE PRELIM(NRECS,                                          &
#include "argppx.h"
     &                 NTIMES,NLEVELS,ErrorStatus,CMESSAGE)
      IMPLICIT NONE

!  Description:
!  Constructs a preliminary STASH list of user requests. Uses interim
!  pointer system, by means of the "extra entry" NELEMP+1 in the LIST_S
!  array. At this stage, the input levels encompass all possible levels.
!  Called by STPROC.
!
!  Method:
!
!  Current code owner:  UM System Team
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Mar. 95    Original code.  S.J.Swarbrick
!   4.1     Apr. 96    Numerous improvements associated with wave model,
!                       correction of output-times table processing,
!                       comprehensive soft-abort system, etc.
!                                             S.J.Swarbrick
!   4.4     Sep. 97    Allow offset for sampling frequency
!                      S.D. Mullerworth
!LL 4.4    21/11/96   Allow daily mean timeseries. R.A.Stratton
!   4.4     Oct. 97    Added checking of error returns from TOTIMP.
!                         Shaun de Witt
!   4.5    18/11/98    Allow new sampling frequencies for vegetation.
!                      Richard Betts
!   5.2  15/11/00      Correct st_end_time_code when using time
!                      processing                    P.Selwood
!   5.0    22/11/99    Correct STASH input code for fields in
!                      secondary space - first step to enabling
!                      STASH output of these fields. R. Rawlins
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!   5.1    15/05/00    S-N ordering consistency correction. R Rawlins
!   5.3    23/10/01    Minor enhancement to error reporting. R Rawlins
!   5.5    17/02/03    Upgrade Wave model from 4.1 to 5.5. D.Holmes-Bell
!   6.0    18/08/03    Check that requested diagnostic exists in
!                      STASHmaster.  E. Leung
!   6.0  15/12/03  Update IOPN inline with
!                  stashmaster 30 digit option codes.  M.Hughes
!   6.2  23/11/05  Removed all references to the wavemodel.
!                  T.Edwards
!   6.2   16/11/05     Fix offsets for radiation/convection diagnostics.
!                      Improve warning reporting. P.Selwood
!
! 6.2  25/11/05 Functionality for improved time stepping, radiative
!               forcing and radiance code added for versions 3C
!               and 3Z of radiation code             (J.-C. Thelen)
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!
!  Global variables:

#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "parparm.h"
#include "typsize.h"
#include "cstash.h"
#include "stextend.h"
#include "model.h"
#include "cntlatm.h"
#include "stparam.h"

! Subroutine arguments


!   Scalar arguments with intent(out):

      INTEGER NRECS
      INTEGER NTIMES
      INTEGER NLEVELS ! Total no. of sets of levs for diags (inpt+outp)
      CHARACTER*80 CMESSAGE

! ErrorStatus:
      INTEGER ErrorStatus

! Local scalars:
      LOGICAL      MODEL_LEV
      LOGICAL      LMASK
      LOGICAL      LMEAN
      LOGICAL      LOFFSET
      INTEGER      TOTIMP
      INTEGER      I
      INTEGER      IBOT1
      INTEGER      IDIAG
      INTEGER      IDOMLEV
      INTEGER      IDOM_L
      LOGICAL      LDUM
      INTEGER      IFIRST
      INTEGER      IFIRST1
      INTEGER      ILAST
      INTEGER      ILAST1
      INTEGER      IM
      INTEGER      IMD
      INTEGER      IPLOF
      INTEGER      MODL_L
      INTEGER      ISEC_L
      INTEGER      ITEM_L
      INTEGER      ITIM_L
      INTEGER      ITIM
      INTEGER      ITOP1
      INTEGER      IUSE_L
      INTEGER      IX1
      INTEGER      IX2
      INTEGER      IY1
      INTEGER      IY2
      INTEGER      JLEV
      INTEGER      LEV_OFFSET
      INTEGER      LBVC
      INTEGER      IMAX          ! to find max of times-table
      INTEGER      ITIMLST       ! column of times-table
      INTEGER      item_chk

      CHARACTER (LEN=256)          :: CMESSAGE2
      CHARACTER (LEN=*), PARAMETER :: RoutineName='PRELIM'

! Function and subroutine calls:
      LOGICAL  DISCT_LEV
      INTEGER  EXPPXI
      EXTERNAL EXPPXI,LEVSRT,TSTMSK,LLTORC,LEVCOD,PSLCOM,PSLIMS

!- End of Header ------------------------------------------------------

! 0.1  Store output-times tables in array ITIM_S

      IF(NTIMES == 0) THEN
      DO I=1,NPROFTP
        IF (IOPT_T(I) == 2.AND.MODL_T(I) >  0) THEN
! Profile has output times list
!  MODL_T(I) labels internal model for times list
          DO ITIM=1,ITIM_T(I)
! DEPENDS ON: totimp
             ITIM_S(ITIM,I)=TOTIMP(ISER_T(ITIM,I),UNT3_T(I),MODL_T(I))
             if (ITIM_S(ITIM,I)  ==  -999) then
                 ErrorStatus = 100
                 write (cmessage,'(a,a,i3)')                            &
     &           'PRELIM:TOTIMP:Error in time period conversion',       &
     &           ' output times table no.=',i
                 write(6,*) cmessage
                 GOTO 9999
              endif
          END DO
          ITIM_S(ITIM_T(I)+1,I)=-1
        ELSE
          ITIM_S(1,I)=-1
        END IF
      END DO
      NTIMES=NPROFTP
      END IF

! 0.2  Store output levels lists in array LEVLST_S

      LEV_OFFSET=NLEVELS ! Initialised to 0 before entering this routine

! Loop over domain profiles in STASH basis file
      DO I=1,NDPROF
        IF (LEVB_D(I) == -1) THEN
! There is a levels list for this dom prof
          IF (IOPL_D(I) == 1.OR.IOPL_D(I) == 2.OR.                      &
     &                          IOPL_D(I) == 6    ) THEN
! Levs list contains model levs - list type is integer
             LLISTTY(I+LEV_OFFSET)='I'
          ELSE
! Not model levs - list type real
             LLISTTY(I+LEV_OFFSET)='R'
          END IF
! LEVT_D(I) = no. of levs in list 'I'
          LEVLST_S(1,I+LEV_OFFSET)=LEVT_D(I)

! Levels list 'I' was read into (R)LEVLST_D(J,I), J=1,LEVT_D(I),
!  by RDBASIS.
!  Transfer this levels list to (R)LEVLST_S(J,I+LEV_OFFSET),
!  J=2,LEVT_D(I)+1.

          DO JLEV=1,LEVT_D(I)
            IF (IOPL_D(I) == 1.OR.IOPL_D(I) == 2.OR.                    &
     &                            IOPL_D(I) == 6    ) THEN
!         Model levels
               LEVLST_S(JLEV+1,I+LEV_OFFSET)= LEVLST_D(JLEV,I)
            ELSE IF (IOPL_D(I) /= 5) THEN
!         Real levels
              RLEVLST_S(JLEV+1,I+LEV_OFFSET)=RLEVLST_D(JLEV,I)
            END IF
          END DO

          IPLOF=I+LEV_OFFSET

!   Sort this levels list into correct order (if not already in order)
! DEPENDS ON: levsrt
          CALL LEVSRT( LLISTTY(  IPLOF), LEVLST_S(1,IPLOF),             &
     &                LEVLST_S(2,IPLOF),RLEVLST_S(2,IPLOF))
        ELSE
! No levels list, i.e., the output from this diag. is on a
!    contiguous range of model levels
          LEVLST_S(1,I+LEV_OFFSET)=0
        END IF
      END DO  !  Domain profiles

      NLEVELS=NDPROF+LEV_OFFSET  ! NDPROF = no. of sets of input levels

      IF(NLEVELS >  NLEVLSTSP) THEN
        WRITE(6,*)                                                      &
     &  'PRELIM: TOO MANY LEVELS LISTS, ARRAYS OVERWRITTEN'
        CMESSAGE=                                                       &
     & 'PRELIM: TOO MANY LEVELS LISTS, ARRAYS OVERWRITTEN'
        GO TO 9999
      END IF

! Section 1. MAIN LOOP - loop over diag requests in STASH basis file

      IF(NDIAG >  0) THEN

      DO IDIAG=1,NDIAG

      MODL_L=MODL_B(IDIAG)
      ISEC_L=ISEC_B(IDIAG)
      ITEM_L=ITEM_B(IDIAG)
      IDOM_L=IDOM_B(IDIAG)
      IUSE_L=IUSE_B(IDIAG)
      ITIM_L=ITIM_B(IDIAG)

      item_chk=0
! DEPENDS ON: exppxi
      item_chk=EXPPXI(MODL_L,ISEC_L,ITEM_L,ppx_item_number,             &
#include "argppx.h"
     &                ErrorStatus,CMESSAGE)
      IF(ITEM_CHK /= ITEM_L)THEN
        WRITE(CMESSAGE2,*)'Diagnostic discarded ',                      &
     &       ' model ',modl_l,' section ',isec_l,' item',item_l,        &
     &       ' No stashmaster record'
        ERRORSTATUS = -10
! DEPENDS ON: ereport
        CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)

        ! Make this item null
        ITEM_B(IDIAG)=0
        GOTO 999
      ENDIF
      IF(ITIM_L /= 0) THEN       ! If the diag is not a null request

! Section 1.0  Extract data required for STASH processing from PPXI

        IF(NRECS == NRECDP) THEN
          WRITE(CMESSAGE2,*)                                            &
     &   'TOO MANY STASH LIST ENTRIES, REQUEST DENIED',                 &
     &   ' (M,S,I)',MODL_L,ISEC_L,ITEM_L

          ERRORSTATUS = -20
! DEPENDS ON: ereport
          CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
          GOTO 999
        END IF

! DEPENDS ON: exppxi
        VMSK    = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_version_mask ,      &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        ISPACE  = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_space_code   ,      &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        ITIMA   = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_timavail_code,      &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IGP     = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_grid_type    ,      &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        ILEV    = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_lv_code      ,      &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IBOT    = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_lb_code      ,      &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        ITOP    = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_lt_code      ,      &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IFLAG   = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_lev_flag     ,      &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
        DO I=1,6


! DEPENDS ON: exppxi
        IOPN(I) = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_opt_code+I-1 ,      &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
        END DO
! DEPENDS ON: exppxi
        IPSEUDO = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_pt_code      ,      &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IPFIRST = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_pf_code      ,      &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IPLAST  = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_pl_code      ,      &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        PTR_PROG= EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_ptr_code     ,      &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        LBVC    = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_lbvc_code    ,      &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)

! Check availability of diagnostic
! DEPENDS ON: tstmsk
        CALL TSTMSK(MODL_L,ISEC_L,LMASK,LDUM,ErrorStatus,CMESSAGE)
        IF(.NOT.LMASK) THEN
          WRITE(CMESSAGE2,*)                                            &
     &   'DIAGNOSTIC NOT AVAILABLE TO THIS VERSION ',                   &
     &   'REQUEST DENIED ',                                             &
     &   '(M,S,I)',                                                     &
     &                MODL_L,ISEC_L,ITEM_L

          ERRORSTATUS = -30
! DEPENDS ON: ereport
          CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
          GOTO 999
        END IF

        NRECS=NRECS+1

        LIST_S(st_model_code  ,NRECS)= MODL_L
        LIST_S(st_sect_no_code,NRECS)= ISEC_L
        LIST_S(st_item_code   ,NRECS)= ITEM_L
! Prelim pointer for 'child' records
        LIST_S(NELEMP+1       ,NRECS)= NRECS
        LIST_S(st_lookup_ptr  ,NRECS)=-1

! Set input code for STASH requests:
!  =0 Use primary or secondary field:       D1(SI(item,section,model))
!  =1 Use field in diagnostic space: STASHwork(SI(item,section,model))
!  =-j Use diagnostic at D1(LIST_S(st_output_addr,j))
        IF( (ISPACE == 2).OR.(ISPACE == 4)                              &
     &  .OR.(ISPACE == 7).OR.(ISPACE == 8) .OR.(ISPACE == 9)) THEN
          LIST_S(st_input_code,NRECS)=0
        ELSE
          LIST_S(st_input_code,NRECS)=1
        END IF

        IF((ITIMA >= 5).AND.(ITIMA <= 12)) THEN
          LMEAN=.TRUE.
        ELSE
          LMEAN=.FALSE.
        END IF


! 1.1   Expand the domain profile ---------------------------

!   Averaging and Weighting
        IM=IMSK_D(IDOM_L)
        IF ((IGP ==  2).OR.(IGP ==  3)    .OR.                          &
     &      (IGP == 12).OR.(IGP == 13))   THEN
! Diags only available over land/sea
          IF((IMSK_D(IDOM_L)  ==  1)      .AND.                         &
     &       (IGP == 3.OR.IGP == 13))     THEN
! Diag requested over land+sea, only available over sea
            IM=3
          ELSE IF((IMSK_D(IDOM_L)  ==  1) .AND.                         &
     &            (IGP == 2.OR.IGP == 12))THEN
! Diag requested over land+sea, only available over land
            IM=2
          ELSE IF((IMSK_D(IDOM_L)  ==  2) .AND.                         &
     &            (IGP == 3.OR.IGP == 13))THEN
! Diag requested over land, only available over sea
            WRITE(6,*)'PRELIM: CHANGED TO SEA DIAG'
            WRITE(6,*) 'MODEL,SECTION,ITEM ',                           &
     &                  MODL_L,ISEC_L,ITEM_L
            IM=3
          ELSE IF((IMSK_D(IDOM_L)  ==  3) .AND.                         &
     &            (IGP == 2.OR.IGP == 12))THEN
! Diag requested over sea, only available over land
            WRITE(6,*)'PRELIM: CHANGED TO LAND DIAG'
            WRITE(6,*) 'MODEL,SECTION,ITEM ',                           &
     &                  MODL_L,ISEC_L,ITEM_L
            IM=2
          END IF
        END IF

        LIST_S(st_gridpoint_code,NRECS)=IM+10*IMN_D(IDOM_L)
        LIST_S(st_weight_code   ,NRECS)=      IWT_D(IDOM_L)

!   Horizontal area
!    - convert lat/long spec to row/column numbers if appropriate;
!    - convert lat/long spec to equatorial lat/long if appropriate.
        IF(IOPA_D(IDOM_L) == 1) THEN
! Full domain
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,90,-90,0,360,                                 &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 2 ) THEN
! N Hemis
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,90,0,0,360,                                   &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 3 ) THEN
! S Hemis
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,0,-90,0,360,                                  &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 4 ) THEN
! 90N-30N
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,90,30,0,360,                                  &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 5 ) THEN
! 30S-90S
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,-30,-90,0,360,                                &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 6 ) THEN
! 30N-00N
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,30,00,0,360,                                  &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 7 ) THEN
! 00S-30S
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,00,-30,0,360,                                 &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 8 ) THEN
! 30N-30S
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,30,-30,0,360,                                 &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 9 ) THEN
! Other lat/long spec
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,INTH_D(IDOM_L),ISTH_D(IDOM_L),                &
     &                    IWST_D(IDOM_L),IEST_D(IDOM_L),                &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 10) THEN
! Grid point spec
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,90,-90,0,360,IY1,IY2,IX1,IX2)
          LIST_S(st_north_code,NRECS)=MIN(INTH_D(IDOM_L),IY2)
          LIST_S(st_south_code,NRECS)=MIN(ISTH_D(IDOM_L),IY2)
          LIST_S(st_west_code ,NRECS)=MIN(IWST_D(IDOM_L),IX2)
          LIST_S(st_east_code ,NRECS)=MIN(IEST_D(IDOM_L),IX2)
        ELSE
          WRITE(CMESSAGE2,*) 'INVALID DOMAIN AREA OPTION=',             &
     &                      IOPA_D(IDOM_L),                             &
     &                      '(M,S,I)',                                  &
     &                      MODL_L,ISEC_L,ITEM_L

          ERRORSTATUS = -35
! DEPENDS ON: ereport
          CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
          NRECS=NRECS-1
          GOTO 999
        END IF

! Input level setting
! DEPENDS ON: disct_lev
        MODEL_LEV=DISCT_LEV(ILEV,ErrorStatus,CMESSAGE)
        IF (MODEL_LEV) THEN
! Model levels
! Set bottom level
! DEPENDS ON: levcod
          CALL LEVCOD(IBOT,IBOT1,ErrorStatus,CMESSAGE)
! Set top level
! DEPENDS ON: levcod
          CALL LEVCOD(ITOP,ITOP1,ErrorStatus,CMESSAGE)
! Contig. range of model levels
          IF(IFLAG == 0) THEN
            LIST_S(st_input_bottom,NRECS)=IBOT1
            LIST_S(st_input_top   ,NRECS)=ITOP1
! Non-contig. levels list
          ELSE IF(IFLAG == 1) THEN
            LIST_S(st_input_bottom,NRECS)=-1
            LIST_S(st_input_top   ,NRECS)= 1
          END IF
        ELSE
! Non-model levels
          IF(ILEV == 3) THEN
!  Pressure levels
            LIST_S(st_input_bottom,NRECS)=-1
            LIST_S(st_input_top   ,NRECS)= 2
          ELSE IF(ILEV == 4) THEN
!  Height levels
            LIST_S(st_input_bottom,NRECS)=-1
            LIST_S(st_input_top   ,NRECS)= 3
          ELSE IF(ILEV == 5) THEN
!  Special levels
            LIST_S(st_input_bottom,NRECS)=100
            LIST_S(st_input_top   ,NRECS)=LBVC
          ELSE IF(ILEV == 7) THEN
!  Theta levels
            LIST_S(st_input_bottom,NRECS)=-1
            LIST_S(st_input_top   ,NRECS)= 4
          ELSE IF(ILEV == 8) THEN
!  PV levels
            LIST_S(st_input_bottom,NRECS)=-1
            LIST_S(st_input_top   ,NRECS)= 5
          ELSE IF(ILEV == 9) THEN
!  Cloud threshold levels
            LIST_S(st_input_bottom,NRECS)=-1
            LIST_S(st_input_top   ,NRECS)= 6
          END IF
        END IF

! Output level specification
! DEPENDS ON: disct_lev
        MODEL_LEV=DISCT_LEV(ILEV,ErrorStatus,CMESSAGE)
        IF (MODEL_LEV) THEN
! Model levels
          IF (LEVB_D(IDOM_L) >= 0) THEN
! Contiguous range of model levels
            LIST_S(st_output_bottom,NRECS)=MAX(LEVB_D(IDOM_L),IBOT1)
            LIST_S(st_output_top   ,NRECS)=MIN(LEVT_D(IDOM_L),ITOP1)
            IF ((LEVB_D(IDOM_L) <  IBOT1).OR.                           &
     &          (LEVT_D(IDOM_L) >  ITOP1)) THEN
            WRITE(CMESSAGE2,*)                                          &
     &     'DIAGNOSTIC HAS LEVEL RANGE OUT OF BOUNDS; CORRECTED ',      &
     &     '(M,S,I) ',                                                  &
     &      MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-40
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            END IF
            IF ( (  TS_D(IDOM_L) ==   'Y').AND.                         &
     &          ((LEVB_D(IDOM_L) <  IBOT1).OR.                          &
     &           (LEVT_D(IDOM_L) >  ITOP1))    ) THEN
              WRITE(CMESSAGE2,*)                                        &
     &       'TIME SERIES DOMAIN',                                      &
     &       'HAS INCONSISTENT LEVELS; DIAGNOSTIC IGNORED',             &
     &       ' (M,S,I) ',                                               &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-50
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
            IF ((LEVT_D(IDOM_L) <  IBOT1).OR.                           &
     &          (LEVB_D(IDOM_L) >  ITOP1)) THEN
              WRITE(CMESSAGE2,*)                                        &
     &       'DIAGNOSTIC HAS TOP/BOT LEVELS INCONSISTENT; DIAG IGNORED',&
     &       ' (M,S,I) ',                                               &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-60
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
          ELSE
! Non-contig. list of model levels
            LIST_S(st_output_bottom,NRECS)=-(IDOM_L+LEV_OFFSET)
            LIST_S(st_output_top   ,NRECS)=1
          END IF
        ELSE
! Non-model levels
          IF(ILEV == 5) THEN
! Special level
            LIST_S(st_output_bottom,NRECS)=100
            LIST_S(st_output_top   ,NRECS)=LBVC
          ELSE IF(ILEV == 3) THEN
! Pressure levels
            LIST_S(st_output_bottom,NRECS)=-(IDOM_L+LEV_OFFSET)
            LIST_S(st_output_top   ,NRECS)=2
          ELSE IF(ILEV == 4) THEN
! Height levels
            LIST_S(st_output_bottom,NRECS)=-(IDOM_L+LEV_OFFSET)
            LIST_S(st_output_top   ,NRECS)=3
          ELSE IF(ILEV == 7 ) THEN
! Theta levels
            LIST_S(st_output_bottom,NRECS)=-(IDOM_L+LEV_OFFSET)
            LIST_S(st_output_top   ,NRECS)=4
          ELSE IF(ILEV == 8 ) THEN
! PV levels
            LIST_S(st_output_bottom,NRECS)=-(IDOM_L+LEV_OFFSET)
            LIST_S(st_output_top   ,NRECS)=5
          ELSE IF(ILEV == 9 ) THEN
! Cloud threshold levels
            LIST_S(st_output_bottom,NRECS)=-(IDOM_L+LEV_OFFSET)
            LIST_S(st_output_top   ,NRECS)=6
          ELSE
            WRITE(CMESSAGE2,*) 'DOMAIN LEVEL OPTION=',IOPL_D(IDOM_L),   &
     &                        'DIAG IGNORED. (M,S,I) ',                 &
     &                        MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-70
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            NRECS=NRECS-1
            GOTO 999
          END IF
        END IF

! Output pseudo-levels level setting
        IF(IPSEUDO /= PLT_D(IDOM_L)) THEN
            WRITE(CMESSAGE2,*)                                          &
     &     'DIAGNOSTIC HAS ',                                           &
     &     'INVALID PSEUDO LEVEL TYPE; IGNORED.',                       &
     &     ' (M,S,I)',                                                  &
     &      MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-80
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            NRECS=NRECS-1
            GOTO 999
        END IF
        LIST_S(st_pseudo_in,NRECS)=0  !(This is set in INPUTL)
        IF(IPSEUDO >  0) THEN
! Pseudo levels list for this diagnostic
            LIST_S(st_pseudo_out,NRECS)=PLPOS_D(IDOM_L)
            LENPLST(PLPOS_D(IDOM_L))   =PLLEN_D(IDOM_L)
            IFIRST=PSLIST_D(1,PLPOS_D(IDOM_L))
            ILAST =PSLIST_D(PLLEN_D(IDOM_L),PLPOS_D(IDOM_L))
! Check pseudo level limits
! DEPENDS ON: pslims
            CALL PSLIMS(IPFIRST,IPLAST,IFIRST1,ILAST1)
            IF(IFIRST <  IFIRST1) THEN
              WRITE(CMESSAGE2,*)                                        &
     &       'DIAGNOSTIC HAS ',                                         &
     &       'FIRST PSEUDO LEVEL TOO LOW; IGNORED.',                    &
     &       '(M,S,I) ',                                                &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-90
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
            IF(ILAST >  ILAST1) THEN
              WRITE(CMESSAGE2,*)                                        &
     &       'DIAGNOSTIC HAS ',                                         &
     &       'LAST PSEUDO LEVEL TOO HIGH; IGNORED',                     &
     &       ' (M,S,I)',                                                &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-95
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
        ELSE
            LIST_S(st_pseudo_out,NRECS)=0
        END IF

! Time-series domain profiles
        IF(TS_D(IDOM_L) == 'Y') THEN
! Pointer for location of time series
            LIST_S(st_series_ptr,NRECS)=NPOS_TS(IDOM_L)
        ELSE
          LIST_S(st_series_ptr,NRECS)=0
        END IF

! 1.2   Expand the useage profile --------------------------

        IF (LOCN_U(IUSE_L) == 5) THEN                    ! PP file

          IF(LMEAN) THEN
            LIST_S(st_output_code,NRECS)=-27
            LIST_S(st_macrotag,NRECS)=0
          ELSE
            WRITE(6,*)                                                  &
     &     'MESSAGE FROM ROUTINE PRELIM: DIAGNOSTIC REQUEST HAS ',      &
     &     'OUTPUT DESTINATION CODE 5 (CLIMATE MEAN PP FILE) ',         &
     &     'BUT DIAGNOSTIC IS NOT A CLIMATE MEAN; REQUEST IGNORED'

            WRITE(CMESSAGE2,*)                                          &
     &     'DIAGNOSTIC IS NOT CLIMATE MEAN, BUT SENT TO MEAN FILE:',    &
     &     'IGNORED. (M,S,I) ',                                         &
     &      MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-100
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            NRECS=NRECS-1
            GOTO 999
          END IF

        ELSE IF (LMEAN) THEN

            WRITE(6,*)                                                  &
     &     'MESSAGE FROM ROUTINE PRELIM: DIAGNOSTIC REQUEST IS A ',     &
     &     'CLIMATE MEAN - SHOULD HAVE OUTPUT DESTINATION CODE 5 ',     &
     &     '(CLIMATE MEAN PP FILE); REQUEST IGNORED'

            WRITE (CMESSAGE2,*)                                         &
     &      'DIAGNOSTIC IS CLIMATE MEAN, INCORRECT DESTINATION:',       &
     &      'IGNORED. (M,S,I) ',                                        &
     &       MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-110
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            NRECS=NRECS-1
            GOTO 999

        ELSE IF (LOCN_U(IUSE_L) == 3) THEN                ! PP file

          LIST_S(st_output_code,NRECS)=-IUNT_U(IUSE_L)
          IF (UNT1_T(ITIM_L)=="DU" .AND. UNT1_T(ITIM_L)=="DU") THEN
             ! Special tag for dump frequency means
             LIST_S(st_macrotag,NRECS)=-999
          ELSE
             LIST_S(st_macrotag,NRECS)=0
          ENDIF

        ELSE IF (LOCN_U(IUSE_L) == 1) THEN ! Dump store: set user tag

          LIST_S(st_output_code,NRECS)=1
          LIST_S(st_macrotag,NRECS)=IUNT_U(IUSE_L)

        ELSE IF (LOCN_U(IUSE_L) == 6) THEN ! Secondary dump store:
                                           !             set user tag
          LIST_S(st_output_code,NRECS)=2
          LIST_S(st_macrotag,NRECS)=IUNT_U(IUSE_L)

        ELSE IF (LOCN_U(IUSE_L) == 2) THEN ! Climate mean: tag set
                                           !   1000*(time mean tag)
          LIST_S(st_output_code,NRECS)=1
          LIST_S(st_macrotag,NRECS)=IUNT_U(IUSE_L)*1000

        ELSE IF (LOCN_U(IUSE_L) == 4)THEN  ! Printed output

          LIST_S(st_output_code,NRECS)=7
          LIST_S(st_macrotag,NRECS)=0

        ELSE

          WRITE(CMESSAGE2,*)                                            &
     &    'INVALID USEAGE OPTION=', LOCN_U(IUSE_L),                     &
     &    ': DIAGNOSTIC IGNORED. (M,S,I) ',                             &
     &     MODL_L,ISEC_L,ITEM_L

          ERRORSTATUS=-120
! DEPENDS ON: ereport
          CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
          NRECS=NRECS-1
          GOTO 999

        END IF

! 1.3   Expand the time profile ------------------------------

! Initialise as single time field

!   Set time processing record

        IF (LMEAN) THEN
          IF (ITYP_T(ITIM_L) /= 1) THEN
            WRITE(CMESSAGE2,*)                                          &
     &     'CLIMATE MEANS MUST NOT BE TIME PROCESSED.',                 &
     &     '(M,S,I) ',                                                  &
     &      MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-130
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
          END IF
          LIST_S(st_proc_no_code,NRECS)=1
        ELSE
          LIST_S(st_proc_no_code,NRECS)=ITYP_T(ITIM_L)
        END IF

! Initialise offset to 0
            LIST_S(st_offset_code,NRECS)=0
!   Set period record

        IF (ITYP_T(ITIM_L) == 1.OR.LMEAN) THEN        ! No period
          LIST_S(st_period_code,NRECS)=0
        ELSE IF((INTV_T(ITIM_L) == -1).AND.                             &
     &          (ITYP_T(ITIM_L) == 2)) THEN
          LIST_S(st_period_code,NRECS)=-1
        ELSE
          LIST_S(st_period_code,NRECS)=                                 &
! DEPENDS ON: totimp
     &           TOTIMP(INTV_T(ITIM_L),UNT1_T(ITIM_L),MODL_L)
          if (LIST_S(st_period_code,NRECS)  ==  -999) then
              ErrorStatus = 101
              write (cmessage,'(a,a,i2,a,i3,a,i3)')                     &
     &        'PRELIM:TOTIMP:Error in time period conversion',          &
     &        ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
              write(6,*) cmessage
              GOTO 9999
           endif
        END IF

        IF (LMEAN.AND.(IOPT_T(ITIM_L) /= 1)) THEN
          WRITE(CMESSAGE2,*)                                            &
     &   'CLIMATE MEANS MUST USE STANDARD FREQUENCY. DIAG IGNORED.',    &
     &   '(M,S,I) ',                                                    &
     &     MODL_L,ISEC_L,ITEM_L

          ERRORSTATUS=-140
! DEPENDS ON: ereport
          CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
          NRECS=NRECS-1
          GOTO 999
        END IF

        IF(IOPT_T(ITIM_L) == 1) THEN
!Regular output times
          LIST_S(st_freq_code,NRECS)=                                   &
! DEPENDS ON: totimp
     &           TOTIMP(IFRE_T(ITIM_L),UNT3_T(ITIM_L),MODL_L)
          if (LIST_S(st_freq_code,NRECS)  ==  -999) then
              ErrorStatus = 102
              write (cmessage,'(a,a,i2,a,i3,a,i3)')                     &
     &        'PRELIM:TOTIMP:Error in time period conversion',          &
     &        ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
              write(6,*) cmessage
              GOTO 9999
           endif
          LIST_S(st_start_time_code,NRECS)=                             &
! DEPENDS ON: totimp
     &           TOTIMP(ISTR_T(ITIM_L),UNT3_T(ITIM_L),MODL_L)
          if (LIST_S(st_start_time_code,NRECS)  ==  -999) then
              ErrorStatus = 103
              write (cmessage,'(a,a,i2,a,i3,a,i3)')                     &
     &        'PRELIM:TOTIMP:Error in time period conversion',          &
     &        ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
              write(6,*) cmessage
              GOTO 9999
           endif
          IF(IEND_T(ITIM_L) == -1) THEN
             LIST_S(st_end_time_code,NRECS)=-1
          ELSE
             LIST_S(st_end_time_code,NRECS)=                            &
! DEPENDS ON: totimp
     &           TOTIMP(IEND_T(ITIM_L),UNT3_T(ITIM_L),MODL_L)
          ENDIF
          if (LIST_S(st_end_time_code,NRECS)  ==  -999) then
              ErrorStatus = 104
              write (cmessage,'(a,a,i2,a,i3,a,i3)')                     &
     &        'PRELIM:TOTIMP:Error in time period conversion',          &
     &        ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
              write(6,*) cmessage
              GOTO 9999
           endif

!   Set end time to -1 if output requested to end of run

!   Correct start time for radiation, periodic convection, leaf
!   phenology and vegetation competition
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
          IF((ITIMA == 2).AND.(A_LW_RADSTEP_DIAG /= 1)) THEN
            IMD=MOD(LIST_S(st_start_time_code,NRECS),A_LW_RADSTEP_DIAG)
            LIST_S(st_start_time_code,NRECS)=                           &
     &      LIST_S(st_start_time_code,NRECS)+1-IMD
            LOFFSET=.TRUE.
          ELSE IF((ITIMA == 3).AND.(A_SW_RADSTEP_DIAG /= 1)) THEN
            IMD=MOD(LIST_S(st_start_time_code,NRECS),A_SW_RADSTEP_DIAG)
            LIST_S(st_start_time_code,NRECS)=                           &
     &      LIST_S(st_start_time_code,NRECS)+1-IMD
            LOFFSET=.TRUE.
#else
          IF((ITIMA == 2).AND.(A_LW_RADSTEP /= 1)) THEN
            IMD=MOD(LIST_S(st_start_time_code,NRECS),A_LW_RADSTEP)
            LIST_S(st_start_time_code,NRECS)=                           &
     &      LIST_S(st_start_time_code,NRECS)+1-IMD
            LOFFSET=.TRUE.
          ELSE IF((ITIMA == 3).AND.(A_SW_RADSTEP /= 1)) THEN
            IMD=MOD(LIST_S(st_start_time_code,NRECS),A_SW_RADSTEP)
            LIST_S(st_start_time_code,NRECS)=                           &
     &      LIST_S(st_start_time_code,NRECS)+1-IMD
            LOFFSET=.TRUE.
#endif
          ELSE IF((ITIMA == 13).AND.(A_CONV_STEP /= 1)) THEN
            IMD=MOD(LIST_S(st_start_time_code,NRECS),A_CONV_STEP)
            LIST_S(st_start_time_code,NRECS)=                           &
     &      LIST_S(st_start_time_code,NRECS)+1-IMD
            LOFFSET=.TRUE.
          ELSE IF((ITIMA == 14).AND.(PHENOL_PERIOD /= 1)) THEN
            IMD=MOD(LIST_S(st_start_time_code,NRECS),PHENOL_PERIOD)
            LIST_S(st_start_time_code,NRECS)=                           &
     &      LIST_S(st_start_time_code,NRECS)+1-IMD
            LOFFSET=.TRUE.
          ELSE IF((ITIMA == 15).AND.(TRIFFID_PERIOD /= 1)) THEN
            IMD=MOD(LIST_S(st_start_time_code,NRECS),TRIFFID_PERIOD)
            LIST_S(st_start_time_code,NRECS)=                           &
     &      LIST_S(st_start_time_code,NRECS)+1-IMD
            LOFFSET=.TRUE.
          ELSE
            LOFFSET=.FALSE.
          END IF
        ELSE IF(IOPT_T(ITIM_L) == 2) THEN
!List of specified output times
            LIST_S(st_freq_code,NRECS)=-ITIM_L
        ELSE
          WRITE(CMESSAGE2,*)                                            &
     &    'INVALID OUTPUT TIMES CODE. DIAG IGNORED.',                   &
     &    '(M,S,I) ',                                                   &
     &     MODL_L,ISEC_L,ITEM_L

          ERRORSTATUS=-150
! DEPENDS ON: ereport
          CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
          NRECS=NRECS-1
          GOTO 999
        END IF

        IF (LMEAN) LIST_S(st_freq_code,NRECS)=1

        IF ((LIST_S(st_proc_no_code,NRECS) >  1).AND.                   &
     &      (LIST_S(st_proc_no_code,NRECS) <= 6)) THEN
! Other than single time field
          IF(NRECS >= NRECDP) THEN
            WRITE(CMESSAGE2,*)                                          &
     &     'TOO MANY S_LIST REQUESTS. REQUEST IGNORED',                 &
     &     '(M,S,I) ',                                                  &
     &      MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-160
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            NRECS=NRECS-1
            GOTO 999
          END IF

          DO I=1,NELEMP+1          ! Copy stash list forward
            LIST_S(I,NRECS+1)=LIST_S(I,NRECS)
          END DO

          IF(LOFFSET) THEN         ! Rad or conv timesteps,
                                   !       1 alresdy added
            LIST_S(st_start_time_code,NRECS+1)=                         &
     &      LIST_S(st_start_time_code,NRECS+1)-1
            IF (LIST_S(st_period_code,NRECS) /= -1) THEN
! Offsets are added to start time
             LIST_S(st_offset_code,NRECS)=                              &
! DEPENDS ON: totimp
     &           TOTIMP(IOFF_T(ITIM_L),UNT2_T(ITIM_L),MODL_L)
             if (LIST_S(st_offset_code,NRECS)  ==  -999) then
               ErrorStatus = 1
               write (cmessage,'(a,a,i2,a,i3,a,i3)')                    &
     &         'PRELIM:TOTIMP:Error in time period conversion',         &
     &         ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
               write(6,*) cmessage
               GOTO 9999
              endif
              LIST_S(st_start_time_code,NRECS)=                         &
     &        LIST_S(st_start_time_code,NRECS)-                         &
     &        LIST_S(st_period_code,NRECS) +                            &
     &        LIST_S(st_offset_code,NRECS)
            ELSE
              LIST_S(st_start_time_code,NRECS)=1
            END IF

          ELSE

            IF (LIST_S(st_period_code,NRECS) /= -1) THEN
! Offsets are added to start time
             LIST_S(st_offset_code,NRECS)=                              &
! DEPENDS ON: totimp
     &           TOTIMP(IOFF_T(ITIM_L),UNT2_T(ITIM_L),MODL_L)
             if (LIST_S(st_offset_code,NRECS)  ==  -999) then
               ErrorStatus = 1
               write (cmessage,'(a,a,i2,a,i3,a,i3)')                    &
     &         'PRELIM:TOTIMP:Error in time period conversion',         &
     &         ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
               write(6,*) cmessage
               GOTO 9999
              endif
              LIST_S(st_start_time_code,NRECS)=                         &
     &        LIST_S(st_start_time_code,NRECS)-                         &
     &        LIST_S(st_period_code,NRECS)+1+                           &
     &        LIST_S(st_offset_code,NRECS)
            ELSE
              LIST_S(st_start_time_code,NRECS)=1
            END IF

          END IF

          IF(LIST_S(st_start_time_code,NRECS) <  1) THEN
            WRITE(CMESSAGE2,*)                                          &
     &     'DIAGNOSTIC START TIME BEFORE PERIOD, SETTING TO 1.',        &
     &     '(M,S,I) ',                                                  &
     &      MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-170
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            LIST_S(st_start_time_code,NRECS)=1
          END IF

! Check if offset corresponds to a valid timestep.
! If not, reject the diagnostic.
          SELECT CASE (ITIMA)
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
            CASE (2)   ! Long-Wave Radiation
              IF (MOD(LIST_S(st_offset_code,NRECS),A_LW_RADSTEP_DIAG)   &
     &            /= 0)                                                 &
     &        THEN
                WRITE(CMESSAGE2,*)                                      &
     &              'OFFSET DOES NOT AGREE WITH LW_RAD STEP.',          &
     &              ' DIAG IGNORED. (M,S,I)',                           &
     &               MODL_L,ISEC_L,ITEM_L

                ERRORSTATUS=-180
! DEPENDS ON: ereport
                CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
                NRECS=NRECS-1
                GOTO 999
              END IF

            CASE (3)   ! Short-Wave Radiation
              IF (MOD(LIST_S(st_offset_code,NRECS),A_SW_RADSTEP_DIAG)   &
     &            /= 0)                                                 &
     &        THEN
                WRITE(CMESSAGE2,*)                                      &
     &              'OFFSET DOES NOT AGREE WITH SW_RAD STEP.',          &
     &              ' DIAG IGNORED. (M,S,I)',                           &
     &               MODL_L,ISEC_L,ITEM_L

                ERRORSTATUS=-190
! DEPENDS ON: ereport
                CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
                NRECS=NRECS-1
                GOTO 999
              END IF
#else
            CASE (2)   ! Long-Wave Radiation
              IF (MOD(LIST_S(st_offset_code,NRECS),A_LW_RADSTEP)        &
     &            /= 0)                                                 &
     &        THEN
                WRITE(CMESSAGE2,*)                                      &
     &              'OFFSET DOES NOT AGREE WITH LW_RAD STEP.',          &
     &              ' DIAG IGNORED. (M,S,I)',                           &
     &               MODL_L,ISEC_L,ITEM_L

                ERRORSTATUS=-180
! DEPENDS ON: ereport
                CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
                NRECS=NRECS-1
                GOTO 999
              END IF

            CASE (3)   ! Short-Wave Radiation
              IF (MOD(LIST_S(st_offset_code,NRECS),A_SW_RADSTEP)        &
     &            /= 0)                                                 &
     &        THEN
                WRITE(CMESSAGE2,*)                                      &
     &              'OFFSET DOES NOT AGREE WITH SW_RAD STEP.',          &
     &              ' DIAG IGNORED. (M,S,I)',                           &
     &               MODL_L,ISEC_L,ITEM_L

                ERRORSTATUS=-190
! DEPENDS ON: ereport
                CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
                NRECS=NRECS-1
                GOTO 999
              END IF
#endif
            CASE (13)   ! Convection
              IF (MOD(LIST_S(st_offset_code,NRECS),A_CONV_STEP)/= 0)    &
     &        THEN
                WRITE(CMESSAGE2,*)                                      &
     &              'OFFSET DOES NOT AGREE WITH CONVECT STEP.',         &
     &              ' DIAG IGNORED. (M,S,I)',                           &
     &               MODL_L,ISEC_L,ITEM_L

                ERRORSTATUS=-200
! DEPENDS ON: ereport
                CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
                NRECS=NRECS-1
                GOTO 999
              END IF

            CASE (14)   ! Leaf Phenology
              IF (MOD(LIST_S(st_offset_code,NRECS),PHENOL_PERIOD)/= 0)  &
     &        THEN
                WRITE(CMESSAGE2,*)                                      &
     &              'OFFSET DOES NOT AGREE WITH PHENOL PERIOD.',        &
     &              ' DIAG IGNORED. (M,S,I)',                           &
     &               MODL_L,ISEC_L,ITEM_L

                ERRORSTATUS=-210
! DEPENDS ON: ereport
                CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
                NRECS=NRECS-1
                GOTO 999
              END IF

            CASE (15)   ! Triffid
              IF (MOD(LIST_S(st_offset_code,NRECS),TRIFFID_PERIOD)/= 0) &
     &        THEN
                WRITE(CMESSAGE2,*)                                      &
     &              'OFFSET DOES NOT AGREE WITH TRIFFID PERIOD.',       &
     &              ' DIAG IGNORED. (M,S,I)',                           &
     &               MODL_L,ISEC_L,ITEM_L

                ERRORSTATUS=-220
! DEPENDS ON: ereport
                CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
                NRECS=NRECS-1
                GOTO 999
              END IF
          END SELECT

          LIST_S(st_proc_no_code ,NRECS+1)=1

          LIST_S(st_input_bottom ,NRECS+1)=                             &
     &    LIST_S(st_output_bottom,NRECS  )

          LIST_S(st_input_top    ,NRECS+1)=                             &
     &    LIST_S(st_output_top   ,NRECS  )

          LIST_S(st_input_code   ,NRECS+1)=-NRECS
          LIST_S(st_output_code  ,NRECS  )=1
          LIST_S(st_series_ptr   ,NRECS+1)=0
          LIST_S(NELEMP+1        ,NRECS+1)=NRECS+1

          LIST_S(st_freq_code,NRECS)=                                   &
                                                    ! Frequency
! DEPENDS ON: totimp
     &    TOTIMP(ISAM_T(ITIM_L),UNT2_T(ITIM_L),MODL_L)
          if (LIST_S(st_freq_code,NRECS)  ==  -999) then
             ErrorStatus = 105
               write (cmessage,'(a,a,i2,a,i3,a,i3)')                    &
     &         'PRELIM:TOTIMP:Error in time period conversion',         &
     &         ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
               write(6,*) cmessage
             GOTO 9999
          endif

!   Correct frequency for radiation, periodic convection, leaf
!   phenology and vegetation competition

          IF (ITIMA == 2) THEN
            IF (LIST_S(st_freq_code,NRECS) == 1) THEN
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
              LIST_S(st_freq_code,NRECS)=A_LW_RADSTEP_DIAG
            ELSE IF                                                     &
     &      (MOD(LIST_S(st_freq_code,NRECS),A_LW_RADSTEP_DIAG) /= 0)    &
     &       THEN
#else
              LIST_S(st_freq_code,NRECS)=A_LW_RADSTEP
            ELSE IF                                                     &
     &      (MOD(LIST_S(st_freq_code,NRECS),A_LW_RADSTEP) /= 0) THEN
#endif
              WRITE(CMESSAGE2,*)                                        &
     &       'INCORRECT SAMPLING FOR LW_RADSTEP. FREQ=',                &
     &        LIST_S(st_freq_code,NRECS), ':IGNORED.',                  &
     &       '(M,S,I) ',                                                &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-225
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
          ELSE IF(ITIMA == 3) THEN
            IF (LIST_S(st_freq_code,NRECS) == 1) THEN
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
              LIST_S(st_freq_code,NRECS)=A_SW_RADSTEP_DIAG
            ELSE IF                                                     &
     &      (MOD(LIST_S(st_freq_code,NRECS),A_SW_RADSTEP_DIAG) /= 0)    &
     &       THEN
#else
              LIST_S(st_freq_code,NRECS)=A_SW_RADSTEP
            ELSE IF                                                     &
     &      (MOD(LIST_S(st_freq_code,NRECS),A_SW_RADSTEP) /= 0) THEN
#endif
              WRITE(CMESSAGE2,*)                                        &
     &       'INCORRECT SAMPLING FOR SW_RADSTEP. FREQ=',                &
     &        LIST_S(st_freq_code,NRECS), ':IGNORED.',                  &
     &       '(M,S,I) ',                                                &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-230
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
          ELSE IF(ITIMA == 13) THEN
            IF (LIST_S(st_freq_code,NRECS) == 1) THEN
              LIST_S(st_freq_code,NRECS)=A_CONV_STEP
            ELSE IF                                                     &
     &      (MOD(LIST_S(st_freq_code,NRECS),A_CONV_STEP) /= 0) THEN
              WRITE(CMESSAGE2,*)                                        &
     &       'INCORRECT SAMPLING FOR CONV_STEP. FREQ=',                 &
     &        LIST_S(st_freq_code,NRECS), ':IGNORED.',                  &
     &       '(M,S,I) ',                                                &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-240
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
          ELSE IF(ITIMA == 14) THEN
            IF (LIST_S(st_freq_code,NRECS) == 1) THEN
              LIST_S(st_freq_code,NRECS)=PHENOL_PERIOD
            ELSE IF                                                     &
     &      (MOD(LIST_S(st_freq_code,NRECS),PHENOL_PERIOD) /= 0) THEN
              WRITE(CMESSAGE2,*)                                        &
     &       'INCORRECT SAMPLING FOR PHENOL_PERIOD. FREQ=',             &
     &        LIST_S(st_freq_code,NRECS), ':IGNORED.',                  &
     &       '(M,S,I) ',                                                &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-250
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
          ELSE IF(ITIMA == 15) THEN
            IF (LIST_S(st_freq_code,NRECS) == 1) THEN
              LIST_S(st_freq_code,NRECS)=TRIFFID_PERIOD
            ELSE IF                                                     &
     &      (MOD(LIST_S(st_freq_code,NRECS),TRIFFID_PERIOD) /= 0) THEN
              WRITE(CMESSAGE2,*)                                        &
     &       'INCORRECT SAMPLING FOR TRIFFID_PERIOD. FREQ=',            &
     &        LIST_S(st_freq_code,NRECS), ':IGNORED.',                  &
     &       '(M,S,I) ',                                                &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-260
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
          END IF

    ! For the NRECS item an end_time_code needs to be set if we
    ! are dealing with a times table rather than  regular diagn.
    ! This should be the maximum timestep in the time list. The list
    ! should be ready sorted (and thus maximum is last member) but
    ! will run through and find maximum to be on the safe side.

          IMAX = 0
          ITIMLST = -LIST_S(st_freq_code,NRECS+1)

          IF (ITIMLST  >   0 .and. ITIMLST /= 9999 ) THEN      ! List *not* regular
            DO I = 1, ITIM_T(ITIMLST)
              IF (IMAX  <   ITIM_S(I, ITIMLST)) THEN
                IMAX = ITIM_S(I, ITIMLST)
              END IF
            END DO

            LIST_S(st_end_time_code,NRECS) = IMAX

          END IF

!   Period

          IF ((INTV_T(ITIM_L) == -1).AND.(ITYP_T(ITIM_L) == 2)) THEN
            LIST_S(st_period_code,NRECS)=-1
          ELSE
            LIST_S(st_period_code,NRECS)=                               &
! DEPENDS ON: totimp
     &      TOTIMP(INTV_T(ITIM_L),UNT1_T(ITIM_L),MODL_L)
            if (LIST_S(st_period_code,NRECS)  ==  -999) then
               ErrorStatus = 106
               write (cmessage,'(a,a,i2,a,i3,a,i3)')                    &
     &         'PRELIM:TOTIMP:Error in time period conversion',         &
     &         ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
               write(6,*) cmessage
               GOTO 9999
            endif
          END IF

!   Add the record - unless the output destination is the dump,
!                      and output at the accumulating period
          IF (    LOCN_U(IUSE_L) >  2                                   &
     &      .OR.                                                        &
     &       ( (LIST_S(st_freq_code  ,NRECS+1) /=                       &
     &          LIST_S(st_period_code,NRECS  ))                         &
     &                                        .AND.                     &
     &         (LIST_S(st_start_time_code,NRECS+1) /=                   &
     &          LIST_S(st_end_time_code  ,NRECS+1))   )                 &
     &        )THEN
! No tag for parent, except for dump frequency means
             IF (UNT1_T(ITIM_L)=="DU" .AND. UNT1_T(ITIM_L)=="DU") THEN
                ! Special tag for dump frequency means
                LIST_S(st_macrotag,NRECS)=-999
             ELSE
                LIST_S(st_macrotag,NRECS)=0
             ENDIF
            NRECS=NRECS+1
          END IF

        ELSE IF (LIST_S(st_proc_no_code,NRECS) == 8) THEN
! Option of "daily" mean timeseries

          IF(NRECS >= NRECDP) THEN
            WRITE(CMESSAGE2,*)                                          &
     &     'TOO MANY S_LIST REQUESTS. REQUEST IGNORED',                 &
     &     '(M,S,I) ',                                                  &
     &      MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-270
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            NRECS=NRECS-1
            GOTO 999
          END IF

! Special case where 2 extra records required
!  Record 1 - time mean only no spatial processing
!  Record 2 - timeseries formed extracting from record 1
!  Record 3 - extract timeseries from dump ie record 2

          DO I=1,NELEMP+1          ! Copy stash list forward
            LIST_S(I,NRECS+1)=LIST_S(I,NRECS)
            LIST_S(I,NRECS+2)=LIST_S(I,NRECS)
          END DO

          IF(LOFFSET) THEN         ! Rad or conv timesteps,
                                   !       1 already added
            LIST_S(st_start_time_code,NRECS+2)=                         &
     &      LIST_S(st_start_time_code,NRECS+2)-1
            IF (LIST_S(st_period_code,NRECS) /= -1) THEN
              LIST_S(st_start_time_code,NRECS)=                         &
     &        LIST_S(st_start_time_code,NRECS)-                         &
     &        LIST_S(st_period_code,NRECS)
            ELSE
              LIST_S(st_start_time_code,NRECS)=1
            END IF

          ELSE

            IF (LIST_S(st_period_code,NRECS) /= -1) THEN
              LIST_S(st_start_time_code,NRECS)=                         &
     &        LIST_S(st_start_time_code,NRECS)-                         &
     &        LIST_S(st_period_code,NRECS)+1
            ELSE
              LIST_S(st_start_time_code,NRECS)=1
            END IF

          END IF

          IF(LIST_S(st_start_time_code,NRECS) <  1) THEN
            WRITE(CMESSAGE2,*)                                          &
     &       'START TIME BEFORE PERIOD, SETTING TO 1',                  &
     &       '(M,S,I) ',MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-280
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            LIST_S(st_start_time_code,NRECS)=1
          END IF

          LIST_S(st_proc_no_code ,NRECS)=3    ! time mean
          LIST_S(st_proc_no_code ,NRECS+1)=8  ! timseries special case
          LIST_S(st_proc_no_code ,NRECS+2)=1  !  extract

! Reset first record to no area weight or spatial processing
! ie first record just controls time meaning

          LIST_S(st_gridpoint_code,NRECS)=1
          LIST_S(st_weight_code,NRECS)=0

          LIST_S(st_input_bottom ,NRECS+1)=                             &
     &    LIST_S(st_output_bottom,NRECS  )
          LIST_S(st_input_bottom ,NRECS+2)=                             &
     &    LIST_S(st_output_bottom,NRECS+1)

          LIST_S(st_input_top    ,NRECS+1)=                             &
     &    LIST_S(st_output_top   ,NRECS  )
          LIST_S(st_input_top    ,NRECS+2)=                             &
     &    LIST_S(st_output_top   ,NRECS+1)

          LIST_S(st_input_code   ,NRECS+1)=-NRECS
          LIST_S(st_input_code   ,NRECS+2)=-NRECS-1
          LIST_S(st_output_code  ,NRECS  )=1
          LIST_S(st_output_code  ,NRECS+1)=1  ! dump
          LIST_S(st_series_ptr   ,NRECS+2)=0
          LIST_S(st_series_ptr   ,NRECS)=0
          LIST_S(NELEMP+1        ,NRECS+1)=NRECS+1
          LIST_S(NELEMP+1        ,NRECS+2)=NRECS+2
!  definition 8 implies frequency of time mean over every timestep
          LIST_S(st_freq_code,NRECS)=1
          LIST_S(st_freq_code,NRECS+1)=                                 &
                                                      ! Frequency
! DEPENDS ON: totimp
     &    TOTIMP(ISAM_T(ITIM_L),UNT2_T(ITIM_L),MODL_L)
          IF (LIST_S(st_freq_code,NRECS)  ==  -999) THEN
             ErrorStatus = 107
             write (cmessage,'(a,a,i2,a,i3,a,i3)')                      &
     &       'PRELIM:TOTIMP:Error in time period conversion',           &
     &       ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
             GOTO 9999
          ENDIF


!   Correct frequency for radiation, periodic convection, leaf
!   phenology and vegetation competition

          IF (ITIMA == 2) THEN
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
              LIST_S(st_freq_code,NRECS)=A_LW_RADSTEP_DIAG
#else
              LIST_S(st_freq_code,NRECS)=A_LW_RADSTEP
#endif
          ELSE IF(ITIMA == 3) THEN
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
              LIST_S(st_freq_code,NRECS)=A_SW_RADSTEP_DIAG
#else
              LIST_S(st_freq_code,NRECS)=A_SW_RADSTEP
#endif
          ELSE IF(ITIMA == 13) THEN
              LIST_S(st_freq_code,NRECS)=A_CONV_STEP
          ELSE IF(ITIMA == 14) THEN
              LIST_S(st_freq_code,NRECS)=PHENOL_PERIOD
          ELSE IF(ITIMA == 15) THEN
              LIST_S(st_freq_code,NRECS)=TRIFFID_PERIOD
          END IF

!   Period
! time mean over sampling period
            LIST_S(st_period_code,NRECS)=                               &
! DEPENDS ON: totimp
     &      TOTIMP(ISAM_T(ITIM_L),UNT2_T(ITIM_L),MODL_L)
            IF (LIST_S(st_period_code,NRECS)  ==  -999) THEN
               ErrorStatus = 108
               write (cmessage,'(a,a,i2,a,i3,a,i3)')                    &
     &         'PRELIM:TOTIMP:Error in time period conversion',         &
     &         ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
               GOTO 9999
            ENDIF

! period for timeseries recycle period
            LIST_S(st_period_code,NRECS+1)=                             &
! DEPENDS ON: totimp
     &      TOTIMP(INTV_T(ITIM_L),UNT1_T(ITIM_L),MODL_L)
            IF (LIST_S(st_period_code,NRECS+1)  ==  -999) THEN
               ErrorStatus = 109
               write (cmessage,'(a,a,i2,a,i3,a,i3)')                    &
     &         'PRELIM:TOTIMP:Error in time period conversion',         &
     &         ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
               GOTO 9999
            ENDIF


! st_start_time for 2 record should be period for first record
! unless offset from start of run. Note value independent of logical
!  OFFSET
!
            IF (LIST_S(st_period_code,NRECS) /= -1) THEN
              IF (LOFFSET) THEN
               LIST_S(st_start_time_code,NRECS+1)=                      &
     &          LIST_S(st_start_time_code,NRECS+1) -                    &
     &          LIST_S(st_period_code,NRECS+1) +                        &
     &          LIST_S(st_freq_code,NRECS+1) - 1
              ELSE
               LIST_S(st_start_time_code,NRECS+1)=                      &
     &          LIST_S(st_start_time_code,NRECS+1) -                    &
     &          LIST_S(st_period_code,NRECS+1) +                        &
     &          LIST_S(st_freq_code,NRECS+1)
              ENDIF
            ELSE
              LIST_S(st_start_time_code,NRECS+1)=1
            END IF


!   Add both record
            LIST_S(st_macrotag,NRECS)=0
            NRECS=NRECS+2

        END IF       ! Other than single time field

      END IF         ! Diag request not null - ITIM_L /= 0
 999  CONTINUE
      END DO         ! Loop over diagnostic requests

      END IF         ! NDIAG >  0

! DEPENDS ON: pslcom
      CALL PSLCOM(NRECS)    ! Compress out unused pseudo levels lists

 9999 RETURN
      END SUBROUTINE PRELIM

!- End of Subroutine code -------------------------------------------


!+ Compress out unused pseudo levels lists

!- End of Subroutine code ------------------------------------------
#endif
