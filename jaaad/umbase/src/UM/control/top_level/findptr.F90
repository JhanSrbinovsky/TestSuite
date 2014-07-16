#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: FINDPTR  -------------------------------------------------
!LL
!LL  Purpose: Locates address within D1 of diagnostic field which may
!LL           be required elsewhere in the model for special-purpose
!LL           diagnostic routine such as zonal mean print, or as an
!LL           internal interfacing field for coupling sub-models.
!LL           The search information is input in STASH format, and the
!LL           STASH list is scanned for a match.  If the specified
!LL           field does not exist in D1 the address is returned as 0.
!LL           NB: Missing data indicators may be supplied if the search
!LL               is to ignore certain elements in the STASH list.
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL
!LL  Author:   T.C.Johns
!LL
!LL  Code version no: 1.3           Date: 04 March 1992
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL 3.4  16/6/94 : Change CHARACTER*(*) to CHARACTER*(80) N.Farnon
!LL 3.5  June 95   Submodels project:
!LL                Added internal_model to subroutine args.
!LL                Altered args in each STINDEX reference.
!LL                Changed hardwired addresses in each STLIST reference
!LL                  to parameter addresses as defined in STPARAM
!LL                Added *CALL STPARAM
!LL                  S.J.Swarbrick
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!LL  5.3  19/02/02  Remove argsize from argument list. D.Robinson
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: C610
!LL
!LL  Project task: C4?
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage Handling and
!LL                                 Diagnostic System.
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE FINDPTR ( internal_model,SECTION,ITEM,                 &
     &                     PROCESS_CODE,FREQ_CODE,START,END,PERIOD,     &
     &                     GRIDPT_CODE,WEIGHT_CODE,                     &
     &                     BOTTOM_LEVEL,TOP_LEVEL,                      &
     &                     GRID_N,GRID_S,GRID_W,GRID_E,                 &
     &                     STASHMACRO_TAG,MDI,ADDRESS,                  &
#include "argsts.h"
     &                     ICODE,CMESSAGE)
!
      IMPLICIT NONE
!
      INTEGER                                                           &
     &       internal_model,                                            &
                              ! IN  - internal_model id.
     &       SECTION,                                                   &
                                ! IN  - STASH section number
     &       ITEM,                                                      &
                                ! IN  - STASH item number
     &       PROCESS_CODE,                                              &
                                ! IN  - STASH processing code
     &       FREQ_CODE,                                                 &
                                ! IN  - STASH frequency code
     &       START,                                                     &
                                ! IN  - STASH start step for processing
     &       END,                                                       &
                                ! IN  - STASH end step for processing
     &       PERIOD,                                                    &
                                ! IN  - STASH processing period
     &       GRIDPT_CODE,                                               &
                                ! IN  - STASH gridpoint code
     &       WEIGHT_CODE,                                               &
                                ! IN  - STASH weighting code
     &       BOTTOM_LEVEL,                                              &
                                ! IN  - STASH input bottom level
     &       TOP_LEVEL,                                                 &
                                ! IN  - STASH input top level
     &       GRID_N,                                                    &
                                ! IN  - STASH N-row grid code
     &       GRID_S,                                                    &
                                ! IN  - STASH S-row grid code
     &       GRID_W,                                                    &
                                ! IN  - STASH W-col grid code
     &       GRID_E,                                                    &
                                ! IN  - STASH E-col grid code
     &       STASHMACRO_TAG,                                            &
                                ! IN  - STASHmacro tag number
     &       MDI,                                                       &
                                ! IN  - Missing Data Indicator
     &       ADDRESS,                                                   &
                                ! OUT - Address in D1
     &       ICODE              ! OUT - Error return code
      CHARACTER*(80)                                                    &
     &       CMESSAGE           ! OUT - Error return message
!
!*----------------------------------------------------------------------
!  Common blocks
!
#include "csubmodl.h"
#include "parparm.h"
#include "typsize.h"
#include "typsts.h"
#include "stparam.h"
!
!  Local variables
!
      INTEGER                                                           &
     &    ISTART,IEND,I,                                                &
                                  ! Start, end + loop index in STASHlist
     &    NMATCH                                                        &
                                  ! Number of matches found
     &    ,im_index               ! Internal model index
      LOGICAL                                                           &
     &    MATCH                   ! TRUE if diagnostic matched
!L----------------------------------------------------------------------
!L 0.  Check that tag field is within the allowed range for user tags
!L
      IF (STASHMACRO_TAG /= MDI .AND.                                   &
     &   (STASHMACRO_TAG <  0 .OR. STASHMACRO_TAG >  999)) THEN
        CMESSAGE="FINDPTR : STASHMACRO_TAG must be in range 0-999"
        ICODE=ABS(STASHMACRO_TAG)
        GOTO 999
      ENDIF
!L----------------------------------------------------------------------
!L 1.  Locate start/end limits within STASHlist for search;
!L     initialise output ADDRESS to zero
!L
      ADDRESS=0
      NMATCH=0
      im_index=internal_model_index(internal_model)
      IF (STINDEX(2,ITEM,SECTION,im_index) >  0) THEN
        ISTART=STINDEX(1,ITEM,SECTION,im_index)
        IEND  =STINDEX(2,ITEM,SECTION,im_index)+ISTART-1
!L
!L 1.1 Loop over STASHlist entries and try to find matches
!L
        DO I=ISTART,IEND
          IF (STLIST(s_modl,I) /= internal_model.OR.                    &
     &        STLIST(s_sect,I) /= SECTION.OR.                           &
     &        STLIST(s_item,I) /= ITEM) THEN
            ICODE=1000*SECTION+ITEM
            CMESSAGE="FINDPTR : Corrupt STASHlist or STASHindex"
            GOTO 999
          ENDIF
          MATCH=((STLIST(s_output,I) == 1).OR.                          &
     &           (STLIST(s_output,I) == 2)).AND.                        &
     &    (PROCESS_CODE == STLIST(s_proc,I)                             &
     &                                .OR.PROCESS_CODE == MDI).AND.     &
     &    (FREQ_CODE   == STLIST( s_freq,I)                             &
     &                                .OR.FREQ_CODE   == MDI) .AND.     &
     &    (START       == STLIST( s_times,I)                            &
     &                                .OR.START       == MDI) .AND.     &
     &    (END         == STLIST( s_timee,I)                            &
     &                                .OR.END         == MDI) .AND.     &
     &    (PERIOD      == STLIST( s_period,I)                           &
     &                                .OR.PERIOD      == MDI) .AND.     &
     &    (GRIDPT_CODE == STLIST( s_grid,I)                             &
     &                                .OR.GRIDPT_CODE == MDI) .AND.     &
     &    (WEIGHT_CODE == STLIST( s_weight,I)                           &
     &                                .OR.WEIGHT_CODE == MDI) .AND.     &
     &    (BOTTOM_LEVEL == STLIST(s_bottom,I)                           &
     &                                .OR.BOTTOM_LEVEL == MDI) .AND.    &
     &    (TOP_LEVEL   == STLIST(s_top,I)                               &
     &                                .OR.TOP_LEVEL   == MDI) .AND.     &
     &    (GRID_N      == STLIST(s_north,I)                             &
     &                                .OR.GRID_N      == MDI) .AND.     &
     &    (GRID_S      == STLIST(s_south,I)                             &
     &                                .OR.GRID_S      == MDI) .AND.     &
     &    (GRID_W      == STLIST(s_west,I)                              &
     &                                .OR.GRID_W      == MDI) .AND.     &
     &    (GRID_E      == STLIST(s_east,I)                              &
     &                                .OR.GRID_E      == MDI) .AND.     &
     &    (STASHMACRO_TAG == MOD(STLIST(st_macrotag,I),1000).OR.        &
     &     STASHMACRO_TAG == MDI)
!
          IF (MATCH) THEN
            ADDRESS=STLIST(st_output_addr,I)
            NMATCH=NMATCH+1
          ENDIF
        ENDDO
!
        IF (NMATCH >  1) THEN
          ICODE=-1000*SECTION-ITEM
          CMESSAGE="FINDPTR : Warning - multiple match for diagnostic"
      WRITE(6,*)"FINDPTR : Warning - multiple match for diagnostic ",   &
     &            SECTION,ITEM
!
        ENDIF
      ENDIF
!
 999  CONTINUE
      RETURN
!L----------------------------------------------------------------------
      END SUBROUTINE FINDPTR
#endif
