#if defined(A31_1A)
#if defined(ATMOS) && !defined(GLOBAL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine CHK_LOOK_BOUNDA
!
! Purpose : Cross checks values in LOOKUP records of boundary data
!           with model run values
!
! Author :  P.Burton
!
! Model             Modification history from model version 5.2
! version    Date
! 5.2      15/09/00 Deck rewritten for 5.x format LBCs      P.Burton
! 5.3      07/11/01 Remove argsize from argument list.  Z.Gardner
!
! ---------------------------------------------------------

      SUBROUTINE CHK_LOOK_BOUNDA(                                       &
     &  ITEM_LIST,FULL_LOOKUP_BOUNDA,                                   &
#include "argbnd.h"
#include "argppx.h"
     &                           ICODE,CMESSAGE)

      IMPLICIT NONE

#include "cmaxsize.h"
#include "parvars.h"
#include "typsize.h"
#include "typbnd.h"
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"

      INTEGER                                                           &
     &  ITEM_LIST(RIM_LOOKUPSA)                                         &
                                   ! IN: STASH codes of expected items
     &, FULL_LOOKUP_BOUNDA(LEN1_LOOKUP,BOUND_LOOKUPSA)                  &
                                   ! IN: Full LOOKUP record for LBCs
     &, ICODE                     ! OUT : Return code

      CHARACTER*(80)                                                    &
     &  CMESSAGE                  ! OUT : Error message

#include "clookadd.h"
#include "stparam.h"

! Functions
      INTEGER                                                           &
     &  EXPPXI                                                          &
     &, GET_FLD_TYPE

! Local variables
      INTEGER                                                           &
     &  variable                                                        &
                          ! Loop counter for variable
     &, code                                                            &
                          ! item_code value for variable
     &, model                                                           &
                          ! model number for variable
     &, section                                                         &
                          ! section number for variable
     &, item                                                            &
                          ! section number for variable
     &, grid_type                                                       &
                          ! grid code for variable
     &, fld_type                                                        &
                          ! P,U or V for variable
     &, halo_type                                                       &
                          ! halo type for variable
     &, level_type                                                      &
                          ! what type of level for variable
     &, bottom_level_code                                               &
                          ! bottom level code
     &, bottom_level                                                    &
                          ! bottom level for variable
     &, top_level_code                                                  &
                          ! top level code
     &, top_level                                                       &
                          ! top level for variable
     &, n_levels_expected                                               &
                          ! number of levels expected
     &, n_levels_lbc                                                    &
                          ! number of levels in file
     &, halo_x_expected                                                 &
                          ! expected size of halo in x
     &, halo_x_lbc                                                      &
                          ! actual size of halo in x
     &, halo_y_expected                                                 &
                          ! expected size of halo in y
     &, halo_y_lbc                                                      &
                          ! actual size of halo in y
     &, size_x_expected                                                 &
                          ! expected size of field in x
     &, size_x                                                          &
                          ! actual size of field in x
     &, size_y_expected                                                 &
                          ! expected size of field in y
     &, size_y                                                          &
                          ! actual size of y
     &, rim_type                                                        &
                          ! Type of RIMWIDTH
     &, rimwidth_expected                                               &
                          ! expected rimwidth
     &, rimwidth_lbc                                                    &
                          ! actual rimwidth
     &, size_expected     ! Expected size

!=====================================================================

! 1.0 Check that the first field in the file is the orography field

      IF (LOOKUP_BOUNDA(ITEM_CODE,1)  /=  31001) THEN ! Orography
        ICODE=1
        CMESSAGE='CHK_LOOK_BOUNDA : No orography in LBC file'
        GOTO 9999
      ENDIF

! 2.0 Check for the expected number of variables for each time.
!     Bear in mind that RIM_LOOKUPSA contains an extra field
!     (orography) which only occurs at the start of the field.
!     So the actual number of fields written out at each LBC output
!     time is actually RIM_LOOKUPSA-1

      IF (FULL_LOOKUP_BOUNDA(ITEM_CODE,2)  /=                           &
     &    FULL_LOOKUP_BOUNDA(ITEM_CODE,2+RIM_LOOKUPSA-1)) THEN
        WRITE(6,*) 'Wrong number of LBC variables found in LBC file'
        WRITE(6,*) 'Expecting record ',2+RIM_LOOKUPSA-1,' to contain ', &
     &             'STASH item code ',FULL_LOOKUP_BOUNDA(ITEM_CODE,2)
        WRITE(6,*) 'But found item code ',                              &
     &             FULL_LOOKUP_BOUNDA(ITEM_CODE,2+RIM_LOOKUPSA-1)
        ICODE=2
        CMESSAGE='CHK_LOOK_BOUNDA : Wrong number of LBC fields'
        GOTO 9999
      ENDIF

! 3.0 Now check the header record for each required variable

      DO variable=1,RIM_LOOKUPSA

        code=ITEM_LIST(variable)  ! item_code value for variable

        IF (LOOKUP_BOUNDA(ITEM_CODE,variable)  /=  code) THEN
          WRITE(6,*) 'Unexpected field in LBC file'
          WRITE(6,*) 'Field ',variable,' was expected to be ',          &
     &               code,' but found ',                                &
     &               LOOKUP_BOUNDA(ITEM_CODE,variable)

! DEPENDS ON: pr_look
          CALL PR_LOOK(                                                 &
#include "argppx.h"
     &                 LOOKUP_BOUNDA,LOOKUP_BOUNDA,LEN1_LOOKUP,variable)

          ICODE=3
          CMESSAGE='CHK_LOOK_BOUNDA : Unexpected field in LBC file'
          GOTO 9999
        ENDIF

        model=LOOKUP_BOUNDA(MODEL_CODE,variable)
        item=MOD(LOOKUP_BOUNDA(ITEM_CODE,variable),1000)
        section=(LOOKUP_BOUNDA(ITEM_CODE,variable)-item)/1000

! DEPENDS ON: exppxi
        grid_type=EXPPXI(model,section,item,ppx_grid_type,              &
#include "argppx.h"
     &                   ICODE, CMESSAGE)
! DEPENDS ON: get_fld_type
        fld_type=GET_FLD_TYPE(grid_type)

! DEPENDS ON: exppxi
        halo_type=EXPPXI(model,section,item,ppx_halo_type,              &
#include "argppx.h"
     &                   ICODE, CMESSAGE)

! DEPENDS ON: exppxi
        level_type=EXPPXI(model,section,item,ppx_lv_code,               &
#include "argppx.h"
     &                    ICODE, CMESSAGE)
        IF (level_type  ==  5) THEN
          n_levels_expected=1
        ELSE
! DEPENDS ON: exppxi
          bottom_level_code=EXPPXI(model,section,item,ppx_lb_code,      &
#include "argppx.h"
     &                          ICODE, CMESSAGE)
! DEPENDS ON: exppxi
          top_level_code=EXPPXI(model,section,item,ppx_lt_code,         &
#include "argppx.h"
     &                          ICODE, CMESSAGE)
! DEPENDS ON: levcod
          CALL LEVCOD(bottom_level_code,bottom_level,ICODE,CMESSAGE)
! DEPENDS ON: levcod
          CALL LEVCOD(top_level_code,top_level,ICODE,CMESSAGE)
          n_levels_expected=top_level-bottom_level+1
        ENDIF

        halo_x_expected=halosize(1,halo_type)
        halo_y_expected=halosize(2,halo_type)
        size_x_expected=glsize(1,fld_type)
        size_y_expected=glsize(2,fld_type)

        IF (fld_type  ==  fld_type_u) THEN
          size_x_expected=size_x_expected-1
        ENDIF

        IF (LOOKUP_BOUNDA(ITEM_CODE,variable)  ==  31001) THEN
          ! Orography
          rim_type=rima_type_orog
        ELSE
          rim_type=rima_type_norm
        ENDIF

        rimwidth_expected=RIMWIDTHA(rim_type)

        halo_x_lbc=MOD(LOOKUP_BOUNDA(LBUSER3,variable),100)
        halo_y_lbc=MOD(LOOKUP_BOUNDA(LBUSER3,variable)-halo_x_lbc,      &
     &                 10000)/100
        rimwidth_lbc=MOD((LOOKUP_BOUNDA(LBUSER3,variable)-              &
     &                 halo_x_lbc-halo_y_lbc*100),1000000)/10000

        n_levels_lbc=LOOKUP_BOUNDA(LBHEM,variable)-100

        size_expected=global_LENRIMA(fld_type,halo_type,rim_type)*      &
     &                n_levels_expected

        IF (n_levels_lbc  /=  n_levels_expected) THEN
          WRITE(6,*) 'Wrong number of levels for LBC field ',variable
          WRITE(6,*) 'Expected ',n_levels_expected,' levels but found ',&
     &               n_levels_lbc

! DEPENDS ON: pr_look
          CALL PR_LOOK(                                                 &
#include "argppx.h"
     &                 LOOKUP_BOUNDA,LOOKUP_BOUNDA,LEN1_LOOKUP,variable)

          ICODE=4
          CMESSAGE='CHK_LOOK_BOUNDA : Wrong number of levels'
          GOTO 9999
        ENDIF

        IF ((halo_x_lbc  /=  halo_x_expected) .OR.                      &
     &      (halo_y_lbc  /=  halo_y_expected)) THEN
          WRITE(6,*) 'Incorrect halos for LBC field ',variable
          WRITE(6,*) 'Expected halo_x= ',halo_x_expected,               &
     &               ' and halo_y= ',halo_y_expected
          WRITE(6,*) 'but found halo_x= ',halo_x_lbc,                   &
     &               ' and halo_y= ',halo_y_lbc

! DEPENDS ON: pr_look
          CALL PR_LOOK(                                                 &
#include "argppx.h"
     &                 LOOKUP_BOUNDA,LOOKUP_BOUNDA,LEN1_LOOKUP,variable)

          ICODE=5
          CMESSAGE='CHK_LOOK_BOUNDA : Incorrect halos'
          GOTO 9999
        ENDIF

        IF (rimwidth_lbc  /=  rimwidth_expected) THEN
          WRITE(6,*) 'Wrong RIMWIDTH for LBC field ',variable
          WRITE(6,*) 'Expected RIMWIDTH= ',rimwidth_expected,           &
     &               'but found RIMWIDTH= ',rimwidth_lbc

! DEPENDS ON: pr_look
          CALL PR_LOOK(                                                 &
#include "argppx.h"
     &                 LOOKUP_BOUNDA,LOOKUP_BOUNDA,LEN1_LOOKUP,variable)

          ICODE=6
          CMESSAGE='CHK_LOOK_BOUNDA : Wrong RIMWIDTH'
          GOTO 9999
        ENDIF

        size_x=LOOKUP_BOUNDA(LBNPT,variable)
        size_y=LOOKUP_BOUNDA(LBROW,variable)

        IF ((size_x  /=                                                 &
     &       size_x_expected) .OR.                                      &
     &      (size_y  /=                                                 &
     &       size_y_expected)) THEN
          WRITE(6,*) 'Incorrect dimensions for LBC field ',variable
          WRITE(6,*) 'Expected ROW_LENGTH= ',size_x_expected,' and ',   &
     &               'ROWS= ',size_y_expected
          WRITE(6,*) 'But found ROW_LENGTH= ',                          &
     &               size_x,' and ROWS= ',size_y

! DEPENDS ON: pr_look
          CALL PR_LOOK(                                                 &
#include "argppx.h"
     &                 LOOKUP_BOUNDA,LOOKUP_BOUNDA,LEN1_LOOKUP,variable)

          ICODE=7
          CMESSAGE='CHK_LOOK_BOUNDA : Wrong dimensions'
          GOTO 9999
        ENDIF

        IF (LOOKUP_BOUNDA(LBLREC,variable)  /=  size_expected) THEN
          WRITE(6,*) 'Wrong size for LBC field ',variable
          WRITE(6,*) 'Expected size was ',size_expected,' but found ',  &
     &               LOOKUP_BOUNDA(LBLREC,variable)

! DEPENDS ON: pr_look
          CALL PR_LOOK(                                                 &
#include "argppx.h"
     &                 LOOKUP_BOUNDA,LOOKUP_BOUNDA,LEN1_LOOKUP,variable)

          ICODE=8
          CMESSAGE='CHK_LOOK_BOUNDA : Wrong size'
          GOTO 9999
        ENDIF

      ENDDO ! variable

 9999 CONTINUE

      RETURN
      END SUBROUTINE CHK_LOOK_BOUNDA
#endif
#endif
