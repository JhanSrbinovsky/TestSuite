#if defined(C70_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine ADDRESS_CHECK --------------------------------------
!LL
!LL Purpose : Check that start addresses of fields read in agree
!LL           with start addresses set up by UI. Called in INITDUMP
!LL           if prognostic fields read in from atmos or ocean dumps.
!LL
!LL  Model   Date     Modification history:
!LL version
!LL   3.2  25/05/93  New routine. D Robinson
!LL 3.4  16/6/94 : Change CHARACTER*(*) to CHARACTER*(80) N.Farnon
!LL   3.5   May 95   Submodels project. Inserted *CALL CSUBMODL,
!LL                  *CALL CPPXREF, *CALL ARGPPX, *CALL PPXLOOK
!LL                  to pass ppxref lookup arrays to F_TYPE.
!LL                  Modified reference to SI to take account of
!LL                    submodel partitioning.
!LL                  S.J.Swarbrick
!LL   4.0   05/01/96 Get Internal identifier from LOOKUP(45) to
!LL                  determine im_index for SI array. D. Robinson
!LL   4.1   21/03/96 MPP code : Added MPP_DUMP_ADDR/LEN arguments for
!LL                  use when checking dump addressing against stash
!LL                  addressing.   P.Burton
!LL   4.1   23/05/96 Remove internal_model from argument list.
!LL                  D. Robinson
!LL  4.4  01/09/97  Add helpful message after consistency checks. RTHB.
!    6.0  19/06/03  Remove non-MPP parts of code. T. White
!LL
!LL Programming Standard : UM documentation paper no. 3
!LL                        version no. 1, dated 15/01/90
!LL
!LL Documentation : None
!LL
!LLEND--------------------------------------------------------------
!
!*L Arguments

      SUBROUTINE ADDRESS_CHECK (LOOKUP,MPP_DUMP_ADDR,MPP_DUMP_LEN,      &
     &                          LEN1_LOOKUP,LEN2_LOOKUP,                &
     &                          SI,NITEMS,NSECTS,LEN_DATA,              &
#include "argppx.h"
     &                          ICODE,CMESSAGE)

      IMPLICIT NONE

#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"

      INTEGER                                                           &
     &    LEN1_LOOKUP                                                   &
                              !  1st dimension of lookup table
     &   ,LEN2_LOOKUP                                                   &
                              !  2nd dimension of lookup table
     &   ,LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                               &
                                           !  Lookup table
     &   ,MPP_DUMP_ADDR(LEN2_LOOKUP)                                    &
                                      ! Addresses of fields as
!                                     ! calculated in READDUMP.
     &   ,MPP_DUMP_LEN(LEN2_LOOKUP)                                     &
                                      ! Lengths of fields as
!                                     ! calculated in READDUMP.

     &   ,LEN_DATA                                                      &
                              !  Expected length of data
     &   ,NITEMS                                                        &
                              !  No of stash items
     &   ,NSECTS                                                        &
                              !  No of stash sections
!  Stash item addresses
     &   ,SI(NITEMS,0:NSECTS,N_INTERNAL_MODEL)                          &
     &   ,ICODE               !  Return code

      CHARACTER*(80)                                                    &
     &          CMESSAGE      !  Error message

!L Dynamic allocation of arrays for F_TYPE
      INTEGER                                                           &
     &    PP_NUM   (LEN2_LOOKUP)                                        &
     &   ,PP_LEN   (LEN2_LOOKUP)                                        &
     &   ,PP_POS   (LEN2_LOOKUP)                                        &
     &   ,PP_LS    (LEN2_LOOKUP)                                        &
     &   ,PP_STASH (LEN2_LOOKUP)                                        &
     &   ,PP_TYPE  (LEN2_LOOKUP)

!L Local array
      INTEGER FIXHD(256)  !  Dummy array until removed from F_TYPE

!L Local variables
      INTEGER                                                           &
     &    ADDRESS_STASH                                                 &
     &   ,ADDRESS_LOOKUP                                                &
     &   ,ITEM_CODE                                                     &
     &   ,J                                                             &
     &   ,LEN                                                           &
     &   ,N_TYPES                                                       &
     &   ,SECT_NO                                                       &
     &   ,im_ident                                                      &
                   ! Internal model identifier
     &   ,im_index                                                      &
                   ! Position of int mod id in INTERNAL_MODEL_LIST
     &,OLD_STASH   ! VALUE OF STASH NUMBER ON PREVIOUS ITERATION OF LOOP

      CHARACTER*80 TITLE

!L Subroutines called
      EXTERNAL F_TYPE

!*---------------------------------------------------------------------
!
!   SET INITIAL VALUE OF PREVIOUS STASH NUMBER
!
      OLD_STASH = -1
!
!L   Internal Structure

      TITLE = 'Prognostic fields'
!     modify f_type later to add len1_lookup and remove fixhd
! DEPENDS ON: f_type
      CALL F_TYPE (LOOKUP,LEN2_LOOKUP,PP_NUM,N_TYPES,PP_LEN,            &
     &             PP_STASH,PP_TYPE,PP_POS,PP_LS,FIXHD,                 &
#include "argppx.h"
     &TITLE)

      DO J=1,N_TYPES

!       Get Stash Section no and Item Code
        ITEM_CODE = MOD ( PP_STASH(J),1000)
        SECT_NO   = (PP_STASH(J)-ITEM_CODE)/1000

!       Get im_ident/index for this field
        im_ident = LOOKUP(45,pp_pos(j))
        im_index = INTERNAL_MODEL_INDEX(im_ident)

!       Get lookup and stash start address
        ADDRESS_LOOKUP = MPP_DUMP_ADDR(PP_POS(J))
        ADDRESS_STASH  = SI(ITEM_CODE,SECT_NO,im_index)

!       Check that they match
!
!            CHECK THAT START ADDRESSES AGREE FOR FIRST OCCURRANCE
!            OF A NEW STASH CODE:  FOR FIXED LENGTH FIELDS THERE
!            IS ONLY ONE ENTRY IN THE PP_STASH ARRAY FOR EACH
!            STASH CODE, BUT FOR PACKED FIELDS (EG OCEAN)
!            EACH LEVEL MIGHT HAVE A DIFFERENT LENGTH AND
!            GENERATE A NEW PP_STASH VALUE.
!
         IF (ADDRESS_STASH  /=  ADDRESS_LOOKUP .AND.                    &
     &       OLD_STASH  /=  PP_STASH(J) ) THEN
          CMESSAGE = 'ADDR_CHK : Mis_match in start addresses'
          WRITE (6,*) ' Stash Sect No ',SECT_NO,' Item No ',ITEM_CODE
          WRITE (6,*) ' Start Address in SI           ',ADDRESS_STASH
          WRITE (6,*) ' Start Address in LOOKUP Table ',ADDRESS_LOOKUP
          WRITE (6,*) ' You probably need to RECONFIGURE the start dump'
          ICODE = 1
          GO TO 999   !  Return
        ENDIF

!     REMEMBER CURRENT VERSOIN OF PP_STASH FOR NEXT TIME THRU LOOP
         OLD_STASH = PP_STASH(J)
!
      ENDDO

!     Check full length
      LEN = 0
      DO J=1,LEN2_LOOKUP
        LEN = LEN + MPP_DUMP_LEN(J)
      ENDDO

      IF (LEN  /=  LEN_DATA) THEN
        CMESSAGE = 'ADDR_CHK : Mismatch in length of data'
        WRITE (6,*) ' Length according to LOOKUP table ',LEN
        WRITE (6,*) ' Length set up in D1 array        ',LEN_DATA
        WRITE (6,*) ' You probably need to RECONFIGURE the start dump'
        ICODE = 2
        GO TO 999   !  Return
      ENDIF

 999  RETURN
      END SUBROUTINE ADDRESS_CHECK
#endif
