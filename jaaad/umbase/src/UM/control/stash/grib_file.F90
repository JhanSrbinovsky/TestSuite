#if defined(C84_1A) || defined(FLDOP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: GRBWRT------------------------------------------------
!LL
!LL  Purpose: This routine acts as an interface between the model and
!LL  GRIB format output routines.
!LL
!LL  Author:   D.M.Goddard        Date:           23 December 1993
!LL  Reviewer:                    Date of review:
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Code version no: 1           Date: 15 October 1993
!LL
!LL  Modification History:
!LL  3.4   11/10/94 : Correct setting of reals in lookup table
!LL                   and add return code and message to PP2GRIB call
!LL                   R A Stratton.
!    4.0   10/03/95 : Allow alternative grib packing to be used and
!                     improve error traping. R A Stratton.
!LL  4.3   06/02/97  Modify I/O calls for MPP use  P.Burton
!LL  5.5   25/04/03  Grib data format not supported on non-CRAY
!LL                  platform                           E.Leung
!    6.1   18/08/04  Re-enable for FLDOP on non-CRAY platforms.
!                                                       P.Dando
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: ...
!LL
!LL  Project task: ...
!LL
!LL  External documentation: On-line UM document ??? - ??????????
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
      SUBROUTINE GRIB_FILE(LEN1_LOOKUP,LEN2_LOOKUP,LOOKUP,RLOOKUP,IENT, &
     &                     FIELD,PPHORIZ_OUT,LENBUF,NUM_CRAY_WORDS,     &
     &                     UNITPP,IWA,GRIB_PACKING,ICODE,CMESSAGE)

      INTEGER                                                           &
     &     LEN1_LOOKUP                                                  &
                       !  IN   first dimension of LOOKUP
     &    ,LEN2_LOOKUP                                                  &
                       !  IN   second dimension of LOOKUP
     &    ,LENBUF                                                       &
                       !  IN   No of points in output field
     &    ,IENT                                                         &
                       !  IN   level indicator for processing LOOKUP.
     &    ,IWA                                                          &
                       !  IN   Record number
     &    ,PPHORIZ_OUT                                                  &
                       !  IN
     &    ,UNITPP                                                       &
                       !  IN   Output PP unit number
     &    ,GRIB_PACKING                                                 &
                        !  IN  Packing profile for grib
     &    ,LEN_FIELD                                                    &
     &    ,ICODE                                                        &
                          !  OUT  Return code
     &    ,NUM_CRAY_WORDS                                               &
                          !  OUT  Number of cray words output in grib
     &    ,LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP) ! Integer lookup headers
      REAL                                                              &
     &     FIELD(PPHORIZ_OUT)                                           &
                               ! IN   Unpacked output array
     &    ,RLOOKUP(LEN1_LOOKUP,LEN2_LOOKUP) ! REAL lookup headers
      CHARACTER                                                         &
     &     CMESSAGE*(*)     ! OUT  Will contain any error messages
!
! LOCAL VARIABLES
!
      INTEGER                                                           &
     &     ILABEL(45)                                                   &
                            ! Integer part of LOOKUP for level IENT
     &    ,LEN_IO                                                       &
     &    ,IX
      REAL                                                              &
     &     RLABEL(19)                                                   &
                            ! Real part of LOOKUP for level IENT
     &    ,WORK_ARRAY(LENBUF)                                           &
                              ! GRIB packed output array
     &    ,BUFOUT(LENBUF)   ! Output PP BUFFER
#include "clookadd.h"
#include "c_mdi.h"

!L
!L 1. Fill arrays ILABEL and RLABEL
!L
      DO J=1,45
        ILABEL(J)=LOOKUP(J,IENT)
      ENDDO
      DO J=1,19
        RLABEL(J)=RLOOKUP(J+45,IENT)
      ENDDO
!L
!L 2. Convert data to GRIB code
!L
! DEPENDS ON: pp2grib
      CALL PP2GRIB(FIELD,WORK_ARRAY,LENBUF,NUM_CRAY_WORDS,GRIB_PACKING, &
     &             ILABEL,RLABEL,ICODE,CMESSAGE)
      IF (ICODE /= 0) THEN
        RETURN
      ENDIF
!     WRITE(6,*) NUM_CRAY_WORDS,LENBUF
!     write(6,*) (ilabel(j),j=1,45)
!     write(6,*) (rlabel(j),j=1,19)
!L
!L 3. Put coded data into BUFOUT for output
!L
      DO I=1,NUM_CRAY_WORDS
        BUFOUT(I)=WORK_ARRAY(I)
      ENDDO
      DO I=NUM_CRAY_WORDS+1,LENBUF
        BUFOUT(I)=0.0
      ENDDO
!L
!L 4. Update lookup for this field
!L
      DO J=1,45
        LOOKUP(J,IENT)=ILABEL(J)
      ENDDO
      DO J=1,19
        RLOOKUP(J+45,IENT)=RLABEL(J)
      ENDDO
      LOOKUP(LBLREC,IENT)=NUM_CRAY_WORDS
      LOOKUP(LBEGIN,IENT)=IWA
      LOOKUP(LBNREC,IENT)=NUM_CRAY_WORDS
      LOOKUP(DATA_TYPE,IENT)=1
      LOOKUP(NADDR,IENT)=IWA
!L
!L 5. Output BUFOUT
!L
      CALL SETPOS_single(UNITPP,IWA,ICODE)
      CALL BUFFOUT_single(UNITPP,BUFOUT(1),NUM_CRAY_WORDS,LEN_IO,IX)
      RETURN
      END SUBROUTINE GRIB_FILE
#endif
