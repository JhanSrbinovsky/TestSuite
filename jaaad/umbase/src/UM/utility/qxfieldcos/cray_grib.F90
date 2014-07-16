#if defined(FLDC) && defined(CRAY)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    subroutine: CRAY_GRIB
!
!    Purpose: To read a   direct access PP file  and convert it to a
!    pure grib file ready to be passed to HDS or workstation.
!
!    Model            Modification history from model version 3.3:
!   version  Date
!    4.0    31/03/95  : Added to FIELDCOS
!
!    Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
!    Logical components covered: C41
!
!    Project task: C4
!
!    External documentation:
!
!
      SUBROUTINE CRAY_GRIB(IDIM,PPUNIT,TOTAL_WORDS,FORMAT_OUT,          &
     &             LEN1_LOOKUP,PP_LEN2_LOOKUP,PP_FIXHD,LOOKUP,          &
     &             ROOKUP,ENTRY_NO,DATA_ADD,MODEL_FLAG,                 &
     &             COS_PPUNIT,IEXTRA,ICODE,CMESSAGE)

      IMPLICIT NONE
!     Arguments
      CHARACTER                                                         &
     &     CMESSAGE*(*)                                                 &
                                  !OUT error messages
     &    ,FORMAT_OUT*6           !IN format required
      LOGICAL                                                           &
     &     MODEL_FLAG             !IN True => dumps, False => fieldsfile
      INTEGER                                                           &
     &     PPUNIT                                                       &
                                  !IN unit no of required fieldsfile
     &    ,TOTAL_WORDS                                                  &
                                  !IN total number of words written
     &    ,COS_PPUNIT                                                   &
                                  !IN unit no of COS output file
     &    ,IDIM                                                         &
                                  !IN NUM_VALUES rounded to an even no
!                                 !  used to dimension The output array
     &    ,DATA_ADD                                                     &
                                  !IN The word address of the data.
     &    ,LEN1_LOOKUP                                                  &
                                  !IN First dimension of the lookup
     &    ,PP_LEN2_LOOKUP                                               &
                                  !IN Size of the LOOKUP on the file
     &    ,IEXTRA(10)                                                   &
                                  !IN Used within READFF
     &    ,ENTRY_NO                                                     &
                                  !IN Lookup entry no of the Field.
     &    ,PP_FIXHD(*)                                                  &
                                  !IN PPfile fixed header
     &    ,LOOKUP(LEN1_LOOKUP,PP_LEN2_LOOKUP)                           &
                                               !IN integer lookup
     &    ,ICODE                  !OUT error code
      REAL                                                              &
     &     ROOKUP(LEN1_LOOKUP,PP_LEN2_LOOKUP)    !IN Real lookup
!----------------------------------------------------------------------
!     Called routines
      EXTERNAL READFF,GBYTES,SBYTES
!----------------------------------------------------------------------
!    LOCAL VARIABLES
      INTEGER                                                           &
     &     MAX_LEN_ILABEL                                               &
                             ! maximum length of INT part of pp header
     &    ,MAX_LEN_RLABEL    ! maximum length of REAL part of pp header
      PARAMETER (MAX_LEN_ILABEL=45,MAX_LEN_RLABEL=32)
      INTEGER                                                           &
     &     LEN_ILABEL                                                   &
                           ! number of values in ILABEL
     &    ,LEN_RLABEL                                                   &
                           ! number of values in RLABEL
     &    ,ILABEL(MAX_LEN_ILABEL)                                       &
                                        ! holds integer part of LOOKUP
     &    ,I                                                            &
                          ! local counter
     &    ,CARRY                                                        &
                          ! local counter
     &    ,NEW_CODE                                                     &
                          ! new field code
     &    ,ITEM                                                         &
                           ! item code
     &    ,SECTION                                                      &
                           ! section code
     &    ,SECTION1(16)    ! UM octet 9 value

      REAL                                                              &
     &     FIELD(IDIM)                                                  &
                                  ! array holding data
     &    ,RLABEL(MAX_LEN_RLABEL) ! holds real part of LOOKUP

#include "clookadd.h"
#include "cgribtab.h"
!----------------------------------------------------------------------

!  access the Fields File.
! DEPENDS ON: readff
      CALL READFF(PPUNIT,FIELD,IDIM,ENTRY_NO,ILABEL,RLABEL,IEXTRA,      &
     &            PP_LEN2_LOOKUP,LEN1_LOOKUP,PP_FIXHD,LOOKUP,ROOKUP,    &
     &            DATA_ADD,MODEL_FLAG,MAX_LEN_ILABEL,MAX_LEN_RLABEL,    &
     &            LEN_ILABEL,LEN_RLABEL,ICODE,CMESSAGE)

      IF(ICODE /= 0) RETURN

!-----------------------------------------------------------------
!  Alter field codes if required
!  FORMAT_OUT
!  GRIB  - UM stash codes - no change
!  GRIB1 - attempt to alter codes to standard table 2 values
!  GRIB2 - attempt to alter codes to other user table 2
!
      IF (FORMAT_OUT == 'GRIB1'.OR.FORMAT_OUT == 'GRIB2') THEN
        SECTION=ILABEL(42)/1000
        ITEM=ILABEL(42) - SECTION*1000
        NEW_CODE=GRIB_TABLE(SECTION,ITEM)
        IF (NEW_CODE == -99) THEN
      WRITE(6,*)' No standard grib code for field ',ilabel(42),' field  &
     & will not be output'
          RETURN
        ELSE
!  Assumes running on a 64 bit word machine
!  Therefore section 0 is field (1) & section 1 starts at field(2)
! Need to alter octets 4 and 9 in section 1
! decode first 16 octets of grib message
          CALL GBYTES(field(2),section1(1),0,8,0,16)
          SECTION1(4)=1
          SECTION1(9)=new_code
! recode first 16 octets of grib message
          CALL SBYTES(field(2),section1(1),0,8,0,16)
        ENDIF
      ENDIF
!-----------------------------------------------------------------
! write out pure grib code

        WRITE(COS_PPUNIT) (FIELD(I),I=1,ILABEL(LBLREC))
        WRITE(6,100) ILABEL(42),ILABEL(LBLREC)
        TOTAL_WORDS=TOTAL_WORDS+ilabel(lblrec)
 100  FORMAT(1x,' written out grib for ',i6,' length of data',i8)

      RETURN
      END SUBROUTINE CRAY_GRIB
! =====================================================================
#endif



