#if defined(FLDC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Program: FIELDCOS 
!LL
!LL  Purpose:
!LL  To read a model dump or direct access fieldsfile and convert it to
!LL  a sequential PP file ready for transfer to a different platform.
!LL
!LL   A general note on fieldcos -- When doing a bit compare on
!LL   the output of fieldcos half words may disagree. This is caused
!LL   by the extra half word after an odd number of words in a field
!LL   and is nothing to worry about. (Simon Tett 13/5/92)
!LL
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 
!LL
!LL  Logical components covered: C41
!LL
!LL  Project task: C4
!LL
!LL  External documentation: UM documentation paper Y8
!LL
!LL  ------------------------------------------------------------------
      PROGRAM FIELDCOS
      IMPLICIT NONE
!*---------------------------------------------------------------------
!     Called routines
      EXTERNAL READFF,SETPOS,IOERROR,READ_WRITE
!*---------------------------------------------------------------------
!     arguments for called routines
      CHARACTER                                                         &
     &     CMESSAGE*80                                                  &
                            ! Error message from lower routines
     &    ,INFILE*80                                                    &
                            ! Pathname of input file.
     &    ,FORMAT_OUT*6     ! IBM/IEEE/VAX for output format
                            ! GRIB - pure binary grib stash codes
                            ! GRIB1 - pure binary grib - standard codes
                            ! GRIB2 - pure binary grib - Other table 2
      CHARACTER*(*) RoutineName
      PARAMETER (RoutineName='FIELDCOS')
      LOGICAL                                                           &
     &     UNPACK                                                       &
                                    ! indicates whether to unpack
     &    ,OPER                                                         &
                                    ! indicates whether operational
     &    ,MASS                     ! indicates whether its for MASS

      NAMELIST /PACK/ UNPACK,FORMAT_OUT
      NAMELIST /TYPE/ OPER, MASS
      INTEGER                                                           &
     &     LEN1_LOOKUP                                                  &
                                    ! First dimension of the lookup
     &    ,PP_LEN2_LOOKUP                                               &
                                    ! Size of the LOOKUP on the file
     &    ,PPUNIT                                                       &
                                    ! unit no of required fieldsfile
     &    ,COS_PPUNIT                                                   &
                                    ! unit no of COS output file
     &    ,IEXTRA(10)                                                   &
                                    ! spare for future use
     &    ,ICODE                                                        &
                                    ! return code
     &    ,DATA_ADD                                                     &
                                    ! The word address of the data.
     &    ,IWA                                                          &
                                    ! Word address in call SETPOS
     &    ,LEN_IO                                                       &
                                    ! Length of IO done
     &    ,LEN_FIXHD                ! Length of fixed length header
      PARAMETER(LEN_FIXHD=256)
      INTEGER                                                           &
     &     PP_FIXHD(LEN_FIXHD)      !  Fixed length header
      REAL                                                              &
     &     A_IO                     ! status returned by BUFFIN
      PARAMETER(LEN1_LOOKUP=64)
      DATA UNPACK/.FALSE./
      DATA FORMAT_OUT/'IBM   '/
      DATA OPER/.FALSE./
      DATA MASS/.FALSE./
! defines NUNITS
#include "chsunits.h"
#include "cgribtab.h"
      LOGICAL FLAG                  ! =T/F file exists/not
      COMMON /FLAG_IO/FLAG(NUNITS)  ! needed for BUFFIN check
!*---------------------------------------------------------------------
!    LOCAL VARIABLES
      INTEGER                                                           &
     &     I                                                            &
                                  ! local counter
     &    ,IX                                                           &
                                  ! used as a dummy variable in UNIT
     &    ,ERR                                                          &
     &    ,DIAG_UNIT
      LOGICAL LCAL360   ! 360 day calendar switch
!  Initialise LCAL360
      DATA LCAL360 /.FALSE./
      CHARACTER*80 DIAGFILE
!=====================================================================
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS LOOKUP GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!=====================================================================
      READ(5,PACK)
      READ(5,TYPE)

      IF (FORMAT_OUT /= 'IEEE' .AND. MASS) THEN
        ! DEPENDS ON: ereport
        ICODE    = 1
        CMESSAGE = 'Can only use IEEE format with MASS.'
        CALL EREPORT(RoutineName,ICODE,CMESSAGE)
      END IF

      WRITE(6,*)'  UNPACK  ',UNPACK
      WRITE(6,*)'  FORMAT  ',FORMAT_OUT
      WRITE(6,*)'  OPER    ',OPER
      WRITE(6,*)'  MASS    ',MASS
      DO I=1,10
        IEXTRA(I)=0
      ENDDO

      DIAG_UNIT = 7
      CALL GET_FILE(DIAG_UNIT,DIAGFILE,80,ICODE)
      OPEN(UNIT=DIAG_UNIT,FILE=DIAGFILE)

      PPUNIT=10
      COS_PPUNIT=11
! -------------------------------------------------------------------
! If FORMAT_OUT is GRIB1 or GRIB2 initialise grib field code
! conversion table
      IF (FORMAT_OUT == 'GRIB1') THEN
! DEPENDS ON: grib_table_init1
        CALL GRIB_TABLE_INIT1
      ELSE IF (FORMAT_OUT == 'GRIB2') THEN
! DEPENDS ON: grib_table_init2
        CALL GRIB_TABLE_INIT2
      ENDIF
!L-------------Read in the FIXED length header------------------------
      CALL GET_FILE(PPUNIT,INFILE,80,ICODE)
! DEPENDS ON: file_open
      CALL FILE_OPEN(PPUNIT,INFILE,80,0,1,ERR)
      FLAG(PPUNIT)=.TRUE.          ! needed for BUFFIN check
! DEPENDS ON: buffin
      CALL BUFFIN(PPUNIT,PP_FIXHD,LEN_FIXHD,LEN_IO,A_IO)
      IF(A_IO /= -1.0.OR.LEN_IO /= LEN_FIXHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('Buffer in fixed length header',A_IO,LEN_IO,     &
     &                  LEN_FIXHD)
          CMESSAGE='FIELDCOS : I/O error reading FIXED LENGTH HEADER'
          ICODE=2
          WRITE(6,*)' I/O error reading FIXED LENGTH HEADER'
! DEPENDS ON: abort_io
          CALL ABORT_IO('FIELDCOS',CMESSAGE,ICODE,PPUNIT)
      ENDIF
      DATA_ADD=PP_FIXHD(160)-1 ! Start address for the data.
      IWA= PP_FIXHD(150)-1     ! Start address for the lookup table.
      PP_LEN2_LOOKUP=PP_FIXHD(152)
      WRITE(6,*)' PP_LEN2_LOOKUP  ',PP_LEN2_LOOKUP
      WRITE(6,*)' dump type=',pp_fixhd(5),                              &
     &       ' 3=fieldsfile,1=dump,2=time mean dump,4=ancil,5=bound'
! DEPENDS ON: read_write
      CALL READ_WRITE(PP_LEN2_LOOKUP,LEN1_LOOKUP,DATA_ADD,              &
     &                PP_FIXHD,                                         &
     &                IWA,UNPACK,FORMAT_OUT,PPUNIT,COS_PPUNIT,          &
     &                IEXTRA,OPER,MASS,ICODE,CMESSAGE,LCAL360)
      IF(ICODE /= 0) THEN
! DEPENDS ON: ereport
        CALL EREPORT(RoutineName,ICODE,CMESSAGE)
      ENDIF
      STOP
      END PROGRAM FIELDCOS
#endif

