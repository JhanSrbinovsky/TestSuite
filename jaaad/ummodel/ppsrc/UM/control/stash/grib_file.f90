
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
!LL  4.3   06/02/97  Modify I/O calls for mpp use  P.Burton
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
!*L------------------ COMDECK LOOKADD ----------------------------------
!LL
!LL Purpose : Contains information about the format
!LL           of the PP header
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   4.0  12/09/95   Change NPERIODS to LBUSER3, BRSVD1 to BULEV,
!LL                   BRSVD2 to BHULEV and definitions for BRLEV and
!LL                   BHRLEV. Corresponding changes made to STWORK1A
!LL                   and PPHEAD1A. (Andrew Brady)
!LL  4.0  12/10/95  Change item 45 from lbuser7 to model_code. RTHBarnes
!LL  5.1  17/04/00    Fixed/Free format. P.Selwood.
!LL  5.2  25/09/00    Add LBCC_xxxx variables for the compressed
!LL                   LBC LOOKUP array                  P.Burton
!LL
!LL Programming standard :
!LL
!LL Logical components covered : F092
!LL
!LL Project task :
!LL
!LL External documentation:
!LL
!LLEND -----------------------------------------------------------------
!

! Validity time
      Integer, Parameter :: LBYR   =1   ! Year
      Integer, Parameter :: LBMON  =2   ! Month
      Integer, Parameter :: LBDAT  =3   ! Day of month
      Integer, Parameter :: LBHR   =4   ! Hour
      Integer, Parameter :: LBMIN  =5   ! Minute
      Integer, Parameter :: LBDAY  =6   ! Day number

! Data time
      Integer, Parameter :: LBYRD  =7   ! Year
      Integer, Parameter :: LBMOND =8   ! Month
      Integer, Parameter :: LBDATD =9   ! Day of month
      Integer, Parameter :: LBHRD  =10  ! Hour
      Integer, Parameter :: LBMIND =11  ! Minute
      Integer, Parameter :: LBDAYD =12  ! Day number

      Integer, Parameter :: LBTIM  =13  ! Time indicator
      Integer, Parameter :: LBFT   =14  ! Forcast period (hours)
      Integer, Parameter :: LBLREC =15  ! Length of data record
      Integer, Parameter :: LBCODE =16  ! Grid type code
      Integer, Parameter :: LBHEM  =17  ! Hemisphere indicator
      Integer, Parameter :: LBROW  =18  ! Number of rows in grid
      Integer, Parameter :: LBNPT  =19  ! Number of points per row
      Integer, Parameter :: LBEXT  =20  ! Length of extra data
      Integer, Parameter :: LBPACK =21  ! Packing method indicator
      Integer, Parameter :: LBREL  =22  ! Header release number
      Integer, Parameter :: LBFC   =23  ! Field code
      Integer, Parameter :: LBCFC  =24  ! Second field code
      Integer, Parameter :: LBPROC =25  ! Processing code
      Integer, Parameter :: LBVC   =26  ! Vertical coordinate type
      Integer, Parameter :: LBRVC  =27  ! Coordinate type for reference
                                        ! level

      Integer, Parameter :: LBEXP  =28  ! Experiment number
      Integer, Parameter :: LBEGIN =29  ! Start record
      Integer, Parameter :: LBNREC =30  ! No of records-Direct access
                                        ! only
      Integer, Parameter :: LBPROJ =31  ! Met-O-8 projection number
      Integer, Parameter :: LBTYP  =32  ! Met-O-8 field type
      Integer, Parameter :: LBLEV  =33  ! Met-O-8 level code
      Integer, Parameter :: LBRSVD1=34  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD2=35  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD3=36  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD4=37  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBSRCE =38  ! =1111 to indicate following
                                        ! apply to UM
      Integer, Parameter :: DATA_TYPE =39  ! Indicator for real/int or
                                           ! timeseries
      Integer, Parameter :: NADDR  =40  ! Start address in DATA_REAL or
                                        ! DATA_INT
      Integer, Parameter :: LBUSER3=41  ! Free for user-defined function
      Integer, Parameter :: ITEM_CODE =42  !Stash code
      Integer, Parameter :: LBPLEV =43  ! Pseudo-level indicator (if
                                        ! defined)
      Integer, Parameter :: LBUSER6=44  ! Free for user-defined function
      Integer, Parameter :: MODEL_CODE =45 ! internal model identifier

      Integer, Parameter :: BULEV  =46  ! Upper level boundary
      Integer, Parameter :: BHULEV =47  ! Upper level boundary
      Integer, Parameter :: BRSVD3 =48  ! Reserved for future PP-package
                                        ! use
      Integer, Parameter :: BRSVD4 =49  ! Reserved for future PP-package
                                        ! use
      Integer, Parameter :: BDATUM =50  ! Datum value
      Integer, Parameter :: BACC   =51  ! (Packed fields) Packing
                                        ! accuracy
      Integer, Parameter :: BLEV   =52  ! Level
      Integer, Parameter :: BRLEV  =53  ! Lower level boundary
      Integer, Parameter :: BHLEV  =54  ! (Hybrid levels) A-level of
                                        ! value
      Integer, Parameter :: BHRLEV =55  ! Lower level boundary
      Integer, Parameter :: BPLAT  =56  ! Real latitude of 'pseudo'
                                        ! N Pole
      Integer, Parameter :: BPLON  =57  ! Real longitude of 'pseudo'
                                        ! N Pole
      Integer, Parameter :: BGOR   =58  ! Grid orientation
      Integer, Parameter :: BZY    =59  ! Zeroth latitude
      Integer, Parameter :: BDY    =60  ! Latitude interval
      Integer, Parameter :: BZX    =61  ! Zeroth longitude
      Integer, Parameter :: BDX    =62  ! Longitude interval
      Integer, Parameter :: BMDI   =63  ! Missing data indicator
      Integer, Parameter :: BMKS   =64  ! M,K,S scaling factor

      Integer, Parameter :: LBCC_LBYR   = 1  ! Year
      Integer, Parameter :: LBCC_LBMON  = 2  ! Month
      Integer, Parameter :: LBCC_LBDAT  = 3  ! Day of the month
      Integer, Parameter :: LBCC_LBHR   = 4  ! Hour
      Integer, Parameter :: LBCC_LBMIN  = 5  ! Minute
      Integer, Parameter :: LBCC_LBDAY  = 6  ! Day number
      Integer, Parameter :: LBCC_LBEGIN = 7  ! Start record
      Integer, Parameter :: LBCC_NADDR  = 8  ! Start address of DATA
! Mapping of MPP_LOOKUP; analogous to mapping in PP header

      Integer, Parameter :: P_NADDR=1    ! Address on local PE
      Integer, Parameter :: P_LBLREC=2   ! Local length of record

!*----------------------------------------------------------------------
! NADDR IS LOCATION IN PP-HEADER (LOOKUP) FOR START POSN OF VARIABLE
! ITEM_CODE is the location in PP header for a code defined as
!           (section number)*1000+item number
! DATA_TYPE is the location in the PP header defining data as REAL or
!           INTEGER.
! LBNPT is the location defining the number of points per row
!
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

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
