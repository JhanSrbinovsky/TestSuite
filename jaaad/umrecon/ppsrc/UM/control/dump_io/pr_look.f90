

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE PR_LOOK----------------------------------------
!LL
!LL  Purpose: Prints out Kth 64-word PP header
!LL
!LL AD, DR      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL   3.1  05/02/93    Trap use of user defined PPXREF file
!LL                    Author: A. Dickinson    Reviewer: R. Stratton
!     3.5  27/06/95  Submodels project. Replace call to RDPPXRF by
!                    call to GETPPX
!                    Author D.M.Goddard    Reviewer S Swarbrick
!LL   4.0  12/09/95    Change NPERIODS to LBUSER, LBUSER to LBPLEV
!LL                    as changes in CLOOKADD and PPHEAD1A.
!LL                    (Andrew Brady)
!LL  4.0  12/10/95  Chg. FORMAT, as last LBUSER is now MODEL_CODE. RTHB
!   4.0  24/08/95    Add dummy argument in call to GETPPX_REC
!                    Author D.M.Goddard    Reviewer S Swarbrick
!     4.1  18/06/96    Changes to cope with changes in STASH addressing
!                      Author D.M. Goddard.
!     5.1  10/04/00  New reconfiguration support. P.Selwood.
!     5.1  02/05/00  Improve layout of output. D. Robinson.
!     5.3  11/03/02  Change format for LBSRCE. D. Robinson
!     6.0  30/12/03  Set Model Code from 10 to 1 for FieldCalc
!                    diagnostics to use Atmos StashMaster file.
!                    D. Robinson
!     6.2  18/08/05  Allow RCF_EXPPX not to find STASHmaster entry.
!                                              R.Sharp
!
!LL  System component: R30/W30
!LL
!LL  System task: F3
!LL
!LL  Programming standard:  Unified Model Documentation Paper No 3
!LL                         Version No 1 15/1/90
!LL
!LL  Documentation:  Unified Model Documentation Paper No F3
!LL                  Version No 5 9/2/90
!LL
!LL------------------------------------------------------------
!*L Arguments:-------------------------------------------------
      SUBROUTINE PR_LOOK(                                               &
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -recon def line to allow for other small
!                execs which had used the recon def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     &  PPXI,PPXC,ppxRecs,                                              &
! End of comdeck
     &                   LOOKUP,RLOOKUP,LEN1_LOOKUP,K)
      Use Rcf_Exppx_Mod, Only : Rcf_Exppx

      Use Rcf_Ppx_Info_Mod, Only : STM_Record_type

      IMPLICIT NONE

      INTEGER LEN1_LOOKUP     ! IN First dimension of Look Up Table
      INTEGER K               ! IN Field number in Look Up Table
      INTEGER                                                           &
     & LOOKUP(LEN1_LOOKUP,*)  ! IN Integer equivalence of PP LOOKUP
      REAL                                                              &
     & RLOOKUP(LEN1_LOOKUP,*) ! IN Real equivalence of PP LOOKUP
! Dummy variables for the argppx.h argument comdeck
      INTEGER PPXI
      INTEGER PPXRECS
      CHARACTER PPXC

      CHARACTER*36 EXPPXC

! -------------------------------------------------------------
! Workspace usage:---------------------------------------------
! None
! -------------------------------------------------------------
!*L External subroutines called:-------------------------------
!*-------------------------------------------------------------
! Comdecks: ------------------------------------------------------------
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

! Local variables:---------------------------------------------
      INTEGER ICODE             !Error code
      INTEGER ITEM              !STASH item number
      INTEGER SECTION           !STASH section number
      INTEGER MODEL             !Internal model number
      INTEGER I                 !Index

      CHARACTER*36 PHRASE       !Character part of PPXREF record
      CHARACTER*80 CMESSAGE     !Error message

      Type (STM_RECORD_TYPE), POINTER ::  RECORD

!--------------------------------------------------------------------


!L Write time and field type
        ITEM=MOD(LOOKUP(42,K),1000)
        SECTION=(LOOKUP(42,K)-ITEM)/1000
        MODEL=LOOKUP(45,K)

!       All diagnostics under model code of 10 are in section 20
!       of Atmos StashMaster file.
        IF ( MODEL == 10 ) Then
          MODEL = 1
        ENDIF
      ICODE = 0
      RECORD => RCF_EXPPX(MODEL,SECTION,ITEM,.TRUE.)
      IF ( ASSOCIATED ( RECORD ) ) THEN
        PHRASE = RECORD % NAME
      ELSE
        PHRASE = "Name not known"
      END IF
        IF(ICODE /= 0)THEN
         PHRASE='NON-STANDARD FIELD'
        ENDIF

      WRITE (6,100) K,PHRASE

  100 FORMAT(' FIELD NO.', I5,4X,A)

      WRITE (6,101) LOOKUP(LBHR,K),LOOKUP(LBMIN,K),LOOKUP(LBDAT,K),     &
     &  LOOKUP(LBMON,K),LOOKUP(LBYR,K),LOOKUP(LBDAY,K),LOOKUP(LBHRD,K), &
     &  LOOKUP(LBMIND,K),LOOKUP(LBDATD,K),LOOKUP(LBMOND,K),             &
     &  LOOKUP(LBYRD,K),LOOKUP(LBDAYD,K)

  101 FORMAT(                                                           &
     &  ' VALID AT: ',  2I2.2,'Z  ',2(I2.2,'/'),I4.4,' DAY',I6,         &
     &  ' DATA TIME: ', 2I2.2,'Z  ',2(I2.2,'/'),I4.4,' DAY',I6)

!L Rest of header
      WRITE(6,200) (LOOKUP(I,K),I=13,45),(RLOOKUP(I,K),I=46,64)

  200 FORMAT(                                                           &
     &  '   LBTIM   LBFT    LBLREC LBCODE  LBHEM  LBROW  LBNPT',        &
     &  '  LBEXT LBPACK',/,                                             &
     &  1X,2I7,I10,6I7,/,                                               &
     &  '   LBREL   LBFC  LBCFC LBPROC   LBVC  LBRVC  LBEXP',           &
     &  '   LBBEGIN    LBNREC',/,                                       &
     &  1X,7I7,2I10,/,                                                  &
     &  '  LBPROJ  LBTYP  LBLEV LBRSVD LBRSVD LBRSVD LBRSVD   LBSRCE',/ &
     &  ,1X,7I7,I9,/,                                                   &
     &  '  DATA_TYPE     NADDR    LBUSER ITEM_CODE    LBPLEV',          &
     &  '    LBUSER MODEL_CODE',/                                       &
     &  1X,6I10,I11,/,                                                  &
     &  9X,'BULEV',7X,'BHULEV',5X,'BRSVD(3)',5X,'BRSVD(4)',             &
     &  7X,'BDATUM',/,                                                  &
     &  1X, 1P, 5E13.4,/,                                               &
     &  10X,'BACC',9X,'BLEV',8X,'BRLEV',8X,'BHLEV',7X,'BHRLEV',/,       &
     &  1X, 1P, 5E13.4, /,                                              &
     &  9X,'BPLAT',8X,'BPLON',9X,'BGOR',10X,'BZY',10X,'BDY',/,          &
     &  1X, 1P, 5E13.4, /,                                              &
     &  11X,'BZX',10X,'BDX',9X,'BMDI',9X,'BMKS',/,                      &
     &  1X, 1P, 4E13.4)

       WRITE(6,'('' '')')
      RETURN
      END SUBROUTINE PR_LOOK

