

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE PR_LFLD----------------------------------------
!LL
!LL  Purpose: Prints out selected values from logical data
!LL           using information from associated PP header.
!LL
!LL  Written by D. Robinson
!LL
!LL  Model            Modification history:
!LL version  date
!LL   3.3  22/11/93  New routine (adapted from deck PRIFLD1A)
!LL   5.1  31/03/00  Added comma missing from write format. D.P.Matthews
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!LL
!LL  System component: R30/W30
!LL
!LL  System task: F3
!LL
!LL  Programming standard:
!LL           Unified Model Documentation Paper No 3
!LL           Version No 1 15/1/90
!LL
!LL  Documentation:
!LL           Unified Model Documentation Paper No F3
!LL           Version No 5 9/2/90
!LL
!LL------------------------------------------------------------
!*L Arguments:-------------------------------------------------
      SUBROUTINE PR_LFLD(LOOKUP,RLOOKUP,LEN1_LOOKUP,LD1,K)

      IMPLICIT NONE

      INTEGER                                                           &
     & K                                                                &
                     !IN Field number ie position in 2nd dim
                     !   of LOOKUP
     &,LEN1_LOOKUP                                                      &
                     !IN First dimension of LOOKUP table
     &,LOOKUP(LEN1_LOOKUP,*)  !IN Integer equivalence of PP LOOKUP

      REAL                                                              &
     & RLOOKUP(LEN1_LOOKUP,*) !IN Real equivalence of PP LOOKUP

      LOGICAL                                                           &
     & LD1(*)        !IN Kth field in data array
! -------------------------------------------------------------
!*L External subroutines called:-------------------------------
! None
!--------------------------------------------------------------
!*L Local control constants:-----------------------------------
      INTEGER                                                           &
     & NS_PTS                                                           &
                     !PARAM No of points down to print
     &,EW_PTS        !PARAM No of points across to print
      PARAMETER(NS_PTS=6,EW_PTS=5)
! -------------------------------------------------------------
! Workspace usage:---------------------------------------------
      REAL LON(EW_PTS)     ! Longitudes printed out
      INTEGER I(EW_PTS)    ! Index of values printed out
      CHARACTER*12 DASH(EW_PTS)  !Stores dashed lines
!*-------------------------------------------------------------
! Local variables:---------------------------------------------
      INTEGER                                                           &
     & N_ROWS                                                           &
                   ! No of rows in field
     &,N_COLS                                                           &
                   ! No of colums in field
     &,ROW                                                              &
                   ! Row number
     &,R_INC,F_INC                                                      &
                   ! No of rows/points between printed lines
     &,J,L                                                              &
                   ! Loop counts
     &,EW_PRINT                                                         &
                   ! No of E-W values printed out
     &,POS_MIN                                                          &
                   ! Position of Minimum value of field
     &,POS_MAX                                                          &
                   ! Position of Maximum value of field
     &,F_MIN                                                            &
                   ! Minimum value of field
     &,F_MAX       ! Maximum value of field

      REAL                                                              &
     & LAT         ! Latitude
!--------------------------------------------------------------

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

!L Internal structure: None

! Initialise string used to create table boundaries
      DO J=1,EW_PTS
        DASH(J)='------------'
      ENDDO

      IF(LOOKUP(LBCODE,K) == IMDI) THEN
!       IF LBCODE IS MISSING DATA, ASSUME THAT THE FIELD IN DUMP
!       HAS NOT BEEN WRITTEN TO BY STASH.
!       THIS SHOULD ONLY OCCUR TO DIAGNOSTIC PARTS OF THE DUMP BEFORE
!       FIRST WRITE BY STASH TO THAT AREA/HEADER.
        WRITE(6,*) 'MESSAGE FROM PR_LFLD'
        WRITE(6,*) 'LBCODE NOT SET; ASSUME DATA NOT SET. NO PRINT'
        RETURN
      END IF

! No of rows and columns in field
      N_ROWS=LOOKUP(LBROW,K)
      N_COLS=LOOKUP(LBNPT,K)


      IF(N_COLS /= 0.AND.N_COLS /= IMDI)THEN

! No of E-W values to be printed
      EW_PRINT=MIN(N_COLS,EW_PTS)

! Calculate longitudes and addresses of values to be printed from 1st ro
      I(1)=1
      LON(1)=RLOOKUP(BZX,K)+RLOOKUP(BDX,K)
      DO 100 J=1,EW_PTS-2
      I(J+1)=I(J)+N_COLS/(EW_PTS-1)
      LON(J+1)=LON(J)+RLOOKUP(BDX,K)*(N_COLS/(EW_PTS-1))
100   CONTINUE
      I(EW_PTS)=N_COLS
      LON(EW_PTS)=RLOOKUP(BZX,K)+RLOOKUP(BDX,K)*N_COLS

! Initialise row and field pointers
      ROW=1
      LAT=RLOOKUP(BZY,K)+RLOOKUP(BDY,K)
      R_INC=N_ROWS/(NS_PTS-1)
      F_INC=R_INC*N_COLS

! Print 1st row
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)
      WRITE(6,'('' FIELD NO'',I4,'':'',9(F10.3,2X))')                   &
     &K,(LON(J),J=1,EW_PRINT)
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)

! Print remaining rows except last
      DO 200 L=1,NS_PTS-1
      WRITE(6,'(1X,I3,'':'',F8.3,'':'',3X,9(L9,3X))')ROW,LAT,           &
     &(LD1(I(J)),J=1,EW_PRINT)
      DO 300 J=1,EW_PTS
      I(J)=I(J)+F_INC
300   CONTINUE
      ROW=ROW+R_INC
      LAT=LAT+R_INC*RLOOKUP(BDY,K)
200   CONTINUE

! Calculate addresses used to print values for last row
      I(1)=1+(N_ROWS-1)*N_COLS
      DO 400 J=1,EW_PTS-2
      I(J+1)=I(J)+N_COLS/(EW_PTS-1)
400   CONTINUE
      I(EW_PTS)=N_ROWS*N_COLS

! Set row pointers to last row
      LAT=RLOOKUP(BZY,K)+RLOOKUP(BDY,K)*N_ROWS
      ROW=N_ROWS

! Print last row
      WRITE(6,'(1X,I3,'':'',F8.3,'':'',3X,9(L9,3X))')ROW,LAT,           &
     &(LD1(I(J)),J=1,EW_PRINT)
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)
      ELSE

! Print out summary of non standard fields

      EW_PRINT=MIN(EW_PTS,LOOKUP(LBLREC,K))
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)
      WRITE(6,'('' FIELD NO'',I4,'':  DATA NOT ON MODEL GRID''          &
     &,'' SO FIRST FEW VALUES PRINTED'')')K
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)
      WRITE(6,'(1X,3X,'':'',8X,'':'',3X,9(L9,3X))')                     &
     &(LD1(J),J=1,EW_PRINT)
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)

      ENDIF

      WRITE(6,'('' '')')

      RETURN
      END SUBROUTINE PR_LFLD

