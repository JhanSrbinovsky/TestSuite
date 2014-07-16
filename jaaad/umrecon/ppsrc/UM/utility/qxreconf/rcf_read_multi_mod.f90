
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Parallel RCF interface to BUFFIN

Module Rcf_Read_Multi_Mod

!  Subroutine Rcf_Read_Multi - parallel interface to buffin
!
! Description:
!  This routine provides an interface to BUFFIN for the parallel
!  reconfiguration. It is used where each process must read in a
!  local section of a global field.
!
! Method:
!  PE 0 reads in the global field, and then distributes the
!  relevant parts of it to each processor.
!  Fields compressed to land points are expanded by PE 0, and
!  recompressed after being received by the relevant processor.
!
! Derived from UM 4.5 code
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.2   15/11/00   Allows packed fieldsfiles for VAR. P.Selwood.
!   5.3   08/10/01   Ensure polar sea-points are averaged. P.Selwood
!   5.5   05/02/03   Portability changes allowing for big_endian
!                    I/O on little_endian platforms.        P.Dando
!   6.1   04/11/04   Allow buffering of IO on big endian systems.
!                    JC Rioual, P Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Contains

SUBROUTINE Rcf_Read_Multi(  NFT, D1, ISIZE, EXPAND_SIZE, LEN_IO,  &
                            LOCAL_LEN, IOSTAT, LOOKUP, FIXHD12,   &
                            Stash_Record )

Use Rcf_Parvars_Mod, Only : &
    mype,                   &
    nproc,                  &
    glsize

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,                &
    PrStatus_Normal

Use Rcf_Average_Polar_Mod, Only : &
    Rcf_Average_Polar

Use Rcf_General_Scatter_Field_Mod, Only : &
    Rcf_General_Scatter_Field

Use Rcf_Ppx_Info_Mod, Only : &
    STM_Record_Type

IMPLICIT NONE
! Subroutine Arguments:
Integer, Intent(In)   :: nft        ! Fortain unit number
Integer, Intent(In)   :: isize      ! no. of words to read (glob)
Integer, Intent(In)   :: expand_size! size of expanded field
Integer, Intent(In)   :: Fixhd12    ! 12th element of fixed header
Integer, Intent(In)   :: Lookup(64) ! lookup table
Integer, Intent(Out)  :: Len_IO     ! no. of words read in (glob)
Integer, Intent(Out)  :: Local_Len  ! no. of local field words

Real,    Intent(Out)  :: IOstat     ! Return Code
Real,    Intent(Out)  :: D1(*)      ! Array to read data into

Type( STM_Record_Type ), Intent(In) :: Stash_Record

! Comdecks
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

! Local variables
Integer      :: info            ! GCOM return code
Integer      :: io_ret_codes(2) ! return codes from I/O
Integer      :: dimx            ! x dimension returned from COEX
Integer      :: dimy            ! y dimension returned from COEX
Integer      :: idum            ! dummy for COEX
Integer      :: ErrorStatus     ! Error code returned from COEX
Logical      :: averaged        ! Have polar rows been averaged?
Real         :: buf( expand_size ) ! Buffer for reading data into
Real         :: buf2( expand_size) ! Temp. buffer for unpacking.
Character (Len=80) :: Cmessage

Integer, Parameter           :: len_full_word = 64 ! Output word size
                                                   ! for COEX
Character (Len=*), Parameter :: RoutineName = 'Rcf_Read_Multi'
!DIR$ CACHE_ALIGN buf

! ------------------------------------------------------------------

IOSTAT=-1.0
LEN_IO=ISIZE
LOCAL_LEN=0
ErrorStatus = 0

! First thing to do is to read the field in to PE 0
IF (mype == 0) THEN
! Check for 32-bit packed data and call the correct buffin version
! accordingly.  This should also handle 32-bit packed lateral boundary
! data.
  IF (MOD((LOOKUP(LBPACK)),10) .EQ. 2) THEN
! Data is packed using CRAY 32 bit method
    IF (LOOKUP(LBHEM) .EQ. 99) THEN
! For packed lateral boundary data ISIZE is the number of 32-bit words
! to be read
      LEN_IO = ISIZE
      CALL BUFFIN32_SINGLE(NFT,buf,ISIZE,LEN_IO,IOSTAT)
    ELSE
! For other packed fields we need to read in 2*ISIZE 32 bit words
! using BUFFIN32
      CALL BUFFIN32_SINGLE(NFT,buf,2*ISIZE,LEN_IO,IOSTAT)
! And then halve len_io to satisfy tests against ISIZE
      LEN_IO = LEN_IO/2
    ENDIF
  ELSE
! For non-packed data
    CALL BUFFIN_SINGLE(NFT,buf,ISIZE,LEN_IO,IOSTAT)
  ENDIF

!       Has the data been read in OK?
!  IF (.NOT.((IOSTAT .NE. -1.0) .OR. (LEN_IO .NE. ISIZE))) THEN

! We must check to see if it is a 32 bit field on disk, and if
! so, expand it before doing anything with it.
    IF (MOD((LOOKUP(LBPACK)),10) .EQ. 2) THEN
      IF (LOOKUP(DATA_TYPE) .EQ. 1) THEN
! For special case of lateral boundary data, the length
! is given by ISIZE.
        IF (LOOKUP(LBHEM) .EQ. 99) THEN
! DEPENDS ON: expand32b
          CALL EXPAND32B( ISIZE , buf, FIXHD12 )
        ELSE
! DEPENDS ON: expand32b
          CALL EXPAND32B( LOOKUP(LBLREC) , buf, FIXHD12 )
        ENDIF
      ELSE
        IOSTAT=100
      ENDIF
! What about WGDOS packing?
    Else If (MOD((LOOKUP(LBPACK)),10) .EQ. 1) THEN
      ! temporary copy
      buf2(:) = buf(:)
! DEPENDS ON: coex
      Call COEX( buf, expand_size, buf2, isize, dimx, dimy,  &
                 idum, idum, .FALSE., rmdi, len_full_word,   &
                 ErrorStatus, cmessage )

      If ( ErrorStatus /= 0 ) Then
        Call Ereport( RoutineName, ErrorStatus, Cmessage )
    ENDIF
  ENDIF

! Lets do a polar row check/averaging (on p  grids) for read data only
  If (LOOKUP(lbhem) == 0 .AND. LOOKUP(lbegin) == 1 .AND.  &
      Stash_Record % grid_type <= 3 ) Then
    Call Rcf_Average_Polar( buf, glsize(2), glsize(1), .TRUE., averaged)

!  Print a warning if we've done averaging
    If (averaged .AND. PrintStatus >= PrStatus_Normal) Then
      Write (6,*) 'Field had had polar rows averaged in Rcf_Read_Multi'
    End If
  End If

  io_ret_codes(1)=LEN_IO
  io_ret_codes(2)=IOSTAT
ENDIF  ! IF (mype .EQ. 0)

! Broadcast the error codes, so we can check if anything's wrong
CALL GC_IBCAST(333,2,0,nproc,info,io_ret_codes)
LEN_IO=io_ret_codes(1)
IOSTAT=io_ret_codes(2)


If ( (IOSTAT == -1.0) .AND. (Len_IO == Isize) ) Then  ! no errors
  ! Now we can distribute it out to the other processes
  CALL Rcf_General_Scatter_Field( D1,             buf,            &
                                  LOCAL_LEN,      LOOKUP(LBLREC), &
                                  Stash_Record,   0 )
End If


RETURN
End Subroutine Rcf_Read_Multi
End Module Rcf_Read_Multi_Mod


