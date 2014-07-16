
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Parallel RCF interface to BUFFOUT

Module Rcf_Write_Multi_Mod

!  Subroutine Rcf_Write_Multi - parallel interface to buffout
!
! Description:
!  This routine provides an interface to BUFFOUT for the parallel
!  Reconfiguration. It is used where each process must write out a
!  local section of a global field.
!
! Method:
!  Each processor sends its local part of the global field to PE 0
!  which assembles all the parts, and then writes them to disk.
!  Fields which are compressed to land points are expanded before
!  sending to PE 0, PE 0 then compresses the global field before
!  writing it to disk.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   20/06/01   Fix to size of compbuf.  P.Selwood
!   5.4    9/09/02   Grib Prep. Add Multi_PE flag to confirm 'gather'
!                                                          Rod.Sharp
!   5.5   28/02/03   Portability changes allowing for big_endian
!                    I/O on little_endian platforms.        P.Dando
!   6.1   28/08/04   Allow buffering of IO on big endian systems.
!                    JC Rioual, P Selwood
!   6.2   21/10/05   Replace GSYNC with SSYNC. P.Selwood.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Write_Multi(  NFT, D1, ISIZE, LEN_IO, LOCAL_LEN, &
                             IOSTAT, LOOKUP, FIXHD12,           &
                             STASH_RECORD, Multi_PE )

Use Rcf_Parvars_Mod, Only : &
    mype,                   &
    nproc

Use Rcf_General_Gather_Field_Mod, Only : &
     &    Rcf_General_Gather_Field

Use Rcf_Ppx_Info_Mod, Only : &
    STM_Record_Type

IMPLICIT NONE

! Subroutine Arguments:
Integer, Intent(In)    :: Nft         ! Fortran unit number
Integer, Intent(In)    :: Isize       ! no. of words to write out
Integer, Intent(In)    :: Fixhd12     ! 12th element of fixed hdr.
Integer, Intent(In)    :: Lookup(64)  ! Lookup table
Integer, Intent(Out)   :: Len_IO      ! no. or words written out
Integer, Intent(Out)   :: Local_Len   ! size of local field written out

Real, Intent(In)       :: D1(*)       ! Array to write out
Real, Intent(Out)      :: IOstat      ! return code

Type( STM_Record_Type ), Intent(In) :: STASH_RECORD

Logical, Intent(In)    :: Multi_PE    ! flag to confirm 'gather' req'd

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

! Local variables
Integer                :: i                 ! looper
Integer                :: info              ! Gcom code
Real                   :: CompBuf(lookup(lbnrec)) ! compression buffer
Real                   :: buf( isize * 2 )  ! buffer for data
!DIR$ CACHE_ALIGN buf

! ------------------------------------------------------------------
IOSTAT=-1.0
LEN_IO=ISIZE
LOCAL_LEN=0

! If data spread across multiple PE's then gather it.
If (Multi_PE) Then
  ! Gather the field from the local D1 array to buf
  CALL Rcf_General_Gather_Field( D1,        buf,            &
                                 LOCAL_LEN, LOOKUP(LBLREC), &
                                 STASH_RECORD,   0 )
Else
! data is already held as a full copy on each PE.
  buf(: LOOKUP(LBLREC)) = D1(: LOOKUP(LBLREC))
End If

IF (mype .EQ. 0) THEN
!       Does this field need to be compressed?
  IF(MOD((LOOKUP(LBPACK)),10) .EQ. 2) THEN
    IF(LOOKUP(DATA_TYPE) .EQ. 1) THEN
! DEPENDS ON: pack21
      CALL PACK21( LOOKUP(LBLREC), buf, COMPBUF )
    ENDIF
  ELSE ! no compression required - just do a copy
    DO i=1, LOOKUP( LBLREC )
      COMPBUF(i)=buf(i)
    ENDDO
  ENDIF

! Now write out the global field

  IF(MOD((LOOKUP(LBPACK)),10) .EQ. 2) THEN
    IF(LOOKUP(DATA_TYPE) .EQ. 1) THEN
! Data is packed using CRAY 32 bit method - note that we need to write
! out 2*ISIZE 32 bit words using BUFFO32
      CALL buffo32_single( NFT,COMPBUF,2*ISIZE,LEN_IO,IOSTAT )
! And then halve LEN_IO to satisfy tests against ISIZE
      LEN_IO = LEN_IO/2
    ENDIF
  ELSE
! For non-packed data
    CALL buffout_single( NFT, COMPBUF, ISIZE, LEN_IO, IOSTAT )
  ENDIF



ENDIF ! am I PE 0 ?

Call GC_SSYNC( nproc, info )

! If the field was compressed for writing on disk, we need to compress
! and expand the field in memory. This ensures the same field exists in
! memory that would exist if this dump was read back in.

If (MOD((LOOKUP(LBPACK)),10) .EQ. 2) Then
  If (LOOKUP(DATA_TYPE) .EQ. 1) Then
! DEPENDS ON: pack21
    Call PACK21(   LOCAL_LEN, D1, COMPBUF )
! DEPENDS ON: expand21
    Call EXPAND21( LOCAL_LEN, COMPBUF, D1 )
  End If
End If

Return
End Subroutine Rcf_Write_Multi
End Module Rcf_Write_Multi_Mod
