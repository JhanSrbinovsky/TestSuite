#if defined(RECON)
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
#include "clookadd.h"
#include "c_mdi.h"

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
#if defined(LITTLE_END)
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
#else
  CALL BUFFIN_SINGLE(NFT,buf,ISIZE,LEN_IO,IOSTAT)
#endif

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


#endif
