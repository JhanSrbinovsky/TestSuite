#if defined(RECON)
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
#include "clookadd.h"

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

#if defined(LITTLE_END)
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
#else
    CALL buffout_single( NFT, COMPBUF, ISIZE, LEN_IO, IOSTAT )
#endif



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
#endif
