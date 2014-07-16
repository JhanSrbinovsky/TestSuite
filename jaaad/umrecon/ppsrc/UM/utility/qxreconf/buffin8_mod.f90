
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

Module BuffIn8_Mod

!+ Parallel UM version of BUFFIN8
!
contains
! Subroutine Interface:
Subroutine BuffIn8(NFT,Array,ArryLen,CharLen,Len_IO,IOStat)

!
! Description:
!  This routine provides a BUFFIN8 routine for the Parallel Unified
!  Model. It is used where the same data has to be read in to
!  all processors. (Contrasting with READ_MULTI where each
!  processor reads its own local data from a larger global field).
!
! Method:
!  The C BUFFIN8 is renamed BUFFIN8_SINGLE under #if defined(mpp).
!  This routine causes BUFFIN8_SINGLE to be called by PE 0 only, and
!  then the data is broadcast to all other processors.
!
! *******************************************************************
! ***                           WARNING                           ***
! *******************************************************************
!  At present this routine is intended to be used with one of
!  ArryLen or CharLen set to 1.
!  It should handle a single multiple charcter string fine or
!  Arrays of single characters but not neccessarily arrays of
!  multiple character strings
! *******************************************************************
!
! History:
!  Model    Date     Comment and Author :
!  version
!    5.4    22/05/02 Code copied and modified from BUFFIN during
!                    development of GRIB reading code - R.Sharp
!
Use Rcf_Parvars_Mod

IMPLICIT NONE

! Subroutine Arguments Intent IN:

Integer, Intent(IN)  :: NFT          !  FORTRAN unit number
Integer, Intent(IN)  :: ArryLen      !  No. of array elements
Integer, Intent(IN)  :: CharLen      !  Length of each element

! Subroutine Arguments Intent OUT:

Integer  , Intent(OUT) :: Len_IO       !  no. of words read in
Integer  , Intent(OUT) :: IOStat       !  Return code
Character(Len=CharLen), Intent(OUT) :: Array(ArryLen)
                                       !  Array to read data in to

! Local variables

Integer              :: Info
Integer              :: ISize        !  no. of words to be read in

!-------------------------------------------------------------------
!calculate length in bytes of buffer to be read
ISize = Charlen * ArryLen   ! assuming 1 char = 1 byte
IOStat  =  -1
Len_IO  =   ISize


IF (mype == 0) THEN
  Call BUFFIN8_SINGLE(NFT,Array,ISize,Len_IO,IOStat)
ENDIF

!--get the broadcast flag
Call find_unit_bcast_flag(NFT, Info)

!--skip the broadcasts if the flag is set
If ( Info == 0 ) Then
  Call GC_CBCAST(1,ISize,0,nproc,Info,Array)    ! data
  Call GC_IBCAST(2,1,0,nproc,Info,IOStat)       ! IO error flag
  Call GC_IBCAST(3,1,0,nproc,Info,Len_IO)       ! Length read
End If


Return
End Subroutine Buffin8

End Module BuffIn8_Mod
