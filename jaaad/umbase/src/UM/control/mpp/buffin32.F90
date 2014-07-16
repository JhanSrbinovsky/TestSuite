#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(RECON) || defined(UTILIO) || defined(FLDIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Parallel UM version of BUFFIN32
!
! Subroutine Interface:
      SUBROUTINE BUFFIN32(NFT,ARRAY,ISIZE,LEN_IO,IOSTAT)
#if defined(RECON)
      Use Rcf_Parvars_Mod
#endif
      IMPLICIT NONE
!
! Description:
!  This routine provides a BUFFIN32 routine for the Parallel Unified
!  Model. It is used where the same data which is 32-bit packed has
!  to be read in to all processors. (Contrast with READ_MULTI where
!  each processor reads its own local data from a larger global
!  field.)
!
! Method:
!  The C BUFFIN32 is renamed BUFFIN32_SINGLE under *DEF,MPP. This
!  routine causes BUFFIN32_SINGLE to be called by PE 0 only, and
!  then the data is broadcast to all other processors.
!
! Current Code Owner: Paul Dando
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    5.5    25/03/03 Code copied and modified from BUFFIN during
!                    development of portable I/O changes allowing
!                    big-endian data to be input and output on
!                    little-endian platforms.          P.Dando
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes
!    6.0    13/10/03 Declare IOSTAT as REAL rather than INTEGER.
!                                                      P.Dando
!
! Code Description:
!
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Subroutine Arguments with INTENT IN:

      INTEGER, INTENT(IN)  :: NFT          ! FORTRAN unit number
      INTEGER, INTENT(IN)  :: ISIZE        ! no. of words to write out

      REAL,    INTENT(IN)  :: ARRAY(ISIZE) ! Array to write out
      REAL,    INTENT(OUT) :: IOSTAT       ! Return code

! Subroutine Arguments with INTENT OUT:

      INTEGER, INTENT(OUT) :: LEN_IO       ! no. of words written out

! Parameters and Common blocks

#if !defined(RECON)
#include "parvars.h"
#endif

! Local variables

      INTEGER :: info                      ! Broadcast information code
      REAL    :: return_codes(2)           ! Array to hold BUFFO32 code


! ------------------------------------------------------------------

      IOSTAT=-1.0
      LEN_IO=ISIZE


      IF (mype  ==  0) THEN
        CALL BUFFIN32_SINGLE(NFT,ARRAY,ISIZE,LEN_IO,IOSTAT)
        return_codes(1)=LEN_IO
        return_codes(2)=IOSTAT
      ENDIF

!--get the broadcast flag
      call find_unit_bcast_flag(nft, info)
!--skip the broadcasts if the flag is set
      IF (info  ==  0) THEN
        CALL GC_RBCAST(1,ISIZE,0,nproc,info,ARRAY)    ! data
        CALL GC_RBCAST(2,2,0,nproc,info,return_codes) ! return codes
      ENDIF

      LEN_IO=return_codes(1)
      IOSTAT=return_codes(2)

      RETURN
      END SUBROUTINE BUFFIN32

#endif
