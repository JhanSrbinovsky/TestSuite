

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Parallel UM version of BUFFO32
!
! Subroutine Interface:
      SUBROUTINE BUFFO32(NFT,ARRAY,ISIZE,LEN_IO,IOSTAT)

      Use Rcf_Parvars_Mod


      IMPLICIT NONE
!
! Description:
!  This routine provides a BUFFO32 routine for the Parallel Unified
!  Model. It is used where the same data which is 32-bit packed has
!  to be written out from all processors. (Contrast with WRITE_MULTI
!  where each processor reads its own local data from a larger global
!  field.)
!
! Method:
!  The C BUFFO32 is renamed BUFFO32_SINGLE under *DEF,mpp. This
!  routine causes BUFFO32_SINGLE to be called by PE 0 only.
!  Note : No check is made that all processors are attempting
!         to write out identical data. It is assumed that the
!         data on PE 0 is the same as that on all other processors.
!
! Current Code Owner: Paul Dando
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    5.5    25/03/03 Code copied and modified from BUFFOUT as
!                    part of the portable I/O changes allowing
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





! Local variables

      INTEGER :: info                      ! Broadcast information code
      REAL    :: stats(2)                  ! Array to hold BUFFO32 code

! ------------------------------------------------------------------

      IOSTAT=-1.0
      LEN_IO=ISIZE

      IF (mype  ==  0) THEN
        CALL BUFFO32_SINGLE(NFT,ARRAY,ISIZE,LEN_IO,IOSTAT)
        stats(1)=LEN_IO
        stats(2)=IOSTAT
      ENDIF

      CALL GC_RBCAST(521,2,0,nproc,info,stats)
      LEN_IO=stats(1)
      IOSTAT=stats(2)

      RETURN
      END SUBROUTINE BUFFO32

