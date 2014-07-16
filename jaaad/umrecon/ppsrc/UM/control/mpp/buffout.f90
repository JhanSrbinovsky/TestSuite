

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Parallel UM version of BUFFOUT
!
! Subroutine Interface:
      SUBROUTINE BUFFOUT(NFT,ARRAY,ISIZE,LEN_IO,IOSTAT)

      Use Rcf_Parvars_Mod


      IMPLICIT NONE
!
! Description:
!  This routine provides a BUFFOUT routine for the Parallel Unified
!  Model. It is used where the same data has to be written out from
!  all processors. (Contrasting with WRITE_MULTI where each
!  processor writes its own local data to a larger global field).
!
! Method:
!  The C BUFFOUT is renamed BUFFOUT_SINGLE under *DEF,mpp. This
!  routine causes BUFFOUT_SINGLE to be called by PE 0 only.
!  Note : No check is made that all processors are attempting
!         to write out identical data. It is assumed that the
!         data on PE 0 is the same as that on all other processors.
!
! Current Code Owner: Paul Burton
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    3.5    5/1/95   New DECK created for the Parallel Unified
!                    Model. P.Burton + D.Salmond
!    4.1    18/3/96  Broadcast return code from I/O
!LL  5.1    10/04/00 New reconfiguration support.
!    5.3    22/11/01 Enable mpp as the only option for
!                    small executables         E.Leung
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes
!
! Subroutine Arguments:

      INTEGER                                                           &
     &  NFT                                                             &
                     !   IN : FORTRAN unit number
     & ,LEN_IO                                                          &
                     !  OUT : no. of words written out
     & ,ISIZE        !   IN : no. of words to write out

      REAL                                                              &
     &  IOSTAT                                                          &
                     !  OUT : Return code
     & ,ARRAY(ISIZE) !   IN : Array to write out

! Parameters and Common blocks





! Local variables

      INTEGER info
      REAL stats(2)


! ------------------------------------------------------------------

      IOSTAT=-1.0
      LEN_IO=ISIZE

      IF (mype  ==  0) THEN
        CALL BUFFOUT_SINGLE(NFT,ARRAY,ISIZE,LEN_IO,IOSTAT)
        stats(1)=LEN_IO
        stats(2)=IOSTAT
      ENDIF

      CALL GC_RBCAST(521,2,0,nproc,info,stats)
      LEN_IO=stats(1)
      IOSTAT=stats(2)

      RETURN
      END SUBROUTINE BUFFOUT

