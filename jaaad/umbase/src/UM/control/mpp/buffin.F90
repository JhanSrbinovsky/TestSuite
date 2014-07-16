#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(RECON) || defined(UTILIO) || defined(FLDIO) \
  || defined(FLUXPROC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Parallel UM version of BUFFIN
!
! Subroutine Interface:
      SUBROUTINE BUFFIN(NFT,ARRAY,ISIZE,LEN_IO,IOSTAT)
#if defined(RECON)
      Use Rcf_Parvars_Mod
#endif
      IMPLICIT NONE
!
! Description:
!  This routine provides a BUFFIN routine for the Parallel Unified
!  Model. It is used where the same data has to be read in to
!  all processors. (Contrasting with READ_MULTI where each
!  processor reads its own local data from a larger global field).
!
! Method:
!  The C BUFFIN is renamed BUFFIN_SINGLE under *DEF,MPP. This
!  routine causes BUFFIN_SINGLE to be called by PE 0 only, and
!  then the data is broadcast to all other processors.
!
! Current Code Owner: Paul Burton
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    3.5    5/1/95   New DECK created for the Parallel Unified
!                    Model. P.Burton + D.Salmond
!    4.1    18/3/96  Broadcast return code from I/O
!    4.2    18/11/96 Add *CALL AMAXSIZE before IOVARS  P.Burton
!    4.2    18/11/96 Ad
!    4.4    24/06/97 Remove use of memory aligned buf array - read
!                    data straight into ARRAY             P.Burton
!    4.5    16/07/98 Added code to get the broadcast flag for
!                    a unit, so that the broadcast of data to
!                    all the PE's can be surpressed if required -
!                    NUMOBS initially.
!                      Authors: Bob Carruthers & Deborah Salmond
!                               Cray Research
!LL  5.1    10/04/00 New reconfiguration support.
!    5.3    22/11/01 Enable MPP as the only option for
!                    small executables         E.Leung
!    5.5    25/04/03 Add defs FLUXPROC        E.Leung
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes
!
! Subroutine Arguments:

      INTEGER                                                           &
     &  NFT                                                             &
                     !   IN : FORTRAN unit number
     & ,LEN_IO                                                          &
                     !  OUT : no. of words read in
     & ,ISIZE        !   IN : no. of words to be read in

      REAL                                                              &
     &  IOSTAT                                                          &
                     !  OUT : Return code
     & ,ARRAY(ISIZE) !  OUT : Array to read data in to

! Parameters and Common blocks

#if !defined(RECON)
#include "parvars.h"
#endif

! Local variables

      INTEGER info
      REAL    return_codes(2)


! ------------------------------------------------------------------

      IOSTAT=-1.0
      LEN_IO=ISIZE


      IF (mype  ==  0) THEN
        CALL BUFFIN_SINGLE(NFT,ARRAY,ISIZE,LEN_IO,IOSTAT)
        return_codes(1)=LEN_IO
        return_codes(2)=IOSTAT
      ENDIF

!--get the broadcast flag
      call find_unit_bcast_flag(nft, info)
!--skip the broadcasts if the flag is set
      if(info == 0) then
      CALL GC_RBCAST(1,ISIZE,0,nproc,info,ARRAY)    ! data
      CALL GC_RBCAST(2,2,0,nproc,info,return_codes) ! return codes
      endif


      LEN_IO=return_codes(1)
      IOSTAT=return_codes(2)

      RETURN
      END SUBROUTINE BUFFIN

#endif
