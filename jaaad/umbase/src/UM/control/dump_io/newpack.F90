#if defined(C80_1A) || defined(UTILIO) || defined(RECON)               \
 || defined(VAROPSVER)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE NEWPACK:--------------------------------------
!LL
!LL  Purpose: Packing codes stored in LOOKUP(21,K) & LOOKUP(39,K)
!LL           are changed from pre vn2.8 values to
!LL           specification required at release 2.8
!LL
!LL  Written by A. Dickinson 28/08/92
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.2    13/07/93 Tidyied up integer declararions
!LL
!LL Programming standard :
!LL
!LL Logical components covered :
!LL
!LL Project task :
!LL
!LL  Documentation: UM Documentation Paper F3
!LL
!LLEND -----------------------------------------------------------------
!
      SUBROUTINE NEWPACK                                                &
     &(LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP)

      IMPLICIT NONE

      INTEGER                                                           &
     & LEN1_LOOKUP                                                      &
     &,LEN2_LOOKUP                                                      &
     &,LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)


      INTEGER                                                           &
     & N1                                                               &
     &,N2                                                               &
     &,N3                                                               &
     &,K

      DO K=1,LEN2_LOOKUP
        N1=0
        N2=0
        N3=0
        IF(LOOKUP(21,K) == -2)N1=2
! Ocean field packed using index array
        IF(LOOKUP(21,K) >  9.AND.LOOKUP(21,K) <  100)THEN
          N2=1
          N3=LOOKUP(21,K)-10
        ENDIF
! Ocean field compressed using bit mask
        IF(LOOKUP(21,K) >  99)THEN
          N2=2
          N3=LOOKUP(21,K)-100
        ENDIF
! Real field stored at land pts
        IF(LOOKUP(39,K) == 4)THEN
          LOOKUP(39,K)=1
          N2=2
          N3=1
        ENDIF
! Integer field stored at land pts
        IF(LOOKUP(39,K) == 5)THEN
          LOOKUP(39,K)=2
          N2=2
          N3=1
        ENDIF

        LOOKUP(21,K)=100*N3+10*N2+N1
        ENDDO

      RETURN
      END SUBROUTINE NEWPACK
#endif
