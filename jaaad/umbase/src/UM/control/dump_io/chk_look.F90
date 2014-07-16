#if defined(C80_1A) || defined(UTILIO) || defined(RECON)               \
 || defined(VAROPSVER)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE CHK_LOOK---------------------------------------
!LL
!LL  Written by A. Dickinson
!LL
!LL  Purpose: Cross checks pointers in PP LOOKUP records with
!LL           model parameters
!LL
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3
!LL                        Version No 1 15/1/90
!LL
!LL  System component: C25
!LL
!LL  System task: F3
!LL
!LL  Documentation: Unified Model Documentation Paper No F3
!LL                 Version No 5 9/2/90
!LL
!LL------------------------------------------------------------
!*L Arguments:-------------------------------------------------
      SUBROUTINE CHK_LOOK(FIXHD,LOOKUP,LEN1_LOOKUP,                     &
                                                        ! Intent (In)
     &                    LEN_DATA,                                     &
                                                        !
#include "argppx.h"
     &                    ICODE,CMESSAGE)               ! Intent (Out)

      IMPLICIT NONE

      INTEGER                                                           &
     & LEN_DATA                                                         &
                       !IN Length of model data
     &,LEN1_LOOKUP                                                      &
                       !IN First dimension for LOOKUP
     &,LOOKUP(LEN1_LOOKUP,*)                                            &
                             !IN Integer equivalence of PP LOOKUP
     &,FIXHD(*)                                                         &
                       !IN Fixed length header
     &,ICODE          !OUT Return code; successful=0
                      !                 error > 0

      CHARACTER*80                                                      &
     & CMESSAGE       !OUT Error message if ICODE > 0

#if defined(RECON)
! Dummy variables for the argppx.h comdeck declarations
      INTEGER PPXI
      INTEGER PPXRECS
      CHARACTER PPXC
#endif
! Comdecks:----------------------------------------------------------
#include "cppxref.h"
#if !defined(RECON)
#include "csubmodl.h"
#include "ppxlook.h"
#endif
#include "clookadd.h"

! -------------------------------------------------------------
! Workspace usage:---------------------------------------------
! None
! -------------------------------------------------------------
!*L External subroutines called:-------------------------------
      EXTERNAL PR_LOOK
!--------------------------------------------------------------
! Local variables:---------------------------------------------
      INTEGER                                                           &
     & K                                                                &
               ! Loop count
     &,LEN_D                                                            &
                 ! Cumulative length of data
     &,N1
!--------------------------------------------------------------

!L Internal structure: None

      ICODE=0
      CMESSAGE=' '

! CHK_LOOK falls over with Boundary Datasets if pre-3.4
      IF (FIXHD(5) == 5 .and. FIXHD(12) <  304) THEN
        write (6,*) ' CHK_LOOK skipped for Boundary Dataset (Pre 3.4)'
      ELSE

      LEN_D=0

      DO 100 K=1,FIXHD(152)

#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)                \
 || defined(CONVPP)
      if(lookup(1,k) == -99) goto 9986
#endif
! Check that data_type is valid no: 1 to 3 or -1 to -3
      IF((LOOKUP(DATA_TYPE,K) >= 1.AND.LOOKUP(DATA_TYPE,K) <= 3) .OR.   &
     &   (LOOKUP(DATA_TYPE,K) <= -1.AND.LOOKUP(DATA_TYPE,K) >= -3))     &
     &   THEN
        LEN_D=LEN_D+LOOKUP(LBLREC,K)
      ELSE
       WRITE(6,'('' *ERROR* Wrong value of'',I9,'' in LOOKUP(DATA_'',   &
     &''TYPE'',I4,'')'')')LOOKUP(DATA_TYPE,K),K
! DEPENDS ON: pr_look
      CALL PR_LOOK(                                                     &
#include "argppx.h"
     &             LOOKUP,LOOKUP,LEN1_LOOKUP,K)
       ICODE=1
       CMESSAGE='CHK_LOOK: Consistency check'
       RETURN
      ENDIF

100   CONTINUE

#if !defined(MPP)
      IF(LEN_DATA /= LEN_D.OR.FIXHD(161) /= LEN_DATA.OR.                &
     &FIXHD(161) /= LEN_D)THEN

       WRITE(6,'('' *ERROR* Length of model data specified wrongly      &
     & : PARAMETER='',I9,''FILE='',I9,''FIXHD(161)'',I9)')              &
     & LEN_DATA,LEN_D,FIXHD(161)
       WRITE(6,*)' Your initial dump may need reconfiguring.'
       ICODE=2
       CMESSAGE='CHK_LOOK: Consistency check - try reconfiguring your in&
     &itial dump'
       RETURN
      ENDIF
#endif

#if defined(CONVIEEE) || defined(CUMF) || defined(PUMF)                \
 || defined(CONVPP)
9986  continue
#endif

      ENDIF  !  Check on pre-3.4 Boundary Datasets
      RETURN
      END SUBROUTINE CHK_LOOK

#endif
