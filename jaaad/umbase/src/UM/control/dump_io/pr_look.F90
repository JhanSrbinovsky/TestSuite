#if defined(C80_1A) || defined(UTILIO) || defined(FLDC)                \
    || defined(FLDOP) || defined(RECON) || defined(VAROPSVER)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE PR_LOOK----------------------------------------
!LL
!LL  Purpose: Prints out Kth 64-word PP header
!LL
!LL AD, DR      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL   3.1  05/02/93    Trap use of user defined PPXREF file
!LL                    Author: A. Dickinson    Reviewer: R. Stratton
!     3.5  27/06/95  Submodels project. Replace call to RDPPXRF by
!                    call to GETPPX
!                    Author D.M.Goddard    Reviewer S Swarbrick
!LL   4.0  12/09/95    Change NPERIODS to LBUSER, LBUSER to LBPLEV
!LL                    as changes in CLOOKADD and PPHEAD1A.
!LL                    (Andrew Brady)
!LL  4.0  12/10/95  Chg. FORMAT, as last LBUSER is now MODEL_CODE. RTHB
!   4.0  24/08/95    Add dummy argument in call to GETPPX_REC
!                    Author D.M.Goddard    Reviewer S Swarbrick
!     4.1  18/06/96    Changes to cope with changes in STASH addressing
!                      Author D.M. Goddard.
!     5.1  10/04/00  New reconfiguration support. P.Selwood.
!     5.1  02/05/00  Improve layout of output. D. Robinson.
!     5.3  11/03/02  Change format for LBSRCE. D. Robinson
!     6.0  30/12/03  Set Model Code from 10 to 1 for FieldCalc
!                    diagnostics to use Atmos StashMaster file.
!                    D. Robinson
!     6.2  18/08/05  Allow RCF_EXPPX not to find STASHmaster entry.
!                                              R.Sharp
!
!LL  System component: R30/W30
!LL
!LL  System task: F3
!LL
!LL  Programming standard:  Unified Model Documentation Paper No 3
!LL                         Version No 1 15/1/90
!LL
!LL  Documentation:  Unified Model Documentation Paper No F3
!LL                  Version No 5 9/2/90
!LL
!LL------------------------------------------------------------
!*L Arguments:-------------------------------------------------
      SUBROUTINE PR_LOOK(                                               &
#include "argppx.h"
     &                   LOOKUP,RLOOKUP,LEN1_LOOKUP,K)
#if defined(RECON) || defined(VAROPSVER)
      Use Rcf_Exppx_Mod, Only : Rcf_Exppx

      Use Rcf_Ppx_Info_Mod, Only : STM_Record_type
#endif

      IMPLICIT NONE

      INTEGER LEN1_LOOKUP     ! IN First dimension of Look Up Table
      INTEGER K               ! IN Field number in Look Up Table
      INTEGER                                                           &
     & LOOKUP(LEN1_LOOKUP,*)  ! IN Integer equivalence of PP LOOKUP
      REAL                                                              &
     & RLOOKUP(LEN1_LOOKUP,*) ! IN Real equivalence of PP LOOKUP
#if defined(RECON) || defined(VAROPSVER)
! Dummy variables for the argppx.h argument comdeck
      INTEGER PPXI
      INTEGER PPXRECS
      CHARACTER PPXC
#endif

      CHARACTER*36 EXPPXC

! -------------------------------------------------------------
! Workspace usage:---------------------------------------------
! None
! -------------------------------------------------------------
!*L External subroutines called:-------------------------------
#if !defined(RECON) && !defined(VAROPSVER)
      EXTERNAL EXPPXC
#endif
!*-------------------------------------------------------------
! Comdecks: ------------------------------------------------------------
#if !defined(RECON) && !defined(VAROPSVER)
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#endif
#include "clookadd.h"

! Local variables:---------------------------------------------
      INTEGER ICODE             !Error code
      INTEGER ITEM              !STASH item number
      INTEGER SECTION           !STASH section number
      INTEGER MODEL             !Internal model number
      INTEGER I                 !Index

      CHARACTER*36 PHRASE       !Character part of PPXREF record
      CHARACTER*80 CMESSAGE     !Error message

#if defined(RECON) || defined(VAROPSVER)
      Type (STM_RECORD_TYPE), POINTER ::  RECORD
#endif

!--------------------------------------------------------------------


!L Write time and field type
        ITEM=MOD(LOOKUP(42,K),1000)
        SECTION=(LOOKUP(42,K)-ITEM)/1000
        MODEL=LOOKUP(45,K)

!       All diagnostics under model code of 10 are in section 20
!       of Atmos StashMaster file.
        IF ( MODEL == 10 ) Then
          MODEL = 1
        ENDIF
      ICODE = 0
#if defined(RECON) || defined(VAROPSVER)
      RECORD => RCF_EXPPX(MODEL,SECTION,ITEM,.TRUE.)
      IF ( ASSOCIATED ( RECORD ) ) THEN
        PHRASE = RECORD % NAME
      ELSE
        PHRASE = "Name not known"
      END IF
#else
! DEPENDS ON: exppxc
      PHRASE=EXPPXC(MODEL,SECTION,ITEM,                                 &
#include "argppx.h"
     &              ICODE,CMESSAGE)
#endif
        IF(ICODE /= 0)THEN
         PHRASE='NON-STANDARD FIELD'
        ENDIF

      WRITE (6,100) K,PHRASE

  100 FORMAT(' FIELD NO.', I5,4X,A)

      WRITE (6,101) LOOKUP(LBHR,K),LOOKUP(LBMIN,K),LOOKUP(LBDAT,K),     &
     &  LOOKUP(LBMON,K),LOOKUP(LBYR,K),LOOKUP(LBDAY,K),LOOKUP(LBHRD,K), &
     &  LOOKUP(LBMIND,K),LOOKUP(LBDATD,K),LOOKUP(LBMOND,K),             &
     &  LOOKUP(LBYRD,K),LOOKUP(LBDAYD,K)

  101 FORMAT(                                                           &
     &  ' VALID AT: ',  2I2.2,'Z  ',2(I2.2,'/'),I4.4,' DAY',I6,         &
     &  ' DATA TIME: ', 2I2.2,'Z  ',2(I2.2,'/'),I4.4,' DAY',I6)

!L Rest of header
      WRITE(6,200) (LOOKUP(I,K),I=13,45),(RLOOKUP(I,K),I=46,64)

  200 FORMAT(                                                           &
     &  '   LBTIM   LBFT    LBLREC LBCODE  LBHEM  LBROW  LBNPT',        &
     &  '  LBEXT LBPACK',/,                                             &
     &  1X,2I7,I10,6I7,/,                                               &
     &  '   LBREL   LBFC  LBCFC LBPROC   LBVC  LBRVC  LBEXP',           &
     &  '   LBBEGIN    LBNREC',/,                                       &
     &  1X,7I7,2I10,/,                                                  &
     &  '  LBPROJ  LBTYP  LBLEV LBRSVD LBRSVD LBRSVD LBRSVD   LBSRCE',/ &
     &  ,1X,7I7,I9,/,                                                   &
     &  '  DATA_TYPE     NADDR    LBUSER ITEM_CODE    LBPLEV',          &
     &  '    LBUSER MODEL_CODE',/                                       &
     &  1X,6I10,I11,/,                                                  &
     &  9X,'BULEV',7X,'BHULEV',5X,'BRSVD(3)',5X,'BRSVD(4)',             &
     &  7X,'BDATUM',/,                                                  &
     &  1X, 1P, 5E13.4,/,                                               &
     &  10X,'BACC',9X,'BLEV',8X,'BRLEV',8X,'BHLEV',7X,'BHRLEV',/,       &
     &  1X, 1P, 5E13.4, /,                                              &
     &  9X,'BPLAT',8X,'BPLON',9X,'BGOR',10X,'BZY',10X,'BDY',/,          &
     &  1X, 1P, 5E13.4, /,                                              &
     &  11X,'BZX',10X,'BDX',9X,'BMDI',9X,'BMKS',/,                      &
     &  1X, 1P, 4E13.4)

       WRITE(6,'('' '')')
      RETURN
      END SUBROUTINE PR_LOOK

#endif
