#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Add inactive records to STASH list, when space is required
!
! Subroutine Interface:



!+Find whether ST_list entry (Im_ident,ISEC,ITEM) is an implied diag
! Subroutine Interface:

      SUBROUTINE IMPLIED                                                &
     &(Im_ident,ISEC,ITEM,LIMPLIED,ErrorStatus,CMESSAGE)
      IMPLICIT NONE
! Description:
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Mar. 95    Original code.  S.J.Swarbrick
!   4.5    09/12/97    Read the Implied data from PE 0 and
!                      distribute it.
!                        Author: Bob Carruthers
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:
#include "csubmodl.h"
#include "version.h"
#include "cstash.h"
#include "stextend.h"
#include "lenfil.h"
#include "chsunits.h"
#include "clfhist.h"
#if defined(MPP) && defined(T3E)
#include "parvars.h"

      integer info, msg
#endif

! Subroutine arguments:
!   Scalar arguments with intent(in):
      INTEGER Im_ident
      INTEGER ISEC
      INTEGER ITEM
      INTEGER ICODE,err ! return code

!   Scalar argument with intent(out):
      LOGICAL LIMPLIED   ! Set to T if diag is implied
      CHARACTER*80 CMESSAGE

! Local scalars:
      LOGICAL LSET ! Set to T when STASH_SET dir name has been obtained
      INTEGER I
      INTEGER J
      INTEGER IHOLD
      INTEGER N_IMPLICATORS
      CHARACTER*256 DIR

! Local arrays:
      INTEGER IMPLICS(100)

! ErrorStatus
      INTEGER ErrorStatus


! External subroutine calls
      INTEGER GETENV

! Function & Subroutine calls:
      External GET_FILE,FORT_GET_ENV
!- End of Header ----------------------------------------------------

      DATA LSET /.FALSE./
#if defined(MPP) && defined(T3E)
      save lset, file, ihold, dir
#else
      SAVE FILE,IHOLD,LSET
#endif

! Construction of file name for "STASH sets"
! (which specify implied diags)
!   On first call: assign directory name STASH_SET to FILE; add '/X'
#if defined(MPP) && defined(T3E)

      if(.not.lset) then
        dir       = ' '
        stash_set = ' '
        call fort_get_env('STASETS_DIR', 11, dir, 256, err)
#else

      DIR             =' '
      STASH_SET       =' '
! Correction for reading in the ~ctldata/stasets directory (N Farnon)
      CALL FORT_GET_ENV('STASETS_DIR',11,DIR,256,err)
#endif
            IF (err  /=  0) THEN
              WRITE(6,*) 'Warning: Environment variable STASETS_DIR has &
     &             not been set. Error code = ',err
#if defined(MPP) && defined(T3E)
              call abort()
#endif
            ENDIF
#if defined(MPP) && defined(T3E)
        stash_set=dir
#else
      STASH_SET=DIR

      IF(.NOT.LSET) THEN
#endif
          FILE=STASH_SET
          LSET=.TRUE.
          DO J=1,256
            IF (FILE(J:J) == ' ') THEN
            I = J
            GOTO 102
            END IF
          END DO
  102     FILE(I:I+1)='/X'
          I=I+2
          IHOLD=I
      END IF

!   Append rest of file name to FILE
      WRITE(FILE(IHOLD:IHOLD+9),501) Im_ident,ISEC,ITEM
  501 FORMAT(I2.2,2I3.3)

! Open STASH sets file; read diags listed in file into IMPLICS
! These diags are implied by the diag Im_ident, ISEC, ITEM
! Error message added (N.Farnon)
#if defined(MPP) && defined(T3E)
      if(mype == 0) then
#endif
      OPEN (3,FILE=FILE,IOSTAT=ICODE)
      IF (ICODE /= 0) THEN
        WRITE(6,*) 'Can not open stash_sets file, ICODE=',ICODE
        call abort()
      ELSE
        WRITE(6,*) 'OPEN: 3: ',FILE,': FILE EXISTS'
      END IF
      READ (3,600) N_IMPLICATORS
  600 FORMAT(I4)
#if defined(MPP) && defined(T3E)
      endif ! read on PE 0
!
!--send the number of implicators to each PE
      msg=7067
      call gc_ibcast(msg, 1, 0, nproc, info, n_implicators)
#endif
!
      if(n_implicators >  100) then
        write(6,*)'IMPLIED: Too Many Implicators',                      &
     &   ' - ',N_IMPLICATORS,' Requested, but there is Space',          &
     &   ' for only 100'
        call abort()
      endif
!
#if defined(MPP) && defined(T3E)
      if(mype == 0) then
#endif
      READ (3,610) (IMPLICS(I),I=1,N_IMPLICATORS)
  610 FORMAT(10I4)
#if defined(MPP) && defined(T3E)
      endif ! read on PE 0
!
      msg=7068
      call gc_ibcast(msg, n_implicators, 0, nproc, info,                &
     & implics)
#endif

!Find out whether any of the diags listed in FILE are present in LIST_S.
! (Any diag in the STASH list has a non-zero entry in SINDX). If one
!  or more of them are present, set LIMPLIED=T - indicating that
! diag Im_ident,ISEC,ITEM is implied.

      DO I=1,N_IMPLICATORS
        IF(INDX_S(2,Im_ident,ISEC,IMPLICS(I)) /= 0) THEN
          LIMPLIED=.TRUE.
          GO TO 9999
        END IF
      END DO

      LIMPLIED=.FALSE.
      CLOSE(UNIT=3)

 9999 RETURN
      END SUBROUTINE IMPLIED


!+Add diagnostic to the STASH list (LIST_S)
! Subroutine Interface:

#endif
