#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE ARCHIVE(filetype_1,filetype_2,month,year)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : Archive files
!-
!-   Inputs  :
!-             MONTH,YEAR,      - time.
!-             FILETYPE_1       - letter for filename construction.
!-             FILETYPE_2       - letter for filename construction.
!-
!-   Outputs : ICODE,CMESSAGE
!-   Controls:
!-
!
! Current Owner of Code: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  4.4    20/07/98  Created.  C.E. Johnson
!  5.3    05/09/00  Filetype 's' added. C.E. Johnson
!  6.1    21/10/04  Calls UM_FORT_FLUSH instead of FLUSH directly.
!                   M.G. Sanderson.
!  6.2    06/04/05  Protected calls to UM_FORT_FLUSH.
!                                               M.G. Sanderson
!
!-
!VVV  V5.2.1  ARCHIVE  16/VIII.01  icode,cmessage not in arguments
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!-----------------------------------------------------------------------

      INTEGER,          INTENT(IN) :: month
      INTEGER,          INTENT(IN) :: year
      CHARACTER(LEN=1), INTENT(IN) :: filetype_1
      CHARACTER(LEN=1), INTENT(IN) :: filetype_2

      INTEGER :: icode
      CHARACTER(LEN=72) :: cmessage
      CHARACTER(LEN=14) :: filename
      CHARACTER(LEN=80) :: pipename

! DEPENDS ON: findname
      CALL FINDNAME(filetype_1,filetype_2,'c',0,0,filename,month,year)
      WRITE(6,*) 'Archiving: ',filename

      CALL GET_FILE(8,pipename,80,icode)
      OPEN(8,FILE=pipename)

      IF (filetype_1 == 'p') THEN      ! mean pp file
        WRITE(8,630) filename          ! archive pp
        WRITE(8,610) filename          ! delete
#if defined(T3E)
! DEPENDS ON: um_fort_flush
        CALL UM_FORT_FLUSH(8, icode)
#else
        CLOSE(8)
        OPEN(8, FILE=pipename)
#endif
      ELSE IF (filetype_1 == 't') THEN ! .tot file
        WRITE(8,620) filename          ! archive
!        WRITE(8,610) filename         ! delete
#if defined(T3E)
! DEPENDS ON: um_fort_flush
        CALL UM_FORT_FLUSH(8, icode)
#else
        CLOSE(8)
        OPEN(8, FILE=pipename)
#endif
      ELSE IF (filetype_1 == 's') THEN ! station output file
        WRITE(8,620) filename          ! archive
        WRITE(8,610) filename          ! delete
#if defined(T3E)
! DEPENDS ON: um_fort_flush
        CALL UM_FORT_FLUSH(8, icode)
#else
        CLOSE(8)
        OPEN(8, FILE=pipename)
#endif
      ELSE IF (filetype_1 == 'f') THEN ! cell following file
        WRITE(8,620) filename          ! archive
        WRITE(8,610) filename          ! delete
#if defined(T3E)
! DEPENDS ON: um_fort_flush
        CALL UM_FORT_FLUSH(8, icode)
#else
        CLOSE(8)
        OPEN(8, FILE=pipename)
#endif
      ELSE
        cmessage='Filetype_1 not recognised'
        WRITE(6,*) 'Filetype_1: ',filetype_1,' ',cmessage
! DEPENDS ON: ereport
        CALL EREPORT('ARCHIVE',1,cmessage)
      END IF

 610  FORMAT('%%% ',A14,' DELETE')
 620  FORMAT('%%% ',A14,' ARCHIVE DUMP')
 630  FORMAT('%%% ',A14,' ARCHIVE PPNOCHART')

      END SUBROUTINE ARCHIVE
#endif
