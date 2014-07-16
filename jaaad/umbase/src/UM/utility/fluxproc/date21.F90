#if defined(PPTOANC) || defined(FLUXPROC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!*********************************************************************
!                                                                     *
!   RECORD OF CHANGES:                                                *
!   ==================                                                *
!                                                                     *
!  0 ORIGINAL VERSION BY JOHN PRINCE, LONG AGO IN THE MISTS OF TIME,  *
!                                     FOR THE IBM WITH STATIC MEMORY  *
!  1 TRANSLATED BY FPP 2.26B16 11/12/89  11:44:56   TDYON=X           *
!  2 RE-WRITE   BY PAUL WHITE  13/12/89                               *
!                        TO MAKE RE-ENTRANT WITH INLINE EXPANSION     *
!                        TRANSLATED BY FPP 2.26B16 13/12/89  11:41:08 *
!                        SWITCHES: LSTOFF=T,OPTON=78,TDYON=FX         *
!  3 RE-COMPILED 29/10/91 TO PRODUCE 31 BIT ADDRESSING MODE VERSION   *
!        BY M. COLLIER - COPIED TO MET.PROGLIB                        *
!  4 Updated 30/1/98 by Edward Jones                                  *
!                        Update ZPDATE subroutine                     *
!                        Added ISALEAP subroutine                     *
!                        Ported to HP, Cray and PC from MET.SRCELIB   *
!  5 Updated 17/2/98 by Edward Jones                                  *
!                        Added DATCHK and MNTHDS Routines             *
!  6 Updated 23/3/98 by Edward Jones                                  *
!                        Added JDAY Routine                           *
!  7fre  Updated 17/4/98 by Stephen Turner                            *
!                        Converted to FREE  format F90, and added     *
!                              Zeller method                          *
!  7fix  Updated 01/9/98 by Stephen Turner                            *
!                        Converted to FIXED format F90 and removed    *
!                        (irrelevant) Zeller method                   *
!  7fix_nomods Updated 6/10/98 by Stephen Turner                      *
!                        Converted to FIXED format F90 without        *
!                        using modules (and without Zeller method)    *
!                                                                     *
!**********************************************************************
!-----------------------------------------------------------------------






!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      SUBROUTINE DATE21 (IDY, IY, ICD)
!
!     DAYS SINCE 1.1.1900, FROM DAY OF YEAR
!
      IMPLICIT NONE
!                        INPUT ARGUMENTS
      INTEGER, INTENT(IN)  :: IDY, IY
!                       OUTPUT ARGUMENTS
      INTEGER, INTENT(OUT) :: ICD
!                       LOCAL VARIABLES
      INTEGER :: IYN

      IYN = IY - 1900
      IF (IYN  >   0) THEN
         ICD = IDY + IYN*365 + (IYN-1)/4 - (IYN-1)/100 + (IYN+299)/400
      ELSE
         ICD = IDY + IYN*365 + IYN/4 - IYN/100
      ENDIF

      END SUBROUTINE DATE21
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! The JDATE Conversion algorithms are based on the algorithm published
! in a letter to the editor of Communications of the ACM (CACM, volume 1
! number 10, October 1968, p.657) by Henry F. Fliegel and
! Thomas Van Flandern
! This algorithm is valid only for dates from
! 1/3/-4900 G onward when converting from a Julian day number to a date,
! or from 1/3/-4800 when converting from a date to a Julian day number.
! It should be noted that these algorithms are valid only in the
! Gregorian Calendar and the Proleptic Gregorian Calendar (after the
! dates given above). They do not handle dates in the Julian Calendar.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!----------------------------------------------------------
#endif
