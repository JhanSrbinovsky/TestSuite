#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(UTILIO) || defined(FLDIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Gathers in atmosphere partial LBCs from boundaries to create
!  global LBC
!
!   5.5   25/04/03   Add defs UTILIO,FLDIO.  Remove MPP.  E.Leung
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes
! Subroutine Interface
      SUBROUTINE GATHER_ATMOS_LBCS()

!     This routine temporarily deleted while compile errors
!     are sorted out
      WRITE(6,*) 'Error: GATHER_ATMOS_LBCS called'
      CALL ABORT()

      RETURN
      END SUBROUTINE GATHER_ATMOS_LBCS
#endif
