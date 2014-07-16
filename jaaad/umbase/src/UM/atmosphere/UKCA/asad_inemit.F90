#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
! Purpose: Dummy routine which the user must supply a new version if they
!     need to use and initialise their emissions scheme.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CINIT
!
! Current code owner: Glenn Carver/Colin Johnson
!                     Oliver Wild
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_INEMIT
        IMPLICIT NONE

!       include any commons here for getting any info
!       from the calling model

!       e.g. common /model/ rates(jpnl)

        INTEGER                       ::  ICODE

        CHARACTER (Len=80)            ::  CMESSAGE
        CHARACTER (Len=* ), Parameter ::  RoutineName='ASAD_INEMIT'

        CMESSAGE = 'Routine should not be callable'
        ICODE = 1
! DEPENDS ON: ereport
        CALL EREPORT(RoutineName,ICODE,CMESSAGE)

        RETURN
        END SUBROUTINE ASAD_INEMIT
#endif
