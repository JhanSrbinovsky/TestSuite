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
! Purpose: Dummy routine for initialising photolysis
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
        SUBROUTINE ASAD_INPHOT
        IMPLICIT NONE

        INTEGER                       ::  ICODE
        CHARACTER (Len=72)            ::  CMESSAGE
        CHARACTER (Len=*), Parameter  ::  RoutineName='ASAD_INPHOT'

        CMESSAGE = 'Dummy routine called'
        ICODE = -1
! DEPENDS ON: ereport
        CALL EReport(RoutineName,ICODE,CMESSAGE)

        RETURN
        END SUBROUTINE ASAD_INPHOT
#endif
