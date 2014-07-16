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
! Purpose: Emissions scheme
!
!     Purpose
!     -------
!     Dummy routine which the user must supply a new version if they
!     need to use an emissions scheme. This subroutine is called
!     by ASAD to set the emission rates.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from UKCA_MAIN1.
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
        SUBROUTINE ASAD_EMISSN

        USE ASAD_MOD
        IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

!     include any commons here for getting any info
!     from the calling model
!
!    e.g. common /model/ rates(jpnl)
!
!     ---------------------------------------------------------------

        RETURN
        END SUBROUTINE ASAD_EMISSN
#endif
