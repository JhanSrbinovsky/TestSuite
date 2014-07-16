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
! Purpose: Dummy heterogeneous chemistry routine.
!
!     The purpose of this routine is to set and return the heterogeneous
!     reaction rates. If the user has heterogeneous chemistry turned on
!     then this subroutine will be called. The user must supply their
!     own version of this routine to compute the heterogeneous rates.
!
!     Note that this subroutine is called repeatedly. It should not
!     therefore be used to do any I/O unless absolutely necessary. The
!     routine inihet is provided to initialise the heterogeneous chemist
!     by reading in files etc.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CDRIVE
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
        SUBROUTINE ASAD_HETERO( n_points )

        USE ASAD_MOD
        IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

        INTEGER, INTENT(IN) :: n_points   ! No of spatial points

!       Local variables

        CHARACTER (LEN=72) :: cmessage    ! Error message

!       include any commons here for getting required info
!       from the calling model

!       e.g. common /model/ het(jpnl)

!       1. Calculate heterogeneous rates
!          --------- ------------- -----

        cmessage='ERROR: YOU SHOULD NOT HAVE CALLED HETERO SUBROUTINE'
! DEPENDS ON: ereport
        CALL EREPORT('ASAD_HETERO',998,cmessage)

        RETURN
        END SUBROUTINE ASAD_HETERO
#endif
