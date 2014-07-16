#if defined(FLUXPROC) || defined(FLXPLPR)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
! History:
! version  date         change
! 4.5      03/09/98     New code
!
! Author:     M. J. Bell
!----------------------------------------------------------------------
! contains routines: add_hours,amend_times
!
! Purpose:  Flux processing routines.
!           add_hours adds an input number of hours
!           to an input date/time (called Clim1)
!           to produce an output date/time (called Valid).
!           amend_times amends date information in Int_Head.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine amend_times (                                          &
#include "avaltim.h"
     &                      Int_head,Len_IntHd )

      implicit none

! declaration of parameters
#include "clookadd.h"

! declaration of globals used
#include "cvaltim.h"
#include "cvaloff.h"

! declaration of argument list
      integer Len_IntHd
      integer Int_Head(Len_IntHd) ! IN/OUT integer part of lookup table

! declarations for validity time

      integer IAddHrs !    used to find validity time
      integer Year1   !    First year in header
      integer Month1  !    First Month in header
      integer Day1    !    First Day in header
      integer Hour1   !    First Hour in header
      integer Min1    !    Always equal to zero
      integer Sec1    !    Always equal to zero

! no other variables used

! declaration of externals
      external add_hours
!----------------------------------------------------------------------

! 1. Set the second time in header
      Int_Head(LBYRD) = ValidYear
      Int_Head(LBMOND) = ValidMonth
      Int_Head(LBDATD) = ValidDay
      Int_Head(LBHRD) = ValidHour

! 2. Set the first time in header

      IAddHrs = 0 - ValidityPeriod

! DEPENDS ON: add_hours
      call add_hours(                                                   &
#include "avaltim.h"
     &       Year1,Month1,Day1,Hour1,Min1,Sec1,                         &
     &       IAddHrs )

      Int_Head(LBYR)  = Year1
      Int_Head(LBMON) = Month1
      Int_Head(LBDAT) = Day1
      Int_Head(LBHR)  = Hour1

      return
      END SUBROUTINE amend_times
!----------------------------------------------------------------------
#endif
