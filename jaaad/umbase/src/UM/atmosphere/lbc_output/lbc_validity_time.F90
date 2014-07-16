#if defined(A32_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculate Validity Time for LBCs.
!
! Subroutine Interface:

      SUBROUTINE LBC_Validity_Time ( LCal360 )

      Implicit None

!
! Description:
!   Calculates the validity time for the LBCs.
!
! Method:
!   Takes the model basis time and derives the validity time
!   through a call to SEC2TIME.
!
! Current Code Owner: Dave Robinson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2    13/11/00   Original code. Dave Robinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables

#include "parvars.h"
#include "cmaxsize.h"
#include "csubmodl.h"
#include "ctime.h"
#include "cntlgen.h"
#include "parlbcs.h"

! Subroutine arguments

      Logical :: LCal360     !  360/365 day calander indicator

! Local scalars:

      Integer  :: yy,mm,dd,hr,mn,ss,day_no,sec

! Function & Subroutine calls:

      External SEC2TIME

!- End of header

! ---------------------------------------------
! Determine Validity Time from Model Basis Time
! ---------------------------------------------

      sec = STEPim(a_im) * SECS_PER_PERIODim(a_im) /                    &
     &          STEPS_PER_PERIODim(a_im)

! DEPENDS ON: sec2time
      Call SEC2TIME(0,SEC,BASIS_TIME_DAYS,BASIS_TIME_SECS,              &
     &                  YY,MM,DD,HR,MN,SS,DAY_NO,LCAL360)

! ------------------------------
! Store validity time in parlbcs
! ------------------------------
      LBC_VT_Year  = YY
      LBC_VT_Month = MM
      LBC_VT_Day   = DD
      LBC_VT_Hour  = HR
      LBC_VT_Min   = MN
      LBC_VT_DayNo = Day_No

      return
      END SUBROUTINE LBC_Validity_Time
#endif
