#if defined(A26_1A) || defined(A26_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine diagnostics_riv

      Subroutine diagnostics_riv(                                       &
     &                       row_length, rows                           &
     &,                       river_row_length, river_rows              &
     &,                      at_extremity                               &
     &,                      at_extremity_riv                           &
     &,                      RIVEROUT                                   &
     &,                      BOX_OUTFLOW, BOX_INFLOW                    &
! Add inland basin outflow to arguments
     &,                      TWATSTOR,INLANDOUT_RIV                     &

     &,                                                                 &
#include "argsts.h"
     & STASHwork                                                        &
     & )

! Description:
!   Calculates river-related diagnostics (held in STASH section 26).
!
! Method:
!   Each diagnostic is simply copied into the STASHwork array
!   to be passed on to STASH for output processing.
!
!   Diagnostics currently available (in order calculated):
!   Item  Description
!    1    River water storage  (river grid)
!    2    gridbox outflow       (   "     )
!    3    gridbox runoff        (   "     )
!    4    coastal outflow       (ATMOS grid)
!         6    inland basin outflow       (river grid)

! History:
! Version   Date     Comment
! ----     -------     -------
! 5.5      28/03/03  Original code. C. Bunton
! 6.0      12/09/03  Change DEFs from A20 to A26. D Robinson
! 6.0      29/07/03  Move diagnostics from section 20 to section 26.
!                    D.Robinson
!   6.2   21/2/06  Re-route outflow from inland basins to soil moisture
!                  P. Falloon

! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
!     &                   ! south, east or west of the processor grid
     &,  at_extremity_riv(4) ! Indicates if this processor is at north,
!     &                   ! south, east or west of the processor grid
! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, river_row_length                                                &
                          ! river row length
     &, river_rows        ! river rows

      LOGICAL                                                           &
     & L_RIVERS           ! IN rivers switched on if .TRUE.

      REAL                                                              &
     & RIVEROUT(row_length, rows)                                       &
     &,BOX_OUTFLOW(river_row_length, river_rows)                        &
     &,BOX_INFLOW(river_row_length, river_rows)                         &
     &,TWATSTOR(river_row_length, river_rows)                           &

! Declare inland basin outflow variable
     &,INLANDOUT_RIV(river_row_length, river_rows)


#include "c_mdi.h"

! Local variables

      Integer                                                           &
     &  i, j, k, l                                                      &
     &,    icode                ! Return code  =0 Normal exit  >1 Error

      Integer, Parameter :: Sect = 26  !  Section No for RR diagnostics

      Character*80  cmessage
      Character(*) RoutineName
      Parameter ( RoutineName='diagnostics_riv')

      Integer                                                           &
     &  im_index        ! internal model index

      Real                                                              &
     &  interp_data(row_length,rows)

#include "csubmodl.h"
#include "typsts.h"

! Diagnostic variables
       Real                                                             &
     &  STASHwork(*)    ! STASH workspace

      External                                                          &
     &  copydiag, copydiag_3d                                           &
     &  ,Ereport

! ------------------------------------------------------------------
! Section 1.  Initialisation.
! ------------------------------------------------------------------


      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! --------------------------------------------------------------------
! River outflow (at seapoints)
! --------------------------------------------------------------------
! Item 26 004 riverout
      IF(sf(004,sect))THEN
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(004,sect,im_index)),riverout,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,004,                                        &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 004)"
            goto 9999
         End if
      ENDIF

! ----------------------------------------------------------------------
! River outflow (at seapoints)
! -------------------------------------------------------------------
! Item 26 001 riverout
      IF(sf(001,sect))THEN
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(001,sect,im_index)),TWATSTOR,       &
     &        river_row_length,river_rows,0,0,0,0, at_extremity_riv,    &
     &        atmos_im,sect,001,                                        &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 001)"
            goto 9999
         End if
      ENDIF
!----------------------------------------------------------------------
! --------------------------------------------------------------------
! River outflow (at seapoints)
! --------------------------------------------------------------------
! Item 26 002 riverout
      IF(sf(002,sect))THEN
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(002,sect,im_index)),BOX_OUTFLOW,    &
     &        river_row_length,river_rows,0,0,0,0, at_extremity_riv,    &
     &        atmos_im,sect,002,                                        &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 002)"
            goto 9999
         End if
      ENDIF
!-------------------------------------------------------------------
! ------------------------------------------------------------------
! River outflow (at seapoints)
! ------------------------------------------------------------------
! Item 26 003 riverout
      IF(sf(003,sect))THEN
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(003,sect,im_index)), BOX_INFLOW,    &
     &        river_row_length,river_rows,0,0,0,0, at_extremity_riv,    &
     &        atmos_im,sect,003,                                        &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 003)"
            goto 9999
         End if
      ENDIF
!---------------------------------------------------------------------

! Output inland basin outflow on TRIP grid

! ------------------------------------------------------------------
! Inland basin outflow
! ------------------------------------------------------------------
! Item 26 006 inlandout_riv
      IF(sf(006,sect))THEN
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(006,sect,im_index)),                &
     &    INLANDOUT_RIV,                                                &
     &        river_row_length,river_rows,0,0,0,0, at_extremity_riv,    &
     &        atmos_im,sect,006,                                        &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 006)"
            goto 9999
         End if
      ENDIF
!---------------------------------------------------------------------

 9999 continue
      If(icode /= 0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,Cmessage)
      Endif

      Return
      END SUBROUTINE diagnostics_riv

#endif
