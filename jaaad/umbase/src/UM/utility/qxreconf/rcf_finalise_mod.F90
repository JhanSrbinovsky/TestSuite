#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ finalise the reconfiguration - ties up loose ends

Module Rcf_Finalise_Mod

!   Subroutine Rcf_Finalise - final tidyup tasks
!
! Description:
! This module performs tidy up tasks including closing files and
! calling Timer for output.
!
! Method:
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Finalise( hdr_in, hdr_out )

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_FortranIO_Mod, Only : &
    Rcf_Free_Unit

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Parvars_mod, Only : &
    mype

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid

Use Rcf_Recon_Mod, Only : &
    Grib

Use Rcf_PrintStatus_Mod, Only : &
    LTimer,                     &
    L_IO_Timer

Implicit None

! Arguments
Type (um_header_type), Intent(InOut)  :: hdr_in
Type (um_header_type), Intent(InOut)  :: hdr_out

! Local variables
Integer                               :: err
Integer                               :: ErrorStatus
Character (Len=*), Parameter          :: RoutineName='Rcf_Finalise'
Character (Len=80)                    :: Cmessage

!-------------------------------------------------------------------
! Close open files
!-------------------------------------------------------------------
err = 0
! Output dump
! DEPENDS ON: file_close
Call File_Close( hdr_out % UnitNum, 'ASTART', 6, 0, 0, err)

If ( err /= 0 ) Then
  Cmessage    = 'Failed to Close Output Dump'
  ErrorStatus = 10
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

Call Rcf_Free_Unit( hdr_out % UnitNum )

! Input dump
If (Input_Grid % Rotated .AND. .NOT. Grib ) Then
! DEPENDS ON: file_close
  Call File_Close( hdr_in % UnitNum, 'RECONTMP', 8, 0, 0, err)
Else
! DEPENDS ON: file_close
  Call File_Close( hdr_in % UnitNum, 'AINITIAL', 8, 0, 0, err)
End If

If ( err /= 0 ) Then
  Cmessage    = 'Failed to Close Input Dump'
  ErrorStatus = 20
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

Call Rcf_Free_Unit( hdr_in % UnitNum )

!-------------------------------------------------------------------
! Produce Timer output
!-------------------------------------------------------------------
If (LTimer) Then
! DEPENDS ON: timer
  Call Timer("Reconfigure",2)
End If

If (L_IO_Timer .AND. mype == 0) Then
  Write (6,*) 'IO timings for PE 0'
  Call io_total_timings()
End If


Return
End Subroutine Rcf_Finalise
End Module Rcf_Finalise_Mod
#endif
