#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Gathers atmosphere LBCs from all processors.

Module Rcf_Gather_Atmos_Lbcs_Mod

!  Subroutine Rcf_Gather_Atmos_LBCs
!
! Description:
!   Gathers a distributed atmosphere LBC field onto a single PE.
!
! Method:
!   Calculates remote sizes and gathers using GC_RSEND/RRECV.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.



Contains

Subroutine Rcf_Gather_Atmos_Lbcs(Full_LBC,     Full_LBC_Size,   &
                                 Part_LBC,     Part_LBC_Size,   &
                                 Stash_Record, Gather_pe)

Use Rcf_Parvars_Mod, Only : &
    mype,               &
    nproc,              &
    glsize

Use Rcf_Recon_Mod, Only : &
    Rimwidtha

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Ppx_Info_Mod, Only : &
    STM_Record_Type

Use Rcf_Level_Code_Mod, Only : &
    Rcf_Level_Code

Use Rcf_Address_Length_Mod, Only : &
    Rcf_Address_Length

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Implicit None

! Arguments
Integer, Intent(In)                 :: Full_LBC_Size
Integer, Intent(Out)                :: Part_LBC_Size
Integer, Intent(In)                 :: Gather_pe

Real, Intent(Out)                   :: Full_LBC( Full_LBC_Size )
Real, Intent(In)                    :: Part_LBC( * )

Type( STM_Record_Type ), Intent(In) :: Stash_Record

! Comdecks
#include "cppxref.h"
#include "gccom.h"

! Local variables
Integer                      :: i
Integer                      :: global_size
Integer                      :: local_size
Integer                      :: local_extra
Integer                      :: part_size
Integer                      :: start_level
Integer                      :: end_level
Integer                      :: istat
Integer                      :: ErrorStatus
Character (Len=*), Parameter :: RoutineName='Rcf_Gather_Atmos_LBCs'
Character (Len=80)           :: Cmessage


! Global size of LBC - found from Rcf_Address_Length( input and
! output grids same size )

Call Rcf_Address_Length( Stash_Record % Grid_Type,             &
                         Stash_Record % Halo_Type, global_size )

! Need to multiply this size by the number of levels
Call Rcf_Level_Code( Stash_Record % lb_code, start_level,  Output_Grid )
Call Rcf_Level_Code( Stash_Record % lt_code, end_level,    Output_Grid )

global_size = global_size * ( (end_level - start_level) + 1)


! Work out how much data per PE now.

local_size = global_size / nproc
local_extra = Mod( global_size, nproc )

! Use gc_rsend and gc_rrecv - could use gcg_ralltoall but this is
! simpler to work out

If (mype /= nproc - 1) Then
  Part_size = local_size
Else
  Part_size = local_size + local_extra
End If

! Now send to Processor gather_pe on all pes
Call Gc_rsend( 200, Part_size, gather_pe, istat, &
               Full_LBC( mype * local_size + 1), Part_LBC )
If (istat /= GC_OK) Then
  ErrorStatus = 30
  Cmessage = 'Failure in Gcom - gc_rsend'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

Call Gc_Ssync( nproc, istat)

! Now recv on processor gather_pe
If (mype == gather_pe) Then
  Do i = 0, nproc - 2
    Call Gc_rrecv( 200, local_size, i, istat, &
                   Full_LBC( i * local_size + 1 ), Part_LBC )
    If (istat /= GC_OK) Then
      ErrorStatus = 10
      Cmessage = 'Failure in Gcom - gc_rrecv'
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If
  End Do

  ! Treat nproc - 1 seperately
  i = nproc - 1
  Call Gc_rrecv( 200, local_size + local_extra, i, istat, &
                 Full_LBC( i * local_size + 1 ), Part_LBC )
  If (istat /= GC_OK) Then
    ErrorStatus = 20
    Cmessage = 'Failure in Gcom - gc_rrecv'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

End If

Call Gc_Ssync( nproc, istat)

! Set the return size
part_LBC_size = local_size
If (mype == nproc - 1) part_LBC_size = part_LBC_size + local_extra

Return
End Subroutine Rcf_Gather_Atmos_Lbcs

End Module Rcf_Gather_Atmos_Lbcs_Mod


#endif
