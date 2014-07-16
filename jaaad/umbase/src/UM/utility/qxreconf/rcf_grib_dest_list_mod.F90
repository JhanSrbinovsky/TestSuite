#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Deallocate memory used for a linked list of GRIB records

Module Rcf_Grib_Dest_List_Mod

! SUBROUTINE Rcf_Grib_Dest_List
!
! Description: Routine to DeAllocate all entries in a linked list
!              and nullify the pointers which locate the 'ends'.
!
! Method:
!         For lists with more than 1 member-
!           Step through the list deallocating the previous entry
!           Deallocate the end
!           Nullify the pointers to the end
!         For lists with only 1 member-
!           Deallocate the member.
!           Nullify the pointers.
!
! Current Code Owner: Roddy Sharp
!
! History:
! Version   Date      Comment
! -------  --------   -------------------------
!  5.4     17/07/02   Original code. Roddy Sharp (frtz)
!  6.2     20/06/05   Added deallocation of vertical
!                     coordinates arrays. Paul Earnshaw (frpe)
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
contains

Subroutine Rcf_Grib_Dest_List(List_Header)

Use Rcf_GRIB_Block_Params_Mod, Only : &
    List_Marker,         &
    Grib_Record

Use EReport_Mod, Only :     &
    EReport

Use Rcf_Parvars_mod, Only : &
    mype

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Diag                   ! =4 Extra Diagnostic output

Implicit None

! Subroutine arguments

!< Array  arguments with intent(InOut):>
Type (List_Marker)               :: List_Header

! Comdecks
#include "c_mdi.h"
#include "clookadd.h"
! contains LBLREC (amongst others)

! Local constants
Character (Len=*), Parameter     :: RoutineName='Rcf_GRIB_Dest_List'

! Local variables

Type (Grib_Record), Pointer      :: Current

Character (Len=80)               :: Cmessage(2)   ! used for EReport
Integer                          :: ErrorStatus   ! used for EReport

Integer                          :: Count

!=======================================================================
!  Routine Code Start :
!=======================================================================

! Check list has _some_ members.
If (Associated (List_Header % Begin )) Then

  Count = 0

  ! Lists with more than 1 member
  If ( List_Header % LstCount  > 1 ) Then
    Current => List_Header % Begin

    Do While (Associated (Current % Next))
      ! Deallocate array attached to pointer
      If (Associated(Current % VertCoords)) Then
        Deallocate(Current % VertCoords)
      End If
      Current => Current % Next
      Deallocate (Current % Prev)
      Count = Count + 1
    End Do
  End If

  ! do for all lists
  Nullify (Current)
  Deallocate (List_Header % End)

  If ( Count < List_Header % LstCount -1 ) Then
      Write(Cmessage(1),'(A,I2)') 'Destroyed less list entries than I&
                                  & thought existed'
      ErrorStatus = 10
      Call EReport( RoutineName, ErrorStatus, Cmessage(1) )
  Else
    If ( PrintStatus >= PrStatus_Diag .And. mype == 0) Then
      Write (6,*) "Destroyed a list containing ",count +1, " members."
    End If
  End If

End If ! list had members

Nullify (List_Header % Begin)
Nullify (List_Header % End)
List_Header % LstCount = 0

Return

End Subroutine Rcf_Grib_Dest_List
End Module Rcf_Grib_Dest_List_Mod
#endif
