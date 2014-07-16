
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Extract a STASHmaster record

Module Rcf_Exppx_Mod

!  Subroutine Rcf_Exppx   - extract a STASHmaster record.
!
! Description:
!   Given a model, section and item this routine returns a pointer
!   to the specified STASHmaster record
!
! Method:
!   The relevant row in the record array (STM_record) is found
!   by the ppxptr index.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   6.2   18/08/05   Add optional argument to allow non-fatal result on
!                    not finding a current STASHmaster entry. R.Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Function Rcf_Exppx( Im_ident, section, item, NoFindArg )

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Ppx_Info_mod, Only : &
    STM_record_type,     &
    STM_record,          &
    ppxptr

IMPLICIT NONE
! Function arguments:
Integer, Intent(In) :: Im_ident    ! Internal model identifier
Integer, Intent(In) :: section     ! STASH section no.
Integer, Intent(In) :: item        ! STASH item no.
Logical, Optional, Intent(In) :: NoFindArg ! Warning with no entry found

! Function value (out)
Type (STM_record_type), Pointer :: Rcf_Exppx


! Local variables
Character (Len=*), Parameter :: RoutineName = 'Rcf_Exppx'
Character (Len=80)           :: Cmessage
Integer                      :: row         ! Row no. in PPXI array
Integer                      :: ErrorStatus !+ve = fatal error
Logical                      :: NoFind  !local value of NoFindArg

If ( Present (NoFindArg)) Then
  NoFind = NoFindArg
Else
  NoFind = .FALSE.
End If

Nullify( Rcf_Exppx )

ErrorStatus = 0
If ( Im_ident <= 0 .OR. section < 0 .OR. item <= 0) Then
  If (Im_ident <= 0) Write(6,*) 'Exppxi: INVALID Im_ident'
  If (section  <  0) Write(6,*) 'Exppxi: INVALID SECTION NO.'
  If (item     <=  0) Write(6,*) 'Exppxi: INVALID ITEM NO.'
  Write(6,*) &
      'Im_ident ',Im_ident,' section ',section,' item ',item
  ErrorStatus=1
  Cmessage='ERROR Exppx: INVALID STASH RECORD ID'

  Call Ereport( RoutineName, ErrorStatus, Cmessage )

Else

  ! Obtain row no. in PPXI array
  row = PPXPTR(Im_ident,section,item)

  ! Obtain required data value
  If ( row > 0 ) Then
    Rcf_Exppx => STM_record( row )
  Else

    If ( NoFind ) Then
      Write (Cmessage,*)  'STASH entry not defined : item ', item,     &
                          ' section ', section,' model ',  Im_ident
      ErrorStatus = -2
    Else
      Write (Cmessage,*)  'Cant find required STASH item ', item,      &
            ' section ', section,' model ',  Im_ident, ' in STASHmaster'
      ErrorStatus = 2
    End If

    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If
End If

Return
End Function Rcf_Exppx

End Module Rcf_Exppx_Mod
