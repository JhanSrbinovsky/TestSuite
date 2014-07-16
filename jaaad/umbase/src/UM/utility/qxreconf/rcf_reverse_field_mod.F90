#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ reverse a N-S grid to give a S-N one

Module Rcf_Reverse_Field_Mod

! Description:
!              Reverse the position of "n_rows" rows for each level of
! the 3D data array passed into the routine such that the first row
! swaps place with the row at position n_rows. The 2nd row swaps place
! with th row at position n_rows-1 and so on. The lookup headers for
! the field are altered to be consistent with the new row ordering.
! Based on the routine - PF_REVERSE at UM4.5
!
! Method:
!
! Current Code Owner: Roddy Sharp
!
! History:
! Version   Date      Comment
! -------  --------   -------------------------
!  5.4     10/07/02   Original code. Roddy Sharp (frtz)
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! -Note- This routine has been overloaded to handle both real and
! logical fields. Please ensure you duplicate 'fixes' performed in one
! routine in the other routine where necessary

Interface Rcf_Reverse_Field
  Module Procedure Rcf_Reverse_Field_Log, Rcf_Reverse_Field_Real
End Interface

Contains

Subroutine Rcf_Reverse_Field_Real(RData,row_length,n_rows,levels,     &
                                  field_no,UM_Hdr)

Use Rcf_UMhead_Mod, Only :  &
    UM_Header_Type            ! Derived containing UM header info

Use Rcf_HeadAddress_Mod, Only  : &
    RC_FirstLat,       &
    RC_LatSpacing

Use EReport_Mod, Only :     &
    EReport

Implicit None
! Subroutine arguments

!< Scalar arguments with intent(in):>
Integer, Intent(In)              :: row_length  ! length of each row
Integer, Intent(In)              :: n_rows      ! no. of rows per level
Integer, Intent(In)              :: levels      ! no. of levels
Integer, Intent(In)              :: field_no    ! pos in lookup of
                                                ! first field

!< Array  arguments with intent(InOut):>
Real,    Intent(InOut)              :: RData(*) ! the data to be reversd
Type (UM_Header_Type),Intent(InOut) :: UM_Hdr   ! UM header info

! Comdecks
#include "clookadd.h"
! defines BDY and BZY

! Local variables

Integer                          :: Int_Val
Integer                          :: top_row
Integer                          :: bottom_row
Integer                          :: point
Integer                          :: lev_shift
Integer                          :: b_row_shift
Integer                          :: t_row_shift
Integer                          :: point_loc
Integer                          :: k
Real                             :: temp
Real                             :: Zeroth_Lat
Real                             :: Lat_Spacing
Real                             :: Real_Val

Integer                          :: ErrorStatus
Character(len=80)                :: cMessage
Character (Len=*), Parameter     :: RoutineName='Rcf_Reverse_Field_Real'

!=======================================================================
! Subroutine Code.
!=======================================================================
!loop over levels
Do k=1,levels

  lev_shift = (k-1) * (row_length * n_rows)

  Do top_row=1,n_rows/2

    bottom_row = n_rows + 1 - top_row

    b_row_shift = lev_shift + ((bottom_row - 1) * row_length)
    t_row_shift = lev_shift + ((top_row    - 1) * row_length)

    Do point=1,row_length

      ! Swap over corresponding point at top and bottom of data.
      temp                      = RData(b_row_shift + point)
      RData(b_row_shift + point) = RData(t_row_shift + point)
      RData(t_row_shift + point) = temp

    End Do ! point

  End Do ! top_row

  ! Now update the Latitude info held in Lookups
  Lat_Spacing = Transfer( UM_Hdr % Lookup( BDY, field_no ) , Real_Val)
  Zeroth_Lat  = Transfer( UM_Hdr % Lookup( BZY, field_no ) , Real_Val)

  ! Transpose Zeroth latitude.
  Zeroth_Lat  = Zeroth_Lat + ( (n_rows + 1) * Lat_Spacing )

  ! Sign of the latitude interval is changed.
  Lat_Spacing = - Lat_Spacing


  ! Put new values back into the Lookups
  UM_Hdr % Lookup( BDY, field_no ) = Transfer(Lat_Spacing , Int_Val)
  UM_Hdr % Lookup( BZY, field_no ) = Transfer(Zeroth_Lat  , Int_Val)

  ! double check the values in the Real Headers
  If (Um_Hdr % RealC (RC_FirstLat) /= ( Zeroth_Lat + Lat_Spacing)) Then
    Um_Hdr % RealC (RC_FirstLat) = ( Zeroth_Lat + Lat_Spacing)
  End If
  If (Um_Hdr % RealC (RC_LatSpacing) /= Abs(Lat_Spacing) ) Then
    !Um_Hdr % RealC (RC_LatSpacing) = Abs(Lat_Spacing)
    Write(cMessage,'(A)') "Latitude Spacing appears to have changed"
    ErrorStatus = 20
    Call EReport( RoutineName, ErrorStatus, Cmessage)

  End If


  ! ***************************************************
  ! At UM4.5 something was done here to re-align u or v
  ! onto the correct spacings.
  ! ***************************************************

End Do ! k over levels


Return

End Subroutine Rcf_Reverse_Field_Real


Subroutine Rcf_Reverse_Field_Log(LData,row_length,n_rows,levels,      &
                                  field_no,UM_Hdr)

Use Rcf_UMhead_Mod, Only :  &
    UM_Header_Type            ! Derived containing UM header info

Use Rcf_HeadAddress_Mod, Only  : &
    RC_FirstLat,       &
    RC_LatSpacing

Use EReport_Mod, Only :     &
    EReport

Implicit None
! Subroutine arguments

!< Scalar arguments with intent(in):>
Integer, Intent(In)              :: row_length  ! length of each row
Integer, Intent(In)              :: n_rows      ! no. of rows per level
Integer, Intent(In)              :: levels      ! no. of levels
Integer, Intent(In)              :: field_no    ! pos in lookup of
                                                ! first field

!< Array  arguments with intent(InOut):>
Logical, Intent(InOut)              :: LData(*) ! the data to be reversd
Type (UM_Header_Type),Intent(InOut) :: UM_Hdr   ! UM header info

! Comdecks
#include "clookadd.h"
! defines BDY and BZY

! Local variables

Integer                          :: Int_Val
Integer                          :: top_row
Integer                          :: bottom_row
Integer                          :: point
Integer                          :: lev_shift
Integer                          :: b_row_shift
Integer                          :: t_row_shift
Integer                          :: point_loc
Integer                          :: k
Real                             :: Zeroth_Lat
Real                             :: Lat_Spacing
Real                             :: Real_Val
Logical                          :: temp

Integer                          :: ErrorStatus
Character(len=80)                :: cMessage
Character (Len=*), Parameter     :: RoutineName='Rcf_Reverse_Field_Log'

!=======================================================================
! Subroutine Code.
!=======================================================================
!loop over levels
Do k=1,levels

  lev_shift = (k-1) * (row_length * n_rows)

  Do top_row=1,n_rows/2

    bottom_row = n_rows + 1 - top_row

    b_row_shift = lev_shift + ((bottom_row - 1) * row_length)
    t_row_shift = lev_shift + ((top_row    - 1) * row_length)

    Do point=1,row_length

      ! Swap over corresponding point at top and bottom of data.
      temp                      = LData(b_row_shift + point)
      LData(b_row_shift + point) = LData(t_row_shift + point)
      LData(t_row_shift + point) = temp

    End Do ! point

  End Do ! top_row

  ! Now update the Latitude info held in Lookups
  Lat_Spacing = Transfer( UM_Hdr % Lookup( BDY, field_no ) , Real_Val)
  Zeroth_Lat  = Transfer( UM_Hdr % Lookup( BZY, field_no ) , Real_Val)

  ! Transpose Zeroth latitude.
  Zeroth_Lat  = Zeroth_Lat + ( (n_rows + 1) * Lat_Spacing )

  ! Sign of the latitude interval is changed.
  Lat_Spacing = - Lat_Spacing


  ! Put new values back into the Lookups
  UM_Hdr % Lookup( BDY, field_no ) = Transfer(Lat_Spacing , Int_Val)
  UM_Hdr % Lookup( BZY, field_no ) = Transfer(Zeroth_Lat  , Int_Val)

  ! double check the values in the Real Headers
  If (Um_Hdr % RealC (RC_FirstLat) /= ( Zeroth_Lat + Lat_Spacing)) Then
    Um_Hdr % RealC (RC_FirstLat) = ( Zeroth_Lat + Lat_Spacing)
  End If
  If (Um_Hdr % RealC (RC_LatSpacing) /= Abs(Lat_Spacing) ) Then
    !Um_Hdr % RealC (RC_LatSpacing) = Abs(Lat_Spacing)
    Write(cMessage,'(A)') "Latitude Spacing appears to have changed"
    ErrorStatus = 20
    Call EReport( RoutineName, ErrorStatus, Cmessage)

  End If


  ! ***************************************************
  ! At UM4.5 something was done here to re-align u or v
  ! onto the correct spacings.
  ! ***************************************************

End Do ! k over levels

Return

End Subroutine Rcf_Reverse_Field_Log

End Module Rcf_Reverse_Field_Mod
#endif
