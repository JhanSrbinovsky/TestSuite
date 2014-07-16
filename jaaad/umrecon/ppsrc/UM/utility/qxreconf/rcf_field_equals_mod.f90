
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Copy one field into another - *IGNORING THE DATA*

Module Rcf_Field_Equals_Mod

!  Subroutine Rcf_Field_Equals - copys a field description into another
!
! Description:
!   Field size descriptions etc are copied from one field to another.
!   No data is copied or data pointers set.
!
! Method:
!   Stashmaster pointers are pointed to!
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.2   15/11/00   Add assignment of interp flag. P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Field_Equals( field_out, field_in )

Use Rcf_Field_Type_Mod, Only : &
    field_type

Implicit None

! Arguments
Type (field_type), Intent(In ) :: field_in
Type (field_type), Intent(Out) :: field_out

field_out % levels          = field_in % levels
field_out % rows            = field_in % rows
field_out % row_len         = field_in % row_len
field_out % level_size      = field_in % level_size
field_out % glob_rows       = field_in % glob_rows
field_out % glob_row_len    = field_in % glob_row_len
field_out % glob_level_size = field_in % glob_level_size
field_out % dump_pos        = field_in % dump_pos
field_out % interp          = field_in % interp
field_out % stashmaster    => field_in % stashmaster

Return
End Subroutine Rcf_Field_Equals
End Module Rcf_Field_Equals_Mod
