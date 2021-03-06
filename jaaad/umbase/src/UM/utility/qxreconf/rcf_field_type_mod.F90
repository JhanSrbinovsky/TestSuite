#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Defines the field_type data-type

Module Rcf_Field_Type_Mod

! Description:
!   Data module defining the field_type data-type. This should contain
!   all information commonly required of a field such as sizes,
!   stashmaster, dump position, etc
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

Use Rcf_Ppx_Info_Mod, Only : &
    STM_record_type

Implicit None

Type field_type

  ! data field - two dimensions x*y and z
  ! For some situations (interpolation of integer field for e.g.),
  ! more than one of these fields may be used at a given time.
  Real, Pointer    :: Data( :, : )
  Integer, Pointer :: Data_Int( :, : )
  Logical, Pointer :: Data_Log( :, : )

  ! Sizes
  Integer          :: levels

  ! Local sizes are default
  Integer          :: rows
  Integer          :: row_len
  Integer          :: level_size           ! Size of single level

  ! Global sizes are also required
  Integer          :: glob_rows
  Integer          :: glob_row_len
  Integer          :: glob_level_size

  ! Information about the sort of field it is
  Integer          :: dump_pos             ! start position in dump
  Integer          :: interp               ! Should it be/has it been
                                           ! interpolated/copied
  ! Pointer to relevant stashmaster field
  Type (STM_record_type), Pointer :: stashmaster    ! stashmaster

End Type field_type

End Module Rcf_Field_Type_Mod

#endif
