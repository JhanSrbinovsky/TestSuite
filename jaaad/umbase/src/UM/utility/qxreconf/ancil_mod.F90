#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
Module Ancil_Mod

!
! Description:
! Define structure to hold Ancillary FIELD information read in
! from ANCILmaster records

! Current Code Owner: D.Robinson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.3   03/10/01   Addition of headers and Implicit none. R.Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
      Implicit none
!
Type ANC_record_type
  Sequence
  Integer                  :: ancil_ref_number
  Integer                  :: model_number
  Integer                  :: section_number
  Integer                  :: item_number
  Integer                  :: anc_file_number
  Integer                  :: anc_field_read

  Character (Len=36)       :: anc_name
  Character (Len=1)        :: anc_flag
End Type ANC_record_type

Integer, Parameter         :: AncRec_IntDataLen  = 5
Integer, Parameter         :: AncRec_CharDataLen = 36

Type (ANC_record_type), Allocatable, Save  :: anc_record(:)

Integer, Save              :: ancRecs     = 0     ! No. of ancil records
Integer, Save              :: max_ancRecs = 0 ! Max no of ancil records


! Define structure to hold Ancillary FILE information read in
! from ANCILmaster records

Type ANC_file_type
  Sequence
  Integer                  :: anc_file_number
  Integer                  :: model_number
  Integer                  :: anc_file_open
  Character (Len=8)        :: anc_env_var
  Character (Len=41)       :: anc_file_title
  Character (Len=1)        :: anc_flag
End Type ANC_file_type

Integer, Parameter         :: AncFile_IntDataLen  = 2
Integer, Parameter         :: AncFile_CharDataLen = 49

Type (ANC_file_type), Allocatable, Save  :: anc_file(:)

Integer , Save :: ancFiles     = 0     ! No of anc files
Integer , Save :: max_ancFiles = 0     ! Max no of anc files
Integer , Save :: AncF_UnitNo          ! Unit No for ancillary files


! Namelist for user ANCILmaster file

Integer, Save              :: n_uancil
Character (Len=80), Save   :: uancfils(20)


Namelist /uancnum/ &
n_uancil, uancfils


! Allocatable work arrays for ancillary processing

Integer, dimension (:), allocatable :: nlookup
Integer, dimension (:), allocatable :: lookup_step
Integer, dimension (:), allocatable :: stashancil
! Integer, dimension (:), allocatable :: fieldcode
Integer, dimension (:), allocatable :: levels
Integer, dimension (:), allocatable :: ancil_add

End Module Ancil_Mod

#endif
