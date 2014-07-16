
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Defines the STASHmaster record

Module Rcf_Ppx_Info_Mod

! Description:
!   Defines the STASHmaster record format, the USTSNUM namelist
!   and other related variables
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   03/10/01   Addition of Implicit None. R.Sharp
!   5.5   10/02/03   Extension of option code to 30 digits. T.White
!   6.1   06/10/04   Increase size of NDiagP - A.A. Dickinson
!   6.2   08/06/06   Increase size of NDiagP again to 2500. R Barnes
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Use Rcf_Submodel_Mod, Only  : &
    N_Internal_Model_Max

Implicit None

! A parameter needed for the following definition.
! Number of packing codes in stashmaster record.
Integer, Parameter        :: Ppxref_pack_profs = 10

! Length of Option and Version Codes.
Integer, Parameter        :: STM_OptCodeLen = 30 !Must be multiple of 5
Integer, Parameter        :: STM_VerCodeLen = 20

! For broadcast: quantity of integer/character data in each STMrecord.
! Must be changed if STMrecord format changes.
Integer, Parameter        :: STM_IntDataLen  = 28 + Ppxref_pack_profs &
                                             + STM_OptCodeLen / 5
Integer, Parameter        :: STM_CharDataLen = 37

! Define STMrecord type - to hold STASHmaster records
! Note that as these are in a sequence and are integer followed by
! character, we can do some fairly efficient comms with them.
! Thus order *is* vital!
Type STM_record_type
  Sequence
  Integer                 :: model
  Integer                 :: section
  Integer                 :: item
  Integer                 :: space_code
  Integer                 :: ptr_code
  Integer                 :: timavail_code
  Integer                 :: grid_type
  Integer                 :: lv_code
  Integer                 :: lb_code
  Integer                 :: lt_code
  Integer                 :: pt_code
  Integer                 :: pf_code
  Integer                 :: pl_code
  Integer                 :: lev_flag
  Integer                 :: opt_code( STM_OptCodeLen / 5 )
  Integer                 :: version_mask
  Integer                 :: halo_type
  Integer                 :: data_type
  Integer                 :: dump_packing
  Integer                 :: packing_acc(Ppxref_pack_profs)
  Integer                 :: rotate_code
  Integer                 :: field_code
  Integer                 :: user_code
  Integer                 :: lbvc_code
  Integer                 :: base_level
  Integer                 :: top_level
  Integer                 :: ref_lbvc_code
  Integer                 :: cf_levelcode
  Integer                 :: cf_fieldcode
  Integer                 :: RowIndex

  Character (Len=36)      :: name
  Character (Len=1)       :: OriginFlag
End Type STM_record_type

!---------------------------------------------------------------
! A few parameters for STASH related sizing
!---------------------------------------------------------------
Integer, Parameter        :: NDiagP            = 2600  ! max size of

                                                       ! STASHmaster
Integer, Parameter        :: NItemP            = 999   ! Max number of
                                                       ! items per sect
Integer, Parameter        :: NSectP            = 99    ! Max # of sects.
Integer, Parameter        :: Ppxref_Sections   = NSectP - 55
Integer, Parameter        :: Ppxref_Items      = NItemP


!-------------------------------------------------------------------
! Guts for the storage of STASHmaster information
!-------------------------------------------------------------------
Integer, Save            :: ppxRecs = 0     ! No. of stash records

Type (STM_record_type), Pointer, Save   :: STM_record(:)

! Referencing for the above
Integer, Save            :: ppxptr( N_Internal_Model_Max,              &
                                 0:Ppxref_Sections , Ppxref_Items ) = 0

! Info for user STASHmaster file
Integer, Save            :: n_ustash
Integer, Save            :: nrecs_ustash
Character (Len=80), Save :: ustsfils(20)

Namelist/ustsnum /n_ustash, nrecs_ustash, ustsfils

End Module Rcf_Ppx_Info_mod
