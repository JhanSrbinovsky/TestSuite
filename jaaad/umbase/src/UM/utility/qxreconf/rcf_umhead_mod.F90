#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Data module defining um_header_type

MODULE Rcf_UMhead_Mod

! Description:
!  Contains parameters used in UM dumps,
!  a structure to hold a UM header
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.2   13/03/01   Use 2 dimensional level dep constants.  P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

!   Global Constants
INTEGER, PARAMETER :: LenFixHd = 256   ! Length of Fixed Length Header.

!   Global Type Definitions:
TYPE UM_header_type
  SEQUENCE

  INTEGER ::      LenIntC         ! Length of Integer Constants array.
  INTEGER ::      LenRealC        ! Length of Real Constants array.
  INTEGER ::      Len1LevDepC     ! 1st dimension \  of Level Dependent
  INTEGER ::      Len2LevDepC     ! 2nd dimension /  Constants array.
  INTEGER ::      Len1RowDepC     ! 1st dimension \  of Row Dependent
  INTEGER ::      Len2RowDepC     ! 2nd dimension /  Constants array.
  INTEGER ::      Len1ColDepC     ! 1st dimension \  of Column Dependent
  INTEGER ::      Len2ColDepC     ! 2nd dimension /  Constants array.
  INTEGER ::      Len1FldsOfC     ! 1st dimension \  of Fields of
  INTEGER ::      Len2FldsOfC     ! 2nd dimension /  Constants array.
  INTEGER ::      LenExtraC       ! Length of Extra Constants array.
  INTEGER ::      LenHistFile     ! Length of History File.
  INTEGER ::      LenCompFldI1    ! Length of Compressed Field Index 1.
  INTEGER ::      LenCompFldI2    ! Length of Compressed Field Index 2.
  INTEGER ::      LenCompFldI3    ! Length of Compressed Field Index 3.
  INTEGER ::      Len1Lookup      ! 1st dimension of Lookup table.
  INTEGER ::      Len2Lookup      ! 2nd dimension of Lookup table.
  INTEGER ::      LenData         ! Length of Data array.
  INTEGER ::      StartData       ! Position of start of Data array.
  INTEGER ::      NumFlds         ! Number of data fields.
  INTEGER ::      MaxFldSize      ! Maximum size of field.
  INTEGER ::      UnitNum         ! Unit number associated with UM dump.

  INTEGER, POINTER ::    FixHd    (:)  ! Fixed length header.
  INTEGER, POINTER ::    IntC     (:)  ! Integer Constants array.
  INTEGER, POINTER ::    CompFldI1(:)  ! Compressed Field Index array 1.
  INTEGER, POINTER ::    CompFldI2(:)  ! Compressed Field Index array 2.
  INTEGER, POINTER ::    CompFldI3(:)  ! Compressed Field Index array 3.
  INTEGER, POINTER ::    Lookup(:,:)   ! Lookup table.

  REAL, POINTER ::       RealC   (:)  ! Real Constants array.
  REAL, POINTER ::       LevDepC (:,:)! Level Dependent Constants array.

!Use 2-dimentional Row and Column dependent consts.
  REAL, POINTER ::       RowDepC (:,:)! Row Dependent Const. array.
  REAL, POINTER ::       ColDepC (:,:)! Column Dependent Const. array.

  REAL, POINTER ::       FldsOfC (:)  ! Field Dependent Constants array.
  REAL, POINTER ::       ExtraC  (:)  ! Extra Constants array.
  REAL, POINTER ::       HistFile(:)  ! History File.

END TYPE UM_header_type

END MODULE Rcf_UMhead_Mod
#endif
