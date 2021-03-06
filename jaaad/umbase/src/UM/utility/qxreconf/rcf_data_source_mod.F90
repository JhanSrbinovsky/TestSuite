#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Definitions for data-type defining source of a field

Module Rcf_data_source_Mod

! Description:
!   Data module defining the data_source_type data-type and
!   related magic numbers.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   03/10/01   Addition of Implicit None. R.Sharp
!   5.3   11/10/01   Alter to use Max_Filename_Len parameter. R.Sharp
!   5.4   29/08/02   Add new magic number -99. D.Robinson
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Use Rcf_FortranIO_Mod, Only : &
    Max_Filename_Len

Implicit None

Type data_source_type

Integer                          :: Source
Integer                          :: Domain
Integer                          :: Ancil_SctnC
Integer                          :: Ancil_ItemC
Real                             :: RConst
Character (Len=Max_Filename_Len) :: Ancil_File

End Type data_source_type

! magic numbers for a source of a field etc
Integer, Parameter     :: Input_Dump     = 1
Integer, Parameter     :: Ancillary_File = 2
Integer, Parameter     :: Set_To_Zero    = 3
Integer, Parameter     :: Set_To_MDI     = 4
Integer, Parameter     :: Tracer_File    = 5
Integer, Parameter     :: Set_To_Const   = 6
Integer, Parameter     :: External_Dump  = 7
Integer, Parameter     :: Field_Calcs    = 8

! magic number for internal use only
Integer, Parameter     :: Already_Processed = -99

! Domain values
Integer, Parameter     :: Whole_Grid     = 1
Integer, Parameter     :: Sub_Grid       = 2

End Module Rcf_data_source_Mod
#endif
