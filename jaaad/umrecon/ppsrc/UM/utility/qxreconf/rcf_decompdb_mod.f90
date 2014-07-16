
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Data module defining the decomp_db_type data-type

Module Rcf_DecompDB_Mod

! Description:
!   This module defines the decomp_db_type data-type - the type used
!   to stored information about various decompositions used in the
!   mpp reconfiguration.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.5   06/01/03   River routing support. P.Selwood.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Use Rcf_Parvars_Mod, Only : &
    Ndim_max,           &
    Maxproc

Use Rcf_DecompTP_Mod, Only : &
    Max_Decomps

Implicit None

Type Decomp_DB_type
  Integer      :: bound( Ndim_max )
  Integer      :: glsize( Ndim_max )
  Integer      :: glsizeu( Ndim_max )
  Integer      :: glsizev( Ndim_max )
  Integer      :: glsizer( Ndim_max )
  Integer      :: gridsize( Ndim_max )
  Integer      :: g_lasize( Ndim_max, 0 : Maxproc )
  Integer      :: g_blsizep( Ndim_max , 0 : Maxproc)
  Integer      :: g_blsizeu( Ndim_max , 0 : Maxproc)
  Integer      :: g_blsizev( Ndim_max , 0 : Maxproc)
  Integer      :: g_blsizer( Ndim_max , 0 : Maxproc)
  Integer      :: g_datastart( Ndim_max, 0 : Maxproc )
  Integer      :: g_datastartr( Ndim_max, 0 : Maxproc )
  Integer      :: g_gridpos( Ndim_max, 0 : Maxproc )
  Integer      :: halosize( Ndim_max )
  Integer      :: neighbour( 4 )
  Integer      :: first_comp_pe
  Integer      :: last_comp_pe
  Integer      :: nproc
  Integer      :: gc_proc_row_group
  Integer      :: gc_proc_col_group
  Integer      :: gc_all_proc_group
  Logical      :: set
End Type Decomp_DB_type

Type (Decomp_DB_type), Save  :: Decomp_DB( Max_Decomps )

Data Decomp_DB(:) % set / Max_Decomps*.FALSE./

End Module Rcf_DecompDB_Mod
