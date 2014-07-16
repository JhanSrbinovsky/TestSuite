#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Top level for reading in reconfiguration namelists.

MODULE Rcf_read_namelists_Mod

!  Subroutine Rcf_Read_Namelists - read the rcf namelists.
!
! Description:
!   Read the namelists, assigning and freeing Fortran units as
!   required.
!
! Method:
!   The namelists are provided in one file - RECONA and RECONO
!   for Atmos and Ocean respectively.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

Subroutine rcf_read_namelists ( )

Use Rcf_Submodel_Mod, Only : &
    Submodel_Ident,      &
    Atmos_IM

Use Rcf_readnl_nsubmodl_Mod,  Only :  Rcf_readnl_nsubmodl
Use Rcf_readnl_ustsnum_Mod,   Only :  Rcf_readnl_ustsnum
Use Rcf_readnl_uancnum_Mod,   Only :  Rcf_readnl_uancnum
Use Rcf_readnl_stshcomp_Mod,  Only :  Rcf_readnl_stshcomp
Use Rcf_readnl_nlstcatm_Mod,  Only :  Rcf_readnl_nlstcatm
Use Rcf_readnl_Recon_Mod,     Only :  Rcf_readnl_recon
Use Rcf_readnl_vertical_Mod,  Only :  Rcf_readnl_vertical
Use Rcf_readnl_horizont_Mod,  Only :  Rcf_readnl_horizont
Use Rcf_readnl_headers_Mod,   Only :  Rcf_readnl_headers
Use Rcf_readnl_items_Mod,     Only :  Rcf_readnl_items
Use Rcf_readnl_trans_Mod,     Only :  Rcf_readnl_trans

Use Ereport_Mod, Only :&
    Ereport

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

Use Rcf_FortranIO_Mod, Only : &
    Rcf_Get_Unit,             &
   Rcf_Free_Unit,             &
   Max_Filename_Len

Use Rcf_Recon_Mod, Only : &
    TRANS

Implicit none
! Arguments

! Local Variables/Paramters
Character (Len=*), Parameter       :: RoutineName = 'Rcf_Read_Namelists'
Character (Len=80)                 :: Cmessage
Character (Len=Max_Filename_Len)   :: FileName
Character (Len=Max_Filename_Len)   :: FileName_vertlevs
Integer                            :: ErrorStatus
Integer                            :: nft
Integer                            :: nft_vertlevs
Integer                            :: nft_horizgrid
Integer                            :: status
Logical                            :: l_exist

! Get a Fortran Unit for the namelists
Call Rcf_Get_Unit( nft )
Call Rcf_Get_Unit( nft_vertlevs )
Call Rcf_Get_Unit( nft_horizgrid )

! First, find the namelist filename from Env Vars
Call Fort_Get_Env( 'RCF_NAMELIST', 12, FileName,                   &
                    Max_Filename_Len, ErrorStatus )

If ( ErrorStatus /= 0 ) Then
  ErrorStatus = 10
  Cmessage = 'Unable to Obtain Namelist Filename from Environment'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

FileName = Trim( FileName )

! Check file exists - Do all this by hand (not using file_open)
! as wish to do this on *all* PEs together
Inquire( file=FileName, exist=l_exist )

If ( .Not. l_exist ) Then
  ErrorStatus = 20
  Cmessage = 'Namelist file does not exist!'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! Open the file containing namelists from RECONA/RECONO file
Open( Unit=nft, File=FileName, IOstat=status )

If ( PrintStatus >= PrStatus_Oper ) Then
  Write (6,*) 'Namelist file: ',FileName
End If

Call rcf_readnl_recon    (nft)
Call rcf_readnl_nsubmodl (nft)
Call rcf_readnl_ustsnum  (nft)
Call rcf_readnl_uancnum  (nft)
Call rcf_readnl_stshcomp (nft)

! Submodel dependent namelists
If (Submodel_Ident == Atmos_IM) Then
  Call rcf_readnl_nlstcatm (nft)
End If

! ------------------------
! Vertical Levels Namelist
! ------------------------

If (Submodel_Ident == Atmos_IM) Then

! Find the filename containing vertical levels from Env Vars
  Call Fort_Get_Env('VERT_LEV', 8, FileName_vertlevs,              &
                     Max_Filename_Len, ErrorStatus)

  If ( ErrorStatus /= 0 ) Then
    ErrorStatus = 30
    Cmessage =  &
    'Unable to Obtain Vertical Levels Filename from Environment'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

  FileName_vertlevs = Trim( FileName_vertlevs )

! Check file exists - Do all this by hand (not using file_open)
! as wish to do this on *all* PEs together
  Inquire( file=FileName_vertlevs, exist=l_exist )

  If ( .Not. l_exist ) Then
    Write (6,*) 'Vertical Levels file: ',FileName_vertlevs
    ErrorStatus = 40
    Cmessage = 'Vertical Levels Namelist file does not exist!'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

! Open the file containing vertical levels
  Open( Unit=nft_vertlevs, File=FileName_vertlevs, IOstat=status )

  If ( PrintStatus >= PrStatus_Oper ) Then
    Write (6,*) 'Vertical Levels file: ',FileName_vertlevs
  End If

End If

! For vertical, two namelists - VERTICAL and VERTLEVS - are read in.
Call rcf_readnl_vertical (nft, nft_vertlevs)
Call rcf_readnl_horizont (nft, nft_horizgrid)
Call rcf_readnl_headers  (nft)
Call rcf_readnl_items    (nft)

If (TRANS) Then
  Call rcf_readnl_trans    (nft)
End If

Close( Unit=nft )
Close( Unit=nft_vertlevs )

Call Rcf_Free_Unit( nft )
Call Rcf_Free_Unit( nft_vertlevs )

Return

END SUBROUTINE  Rcf_read_namelists

END MODULE Rcf_read_namelists_Mod
#endif
