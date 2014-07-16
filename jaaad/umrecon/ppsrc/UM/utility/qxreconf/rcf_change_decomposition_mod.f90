
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Parallel Reconfiguration: Select a new decomposition

Module Rcf_Change_Decomposition_Mod

!  Subroutine  Rcf_Change_Decompostions - swap decomposition variables
!
! Description:
! Sets up the Rcf_Parvars module with the correct information for
! decomposition new_decomp
!
! Method:
! If new_decomp is already the current decomposition, exit and do
! nothing.
! If decomposition new_decomp has not been initialised, print a
! message and exit
! Otherwise, copy the information from the decomp_db arrays in the
! DECOMPDB module into the PARVARS module. Also swaps land-sea mask
! pointers appropriately.
!
! Adapted from UM 4.5 code.
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

Contains

SUBROUTINE Rcf_Change_Decomposition( new_decomp )

Use Rcf_Parvars_Mod

Use Rcf_DecompTP_Mod

Use Rcf_DecompDB_Mod, Only : &
    Decomp_DB

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Lsm_Mod

IMPLICIT NONE

! Subroutine arguments:
Integer, Intent( In)         :: new_decomp  ! new decomposition to use

! Local variables
Integer                      :: ErrorStatus
Integer                      :: ineb
Integer                      :: idim
Integer                      :: iproc
Character (Len=*), Parameter :: RoutineName='Rcf_Change_Decomposition'
Character (Len=80)           :: Cmessage


! ------------------------------------------------------------------

ErrorStatus=0

! Check that the new_decomp argument is sensible
IF ((new_decomp .GT. max_decomps) .OR. &
   ((new_decomp .LT. 1) .AND. (new_decomp .NE. decomp_unset)))  THEN
  IF (mype .EQ. 0) THEN
    WRITE(6,*) 'Error: Cannot change to decomposition ', &
                new_decomp
    WRITE(6,*) 'This decomposition does not exist'
    WRITE(6,*) 'Exiting.'
  ENDIF
  ErrorStatus=10
  Cmessage = 'Required decomposition does not exist'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
ENDIF

! Check if this is already the current decomposition

IF (new_decomp .EQ. current_decomp_type) GOTO 999

! Check to see if setting decomposition to unset

IF (new_decomp .EQ. decomp_unset) THEN
  current_decomp_type = decomp_unset
  GOTO 999
ENDIF

! Check if this decomposition has been initialised

IF ( .NOT. Decomp_DB( new_decomp ) % set ) THEN
  IF (mype .EQ. 0) THEN
    WRITE(6,*) 'Error : Attempt to select uninitialised ', &
               'decomposition ',new_decomp
    WRITE(6,*) 'Exiting.'
  ENDIF
  ErrorStatus=20
  Cmessage = 'Cannot change to uninitialise decomposition'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
ENDIF

! Now we can copy the information into PARVARS

first_comp_pe=Decomp_DB( new_decomp ) % first_comp_pe
last_comp_pe=Decomp_DB( new_decomp ) % last_comp_pe

nproc=Decomp_DB( new_decomp ) % nproc
nproc_x=Decomp_DB( new_decomp ) % gridsize(1)
nproc_y=Decomp_DB( new_decomp ) % gridsize(2)

Offx=Decomp_DB( new_decomp ) % halosize(1)
Offy=Decomp_DB( new_decomp ) % halosize(2)

gc_proc_row_group=Decomp_DB( new_decomp ) % gc_proc_row_group
gc_proc_col_group=Decomp_DB( new_decomp ) % gc_proc_col_group
gc_all_proc_group=Decomp_DB( new_decomp ) % gc_all_proc_group

DO ineb=1,4
  neighbour(ineb)=Decomp_DB( new_decomp ) % neighbour(ineb)
ENDDO

DO idim=1,Ndim_max
  bound(idim)=Decomp_DB( new_decomp ) % bound(idim)
  glsize(idim)=Decomp_DB( new_decomp ) % glsize(idim)
  glsizeu(idim)=Decomp_DB( new_decomp ) % glsizeu(idim)
  glsizev(idim)=Decomp_DB( new_decomp ) % glsizev(idim)
  glsizer(idim)=Decomp_DB( new_decomp ) % glsizer(idim)
  gridsize(idim)=Decomp_DB( new_decomp ) % gridsize(idim)
ENDDO

DO iproc=first_comp_pe,last_comp_pe
  DO idim=1,Ndim_max
    g_lasize(idim,iproc)=  &
                  Decomp_DB( new_decomp ) % g_lasize(idim,iproc)
    g_blsizep(idim,iproc)= &
                   Decomp_DB( new_decomp ) % g_blsizep(idim,iproc)
    g_blsizeu(idim,iproc)= &
                   Decomp_DB( new_decomp ) % g_blsizeu(idim,iproc)
    g_blsizev(idim,iproc)= &
                   Decomp_DB( new_decomp ) % g_blsizev(idim,iproc)
    g_blsizer(idim,iproc)= &
                   Decomp_DB( new_decomp ) % g_blsizer(idim,iproc)
    g_datastart(idim,iproc)= &
                     Decomp_DB( new_decomp ) % g_datastart(idim,iproc)
    g_datastartr(idim,iproc)= &
                     Decomp_DB( new_decomp ) % g_datastartr(idim,iproc)
    g_gridpos(idim,iproc)=  &
                   Decomp_DB( new_decomp ) % g_gridpos(idim,iproc)
  ENDDO
ENDDO

DO idim=1,Ndim_max
  lasize(idim)=g_lasize(idim,mype)
  blsizep(idim)=g_blsizep(idim,mype)
  blsizeu(idim)=g_blsizeu(idim,mype)
  blsizev(idim)=g_blsizev(idim,mype)
  blsizer(idim)=g_blsizer(idim,mype)
  datastart(idim)=g_datastart(idim,mype)
  datastartr(idim)=g_datastartr(idim,mype)
  gridpos(idim)=g_gridpos(idim,mype)
ENDDO

atSouth = ( gridpos(2) .EQ. 0)
atNorth = ( gridpos(2) .EQ. (gridsize(2)-1))
atEast  = ( gridpos(1) .EQ. (gridsize(1)-1))
atWest  = (gridpos(1) .EQ. 0)

! Additional info for reconfiguration - for LandSeaMask
If (new_decomp == decomp_rcf_input ) Then
  glob_atmos_landmask       => glob_lsm_in
  local_atmos_landmask      => local_lsm_in
  local_land_field          => local_land_in
  glob_land_field           => glob_land_in
Else If (new_decomp == decomp_rcf_output) Then
  glob_atmos_landmask       => glob_lsm_out
  local_atmos_landmask      => local_lsm_out
  local_land_field          => local_land_out
  glob_land_field           => glob_land_out
End If

current_decomp_type=new_decomp

999  CONTINUE

RETURN
END Subroutine Rcf_Change_Decomposition

End Module Rcf_Change_Decomposition_Mod



