
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Parallel RCF: sets array of tids of adjacent processes.

Module Rcf_Set_Neighbour_Mod

!  Subroutine Rcf_Set_Neighbour - What are my neighbouring processors?
!
! Description:
! This routine finds the tids of the North, South, East and West
! neighbouring processes. It takes account of the boundary
! condition in each dimension (X and Y) which can be either:
! cyclic : wrap around
! static : no wrap around
!
! Method:
! The tid of each neighbouring process is calculated (taking into
! account the relevant boundary conditions) and placed in the
! neighbour array.
!
! Derived from UM4.5 code.
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


Contains

SUBROUTINE Rcf_Set_Neighbour(decomp_type)

Use Rcf_Parvars_Mod
Use Rcf_DecompTP_Mod
Use Rcf_DecompDB_Mod

IMPLICIT NONE

! Subroutine Arguments:
INTEGER decomp_type  ! Decomposition type to update neighbours for

! ------------------------------------------------------------------

! Set Northern neighbour

IF ( decomp_DB (decomp_type) % g_gridpos (2,mype) .GT. 0) THEN
!       This processor is not at the top of the LPG
  decomp_DB (decomp_type) % neighbour (PNorth) = &
    mype - decomp_DB (decomp_type) % gridsize (1)
ELSEIF (decomp_DB (decomp_type) % bound (2) .EQ. BC_CYCLIC) THEN
!       This processor at the top of LPG, and has cyclic BCs.
  decomp_DB (decomp_type) % neighbour (PNorth) = &
    mype + decomp_DB (decomp_type) % nproc  - &
    decomp_DB (decomp_type) % gridsize (1)
ELSE
!       This processor at top of LPG and has static BCs
  decomp_DB (decomp_type) % neighbour (PNorth) =NoDomain
ENDIF

! Set Southern neighbour

IF ( decomp_DB (decomp_type) % g_gridpos (2,mype) .LT. &
    (decomp_DB (decomp_type) % gridsize (2)-1) ) THEN
!       This processor is not at the bottom of the LPG
  decomp_DB (decomp_type) % neighbour (PSouth) = &
    mype + decomp_DB (decomp_type) % gridsize (1)
ELSEIF (decomp_DB (decomp_type) % bound (2) .EQ. BC_CYCLIC) THEN
!       This processor at the bottom of LPG, and has cyclic BCs.
  decomp_DB (decomp_type) % neighbour (PSouth) = &
    mype - decomp_DB (decomp_type) % nproc  + &
    decomp_DB (decomp_type) % gridsize (1)
ELSE
!       This processor at top of LPG and has static BCs
  decomp_DB (decomp_type) % neighbour (PSouth) =NoDomain
ENDIF

! Set Western neighbour

IF ( decomp_DB (decomp_type) % g_gridpos (1,mype) .GT. 0) THEN
!       This processor is not at the left of the LPG
  decomp_DB (decomp_type) % neighbour (PWest) = &
    mype - 1
ELSEIF (decomp_DB (decomp_type) % bound (1) .EQ. BC_CYCLIC) THEN
!       This processor at the left of the LPG, and has cyclic BCs.
  decomp_DB (decomp_type) % neighbour (PWest) = &
    mype + decomp_DB (decomp_type) % gridsize (1) - 1
ELSE
!       This processor at top of LPG and has static BCs
  decomp_DB (decomp_type) % neighbour (PWest) =NoDomain
ENDIF

! Set Eastern neighbour
IF ( decomp_DB (decomp_type) % g_gridpos (1,mype) .LT. &
    (decomp_DB (decomp_type) % gridsize (1)-1) ) THEN
!       This processor is not at the right of the LPG
  decomp_DB (decomp_type) % neighbour (PEast) = &
    mype + 1
ELSEIF (decomp_DB (decomp_type) % bound (1) .EQ. BC_CYCLIC) THEN
!       This processor at the left of the LPG, and has cyclic BCs.
  decomp_DB (decomp_type) % neighbour (PEast) = &
    mype - decomp_DB (decomp_type) % gridsize (1) + 1
ELSE
!       This processor at top of LPG and has static BCs
  decomp_DB (decomp_type) % neighbour (PEast) =NoDomain
ENDIF

RETURN
END Subroutine Rcf_Set_Neighbour
End Module Rcf_Set_Neighbour_Mod


