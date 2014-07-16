#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(UTILIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Parallel UM: Sets up array of tids of adjacent processes.
!
! Subroutine interface:
      SUBROUTINE SET_NEIGHBOUR(decomp_type)

      IMPLICIT NONE

!
! Description:
! This routine finds the tids of the North, South, East and West
! neighbouring processes. It takes account of the boundary
! condition in each dimension (X and Y) which can be either:
! cyclic : wrap around
! overpole : data moves from processor on opposite side of pole
! static : no wrap around
!
! Method:
! The tid of each neighbouring process is calculated (taking into
! account the relevant boundary conditions) and placed in the
! neighbour array.
!
! Current Code Owner: Paul Burton
!
! History:
!  Model    Date     Modification history from model version 3.5
!  version
!    3.5    4/1/95   New DECK created for the Parallel Unified
!                    Model. P.Burton + R.Skaalin
!    4.2    21/08/96 Updated to use DECOMPDB variables rather than
!                    PARVARS. Only argument required is the
!                    decomposition type.              P.Burton
!    5.0    15/04/99 Added "BC_OVERPOLE" boundary condition.
!                                                      P.Burton
!    5.2    06/02/01  fix bug for cyclic LAM              A. Malcolm
!    5.3    22/11/01  Enable MPP as the only option for
!                     small executables         E.Leung
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes

!
! Subroutine Arguments:

      INTEGER decomp_type  ! Decomposition type to update neighbours
!                          ! for

! Parameters and Common blocks

#include "parvars.h"
#include "decomptp.h"
#include "decompdb.h"

! ------------------------------------------------------------------


! Set Northen neighbour

      IF ( decomp_db_g_gridpos(2,mype,decomp_type)  <                   &
     &    (decomp_db_gridsize(2,decomp_type)-1) ) THEN
!       This processor is not at the North of the LPG
        decomp_db_neighbour(PNorth,decomp_type) =                       &
     &    mype + decomp_db_gridsize(1,decomp_type)
      ELSEIF (decomp_db_bound(2,decomp_type)  ==  BC_CYCLIC) THEN
!       This processor at the top of LPG, and has cyclic BCs.
        decomp_db_neighbour(PNorth,decomp_type) =                       &
     &    mype - decomp_db_nproc(decomp_type) +                         &
     &    decomp_db_gridsize(1,decomp_type)
      ELSEIF (decomp_db_bound(2,decomp_type)  ==  BC_OVERPOLE) THEN
! This processor passes data over the pole
         IF ((decomp_db_g_gridpos(1,mype,decomp_type)+1)  <=            &
     &       (decomp_db_gridsize(1,decomp_type)/2)) THEN
           decomp_db_neighbour(PNorth,decomp_type) =                    &
     &       mype+decomp_db_gridsize(1,decomp_type)/2
         ELSE
           decomp_db_neighbour(PNorth,decomp_type) =                    &
     &       mype-decomp_db_gridsize(1,decomp_type)/2
         ENDIF
       ELSE
!       This processor at top of LPG and has static BCs
        decomp_db_neighbour(PNorth,decomp_type) =NoDomain
      ENDIF

! Set Southern neighbour

      IF ( decomp_db_g_gridpos(2,mype,decomp_type)  >   0) THEN
!       This processor is not at the South of the LPG
        decomp_db_neighbour(PSouth,decomp_type) =                       &
     &    mype - decomp_db_gridsize(1,decomp_type)
      ELSEIF (decomp_db_bound(2,decomp_type)  ==  BC_CYCLIC) THEN
!       This processor at the bottom of LPG, and has cyclic BCs.
        decomp_db_neighbour(PSouth,decomp_type) =                       &
     &    mype + decomp_db_nproc(decomp_type) -                         &
     &    decomp_db_gridsize(1,decomp_type)
      ELSEIF (decomp_db_bound(2,decomp_type)  ==  BC_OVERPOLE) THEN
! This processor passes data over the pole
         IF ((decomp_db_g_gridpos(1,mype,decomp_type)+1)  <=            &
     &       (decomp_db_gridsize(1,decomp_type)/2)) THEN
           decomp_db_neighbour(PSouth,decomp_type) =                    &
     &       mype+decomp_db_gridsize(1,decomp_type)/2
         ELSE
           decomp_db_neighbour(PSouth,decomp_type) =                    &
     &       mype-decomp_db_gridsize(1,decomp_type)/2
         ENDIF
       ELSE
!       This processor at top of LPG and has static BCs
        decomp_db_neighbour(PSouth,decomp_type) =NoDomain
      ENDIF

! Set Western neighbour

      IF ( decomp_db_g_gridpos(1,mype,decomp_type)  >   0) THEN
!       This processor is not at the left of the LPG
        decomp_db_neighbour(PWest,decomp_type) =                        &
     &    mype - 1
      ELSEIF (decomp_db_bound(1,decomp_type)  ==  BC_CYCLIC) THEN
!       This processor at the left of the LPG, and has cyclic BCs.
        decomp_db_neighbour(PWest,decomp_type) =                        &
     &    mype + decomp_db_gridsize(1,decomp_type) - 1
      ELSE
!       This processor at top of LPG and has static BCs
        decomp_db_neighbour(PWest,decomp_type) =NoDomain
      ENDIF

! Set Eastern neighbour
      IF ( decomp_db_g_gridpos(1,mype,decomp_type)  <                   &
     &    (decomp_db_gridsize(1,decomp_type)-1) ) THEN
!       This processor is not at the right of the LPG
        decomp_db_neighbour(PEast,decomp_type) =                        &
     &    mype + 1
      ELSEIF (decomp_db_bound(1,decomp_type)  ==  BC_CYCLIC) THEN
!       This processor at the left of the LPG, and has cyclic BCs.
        decomp_db_neighbour(PEast,decomp_type) =                        &
     &    mype - decomp_db_gridsize(1,decomp_type) + 1
      ELSE
!       This processor at top of LPG and has static BCs
        decomp_db_neighbour(PEast,decomp_type) =NoDomain
      ENDIF

      RETURN
      END SUBROUTINE SET_NEIGHBOUR
#endif
