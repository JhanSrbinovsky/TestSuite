#if defined(ATMOS)
#if defined(CONTROL) || defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Top-level control routine for vegetation section
!
! Subroutine Interface:
      SUBROUTINE VEG_CTL(                                               &
     &                   row_length, rows, n_rows                       &
     &,                  global_row_length, global_rows                 &
     &,                  DIM_CS1, DIM_CS2                               &
     &,                  halo_i, halo_j, off_x, off_y, me               &
     &,                  n_proc, n_procx, n_procy                       &
     &,                  g_rows, g_row_length, g_datastart              &
     &,                  at_extremity                                   &
     &,                  LAND_PTS,LAND_INDEX,NTILES,CAN_MODEL           &
     &,                  A_STEP,ASTEPS_SINCE_TRIFFID                    &
     &,                  PHENOL_PERIOD,TRIFFID_PERIOD                   &
     &,                  L_PHENOL,L_TRIFFID,L_TRIF_EQ                   &
     &,                  ATIMESTEP,FRAC_DISTURB,SATCON                  &
     &,                  G_LEAF_AC,G_LEAF_PHEN_AC,NPP_AC                &
     &,                  RESP_S_AC,RESP_W_AC                            &
     &,                  CS,FRAC,LAI,CLAY_FRAC,HT                       &
     &,                  CATCH_S,CATCH_T,INFIL_T,Z0_T                   &
     &,                                                                 &
#include "argsts.h"
     &                   STASHwork                                      &
     &                   )


      IMPLICIT NONE
!
! Description:  Calls interim control routine VEG_INTCTL
!
! Current Code Owner:  Richard Betts
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.4    6/10/97   Original code.  Richard Betts
!   4.5    5/8/98    Pass info on grid and halo dimensions into VEG_IC.
!                    Richard Betts
!   4.5   23/11/98   Write diagnostics 19001-19016 to STASH.
!                    Richard Betts
!   5.2   15/11/00   Re-Written for New Dynamics, with an additional
!                    subroutine call to diagnostics_veg       M. Best
! 5.3  17/10/01 Changes required for Single Column Model
!                                             Z. Gardner
!   5.3   26/09/01   Add code required for additional vegetation
!                    diagnostics.    M. Best
!   5.4   28/08/02   Arguement list changed for canopy snow
!                    R. Essery
!   6.2   22/06/06   Fix dependency on diagnostics_veg for SCM.
!                    P.Selwood
!  6.2  01/03/06  increase CS and RESP_S to dimension 4, and pass
!                 down to VEG routines.                     C.D. Jones
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!


      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows           ! number of rows in a v field

      Integer                                                           &
     &  global_row_length                                               &
                            !IN. NUMBER OF points on a global row
     &, global_rows                                                     &
                            !IN. NUMBER OF global rows
     &, me                                                              &
                            !IN. Processor number
     &, halo_i                                                          &
                            !IN. size of large halo in x direction
     &, halo_j                                                          &
                            !IN. size of large halo in y direction
     &, off_x                                                           &
                            !IN. size of small halo in x direction
     &, off_y                                                           &
                            !IN. size of small halo in y direction
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, g_rows(0:n_proc-1)                                              &
     &, g_row_length(0:n_proc-1)                                        &
     &, g_datastart(3,0:n_proc-1)


      INTEGER                                                           &
     & LAND_PTS                                                         &
                                    ! IN Number of land points.
     &,NTILES                                                           &
                                    ! IN Number of land-surface tiles.
     &,CAN_MODEL                                                        &
                                    ! IN Swith for thermal vegetation
     &,A_STEP                                                           &
                                    ! IN Atmospheric timestep number.
     &,ASTEPS_SINCE_TRIFFID                                             &
                                    ! INOUT Number of atmospheric
!                                   !       timesteps since last call
!                                   !       to TRIFFID.
     &,PHENOL_PERIOD                                                    &
                                    ! IN Phenology period (days).
     &,TRIFFID_PERIOD                                                   &
                                    ! IN TRIFFID period (days).
     &,DIM_CS1, DIM_CS2             ! IN soil carbon dimensions

#include "nstypes.h"

      INTEGER                                                           &
     & LAND_INDEX(LAND_PTS)         ! IN I=LAND_INDEX(L) => the Ith
!                                   !    point in P_FIELD is the Lth
!                                   !    land point.
      LOGICAL                                                           &
     & L_PHENOL                                                         &
                                    ! IN .T. for interactive leaf
!                                   !    phenology.
     &,L_TRIFFID                                                        &
                                    ! IN .T. for interactive vegetation.
     &,L_TRIF_EQ                    ! IN .T. for vegetation equilibrium.

      REAL                                                              &
     & ATIMESTEP                                                        &
                                    ! IN Atmospheric timestep (s).
     &,FRAC_DISTURB(LAND_PTS)                                           &
                                    ! IN Fraction of gridbox in which
!                                   !    vegetation is disturbed.
     &,SATCON(LAND_PTS)                                                 &
                                    ! IN Saturated hydraulic
!                                   !    conductivity (kg/m2/s).
     &,G_LEAF_AC(LAND_PTS,NPFT)                                         &
                                    ! INOUT Accumulated leaf turnover
!                                   !       rate.
     &,G_LEAF_PHEN_AC(LAND_PTS,NPFT)                                    &
                                    ! INOUT Accumulated leaf turnover
!                                   !       rate including phenology.
     &,NPP_AC(LAND_PTS,NPFT)                                            &
                                    ! INOUT Accumulated NPP (kg C/m2).
     &,RESP_W_AC(LAND_PTS,NPFT)                                         &
                                    ! INOUT Accumulated wood respiration
!                                   !       (kg C/m2).
     &,RESP_S_AC(LAND_PTS,DIM_CS1)                                      &
                                    ! INOUT Accumulated soil respiration
!                                   !       (kg C/m2).
     &,CS(LAND_PTS,DIM_CS1)                                             &
                               ! INOUT Soil carbon content
!                                   !       (kg C/m2).
     &,FRAC(LAND_PTS,NTYPE)                                             &
                                    ! INOUT Fractions of surface types.
     &,LAI(LAND_PTS,NPFT)                                               &
                                    ! INOUT LAI of plant functional
!                                   !       types.
     &,HT(LAND_PTS,NPFT)                                                &
                                    ! INOUT Height of plant functional
!                                   !       types (m).
     &,ALBSNC(LAND_PTS)                                                 &
                                    ! OUT Snow-covered albedo.
     &,ALBSNF(LAND_PTS)                                                 &
                                    ! OUT Snow-free albedo.
     &,CATCH_S(LAND_PTS)                                                &
                                    ! OUT Canopy snow capacity for NLT
!                                   !     tile (kg/m2).
     &,CATCH_T(LAND_PTS,NTILES)                                         &
                                    ! OUT Canopy capacity for tiles
!                                   !     (kg/m2).
     &,INFIL_T(LAND_PTS,NTILES)                                         &
                                    ! OUT Maximum surface infiltration
!                                   !     rate for tiles (kg/m2/s).
     &,Z0_T(LAND_PTS,NTILES)        ! OUT Roughness length for tiles (m)


#include "csubmodl.h"
#include "typsts.h"

! Diagnostic variables
       Real                                                             &
     &  STASHwork(*)    ! STASH workspace


! Local Variables

      REAL                                                              &
     & C_VEG(LAND_PTS,NPFT)                                             &
                                    ! LOCAL Total carbon content of
!                                   !     the vegetation (kg C/m2).
     &,CV(LAND_PTS)                                                     &
                                    ! LOCAL Gridbox mean vegetation
!                                   !     carbon (kg C/m2).
     &,CS_TOT(DIM_CS2)                                                  &
                                    ! total soil carbon (kg C/m2).
     &,G_LEAF_PHEN(LAND_PTS,NPFT)                                       &
                                    ! LOCAL Mean leaf turnover rate over
!                                   !     phenology period (/360days).
     &,LIT_C(LAND_PTS,NPFT)                                             &
                                    ! LOCAL Carbon Litter
!                                   !     (kg C/m2/360days).
     &,LIT_C_MN(LAND_PTS)                                               &
                                    ! LOCAL Gridbox mean carbon litter
!                                   !     (kg C/m2/360days).
     &,G_LEAF_DAY(LAND_PTS,NPFT)                                        &
                                    ! LOCAL Mean leaf turnover rate for
!                                   !     input to PHENOL (/360days).
     &,LAI_PHEN(LAND_PTS,NPFT)                                          &
                                    ! LOCAL LAI of PFTs after phenology.
     &,G_LEAF_DR_OUT(LAND_PTS,NPFT)                                     &
                                    ! LOCAL Mean leaf turnover rate for
!                                   !       driving TRIFFID (/360days).
     &,NPP_DR_OUT(LAND_PTS,NPFT)                                        &
                                    ! LOCAL Mean NPP for driving TRIFFID
!                                   !     (kg C/m2/360days).
     &,RESP_W_DR_OUT(LAND_PTS,NPFT)                                     &
                                    ! LOCAL Mean wood respiration for
!                                   !       driving TRIFFID
!                                   !       (kg C/m2/360days).
     &,RESP_S_DR_OUT(LAND_PTS,DIM_CS1+1)                                &
                                        ! Mean soil respiration for
                                        ! driving TRIFFID
                                        ! (kg C/m2/360days).
     &,CLAY_FRAC(ROW_LENGTH, ROWS)      ! IN  Clay fraction of soil


! DEPENDS ON: veg_ic
      CALL VEG_IC(LAND_PTS,LAND_INDEX,NTILES,CAN_MODEL                  &
     &,           row_length, rows                                      &
     &,           A_STEP,ASTEPS_SINCE_TRIFFID                           &
     &,           PHENOL_PERIOD,TRIFFID_PERIOD                          &
     &,           L_PHENOL,L_TRIFFID,L_TRIF_EQ                          &
     &,           ATIMESTEP,FRAC_DISTURB,SATCON                         &
     &,           G_LEAF_AC,G_LEAF_PHEN_AC,NPP_AC                       &
     &,           RESP_S_AC,RESP_W_AC                                   &
     &,           CS,FRAC,LAI,CLAY_FRAC,HT                              &
     &,           CATCH_S,CATCH_T,INFIL_T,Z0_T                          &
     &,           C_VEG,CV,LIT_C,LIT_C_MN,G_LEAF_DAY,G_LEAF_PHEN        &
     &,           LAI_PHEN,G_LEAF_DR_OUT,NPP_DR_OUT,RESP_W_DR_OUT       &
     &,           RESP_S_DR_OUT                                         &
     &                )



! ----------------------------------------------------------------------
! Vegetation output diagnostics
! ----------------------------------------------------------------------

#if !defined(SCMA)
! Check that vegetation diagnostics requested this timestep
      if (sf(0,19)) then

! DEPENDS ON: diagnostics_veg
        Call diagnostics_veg(                                           &
     &                     row_length, rows, n_rows                     &
     &,                    global_row_length, global_rows               &
     &,                    DIM_CS1, DIM_CS2                             &
     &,                    halo_i, halo_j, off_x, off_y, me             &
     &,                    n_proc, n_procx, n_procy                     &
     &,                    g_rows, g_row_length, g_datastart            &
     &,                    at_extremity                                 &
     &,                    land_pts                                     &
     &,                    land_index                                   &
     &,                    ntype,npft                                   &
     &,                    c_veg,cv,g_leaf_phen                         &
     &,                    lit_c,lit_c_mn,g_leaf_day                    &
     &,                    lai_phen,g_leaf_dr_out,npp_dr_out            &
     &,                    resp_w_dr_out,resp_s_dr_out,frac_disturb     &
     &,                    frac,lai,ht,cs                               &
     &,                                                                 &
#include "argsts.h"
     &   STASHwork                                                      &
     &     )

      endif
#endif



      RETURN
      END SUBROUTINE VEG_CTL
#endif
#endif
