#if defined(A19_0A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine Interface:
      SUBROUTINE VEG_IC(LAND_PTS,LAND_INDEX,NTILES,CAN_MODEL            &
     &,                 row_length, rows                                &
     &,                 A_STEP,ASTEPS_SINCE_TRIFFID                     &
     &,                 PHENOL_PERIOD,TRIFFID_PERIOD                    &
     &,                 L_PHENOL,L_TRIFFID,L_TRIF_EQ                    &
     &,                 ATIMESTEP,FRAC_DISTURB,SATCON                   &
     &,                 G_LEAF_AC,G_LEAF_PHEN_AC,NPP_AC                 &
     &,                 RESP_S_AC,RESP_W_AC                             &
     &,                 CS,FRAC,LAI,CLAY_FRAC,HT                        &
     &,                 CATCH_S,CATCH_T,INFIL_T,Z0_T                    &
     &,                 C_VEG,CV,LIT_C,LIT_C_MN,G_LEAF_DAY,G_LEAF_PHEN  &
     &,                 LAI_PHEN,G_LEAF_DR_OUT,NPP_DR_OUT,RESP_W_DR_OUT &
     &,                 RESP_S_DR_OUT                                   &
     &                  )


      IMPLICIT NONE
!
! Description:  Version 0A of vegetation section: Dummy routine to
!               ensure compilation
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2   15/11/00   New Deck         M. Best
!   5.4   28/08/02   Arguement list changed for canopy snow
!                    R. Essery
!   6.0   11/09/03   Replace ABORT with call to EREPORT.   P.Dando
!   6.2   01/03/06   Arguement list changed for CLAY_FRAC
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.

      INTEGER                                                           &
     & LAND_PTS                                                         &
                                    ! Number of land points.
     &,NTILES                                                           &
                                    ! Number of land-surface tiles.
     &,CAN_MODEL                                                        &
                                    ! Swith for thermal vegetation
     &,row_length                                                       &
                        ! number of points on a row
     &,rows                                                             &
                        ! number of rows in a theta field
     &,A_STEP                                                           &
                                    ! Atmospheric timestep number.
     &,ASTEPS_SINCE_TRIFFID                                             &
                                    ! Number of atmospheric
!                                   !       timesteps since last call
!                                   !       to TRIFFID.
     &,PHENOL_PERIOD                                                    &
                                    ! Phenology period (days).
     &,TRIFFID_PERIOD               ! TRIFFID period (days).

#include "nstypes.h"

      INTEGER                                                           &
     & LAND_INDEX(LAND_PTS)         ! I=LAND_INDEX(L) => the Ith
!                                   !    point in P_FIELD is the Lth
!                                   !    land point.
      LOGICAL                                                           &
     & L_PHENOL                                                         &
                                    ! .T. for interactive leaf
!                                   !    phenology.
     &,L_TRIFFID                                                        &
                                    ! .T. for interactive vegetation.
     &,L_TRIF_EQ                    ! .T. for vegetation equilibrium.

      REAL                                                              &
     & ATIMESTEP                                                        &
                                    ! Atmospheric timestep (s).
     &,FRAC_DISTURB(LAND_PTS)                                           &
                                    ! Fraction of gridbox in which
!                                   !    vegetation is disturbed.
     &,SATCON(LAND_PTS)                                                 &
                                    ! Saturated hydraulic
!                                   !    conductivity (kg/m2/s).
     &,G_LEAF_AC(LAND_PTS,NPFT)                                         &
                                    ! Accumulated leaf turnover
!                                   !       rate.
     &,G_LEAF_PHEN_AC(LAND_PTS,NPFT)                                    &
                                    ! Accumulated leaf turnover
!                                   !       rate including phenology.
     &,NPP_AC(LAND_PTS,NPFT)                                            &
                                    ! Accumulated NPP (kg C/m2).
     &,RESP_W_AC(LAND_PTS,NPFT)                                         &
                                    ! Accumulated wood respiration
!                                   !       (kg C/m2).
     &,RESP_S_AC(LAND_PTS)                                              &
                                    ! Accumulated soil respiration
!                                   !       (kg C/m2).
     &,CLAY_FRAC(ROW_LENGTH, ROWS)                                      &
                                        ! IN  Clay fraction of soil
     &,CS(LAND_PTS)                                                     &
                                    ! Soil carbon content
!                                   !       (kg C/m2).
     &,FRAC(LAND_PTS,NTYPE)                                             &
                                    ! Fractions of surface types.
     &,LAI(LAND_PTS,NPFT)                                               &
                                    ! LAI of plant functional
!                                   !       types.
     &,HT(LAND_PTS,NPFT)                                                &
                                    ! Height of plant functional
!                                   !       types (m).
     &,ALBSNC(LAND_PTS)                                                 &
                                    ! Snow-covered albedo.
     &,ALBSNF(LAND_PTS)                                                 &
                                    ! Snow-free albedo.
     &,CATCH_S(LAND_PTS)                                                &
                                    ! Canopy snow capacity for NLT
     &,CATCH_T(LAND_PTS,NTILES)                                         &
                                    ! Canopy capacity for tiles
!                                   !     (kg/m2).
     &,INFIL_T(LAND_PTS,NTILES)                                         &
                                    ! Maximum surface infiltration
!                                   !     rate for tiles (kg/m2/s).
     &,G_LEAF_DAY(LAND_PTS,NPFT)                                        &
                                    ! Mean leaf turnover rate for
!                                   !     input to PHENOL (/360days).
     &,G_LEAF_PHEN(LAND_PTS,NPFT)                                       &
                                    ! Mean leaf turnover rate over
!                                   !     phenology period (/360days).
     &,G_LEAF_DR_OUT(LAND_PTS,NPFT)                                     &
                                    ! Mean leaf turnover rate for
!                                   !       driving TRIFFID (/360days).
     &,LAI_PHEN(LAND_PTS,NPFT)                                          &
                                    ! LAI of PFTs after phenology.
     &,NPP_DR_OUT(LAND_PTS,NPFT)                                        &
                                    ! Mean NPP for driving TRIFFID
!                                   !     (kg C/m2/360days).
     &,RESP_W_DR_OUT(LAND_PTS,NPFT)                                     &
                                    ! Mean wood respiration for
!                                   !       driving TRIFFID
!                                   !       (kg C/m2/360days).
     &,RESP_S_DR_OUT(LAND_PTS)                                          &
                                    ! Mean soil respiration for
!                                   !     driving TRIFFID
!                                   !     (kg C/m2/360days).
     &,Z0_T(LAND_PTS,NTILES)                                            &
                                    ! Roughness length for tiles (m)
     &,C_VEG(LAND_PTS,NPFT)                                             &
                                    ! Total carbon content of
!                                   !     the vegetation (kg C/m2).
     &,CV(LAND_PTS)                                                     &
                                    ! Gridbox mean vegetation
!                                   !     carbon (kg C/m2).
     &,LIT_C(LAND_PTS,NPFT)                                             &
                                    ! Carbon Litter
!                                   !     (kg C/m2/360days).
     &,LIT_C_MN(LAND_PTS)           ! Gridbox mean carbon litter
!                                   !     (kg C/m2/360days).

! Local variables:
      INTEGER             :: ErrorStatus = 0     ! Error return code
      CHARACTER (Len=256) :: CMessage = ' '      ! Error message if
                                                 ! return code >0
! Local parameters:
      CHARACTER (Len=*), Parameter ::  RoutineName = 'VEG_IC'


      WRITE (6,*) '**ERROR**: VEG_IC has been called but is'
      WRITE (6,*) 'unavailable.  Either set L_PHENOL and'
      WRITE (6,*) 'L_TRIFFID to .FALSE. or'
      WRITE (6,*) 'select section A19_1A or A19_2A.'

      CMessage='Routine unavailable - see output for details'
      ErrorStatus = 1

! DEPENDS ON: ereport
      CALL EREPORT(RoutineName,ErrorStatus,CMessage)

      RETURN
      END SUBROUTINE VEG_IC
#endif
