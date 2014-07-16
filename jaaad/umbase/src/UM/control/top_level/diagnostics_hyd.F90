#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine diagnostics_hyd

      Subroutine diagnostics_hyd(                                       &
     &                       row_length, rows, model_levels             &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      at_extremity                               &
     &,                      land_points, dsm_levels                    &
!  Add inland basin outflow to arguments
     &,                      land_index,inlandout_atm                   &
     &,                      smc, surf_roff, sub_surf_roff              &
     &,                      snow_depth_land, snow_melt, canopy, t_soil &
     &,                      soil_layer_moisture                        &
     &,                      ntiles, snomlt_surf_htf, sthu, sthf        &
     &,                      tot_tfall, snow_tile, melt_tile            &
     &,                      rgrain, land_sea_mask                      &
     &,                      dun_roff, drain, qbase, qbase_zw           &
     &,                      fch4_wetl, fexp, gamtot, ti_mean, ti_sig   &
     &,                      fsat, fwetl, zw, sthzw                     &
     &,                      timestep                                   &
     &,                                                                 &
#include "argsts.h"
     & STASHwork                                                        &
     &     )

! Description:
!   Calculates hydrology-related diagnostics (held in STASH section 8).
!
! Method:
!   Required level lists and logical switches are determined by the
!   calling routine from STASH requests and STASHflags.
!   Intercepted arrays and diagnostic arrays are input from the
!   hydrology routines called previously. Each diagnostic is simply
!   copied into the STASHwork array to be passed on to STASH for
!   output processing.
!
!   Diagnostics currently available (in order calculated):
!   Item  Description
!    208  Soil moisture content
!     23  Snow depth
!    201  Snow melt
!    234  Surface run-off
!    235  Sub-surface run-off
!    225  Soil temperature
!    223  Soil layer moisture
!    204  Surface run-off (accumulated)
!    205  Sub-surface run-off (accumulated)
!    209  Canopy water content
!    202  Land snow melt heat flux
!    231  Land snow melt rate
!    233  Canopy throughfall rate
!    236  Snow amount on tiles
!    237  Snow melt rate on tiles
!    238  Snow grain size on tiles
!    229  Unfrozen soil moisture fraction
!    230  Frozen soil moisture fraction
!    239  Baseflow
!    240  Dunne Runoff
!    241  Baseflow from water table layer
!    242  Wetland methane flux
!    252  Drainage out of "nshyd"th model layer
!    245  Inland basin outflow on atmos grid
!    252  Drainage out of bottom "nshyd"th soil layer (currently layer 4).

! History:
! Version   Date     Comment
! ----     -------     -------
! 5.0  30/11/99 Original version. J-C Thil.
! 5.1  09/12/99 Add error trapping. Replace rmdi_pp by rmdi.
!               Rick Rawlins
! 5.1  27/03/00 Add new diagnostics. J-C Thil.
!  5.1  28/04/00   Access land points via proper 2D indexes instead
!                  of the much objectionable out of bounds first
!                  index & second index set to 1. JC Thil
! 5.2  12/10/00 Add diagnostics 201,209. D Matthews
! 5.3  26/09/01 Add additional hydrology diagnostics.    M. Best
! 5.3  27/07/01 Correct extraction to be over soil not model levels
!               for (8,223),(8,225). R Rawlins
! 5.5  04/02/03 Add additional hydrology diagnostics.    Nic Gedney
!   6.2   21/2/06  Re-route outflow from inland basins to soil moisture
!                  P. Falloon

! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels
     &, number_format                                                   &
                         ! switch controlling number format diagnostics
                         ! are written out in. See PP_WRITE for details.
     &, model_domain     ! indicator as to model type, ie global, lam

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
     &, g_datastart(3,0:n_proc-1)                                       &
     &, land_points                                                     &
                    ! No.of land points being processed, can be 0.
     &, dsm_levels                                                      &
     &, ntiles      ! No. of land-surface tiles ( MOSES II )

      Real                                                              &
     &  lat_rot_NP                                                      &
     &, long_rot_NP


      Real                                                              &
     &  timestep

! Primary Arrays used in all models
      Integer                                                           &
     &  land_index(land_points)      ! set from land_sea_mask

      Real                                                              &
     &  snow_depth_land(land_points)                                    &
                                     !
     &, snow_melt(land_points)                                          &
                                ! snowmelt (kg/m2/s)
     &, canopy (land_points)                                            &
                              ! canopy water content (kg/m2)
     &, smc(land_points)                                                &
                            ! available soil moisture in the
!                                 rootzone (kg/m2).
     &, surf_roff(land_points)                                          &
                                ! surface runoff (kg/m2/s).
     &, sub_surf_roff(land_points)                                      &
                                    ! sub-surface runoff
! Declare inland basin outflow variable
     &, inlandout_atm(land_points)                                      &
                                    !inland basin outflow

     &, t_soil(land_points,dsm_levels)                                  &
     &, soil_layer_moisture(land_points,dsm_levels)                     &
     &, snomlt_surf_htf(row_length, rows)                               &
     &, sthu(land_points,dsm_levels)                                    &
                                       ! Unfrozen soil moisture
!                content of each layer as a fraction of saturation
     &, sthf(land_points,dsm_levels)                                    &
                                       ! Frozen soil moisture content
!                of each layer as a fraction of saturation
     &, tot_tfall(land_points)                                          &
                                       ! total throughfall (kg/m2/s)
     &, snow_tile(land_points,ntiles)                                   &
                                       ! Lying snow on tiles (kg/m2)
     &, melt_tile(land_points,ntiles)                                   &
                                       ! Snowmelt on tiles (kg/m2/s)
     &, rgrain(land_points,ntiles)                                      &
                                       ! Snow grain size (microns)

! Additional variables required for large-scale hydrology:
     &, qbase(land_points)                                              &
                                    ! Baseflow (kg/m2/s).
     &, dun_roff(land_points)                                           &
                                    ! Dunne runoff (kg/m2/s).
     &, qbase_zw(land_points)                                           &
                                    ! Baseflow from Zw layer (kg/m2/s).
     &, drain(land_points)                                              &
                                    ! Drainage out of "nshyd" later
!                                   ! (kg/m2/s).
     &, fch4_wetl(land_points)                                          &
                                    ! Wetland methane flux (kg C/m2/s).
     &, TI_MEAN(LAND_POINTS)                                            &
                                    ! Mean topographic index.
     &, TI_SIG(LAND_POINTS)                                             &
                                    ! Std. dev. of topographic index.
     &, FEXP(LAND_POINTS)                                               &
                                    ! Decay factor in Sat. Conductivity
!                                   !   in water table layer.
     &, GAMTOT(LAND_POINTS)                                             &
                                    ! Integrated complete Gamma
!                                   !   function.
     &, FSAT(LAND_POINTS)                                               &
                                    ! Surface saturation fraction.
     &, FWETL(LAND_POINTS)                                              &
                                    ! Wetland fraction.
     &, ZW(LAND_POINTS)                                                 &
                                    ! Water table depth (m).
     &, STHZW(LAND_POINTS)          ! Saturation fraction in water
!                                   !   table layer.


      Logical                                                           &
     &  land_sea_mask(row_length, rows)

      Integer                                                           &
     & PSLEVEL                                                          &
                     !  loop counter for pseudolevels
     &,PSLEVEL_OUT                                                      &
                     !  index for pseudolevels sent to STASH
     &,LEVEL                                                            &
                     !  loop counter for levels
     &,LEVEL_OUT     !  index for levels sent to STASH

      Logical                                                           &
     & PLLTILE(NTILES)                                                  &
                          ! pseudolevel list for surface types
     &,LIST(DSM_LEVELS)   ! level list for soil levels


#include "c_mdi.h"

! Local variables

      Integer                                                           &
     &  i, j, k, l                                                      &
     &,    icode                ! Return code  =0 Normal exit  >1 Error

      Character*80  cmessage
      Character(*) RoutineName
      Parameter ( RoutineName='diagnostics_hyd')

      Integer                                                           &
     &  im_index        ! internal model index

      Real                                                              &
     &  interp_data(row_length,rows)                                    &
     &, interp_data_3(row_length,rows,dsm_levels)

#include "csubmodl.h"
#include "typsts.h"

! Diagnostic variables
       Real                                                             &
     &  STASHwork(*)    ! STASH workspace

      External                                                          &
     &  copydiag, copydiag_3d                                           &
     &  ,Ereport

! ----------------------------------------------------------------------
! Section 1.  Initialisation.
! ----------------------------------------------------------------------


      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! ----------------------------------------------------------------------
! Soil moisture content
! ----------------------------------------------------------------------
! Item 8 208  smc

      If (sf(208,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = smc(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(208,8,im_index)),interp_data,      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,208,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 208)"
            goto 9999
         End if

      End if


! ----------------------------------------------------------------------
! Snow depth
! ----------------------------------------------------------------------
! Item 8 023 snow_depth_land

      If (sf(023,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = snow_depth_land(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(023,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,023,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 023)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Snow melt
! ----------------------------------------------------------------------
! Item 201

      If (sf(201,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = snow_melt(l) * timestep
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(201,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,201,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 201)"
            goto 9999
         End if

      End if


! ----------------------------------------------------------------------
! Surface run-off.
! ----------------------------------------------------------------------
! Item 234  surf_roff

      If (sf(234,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = surf_roff(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(234,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,234,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 234)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Sub-surface run-off.
! ----------------------------------------------------------------------
! Item 235  sub_surf_roff

      If (sf(235,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = sub_surf_roff(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(235,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,235,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 235)"
            goto 9999
         End if

      End if


! ----------------------------------------------------------------------
!  Soil temperature
! ----------------------------------------------------------------------
! Item 8 225 t_soil

      If (sf(225,8)) Then

         Do k = 1, dsm_levels
            Do j= 1, rows
               Do i = 1, row_length
                  interp_data_3(i,j,k) = rmdi
               End Do
            End Do

            Do l = 1, land_points
               j=(land_index(l)-1)/row_length + 1
               i=land_index(l) - (j-1)*row_length
               interp_data_3(i,j,k) = t_soil(l,k)
            End Do

         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(225,8,im_index)),interp_data_3,  &
     &        row_length,rows,dsm_levels,0,0,0,0,                       &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,225,8,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,8,225,                                           &
     &        icode,cmessage)



         If (icode  >   0) then
            cmessage="Error in copydiag_3d( item 225)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Soil layer moisture
! ----------------------------------------------------------------------
! Item 8 223 soil_layer_moisture

      If (sf(223,8)) Then

         Do k = 1, dsm_levels
            Do j= 1, rows
               Do i = 1, row_length
                  interp_data_3(i,j,k) = rmdi
               End Do
            End Do
             Do l = 1, land_points
               j=(land_index(l)-1)/row_length + 1
               i=land_index(l) - (j-1)*row_length
               interp_data_3(i,j,k) =                                   &
     &              soil_layer_moisture(l,k)
             End Do
         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(223,8,im_index)),interp_data_3,  &
     &        row_length,rows,dsm_levels,0,0,0,0,                       &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,223,8,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,8,223,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 223)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Surface run-off (accumulated)
! ----------------------------------------------------------------------
! Item 204  surf_roff

      If (sf(204,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do i = 1, land_points
            interp_data(land_index(i),1) = surf_roff(i) * timestep
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(204,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,204,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 204 )"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Sub-surface run-off (accumulated)
! ----------------------------------------------------------------------
! Item 205  sub_surf_roff

      If (sf(205,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do i = 1, land_points
            interp_data(land_index(i),1) = sub_surf_roff(i)             &
     &     *  timestep
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(205,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,205,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 205)"
            goto 9999
         endif
      endif

! ----------------------------------------------------------------------
! Canopy water content
! ----------------------------------------------------------------------
! Item 209

      If (sf(209,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do i = 1, land_points
            interp_data(land_index(i),1) = canopy(i)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(209,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,209,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 209)"
            goto 9999
         endif
      endif


! ----------------------------------------------------------------------
! Land snow melt heat flux (W/m2)
! ----------------------------------------------------------------------
! Item 202 SNOMLT_SURF_HTF

      IF (SF(202,8)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(202,8,im_index)),                              &
     &       snomlt_surf_htf,                                           &
     &       land_sea_mask,row_length*rows,land_points)
      END IF


! ----------------------------------------------------------------------
! Land snow melt heat rate (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 231 snow_melt

      IF (SF(231,8)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(231,8,im_index)),                              &
     &       snow_melt,                                                 &
     &       land_sea_mask,row_length*rows,land_points)
      END IF


! ----------------------------------------------------------------------
! Canopy throughfall rate (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 233 tot_tfall

      IF (SF(233,8)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(233,8,im_index)),                              &
     &       tot_tfall,                                                 &
     &       land_sea_mask,row_length*rows,land_points)
      END IF


! ----------------------------------------------------------------------
! Snow amount on tiles (Kg/m2)
! ----------------------------------------------------------------------
! Item 236 snow_tile

      IF (SF(236,8)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,236,8,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "daghyd  : error in set_pseudo_list(item 236 = snow_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(236,8,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),snow_tile(1,PSLEVEL_OUT),            &
     &           land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF



! ----------------------------------------------------------------------
! Snow melt rate on tiles (Kg/m2)
! ----------------------------------------------------------------------
! Item 237 melt_tile

      IF (SF(237,8)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,237,8,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "daghyd  : error in set_pseudo_list(item 237 = melt_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(237,8,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),melt_tile(1,PSLEVEL_OUT),            &
     &           land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF



! ----------------------------------------------------------------------
! Snow grain size on tiles
! ----------------------------------------------------------------------
! Item 238 rgrain

      IF (SF(238,8)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,238,8,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "daghyd  : error in set_pseudo_list(item 238 = rgrain)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(238,8,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),rgrain(1,PSLEVEL_OUT),               &
     &           land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF



! ----------------------------------------------------------------------
! Unfrozen soil moisture fraction
! ----------------------------------------------------------------------
! Item 229 sthu

      IF (SF(229,8)) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(DSM_LEVELS,LEN_STLIST,                     &
     &       STLIST(1,STINDEX(1,229,8,im_index)),                       &
     &       LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "daghyd  : error in set_levels_list(item 229 = sthu)"
            goto 9999
        END IF
        LEVEL_OUT=0
        DO LEVEL=1,DSM_LEVELS
          IF(LIST(LEVEL)) THEN
            LEVEL_OUT=LEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(229,8,im_index)+(LEVEL_OUT-1)              &
     &           *row_length*rows),sthu(1,LEVEL_OUT),                   &
     &           land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF


! ----------------------------------------------------------------------
! Frozen soil moisture fraction
! ----------------------------------------------------------------------
! Item 230 sthf

      IF (SF(230,8)) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(DSM_LEVELS,LEN_STLIST,                     &
     &       STLIST(1,STINDEX(1,230,8,im_index)),                       &
     &       LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "daghyd  : error in set_levels_list(item 230 = sthf)"
            goto 9999
        END IF
        LEVEL_OUT=0
        DO LEVEL=1,DSM_LEVELS
          IF(LIST(LEVEL)) THEN
            LEVEL_OUT=LEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(230,8,im_index)+(LEVEL_OUT-1)              &
     &           *row_length*rows),sthf(1,LEVEL_OUT),                   &
     &           land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF


! ----------------------------------------------------------------------
! Baseflow
! ----------------------------------------------------------------------
! Item 239  baseflow

      If (sf(239,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = qbase(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(239,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,239,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 239)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Dunne Runoff
! ----------------------------------------------------------------------
! Item 240 Dunne Runoff

      If (sf(240,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = dun_roff(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(240,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,240,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 240)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Baseflow from water table layer
! ----------------------------------------------------------------------
! Item 241  baseflow from water table layer

      If (sf(241,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = qbase_zw(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(241,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,241,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 241)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Wetland methane flux
! ----------------------------------------------------------------------
! Item 242  Wetland methane flux

      If (sf(242,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fch4_wetl(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(242,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,242,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 242)"
            goto 9999
         End if

      End if
      If (sf(243,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = ti_mean(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(243,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,243,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 243)"
            goto 9999
         End if
      End if
      If (sf(244,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = ti_sig(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(244,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,244,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 244)"
            goto 9999
         End if
      End if
      If (sf(251,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fexp(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(251,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,251,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 251)"
            goto 9999
         End if
      End if
      If (sf(246,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = gamtot(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(246,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,246,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 246)"
            goto 9999
         End if
      End if
      If (sf(247,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fsat(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(247,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,247,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 247)"
            goto 9999
         End if
      End if
      If (sf(248,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fwetl(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(248,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,248,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 248)"
            goto 9999
         End if
      End if
      If (sf(249,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = zw(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(249,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,249,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 249)"
            goto 9999
         End if
      End if
      If (sf(250,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = sthzw(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(250,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,250,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 250)"
            goto 9999
         End if
      End if
! ----------------------------------------------------------------------
! Drainage out of  "nshyd"th model layer
! ----------------------------------------------------------------------
! Item 252  drainage

      If (sf(252,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = drain(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(252,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,252,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 252)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Output inland basin outflow on atmosphere grid

! --------------------------------------------------------------------
! Inland basin outflow (atmos grid)
! --------------------------------------------------------------------
! Item 245  Inland basin outflow

      If (sf(245,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = inlandout_atm(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(245,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,245,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 245)"
            goto 9999
         End if

      End if
! -------------------------------------------
 9999 continue
      If(icode /= 0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,Cmessage)
      Endif

      Return
      END SUBROUTINE diagnostics_hyd

#endif
