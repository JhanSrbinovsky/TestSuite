#if defined(A70_1B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine diagnostics_sw

      Subroutine diagnostics_sw(                                        &
     &                       row_length, rows, model_levels             &
     &,                      wet_model_levels, cloud_levels             &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      at_extremity                               &
     &,                      timestep                                   &
     &,                      T_n, T_inc                                 &
     &,                      q_n, qcl_n, cf_n, cfl_n                    &
     &,                      T_latest, q_latest, qcl_latest             &
     &,                      cf_latest, cfl_latest                      &
     &,                      surfsw, itoasw, solar_out_toa              &
     &,                      solar_out_clear, surface_down_flux         &
     &,                      surf_down_clr, surf_up_clr                 &
     &,                      SWsea, flux_below_690nm_surf               &
     &,                      photosynth_act_rad, dirpar_flux            &
     &,                      fl_solid_below_690nm_surf                  &
     &,                      fl_sea_below_690nm_surf                    &
     &,                      orog_corr, sol_bearing                     &
     &,                      f_orog                                     &
     &,                      slope_aspect, slope_angle                  &
     &,                      sw_net_land,sw_net_sice                    &
     &,                      T_incr_diagnostic                          &
     &,                      clear_hr                                   &
     &,                      flux_direct, flux_diffuse                  &
     &,                      net_flux_trop, up_flux_trop                &
     &,                      re_strat, wgt_strat, lwp_strat             &
     &,                      re_conv, wgt_conv                          &
     &,                      ntot_diag, strat_lwc_diag                  &
     &,                      so4_ccn_diag, cond_samp_wgt                &
     &,                      weighted_re, sum_weight_re                 &
     &,                      weighted_warm_re, sum_weight_warm_re       &
     &,                      Nc_diag, Nc_weight                         &
     &,                      sea_salt_film, sea_salt_jet                &
     &,                      salt_dim1, salt_dim2, salt_dim3            &
     &,                      cloud_extinction                           &
     &,                      cloud_weight_extinction                    &
     &,                      ls_cloud_extinction                        &
     &,                      ls_cloud_weight_extinction                 &
     &,                      cnv_cloud_extinction                       &
     &,                      cnv_cloud_weight_extinction                &
     &,                                                                 &
#include "argsts.h"
     & STASHwork                                                        &
     &     )

! Purpose:
!          Calculates diagnostics and outputs them.
!
! Method:
!
! History:
! Version   Date     Comment
! ---      -----     -------
! 5.0  30/11/99 Original version in UM. J-C Thil.
! 5.1  09/12/99 Add error trapping. Rick Rawlins
! ----     -------     -------
! 5.1  27/03/00 Add new diagnostics. J-C Thil.
! 5.1  28/02/00 Provide explicit model increments as STASH output
!               diagnostics. R Rawlins
! 5.2  14/11/00 Provide code for reactivated diagnostics.
!                                 (J. M. Edwards)
! 5.2  15/11/00 Add sea-salt aerosol diagnostics. A. Jones
! 5.3  25/04/01  Add diagnostics for coastal tiling.
!                                                   N. Gedney
! 5.4  03/08/01 Add PC2 cloud scheme diagnostics. D.R. Wilson
! 5.4  29/05/02 Add column-integrated cloud droplet and warm-cloud-only
!               satellite-view rE diagnostics.
!                                                 A. Jones
! 5.5  09/01/02 Add Extinction diagnostics. A.B.Keen/K.D.Williams
! 5.5  15/05/03 Correct cmessage for diagnostic 248.  A. Jones
! 6.2  02/03/06 Added diagnostics for total and direct component
!               of surface PAR flux.  M.G. Sanderson
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels
     &, cloud_levels                                                    &
                         ! number of cloudy levels
     &, wet_model_levels                                                &
                         ! number of model levels where moisture
                         ! variables are held
     &, number_format                                                   &
                         ! switch controlling number format diagnostics
                         ! are written out in. See PP_WRITE for details.
     &, model_domain                                                    &
                         ! indicator as to model type, ie global, lam
     &, salt_dim1                                                       &
                         !
     &, salt_dim2                                                       &
                         ! Dimensions for sea-salt aerosol diagnostics.
     &, salt_dim3        !

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
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       &
     &, g_datastart (3,0:n_proc-1)

      Real                                                              &
     &  timestep

! Primary Arrays used in all models

      Real                                                              &
     &  p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, T_n(row_length, rows, model_levels)                             &
     &, T_inc(row_length, rows, model_levels)                           &
     &, q_n(row_length, rows, wet_model_levels)                         &
     &, qcl_n(row_length, rows, wet_model_levels)                       &
     &, cf_n(row_length, rows, wet_model_levels)                        &
     &, cfl_n(row_length, rows, wet_model_levels)                       &
     &, T_latest(row_length, rows, model_levels)                        &
     &, q_latest(row_length, rows, wet_model_levels)                    &
     &, qcl_latest(row_length, rows, wet_model_levels)                  &
     &, cf_latest(row_length, rows, wet_model_levels)                   &
     &, cfl_latest(row_length, rows, wet_model_levels)                  &

     &, exner_rho_levels(1-off_x:row_length+off_x,                      &
     &                     1-off_y:rows+off_y, model_levels)            &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                       1-off_y:rows+off_y, model_levels)          &
     &, p_theta_levels(1-off_x:row_length+off_x,                        &
     &                   1-off_y:rows+off_y, model_levels)              &
     &, q(row_length,rows, wet_model_levels)                            &
     &, qcl(row_length, rows, wet_model_levels)                         &
     &, qcf(row_length, rows, wet_model_levels)

      Real                                                              &
     &  itoasw (row_length, rows)                                       &
     &, surfsw (row_length, rows)                                       &
     &, solar_out_toa (row_length, rows)                                &
     &, solar_out_clear (row_length, rows)                              &
     &, surface_down_flux (row_length, rows)                            &
     &, surf_down_clr (row_length, rows)                                &
     &, surf_up_clr (row_length, rows)                                  &
     &, SWsea(row_length, rows)                                         &
                                  ! Net short-wave absorbed by planet
     &, flux_below_690nm_surf(row_length, rows)                         &
     &, photosynth_act_rad(row_length, rows)                            &
     &, dirpar_flux(row_length, rows)                                   &
     &, fl_solid_below_690nm_surf(row_length, rows)                     &
     &, fl_sea_below_690nm_surf(row_length, rows)                       &
     &, sw_net_land(row_length, rows)                                   &
                                       !SW net local flux over land
     &, sw_net_sice(row_length, rows)                                   &
                                       !SW net local flux over sea-ice
     &, T_incr_diagnostic(row_length,rows,model_levels)                 &
     &, clear_hr(row_length, rows, model_levels)                        &
     &, flux_direct(row_length, rows, model_levels+1)                   &
     &, flux_diffuse(row_length, rows, model_levels+1)                  &
     &, net_flux_trop(row_length, rows)                                 &
     &, up_flux_trop(row_length, rows)                                  &
     &, re_strat(row_length, rows, cloud_levels)                        &
     &, wgt_strat(row_length, rows, cloud_levels)                       &
     &, lwp_strat(row_length, rows, cloud_levels)                       &
     &, re_conv(row_length, rows, cloud_levels)                         &
     &, wgt_conv(row_length, rows, cloud_levels)                        &
     &, ntot_diag(row_length, rows, cloud_levels)                       &
     &, strat_lwc_diag(row_length, rows, cloud_levels)                  &
     &, so4_ccn_diag(row_length, rows, cloud_levels)                    &
     &, cond_samp_wgt(row_length, rows, cloud_levels)                   &
     &, weighted_re(row_length, rows)                                   &
     &, sum_weight_re(row_length, rows)                                 &
     &, weighted_warm_re(row_length, rows)                              &
     &, sum_weight_warm_re(row_length, rows)                            &
     &, Nc_diag(row_length, rows)                                       &
     &, Nc_weight(row_length, rows)                                     &
     &, sea_salt_film(salt_dim1, salt_dim2, salt_dim3)                  &
     &, sea_salt_jet(salt_dim1, salt_dim2, salt_dim3)                   &
     &, cloud_extinction(row_length, rows, cloud_levels)                &
     &, cloud_weight_extinction(row_length, rows, cloud_levels)         &
     &, ls_cloud_extinction(row_length, rows, cloud_levels)             &
     &, ls_cloud_weight_extinction(row_length, rows, cloud_levels)      &
     &, cnv_cloud_extinction(row_length, rows, cloud_levels)            &
     &, cnv_cloud_weight_extinction(row_length, rows, cloud_levels)

! Orography variables

      Real orog_corr(row_length, rows)   ! Orography correction factor
      Real sol_bearing(row_length,rows)  ! Local solar bearing
      Real f_orog(row_length,rows)       ! Extra SW surf flux
      Real slope_aspect(row_length,rows) ! Gridbox mean slope aspect
      Real slope_angle(row_length,rows)  ! Gridbox mean slope angle

#include "csubmodl.h"
#include "typsts.h"

! Diagnostic variables
      Real                                                              &
     & STASHwork(*)    ! STASH workspace

! Local array & variables
      Real                                                              &
     &  T_plus_T_inc(row_length, rows, model_levels)                    &
     &, heating_rate(row_length,rows,model_levels)                      &
     &, work_3d(row_length,rows,model_levels)

      Integer                                                           &
     &  i, j, k                                                         &
     &,    icode                ! Return code  =0 Normal exit  >1 Error

      Character*80  cmessage
      Character(*) RoutineName
      Parameter ( RoutineName='diagnostics_sw')

      Integer                                                           &
     &  im_index                                                        &
                        ! internal model index
     & ,item                                                            &
                        ! STASH item
     & ,sect            ! STASH section
      Parameter( sect = 1 ) ! for sw radiation

! External routines
      External                                                          &
     &  copydiag, copydiag_3d                                           &
     &  ,Ereport

      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! ----------------------------------------------------------------------
! Section 1.  Diagnostic Calculation and output.
! ----------------------------------------------------------------------

      If (sf(004,1)) Then

!L   T+T_inc

         Do k = 1, model_levels
            Do j = 1, rows
               Do i = 1, row_length
                  T_plus_T_inc(i,j,k) = T_n(i,j,k) + T_inc(i,j,k)
               End Do
            End Do
         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(004,1,im_index)),                &
     &        T_plus_T_inc,                                             &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,004,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,004,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage="Error in copydiag_3d( item 004)"
            goto 9999
         End if

      End if

! increment diagnostics= modified - previous

      item = 161  ! temperature increment
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        T_incr_diagnostic,                                        &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 161)"//cmessage
         End if

      Endif  !  sf(item,sect)
      item = 181  ! temperature increment including condensation
      If (icode <= 0 .and. sf(item,sect)) Then
!
      Do k = 1,model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = T_latest(i,j,k) - T_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 181)"//cmessage
         End if

      Endif  !  sf(item,sect)

      item = 182  ! Vapour increment
      If (icode <= 0 .and. sf(item,sect)) Then
!
      Do k = 1,wet_model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = q_latest(i,j,k) - q_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 182)"//cmessage
         End if

      Endif  !  sf(item,sect)

      item = 183  ! liquid water content increment
      If (icode <= 0 .and. sf(item,sect)) Then
!
      Do k = 1,wet_model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = qcl_latest(i,j,k) - qcl_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 183)"//cmessage
         End if

      Endif  !  sf(item,sect)

      item = 192  ! total cloud fraction increment
      If (icode <= 0 .and. sf(item,sect)) Then
!
      Do k = 1,wet_model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = cf_latest(i,j,k) - cf_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 192)"//cmessage
         End if

      Endif  !  sf(item,sect)
!
      item = 193  ! liquid cloud fraction increment
      If (icode <= 0 .and. sf(item,sect)) Then
!
      Do k = 1,wet_model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = cfl_latest(i,j,k) - cfl_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 193)"//cmessage
         End if

      Endif  !  sf(item,sect)


      If (sf(201,1)) Then

!L   surfsw

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(201,1,im_index)),surfsw,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,201,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 201)"
            goto 9999
         End if

      End if



      If (sf(207,1)) Then

!L   itoasw

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(207,1,im_index)),itoasw,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,207,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 207)"
            goto 9999
         End if

      End if



      If (sf(208,1)) Then

!L   solar_out_toa

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(208,1,im_index)),solar_out_toa,    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,208,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 208)"
            goto 9999
         End if

      End if


      If (sf(209,1)) Then

!L   solar_out_clear

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(209,1,im_index)),solar_out_clear,  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,209,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 209)"
            goto 9999
         End if

      End if



      If (sf(210,1)) Then

!L   surf_down_clr

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(210,1,im_index)),surf_down_clr,    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,210,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 210)"
            goto 9999
         End if

      End if



      If (sf(211,1)) Then

!L   surf_up_clr

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(211,1,im_index)),surf_up_clr,      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,211,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 211)"
            goto 9999
         End if

      End if

      If (sf(235,1)) Then

!L   surface_down_flux

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(235,1,im_index)),surface_down_flux,&
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,235,                                           &
!#include <argppx/argppx.h>
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 235)"
            goto 9999
         End if

      End if


      If (sf(203,1)) Then

!L   SWsea

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(203,1,im_index)),SWsea,            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,203,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 203)"
            goto 9999
         End if

      End if


      If (sf(204,1)) Then

!L flux_below_690nm_surf

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(204,1,im_index)),                   &
     &        flux_below_690nm_surf,                                    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,204,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 204)"
            goto 9999
         End if

      End if
!
      If (sf(257,1)) Then

!L sw_net_land

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(257,1,im_index)),                   &
     &        sw_net_land,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,257,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 257)"
            goto 9999
         End if

      End if
!
      If (sf(258,1)) Then

!L sw_net_sice

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(258,1,im_index)),                   &
     &        sw_net_sice,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,258,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 258)"
            goto 9999
         End if

      End if
!
      If (sf(259,1)) Then

!L fl_solid_below_690nm_surf

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(259,1,im_index)),                   &
     &        fl_solid_below_690nm_surf,                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,259,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 259)"
            goto 9999
         End if

      End if
!
      If (sf(260,1)) Then

!L fl_sea_below_690nm_surf

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(260,1,im_index)),                   &
     &        fl_sea_below_690nm_surf,                                  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,260,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 260)"
            goto 9999
         End if

      End if
!


      If (sf(221,1)) Then

!L re_strat

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(221,1,im_index)),                &
     &        re_strat,                                                 &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,221,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,221,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 221)"
            goto 9999
         End if

      End if



      If (sf(223,1)) Then

!L wgt_strat

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(223,1,im_index)),                &
     &        wgt_strat,                                                &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,223,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,223,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 223)"
            goto 9999
         End if

      End if



      If (sf(224,1)) Then

!L lwp_strat

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(224,1,im_index)),                &
     &        lwp_strat,                                                &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,224,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,224,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 224)"
            goto 9999
         End if

      End if



      If (sf(225,1)) Then

!L re_conv

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(225,1,im_index)),                &
     &        re_conv,                                                  &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,225,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,225,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 225)"
            goto 9999
         End if

      End if



      If (sf(226,1)) Then

!L wgt_conv

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(226,1,im_index)),                &
     &        wgt_conv,                                                 &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,226,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,226,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 226)"
            goto 9999
         End if

      End if

      If (sf(230,1)) Then
      
! flux_direct

         Call copydiag_3d(STASHwork(si(230,1,im_index)),                &
     &        flux_direct,                                              &
     &        row_length,rows,model_levels+1,0,0,0,0, at_extremity,     &
     &        stlist(1,stindex(1,230,1,im_index)),len_stlist,           &  
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,230,                                           &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 230)"
            goto 9999
         End if

      End if

      If (sf(231,1)) Then

! flux_diffuse

         Call copydiag_3d(STASHwork(si(231,1,im_index)),                &
     &        flux_diffuse,                                             & 
     &        row_length,rows,model_levels+1,0,0,0,0, at_extremity,     &
     &        stlist(1,stindex(1,231,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,231,                                           &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 231)"
            goto 9999
         End if

      End if

!
! SW heating =
!            sw radiation temperature increment per timestep / timestep

      If (icode <= 0 .and. sf(232,1)) Then

        Do k=1,model_levels
         Do j=1,rows
          Do i=1,row_length
            heating_rate(i,j,k) =  T_incr_diagnostic(i,j,k) /           &
     &                                  timestep
          Enddo ! i
         Enddo ! j
        Enddo ! k

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(232,1,im_index)),                &
     &        heating_rate,                                             &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,232,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,232,                                           &
     &        icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 232)"//cmessage
        End if

      Endif

      If (sf(233,1)) Then

!L clear_hr

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(233,1,im_index)),                &
     &        clear_hr,                                                 &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,233,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,233,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 233)"
            goto 9999
         End if

      End if


      If (sf(237,1)) Then

!L   net_flux_trop

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(237,1,im_index)),net_flux_trop,    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,237,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 237)"
            goto 9999
         End if

      End if


      If (sf(238,1)) Then

!L   up_flux_trop

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(238,1,im_index)),up_flux_trop,     &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,238,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 238)"
            goto 9999
         End if

      End if



      If (sf(241,1)) Then

!L ntot_diag

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(241,1,im_index)),                &
     &        ntot_diag,                                                &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,241,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,241,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 241)"
            goto 9999
         End if

      End if



      If (sf(242,1)) Then

!L strat_lwc_diag

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(242,1,im_index)),                &
     &        strat_lwc_diag,                                           &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,242,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,242,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 242)"
            goto 9999
         End if

      End if



      If (sf(243,1)) Then

!L so4_ccn_diag

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(243,1,im_index)),                &
     &        so4_ccn_diag,                                             &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,243,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,243,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 243)"
            goto 9999
         End if

      End if



      If (sf(244,1)) Then

!L cond_samp_wgt

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(244,1,im_index)),                &
     &        cond_samp_wgt,                                            &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,244,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,244,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 244)"
            goto 9999
         End if

      End if



      If (sf(245,1)) Then

!L weighted_re

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(245,1,im_index)),                  &
     &        weighted_re,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,245,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 245)"
            goto 9999
         End if

      End if



      If (sf(246,1)) Then

!L sum_weight_re

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(246,1,im_index)),                  &
     &        sum_weight_re,                                            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,246,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 246)"
            goto 9999
         End if

      End if



      If (sf(254,1)) Then

!L weighted_warm_re

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(254,1,im_index)),                  &
     &        weighted_warm_re,                                         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,254,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 254)"
            goto 9999
         End if

      End if



      If (sf(255,1)) Then

!L sum_weight_warm_re

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(255,1,im_index)),                  &
     &        sum_weight_warm_re,                                       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,255,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 255)"
            goto 9999
         End if

      End if



      If (sf(247,1)) Then

!L sea_salt_film

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(247,1,im_index)),                &
     &        sea_salt_film,                                            &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,247,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,247,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 247)"
            goto 9999
         End if

      End if



      If (sf(248,1)) Then

!L sea_salt_jet

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(248,1,im_index)),                &
     &        sea_salt_jet,                                             &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,248,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,248,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 248)"
            goto 9999
         End if

      End if



      If (sf(280,1)) Then

!L Nc_diag

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(280,1,im_index)),                  &
     &        Nc_diag,                                                  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,280,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 280)"
            goto 9999
         End if

      End if



      If (sf(281,1)) Then

!L Nc_weight

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(281,1,im_index)),                  &
     &        Nc_weight,                                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,281,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 281)"
            goto 9999
         End if

      End if


      If (sf(262,1)) Then

!   cloud_extinction

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(262,1,im_index)),                &
     &        cloud_extinction,                                         &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,262,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,262,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 262)"
            goto 9999
         End if

      End if

      If (sf(263,1)) Then

!   cloud_weight_extinction

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(263,1,im_index)),                &
     &        cloud_weight_extinction,                                  &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,263,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,263,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 263)"
            goto 9999
         End if

      End if

      If (sf(264,1)) Then

!   ls_cloud_extinction

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(264,1,im_index)),                &
     &        ls_cloud_extinction,                                      &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,264,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,264,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 264)"
            goto 9999
         End if

      End if

      If (sf(265,1)) Then

!   ls_cloud_weight_extinction

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(265,1,im_index)),                &
     &        ls_cloud_weight_extinction,                               &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,265,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,265,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 265)"
            goto 9999
         End if

      End if

      If (sf(266,1)) Then

!   cnv_cloud_extinction

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(266,1,im_index)),                &
     &        cnv_cloud_extinction,                                     &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,266,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,266,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 266)"
            goto 9999
         End if

      End if

      If (sf(267,1)) Then

!   cnv_cloud_weight_extinction

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(267,1,im_index)),                &
     &        cnv_cloud_weight_extinction,                              &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,267,1,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,267,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 267)"
            goto 9999
         End if

      End if


      If (sf(290,1)) Then

!L photosynth_act_rad (total PAR flux at surface)

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(290,1,im_index)),                   &
     &        photosynth_act_rad,                                       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,290,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 290)"
            goto 9999
         End if

      End if

      If (sf(291,1)) Then

!L flux_direct_par (direct component of PAR flux at surface)

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(291,1,im_index)),                   &
     &        dirpar_flux,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,291,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 291)"
            goto 9999
         End if

      End if


      If (sf(292,1)) Then

         !  sol_bearing

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(292,1,im_index)),sol_bearing,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,292,                                           &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 292)"
            goto 9999
         End if

      End if


      If (sf(293,1)) Then

         !  slope_aspect

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(293,1,im_index)),slope_aspect,      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,293,                                           &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 293)"
            goto 9999
         End if

      End if


      If (sf(294,1)) Then

         !  slope_angle

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(294,1,im_index)),slope_angle,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,294,                                           &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 294)"
            goto 9999
         End if

      End if


      If (sf(295,1)) Then

         !  orog_corr

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(295,1,im_index)),orog_corr,         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,295,                                           &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 295)"
            goto 9999
         End if

      End if


      If (sf(296,1)) Then

         !  f_orog

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(296,1,im_index)),f_orog,            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,296,                                           &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 296)"
            goto 9999
         End if

      End if


 9999 continue  ! exit point on error
      If(icode /= 0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,Cmessage)
      Endif

      Return
      END SUBROUTINE diagnostics_sw

#endif
