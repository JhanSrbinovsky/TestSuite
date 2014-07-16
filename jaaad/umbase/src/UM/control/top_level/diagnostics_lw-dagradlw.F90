#if defined(A70_1B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine diagnostics_lw
!***********************************************************************
!*  This section has the following development path                   **
!*  vn2.6  ...  Basic code                                            **
!*  vn2.6a ...  code fixed/developed                                  **
!***********************************************************************

      Subroutine diagnostics_lw(                                        &
     &                       row_length, rows, model_levels             &
     &,                      wet_model_levels, ozone_levels             &
     &,                      cloud_levels                               &
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
     &,                      surflw, OLR, clear_olr                     &
     &,                      total_cloud_cover_lw                       &
     &,                      surface_down_flux_lw, surf_down_clr_lw     &
     &,                      T_incr_diagnostic                          &
     &,                      clear_hr_lw                                &
     &,                      net_flux_trop_lw, down_flux_trop_lw        &
     &,                      total_cloud_on_levels                      &
     &,                      cloud_absorptivity                         &
     &,                      cloud_weight_absorptivity                  &
     &,                      ls_cloud_absorptivity                      &
     &,                      ls_cloud_weight_absorptivity               &
     &,                      cnv_cloud_absorptivity                     &
     &,                      cnv_cloud_weight_absorptivity              &
     &,                      isccp_weights                              &
     &,                      isccp_cf                                   &
     &,                      isccp_cf_tau_0_to_p3                       &
     &,                      isccp_cf_tau_p3_to_1p3                     &
     &,                      isccp_cf_tau_1p3_to_3p6                    &
     &,                      isccp_cf_tau_3p6_to_9p4                    &
     &,                      isccp_cf_tau_9p4_to_23                     &
     &,                      isccp_cf_tau_23_to_60                      &
     &,                      isccp_cf_tau_ge_60                         &
     &,                      meanalbedocld                              &
     &,                      meantaucld                                 &
     &,                      meanptop                                   &
     &,                      totalcldarea                               &
     &,                      ls_qcl_rad                                 &
     &,                      ls_qcf_rad                                 &
     &,                      cc_qcl_rad                                 &
     &,                      cc_qcf_rad                                 &
     &,                      ls_cl_rad                                  &
     &,                      ls_cf_rad                                  &
     &,                      cc_cl_rad                                  &
     &,                      cc_cf_rad                                  &
     &,                      ozone                                      &
     &,                      O3_trop_level                              &
     &,                      O3_trop_height                             &
     &,                      T_trop_level                               &
     &,                      T_trop_height                              &
     &,                      LW_incs, LWsea                             &
     &,                      n_aod_wavel                                &
     &,                      aod_sulphate                               &
     &,                      aod_dust                                   &
     &,                      aod_seasalt                                &
     &,                      aod_soot                                   &
     &,                      aod_biomass                                &
     &,                      aod_biogenic                               &
     &,                      aod_ocff                                   &
     &,                      aod_delta                                  &
     &     ,                                                            &
#include "argsts.h"
     & STASHwork                                                        &
     &     )


! Purpose:
!          Calculates diagnostics and outputs them.
!
! Method:
!
! History:
! Date     Version     Comment
! ----     -------     -------
! 22/04/98   2.4   Distinguish between area and bulk cloud fractions in
!                  radiation calculations: bulk for Re, area otherwise.
!                  Add Total Cloud in LW diagnostic. Andrew C. Bushell
! Version   Date     Comment
! ---      -----     -------
! 5.0  30/11/99 Original version in UM. J-C Thil.
! 5.1  09/12/99 Add error trapping. Rick Rawlins
! 5.1  28/02/00 Provide explicit model increments as STASH output
!               diagnostics. R Rawlins
! 5.2  14/11/00 Provide code for reactivated diagnostics.
!                                 (J. M. Edwards)
! 5.4  03/08/01 Diagnostics for PC2 cloud scheme. D.R. Wilson
! 5.5  09/01/02 Add absorptivity and isccp diagnostics.
!                                             A.B.Keen/K.D.Williams
!
! 5.5  24/02/03 Output O3_trop_level, O3_trop_height,
!               T_trop_level, T_trop_height.
!                                      (J.-C. Thelen)
! 6.1  05/05/04 Add extra ISCCP diagnostics.     A. Jones
! 6.2  05/11/05 Add aerosol optical depth diagnostics.    N. Bellouin
! 6.4  16/01/07 Add new grid-box mean cloud diagnostics.  Adrian Lock
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
      Real                                                              &
     &  timestep         ! atmosphere timestep


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
     &, ozone_levels                                                    &
                         ! number of levels where ozone is held
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
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       &
     &, g_datastart (3,0:n_proc-1)

! Primary Arrays used in all models
      Real                                                              &
     &  T_n(row_length, rows, model_levels)                             &
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
     &, LW_incs(row_length, rows, 0:model_levels)

      Real                                                              &
     &  OLR (row_length, rows)                                          &
     &, surflw (row_length, rows)                                       &
     &, clear_olr (row_length, rows)                                    &
     &, total_cloud_cover_lw (row_length, rows)                         &
     &, surface_down_flux_lw (row_length, rows)                         &
     &, surf_down_clr_lw (row_length, rows)                             &
     &, LWsea(row_length, rows)                                         &
     &, T_incr_diagnostic(row_length,rows,model_levels)                 &
     &, clear_hr_lw(row_length,rows,model_levels)                       &
     &, net_flux_trop_lw(row_length,rows)                               &
     &, down_flux_trop_lw(row_length,rows)                              &
     &, ozone(row_length,rows,ozone_levels)                             &
     &, total_cloud_on_levels(row_length, rows, cloud_levels)           &
     &, cloud_absorptivity(row_length, rows, cloud_levels)              &
     &, cloud_weight_absorptivity(row_length, rows, cloud_levels)       &
     &, ls_cloud_absorptivity(row_length, rows, cloud_levels)           &
     &, ls_cloud_weight_absorptivity(row_length, rows, cloud_levels)    &
     &, cnv_cloud_absorptivity(row_length, rows, cloud_levels)          &
     &, cnv_cloud_weight_absorptivity(row_length, rows, cloud_levels)   &
     &, isccp_weights(row_length, rows)                                 &
     &, isccp_cf(row_length, rows, 7)                                   &
     &, isccp_cf_tau_0_to_p3(row_length, rows, 7)                       &
     &, isccp_cf_tau_p3_to_1p3(row_length, rows, 7)                     &
     &, isccp_cf_tau_1p3_to_3p6(row_length, rows, 7)                    &
     &, isccp_cf_tau_3p6_to_9p4(row_length, rows, 7)                    &
     &, isccp_cf_tau_9p4_to_23(row_length, rows, 7)                     &
     &, isccp_cf_tau_23_to_60(row_length, rows, 7)                      &
     &, isccp_cf_tau_ge_60(row_length, rows, 7)                         &
     &, meanalbedocld(row_length,rows)                                  &
     &, meantaucld(row_length,rows)                                     &
     &, meanptop(row_length,rows)                                       &
     &, totalcldarea(row_length,rows)                                   &
     &, ls_qcl_rad(row_length,rows,model_levels)                        &
     &, ls_qcf_rad(row_length,rows,model_levels)                        &
     &, cc_qcl_rad(row_length,rows,model_levels)                        &
     &, cc_qcf_rad(row_length,rows,model_levels)                        &
     &, ls_cl_rad(row_length,rows,model_levels)                         &
     &, ls_cf_rad(row_length,rows,model_levels)                         &
     &, cc_cl_rad(row_length,rows,model_levels)                         &
     &, cc_cf_rad(row_length,rows,model_levels)                         &
     &, O3_trop_level(row_length,rows)                                  &
     &, O3_trop_height(row_length,rows)                                 &
     &, T_trop_level(row_length,rows)                                   &
     &, T_trop_height(row_length,rows)

! Aerosol optical depth diagnostics
      Integer                                                           &
     &  n_aod_wavel
      Real                                                              &
     &  aod_sulphate(row_length, rows, n_aod_wavel)                     &
     &, aod_dust(row_length, rows, n_aod_wavel)                         &
     &, aod_seasalt(row_length, rows, n_aod_wavel)                      &
     &, aod_soot(row_length, rows, n_aod_wavel)                         &
     &, aod_biomass(row_length, rows, n_aod_wavel)                      &
     &, aod_biogenic(row_length, rows, n_aod_wavel)                     &
     &, aod_ocff(row_length, rows, n_aod_wavel)                         &
     &, aod_delta(row_length, rows, n_aod_wavel)

#include "csubmodl.h"
#include "typsts.h"

! Diagnostic variables
      Real                                                              &
     & STASHwork(*)    ! STASH workspace

! Local array & variables
      Real                                                              &
     &  work_3d(row_length, rows, model_levels)                         &
     &  ,isccp_dummy_3d(row_length,rows,cloud_levels)

      Integer                                                           &
     &  i, j, k                                                         &
     &,    icode                                                        &
                                ! Return code  =0 Normal exit  >1 Error
     & ,item                                                            &
                        ! STASH item
     & ,sect            ! STASH section
      Parameter( sect = 2 ) ! for lw radiation

      Character*80  cmessage
      Character(*) RoutineName
      Parameter ( RoutineName='diagnostics_lw')

      Integer                                                           &
     &  im_index        ! internal model index

! External routines
      External                                                          &
     &  copydiag, copydiag_3d                                           &
     &  ,Ereport

      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! ----------------------------------------------------------------------
! Section 1.  Diagnostic Calculation and output.
! ----------------------------------------------------------------------

! pseudo temperature after lw radiation (diagnostic only)

      item =   4  ! T + LW T_increment
      If (icode <= 0 .and. sf(item,sect)) Then

        Do k = 1,model_levels
          Do j = 1,rows
            Do i = 1,row_length
               work_3d(i,j,k) = T_n(i,j,k) + T_inc(i,j,k)
            Enddo ! i
          Enddo   ! j
        Enddo     ! k

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 004)"//cmessage
        End if

      Endif  !  sf(item,sect)

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

      If (sf(201,2)) Then

!L   surflw

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(201,2,im_index)),surflw,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,201,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 201)"
            goto 9999
         End if

      End if



      If (sf(205,2)) Then

!L   OLR

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(205,2,im_index)),OLR,              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,205,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 205)"
            goto 9999
         End if

      End if


      If (sf(206,2)) Then

!L  clear_olr

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(206,2,im_index)),clear_olr,        &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,206,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 206)"
            goto 9999
         End if

      End if




      If (sf(203,2)) Then

!L    LWsea : 'net down lw rad flux: open sea'

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(203,2,im_index)),                  &
     &        LWsea,                                                    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,203,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 203 = LWsea)"
            goto 9999
         End if

      End if


      If (sf(204,2)) Then

!L  total_cloud_cover_lw

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(204,2,im_index)),                  &
     &        total_cloud_cover_lw,                                     &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,204,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 204)"
            goto 9999
         End if

      End if



      item = 207 ! surface downward LW flux
      If (icode <= 0 .and. sf(item,sect)) Then

!L    surface_down_flux_lw

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(item,sect,im_index)),              &
     &        surface_down_flux_lw,                                     &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 207)"
            goto 9999
         End if

      End if



      item = 208 ! surface downward clear skies LW flux
      If (icode <= 0 .and. sf(item,sect)) Then

!L    surf_down_clr_lw,

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(item,sect,im_index)),              &
     &        surf_down_clr_lw,                                         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 208)"
            goto 9999
         End if

      End if




      If (sf(233,2)) Then

!L  clear_hr_lw

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(233,2,im_index)),               &
     &        clear_hr_lw,                                              &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,233,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,233,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 233)"//cmessage
            goto 9999
         End if

      End if



      If (sf(237,2)) Then

!L  net_flux_trop_lw

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(237,2,im_index)),                  &
     &        net_flux_trop_lw,                                         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,237,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 237)"
            goto 9999
         End if

      End if



      If (sf(238,2)) Then

!L  down_flux_trop_lw

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(238,2,im_index)),                  &
     &        down_flux_trop_lw,                                        &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,238,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 238)"
            goto 9999
         End if

      End if


      If (sf(260,2)) Then

!L  Ozone

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(260,2,im_index)),               &
     &        ozone,                                                    &
     &        row_length,rows,ozone_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,260,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,260,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 260)"//cmessage
            goto 9999
         End if

      End if

      If (sf(280,2)) then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(280,2,im_index)),                  &
     &        O3_trop_level,                                            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,280,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag(item 280)"//cmessage
         End if
      End if

      If (sf(281,2)) then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(281,2,im_index)),                  &
     &        O3_trop_height,                                           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,281,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag(item 281)"//cmessage
         End if
      End if

      If (sf(282,2)) then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(282,2,im_index)),                  &
     &        T_trop_level,                                             &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,282,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag(item 282)"//cmessage
         End if
      End if

      If (sf(283,2)) then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(283,2,im_index)),                  &
     &        T_trop_height,                                            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,283,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag(item 283)"//cmessage
         End if
      End if
! LW heating =
!            lw radiation temperature increment per timestep / timestep

      item = 232  ! LW heating
      If (icode <= 0 .and. sf(item,sect)) Then

        Do k = 1,model_levels
          Do j = 1,rows
            Do i = 1,row_length
               work_3d(i,j,k) = LW_incs(i,j,k)/timestep
            Enddo ! i
          Enddo   ! j
        Enddo     ! k

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 232)"//cmessage
        End if

      Endif  !  sf(item,sect)


      If (sf(261,2)) Then

!   total_cloud_on_levels

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(261,2,im_index)),               &
     &        total_cloud_on_levels,                                    &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,261,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,261,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 261)"//cmessage
            goto 9999
         End if

      End if

       If (sf(262,2)) Then

!   cloud_absorptivity

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(262,2,im_index)),               &
     &        cloud_absorptivity,                                       &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,262,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,262,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 262)"//cmessage
            goto 9999
         End if

      End if

      If (sf(263,2)) Then

!   cloud_weight_absorptivity

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(263,2,im_index)),               &
     &        cloud_weight_absorptivity,                                &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,263,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,263,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 263)"//cmessage
            goto 9999
         End if

      End if

      If (sf(264,2)) Then

!   ls_cloud_absorptivity

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(264,2,im_index)),               &
     &        ls_cloud_absorptivity,                                    &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,264,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,264,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 264)"//cmessage
            goto 9999
         End if

      End if

      If (sf(265,2)) Then

!   ls_cloud_weight_absorptivity

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(265,2,im_index)),               &
     &        ls_cloud_weight_absorptivity,                             &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,265,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,265,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 265)"//cmessage
            goto 9999
         End if

      End if

      If (sf(266,2)) Then

!   cnv_cloud_absorptivity

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(266,2,im_index)),               &
     &        cnv_cloud_absorptivity,                                   &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,266,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,266,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 266)"//cmessage
            goto 9999
         End if

      End if

      If (sf(267,2)) Then

!   cnv_cloud_weight_absorptivity

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(267,2,im_index)),               &
     &        cnv_cloud_weight_absorptivity,                            &
     &        row_length,rows,cloud_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,267,2,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,2,267,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 267)"//cmessage
            goto 9999
         End if

      End if

      If (sf(269,2)) Then

!   isccp_weights

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(269,2,im_index)),                  &
     &        isccp_weights,                                            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,269,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 269)"
            goto 9999
         End if

      End if

      If (sf(290,2)) Then

!   meanalbedocld

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(290,2,im_index)),                  &
     &        meanalbedocld,                                            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,290,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 290)"
            goto 9999
         End if

      End if

      If (sf(291,2)) Then

!   meantaucld

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(291,2,im_index)),                  &
     &        meantaucld,                                               &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,291,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 291)"
            goto 9999
         End if

      End if

      If (sf(292,2)) Then

!   meanptop

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(292,2,im_index)),                  &
     &        meanptop,                                                 &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,292,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 292)"
            goto 9999
         End if

      End if

      If (sf(293,2)) Then

!   totalcldarea

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(293,2,im_index)),                  &
     &        totalcldarea,                                             &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,293,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 293)"
            goto 9999
         End if

      End if

!   Cloud water mixing ratios
      item = 308  ! LS cloud liquid water mixing ratio
      If (icode.le.0 .and. sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        ls_qcl_rad,                                               &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 308)"//cmessage
         End if
      Endif
      item = 309  ! LS cloud ice water mixing ratio
      If (icode.le.0 .and. sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        ls_qcf_rad,                                               &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 309)"//cmessage
         End if
      Endif
      item = 310  ! Convective cloud liquid water mixing ratio
      If (icode.le.0 .and. sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        cc_qcl_rad,                                               &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 310)"//cmessage
         End if
      Endif
      item = 311  ! Convective cloud ice water mixing ratio
      If (icode.le.0 .and. sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        cc_qcf_rad,                                               &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 311)"//cmessage
         End if
      Endif
!   Cloud amounts
      item = 312  ! LS cloud fraction of grdbox seen by radiation.
                  ! Liquid
      If (icode.le.0 .and. sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        ls_cl_rad,                                                &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 312)"//cmessage
         End if
      Endif
      item = 313  ! LS cloud fraction of grdbox seen by radiation.
                  ! Ice
      If (icode.le.0 .and. sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        ls_cf_rad,                                                &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 313)"//cmessage
         End if
      Endif
      item = 314  ! CONV cloud fraction of grdbox seen by radiation.
                  ! Liquid
      If (icode.le.0 .and. sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        cc_cl_rad,                                                &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 314)"//cmessage
         End if
      Endif
      item = 315  ! CONV cloud fraction of grdbox seen by radiation.
                  ! Ice
      If (icode.le.0 .and. sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        cc_cf_rad,                                                &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 315)"//cmessage
         End if
      Endif

! Copy ISCCP diagnostics by looping over 7 levels in call to copydiag.
! This is because copydiag_3d cannot handle ISCCP levels.

      If (sf(270,2)) Then

!   isccp_cf

        Do k = 1,7
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(270,2,im_index)                    &
     &         +(row_length*rows*(k-1))),                               &
     &        isccp_cf(1,1,k),                                          &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,270,                                           &
     &        icode,cmessage)
        Enddo


         If (icode  >   0) then
            cmessage=": error in copydiag( item 270)"//cmessage
            goto 9999
         End if

      End if

      If (sf(271,2)) Then

!   isccp_cf_tau_0_to_p3

        Do k = 1,7
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(271,2,im_index)                    &
     &         +(row_length*rows*(k-1))),                               &
     &        isccp_cf_tau_0_to_p3(1,1,k),                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,271,                                           &
     &        icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 271)"//cmessage
            goto 9999
         End if

      End if

      If (sf(272,2)) Then

!   isccp_cf_tau_p3_to_1p3

        Do k = 1,7
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(272,2,im_index)                    &
     &         +(row_length*rows*(k-1))),                               &
     &        isccp_cf_tau_p3_to_1p3(1,1,k),                            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,272,                                           &
     &        icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 272)"//cmessage
            goto 9999
         End if

      End if

      If (sf(273,2)) Then

!   isccp_cf_tau_1p3_to_3p6

        Do k = 1,7
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(273,2,im_index)                    &
     &         +(row_length*rows*(k-1))),                               &
     &        isccp_cf_tau_1p3_to_3p6(1,1,k),                           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,273,                                           &
     &        icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 273)"//cmessage
            goto 9999
         End if

      End if

      If (sf(274,2)) Then

!   isccp_cf_tau_3p6_to_9p4

        Do k = 1,7
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(274,2,im_index)                    &
     &         +(row_length*rows*(k-1))),                               &
     &        isccp_cf_tau_3p6_to_9p4(1,1,k),                           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,274,                                           &
     &        icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 274)"//cmessage
            goto 9999
         End if

      End if

      If (sf(275,2)) Then

!   isccp_cf_tau_9p4_to_23

        Do k = 1,7
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(275,2,im_index)                    &
     &         +(row_length*rows*(k-1))),                               &
     &        isccp_cf_tau_9p4_to_23(1,1,k),                            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,275,                                           &
     &        icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 275)"//cmessage
            goto 9999
         End if

      End if

      If (sf(276,2)) Then

!   isccp_cf_tau_23_to_60

        Do k = 1,7
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(276,2,im_index)                    &
     &         +(row_length*rows*(k-1))),                               &
     &        isccp_cf_tau_23_to_60(1,1,k),                             &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,276,                                           &
     &        icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 276)"//cmessage
            goto 9999
         End if

      End if

      If (sf(277,2)) Then

!   isccp_cf_tau_ge_60

        Do k = 1,7
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(277,2,im_index)                    &
     &         +(row_length*rows*(k-1))),                               &
     &        isccp_cf_tau_ge_60(1,1,k),                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,277,                                           &
     &        icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 277)"//cmessage
            goto 9999
         End if

      End if

!   Aerosol optical depth diagnostics
!   (loop on wavelength)

      If (sf(284,2)) Then
        Do k = 1, n_aod_wavel
! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(284,2,im_index)                   &
     &        +(row_length*rows*(k-1))),                                &
     &        aod_sulphate(1,1,k),                                      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,284,                                           &
     &        icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 284)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

      If (sf(285,2)) Then
        Do k = 1, n_aod_wavel
! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(285,2,im_index)                   &
     &        +(row_length*rows*(k-1))),                                &
     &        aod_dust(1,1,k),                                          &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,285,                                           &
     &        icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 285)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

      If (sf(286,2)) Then
        Do k = 1, n_aod_wavel
! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(286,2,im_index)                   &
     &        +(row_length*rows*(k-1))),                                &
     &        aod_seasalt(1,1,k),                                       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,286,                                           &
     &        icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 286)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

      If (sf(287,2)) Then
        Do k = 1, n_aod_wavel
! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(287,2,im_index)                   &
     &        +(row_length*rows*(k-1))),                                &
     &        aod_soot(1,1,k),                                          &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,287,                                           &
     &        icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 287)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

      If (sf(288,2)) Then
        Do k = 1, n_aod_wavel
! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(288,2,im_index)                   &
     &        +(row_length*rows*(k-1))),                                &
     &        aod_biomass(1,1,k),                                       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,288,                                           &
     &        icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 288)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

      If (sf(289,2)) Then
        Do k = 1, n_aod_wavel
! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(289,2,im_index)                   &
     &        +(row_length*rows*(k-1))),                                &
     &        aod_biogenic(1,1,k),                                      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,289,                                           &
     &        icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 289)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

      If (sf(295,2)) Then
        Do k = 1, n_aod_wavel
! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(295,2,im_index)                   &
     &        +(row_length*rows*(k-1))),                                &
     &        aod_ocff(1,1,k),                                          &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,295,                                           &
     &        icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 295)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

      If (sf(296,2)) Then
        Do k = 1, n_aod_wavel
! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(296,2,im_index)                   &
     &        +(row_length*rows*(k-1))),                                &
     &        aod_delta(1,1,k),                                         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,296,                                           &
     &        icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 296)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

 9999 continue  ! exit point on error
      If(icode /= 0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,Cmessage)
      Endif

      Return
      END SUBROUTINE diagnostics_lw
#endif
