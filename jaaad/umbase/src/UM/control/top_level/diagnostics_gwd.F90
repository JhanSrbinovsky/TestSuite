#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine diagnostics_gwd

      Subroutine diagnostics_gwd(                                       &
     &                       row_length, rows, model_levels             &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      at_extremity                               &
     &,                      u, v, R_u, R_v                             &
     &,                      u_incr_diagnostic,v_incr_diagnostic        &
     &,                      stress_ud     ,  stress_vd                 &
     &,                      stress_ud_satn,  stress_vd_satn            &
     &,                      stress_ud_wake,  stress_vd_wake            &
     &,                      du_dt_satn    ,  dv_dt_satn                &
     &,                      du_dt_wake   ,  dv_dt_wake                 &
     &,                      u_s_d, v_s_d, nsq_s_d                      &
     &,                      num_lim_d, num_fac_d                       &
     &,                      fr_d, bld_d ,bldt_d                        &
     &, tausx_d, tausy_d, taus_scale_d                                  &
     &, sd_orog_land, land_sea_mask, land_points                        &
     &, orog_grad_xx_land, orog_grad_xy_land, orog_grad_yy_land         &
     &, GWSPEC_EFLUX, GWSPEC_SFLUX, GWSPEC_WFLUX                        &
     &, GWSPEC_NFLUX, GWSPEC_EWACC, GWSPEC_NSACC,                       &
#include "argsts.h"
     &     STASHwork                                                    &
     &     )


! Purpose:
!          Calculates diagnostics and outputs them.
!
! Method:
!   Item by item copying of diagnostic arrays made available from gwd
! scheme into the corresponding STASHwork (master array for output)
! diagnostic space. All sizes and address pointers have been
! calculated by the STASH initialisation routines previously. Each
! diagnostic is protected by a STASHflag logical switch which is
! re-set every timestep.
!   Note that diagnostics for du_dt,dv_dt are available on all
! model rho levels. Stress diagnostics are now available on all model
! theta_levels - running from level 0 at the ground to
! level=model_levels at the model lid.
! New single level diagnostics (6214-6222) are output at theta points
! in the horizontal, i.e. they are output at the points on which the
! GWD scheme operates. 6233 and 6234 are the same.
!
!
! History:
! Version   Date     Comment
! ----     -------     -------
! 5.0  30/11/99 Original version. J-C Thil.
! 5.1  09/12/99 Correct interfaces into copydiag_3D for 003 and 201.
!                  Rick Rawlins
! 5.1  15/03/00 Correct no. of vertical levels for most gwd
!               diagnostics. R Rawlins
! 5.1  28/02/00 Provide explicit model increments as STASH output
!               diagnostics. R Rawlins
! 5.2  15/11/00 Calculate stress diagnostics 6201,6202 for all
!               model levels (previously start_level to
!               model_levels-1). Add stress diagnostics for individual
!               components of the GWD scheme (STASH nos 6223-6230)
!               Add blocked flow wake drag wind increments
!               diagnostics (6231,6232).
!               Add 9 single level diagnostics of the GWD surface
!               calculations (6214-6222).
!                                               Stuart Webster.
! 5.3  11/10/01 Calculate orographic standard deviation diagnostic 06203
!               from compressed field. D.M. Goddard
! 5.4  28/08/02 Add numerical limiter diagnostics 6233,6234. S. Webster
! 5.4  05/09/02 Add diagnostics for spectral (non-orographic) gravity
!               wave forcing.       Adam Scaife
! 5.5  25/02/03 Remove 3B GWD diagnostics.   Stuart Webster.
! 6.2  21/02/06 Add orog surface stress and sigma xx, xy and yy
!               diagnostics.                 Stuart Webster.
! 6.4  19/05/06 Correction of USSP diagnostic flux output.
!                                                   A.C. Bushell
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

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
     &, number_format                                                   &
                         ! switch controlling number format diagnostics
                         ! are written out in. See PP_WRITE for details.
     &, model_domain                                                    &
                         ! indicator as to model type, ie global, lam
     &, land_points               ! Number of land points

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
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, R_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &      model_levels)                                               &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &      model_levels)

      Real                                                              &
     &  stress_ud      (row_length,  rows,0:model_levels)               &
     &, stress_vd      (row_length,n_rows,0:model_levels)               &
     &, stress_ud_satn (row_length,  rows,0:model_levels)               &
     &, stress_vd_satn (row_length,n_rows,0:model_levels)               &
     &, stress_ud_wake (row_length,  rows,0:model_levels)               &
     &, stress_vd_wake (row_length,n_rows,0:model_levels)               &
     &, du_dt_satn     (row_length,  rows,  model_levels)               &
     &, dv_dt_satn     (row_length,n_rows,  model_levels)               &
     &, du_dt_wake     (row_length, rows,   model_levels)               &
     &, dv_dt_wake     (row_length,n_rows,  model_levels)               &
     &, GWSPEC_EFLUX(row_length,rows,model_levels)                      &
     &, GWSPEC_SFLUX(row_length,n_rows,model_levels)                    &
     &, GWSPEC_WFLUX(row_length,rows,model_levels)                      &
     &, GWSPEC_NFLUX(row_length,n_rows,model_levels)                    &
     &,GWSPEC_EWACC(row_length,rows,model_levels)                       &
     &,GWSPEC_NSACC(row_length,n_rows,model_levels)                     &
! Note: space is only allocated for _incr_diagnostic arrays in calling
!       routine if coresponding STASHflags are set .true.
     &, u_incr_diagnostic(row_length,  rows,model_levels)               &
     &, v_incr_diagnostic(row_length,n_rows,model_levels)               &
     &, num_lim_d(row_length, rows)                                     &
     &, num_fac_d(row_length, rows)                                     &
     &, u_s_d  (row_length, rows)                                       &
     &, v_s_d  (row_length, rows)                                       &
     &, nsq_s_d  (row_length, rows)                                     &
     &, fr_d   (row_length, rows)                                       &
     &, bld_d  (row_length, rows)                                       &
     &, bldt_d (row_length, rows)                                       &
     &, tausx_d (row_length, rows)                                      &
     &, tausy_d (row_length, n_rows)                                    &
     &, taus_scale_d (row_length, rows)                                 &
     &, sd_orog_land(land_points)                                       &
                                  ! orographic std dev on land points
     &, orog_grad_xx_land(land_points)                                  &
                                       ! sigma_xx on land points
     &, orog_grad_xy_land(land_points)                                  &
                                       ! sigma_xy on land points
     &, orog_grad_yy_land(land_points) ! sigma_yy on land points

      LOGICAL, INTENT(IN) ::                                            &
     &  land_sea_mask (land_points)  ! land_sea_mask

#include "csubmodl.h"
#include "typsts.h"

! Diagnostic variables
      Real                                                              &
     & STASHwork(*)     ! STASH workspace

! Local variables
      Integer                                                           &
     &  i, j, k, l                                                      &
                      ! loop indices
     &,    icode                                                        &
                                ! Return code  =0 Normal exit  >1 Error
     & ,item                                                            &
                        ! STASH item
     & ,sect            ! STASH section
      Parameter( sect = 6 ) ! for gravity wave drag

      Character*80  cmessage
      Character(*) RoutineName
      Parameter ( RoutineName='diagnostics_gwd')

      Integer                                                           &
     &  im_index        ! internal model index

      Real                                                              &
     & interp_data_3(row_length*rows*model_levels)                      &
     &, sd_orog     (row_length,rows)                                   &
                                       ! orog std dev full field
     &, orog_grad_xx(row_length,rows)                                   &
                                       !     sigma_xx full field
     &, orog_grad_xy(row_length,rows)                                   &
                                       !     sigma_xy full field
     &, orog_grad_yy(row_length,rows)  !     sigma_yy full field

! External routines
      External                                                          &
     &   copydiag, copydiag_3d                                          &
     &  ,Ereport                                                        &
     &, From_land_points

! ND pp codes temporarily kept around :
!
!      PP_code_r                =   1
!      PP_code_p                =   8
!      PP_code_eta              =  10   !Note: currently no PP code
!!                             for eta, so set to 10 which is sigma.
!      PP_code_T                =  16
!      PP_code_msl              = 128
!      PP_code_surf             = 129
!      PP_code_du_dt_satn       =  68
!      PP_code_dv_dt_satn       =  69
!      PP_code_gwd_stress_u     = 161
!      PP_code_gwd_stress_v     = 162
!      PP_code_u_tend           = 975
!      PP_code_v_tend           = 976
!
      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! ----------------------------------------------------------------------
! Section 1.  Diagnostic Calculation and output.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 1.1 Diagnostics for GWD
! ----------------------------------------------------------------------

!-----------------------------------------------------------------------
! u,v increment diagnostics= modified - previous
!-----------------------------------------------------------------------

      item = 185  ! u increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k=1,model_levels
          Do j=1,rows
           Do i=1,row_length
            u_incr_diagnostic(i,j,k) = R_u(i,j,k) -                     &
     &                                      u_incr_diagnostic(i,j,k)
           End Do !i
          End Do  !j
         End Do   !k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        u_incr_diagnostic,                                        &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 185)"//cmessage
         End if

      Endif  !  sf(item,sect)

      item = 186  ! v increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k=1,model_levels
          Do j=1,n_rows
           Do i=1,row_length
            v_incr_diagnostic(i,j,k) = R_v(i,j,k) -                     &
     &                                      v_incr_diagnostic(i,j,k)
           End Do !i
          End Do  !j
         End Do   !k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        v_incr_diagnostic,                                        &
     &        row_length,n_rows,model_levels,0,0,0,0, at_extremity,     &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 186)"//cmessage
         End if

      Endif  !  sf(item,sect)

! ----------------------------------------------------------------------
!  U component of stress
! ----------------------------------------------------------------------
! Item 201 stress_ud

      If (sf(201,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(201,6,im_index)),stress_ud,      &
     &        row_length,rows,model_levels+1,0,0,0,0,                   &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,201,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,201,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(stress_ud)"
            goto 9999
         Endif

      End if


! ----------------------------------------------------------------------
!  V component of stress
! ----------------------------------------------------------------------
! Item 202 stress_vd

      If (sf(202,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(202,6,im_index)),stress_vd,      &
     &        row_length,n_rows,model_levels+1,0,0,0,0,                 &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,202,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,202,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(stress_vd)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  Orographic standard deviation
! ----------------------------------------------------------------------
! Item 203 Orographic standard deviation

      IF (sf(203,6)) THEN

! DEPENDS ON: from_land_points
         CALL From_land_points(sd_orog,sd_orog_land, land_sea_mask,     &
     &        row_length*rows,land_points)

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(203,6,im_index)),                   &
     &        sd_orog,                                                  &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,6,203,                                           &
     &        icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag(sd_orog)"
            goto 9999
         END IF

      END IF

! ----------------------------------------------------------------------
!  Orographic sigma_xx field
! ----------------------------------------------------------------------
! Item 204 Orographic sigma_xx field

      IF (sf(204,6)) THEN

! DEPENDS ON: from_land_points
         CALL From_land_points(orog_grad_xx,orog_grad_xx_land,          &
     &                land_sea_mask,row_length*rows,land_points)

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(204,6,im_index)),                   &
     &        orog_grad_xx,                                             &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,6,204,                                           &
     &        icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag(orog_grad_xx)"
            goto 9999
         END IF

      END IF

! ----------------------------------------------------------------------
!  Orographic sigma_xy field
! ----------------------------------------------------------------------
! Item 205 Orographic sigma_xy field

      IF (sf(205,6)) THEN

! DEPENDS ON: from_land_points
         CALL From_land_points(orog_grad_xy,orog_grad_xy_land,          &
     &                land_sea_mask,row_length*rows,land_points)

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(205,6,im_index)),                   &
     &        orog_grad_xy,                                             &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,6,205,                                           &
     &        icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag(orog_grad_xy)"
            goto 9999
         END IF

      END IF

! ----------------------------------------------------------------------
!  Orographic sigma_yy field
! ----------------------------------------------------------------------
! Item 206 Orographic sigma_yy field

      IF (sf(206,6)) THEN

! DEPENDS ON: from_land_points
         CALL From_land_points(orog_grad_yy,orog_grad_yy_land,          &
     &                land_sea_mask,row_length*rows,land_points)

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(206,6,im_index)),                   &
     &        orog_grad_yy,                                             &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,6,206,                                           &
     &        icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag(orog_grad_yy)"
            goto 9999
         END IF

      END IF

! ----------------------------------------------------------------------
!  U component of satn stress
! ----------------------------------------------------------------------
! Item 223 stress_ud_satn

      If (sf(223,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(223,6,im_index)),stress_ud_satn, &
     &        row_length,rows,model_levels+1,0,0,0,0,                   &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,223,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,223,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(stress_ud_satn)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  V component of satn stress
! ----------------------------------------------------------------------
! Item 224 stress_vd_satn

      If (sf(224,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(224,6,im_index)),stress_vd_satn, &
     &        row_length,n_rows,model_levels+1,0,0,0,0,                 &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,224,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,224,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(stress_vd_satn)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  U component  of wake stress
! ----------------------------------------------------------------------
! Item 227 stress_ud_wake

      If (sf(227,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(227,6,im_index)),stress_ud_wake, &
     &        row_length,rows,model_levels+1,0,0,0,0,                   &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,227,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,227,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(stress_ud_wake)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  V component of wake stress
! ----------------------------------------------------------------------
! Item 228 stress_vd_wake

      If (sf(228,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(228,6,im_index)),stress_vd_wake, &
     &        row_length,n_rows,model_levels+1,0,0,0,0,                 &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,228,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,228,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(stress_vd_wake)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
! du/dt saturation component
! ----------------------------------------------------------------------
! Item 207 du_dt_satn

      If (sf(207,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(207,6,im_index)),du_dt_satn,     &
     &        row_length,rows,model_levels,0,0,0,0,                     &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,207,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,207,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(du_dt_satn)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  dv/dt saturation component
! ----------------------------------------------------------------------
! Item 208 dv_dt_satn

      If (sf(208,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(208,6,im_index)),dv_dt_satn,     &
     &        row_length,n_rows,model_levels,0,0,0,0,                   &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,208,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,208,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(dv_dt_satn)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  du/dt wake component
! ----------------------------------------------------------------------
! Item 231 du_dt_wake

      If (sf(231,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(231,6,im_index)),du_dt_wake,     &
     &        row_length,rows,model_levels,0,0,0,0,                     &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,231,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,231,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(du_dt_wake)"
            goto 9999
         Endif

      End if


! ----------------------------------------------------------------------
!  dv/dt wake component
! ----------------------------------------------------------------------
! Item 232 dv_dt_wake

      If (sf(232,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(232,6,im_index)),dv_dt_wake,     &
     &        row_length,n_rows,model_levels,0,0,0,0,                   &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,232,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,232,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(dv_dt_wake)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  u latest.
! ----------------------------------------------------------------------
! Item 002 u

      If (sf(002,6)) Then

         Do k = 1, model_levels
            Do j = 1, rows
               Do i = 1, row_length
                 l = i + (j-1)*row_length + (k-1)*rows*row_length
                 interp_data_3(l) = u(i,j,k) + R_u(i,j,k)
               End Do
            End Do
         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(002,6,im_index)),interp_data_3,  &
     &        row_length,rows,model_levels,0,0,0,0,                     &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,002,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,002,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(u+R_u)"
            goto 9999
         Endif

      End if



! ----------------------------------------------------------------------
!  v latest
! ----------------------------------------------------------------------
! Item 003 v

      If (sf(003,6)) Then

         Do k = 1, model_levels
            Do j = 1, n_rows
               Do i = 1, row_length
                 l = i + (j-1)*row_length + (k-1)*n_rows*row_length
                 interp_data_3(l) = v(i,j,k) + R_v(i,j,k)
               End Do
            End Do
         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(003,6,im_index)),interp_data_3,  &
     &        row_length,n_rows,model_levels,0,0,0,0,                   &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,003,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,003,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(v+R_v)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  U_S diag
! ----------------------------------------------------------------------
! Item 214 U_S_D

      If (sf(214,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(214,6,im_index)),                   &
     &        U_S_D,                                                    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,6,214,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(u_s)"
            goto 9999
         Endif

      End if
! ----------------------------------------------------------------------
!  V_S diag
! ----------------------------------------------------------------------
! Item 215 V_S_d

      If (sf(215,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(215,6,im_index)),                   &
     &        V_S_D,                                                    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,6,215,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(v_s)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  NSQ_S diag
! ----------------------------------------------------------------------
! Item 216 nsq_s_d

      If (sf(216,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(216,6,im_index)),                   &
     &        NSQ_S_D,                                                  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,6,216,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(nsq_s)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  Froude number diag
! ----------------------------------------------------------------------
! Item 217 FR_d

      If (sf(217,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(217,6,im_index)),                   &
     &        FR_D,                                                     &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,6,217,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(Fr)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  Blocked Layer Depth diag
! ----------------------------------------------------------------------
! Item 218 BLD_d

      If (sf(218,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(218,6,im_index)),                   &
     &        BLD_D,                                                    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,6,218,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(BLD)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  % of time of blocked layer
! ----------------------------------------------------------------------
! Item 222 BLDT_D

      If (sf(222,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(222,6,im_index)),                   &
     &        BLDT_D,                                                   &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,6,222,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(BLDT)"

            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  % of time numerical limiter invoked
! ----------------------------------------------------------------------
! Item 233 NUM_LIM_D

      If (sf(233,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(233,6,im_index)),                   &
     &        num_lim_d,                                                &
     &    row_length,rows,0,0,0,0, at_extremity,                        &
     &        atmos_im,6,233,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(num_lim)"

            goto 9999
         Endif

      End if



! ----------------------------------------------------------------------
!  % redn. of flow-blocking stress after numerical limiter invoked
! ----------------------------------------------------------------------
! Item 234 NUM_FAC_D

      If (sf(234,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(234,6,im_index)),                   &
     &        num_fac_d,                                                &
     &    row_length,rows,0,0,0,0, at_extremity,                        &
     &        atmos_im,6,234,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(num_fac)"

            goto 9999
         Endif

      End if


! ----------------------------------------------------------------------
!  Orographic X-component of the total surface stress
!  Same as 6201, level 1 output, but here the pp headers are correctly
!  labelled
! ----------------------------------------------------------------------
! Item 235 TAUSX_D

      If (sf(235,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(235,6,im_index)),                   &
     &        tausx_d,                                                  &
     &    row_length,rows,0,0,0,0, at_extremity,                        &
     &        atmos_im,6,235,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(tausx)"

            goto 9999
         Endif

      End if


! ----------------------------------------------------------------------
!  Orographic Y-component of the total surface stress
!  Same as 6202, level 1 output, but here the pp headers are correctly
!  labelled
! ----------------------------------------------------------------------
! Item 236 TAUSY_D

      If (sf(236,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(236,6,im_index)),                   &
     &        tausy_d,                                                  &
     &    row_length,n_rows,0,0,0,0, at_extremity,                      &
     &        atmos_im,6,236,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(tausy)"

            goto 9999
         Endif

      End if


! ----------------------------------------------------------------------
!  Orographic surface stress Froude number dependent scaling factor
! ----------------------------------------------------------------------
! Item 237 TAUS_SCALE_D

      If (sf(237,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(237,6,im_index)),                   &
     &        taus_scale_d,                                             &
     &    row_length,rows,0,0,0,0, at_extremity,                        &
     &        atmos_im,6,237,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(taus_scale)"

            goto 9999
         Endif

      End if


!---------------------------------------------------------------------
! Section 2.
! SPECTRAL (NON-OROGRAPHIC) GRAVITY WAVE SCHEME DIAGNOSTICS
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!  USSP Eastward flux component
! --------------------------------------------------------------------
! Item 101 GWSPEC_EFLUX(i,j,k)
      If (sf(101,6)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(101,6,im_index)),GWSPEC_EFLUX,   &
     &        row_length,rows,model_levels,0,0,0,0,                     &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,101,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,101,                                           &
     &        icode,cmessage)
         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(GWSPEC_EFLUX)"
            goto 9999
         Endif
      Endif
! --------------------------------------------------------------------
!  USSP Southward flux component
! --------------------------------------------------------------------
! Item 102 GWSPEC_SFLUX(i,j,k)
      If (sf(102,6)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(102,6,im_index)),GWSPEC_SFLUX,   &
     &        row_length,n_rows,model_levels,0,0,0,0,                   &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,102,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,102,                                           &
     &        icode,cmessage)
         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(GWSPEC_SFLUX)"
            goto 9999
         Endif
      Endif
! --------------------------------------------------------------------
!  USSP Westward flux component
! --------------------------------------------------------------------
! Item 103 GWSPEC_WFLUX(i,j,k)
      If (sf(103,6)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(103,6,im_index)),GWSPEC_WFLUX,   &
     &        row_length,rows,model_levels,0,0,0,0,                     &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,103,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,103,                                           &
     &        icode,cmessage)
         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(GWSPEC_WFLUX)"
            goto 9999
         Endif
      Endif
! --------------------------------------------------------------------
!  USSP Northward flux component
! --------------------------------------------------------------------
! Item 104 GWSPEC_NFLUX(i,j,k)
      If (sf(104,6)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(104,6,im_index)),GWSPEC_NFLUX,   &
     &        row_length,n_rows,model_levels,0,0,0,0,                   &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,104,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,104,                                           &
     &        icode,cmessage)
         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(GWSPEC_NFLUX)"
            goto 9999
         Endif
      Endif
! --------------------------------------------------------------------
!  USSP EW acceleration component
! --------------------------------------------------------------------
! Item 105 GWSPEC_EWACC(i,j,k)
      If (sf(105,6)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(105,6,im_index)),GWSPEC_EWACC,   &
     &        row_length,rows,model_levels,0,0,0,0,                     &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,105,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,105,                                           &
     &        icode,cmessage)
         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(GWSPEC_EWACC)"
            goto 9999
         Endif
      Endif
! --------------------------------------------------------------------
!  USSP NS acceleration component
! --------------------------------------------------------------------
! Item 106 GWSPEC_NSACC(i,j,k)
      If (sf(106,6)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(106,6,im_index)),GWSPEC_NSACC,   &
     &        row_length,n_rows,model_levels,0,0,0,0,                   &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,106,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,106,                                           &
     &        icode,cmessage)
         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(GWSPEC_NSACC)"
            goto 9999
         Endif
      Endif

 9999 continue                  !  single point error handling.
      If(icode /= 0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,Cmessage)
      Endif

      Return
      END SUBROUTINE diagnostics_gwd
#endif
