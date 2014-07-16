#if defined(A70_1C) || defined(A70_1Z)
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT******************************

      Subroutine diagnostics_lw(                                        &
     &      row_length, rows, model_levels                              & 
     &,     wet_model_levels, ozone_levels                              &
     &,     cloud_levels                                                & 
     &,     n_rows, global_row_length, global_rows                      &
     &,     halo_i, halo_j, off_x, off_y, me                            &
     &,     n_proc, n_procx, n_procy                                    &
     &,     g_rows, g_row_length, g_datastart                           &
     &,     at_extremity                                                &
     &,     timestep,i_off                                              & 
     &,     T_n, T_inc                                                  & 
     &,     q_n, qcl_n, cf_n, cfl_n                                     &
     &,     T_latest, q_latest, qcl_latest                              &
     &,     cf_latest, cfl_latest                                       &
     &,     surflw, OLR                                                 &
     &,     T_incr_diagnostic                                           &
     &,     n_channel                                                   &
     &,     ozone                                                       &
     &,     O3_trop_level                                               &
     &,     O3_trop_height                                              &
     &,     T_trop_level                                                &
     &,     T_trop_height                                               &
     &,     LW_incs, LWsea                                              &
     &,     n_aod_wavel                                                 &
     &,     l_out_lwdiag, j_lw                                          &
     &,                                                                 &
#include "argsts.h"
     &      STASHwork)


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
! 6.2  05/11/05 Add aerosol optical depth diagnostics.    N. Bellouin
! 6.2  19/01/06 Add radiative forcing and radiance diagnostics
!                                      (J.-C. Thelen)
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Use lw_diag_mod, Only: LW_diag

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Logical ::  at_extremity(4)  
!        Indicates if this processor is at north,
!        south, east or west of the processor grid.

      Logical :: l_out_lwdiag
!        Flag to output radiation diagnostics (mainly
!         required for time-stepping scheme.

      Real :: timestep         
!         Atmosphere timestep

      Integer :: row_length       ! number of points on a row
      Integer :: rows             ! number of rows in a theta field
      Integer :: n_rows           ! number of rows in a v field
      Integer :: model_levels     ! number of model levels
      Integer :: cloud_levels     ! number of cloudy levels
      Integer :: wet_model_levels ! number of model levels where moisture
                                  ! variables are held
      Integer :: ozone_levels     ! number of levels where ozone is held
      Integer :: number_format    ! switch controlling number format diagnostics
                                  ! are written out in. See PP_WRITE for details.
      Integer :: model_domain     ! indicator as to model type, ie global, lam
      Integer :: n_channel        ! Number of satellite channels used

      Integer :: global_row_length   !IN. NUMBER OF points on a global row
      Integer :: global_rows         !IN. NUMBER OF global rows
      Integer :: me                  !IN. Processor number
      Integer :: halo_i              !IN. size of large halo in x direction
      Integer :: halo_j              !IN. size of large halo in y direction
      Integer :: off_x               !IN. size of small halo in x direction
      Integer :: off_y               !IN. size of small halo in y direction
      Integer :: n_proc
      Integer :: n_procx
      Integer :: n_procy
      Integer :: g_rows (0:n_proc-1)
      Integer :: g_row_length (0:n_proc-1)
      Integer :: g_datastart (3,0:n_proc-1)

      Integer, Intent(IN) :: i_off
!          Offset for diagnostics
      Integer, Intent(IN) ::  n_aod_wavel
!          Aerosol optical depth diagnostics
      Integer, Intent(IN) :: j_lw
!       Call to LW radiation

!
! Primary Arrays used in all models
!
      Real :: T_n(row_length, rows, model_levels)
      Real :: T_inc(row_length, rows, model_levels)
      Real :: q_n(row_length, rows, wet_model_levels)
      Real :: qcl_n(row_length, rows, wet_model_levels)
      Real :: cf_n(row_length, rows, wet_model_levels)
      Real :: cfl_n(row_length, rows, wet_model_levels)
      Real :: T_latest(row_length, rows, model_levels)
      Real :: q_latest(row_length, rows, wet_model_levels)
      Real :: qcl_latest(row_length, rows, wet_model_levels)
      Real :: cf_latest(row_length, rows, wet_model_levels)
      Real :: cfl_latest(row_length, rows, wet_model_levels)

      Real :: OLR (row_length, rows)
      Real :: surflw (row_length, rows)
      Real :: LWsea(row_length, rows)
      Real :: T_incr_diagnostic(row_length,rows,model_levels)
      Real :: ozone(row_length,rows,ozone_levels)
      Real :: O3_trop_level(row_length,rows)
      Real :: O3_trop_height(row_length,rows)
      Real :: T_trop_level(row_length,rows)
      Real :: T_trop_height(row_length,rows)
      Real :: LW_incs(row_length, rows, 0:model_levels)

#include "csubmodl.h"
#include "typsts.h"

! Diagnostic variables
      
      Real :: STASHwork(*)    ! STASH workspace

! Local array & variables
      
      Real :: work_3d(row_length, rows, model_levels)
      Real :: isccp_dummy_3d(row_length,rows,cloud_levels)

      Integer :: i, j, k
      Integer :: icode                ! Return code  =0 Normal exit  >1 Error
      Integer :: item                 ! STASH item
      Integer, Parameter :: sect=2    ! STASH LW Radiation section

      CHARACTER(LEN=1024)  :: cmessage
      Character(*) RoutineName
      Parameter ( RoutineName='diagnostics_lw')

      Integer :: im_index        ! internal model index
      Integer :: ptr_stash       ! Pointer to position in STASH


      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! ----------------------------------------------------------------------
! Section 1.  Diagnostic Calculation and output.
! ----------------------------------------------------------------------

! First output diagnostics that are not contained in the LW_DIAG
! structure.

!
! Pseudo temperature after lw radiation (diagnostic only)
!
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
     &       work_3d,                                                   &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 004)"//cmessage
        End if

      Endif  !  sf(item,sect)

!
! Increment diagnostics= modified - previous
!
      item = 161  ! temperature increment
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &       T_incr_diagnostic,                                         &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 161)"//cmessage
         End if

      Endif  !  sf(item,sect)

!
! Temperature increment including condensation
!
      item = 181  
      If (icode <= 0 .and. sf(item,sect)) Then

      Do k = 1,model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = T_latest(i,j,k) - T_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &       work_3d,                                                   &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        & 
     &       icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 181)"//cmessage
         End if

      Endif  !  sf(item,sect)

!
! Vapour increment
!
      item = 182  
      If (icode <= 0 .and. sf(item,sect)) Then

      Do k = 1,wet_model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = q_latest(i,j,k) - q_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &      work_3d,                                                    &
     &      row_length,rows,wet_model_levels,0,0,0,0, at_extremity,     &
     &      stlist(1,stindex(1,item,sect,im_index)),len_stlist,         &
     &      stash_levels,num_stash_levels+1,                            &
     &      atmos_im,sect,item,                                         &
     &      icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 182)"//cmessage
         End if

      Endif  !  sf(item,sect)

!
! Liquid water content increment
!
      item = 183  
      If (icode <= 0 .and. sf(item,sect)) Then

      Do k = 1,wet_model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = qcl_latest(i,j,k) - qcl_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &       work_3d,                                                   &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 183)"//cmessage
         End if

      Endif  !  sf(item,sect)

!
! Total cloud fraction increment
!
      item = 192  
      If (icode <= 0 .and. sf(item,sect)) Then

      Do k = 1,wet_model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = cf_latest(i,j,k) - cf_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           & 
     &       work_3d,                                                   &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 192)"//cmessage
         End if

      Endif  !  sf(item,sect)

!
! Liquid cloud fraction increment
!
      item = 193  
      If (icode <= 0 .and. sf(item,sect)) Then

      Do k = 1,wet_model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = cfl_latest(i,j,k) - cfl_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &       work_3d,                                                   &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        & 
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 193)"//cmessage
         End if

      Endif  !  sf(item,sect)

!
! Surflw
!
      If (sf(201,2)) Then
      
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(201,2,im_index)),surflw,           &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,201,                                            & 
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 201)"
            goto 9999
         End if

      End if

!
! OLR
!
      If (sf(205,2)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(205,2,im_index)),OLR,              &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,205,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 205)"
            goto 9999
         End if

      End if

!
! LWsea : 'net down lw rad flux: open sea'
!
      If (sf(203,2)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(203,2,im_index)),                  &
     &       LWsea,                                                     &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,203,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 203 = LWsea)"
            goto 9999
         End if

      End if

!
! Ozone & ozone Troposphere Diagnostics
!
      If (sf(260,2)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(260,2,im_index)),               &
     &       ozone,                                                     &
     &       row_length,rows,ozone_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,260,2,im_index)),len_stlist,            & 
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,260,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 260)"//cmessage
            goto 9999
         End if

      End if

      If (sf(280,2)) then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(280,2,im_index)),                  &
     &       O3_trop_level,                                             & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,280,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag(item 280)"//cmessage
         End if
      End if

      If (sf(281,2)) then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(281,2,im_index)),                  & 
     &       O3_trop_height,                                            &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,281,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag(item 281)"//cmessage
         End if
      End if

      If (sf(282,2)) then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(282,2,im_index)),                  &
     &       T_trop_level,                                              &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,282,                                            &
     &       icode,cmessage)      

         If (icode  >   0) then
            cmessage=": error in copydiag(item 282)"//cmessage
         End if
      End if

      If (sf(283,2)) then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(283,2,im_index)),                  &
     &       T_trop_height,                                             &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,283,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag(item 283)"//cmessage
         End if
      End if

!
! LW Heating = Lw radiation temp. incr. per timestep / timestep
!
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
     &       work_3d,                                                   &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        & 
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        & 
     &       icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 232)"//cmessage
        End if

      Endif  !  sf(item,sect)


! Now treat diagnostics contained in structure LW_DIAG.

      If (l_out_lwdiag) then
!
! Clear Olr. Stash: sf(206,2)
!
        If (LW_diag(j_lw)%L_clear_olr) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(206+i_off,2,im_index)),           &
     &       LW_diag(j_lw)%clear_olr,                                   &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,206+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage="Error in copydiag( item 206)"
            goto 9999
          End if

        End if

!
! Surface Down Flux. Stash: sf(207,2) 
!
        If (icode <= 0 .and. LW_diag(j_lw)%L_surface_down_flux) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(207+i_off,2,im_index)),           &
     &       LW_diag(j_lw)%surface_down_flux,                           & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,207+i_off,                                   &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage="Error in copydiag( item 207)"
            goto 9999
          End if

        End if

!
! Clear Sky Surface Down Flux. Stash: sf(208,2) 
!
        If (icode <= 0 .and.LW_diag(j_lw)%L_surf_down_clr) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(208+i_off,2,im_index)),           &
     &       LW_diag(j_lw)%surf_down_clr,                               & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,208+i_off,                                   &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage="Error in copydiag( item 208)"
            goto 9999
          End if

        End if

!
! Clear-Sky heating Rates. Stash: sf(233,2)
!
        If (LW_diag(j_lw)%L_clear_hr) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(233+i_off,2,im_index)),        &
     &       LW_diag(j_lw)%clear_hr,                                    & 
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,233+i_off,2,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,233+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 233)"//cmessage
            goto 9999
          End if

        End if

!
! Net LW Flux at Tropopause. Stash: sf(237,2)
!
        If (LW_diag(j_lw)%L_net_flux_trop) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(237+i_off,2,im_index)),           &
     &       LW_diag(j_lw)%net_flux_trop,                               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,237+i_off,                                      &
     &       icode,cmessage) 

          If (icode  >   0) then
            cmessage="Error in copydiag( item 237)"
            goto 9999
          End if

        End if

!
! LW Down Flux at Tropopause. Stash: sf(238,2)
!
        If (LW_diag(j_lw)%L_down_flux_trop) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(238+i_off,2,im_index)),           &
     &       LW_diag(j_lw)%down_flux_trop,                              &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,238+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage="Error in copydiag( item 238)"
            goto 9999
          End if
        End if

!
! Cloud Absorptivity. Stash: sf(262,2)
!
        If (LW_diag(j_lw)%L_cloud_absorptivity) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(262+i_off,2,im_index)),        &
     &       LW_diag(j_lw)%cloud_absorptivity,                          &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,262+i_off,2,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,262+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 262)"//cmessage
            goto 9999
          End if

        End if

!
! Cloud Weight Absorptivity. Stash: sf(263,2)
!
        If (LW_diag(j_lw)%L_cloud_weight_absorptivity) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(263+i_off,2,im_index)),        &
     &       LW_diag(j_lw)%cloud_weight_absorptivity,                   &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,263+i_off,2,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,263+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 263)"//cmessage
            goto 9999
          End if

        End if

!
! LS. Cloud Absorptivity. Stash: sf(264,2)
!
        If (LW_diag(j_lw)%L_ls_cloud_absorptivity) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(264+i_off,2,im_index)),        & 
     &       LW_diag(j_lw)%ls_cloud_absorptivity,                       &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,264+i_off,2,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,264+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 264)"//cmessage
            goto 9999
          End if

        End if

!
! Ls. Cloud Weight Absorptivity. Stash: sf(265,2)
!
        If (LW_diag(j_lw)%L_ls_cloud_weight_absorptivity) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(265+i_off,2,im_index)),        & 
     &       LW_diag(j_lw)%ls_cloud_weight_absorptivity,                &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,265+i_off,2,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,265+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 265)"//cmessage
            goto 9999
          End if

        End if

!
! Convective Cloud Absorptivity. Stash: sf(266,2)
!
        If (LW_diag(j_lw)%L_cnv_cloud_absorptivity) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(266+i_off,2,im_index)),        &
     &       LW_diag(j_lw)%cnv_cloud_absorptivity,                      &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,266+i_off,2,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,266+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 266)"//cmessage
            goto 9999
          End if

        End if

!
! Convective Cloud Weight Absorptivity. Stash: sf(267,2)
!
        If (LW_diag(j_lw)%L_cnv_cloud_weight_absorptivity) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(267+i_off,2,im_index)),        & 
     &       LW_diag(j_lw)%cnv_cloud_weight_absorptivity,               &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,267+i_off,2,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,267+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 267)"//cmessage
            goto 9999
          End if

        End if

!
! Top of the Atmosphere Radiances. Stash: sf(297+i,2)
!
        If (LW_diag(j_lw)%L_toa_radiance) Then

           Do i=1, n_channel

             ptr_stash = si(297+i_off,2,im_index)                      &
     &         + (i-1) * row_length * rows
     
! DEPENDS ON: copydiag     
             Call copydiag (STASHwork(ptr_stash),                      &
     &           LW_diag(j_lw)%toa_radiance(1, 1, i),                  &
     &           row_length,rows,0,0,0,0, at_extremity,                &
     &           atmos_im,2,297+i_off,                                 &
     &           icode,cmessage)

             If (icode  >   0) then
                cmessage="Error in copydiag( item 297)"
                goto 9999
             End if
           Enddo

        End if

!
! Total Cloud Cover. Stash: sf(204,2)
!
      If (LW_diag(j_lw)%L_total_cloud_cover) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(204+i_off,2,im_index)),            &
     &       LW_diag(j_lw)%total_cloud_cover,                           & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,204+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 204)"
            goto 9999
         End if

      End if

!
! Total Cloud on Levels. Stash: sf(261,2)
!
      If (LW_diag(j_lw)%L_total_cloud_on_levels) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(261+i_off,2,im_index)),         &
     &       LW_diag(j_lw)%total_cloud_on_levels,                       &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,261+i_off,2,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,261+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 261)"//cmessage
            goto 9999
         End if

      End if

!
! ISCCP Weights: Stash: sf(269,2)
!
      If (LW_diag(j_lw)%L_isccp_weights) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(269+i_off,2,im_index)),            &
     &       LW_diag(j_lw)%isccp_weights,                               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,269+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 269)"
            goto 9999
         End if

      End if

! Copy ISCCP diagnostics by looping over 7 levels in call to copydiag.
! This is because copydiag_3d cannot handle ISCCP levels.

!
! ISCCP CF: Stash: sf(270,2)
!
      If (LW_diag(j_lw)%L_isccp_cf) Then

        Do k = 1,7
        
! DEPENDS ON: copydiag  
         Call copydiag (STASHwork(si(270+i_off,2,im_index)              &
     &        +(row_length*rows*(k-1))),                                &
     &       LW_diag(j_lw)%isccp_cf(1,1,k),                             &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,270+i_off,                                      &
     &       icode,cmessage)
        Enddo


         If (icode  >   0) then
            cmessage=": error in copydiag( item 270)"//cmessage
            goto 9999
         End if

      End if

!
! isccp_cf_tau_0_to_p3. Stash: sf(271,2)
!
      If (LW_diag(j_lw)%L_isccp_cf_tau_0_to_p3) Then

        Do k = 1,7
        
! DEPENDS ON: copydiag  
         Call copydiag (STASHwork(si(271+i_off,2,im_index)              &
     &        +(row_length*rows*(k-1))),                                & 
     &       LW_diag(j_lw)%isccp_cf_tau_0_to_p3(1,1,k),                 & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,271+i_off,                                      &
     &       icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 271)"//cmessage
            goto 9999
         End if

      End if

!
! Isccp_cf_tau_p3_to_1p3. Stash: sf(272,2)
!
      If (LW_diag(j_lw)%L_isccp_cf_tau_p3_to_1p3) Then

        Do k = 1,7
        
! DEPENDS ON: copydiag  
         Call copydiag (STASHwork(si(272+i_off,2,im_index)              &
     &        +(row_length*rows*(k-1))),                                & 
     &       LW_diag(j_lw)%isccp_cf_tau_p3_to_1p3(1,1,k),               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,272+i_off,                                      & 
     &       icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 272)"//cmessage
            goto 9999
         End if

      End if

!
! Isccp_cf_tau_1p3_to_3p6. Stash: sf(273,2)
!
      If (LW_diag(j_lw)%L_isccp_cf_tau_1p3_to_3p6) Then

        Do k = 1,7
        
! DEPENDS ON: copydiag  
         Call copydiag (STASHwork(si(273+i_off,2,im_index)              &
     &        +(row_length*rows*(k-1))),                                & 
     &       LW_diag(j_lw)%isccp_cf_tau_1p3_to_3p6(1,1,k),              &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,273+i_off,                                      &
     &       icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 273)"//cmessage
            goto 9999
         End if

      End if

!
! Isccp_cf_tau_3p6_to_9p4. Stash: sf(274,2)
!
      If (LW_diag(j_lw)%L_isccp_cf_tau_3p6_to_9p4) Then

        Do k = 1,7
        
! DEPENDS ON: copydiag  
         Call copydiag (STASHwork(si(274+i_off,2,im_index)              &
     &        +(row_length*rows*(k-1))),                                & 
     &       LW_diag(j_lw)%isccp_cf_tau_3p6_to_9p4(1,1,k),              &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,274+i_off,                                      &
     &       icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 274)"//cmessage
            goto 9999
         End if

      End if

!
! Isccp_cf_tau_9p4_to_23. Stash: sf(275,2)
!
      If (LW_diag(j_lw)%L_isccp_cf_tau_9p4_to_23) Then

        Do k = 1,7
        
! DEPENDS ON: copydiag  
         Call copydiag (STASHwork(si(275+i_off,2,im_index)              &
     &        +(row_length*rows*(k-1))),                                & 
     &       LW_diag(j_lw)%isccp_cf_tau_9p4_to_23(1,1,k),               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,275+i_off,                                      &
     &       icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 275)"//cmessage
            goto 9999
         End if

      End if

!
!   isccp_cf_tau_23_to_60. Stash: sf(276,2)
!
      If (LW_diag(j_lw)%L_isccp_cf_tau_23_to_60) Then

        Do k = 1,7
        
! DEPENDS ON: copydiag  
         Call copydiag (STASHwork(si(276+i_off,2,im_index)              &
     &        +(row_length*rows*(k-1))),                                &
     &       LW_diag(j_lw)%isccp_cf_tau_23_to_60(1,1,k),                &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,276+i_off,                                      &
     &       icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 276)"//cmessage
            goto 9999
         End if

      End if

!
! Isccp_cf_tau_ge_60. Stash: sf(277,2)
!
      If (LW_diag(j_lw)%L_isccp_cf_tau_ge_60) Then

        Do k = 1,7
        
! DEPENDS ON: copydiag  
         Call copydiag (STASHwork(si(277+i_off,2,im_index)              &
     &        +(row_length*rows*(k-1))),                                &
     &       LW_diag(j_lw)%isccp_cf_tau_ge_60(1,1,k),                   &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,277+i_off,                                      &
     &       icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 277)"//cmessage
            goto 9999
         End if

      End if

!
! Mean Albedo Cloud. Stash: sf(290,2)
!
      if (LW_diag(j_lw)%L_meanalbedocld) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(290+i_off,2,im_index)),            &
     &        LW_diag(j_lw)%meanalbedocld,                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,290+i_off,                                     &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 290)"
            goto 9999
         End if

      End if

!
! Mean Tau Cloud. Stash: sf(291,2)
!
      if (LW_diag(j_lw)%L_meantaucld) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(291+i_off,2,im_index)),            &
     &       LW_diag(j_lw)%meantaucld,                                  &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,291+i_off,                                      & 
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 291)"
            goto 9999
         End if

      End if

!
! Mean Top: Stash: sf(292,2)
!
      if (LW_diag(j_lw)%L_meanptop) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(292+i_off,2,im_index)),            & 
     &       LW_diag(j_lw)%meanptop,                                    &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,292+i_off,                                      & 
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 292)"
            goto 9999
         End if

      End if

!
! Total Cloud Area
!
      if (LW_diag(j_lw)%L_totalcldarea) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(293+i_off,2,im_index)),            &
     &       LW_diag(j_lw)%totalcldarea,                                &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,293+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 293)"
            goto 9999
         End if

      End if

!   Aerosol optical depth diagnostics
!   (loop on wavelength)

!
! AOD Sulphate. Stash: sf(284,2)
!
      If (LW_diag(j_lw)%L_aod_sulphate) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(284+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(j_lw)%aod_sulphate(1,1,k),                         & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,284+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 284)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Dust. Stash: sf(285,2)
!
      If (LW_diag(j_lw)%L_aod_dust) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(285+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(j_lw)%aod_dust(1,1,k),                             & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,285+i_off,                                      & 
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 285)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Seasalt. Stash: sf(286,2)
!
      If (LW_diag(j_lw)%L_aod_seasalt) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(286+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(j_lw)%aod_seasalt(1,1,k),                          & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,286+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 286)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Soot. Stash: sf(287,2)
!
      If (LW_diag(j_lw)%L_aod_soot) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(287+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(j_lw)%aod_soot(1,1,k),                             &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,287+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 287)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Biomass. Stash: sf(288,2)
!
      If (LW_diag(j_lw)%L_aod_biomass) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(288+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(j_lw)%aod_biomass(1,1,k),                          & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,288+i_off,                                      & 
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 288)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if
      
!
! AOD Biogenic. Stash: sf(289,2)
!
      If (LW_diag(j_lw)%L_aod_biogenic) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(289+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(j_lw)%aod_biogenic(1,1,k),                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     & 
     &       atmos_im,2,289+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 289)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Fossil-fuel organic carbon. Stash: sf(295,2)
!
      If (LW_diag(j_lw)%L_aod_ocff) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(295+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(j_lw)%aod_ocff(1,1,k),                             &
     &       row_length,rows,0,0,0,0, at_extremity,                     & 
     &       atmos_im,2,295+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 295)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Delta aerosol. Stash: sf(296,2)
!
      If (LW_diag(j_lw)%L_aod_delta) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(296+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(j_lw)%aod_delta(1,1,k),                            &
     &       row_length,rows,0,0,0,0, at_extremity,                     & 
     &       atmos_im,2,296+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 296)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if
!
!   Cloud water mixing ratios
!
      item = 308+i_off  ! LS cloud liquid water mixing ratio
      If (LW_diag(j_lw)%L_ls_qcl_rad) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              LW_diag(j_lw)%ls_qcl_rad,                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 308)"//cmessage
            goto 9999
         End if
      Endif

      item = 309+i_off  ! LS cloud ice water mixing ratio
      If (LW_diag(j_lw)%L_ls_qcf_rad) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              LW_diag(j_lw)%ls_qcf_rad,                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 309)"//cmessage
            goto 9999
         End if
      Endif

      item = 310+i_off  ! Convective cloud liquid water mixing ratio
      If (LW_diag(j_lw)%L_cc_qcl_rad) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              LW_diag(j_lw)%cc_qcl_rad,                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 310)"//cmessage
            goto 9999
         End if
      Endif

      item = 311+i_off  ! Convective cloud ice water mixing ratio
      If (LW_diag(j_lw)%L_cc_qcf_rad) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              LW_diag(j_lw)%cc_qcf_rad,                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 311)"//cmessage
            goto 9999
         End if
      Endif

!
!   Cloud amounts
!
      item = 312+i_off  ! LS cloud fraction of grdbox seen by radiation.
                        ! Liquid
      If (LW_diag(j_lw)%L_ls_cl_rad) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              LW_diag(j_lw)%ls_cl_rad,                                  &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 312)"//cmessage
            goto 9999
         End if
      Endif

      item = 313+i_off  ! LS cloud fraction of grdbox seen by radiation.
                        ! Ice
      If (LW_diag(j_lw)%L_ls_cf_rad) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              LW_diag(j_lw)%ls_cf_rad,                                  &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 313)"//cmessage
            goto 9999
         End if
      Endif

      item = 314+i_off  ! CONV cloud fraction of grdbox seen by radiation.
                        ! Liquid
      If (LW_diag(j_lw)%L_cc_cl_rad) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              LW_diag(j_lw)%cc_cl_rad,                                  &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 314)"//cmessage
            goto 9999
         End if
      Endif

      item = 315+i_off  ! CONV cloud fraction of grdbox seen by radiation.
                        ! Ice
      If (LW_diag(j_lw)%L_cc_cf_rad) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              LW_diag(j_lw)%cc_cf_rad,                                  &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 315)"//cmessage
            goto 9999
         End if
      Endif

      End If  ! l_out_lwdiag


 9999 continue  ! exit point on error
      If(icode /= 0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,Cmessage)
      Endif

      Return
    End Subroutine diagnostics_lw
#endif
