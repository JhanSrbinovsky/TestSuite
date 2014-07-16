#if defined(A70_1C) || defined(A70_1Z)
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT******************************

! Subroutine diagnostics_sw
      Subroutine diagnostics_sw(                                        &
     &      row_length, rows, model_levels                              &
     &,     wet_model_levels, cloud_levels                              &
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
     &,     surfsw, itoasw                                              &
     &,     SWsea, flux_below_690nm_surf                                &
     &,     photosynth_act_rad                                          &
     &,     f_orog                                                      &
     &,     slope_aspect, slope_angle                                   &
     &,     sw_net_land,sw_net_sice                                     &
     &,     T_incr_diagnostic                                           &
     &,     n_channel                                                   &
     &,     sea_salt_film, sea_salt_jet                                 &
     &,     salt_dim1, salt_dim2, salt_dim3                             &
     &,     l_out_swdiag, j_sw                                          &
     &,                                                                 &
#include "argsts.h"
     &      STASHwork )
     
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
! 6.2  13/06/05 Add orography correction diagnostics. J. Manners
! 6.2  19/01/06 Add radiative forcing, radiance, UV-flux diagnostics.
!                                                J.C. Thelen
! 6.2  02/03/06 Added diagnostics for total and direct component
!               of surface PAR flux.             M.G. Sanderson
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Use sw_diag_mod, Only: SW_diag

      Implicit None
      
!
! Arguments with Intent IN. ie: Input variables.
!

      Logical :: at_extremity(4)  
!       Indicates if this processor is at north,
!       south, east or west of the processor grid
      Logical :: l_out_swdiag
!        Flag to output radiation diagnostics (mainly
!         required for time-stepping scheme.

      Integer :: row_length       ! number of points on a row
      Integer :: rows             ! number of rows in a theta field
      Integer :: n_rows           ! number of rows in a v field
      Integer :: model_levels     ! number of model levels
      Integer :: cloud_levels     ! number of cloudy levels
      Integer :: wet_model_levels ! number of model levels where moisture
                                  ! variables are held
      Integer :: number_format    ! switch controlling number format diagnostics
                                  ! are written out in. See PP_WRITE for details.
      Integer :: model_domain     ! indicator as to model type, ie global, lam
      Integer :: n_channel        ! Number of satellite channels used
      Integer :: salt_dim1        !
      Integer :: salt_dim2        ! Dimensions for sea-salt aerosol diagnostics.
      Integer :: salt_dim3        !

      
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

      Real :: timestep

      Integer, Intent(IN) :: i_off
!       Offset to diagnostics in multiple calls to radiation
      Integer, Intent(IN) :: j_sw
!       Call to SW radiation

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

      Real :: itoasw (row_length, rows)
      Real :: surfsw (row_length, rows)
      Real :: SWsea(row_length, rows)   ! Net short-wave absorbed by planet
      Real :: flux_below_690nm_surf(row_length, rows)
      Real :: sw_net_land(row_length, rows)
      Real :: sw_net_sice(row_length, rows)
      Real :: photosynth_act_rad(row_length, rows)
      Real :: T_incr_diagnostic(row_length,rows,model_levels)
      Real :: sea_salt_film(salt_dim1, salt_dim2, salt_dim3)
      Real :: sea_salt_jet(salt_dim1, salt_dim2, salt_dim3)
!
! Orography variables
!
      Real :: f_orog(row_length,rows)      ! Extra SW surf flux
      Real :: slope_aspect(row_length,rows)! Gridbox mean slope aspect
      Real :: slope_angle(row_length,rows) ! Gridbox mean slope angle

#include "csubmodl.h"
#include "typsts.h"

!
! Diagnostic variables
!
      Real ::  STASHwork(*)    ! STASH workspace
!
! Local array & variables
!
      
      Real :: T_plus_T_inc(row_length, rows, model_levels)
      Real :: heating_rate(row_length,rows,model_levels)
      Real :: work_3d(row_length,rows,model_levels)

      Integer :: i, j, k
      Integer :: icode           ! Return code  =0 Normal exit  >1 Error
      Integer :: ptr_stash       ! Pointer to position in STASH

      CHARACTER(LEN=1024)  :: cmessage
      Character(*) RoutineName
      Parameter ( RoutineName='diagnostics_sw')

      Integer :: im_index                ! internal model index
      Integer :: item                    ! STASH item
      Integer, Parameter :: sect=1       ! STASH SW radiation section


      icode    = 0                       ! Initialise error status
      im_index = internal_model_index(atmos_im)

! ----------------------------------------------------------------------
! Section 1.  Diagnostic Calculation and output.
! ----------------------------------------------------------------------

!
! First treat diagnostics that are not contained in the SW_Diag Structure.
!

!
! Temperature
!
     If (sf(004,1)) Then

         Do k = 1, model_levels
            Do j = 1, rows
               Do i = 1, row_length
                  T_plus_T_inc(i,j,k) = T_n(i,j,k) + T_inc(i,j,k)
               End Do
            End Do
         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(004,1,im_index)),                &
     &       T_plus_T_inc,                                              &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,004,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,004,                                            &
     &       icode,cmessage)

         If (icode >  0) Then
            cmessage="Error in copydiag_3d( item 004)"
            goto 9999
         End if

      End if

!
! Temperature Increment: increment diagnostics= modified - previous
!
      item = 161
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
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            & 
     &       work_3d,                                                   &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        & 
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)

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
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &       work_3d,                                                   &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)

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
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
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
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
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
! Surfsw
!
      If (sf(201,1)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(201,1,im_index)),surfsw,           &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,201,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 201)"
            goto 9999
         End if

      End if

!
! Incoming SW at TOA
!
      If (sf(207,1)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(207,1,im_index)),itoasw,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,207,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 207)"
            goto 9999
         End if

      End if

!
! SWSea
!
      If (sf(203,1)) Then

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

!
! Flux Below 690nm at Surface
!
      If (sf(204,1)) Then

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
! Surface Net Land
!
      If (sf(257,1)) Then

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
! SW Net Sice
!
      If (sf(258,1)) Then

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
! SW Heating:  SW heating =
! sw radiation temperature increment per timestep / timestep
!
      If (icode <= 0 .and. sf(232,1)) Then

        Do k=1,model_levels
         Do j=1,rows
          Do i=1,row_length
            heating_rate(i,j,k) =  T_incr_diagnostic(i,j,k) /           &
     &                                       timestep
          Enddo ! i
         Enddo ! j
        Enddo ! k

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(232,1,im_index)),                &
     &       heating_rate,                                              &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,232,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,232,                                            &
     &       icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 232)"//cmessage
        End if

      Endif

!
! Photosynth_act_rad (total PAR flux at surface)
! 
      If (sf(290,1)) Then

! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(290,1,im_index)),                  &
     &       photosynth_act_rad,                                        &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,290,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
             cmessage="Error in copydiag( item 290)"
             goto 9999
          End if

      End if

!
! Slope Aspect
!      
      If (sf(293,1)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(293,1,im_index)),                   &
     &       slope_aspect,                                              &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,293,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 293)"
            goto 9999
         End if

      End if

!
! Slope Angle
!
      If (sf(294,1)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(294,1,im_index)),                   &
     &       slope_angle,                                               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,294,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 294)"
            goto 9999
         End if
      End if

!      
! F_orog
!
      If (sf(296,1)) Then
      
! DEPENDS ON: copydiag      
         Call copydiag(STASHwork(si(296,1,im_index)),                   &
     &       f_orog,                                                    &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,296,                                            & 
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 296)"
            goto 9999
         End if

      End if

!
! Sea Salt film
!      
      If (sf(247+i_off,1)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(247+i_off,1,im_index)),          &
     &       sea_salt_film,                                             &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,247,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,247+i_off,                                      & 
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 247)"
            goto 9999
         End if

      End if

!
! Sea Salt Jet
!
      If (sf(248+i_off,1)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(248+i_off,1,im_index)),          &
     &       sea_salt_jet,                                              &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,248+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,248+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 248)"
            goto 9999
         End if

      End if      


! Now treat diagnostics contained in structure SW_DIAG.

      If (l_out_swdiag) then     

!
! Outwards Solar Flux at TOA. Stash: Sf(208,1)
!           
        If (SW_diag(j_sw)%L_solar_out_toa) Then

! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(208+i_off,1,im_index)),            &
     &       SW_diag(j_sw)%solar_out_toa,                               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,208+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage="Error in copydiag( item 208)"
            goto 9999
          End if

        End if

!
! Outwards Solar Clear Flux at TOA. Stash: Sf(209,1)
! 
        If (SW_diag(j_sw)%L_solar_out_clear) Then

! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(209+i_off,1,im_index)),            &
     &       SW_diag(j_sw)%solar_out_clear,                             &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,209+i_off,                                      &
     &       icode,cmessage)                     

           If (icode  >   0) then
             cmessage="Error in copydiag( item 209)"
             goto 9999
           End if

        End if

!
! Surface Down Clear: Stash: sf(210,1)
!
        If (SW_diag(j_sw)%L_surf_down_clr) Then

! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(210+i_off,1,im_index)),            &
     &       SW_diag(j_sw)%surf_down_clr,                               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,210+i_off,                                      &
     &       icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 210)"
             goto 9999
           End if

        End if

!
! Surface up Clear. Stash: sf(211,1)
!
        If (SW_diag(j_sw)%L_surf_up_clr) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(211+i_off,1,im_index)),           &
     &       SW_diag(j_sw)%surf_up_clr,                                 &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,211+i_off,                                      &
     &       icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 211)"
             goto 9999
           End if

        End if

!
! Surface Down Flux. Stash: sf(235,1)
!
        If (SW_diag(j_sw)%L_surface_down_flux) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(235+i_off,1,im_index)),           &
     &       SW_diag(j_sw)%surface_down_flux,                           &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,235+i_off,                                      &
     &       icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 235)"
             goto 9999
           End if

        End if

!
! Direct UV-Flux. Stash: sf(212,1)
!
      If (SW_diag(j_sw)%L_uvflux_direct) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(212+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%uvflux_direct,                               &
     &       row_length,rows,model_levels+1,0,0,0,0, at_extremity,      &
     &       stlist(1,stindex(1,212+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,212+i_off,                                      &
     &       icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 212)"
            goto 9999
         End if
      endif
!
! UV Flux Up. Stash: sf(213,1)
!
      If (SW_diag(j_sw)%L_uvflux_up) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(213+i_off,1,im_index)),          & 
     &       SW_diag(j_sw)%uvflux_up,                                   &
     &       row_length,rows,model_levels+1,0,0,0,0, at_extremity,      &
     &       stlist(1,stindex(1,213+i_off,1,im_index)),len_stlist,      & 
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,213+i_off,                                      &
     &       icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 213)"
            goto 9999
         End if

      End if
!
! Net UV Flux. Stash: sf(214,1)
!
      If (SW_diag(j_sw)%L_uvflux_net) then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(214+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%uvflux_net,                                  &
     &       row_length,rows,model_levels+1,0,0,0,0, at_extremity,      & 
     &       stlist(1,stindex(1,214+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,214+i_off,                                      &
     &       icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 214)"
            goto 9999
         End if

      End if

!
! Direct Flux. Stash: sf(230,1)
!
      If (SW_diag(j_sw)%L_flux_direct) Then
      
         Call copydiag_3d(STASHwork(si(230+i_off,1,im_index)),          &
     &        SW_diag(j_sw)%flux_direct,                                &
     &        row_length,rows,model_levels+1,0,0,0,0, at_extremity,     & 
     &        stlist(1,stindex(1,230+i_off,1,im_index)),len_stlist,     &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,230+i_off,                                     &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 230)"
            goto 9999
         End if

      End if
!
! Diffuse Flux. Stash: sf(231,1)
!
      If (SW_diag(j_sw)%L_flux_diffuse) Then

         Call copydiag_3d(STASHwork(si(231+i_off,1,im_index)),          &
     &        SW_diag(j_sw)%flux_diffuse,                               &
     &        row_length,rows,model_levels+1,0,0,0,0, at_extremity,     &
     &        stlist(1,stindex(1,231+i_off,1,im_index)),len_stlist,     &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,231+i_off,                                     &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 231)"
            goto 9999
         End if

      End if    
        
!
! Clear Sky Heating Rates. Stash: sf(233,1)
!      
        If (SW_diag(j_sw)%L_clear_hr) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d(STASHwork(si(233+i_off,1,im_index)),         &
     &       SW_diag(j_sw)%clear_hr,                                    & 
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,233+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,233+i_off,                                      &
             icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 233)"
             goto 9999
           End if

        End if

!
! Net SW Flux at Tropopause: Stash: sf(237,1)
!
        If (SW_diag(j_sw)%L_net_flux_trop) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(237+i_off,1,im_index)),           &
     &       SW_diag(j_sw)%net_flux_trop,                               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,237+i_off,                                      &
     &       icode,cmessage)

           If (icode  >   0) then
              cmessage="Error in copydiag( item 237)"
              goto 9999
           End if

        End if

!
! SW Up Flux at Tropopause: Stash: sf(238,i)
!
        If (sf(238+i_off,1)) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(238+i_off,1,im_index)),           &
     &       SW_diag(j_sw)%up_flux_trop,                                &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,238+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage="Error in copydiag( item 238)"
            goto 9999
          End if

        End if

!
! Cloud Extinction. Stash: sf(262,1)
!
        If (SW_diag(j_sw)%L_cloud_extinction) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d(STASHwork(si(262+i_off,1,im_index)),         &
     &       SW_diag(j_sw)%cloud_extinction,                            &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,262+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,262+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage="Error in copydiag( item 262)"
            goto 9999
          End if

        End if

!
! Cloud Weight Extinction. Stash: sf(263,1)
!
        If (SW_diag(j_sw)%L_cloud_weight_extinction) Then
 
 ! DEPENDS ON: copydiag_3d
          Call copydiag_3d(STASHwork(si(263+i_off,1,im_index)),         & 
     &       SW_diag(j_sw)%cloud_weight_extinction,                     &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,263+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,263+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
             cmessage="Error in copydiag( item 263)"
             goto 9999
          End if

        End if

!
! Large-Scale Cloud Extinction. Stash: sf(264,1)
!
        If (SW_diag(j_sw)%L_ls_cloud_extinction) Then
 
! DEPENDS ON: copydiag_3d
           Call copydiag_3d(STASHwork(si(264+i_off,1,im_index)),        &
     &       SW_diag(j_sw)%ls_cloud_extinction,                         &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,264+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,264+i_off,                                      &
     &       icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 264)"
             goto 9999
           End if
      
        End if

!
! Large-Scale Cloud Weight Extinction. Stash: sf(265,1)
!
        If (SW_diag(j_sw)%L_ls_cloud_weight_extinction) Then

! DEPENDS ON: copydiag_3d
           Call copydiag_3d(STASHwork(si(265+i_off,1,im_index)),        &
     &       SW_diag(j_sw)%ls_cloud_weight_extinction,                  &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,265+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,265+i_off,                                      &
     &       icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 265)"
             goto 9999
           End if

        End if

!
! Convective Cloud Extinction. Stash: sf(266,1)
!
        If (SW_diag(j_sw)%L_cnv_cloud_extinction) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d(STASHwork(si(266+i_off,1,im_index)),         &
     &       SW_diag(j_sw)%cnv_cloud_extinction,                        &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,266+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,266+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
             cmessage="Error in copydiag( item 266)"
             goto 9999
          End if

        End if

!
! Convective Cloud weight Extinction. Stash: sf(267,1)
!
        If (SW_diag(j_sw)%L_cnv_cloud_weight_extinction) Then

! DEPENDS ON: copydiag_3d
           Call copydiag_3d(STASHwork(si(267+i_off,1,im_index)),        &
     &       SW_diag(j_sw)%cnv_cloud_weight_extinction,                 &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,267+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,267+i_off,                                      &
     &       icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 267)"
             goto 9999
           End if

        End if

!
! Radiances at TOA: Stash: sf(297,1)
!
        If (SW_diag(j_sw)%L_toa_radiance) Then

           Do i=1, n_channel

             ptr_stash = si(297+i_off,1,im_index)                       &
     &                 + (i-1) * row_length * rows
     
! DEPENDS ON: copydiag     
             Call copydiag (STASHwork(ptr_stash),                       &
     &           SW_diag(j_sw)%toa_radiance(1, 1, i),                   &
     &           row_length,rows,0,0,0,0, at_extremity,                 &
     &           atmos_im,1,297+i_off,                                  &
     &           icode,cmessage)

             If (icode  >   0) then
                cmessage="Error in copydiag( item 297)"
                goto 9999
             End if
           End do

        End if

!
! Fl_solid_below_690nm_surf. Stash: sf(259,1)
!
      If (SW_diag(j_sw)%L_FlxSolBelow690nmSurf) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(259+i_off,1,im_index)),             &
     &       SW_diag(j_sw)%FlxSolBelow690nmSurf,                        &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,259+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 259)"
            goto 9999
         End if

      End if

!
! Fl_sea_below_690nm_surf. Stash: sf(260,i)
! 
      If (SW_diag(j_sw)%L_FlxSeaBelow690nmSurf) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(260+i_off,1,im_index)),             &
     &       SW_diag(j_sw)%FlxSeaBelow690nmSurf,                        &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,260+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 260)"
            goto 9999
         End if

      End if

!
! Re. Strat. Stash: sf(221,1)
!
      If (SW_diag(j_sw)%re_strat_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(221+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%re_strat,                                    &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,221+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,221+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 221)"
            goto 9999
         End if

      End if

!
! Wgt. Strat. Stash: sf(223,1)
!
      If (SW_diag(j_sw)%wgt_strat_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(223+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%wgt_strat,                                   &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,223+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,223+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 223)"
            goto 9999
         End if

      End if

!
! LWP. Strat. Stash: sf(224,i)
!
      If (SW_diag(j_sw)%lwp_strat_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(224+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%lwp_strat,                                   &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,224+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,224+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 224)"
            goto 9999
         End if

      End if

!
! Re. Conv. Stash: sf(225,i)
!
      If (SW_diag(j_sw)%re_conv_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(225+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%re_conv,                                     &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,225+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           & 
     &       atmos_im,1,225+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 225)"
            goto 9999
         End if

      End if

!
! Wgt. Conv. Stash: sf(226,1)
!
      If (SW_diag(j_sw)%wgt_conv_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(226+i_off,1,im_index)),          & 
     &       SW_diag(j_sw)%wgt_conv,                                    &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,226+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,226+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 226)"
            goto 9999
         End if

      End if

!
! Ntot. Diag. Stash: sf(241,1)
!
      If (SW_diag(j_sw)%ntot_diag_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(241+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%ntot_diag,                                   &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,241+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,241+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 241)"
            goto 9999
         End if

      End if

!
! Strat. LWC Diag. Stash: sf(242,1)
!
      If (SW_diag(j_sw)%strat_lwc_diag_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(242+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%strat_lwc_diag,                              &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,242+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,242+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 242)"
            goto 9999
         End if

      End if

!
! SO4 Cloud Condensation Nuclei. Stash: sf(243,1)
!
      If (SW_diag(j_sw)%so4_ccn_diag_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(243+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%so4_ccn_diag,                                & 
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,243+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,243+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 243)"
            goto 9999
         End if

      End if

!
! Cond. Samp. Wgt. Stash(244,1)
!
      If (SW_diag(j_sw)%cond_samp_wgt_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(244+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%cond_samp_wgt,                               &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,244+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,244+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 244)"
            goto 9999
         End if

      End if

!
!Weighted Re. Stash: sf(245,1)
! 
      If (SW_diag(j_sw)%weighted_re_flag) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(245+i_off,1,im_index)),            & 
     &       SW_diag(j_sw)%weighted_re,                                 &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,245+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 245)"
            goto 9999
         End if

      End if

!
! Sum Weighted Re. Stash: sf(246,1)
!
      If (SW_diag(j_sw)%sum_weight_re_flag) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(246+i_off,1,im_index)),            &
     &       SW_diag(j_sw)%sum_weight_re,                               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,246+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 246)"
            goto 9999
         End if

      End if

!
! Flux_direct_par (direct component of PAR flux at surface)
!
      if (SW_diag(j_sw)%L_direct_par) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(291+i_off,1,im_index)),             &
     &       SW_diag(j_sw)%flxdirparsurf,                               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,291+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 291)"
            goto 9999
         End if

      End if

!
! Weighted Warm Re. Stash: sf(254,1)
!
      If (SW_diag(j_sw)%wgtd_warm_re_flag) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(254+i_off,1,im_index)),            &
     &       SW_diag(j_sw)%weighted_warm_re,                            &
     &       row_length,rows,0,0,0,0, at_extremity,                     & 
     &       atmos_im,1,254+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 254)"
            goto 9999
         End if

      End if

!
! Sum Weighted Warm Re. Stash: sf(255,1)
!
      If (SW_diag(j_sw)%sum_wgt_warm_re_flag) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(255+i_off,1,im_index)),            &
     &       SW_diag(j_sw)%sum_weight_warm_re,                          & 
     &       row_length,rows,0,0,0,0, at_extremity,                     & 
     &       atmos_im,1,255+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 255)"
            goto 9999
         End if

      End if

!
! Nc. Diag. Stash: sf(280,1)
!
      If (SW_diag(j_sw)%Nc_diag_flag) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(280+i_off,1,im_index)),            &
     &       SW_diag(j_sw)%Nc_diag,                                     &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,280+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 280)"
            goto 9999
         End if

      End if

!
! Nc. Weight. Stash: sf(281,1)
!
      If (SW_diag(j_sw)%Nc_weight_flag) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(281+i_off,1,im_index)),            &
     &       SW_diag(j_sw)%Nc_weight,                                   &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,281+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 281)"
            goto 9999
         End if

      End if

!
! Solar Bearing. Stash: sf(292,1)
!
      If (SW_diag(j_sw)%L_sol_bearing) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(292+i_off,1,im_index)),             &
     &       SW_diag(j_sw)%sol_bearing,                                 &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,292+i_off,                                      &
     &       icode,cmessage)              

         If (icode  >   0) then
            cmessage="Error in copydiag( item 292)"
            goto 9999
         End if

      End if

!
! Orographic Correction. Stash: sf(295,1)
!
      If (SW_diag(j_sw)%L_orog_corr) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(295+i_off,1,im_index)),             &
     &       SW_diag(j_sw)%orog_corr,                                   & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,295+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 295)"
            goto 9999
         End if

      End if

      End if  ! l_out_swdiag


 9999 continue  ! exit point on error
      If(icode /= 0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,Cmessage)
      Endif

      Return
    End Subroutine diagnostics_sw
#endif
