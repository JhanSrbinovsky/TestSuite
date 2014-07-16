#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Diagnostic output area for Large Scale Cloud (Section 9) Diagnostics.
! Subroutine Interface:
      Subroutine diagnostics_lscld(                                     &
     &                       row_length, rows, model_levels             &
     &,                      rhc_row_length, rhc_rows                   &
     &,                      wet_model_levels, boundary_layer_levels    &
     &,                      cloud_levels                               &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      at_extremity, p_theta_levels               &
     &,                      p, r_theta_levels, r_rho_levels            &
     &,                      T_incr_diagnostic, q_incr_diagnostic       &
     &,                      qcl_incr_diagnostic                        &
     &,                      T, q, qcl, qcf                             &
     &,                      area_cloud_fraction, bulk_cloud_fraction   &
     &,                      p_star, rhcpt                              &
     &,                      combined_cloud, cca, ccb, cct              &
     &,                      n_cca_levels, L_murk, Aerosol, RHcrit      &
     &,                      l_mixing_ratio                             &
     &,                                                                 &
#include "argsts.h"
     & STASHwork                                                        &
     &  )
!
Use ac_diagnostics_mod, Only :  &
    cf_lsc, qcl_lsc

      Implicit None
!
! Purpose:
!          Calculates diagnostics associated with large scale cloud
! (section 9) for output through STASH system. This routine is
! positioned at the end of moist processes calculation, ie the end
! of the moist timestep, and is therefore the appropriate location
! for output of moisture related diagnostics.
!
! Method:
!  Required level lists and logical switches are determined by the
! calling routine from STASH requests and STASHflags.
! Intercepted arrays - calculated unconditionally - and diagnostic
! arrays - dependent on STASHflags - are input from the calling
! routine. After minor calculations, each diagnostic is simply
! copied into the STASHwork array to be passed on to STASH for
! output processing.
!
!  Diagnostics currently available: (in order calculated)
!
! STASH item (all section 9 )
!   4 Temperature on model levels after large scale cloud
! 181 T   increment over boundary layer + large scale cloud routines
! 182 q   increment over boundary layer + large scale cloud routines
! 183 qcl increment over boundary layer + large scale cloud routines
! 201 bulk cloud fraction
! 202 very low cloud amount
! 203 low      cloud amount
! 204 medium   cloud amount
! 205 high     cloud amount
! 206 qcl after large scale cloud. Note this is only needed to provide
!     assimilation increments for latent heat nudging. Otherwise should
!     use prognostic (0,254).
! 208 cloud base for cover >  0.1 octa kft
! 209 cloud base for cover >  1.5 octa kft
! 210 cloud base for cover >  2.5 octa kft
! 211 cloud base for cover >  3.5 octa kft
! 212 cloud base for cover >  4.5 octa kft
! 213 cloud base for cover >  5.5 octa kft
! 214 cloud base for cover >  6.5 octa kft
! 215 cloud base for cover >  7.9 octa kft
! 216 total cloud ; random overlap
! 217 total cloud ; max/random overlap
! 218 cloud fraction below 1000 ft asl
! 219 low cloud base ft asl
! 220 low cloud top ft asl
! 221 wet bulb freezing level
! 222 wet bulb temperature
! 223 total cloud top height kft
! 229 relative humidity on model levels. Note: percentage units.
! 226 layer cloud frequency
! 228 critical relative humidity on model levels
! 230 Visibility on model levels. Note that if murk aerosol is included
!     then vis depends on aerosol as well as humidity. Visibility is
!     available over all model levels but aerosol is only actively
!     mixed over boundary levels.
! 231 combined cloud amount on model levels
!
! Current Owner of Code: Owner of Section 9 (Large-Scale Cloud).
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
!  Global Variables:----------------------------------------------------
!     None.
!
! Arguments with Intent IN. ie: Input variables.

      LOGICAL ,INTENT(IN) ::                                            &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

      INTEGER ,INTENT(IN) ::                                            &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, rhc_row_length                                                  &
                         ! row_length if diagnostic RHcrit ON, else 1.
     &, rhc_rows                                                        &
                         ! rows       if diagnostic RHcrit ON, else 1.
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels
     &, cloud_levels                                                    &
                         ! number of cloud model levels
     &, wet_model_levels                                                &
                         ! number of model levels where moisture
                         ! variables are held
     &, boundary_layer_levels                                           &
                              ! number of boundary layer levels
     &, n_cca_levels     ! Number of levels for conv cloud
                         ! amount: 1 for 2D, nlevs for 3D.

      INTEGER ,INTENT(IN) ::                                            &
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

      LOGICAL ,INTENT(IN) ::                                            &
     &  L_murk                                                          &
                            ! murk aerosol included if T
     &, L_mixing_ratio      ! Use mixing ratio formulation

! Primary Arrays used in all models
      Real ,INTENT(IN) ::                                               &
     &  p_theta_levels(row_length, rows, model_levels)                  &
     &, p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, model_levels)              &
     &, T(row_length, rows, model_levels)                               &
     &, q(row_length,rows, wet_model_levels)                            &
     &, qcl(row_length, rows, wet_model_levels)                         &
     &, qcf(row_length, rows, wet_model_levels)                         &
     &, area_cloud_fraction(row_length, rows, wet_model_levels)         &
     &, bulk_cloud_fraction(row_length, rows, wet_model_levels)         &
     &, p_star(row_length, rows)                                        &
     &, RHCPT(rhc_row_length, rhc_rows, wet_model_levels)               &
     &, cca(row_length, rows, n_cca_levels)                             &
     &, T_incr_diagnostic(row_length, rows, model_levels)               &
     &, q_incr_diagnostic(row_length,rows, wet_model_levels)            &
     &, qcl_incr_diagnostic(row_length, rows, wet_model_levels)         &
     &, rhcrit(wet_model_levels)                                        &
                                    ! critical RH for cloud formation
     &, Aerosol(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          model_levels)                                           &
                                    ! murk aerosol
     &, combined_cloud(row_length, rows, wet_model_levels)
!                                        Mixed CCA and CF per gridbox.
!
      INTEGER  ,INTENT(IN) ::                                           &
     &  ccb(row_length, rows)                                           &
     &, cct(row_length, rows)

#include "csubmodl.h"
#include "typsts.h"
#include "cmaxsize.h"
#include "cconsts.h"
!         cconsts defines LOW_BOT_LEVEL to HIGH_TOP_LEVEL cloud splits.
#include "clschm3a.h"
!         clschm3a defines reference numbers for cloud overlap choices.
#include "c_0_dg_c.h"
#include "c_kt_ft.h"
#include "c_lowcld.h"
#include "c_mdi.h"
#include "c_a.h"
!
! arguments with intent in/out. ie: input variables changed on output.
!
! Diagnostics info
      Real                                                              &
     & STASHwork(*)     ! STASH workspace
!
!  Local parameters and other physical constants------------------------
!
      Integer sect ! STASH section for diagnostics
      Parameter ( sect = 9 )
!      ( LS Cloud diagnostics are output under STASH Section 09 )
!
      Character(*) RoutineName
      Parameter  ( RoutineName='diagnostics_lscld')
!
!  Local scalars--------------------------------------------------------
!
      Integer                                                           &
     &  i, j, k, l, ji                                                  &
     &, kinvert                                                         &
                                ! vertical index for inverted arrays.
     &, item                                                            &
                                ! STASH item of individual diagnostics.
     &,    icode                ! Return code  =0 Normal exit  >1 Error

      REAL, PARAMETER:: vis_probability=0.5 ! median visibility diag

      Character*80  cmessage
!
      Integer                                                           &
     &  im_index                                                        &
                        ! internal model index
!                  ( LS Cloud is part of the Atmosphere internal model )
     &, nclds      ! Number of radiation cloud levels ( <=  wet levels)
!  Local scalars required to calculate Diagnostics 09208 - 09215
      INTEGER, PARAMETER :: noktas=8
      REAL,    PARAMETER :: roktas=8.0
      REAL    :: cloud_base_height         ! Height of cloud base (m)
      REAL    :: combined_cloud_in_oktas   ! Combined cloud in oktas

!  Local dynamic arrays required to calculate Diagnostics 09208 - 09215
      INTEGER ::                                                        &
     & cloud_base_level(row_length, rows, noktas)
      REAL ::                                                           &
     & cloud_base(row_length, rows, noktas)
      REAL :: oktas(noktas)
      DATA oktas/0.1,1.5,2.5,3.5,4.5,5.5,6.5,7.9/

!  Local scalars required to calculate Diagnostics 09218, 09219, 09220
      REAL ::                                                           &
     & pu,pl                                                            &
                                       ! Upper and lower half
                                       !       level pressure
     &,pt                                                               &
                                       ! Pressure thickness accumulator
     &,ft                                                               &
                                       ! Cloud fraction accumulator
     &,dp                                                               &
                                       ! Later pressure thickness
     &,h_asl                                                            &
                                       ! Layer base height asl
     &,h_asln                                                           &
                                       ! Layer top height asl
     &,fr                              ! Layer fraction below ceiling

      REAL, PARAMETER ::                                                &
     & str_ceilm = str_ceil * ft2m     ! STR_CEIL in metres

!  Local scalars required to calculate Diagnostic 09221 and 09222
      LOGICAL, PARAMETER ::                                             &
     & l_potential=.false.             ! Wet bulb temperature required
                                       ! from subroutine ThetaW
      REAL ::                                                           &
     & frac

!  Local scalars required to calculate Diagnostic 09223
      INTEGER ::                                                        &
     & cld_top_lev                     ! Top layer in cloud
      REAL, PARAMETER::                                                 &
     & thresh = 0.0627,                                                 &
     & m_to_kft = (1./ft2m)*0.001
!
!  Local dynamic arrays-------------------------------------------------
       Real                                                             &
     &  work2d_1(row_length, rows)                                      &
                                       ! Single-level work array (cloud)
     &, work2d_2(row_length, rows)                                      &
     &, work2d_3(row_length, rows)                                      &
     &, work3d_1(row_length, rows, wet_model_levels)
!                                        Full work array (cloud/moist).
!
!  External subroutine calls: ------------------------------------------
      External                                                          &
     &  R2_calc_total_cloud_cover                                       &
     &, QSAT                                                            &
     &, copydiag, copydiag_3d                                           &
     &, Ereport                                                         &
     &, Thetaw
!- End of Header
!
! ==Main Block==--------------------------------------------------------
      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)
!
! ----------------------------------------------------------------------
! DIAG.09004 Copy T from Main Cloud to stashwork
! ----------------------------------------------------------------------
      item = 4
! Diag09004_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),T,         &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
         If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(T Main Cloud)"
         End if
      End if  ! Diag09004_if1
!
! ----------------------------------------------------------------------
! DIAG.09010 Copy q to stashwork
! ----------------------------------------------------------------------
      If (.false.) Then
! Prevent access to diagnostic which duplicates 00010 and thus redundant
! Pending complete removal of code at UM Version 5.3.
      item = 10
! Diag09010_if1:
!     If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),q,         &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
!#include <argppx/argppx.h>
     &      icode,cmessage)

         If (icode >  0) Then
            cmessage="cld_ctl  : error in copydiag_3d(q)"
         End if
      End if  ! Diag09010_if1
!
! ----------------------------------------------------------------------
! DIAG.09181 Copy T INC: bdy layer + ls cld to stashwork
! ----------------------------------------------------------------------
      item = 181
! Diag09181_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &       T_incr_diagnostic,                                         &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(T INC: bl + ls cld)"
        End if
      End if  ! Diag09181_if1
!
! ----------------------------------------------------------------------
! DIAG.09182 Copy q INC: bdy layer + ls cld to stashwork
! ----------------------------------------------------------------------
      item = 182
! Diag09182_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &       q_incr_diagnostic,                                         &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(Q INC: bl + ls cld)"
        End if
      End if  ! Diag09182_if1
!
! ----------------------------------------------------------------------
! DIAG.09183 Copy qCL INC: bdy layer + ls cld to stashwork
! ----------------------------------------------------------------------
      item = 183
! Diag09183_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
! DEPENDS ON: copydiag_3d
        Call copydiag_3d ( stashwork(si(item,sect,im_index)),           &
     &       qcl_incr_diagnostic,                                       &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
        cmessage="cld_ctl  : error in copydiag_3d(Qcl INC: bl + ls cld)"
        End if
      End if  ! Diag09183_if1
!
! ----------------------------------------------------------------------
! DIAG.09201 Copy Bulk cloud fraction to stashwork
! ----------------------------------------------------------------------
! Needed for Assimilation, otherwise duplicates 00266 output (preferred)
      item = 201
! Diag09201_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then

      If (.not.allocated(cf_lsc)) then
        Allocate ( cf_lsc(row_length*rows,wet_model_levels) )
      End If
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        bulk_cloud_fraction,                                      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
! Copy bulk_cloud_fraction into cf_lsc for ac_diagnostics_mod

      do k = 1,wet_model_levels
        do j = 1,rows
          do i = 1,row_length
            ji = (j-1)*row_length+i
            cf_lsc(ji,k) = bulk_cloud_fraction(i,j,k)
          end do
        end do
      end do

        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(Bulk CF Main Cloud)"
        End if
      End if  ! Diag09201_if1
!
! ----------------------------------------------------------------------
! DIAG.09202 Copy Very Low Cloud Amount to stashwork
! ----------------------------------------------------------------------
      item = 202
! Diag09202_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Initialize cloud amount to lowest level value.
! Diag09202_do1:
        Do j = 1, rows
          Do i = 1, row_length
            work2d_1(i, j) = area_cloud_fraction(i, j, 1)
          End Do
        End Do  ! Diag09203_do1
!
! Cloud amount is calculated under maximum overlap assumption.
! Diag09202_do2:
        Do k = 1, MAX((LOW_BOT_LEVEL - 1), 1)
          Do j = 1, rows
            Do i = 1, row_length
              work2d_1(i, j) =                                          &
     &             MAX(work2d_1(i, j), area_cloud_fraction(i, j, k))
            End Do
          End Do
        End Do  ! Diag09202_do2
!
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode  >   0) then
          cmessage="cld_ctl  : error in copydiag(LOW Cloud Amount)"
         End if
      End if  ! Diag09203_if1
!
! ----------------------------------------------------------------------
! DIAG.09203 Copy Low Cloud Amount to stashwork
! ----------------------------------------------------------------------
      item = 203
! Diag09203_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Initialize cloud amount to lowest level value.
! Diag09203_do1:
        Do j = 1, rows
          Do i = 1, row_length
            work2d_1(i, j) = area_cloud_fraction(i, j, LOW_BOT_LEVEL)
          End Do
        End Do  ! Diag09203_do1
!
! Cloud amount is calculated under maximum overlap assumption.
! Diag09203_do2:
        Do k = LOW_BOT_LEVEL + 1, LOW_TOP_LEVEL
          Do j = 1, rows
            Do i = 1, row_length
              work2d_1(i, j) =                                          &
     &             MAX(work2d_1(i, j), area_cloud_fraction(i, j, k))
            End Do
          End Do
        End Do  ! Diag09203_do2
!
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode  >   0) then
          cmessage="cld_ctl  : error in copydiag(LOW Cloud Amount)"
         End if
      End if  ! Diag09203_if1
!
! ----------------------------------------------------------------------
! DIAG.09204 Copy Medium Cloud Amount to stashwork
! ----------------------------------------------------------------------
      item = 204
! Diag09204_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Initialize cloud amount to lowest level value.
! Diag09204_do1:
        Do j = 1, rows
          Do i = 1, row_length
            work2d_1(i, j) = area_cloud_fraction(i, j, MED_BOT_LEVEL)
          End Do
        End Do  ! Diag09204_do1
!
! Cloud amount is calculated under maximum overlap assumption.
! Diag09204_do2:
        Do k = MED_BOT_LEVEL + 1, MED_TOP_LEVEL
          Do j = 1, rows
            Do i = 1, row_length
              work2d_1(i, j) =                                          &
     &             MAX(work2d_1(i, j), area_cloud_fraction(i, j, k))
            End Do
          End Do
        End Do  ! Diag09204_do2
!
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode  >   0) then
          cmessage="cld_ctl  : error in copydiag(MEDIUM Cloud Amount)"
        End if
      End if  ! Diag09204_if1
!
! ----------------------------------------------------------------------
! DIAG.09205 Copy High Cloud Amount to stashwork
! ----------------------------------------------------------------------
      item = 205
! Diag09205_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Initialize cloud amount to lowest level value.
! Diag09205_do1:
        Do j = 1, rows
          Do i = 1, row_length
            work2d_1(i, j) = area_cloud_fraction(i, j, HIGH_BOT_LEVEL)
          End Do
        End Do  ! Diag09205_do1
!
! Cloud amount is calculated under maximum overlap assumption.
! Diag09205_do2:
        Do k = HIGH_BOT_LEVEL + 1, HIGH_TOP_LEVEL
          Do j = 1, rows
            Do i = 1, row_length
              work2d_1(i, j) =                                          &
     &             MAX(work2d_1(i, j), area_cloud_fraction(i, j, k))
            End Do
          End Do
        End Do  ! Diag09205_do2
!
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode  >   0) then
          cmessage="cld_ctl  : error in copydiag(HIGH Cloud Amount)"
        End if
      End if  ! Diag09205_if1
!
! ----------------------------------------------------------------------
! DIAG.09206 Copy cloud liquid condensate to stashwork
! ----------------------------------------------------------------------
! Needed for Assimilation, otherwise duplicates 00254 output (preferred)
      item = 206
! Diag09206_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
      If (.not.allocated(qcl_lsc)) Then
        Allocate ( qcl_lsc(row_length*rows,wet_model_levels) )
      End If

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(stashwork(si(item,sect,im_index)),qcl,        &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

! Copy qcl to qcl_lsc for ac_diagnostics_mod

      do k = 1,wet_model_levels
        do j = 1,rows
          do i = 1,row_length
            ji = (j-1)*row_length+i
            qcl_lsc(ji,k) = qcl(i,j,k)
          end do
        end do
      end do

         If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(Qcl Main Cloud)"
         End if
      End if  ! Diag09206_if1
!
! ----------------------------------------------------------------------
! DIAG.09207 Copy cloud frozen condensate to stashwork
! ----------------------------------------------------------------------
      If (.false.) Then
! Prevent access to diagnostic after repositioning of routine until a
! final decision is made on if/how to ingest it. Timetable UM5.1 -> 5.2
      item = 207
! Diag09207_if1:
!     If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(stashwork(si(item,sect,im_index)),qcf,        &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage="cld_ctl  : error in copydiag_3d(cloud ice)"
         End if
      End if  ! Diag09207_if1
!
! ----------------------------------------------------------------------
! Find cloud base for pre defined cloud cover threshholds
! for Diagnostics 09208 - 09215.
! ----------------------------------------------------------------------
      IF (icode  <=  0  .AND. ( sf(208,sect) .OR.                       &
     &                          sf(209,sect) .OR.                       &
     &                          sf(210,sect) .OR.                       &
     &                          sf(211,sect) .OR.                       &
     &                          sf(212,sect) .OR.                       &
     &                          sf(213,sect) .OR.                       &
     &                          sf(214,sect) .OR.                       &
     &                          sf(215,sect))) THEN
! Initialise output arrays
        cloud_base(:,:,:)=RMDI
        cloud_base_level(:,:,:)=IMDI

! Set cloud base to model levels
        DO l=1,noktas                 ! Loop over threshholds
          DO k=1,wet_model_levels     ! Loop over levels
            DO j=1,rows               ! Loop over rows
              DO i=1,row_length       ! Loop over points
                kinvert = wet_model_levels + 1 - k
  ! Convert to oktas
                combined_cloud_in_oktas=combined_cloud(i,j,kinvert)     &
     &                                  *roktas
  ! Calculate cloud_base_level
                IF(combined_cloud_in_oktas >  oktas(l))THEN
                  IF(cloud_base_level(i,j,l) <  0)THEN
                    cloud_base_level(i,j,l) = k
                  END IF
                END IF
              END DO
            END DO
          END DO
        END DO

! Convert level numbers to heights (M converted to Kft)
        DO l=1,noktas                 ! Loop over threshholds
          DO j=1,rows                 ! Loop over rows
            DO i=1,row_length         ! Loop over points
              IF(cloud_base_level(i,j,l) > 0)THEN
                cloud_base_height =                                     &
     &                 r_rho_levels(i,j,cloud_base_level(i,j,l))        &
     &                 - earth_radius
                cloud_base(i,j,l)=cloud_base_height*m_to_kft
              END IF
            END DO
          END DO
        END DO
      END IF
!
! ----------------------------------------------------------------------
! DIAG.09208 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 0.1 okta
! ----------------------------------------------------------------------
      item = 208
! Diag09208_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,1),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 0.1 okta)"
        END IF
      END IF  ! Diag09208_if1
!
! ----------------------------------------------------------------------
! DIAG.09209 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 1.5 okta
! ----------------------------------------------------------------------
      item = 209
! Diag09209_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,2),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 1.5 okta)"
        END IF
      END IF  ! Diag09209_if1
!
! ----------------------------------------------------------------------
! DIAG.09210 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 2.5 okta
! ----------------------------------------------------------------------
      item = 210
! Diag09210_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,3),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 2.5 okta)"
        END IF
      END IF  ! Diag09210_if1
!
! ----------------------------------------------------------------------
! DIAG.09211 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 3.5 okta
! ----------------------------------------------------------------------
      item = 211
! Diag09211_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,4),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 3.5 okta)"
        END IF
      END IF  ! Diag09211_if1
!
! ----------------------------------------------------------------------
! DIAG.09212 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 4.5 okta
! ----------------------------------------------------------------------
      item = 212
! Diag09212_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,5),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 4.5 okta)"
        END IF
      END IF  ! Diag09212_if1
!
! ----------------------------------------------------------------------
! DIAG.09213 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 5.5 okta
! ----------------------------------------------------------------------
      item = 213
! Diag09213_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,6),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 5.5 okta)"
        END IF
      END IF  ! Diag09213_if1
!
! ----------------------------------------------------------------------
! DIAG.09214 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 6.5 okta
! ----------------------------------------------------------------------
      item = 214
! Diag09214_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,7),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 6.5 okta)"
        END IF
      END IF  ! Diag09214_if1
!
! ----------------------------------------------------------------------
! DIAG.09215 Copy Height of Lowest Cloud Base for
!                 Cloud Cover &gt; 7.9 okta
! ----------------------------------------------------------------------
      item = 215
! Diag09215_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),                &
     &       cloud_base(1,1,8),                                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Ht lowest cld base > 7.9 okta)"
        END IF
      END IF  ! Diag09215_if1
!
! ***** Code matches definition in COMMON RADIATION Section 70. ****
      nclds = MIN(cloud_levels, wet_model_levels)
!
! ----------------------------------------------------------------------
! DIAG.09216 Copy Total Cloud Amount RANDOM Overlap to stashwork
! ----------------------------------------------------------------------
      item = 216
! NB: Convention in Sect 70 (Radiation) is to invert levels, 1 at top.
! Diag09216_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: r2_calc_total_cloud_cover
        Call R2_calc_total_cloud_cover(                                 &
     &         row_length*rows, wet_model_levels, nclds                 &
     &       , IP_CLOUD_MIX_RANDOM, combined_cloud(1,1,1), work2d_1     &
     &       , row_length*rows, wet_model_levels                        &
     &        )
!
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="cld_ctl  : error in copydiag(total cloud Random)"
         Endif
      Endif  ! Diag09216_if1
!
! ----------------------------------------------------------------------
! DIAG.09217 Copy Total Cloud Amount MAX/RANDOM Overlap to stashwork
! ----------------------------------------------------------------------
      item = 217
! NB: Convention in Sect 70 (Radiation) is to invert levels, 1 at top.
! Diag09217_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: r2_calc_total_cloud_cover
        Call R2_calc_total_cloud_cover(                                 &
     &         row_length*rows, wet_model_levels, nclds                 &
     &       , IP_CLOUD_MIX_MAX, combined_cloud(1,1,1), work2d_1        &
     &       , row_length*rows, wet_model_levels                        &
     &        )
!
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)),work2d_1,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

        If (icode  >   0) then
          cmessage="cld_ctl  : error in copydiag(total cloud Max/Rand)"
        Endif
      Endif  ! Diag09217_if1
! ----------------------------------------------------------------------
! Find Cloud Fraction in air &lt; STR_CEIL, Low Cloud Base and Top
! for Diagnostics 09218 and/or 09219 and/or 09220
! ----------------------------------------------------------------------
      IF (icode  <=  0  .AND. ( sf(218,sect) .OR.                       &
     &                          sf(219,sect) .OR.                       &
     &                          sf(220,sect))) THEN
        DO j=1,rows               ! Loop over rows
          DO i=1,row_length       ! Loop over points
            pt=0.0
            ft=0.0
            work2d_1(i,j)=RMDI    ! Initialise work arrays
            work2d_2(i,j)=RMDI
            work2d_3(i,j)=RMDI
            h_asln=r_theta_levels(i,j,0) - earth_radius
            pu=p_star(i,j)
            DO k=1,min(wet_model_levels,model_levels-1)
              h_asl=h_asln
              h_asln=r_rho_levels(i,j,k+1) - earth_radius

! Check if have not already found low cloud base
              IF(work2d_2(i,j) == RMDI) THEN
! IF not, and cloud cover in this layer &gt; threshhold
                IF(area_cloud_fraction(i,j,k) >= cloud_threshold) THEN
! THEN call the bottom of this layer the base of low cloud
                  work2d_2(i,j) = h_asl / ft2m
                END IF
              END IF

! Check if already found low cloud base but not top
              IF((work2d_2(i,j) /= RMDI).AND.(work2d_3(i,j) == RMDI))   &
     &          THEN
! IF not, and cloud cover in this layer &lt; threshhold
                IF(area_cloud_fraction(i,j,k) <= cloud_threshold) THEN
! THEN call the bottom of this layer the top of low cloud
                  work2d_3(i,j) = h_asl / ft2m
                END IF
              END IF

! IF bottom of layer is below low cloud ceiling (1000ft)
               IF(h_asl  <   str_ceilm) THEN
! Calculate top and bottom layer pressures
                 pl = pu
                 pu = p(i,j,k+1)
! And accumulate presure thickness and pressure weighted cloud amount
                 dp = pu - pl
! IF whole layer below ceiling, simply accumulate whole layer.
                 IF(h_asln  <   str_ceilm) THEN
                   pt = pt + dp
                   ft = ft + dp * area_cloud_fraction(i,j,k)
                 ELSE
! Otherwise height interpolate
                  fr = (str_ceilm - h_asl ) / (h_asln - h_asl )
                  pt = pt + dp * fr
                  ft = ft + dp * area_cloud_fraction(i,j,k) * fr
! And set results
                  work2d_1(i,j) = ft / pt

                 END IF
               END IF
             END DO               ! k over max
           END DO                 ! END loop over points
         END DO                   ! END loop over rows
       END IF
!
! ----------------------------------------------------------------------
! DIAG.09218 Copy Cloud fraction below 1000ft ASL to stashwork.
! ----------------------------------------------------------------------
      item = 218
! Diag09218_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)), work2d_1,      &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     & "cld_ctl  : error in copydiag(Cloud fraction below 1000ft ASL)"
        END IF
      END IF  ! Diag09218_if1
!
! ----------------------------------------------------------------------
! DIAG.09219 Copy Low Cloud Base to stashwork.
! ----------------------------------------------------------------------
      item = 219
! Diag09219_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)), work2d_2,      &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Low cloud base)"
        END IF
      END IF  ! Diag09219_if1
!
! ----------------------------------------------------------------------
! DIAG.09220 Copy Low Cloud Top to stashwork.
! ----------------------------------------------------------------------
      item = 220
! Diag09220_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)), work2d_3,      &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Low cloud top)"
        END IF
      END IF  ! Diag09220_if1
!
! ----------------------------------------------------------------------
! Calculate wet bulb temperature for Diagnostics 09221 and/or 09222
! ----------------------------------------------------------------------
      IF (icode  <=  0  .AND. ( sf(221,sect) .OR.                       &
     &                          sf(222,sect))) THEN

        DO k = 1, wet_model_levels

! Initialize 2d work arrays
          DO j = 1, rows
            DO i = 1, row_length
              work2d_1(i,j) = t(i,j,k)
              work2d_2(i,j) = q(i,j,k)
              work2d_3(i,j) = p_theta_levels(i,j,k)
            END DO
          END DO

! Calculate tw
! DEPENDS ON: thetaw
          Call ThetaW(row_length*rows,                                  &
     &                work2d_1, work2d_2, work2d_3, l_potential,        &
     &                work3d_1(1,1,k))
        END DO
      END IF
!
! ----------------------------------------------------------------------
! DIAG.09222 Copy Wet Bulb Temperature to stashwork.
! ----------------------------------------------------------------------
      item = 222
! Diag09222_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d(stashwork(si(item,sect,im_index)), work3d_1,   &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag_3d(Wet Bulb Temperature)"
        END IF
      END IF  ! Diag09222_if1
!
! ----------------------------------------------------------------------
! Calculate wet bulb freezing level for Diagnostic 09221
! ----------------------------------------------------------------------
      IF (icode  <=  0  .AND.  sf(221,sect)) THEN
        DO j = 1, rows                 ! Loop over rows
          DO i = 1, row_length         ! Loop over points
            DO k = 1, wet_model_levels ! Loop over all wet-levels
              IF (work3d_1(i,j,k)  /=  rmdi) THEN
                IF (work3d_1(i,j,k)  <=  zerodegc) THEN
                  IF ( k  ==  1) THEN
                    work2d_1(i,j)=r_theta_levels(i,j,k) -               &
     &                            earth_radius
                  ELSE
                    frac = (zerodegc - work3d_1(i,j,k-1))/              &
     &                     (work3d_1(i,j,k)-work3d_1(i,j,k-1))
                    work2d_1(i,j)=r_theta_levels(i,j,k)*frac +          &
     &                            r_theta_levels(i,j,k-1)*(1.0-frac)    &
     &                            - earth_radius
                  END IF
                  EXIT
                END IF
              ELSE
                work2d_1(i,j) = rmdi
              END IF
            END DO
          END DO
        END DO
      END IF
!
! ----------------------------------------------------------------------
! DIAG.09221 Copy Wet Bulb Freezing Level to stashwork.
! ----------------------------------------------------------------------
      item = 221
! Diag09221_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)), work2d_1,      &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Wet Bulb Freezing Level)"
        END IF
      END IF  ! Diag09221_if1

!
! ----------------------------------------------------------------------
! Calculate Total Cloud Top Height for Diagnostic 09223
! ----------------------------------------------------------------------
      IF (icode  <=  0  .AND. sf(223,sect)) THEN
        DO j = 1, rows                   ! Loop over rows
          DO i = 1, row_length           ! Loop over points
            cld_top_lev = IMDI
            DO k = 1, wet_model_levels   ! Loop over all wet-levels
              IF (combined_cloud(i,j,k)  >   thresh) THEN
                cld_top_lev = wet_model_levels + 1 - k
              END IF
              IF (cld_top_lev  >   0) EXIT
            END DO
            IF (cld_top_lev  >   0) THEN
              work2d_1(i,j) = r_theta_levels(i,j,cld_top_lev)           &
     &                        - earth_radius
! Convert to Kiloft
              work2d_1(i,j) = work2d_1(i,j) * m_to_kft
            ELSE
              work2d_1(i,j) = RMDI
            END IF
          END DO
        END DO
      END IF
!
! ----------------------------------------------------------------------
! DIAG.09223 Copy Total Cloud Top Height to stashwork.
! ----------------------------------------------------------------------
      item = 223
! Diag09223_if1:
      IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag
        Call copydiag(stashwork(si(item,sect,im_index)), work2d_1,      &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        IF (icode >  0) THEN
          cmessage=                                                     &
     &     "cld_ctl  : error in copydiag(Total cloud top height)"
        END IF
      END IF  ! Diag09223_if1
!
! ----------------------------------------------------------------------
! DIAG.09229 Copy Relative Humidity to stashwork
! ----------------------------------------------------------------------
      item = 229
! Diag09229_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Calculate Qsat then work out Relative Humidity.
! Diag09229_do1:
        Do k = 1, wet_model_levels
! DEPENDS ON: qsat_mix
          Call qsat_mix(work2d_1,T(1,1,K),p_theta_levels(1,1,k),        &
     &              row_length*rows,l_mixing_ratio)
!
          Do j = 1, rows
            Do i = 1, row_length
              ! relative humidity in per cent
              work3d_1(i, j, k) = q(i, j, k) / work2d_1(i, j) *100.0

!  Supersaturation (>100%) can occur with mixed phase scheme but
!  negative humidity is removed from the diagnostic:
              If ( work3d_1(i, j, k) < 0.0) Then
                  work3d_1(i, j, k) = 0.0
              Endif

            End Do
          End Do
!
        End Do  ! Diag09229_do1
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(RH Main Cloud)"
        End if
      End if  ! Diag09229_if1
!
! ----------------------------------------------------------------------
! DIAG.09226 Copy Layer Cloud Frequency to stashwork.
! ----------------------------------------------------------------------
      item = 226
! Diag09226_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Freq_k_do1:
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              If (area_cloud_fraction(i, j, k)  <=  0.)  Then
                work3d_1(i, j, k) = 0.
              Else
                work3d_1(i, j, k) = 1.
              End if
            End Do
          End Do
        End Do  ! Freq_k_do1
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),work3d_1,   &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(Layer Cloud Freq)"
        End if
      End if  ! Diag09226_if1
!
! ----------------------------------------------------------------------
! DIAG.09228 Copy Critical Relative Humidity to stashwork.
! ----------------------------------------------------------------------
      item = 228
! Diag09228_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!     STASH will only permit this diagnostic to be chosen if 3D RHcrit
!     diagnosis is active. So RHCPT should be (row_length,rows,wet_ml).
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),RHCPT,      &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)
!
        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(Critical RH)"
        End if
      End if  ! Diag09228_if1
!
! ----------------------------------------------------------------------
! DIAG.09230 Calculate visibility on model levels and copy to stashwork
! ----------------------------------------------------------------------
      item = 230
! Note that the following calculates the visibility for all model levels
! and then extracts just those levels that have been requested. If the
! diagnostic is normally called for a small subset of levels it will be
! more efficient to introduce a call to set_levels_list and calculate
! only those levels that are needed.
! Diag09230_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then

        do k=1,wet_model_levels

          ! visibility is calculated using prescribed critical RH
          ! (not diagnostic RH)
! DEPENDS ON: visbty
          Call Visbty(                                                  &
     &     p_theta_levels(1,1,k), T(1,1,k), q(1,1,k)                    &
                                                           !INPUT
     &     ,qcl(1,1,k), qcf(1,1,k)                                      &
                                                           !INPUT
     &     ,Aerosol(1:row_length,1:rows,k)                              &
                                                           !INPUT
     &     ,vis_probability, RHcrit(k), L_murk                          &
                                                           !INPUT
     &     ,row_length * rows                                           &
                                                           !INPUT
     &     ,work3d_1(1,1,k))                               !OUTPUT

        End do ! k   wet_model_levels

! DEPENDS ON: copydiag_3d
        Call copydiag_3d(stashwork(si(item,sect,im_index)),work3d_1,    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)


      End if  ! Diag09230_if1

! ----------------------------------------------------------------------
! DIAG.09231 Copy Combined Cloud fraction to stashwork
! ----------------------------------------------------------------------
      item = 231
! Diag09231_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! Diag09231_do1:
        Do k = 1, wet_model_levels
! NB: Convention in Sect 70 (Radiation) is to invert levels, 1 at top.
!     Combined_cloud is calculated this way but re-inverted for STASH.
!
          kinvert = wet_model_levels+1-k
          Do j = 1, rows
            Do i = 1, row_length
              work3d_1(i, j, k) = combined_cloud(i,j,kinvert)
            End Do
          End Do
        End Do  ! Diag09231_do1
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)), work3d_1,  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
          cmessage="cld_ctl  : error in copydiag_3d(Combined Cld On Lv)"
        End if
      End if  ! Diag09231_if1
!
! ----------------------------------------------------------------------
      If (icode  /=  0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,Cmessage)
      End if
!
      Return
      END SUBROUTINE diagnostics_lscld
#endif
