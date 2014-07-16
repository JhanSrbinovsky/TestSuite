#if defined(A09_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine qt_bal_cld
       subroutine qt_bal_cld                                            &
     &     (p_star,p_theta_levels,p,                                    &
     &      theta,exner_theta_levels,                                   &
     &      q,qcl,qcf,qcf2,                                             &
     &      rhcpt, rhc_row_length, rhc_rows, bl_levels,                 &
     &      cloud_fraction_method,overlap_ice_liquid,                   &
     &      ice_fraction_method,ctt_weight,t_weight,                    &
     &      qsat_fixed,sub_cld,                                         &
     &      row_length,rows,model_levels,wet_levels,                    &
     &      offx,offy,halo_i,halo_j,                                    &
     &      delta_lambda, delta_phi,                                    &
     &      r_theta_levels, FV_cos_theta_latitude,                      &
     &      lc,cp, L_cld_area, L_ACF_Cusack, L_ACF_Brooks, L_eacf,      &
     &      L_mcr_qcf2, l_mixing_ratio, ntml, cumulus,                  &
     &      area_cloud_fraction,  bulk_cloud_fraction,                  &
     &      cloud_fraction_liquid,  cloud_fraction_frozen,              &
     &      mype)

! Purpose:
!        reset q, t and the cloud fields to be consistent at the
!        end of the timestep
!
! Method:
!          Is described in ;
!
! Original Progammer: Andy Malcolm
! Current code owner: Andy Malcolm
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      implicit none

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, model_levels                                                    &
                         ! number of model levels
     &, wet_levels                                                      &
                   ! number of model levels where moisture
                         ! variables are held
     &, halo_i                                                          &
                             ! Size of halo in i direction.
     &, halo_j                                                          &
                             ! Size of halo in j direction.
     &, offx                                                            &
                            ! Size of small halo in i
     &, offy                                                            &
                            ! Size of small halo in j.
     &, mype                                                            &
                           ! My processor number
     &, rhc_row_length                                                  &
     &, rhc_rows                                                        &
     & , bl_levels
! physical constants
      Real                                                              &
     &  lc, cp

      Logical                                                           &
     &  L_CLD_AREA                                                      &
                           ! true if using area cloud fraction (ACF)
     &, L_ACF_Cusack                                                    &
                           ! ... and selected Cusack and PC2 off
     &, L_ACF_Brooks                                                    &
                           ! ... and selected Brooks
     &, L_eacf                                                          &
                           ! true if using empirically adjusted
                           !               cloud fraction
     &, L_mcr_qcf2                                                      &
                           ! true if second cloud ice variable in use
     &, l_mixing_ratio     ! true if using mixing ratio formulation

      Real                                                              &
     &  p(1-offx:row_length+offx, 1-offy:rows+offy,                     &
     &    model_levels+1)                                               &
     &, p_theta_levels(1-offx:row_length+offx, 1-offy:rows+offy,        &
     &    model_levels)                                                 &
     &, p_star(row_length, rows)                                        &
     &, theta(1-offx:row_length+offx, 1-offy:rows+offy,                 &
     &        model_levels)                                             &
     &, exner_theta_levels(1-offx:row_length+offx, 1-offy:rows+offy,    &
     &    model_levels)                                                 &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &    wet_levels)                                                   &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &      wet_levels)                                                 &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &      wet_levels)                                                 &
     &, qcf2(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &      wet_levels)                                                 &
     &, rhcpt(rhc_row_length, rhc_rows, wet_levels)                     &
! coordinate arrays
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
     &, delta_lambda                                                    &
     &, delta_phi                                                       &
! trig arrays
     &, FV_cos_theta_latitude (1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy)
!

      Integer                                                           &
                        !, Intent(IN)
     &  cloud_fraction_method                                           &
                               ! Method for calculating
                               ! total cloud fraction
     &, ice_fraction_method    ! Method for calculating ice cloud frac.

      Real                                                              &
                        !, Intent(IN)
     &  overlap_ice_liquid                                              &
                               ! Overlap between ice and liquid phases
     &, ctt_weight                                                      &
                               ! Weighting of cloud top temperature
     &, t_weight                                                        &
                               ! Weighting of local temperature
     &, qsat_fixed                                                      &
                               ! Fixed value of saturation humidity
     &, sub_cld                ! Scaling parameter

! Diagnostic variables

      Real                                                              &
     &  area_cloud_fraction(row_length, rows, wet_levels)               &
     &, bulk_cloud_fraction(1-halo_i:row_length+halo_i,                 &
     &                      1-halo_j:rows+halo_j, wet_levels)           &
     &, cloud_fraction_liquid(1-halo_i:row_length+halo_i,               &
     &                        1-halo_j:rows+halo_j, wet_levels)         &
     &, cloud_fraction_frozen(1-halo_i:row_length+halo_i,               &
     &                        1-halo_j:rows+halo_j, wet_levels)


      Integer                                                           &
     &  ntml (row_length, rows)

      Logical                                                           &
     &  cumulus (row_length, rows) ! bl convection flag

! Local variables
      Real                                                              &
     & p_layer_centres(row_length, rows, 0:model_levels) ,              &
              ! pressure at layer centres. Same as p_theta_levels
              !except bottom level = p_star, and at top = 0.
     & p_layer_boundaries(row_length, rows, 0:model_levels) ,           &
     & tl(row_length, rows, wet_levels) ,                               &
     & qt(row_length, rows, wet_levels) ,                               &
     & qcf_in(row_length, rows, wet_levels) ,                           &
     & qcl_out(row_length, rows, wet_levels),                           &
     & cf_inout(row_length, rows, wet_levels),                          &
     & cfl_inout(row_length, rows, wet_levels),                         &
     & cff_inout(row_length, rows, wet_levels)


      Integer                                                           &
     &  large_levels                                                    &
     &, levels_per_level

      Integer                                                           &
     &  i,j,k,errorstatus     ! loop variables

      Do j = 1, rows
        Do i = 1, row_length
          p_layer_centres(i,j,0) = p_star(i,j)
          p_layer_boundaries(i,j,0) = p_star(i,j)
        End Do
      End Do
      Do k = 1, model_levels-1
        Do j = 1, rows
          Do i = 1, row_length
            p_layer_centres(i,j,k) = p_theta_levels(i,j,k)
            p_layer_boundaries(i,j,k) = p(i,j,k+1)
          End Do
        End Do
      End Do
      k=model_levels
      Do j = 1, rows
        Do i = 1, row_length
          p_layer_centres(i,j,k) = p_theta_levels(i,j,k)
          p_layer_boundaries(i,j,k) = 0.0
        End Do
      End Do

! ----------------------------------------------------------------------
! Section  Convert qT and Tl for input to cloud scheme.
! ----------------------------------------------------------------------

! Create Tl and qT
      Do k = 1, wet_levels
        Do j = 1, rows
          Do i = 1, row_length
            Tl(i,j,k) =                                                 &
     &             (theta(i,j,k) *exner_theta_levels(i,j,k))            &
     &                     - (lc * qcl(i,j,k)) / cp
            qt(i,j,k) = q(i,j,k) + qcl(i,j,k)
            qcf_in(i,j,k)=qcf(i,j,k)
            cf_inout(i,j,k)= bulk_cloud_fraction(i,j,k)
            cfl_inout(i,j,k)= cloud_fraction_liquid(i,j,k)
            cff_inout(i,j,k)= cloud_fraction_frozen(i,j,k)
          End Do
        End Do
      End Do

      ! If second cloud ice variable in use then add to qcf
      ! for the cloud scheme call
      If (L_mcr_qcf2)                                                   &
     &  qcf_in(:,:,:) = qcf_in(:,:,:) + qcf2(1:row_length,1:rows,:)

! ----------------------------------------------------------------------
! Section BL.4b Call cloud scheme to convert Tl and qT to T, q and qcl
!              calculate bulk_cloud fields from qT and qcf
!               and calculate area_cloud fields.
! ----------------------------------------------------------------------
!

! Determine number of sublevels for vertical gradient area cloud
! Want an odd number of sublevels per level: 3 is hardwired in do loops
      levels_per_level = 3
      large_levels = ((wet_levels - 2)*levels_per_level) + 2
!
! DEPENDS ON: ls_arcld
      CALL ls_arcld( p_layer_centres, RHCPT, p_layer_boundaries,        &
     &                 model_levels, wet_levels, row_length, rows,      &
     &                 rhc_row_length, rhc_rows, bl_levels,             &
     &                 cloud_fraction_method,overlap_ice_liquid,        &
     &                 ice_fraction_method,ctt_weight,t_weight,         &
     &                 qsat_fixed,sub_cld,                              &
     &                 levels_per_level, large_levels,                  &
     &                 L_cld_area, L_ACF_Cusack,L_ACF_Brooks, L_eacf,   &
     &                 halo_i, halo_j, offx, offy,                      &
     &                 delta_lambda, delta_phi,                         &
     &                 r_theta_levels, FV_cos_theta_latitude,           &
     &                 ntml, cumulus, l_mixing_ratio, qcf_in,           &
     &                 Tl, qt, qcl_out,                                 &
     &                 area_cloud_fraction,  cf_inout,                  &
     &                 cfl_inout, cff_inout ,                           &
     &                 errorstatus, mype)

! qt holds q (no halos), tl holds t(no halos),
! qcl_out holds qcl(no halos)
      Do k = 1, wet_levels
        Do j = 1, rows
          Do i = 1, row_length
              theta(i,j,k) = tl(i,j,k)/exner_theta_levels(i,j,k)
              q(i,j,k) = qt(i,j,k)
              qcl(i,j,k) =qcl_out(i,j,k)
              bulk_cloud_fraction(i,j,k) = cf_inout(i,j,k)
              cloud_fraction_liquid(i,j,k) = cfl_inout(i,j,k)
              cloud_fraction_frozen(i,j,k) = cff_inout(i,j,k)
          End Do
        End Do
      End Do

      return
      END SUBROUTINE qt_bal_cld
#endif
