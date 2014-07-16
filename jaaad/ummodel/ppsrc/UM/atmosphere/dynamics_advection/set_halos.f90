
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine set_halos
      subroutine set_halos(                                             &
     &                     u, v, w, u_adv, v_adv, w_adv, theta,         &
     &                     q, qcl, qcf, qcf2, qrain, qgraup,            &
     &                     cf, cfl, cff,                                &
     &                     rho, p, p_theta_levels,                      &
     &                     exner_rho_levels, exner_theta_levels,        &
     &                     MURK,                                        &
     &                     DUST_DIV1,DUST_DIV2,DUST_DIV3,               &
     &                     DUST_DIV4,DUST_DIV5,DUST_DIV6,               &
     &                     SO2, SO4_AITKEN, SO4_ACCU, SO4_DISS,         &
     &                     dms, nh3, soot_new, soot_agd, soot_cld,      &
     &                     bmass_new, bmass_agd, bmass_cld,             &
     &                     ocff_new, ocff_agd, ocff_cld,                &
     &                     co2, free_tracers, ukca_tracers,             &
     &                     row_length, rows, n_rows, model_levels,      &
     &                     wet_levels, offx, offy, halo_i, halo_j,      &
     &                     TR_LEVELS, TR_VARS, tr_ukca, L_MURK,         &
     &                     L_DUST, L_SULPC_SO2,                         &
     &                     l_sulpc_nh3, l_sulpc_dms, l_soot,            &
     &                     l_biomass, l_ocff, l_co2_interactive,        &
     &                     L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,       &
     &                     OZONE_TRACER,                                &
     &                     L_USE_CARIOLLE                               &
     &                     )

! Purpose:
!          set halos of arrays in a subroutine to tidy ATM_STEP.
!
! Method:
!          Is described in ;
!
!
! Original Programmer: Andrew J. Malcolm
! Current code owner: Andrew J. Malcolm
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Use swapable_field_mod, Only: &
          swapable_field_pointer_type

      Implicit None

! Arguments with Intent IN. ie: Input variables.
      Integer, Intent(In) ::                                            &
     &       ROW_LENGTH,                                                &
                                  ! IN: No of points per local row
     &       ROWS,                                                      &
                                  ! IN: No of local (theta) rows
     &       n_ROWS,                                                    &
                                  ! IN: No of local (v) rows
     &       MODEL_LEVELS,                                              &
                                  ! IN: No of model levels
     &       WET_LEVELS,                                                &
                                  ! IN: No of moist-levels
     &       Offx,                                                      &
                                  ! standard halo size in East-West
     &       Offy,                                                      &
                                  ! standard halo size in North-South
     &       halo_i,                                                    &
                                  ! extended halo size in East-West
     &       halo_j,                                                    &
                                  ! extended halo size in North-South
     &       tr_vars,                                                   &
                                  ! number of free tracer variables
     &       tr_ukca,                                                   &
                                  ! number of ukca tracer variables
     &       tr_levels            ! number of free tracer levels

      Logical, Intent(In) ::                                            &
                                  ! Inclusion switches.
     &       l_murk,                                                    &
                                  ! aerosol
     &       L_DUST,                                                    &
                                  ! mineral dust
     &       l_sulpc_so2,                                               &
                                  ! sulphur cycle - main switch
     &       l_sulpc_dms,                                               &
                                  ! dms in sulphur cycle
     &       l_sulpc_nh3,                                               &
                                  ! nh3 in sulphur cycle
     &       l_soot,                                                    &
                                  ! Soot cycle
     &       l_biomass,                                                 &
                                  ! Biomass aerosol
     &       l_ocff,                                                    &
                                  ! Fossil-fuel organic carbon aerosol
     &       l_co2_interactive                                          &
                                  ! Active carbon cycle
     &,      L_mcr_qcf2                                                 &
                                  ! T => Use second ice prognostic
     &,      L_mcr_qrain                                                &
                                  ! T => Use rain prognostic
     &,      L_mcr_qgraup,                                              &
                                  ! T => Use graupel prognostic
     &       L_USE_CARIOLLE
                                  ! Use cariolle tracer scheme

      Real, Target, Intent (InOut) ::                                   &
     &  u(1-offx:row_length+offx, 1-offy:rows+offy,                     &
     &      model_levels)                                               &
     &, v(1-offx:row_length+offx, 1-offy:n_rows+offy,                   &
     &      model_levels)                                               &
     &, w(1-offx:row_length+offx, 1-offy:rows+offy,                     &
     &      0:model_levels)                                             &
     &, u_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &          model_levels)                                           &
     &, v_adv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,       &
     &          model_levels)                                           &
     &, w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &          0:model_levels)                                         &
     &, theta(1-offx:row_length+offx, 1-offy:rows+offy,                 &
     &          model_levels)                                           &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &      wet_levels)                                                 &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_levels)                                               &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_levels)                                               &
     &, qcf2(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &        wet_levels)                                               &
     &, qrain(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        wet_levels)                                               &
     &, qgraup(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &        wet_levels)                                               &
     &, cf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
     &      wet_levels)                                                 &
     &, cfl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_levels)                                               &
     &, cff(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_levels)                                               &
     &, rho(1-offx:row_length+offx, 1-offy:rows+offy,                   &
     &        model_levels)                                             &
     &, p(1-offx:row_length+offx, 1-offy:rows+offy,                     &
     &      model_levels+1)                                             &
     &, p_theta_levels(1-offx:row_length+offx,                          &
     &                   1-offy:rows+offy, model_levels)                &
     &, exner_rho_levels(1-offx:row_length+offx,                        &
     &                   1-offy:rows+offy, model_levels+1)              &
     &, exner_theta_levels(1-offx:row_length+offx,                      &
     &                     1-offy:rows+offy, model_levels)              &
     &, murk(1-offx:row_length+offx, 1-offy:rows+offy,                  &
     &       model_levels)                                              &
     &, DUST_DIV1(1-OFFX:ROW_LENGTH+OFFX, 1-OFFY:ROWS+OFFY,             &
     &       MODEL_LEVELS)                                              &
     &, DUST_DIV2(1-OFFX:ROW_LENGTH+OFFX, 1-OFFY:ROWS+OFFY,             &
     &       MODEL_LEVELS)                                              &
     &, DUST_DIV3(1-OFFX:ROW_LENGTH+OFFX, 1-OFFY:ROWS+OFFY,             &
     &       MODEL_LEVELS)                                              &
     &, DUST_DIV4(1-OFFX:ROW_LENGTH+OFFX, 1-OFFY:ROWS+OFFY,             &
     &       MODEL_LEVELS)                                              &
     &, DUST_DIV5(1-OFFX:ROW_LENGTH+OFFX, 1-OFFY:ROWS+OFFY,             &
     &       MODEL_LEVELS)                                              &
     &, DUST_DIV6(1-OFFX:ROW_LENGTH+OFFX, 1-OFFY:ROWS+OFFY,             &
     &       MODEL_LEVELS)                                              &
     &, so2(1-offx:row_length+offx, 1-offy:rows+offy,                   &
     &      model_levels)                                               &
     &, so4_aitken(1-offx:row_length+offx, 1-offy:rows+offy,            &
     &             model_levels)                                        &
     &, so4_accu(1-offx:row_length+offx, 1-offy:rows+offy,              &
     &           model_levels)                                          &
     &, so4_diss(1-offx:row_length+offx, 1-offy:rows+offy,              &
     &           model_levels)                                          &
     &, dms(1-offx:row_length+offx, 1-offy:rows+offy,                   &
     &      model_levels)                                               &
     &, nh3(1-offx:row_length+offx, 1-offy:rows+offy,                   &
     &      model_levels)                                               &
     &, soot_new(1-offx:row_length+offx, 1-offy:rows+offy,              &
     &           model_levels)                                          &
     &, soot_agd(1-offx:row_length+offx, 1-offy:rows+offy,              &
     &           model_levels)                                          &
     &, soot_cld(1-offx:row_length+offx, 1-offy:rows+offy,              &
     &           model_levels)                                          &
     &, bmass_new(1-offx:row_length+offx, 1-offy:rows+offy,             &
     &            model_levels)                                         &
     &, bmass_agd(1-offx:row_length+offx, 1-offy:rows+offy,             &
     &            model_levels)                                         &
     &, bmass_cld(1-offx:row_length+offx, 1-offy:rows+offy,             &
     &            model_levels)                                         &
     &, ocff_new(1-offx:row_length+offx, 1-offy:rows+offy,              &
     &           model_levels)                                          &
     &, ocff_agd(1-offx:row_length+offx, 1-offy:rows+offy,              &
     &           model_levels)                                          &
     &, ocff_cld(1-offx:row_length+offx, 1-offy:rows+offy,              &
     &           model_levels)                                          &
     &, co2(1-offx:row_length+offx, 1-offy:rows+offy,                   &
     &            model_levels)                                         &
     &, free_tracers(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &               tr_levels, tr_vars)                                &
     &, ukca_tracers(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &               tr_levels, tr_ukca)                                &
     &, OZONE_TRACER(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &            model_levels)  

! FLDTYPE definitions for the different field types recognised on the
! decomposition
      INTEGER,PARAMETER:: Nfld_max=7 ! maximum number of field types
      INTEGER,PARAMETER:: fld_type_p=1       ! grid on P points
      INTEGER,PARAMETER:: fld_type_u=2       ! grid on U points
      INTEGER,PARAMETER:: fld_type_v=3       ! grid on V points
      INTEGER,PARAMETER:: fld_type_comp_wave  = 4
                              ! Compressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_full_wave  = 5
                              ! Uncompressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_rim_wave   = 6
                              ! Boundary data for WAM Wave Field
      INTEGER,PARAMETER:: fld_type_r=7       ! grid on river points
      INTEGER,PARAMETER:: fld_type_unknown=-1! non-standard grid
! FLDTYPE end

! Local Variables
      Integer :: i    ! looper
      Integer :: i_field

! Note need to set the maximum number of fields in a single swap.
! 21 is tracer species (sulphates, dust, soot, co2, murk, cariolle 
! ozone) with all options true. tr_vars is number of free tracers.
! tr_ukca is the number of ukca tracers.

      Type(swapable_field_pointer_type) ::                              &
     &            fields_to_swap(max(24,tr_vars,tr_ukca))
!----------------------------------------------------------------------

        i_field = 0

! Multivariate swapbounds for fields with rows/row_length and
! extended halos

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => u_adv(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_u
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .TRUE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => v_adv(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_v
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  n_rows
        fields_to_swap(i_field) % vector      =  .TRUE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => w_adv(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  model_levels+1
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => q(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  wet_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => qcl(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  wet_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => qcf(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  wet_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.
        
        If (L_mcr_qcf2) Then
          i_field = i_field + 1
          fields_to_swap(i_field) % field       => qcf2(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  wet_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.
        End If

        If (L_mcr_qrain) Then        
          i_field = i_field + 1
          fields_to_swap(i_field) % field       => qrain(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  wet_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.
        End If

        If (L_mcr_qgraup) Then
          i_field = i_field + 1
          fields_to_swap(i_field) % field       => qgraup(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  wet_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.
        End If

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => cf(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  wet_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.
        
        i_field = i_field + 1
        fields_to_swap(i_field) % field       => cfl(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  wet_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.
        
        i_field = i_field + 1
        fields_to_swap(i_field) % field       => cff(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  wet_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.
        
! DEPENDS ON: swap_bounds_mv
        call Swap_Bounds_mv( fields_to_swap, i_field,                   &
                             row_length, halo_i, halo_j ) 

! Multivariate swapbounds for rows/row_length with single point halos
        i_field = 0
        i_field = i_field + 1
        fields_to_swap(i_field) % field       => u(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_u
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .TRUE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => v(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_v
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  n_rows
        fields_to_swap(i_field) % vector      =  .TRUE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => w(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  model_levels+1
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => theta(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.
        
        i_field = i_field + 1
        fields_to_swap(i_field) % field       => rho(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.
        
        i_field = i_field + 1
        fields_to_swap(i_field) % field       => p(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  model_levels + 1
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.
        
        i_field = i_field + 1
        fields_to_swap(i_field) % field       => p_theta_levels(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.
        
        i_field = i_field + 1
        fields_to_swap(i_field) % field       => exner_rho_levels(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  model_levels + 1
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => exner_theta_levels(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.

        
! DEPENDS ON: swap_bounds_mv
        call Swap_Bounds_mv( fields_to_swap, i_field,                   &
                             row_length, offx, offy) 

! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(p,row_length, rows,                    &
     &                     model_levels+1,offx, offy)

! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(p_theta_levels,row_length, rows,       &
     &                     model_levels,offx, offy)

! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(exner_theta_levels,row_length, rows,   &
     &                     model_levels,offx, offy)

! All the tracer fields should be updated only if they are used.
! Will check on the logicals for this.

        i_field = 0
        If (L_Murk) Then
          i_field = i_field + 1
          fields_to_swap(i_field) % field       => murk(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.
        End If

        IF (L_DUST) THEN        ! mineral dust

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => dust_div1(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => dust_div2(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => dust_div3(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => dust_div4(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => dust_div5(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => dust_div6(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

        ENDIF  ! mineral dust


        If (L_Sulpc_SO2) Then    ! sulphur cycle
          i_field = i_field + 1
          fields_to_swap(i_field) % field       => so2(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => so4_aitken(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => so4_accu(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => so4_diss(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          If (l_sulpc_dms) Then
            i_field = i_field + 1
            fields_to_swap(i_field) % field       => dms(:,:,:)
            fields_to_swap(i_field) % field_type  =  fld_type_p
            fields_to_swap(i_field) % levels      =  model_levels
            fields_to_swap(i_field) % rows        =  rows
            fields_to_swap(i_field) % vector      =  .FALSE.
          End If

          If (l_sulpc_nh3) Then
            i_field = i_field + 1
            fields_to_swap(i_field) % field       => nh3(:,:,:)
            fields_to_swap(i_field) % field_type  =  fld_type_p
            fields_to_swap(i_field) % levels      =  model_levels
            fields_to_swap(i_field) % rows        =  rows
            fields_to_swap(i_field) % vector      =  .FALSE.
          End If

        End If     ! sulphur cycle



        If (L_Soot) Then    ! soot cycle
          i_field = i_field + 1
          fields_to_swap(i_field) % field       => soot_new(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => soot_agd(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => soot_cld(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.
        End If    ! soot cycle

        If (L_Biomass) Then    ! biomass aerosol
          i_field = i_field + 1
          fields_to_swap(i_field) % field       => bmass_new(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => bmass_agd(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => bmass_cld(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.
        End If    ! biomass aerosol


        If (L_CO2_Interactive) Then  ! Active CO2 cycle
          i_field = i_field + 1
          fields_to_swap(i_field) % field       => co2(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.
        End If    ! active co2 cycle

        If (L_OCFF) Then       ! Fossil-fuel organic carbon aerosol
          i_field = i_field + 1
          fields_to_swap(i_field) % field       => ocff_new(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => ocff_agd(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.

          i_field = i_field + 1
          fields_to_swap(i_field) % field       => ocff_cld(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.
        End If    ! OCFF cycle     
        
        If (L_USE_CARIOLLE) Then  ! Active Ozone Cariolle Scheme
          i_field = i_field + 1
          fields_to_swap(i_field) % field       => ozone_tracer(:,:,:)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.
        End If    ! active ozone cariolle scheme

! All tracer fields up to now
        If (i_field > 0) Then
! DEPENDS ON: swap_bounds_mv
          Call Swap_Bounds_mv( fields_to_swap, i_field,                 &
                               row_length, offx, offy) 
        End If

        i_field = 0
! Need to loop over free tracers.
        Do i = 1, tr_vars
          i_field = i_field + 1
          fields_to_swap(i_field) % field       => free_tracers(:,:,:,i)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.
        End Do

        If (i_field > 0) Then
! DEPENDS ON: swap_bounds_mv
          Call Swap_Bounds_mv( fields_to_swap, i_field,                 &
                               row_length, offx, offy) 
        End If

        i_field = 0
! And ukca tracers.
        Do i = 1, tr_ukca
          i_field = i_field + 1
          fields_to_swap(i_field) % field       => ukca_tracers(:,:,:,i)
          fields_to_swap(i_field) % field_type  =  fld_type_p
          fields_to_swap(i_field) % levels      =  model_levels
          fields_to_swap(i_field) % rows        =  rows
          fields_to_swap(i_field) % vector      =  .FALSE.
        End Do

        If (i_field > 0) Then
! DEPENDS ON: swap_bounds_mv
          Call Swap_Bounds_mv( fields_to_swap, i_field,                 &
                               row_length, offx, offy) 
        End If

      Return
      END SUBROUTINE set_halos

