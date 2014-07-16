


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
Module swrdiag_mod
!+ ---------------------------------------------------------------------
!  Module to define diagnostic structures for SW diagnostics.
!
!  Current Code Owner: J. M. Edwards
!
!- ----------------------------------------------------------------------
!
!     Note for coders: Adding new diagnostics to the radiation
!                      code from 5.3
!
!     As the number of diagnostics required from the model has grown
!     the old practice of returning diagnostics from the radiation
!     code has resulted in our reaching the limit on the number of
!     continuation lines allowed. To obviate this problem diagnostics
!     are now bundled into a structure. For convenience, those, such
!     as the SW flux below 690 nm, that may be used to advance the model
!     are not included; those not calculated within the radiation code
!     itself are also omitted. Neither of these points is an absolute
!     requirement.
!
!     New diagnostics should be added within this structure by creating
!     a flag and a pointer array for the diagnostic. Under control of
!     STASH, allocate space for the diagnostic and initialize to 0
!     in the SW (to fill unlit points) in RAD_CTL and remember to
!     deallocate the space later. Below RAD_CTL, the element of the
!     structure can be treated as a normal array. 2-D indexing is
!     used to match the structure of the main fields in the model.
!
!
      TYPE StrSWDiag
!
        Logical ::  L_solar_out_toa
        Logical ::  L_solar_out_clear
        Logical ::  L_surface_down_flux
        Logical ::  L_surf_down_clr
        Logical ::  L_surf_up_clr
        Logical ::  L_clear_hr
        Logical ::  L_net_flux_trop
        Logical ::  L_up_flux_trop
        Logical ::  L_flux_direct
        Logical ::  L_flux_diffuse
        Logical ::  re_conv_flag
        Logical ::  re_strat_flag
        Logical ::  wgt_conv_flag
        Logical ::  wgt_strat_flag
        Logical ::  lwp_strat_flag
        Logical ::  weighted_re_flag
        Logical ::  sum_weight_re_flag
        Logical ::  wgtd_warm_re_flag
        Logical ::  sum_wgt_warm_re_flag
        Logical ::  ntot_diag_flag
        Logical ::  strat_lwc_diag_flag
        Logical ::  so4_ccn_diag_flag
        Logical ::  cond_samp_wgt_flag
        Logical ::  seasalt_film_flag
        Logical ::  seasalt_jet_flag
        Logical ::  Nc_diag_flag
        Logical ::  Nc_weight_flag
        Logical ::  L_FlxSolBelow690nmSurf
        Logical ::  L_FlxSeaBelow690nmSurf
        Logical ::  L_cloud_extinction
        Logical ::  L_cloud_weight_extinction
        Logical ::  L_ls_cloud_extinction
        Logical ::  L_ls_cloud_weight_extinction
        Logical ::  L_cnv_cloud_extinction
        Logical ::  L_cnv_cloud_weight_extinction
        Logical ::  L_direct_par
! Logicals for Orography
        Logical ::  L_orog_corr                    = .false.
        Logical ::  L_sol_bearing                  = .false.
!
!
        Real, Pointer :: solar_out_toa(:, :)                  => NULL()
!                          Reflected SW flux at the top
!                          of the atmosphere
        Real, Pointer :: solar_out_clear(:, :)                => NULL()
!                          Reflected clear-sky SW flux
!                          at the top of the atmosphere:
!                          this is calculated at all points
!                          omitting cloud (Method II)
        Real, Pointer :: surface_down_flux(:, :)              => NULL()
!                          Downward SW flux at the surface
!                          (not net)
        Real, Pointer :: surf_down_clr(:, :)                  => NULL()
!                          Clear-sky downward SW flux at the
!                          surface (Method II)
        Real, Pointer :: surf_up_clr(:, :)                    => NULL()
!                          Clear-sky upward SW flux at the
!                          surface (Method II)
        Real, Pointer :: net_flux_trop(:, :)                  => NULL()
!                          Net downward flux at the tropopause
        Real, Pointer :: up_flux_trop(:, :)                   => NULL()
!                          Actual upward flux at the tropopause
        Real, Pointer :: flux_direct(:, :, :)
!                          Direct downward SW flux
        Real, Pointer :: flux_diffuse(:, :, :)
!                          Diffuse downward SW flux
        Real, Pointer :: clear_hr(:, :, :)
!                          Clear-sky heating rates, calculated
!                          by ignoring all clouds (Method II):
!                          these are not necessarily the same
!                          as the heating rates in the cloud-
!                          free parts of a grid-box
        Real, Pointer :: re_strat(:, :, :)                    => NULL()
!                          The weighted effective radius in
!                          stratiform clouds multiplied by 10^6
!                          to avoid packing problems
        Real, Pointer :: wgt_strat(:, :, :)                   => NULL()
!                          The weighting factor for re_strat
!                          and lwp_strat
        Real, Pointer :: lwp_strat(:, :, :)                   => NULL()
!                          The liquid water path in stratiform
!                          cloud weighted by wgt_strat
        Real, Pointer :: re_conv(:, :, :)                     => NULL()
!                          The weighted effective radius in
!                          Convective clouds multiplied by 10^6
!                          to avoid packing problems
        Real, Pointer :: wgt_conv(:, :, :)                    => NULL()
!                          The weighting factor for re_conv
        Real, Pointer :: ntot_diag(:, :, :)                   => NULL()
!                          The number concentration of droplets
!                          multiplied by stratiform weighting factor
        Real, Pointer :: strat_lwc_diag(:, :, :)              => NULL()
!                          The liquid water content of stratiform
!                          clouds multiplied by stratiform weighting
!                          factor
        Real, Pointer :: so4_ccn_diag(:, :, :)                => NULL()
!                          The mass concentration of SO4 CCN multiplied
!                          by the conditional sampling weight
        Real, Pointer :: cond_samp_wgt(:, :, :)               => NULL()
!                          The conditional sampling weight for
!                          so4_ccn_diag
        Real, Pointer :: weighted_re(:, :)                    => NULL()
!                          The effective radius as seen from space
!                          multiplied by an appropriate weight
        Real, Pointer :: sum_weight_re(:, :)                  => NULL()
!                          The weighting factor for the effective
!                          radius as viewed from space
        Real, Pointer :: weighted_warm_re(:, :)               => NULL()
!                          The effective radius as seen from space
!                          for warm clouds only (T>273K) multiplied
!                          by an appropriate weight
        Real, Pointer :: sum_weight_warm_re(:, :)             => NULL()
!                          The weighting factor for the warm-cloud-
!                          only effective radius as viewed from space
        Real, Pointer :: Nc_diag(:, :)                        => NULL()
!                          The column-integrated cloud droplet number
!                          multiplied by an appropriate weight
        Real, Pointer :: Nc_weight(:, :)                      => NULL()
!                          The weighting factor for the column-
!                          integrated cloud droplet number
        Real, Pointer :: FlxSolBelow690nmSurf(:, :)           => NULL()
!                          The grid-box mean flux below 690 nm
!                          into the solid surface
        Real, Pointer :: FlxSeaBelow690nmSurf(:, :)           => NULL()
!                          The grid-box mean flux below 690 nm
!                          into the sea surface
        Real, pointer :: cloud_extinction(:, :, :)            => NULL()
!                           Mean extinction coefficient in clouds,
!                           weighted by the cloud amount and the
!                           clear sky flux
        Real, pointer :: inhom_cld(:, :, :)
! sza                       Cloud inhomgeneous parameter
        Real, pointer :: cloud_weight_extinction(:, :, :)     => NULL()
!                           Weighting factor for extinction in clouds :
!                           the product of the cloud amount and the
!                           clear-sky direct flux.
        Real, pointer :: ls_cloud_extinction(:, :, :)         => NULL()
!                           Mean extinction coefficient in layer clouds,
!                           weighted by the cloud amount and the
!                           clear sky flux
        Real, pointer :: ls_cloud_weight_extinction(:, :, :)  => NULL()
!                           Weighting factor for extinction in layer
!                           clouds : the product of the cloud amount
!                           and the clear-sky direct flux.
        Real, pointer :: cnv_cloud_extinction(:, :, :)        => NULL()
!                           Mean extinction coefficient in conv. clouds,
!                           weighted by the cloud amount and the
!                           clear sky flux
        Real, pointer :: cnv_cloud_weight_extinction(:, :, :) => NULL()
!                           Weighting factor for extinction in conv.
!                           clouds : the product of the cloud amount
!                           and the clear-sky direct flux.
        Real, Pointer :: orog_corr(:, :)                      => NULL()
!                          Correction factor for the direct solar flux
!                          reaching the surface for sloping terrain.
        Real, Pointer :: sol_bearing(:, :)                    => NULL()
!                          Mean local bearing of the sun over the time
!                          step, in radians clockwise from grid north.
        Real, pointer :: flxdirparsurf(:, :)                  => NULL()
!                           The direct component of the grid box mean
!                           flux below 690nm (PAR) at the surface
!
      END TYPE StrSWDiag
! ----------------------------------------------------------------------
End Module swrdiag_mod
! ----------------------------------------------------------------------
