#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Calculate and output some SCM diagnostics

      subroutine dgnstcs_imp_ctl2(                                      &
      ! IN
     &     row_length,rows,rhc_row_length,rhc_rows,                     &
     &     cloud_fraction_method,overlap_ice_liquid,                    &
     &     ice_fraction_method,ctt_weight,t_weight,                     &
     &     qsat_fixed,sub_cld,                                          &
     &     x1i,x1ic,x1r,x2r,x4r,l_psd,ai,bi,aic,bic,                    &
     &     lsp_ei,lsp_fi,lsp_eic,lsp_fic,                               &
     &     model_levels,wet_model_levels,bl_levels,                     &
     &     cloud_levels,off_x,off_y,eta_theta_levels,                   &
     &     rhcrit,L_murk,combined_cloud,nclds,                          &
     &     cumulus,ntml,L_eacf,plsp,                                    &
     &     conv_rain,conv_snow,ls_rain,ls_snow,                         &
     &     n_rows,R_u,R_v,u_incr_diagnostic,v_incr_diagnostic,          &
     &     bl_depth,bl_top,bl_type_1,bl_type_2,bl_type_3,bl_type_4,     &
     &     bl_type_5,bl_type_6,bl_type_7,bl_alltypes,                   &
     &     fqt,ftl,T,q,cf,cfl,cff,                                      &
     &     T_earliest,q_earliest,qcl_earliest,qcf_earliest,             &
     &     cf_earliest,cfl_earliest,cff_earliest,                       &
     &     p_theta_levels,LWP,IWP,z0m,z0h_eff_gb,z0m_eff_gb,            &
     &     sea_ice_htf,sice_mlt_htf,taux,tauy,area_cloud_fraction,      &
     &     p_star,u10m,v10m,latent_heat,cca_2d,rho1,qcl,qcf,            &
     &     surf_ht_flux_sice,rhcpt,surf_ht_flux_gb,aerosol,             &
     &     nSCMDpkgs,L_SCMDiags,                                        &
      ! INOUT
     &     q1p5m,t1p5m)
      implicit none

! Description:
      ! Much of what is in here is (regrettably) duplication of code
      ! in diagnostics_lscld and diagnostics_bl. These routines can't
      ! be called in the SCM because of explicit STASH calls, but we
      ! still want the diagnostics.

! Owner: Luke Jones

! Code Description:
      ! Language: Fortran77 with some bits of Fortran90

      ! ALL PARAMETERS ARE INTENT IN APART FROM T1P5M AND Q1P5M
      ! WHICH ARE INOUT.

      integer row_length,rows,rhc_row_length,rhc_rows,model_levels      &
     &     ,wet_model_levels,bl_levels,cloud_levels                     &
     &     ,off_x,off_y                                                 &
                         ! MPP associated variables, expect both to be
                         ! zero since this is the SCM
     &     ,nclds                                                       &
     &     ,n_rows

      logical                                                           &
     &     L_murk                                                       &
                           ! Switch for (visibility) aerosol
     &     ,L_eacf                                                      &

     &     ,cumulus(row_length,rows)                                    &
                                     ! Logical indicator of convection
     &     ,ntml(row_length,rows) ! Height of diagnosed BL top

      real eta_theta_levels(0:model_levels),rhcrit(wet_model_levels)    &
     &     ,plsp(row_length,rows)                                       &
     &     ,conv_rain(row_length,rows)                                  &
     &     ,conv_snow(row_length,rows)                                  &
     &     ,ls_rain(row_length,rows)                                    &
     &     ,ls_snow(row_length,rows)                                    &
     &     ,R_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          model_levels)                                           &
     &     ,R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,         &
     &          model_levels)                                           &
     &     ,u_incr_diagnostic(row_length, rows, model_levels)           &
     &     ,v_incr_diagnostic(row_length, n_rows, model_levels)         &
     &     ,bl_depth(row_length,rows)                                   &
     &     ,bl_top(row_length,rows)                                     &
     &     ,bl_type_1(row_length,rows)                                  &
                                      ! IN Indicator set to 1.0 if
!                                     ! stable b.l. diagnosed, 0.0
!                                     ! otherwise.
     &     ,bl_type_2(row_length,rows)                                  &
                                      ! IN Indicator set to 1.0 if Sc
!                                     ! over stable surface layer
!                                     ! diagnosed, 0.0 otherwise.
     &     ,bl_type_3(row_length,rows)                                  &
                                      ! IN Indicator set to 1.0 if
!                                     ! well mixed b.l. diagnosed,
!                                     ! 0.0 otherwise.
     &     ,bl_type_4(row_length,rows)                                  &
                                      ! IN Indicator set to 1.0 if
!                                     ! decoupled Sc layer (not over
!                                     ! cumulus) diagnosed,
!                                     ! 0.0 otherwise.
     &     ,bl_type_5(row_length,rows)                                  &
                                      ! IN Indicator set to 1.0 if
!                                     ! decoupled Sc layer over cumulus
!                                     ! diagnosed, 0.0 otherwise.
     &     ,bl_type_6(row_length,rows)                                  &
                                      ! IN Indicator set to 1.0 if a
!                                     ! cumulus capped b.l. diagnosed,
!                                     ! 0.0 otherwise.
     &     ,bl_type_7(row_length,rows)                                  &
                                      ! IN Indicator set to 1.0 if a
!                                     ! shear-dominated b.l. diagnosed,
!                                     ! 0.0 otherwise.
     &     ,bl_alltypes(row_length,rows)                                &
     &     ,fqt(row_length,rows,bl_levels)                              &
     &     ,ftl(row_length,rows,bl_levels)                              &
     &     ,z0m(row_length,rows)                                        &
                                       ! IN Roughness length for mom
     &     ,z0m_eff_gb(row_length,rows)                                 &
                                       ! IN Orographic roughness length
     &     ,z0h_eff_gb(row_length,rows)                                 &
                                       ! IN Roughness length for heat
! variables needed for PC2 diagnostics
     &     ,T(row_length, rows, model_levels)                           &
     &     ,q(row_length, rows, wet_model_levels)                       &
     &     ,cf(row_length, rows, wet_model_levels)                      &
     &     ,cfl(row_length, rows, wet_model_levels)                     &
     &     ,cff(row_length, rows, wet_model_levels)                     &
! _earliest arrays contain fields at start of imp_ctl
     &     ,T_earliest(row_length, rows, model_levels)                  &
     &     ,q_earliest(row_length, rows, wet_model_levels)              &
     &     ,qcl_earliest(row_length, rows, wet_model_levels)            &
     &     ,qcf_earliest(row_length, rows, wet_model_levels)            &
! _earliest values contain the values of temperature, water contents
! and cloud fractions before the boundary layer call.
     &     ,cf_earliest(row_length, rows, wet_model_levels)             &
     &     ,cfl_earliest(row_length, rows, wet_model_levels)            &
     &     ,cff_earliest(row_length, rows, wet_model_levels)            &
     &     ,p_theta_levels(row_length, rows, model_levels)              &
     &     ,LWP(row_length, rows)                                       &
                                   ! Liquid water path  (kg/m2)
     &     ,IWP(row_length, rows)                                       &
                                   ! Ice water path     (kg/m2)
     &     ,sea_ice_htf(row_length,rows)                                &
     &     ,sice_mlt_htf(row_length,rows)                               &
     &     ,taux(row_length,rows,bl_levels)                             &
     &     ,tauy(row_length,rows,bl_levels)                             &
     &     ,area_cloud_fraction(row_length,rows,wet_model_levels)       &
     &     ,p_star(row_length,rows)                                     &
     &     ,u10m(row_length,rows)                                       &
     &     ,v10m(row_length,rows)                                       &
     &     ,latent_heat(row_length,rows)                                &
     &     ,cca_2d(row_length,rows)                                     &
     &     ,rho1(row_length,rows)                                       &
                                     ! Air density at level 1 / kg m-3
     &     ,qcl(row_length,rows,wet_model_levels)                       &
     &     ,qcf(row_length,rows,wet_model_levels)                       &
     &     ,q1p5m(row_length,rows)                                      &
                                   ! INOUT
     &     ,t1p5m(row_length,rows)                                      &
                                   ! INOUT
     &     ,surf_ht_flux_sice(row_length,rows)                          &
     &     ,surf_ht_flux_gb(row_length,rows)                            &
     &     ,rhcpt(rhc_row_length,rhc_rows,wet_model_levels)             &
     &     ,aerosol(1-off_x:row_length+off_x,1-off_y:rows+off_y,        &
     &                                                    model_levels) &
     &     ,combined_cloud(row_length,rows,wet_model_levels)

      Integer                                                           &
                        !, Intent(IN)
     &      cloud_fraction_method                                       &
                                  ! Method for calculating
                                  ! total cloud fraction
     &     ,ice_fraction_method ! Method for calculating ice cloud frac.

      Real                                                              &
                        !, Intent(IN)
     &      overlap_ice_liquid                                          &
                                ! Overlap between ice and liquid phases
     &     ,ctt_weight                                                  &
                                ! Weighting of cloud top temperature
     &     ,t_weight                                                    &
                                ! Weighting of local temperature
     &     ,qsat_fixed                                                  &
                                ! Fixed value of saturation humidity
     &     ,sub_cld             ! Scaling parameter

      Real                                                              &
                  !, Intent(IN)
     &  x1i                                                             &
                  ! Intercept of aggregate size distribution
     &, x1ic                                                            &
                  ! Intercept of crystal size distribution
     &, x1r                                                             &
                  ! Intercept of raindrop size distribution
     &, x2r                                                             &
                  ! Scaling parameter of raindrop size distribution
     &, x4r                                                             &
                  ! Shape parameter of raindrop size distribution
     &, ai, bi, aic, bic                                                &
                  ! Ice mass diameter relationships m(D) = ai D^bi
     &, lsp_ei, lsp_fi, lsp_eic, lsp_fic                                
                  ! Ice particle Best-Reynolds number relationships
                  ! Re(D) =LSP_EI Be^LSP_FI

! Logicals for SCM Diagnostics packages
      Integer                                                           &
     &  nSCMDpkgs             ! No of SCM diagnostics packages

      Logical                                                           &
     &  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages
      Logical  L_psd          ! Use generic ice particle size distn.
!
! Include parameters necessary for calls to SCMoutput...
#include "s_scmop.h"

      ! Local variables...
      Character*(*), Parameter :: RoutineName = 'dgnstcs_imp_ctl2'

      integer i,j,k,error_code_ignored                                  &
     &, kinvert               ! Vertical index for inverted arrays

      ! Work array
      real interp_data(row_length*rows*model_levels)                    &
     &    ,work_2d(row_length,rows)                                     &
     &    ,work_3d_w(row_length,rows,wet_model_levels)                  &
     &    ,work_3d(row_length,rows,model_levels)

      ! Reference numbers for cloud schemes
#include "clschm3a.h"

      ! Global variables
#include "c_r_cp.h"
#include "c_lheat.h"

      ! Local parameters and other physical constants
      Real LCRCP                        ! Derived parameter.
      Parameter ( LCRCP=LC/CP )         ! Lat ht of condensation/Cp.

      ! Variables used to calculate visibility stuff
#include "c_pi.h"
! vis_thresh and n_vis_thresh are in here:
#include "c_visbty.h"
      real                                                              &
     &      Beta_LS_Rain(row_length,rows)                               &
                                          ! Scattering in LS Rain.
     &     ,Beta_LS_Snow(row_length,rows)                               &
                                          ! Scattering in LS Snow.
     &     ,Beta_C_Rain(row_length,rows)                                &
                                         ! Scattering in Conv Rain
     &     ,Beta_C_Snow(row_length,rows)                                &
                                         ! Scattering in Conv Snow
     &     ,Vis_Threshold(row_length,rows,1,n_vis_thresh)               &
                                                          ! FOG_FR works
                                ! for n levels, we want 1
     &     ,PVis(row_length,rows,n_vis_thresh)

      real  v_1p5m(row_length,rows)                                     &
                                        ! Visibility overall at 1.5m
     &     ,v_no_precip(row_length,rows)                                &
                                        ! Visibility without precip.
     &     ,v_ls_precip(row_length,rows)                                &
                                        ! Visibility in LS Precip.
     &     ,v_c_precip(row_length,rows)                                 &
                                        ! Visibility in Conv Precip.
     &     ,v_probs(row_length,rows,1,n_vis_thresh) ! Vis probs at
                                      ! first level (decimal fraction)

      ! Total cloud amounts.
      ! (X,X,1) - max. overlap. (X,X,2) - random overlap.
      real tot_cloud(row_length,rows,2)

      ! More output diagnostics
      real layer_cloud1p5m(row_length,rows)                             &
                                            ! Layer cloud at 1.5m
                                            ! (decimal fraction)
     &     ,qcl1p5m(row_length,rows) ! Cloud liquid water content at
                                     ! 1.5m (kg per kg of air)

      ! Variables used in the calculation of the low_cloud, med_cloud
      ! and high_cloud diagnostics.
      integer, parameter :: num_cloud_types=3
      real h_split(num_cloud_types+1),cloud_bound(num_cloud_types+1)
      integer kk,level,low_bot_level,low_top_level,                     &
     &     med_bot_level,med_top_level,high_bot_level,high_top_level
! max_model_levels is in cmaxsize.h which is only required because it
! is referenced in s_vcoord.h:
#include "cmaxsize.h"
! z_top_of_model is in here:
#include "s_vcoord.h"
      real  low_cloud(row_length,rows)                                  &
     &     ,med_cloud(row_length,rows)                                  &
     &     ,high_cloud(row_length,rows)                                 &
     &     ,td1p5m(row_length,rows)                                     &
                                       ! 1.5m dewpoint temperature
     &     ,rh1p5m(row_length,rows)                                     &
                                       ! 1.5m relative humidity
     &     ,rhw1p5m(row_length,rows)                                    &
                                       ! 1.5m relative humidity over
                                       ! liquid water
     &     ,wspd10m(row_length,rows)                                    &
                                       ! 10m wind speed
     &     ,wdrn10m(row_length,rows)                                    &
                                       ! 10m wind direction
     &     ,bl_gust(row_length,rows)                                    &
                                       ! Shear wind gust
     &     ,alpha                                                       & 
                                       ! 'Angle' of wind
     &     ,qst                        ! A saturation mixing ratio

      ! Parameters
      Logical, Parameter :: l_mixing_ratio=.false. ! Use mixing rations
!
!
!-----------------------------------------------------------------------
!     SCM Large Scale Cloud Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_lscld)) Then

! To calculate low, medium and high cloud amounts, first need to set up
! variables low_bot_level, low_top_level, med_bot_level, etc. In global
! and limited-area models this is done during initialisation.
         ! Default settings for h_split as set in routine READLSTA
         ! (which is not called in the SCM)
         !   low   : middle :  high   cloud model levels
         ! (1)->(2):(2)->(3):(3)->(4)
         h_split(1) =   111.    ! ICAO 1000mb height (m)
         h_split(2) =  1949.    ! ICAO  800mb height (m)
         h_split(3) =  5574.    ! ICAO  500mb height (m)
         h_split(4) = 13608.    ! ICAO  150mb height (m)
         ! Set up cloud type boundaries for low/medium/high cloud as in
         ! routine SETCONA (which is not called in the SCM). See this
         ! routine for more comments on what's going on here.
         do kk=1,num_cloud_types+1
            level=1
            do while ((eta_theta_levels(level)*z_top_of_model           &
     &            <= h_split(kk))                                       &
     &           .and.                                                  &
     &           (level <= model_levels))
               level=level+1
            enddo
            if (level >  model_levels) then
               print*,'ERROR in ni_imp_ctl: level>model_levels',        &
     &              kk,level,model_levels
               level=model_levels
            endif
            cloud_bound(kk) = level
         enddo
         low_bot_level  = cloud_bound(1)
         low_top_level  = cloud_bound(2) - 1
         med_bot_level  = cloud_bound(2)
         med_top_level  = cloud_bound(3) - 1
         high_bot_level = cloud_bound(3)
         high_top_level = cloud_bound(4) - 1

         if (low_top_level >  cloud_levels) then
            print*,'ni_imp_ctl ERROR: no of cloud levels less than '//  &
     &           'top of low',low_top_level,cloud_levels
         endif
         if (med_top_level >  cloud_levels) then
            print*,'ni_imp_ctl ERROR: no of cloud levels less than '//  &
     &           'top of med',med_top_level,cloud_levels
         endif
         if (high_top_level >  cloud_levels) then
            print*,'ni_imp_ctl ERROR: no of cloud levels less than '//  &
     &           'top of high',high_top_level,cloud_levels
         endif

!---------------------------------------------------------------------
! In non-SCM runs the following stuff would be done in diagnostics_lscld

         ! Diagnostic 203: Low Cloud Amount
         ! Initialize cloud amount to lowest level value.
         do j=1,rows
            do i=1,row_length
               low_cloud(i,j)=area_cloud_fraction(i,j,LOW_BOT_LEVEL)
            enddo
         enddo
         ! Cloud amount is calculated under maximum overlap assumption.
         ! Limit cloud amount to maximum of 1.
         do k=low_bot_level+1,low_top_level
            do j=1,rows
               do i=1,row_length
                  low_cloud(i,j) =                                      &
     &                 max(low_cloud(i,j),area_cloud_fraction(i,j,k))
                  low_cloud(i,j) = min(low_cloud(i,j),1.0)
               enddo
            enddo
         enddo
! DEPENDS ON: scmoutput
         Call SCMoutput(low_cloud,                                      &
              'lowcld','Low cloud fraction','Fraction',                 &
              t_avg,d_sl,default_streams,'',RoutineName)

         ! Diagnostic 204: Medium Cloud Amount
         ! Initialize cloud amount to lowest level value.
         do j=1,rows
            do i=1,row_length
               med_cloud(i,j)=area_cloud_fraction(i,j,med_bot_level)
            enddo
         enddo
         ! Cloud amount is calculated under maximum overlap assumption.
         ! Limit cloud amount to maximum of 1.
         do k=med_bot_level+1,med_top_level
            do j=1,rows
               do i=1,row_length
                  med_cloud(i,j) =                                      &
     &                 max(med_cloud(i,j),area_cloud_fraction(i,j,k))
                  med_cloud(i,j) = min(med_cloud(i,j),1.0)
               enddo
            enddo
         enddo
! DEPENDS ON: scmoutput
         Call SCMoutput(med_cloud,                                      &
              'medcld','Medium cloud fraction','Fraction',              &
              t_avg,d_sl,default_streams,'',RoutineName)

         ! Diagnostic 205: High Cloud Amount
         ! Initialize cloud amount to lowest level value.
         do j=1,rows
            do i=1,row_length
               high_cloud(i,j)=area_cloud_fraction(i,j,high_bot_level)
            enddo
         enddo
         ! Cloud amount is calculated under maximum overlap assumption.
         ! Limit cloud amount to maximum of 1.
         do k=high_bot_level+1,high_top_level
            do j=1,rows
               do i=1,row_length
                  high_cloud(i,j) =                                     &
     &                 max(high_cloud(i,j),area_cloud_fraction(i,j,k))
                  high_cloud(i,j) = min(high_cloud(i,j),1.0)
               enddo
            enddo
         enddo
! DEPENDS ON: scmoutput
         Call SCMoutput(high_cloud,                                     &
              'highcld','High cloud fraction','Fraction',               &
              t_avg,d_sl,default_streams,'',RoutineName)

         ! Total Cloud Amount RANDOM Overlap
! DEPENDS ON: r2_calc_total_cloud_cover
         CALL r2_calc_total_cloud_cover(                                &
     &        row_length*rows,wet_model_levels,nclds                    &
     &        ,IP_CLOUD_MIX_RANDOM,combined_cloud(1,1,1)                &
     &        ,tot_cloud(1,1,2)                                         &
     &        ,row_length*rows,wet_model_levels)

         ! Total Cloud Amount MAX/RANDOM Overlap
! DEPENDS ON: r2_calc_total_cloud_cover
         CALL r2_calc_total_cloud_cover(                                &
     &        row_length*rows,wet_model_levels,nclds                    &
     &        ,IP_CLOUD_MIX_MAX,combined_cloud(1,1,1)                   &
     &        ,tot_cloud(1,1,1)                                         &
     &        ,row_length*rows,wet_model_levels)

! DEPENDS ON: scmoutput
         call SCMoutput(tot_cloud(1,1,2),                               &
              'tcarndm','Total cloud amount (random overlap)',          &
              'Fraction',                                               &
              t_avg,d_sl,default_streams,'',RoutineName)
! DEPENDS ON: scmoutput
         call SCMoutput(tot_cloud(1,1,1),                               &
              'tcamxrn','Total cloud amount (max. random)',             &
              'Fraction',                                               &
              t_avg,d_sl,default_streams,'',RoutineName)

! Note that stash equivalents 9,18[123] are identical to 3,18[123]
! and so done below stash 9,226
           Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                If (area_cloud_fraction(i, j, k) <= 0.0)  Then
                  work_3d_w(i, j, k) = 0.0
                Else
                  work_3d_w(i, j, k) = 1.0
                End If
              End Do ! i
            End Do ! j
          End Do ! k
!
! DEPENDS ON: scmoutput
          Call SCMoutput(work_3d_w,                                     &
               'lyrcldfreq','Layer cloud frequency indicator',          &
               'Indicator',                                             &
               t_avg,d_wet,default_streams,'',RoutineName)

!         Stash 9,228
! DEPENDS ON: scmoutput
          Call SCMoutput(rhcpt,                                         &
               'rhcpt','Critical relative humidity','%',                &
               t_avg,d_wet,default_streams,'',RoutineName)

!         Stash 9,229
          Do k = 1, wet_model_levels
! DEPENDS ON: qsat
            Call QSAT(work_2d,T(1,1,K),p_theta_levels(1,1,k),           &
     &                row_length*rows)
            Do j = 1, rows
              Do i = 1, row_length
!               ! relative humidity in per cent
                work_3d_w(i, j, k) = q(i, j, k) / work_2d(i, j) *100.0
!               ! Supersaturation (>100%) can occur with mixed phase
!               ! scheme but negative humidity is removed from the
!               ! diagnostic:
                If ( work_3d_w(i, j, k) < 0.0) Then
                  work_3d_w(i, j, k) = 0.0
                End If
              End Do ! i
            End Do ! j
          End Do ! k
!
! DEPENDS ON: scmoutput
          Call SCMoutput(work_3d_w,                                     &
               'rhaftcld','Relative humidity after main cloud','%',     &
               t_avg,d_wet,default_streams,'',RoutineName)

!         Stash 9,231
          Do k = 1, wet_model_levels
!           NB: Convention in Sect 70 (Radiation) is to invert levels,
!           1 at top. Combined_cloud is calculated this way but
!           re-inverted for STASH.

            kinvert = wet_model_levels+1-k

            Do j = 1, rows
              Do i = 1, row_length
                work_3d_w(i, j, k) = combined_cloud(i,j,kinvert)
              End Do ! i
            End Do ! j
          End Do ! k
!
! DEPENDS ON: scmoutput
          Call SCMoutput(work_3d_w,                                     &
               'combca_ls','Combined cloud amount in each layer',       &
               'Fraction',                                              &
               t_avg,d_wet,default_streams,'',RoutineName)
!

! End of stuff that would be done in diagnostics_lscld               !
!---------------------------------------------------------------------


      End If ! L_SCMDiags(SCMDiag_lscld)

!
!-----------------------------------------------------------------------
!     SCM Surface Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_surf)) Then

!       Stash 3,255
! DEPENDS ON: scmoutput
        Call SCMoutput(q1p5m,                                           &
             'qt1p5m','1.5m total water kg water/kg air',               &
             'kgH20/kgAIR',                                             &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 3,254
! DEPENDS ON: scmoutput
        Call SCMoutput(t1p5m,                                           &
             'tl1p5m','1.5m liquid temperature','K',                    &
             t_avg,d_sl,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_surf)
!
!-----------------------------------------------------------------------
!     SCM Boundary Layer OR Large Scale Cloud Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_bl) .OR. L_SCMDiags(SCMDiag_lscld)) Then

! DEPENDS ON: scmoutput
        Call SCMoutput(lwp,                                             &
             'LWP','Liquid water path','kg/m2',                         &
             t_avg,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(iwp,                                             &
             'IWP','Ice water path','kg/m2',                            &
             t_avg,d_sl,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_bl) .OR. L_SCMDiags(SCMDiag_lscld)


!---------------------------------------------------------------------
! In non-SCM runs the following stuff would be done in diagnostics_bl

      ! This call to ls_cld converts 1.5m TL and QT to T and q,
      ! assuming p_star and rhcrit at model level 1 approximates to
      ! 1.5m. This should be done even if ssfm=.false. It also
      ! calculates layer_cloud1p5m and qcl1p5m.
! DEPENDS ON: ls_cld
      CALL ls_cld(                                                      &
     &     p_star,rhcpt,1,bl_levels,row_length,rows                     &
     &     ,rhc_row_length,rhc_rows                                     &
     &     ,cloud_fraction_method,overlap_ice_liquid                    &
     &     ,ice_fraction_method,ctt_weight,t_weight,qsat_fixed,sub_cld  &
     &     ,ntml,cumulus,L_eacf                                         &
                                            ! in
     &     ,l_mixing_ratio                                              &
                                            ! in
     &     ,t1p5m                                                       &
                                            ! in/out
     &     ,layer_cloud1p5m                                             &
                                            ! out: CF Cld frctn
     &     ,q1p5m                                                       &
                                            ! in/out
     &     ,qcf                                                         &
                                            ! in
     &     ,qcl1p5m                                                     &
                                            ! out: QCL Cld liq H2O cntnt
     &     ,interp_data(2*row_length*rows+1)                            &
                                            ! out: CFL Liq cld frctn
     &     ,interp_data(3*row_length*rows+1)                            &
                                            ! out: CFF Frzn cld frctn
     &     ,error_code_ignored)
!
!-----------------------------------------------------------------------
!     SCM Large Scale Cloud OR Surface Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_lscld) .OR. L_SCMDiags(SCMDiag_surf)) Then

! DEPENDS ON: scmoutput
      call SCMoutput(layer_cloud1p5m,                                   &
           'lca1p5m','1.5m layer cloud amount','Fraction',              &
           t_avg,d_sl,default_streams,'',RoutineName)
! DEPENDS ON: scmoutput
      call SCMoutput(qcl1p5m,                                           &
           'qcl1p5m','1.5m cloud water','kg/kg',                        &
           t_avg,d_sl,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_lscld) .OR. L_SCMDiags(SCMDiag_surf)

!
!-----------------------------------------------------------------------
!     SCM Surface Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_surf)) Then

         ! 1.5m Fog Fraction and Mist Fraction
         do i = 1,row_length
            do j = 1,rows
               do k=1,n_vis_thresh
                  Vis_Threshold(i,j,1,k)=vis_thresh(k)
               enddo
            enddo
         enddo
! DEPENDS ON: fog_fr
         CALL fog_fr(p_star,rhcrit,1,                                   &
     &        row_length*rows,                                          &
     &        t1p5m,aerosol,L_murk,                                     &
     &        q1p5m,qcl,qcf,                                            &
     &        vis_threshold,                                            &
     &        PVis,n_vis_thresh)

         ! Calculate scattering coefficients due to precipitation
! DEPENDS ON: beta_precip
         CALL beta_precip(ls_rain,ls_snow                               &
     &        ,conv_rain,conv_snow,qcf(1,1,1)                           &
     &        ,rho1,T1p5m,p_star                                        &
     &        ,plsp,cca_2d,.false.,.true.                               &
     &        ,row_length*rows,row_length*rows,1                        &
     &        ,x1i,x1ic,x1r,x2r,x4r,l_psd,ai,bi,aic,bic                 &
     &        ,lsp_ei,lsp_fi,lsp_eic,lsp_fic                            &
     &        ,beta_ls_rain,beta_ls_snow                                &
     &        ,beta_c_rain,beta_c_snow,error_code_ignored)

         ! Calculate screen level probability of visibility
         ! less than thresholds
! DEPENDS ON: calc_vis_prob
         CALL calc_vis_prob(p_star,                                     &
                                               ! input
     &        rhcrit,1,                                                 &
                                               ! input
     &        row_length*rows,row_length*rows,                          &
                                               ! input
     &        t1p5m,aerosol,l_murk,                                     &
                                               ! input
     &        q1p5m,qcl,qcf,                                            &
                                               ! input
     &        vis_thresh,n_vis_thresh,                                  &
                                               ! input
     &        plsp,cca_2d,.false.,                                      &
                                               ! input
     &        beta_ls_rain,beta_ls_snow,                                &
                                               ! input
     &        beta_c_rain,beta_c_snow,                                  &
                                               ! input
     &        v_probs,                                                  &
                                               ! OUTput
     &        error_code_ignored)              ! OUTput

         ! Visibility at 1.5 m including precipitation
! DEPENDS ON: visbty
         CALL visbty(                                                   &
     &         p_star,T1p5m,q1p5m,Qcl,Qcf                               &
                                                        ! input
     &        ,Aerosol,calc_prob_of_vis,RHcrit,L_murk                   &
                                                        ! input
     &        ,row_length * rows                                        &
                                                        ! input
     &        ,v_no_precip)                             ! OUTput
! DEPENDS ON: vis_precip
         CALL vis_precip(v_no_precip                                    &
                                                        ! input
     &        ,plsp,cca_2d,.false.                                      &
                                                      ! input
     &        ,Beta_ls_Rain,Beta_ls_Snow                                &
                                                      ! input
     &        ,Beta_C_Rain,Beta_C_Snow                                  &
                                                      ! input
     &        ,row_length*rows,row_length*rows,1                        &
                                                      ! input
     &        ,v_1p5m,v_ls_precip,v_C_Precip                            &
                                                      ! OUTput
     &        ,error_code_ignored)                    ! OUTput

! DEPENDS ON: scmoutput
         call SCMoutput(v_probs(1,1,1,1),                               &
              'pfog1p5m','Probability of fog at 1.5m','-',              &
              t_avg,d_sl,default_streams,'',RoutineName)
! DEPENDS ON: scmoutput
         call SCMoutput(v_probs(1,1,1,2),                               &
              'pmist1p5m','Probability of mist at 1.5m','-',            &
              t_avg,d_sl,default_streams,'',RoutineName)
! DEPENDS ON: scmoutput
         call SCMoutput(v_probs,                                        &
              'pvisthresh','Probability of vis<=threshold','-',         &
              t_avg,d_vis,default_streams,'',RoutineName)
! DEPENDS ON: scmoutput
         call SCMoutput(v_no_precip,                                    &
              'visnop1p5m','1.5m visibility outside precip','m',        &
              t_avg,d_sl,default_streams,'',RoutineName)
! DEPENDS ON: scmoutput
         call SCMoutput(v_ls_precip,                                    &
              'vislsp1p5m','1.5m visibility in LS precip','m',          &
              t_avg,d_sl,default_streams,'',RoutineName)
! DEPENDS ON: scmoutput
         call SCMoutput(v_c_precip,                                     &
              'viscp1pm5','1.5m visibility in conv precip','m',         &
              t_avg,d_sl,default_streams,'',RoutineName)
! DEPENDS ON: scmoutput
         call SCMoutput(v_1p5m,'vis1p5m','1.5m visibility','m',         &
              t_avg,d_sl,default_streams,'',RoutineName)

         ! Calculate relative humidities at screen level
         do j=1,rows
            do i=1,row_length
! DEPENDS ON: qsat
               call QSAT(qst,t1p5m(i,j),p_star(i,j),1)
               rh1p5m(i,j)=q1p5m(i,j)/qst*100.0
! DEPENDS ON: qsat_wat
               call QSAT_WAT(qst,t1p5m(i,j),p_star(i,j),1)
               rhw1p5m(i,j)=q1p5m(i,j)/qst*100.0
            enddo
         enddo
! DEPENDS ON: scmoutput
         call SCMoutput(rh1p5m,                                         &
              'rh1p5m','Relative humidity at 1.5m','%',                 &
              t_avg,d_sl,default_streams,'',RoutineName)
! DEPENDS ON: scmoutput
         call SCMoutput(rhw1p5m,                                        &
              'rhw1p5m','Relative humidity wrt H2O at 1.5m','W/m2',     &
              t_avg,d_sl,default_streams,'',RoutineName)

         ! Calculate dewpoint
! DEPENDS ON: dewpnt
         call dewpnt(q1p5m,p_star,t1p5m,row_length*rows,td1p5m)
! DEPENDS ON: scmoutput
         call SCMoutput(td1p5m,                                         &
              'td1p5m','1.5m dewpoint temperature','W/m2',              &
              t_avg,d_sl,default_streams,'',RoutineName)

         ! Calculate the magnitude and direction (by meteorological
         ! convention) of the wind at 10m from the 10m u and v
         ! components
         do j=1,rows
            do i=1,row_length
               wspd10m(i,j)=sqrt(u10m(i,j)**2+v10m(i,j)**2)
               if (wspd10m(i,j) == 0.0) then
                  wdrn10m(i,j)=0.0
               else
                  alpha=asin(v10m(i,j)/wspd10m(i,j))*                   &
     &                 recip_pi_over_180
                  if (u10m(i,j) >  0.0) then
                     wdrn10m(i,j)=270-alpha
                  else
                     wdrn10m(i,j)=90+alpha
                  endif
               endif
            enddo
         enddo
! DEPENDS ON: scmoutput
         call SCMoutput(wspd10m,'wspd10m','10m wind speed','m/s',       &
              t_avg,d_sl,default_streams,'',RoutineName)
! DEPENDS ON: scmoutput
         call SCMoutput(wdrn10m,                                        &
              'wdrn10m','10m wind direction','degs',                    &
              t_avg,d_sl,default_streams,'',RoutineName)

         ! Boundary layer gust diagnostic
         do j=1,rows
            do i=1,row_length
               bl_gust(i,j)=wspd10m(i,j)+2.0*2.5                        &
     &              *sqrt(sqrt(taux(i,j,1)**2+tauy(i,j,1)**2))
            enddo
         enddo
! DEPENDS ON: scmoutput
         Call SCMoutput(bl_gust,                                        &
              'gust10m','10m Gust','m/s',                               &
              t_avg,d_sl,default_streams,'',RoutineName)


      End If ! L_SCMDiags(SCMDiag_surf)

!
!-----------------------------------------------------------------------
!     SCM Sea Points Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_sea)) Then

! DEPENDS ON: scmoutput
        Call SCMoutput(sice_mlt_htf,                                    &
             'sice_mlt_htf','Heat flux due to melting sea ice','W/m2',  &
             t_avg,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(sea_ice_htf,                                     &
             'sea_ice_htf','Heat flux through sea ice','W/m2',          &
             t_avg,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(surf_ht_flux_sice,                               &
             'surf_ht_flux_si','Net downward heat flux into sea frctn', &
             'W/m2',                                                    &
             t_avg,d_sl,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_sea)
!
!-----------------------------------------------------------------------
!     SCM Surface Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_surf)) Then

! DEPENDS ON: scmoutput
        Call SCMoutput(latent_heat,                                     &
             'lat_ht','Surface latent heat flux','W/m2',                &
             t_avg,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(surf_ht_flux_gb,                                 &
             'surf_ht_flux','Net downward heat flux at surface','W/m2', &
             t_avg,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(z0m,                                             &
             'z0m','Roughness length for momentum','m',                 &
             t_avg,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(z0h_eff_gb,                                      &
             'z0h','Effective roughness length for heat','m',           &
             t_avg,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(z0m_eff_gb,                                      &
             'z0m_eff','Effective roughness length for momentum','m',   &
             t_avg,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(q1p5m,                                           &
             'q1p5m','1.5m specific humidity','kg/kg',                  &
             t_avg,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(t1p5m,                                           &
             't1p5m','1.5m temperature','K',                            &
             t_avg,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(t1p5m,                                           &
             't1p5m_max','Max 1.5m temperature','K',                    &
             t_max,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(t1p5m,                                           &
             't1p5m_min','Min 1.5m temperature','K',                    &
             t_min,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(u10m,                                            &
             'u10m','Zonal 10m wind','m/s',                             &
             t_avg,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(v10m,                                            &
             'v10m','Meridional 10m wind','m/s',                        &
             t_avg,d_sl,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_surf)
!
!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_bl)) Then

!       Stash 3,025
! DEPENDS ON: scmoutput
        Call SCMoutput(bl_depth,                                        &
             'bl_depth','Boundary layer depth from B.layer','m',        &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Since model level fields are on rho_levels and the surface
!       is a theta_level, these are output as separate diagnostics

!       Stash 3,217
        Do j = 1,rows
          Do i = 1,row_length
            work_2d(i,j) = ftl(i,j,1)
          End Do ! i
        End Do ! j

! DEPENDS ON: scmoutput
        Call SCMoutput(work_2d,                                         &
             'ftl_surf',                                                &
             'Surface sensible heat flux from B.layer','W/m2',          &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 3,223
        Do j = 1,rows
          Do i = 1,row_length
            work_2d(i,j) = fqt(i,j,1)
          End Do ! i
        End Do ! j

! DEPENDS ON: scmoutput
        Call SCMoutput(work_2d,                                         &
             'fqt_surf',                                                &
             'Surface sensible moisture flux from B.layer','kg/m2/s',   &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 3,222
! DEPENDS ON: scmoutput
        Call SCMoutput(fqt,                                             &
             'fqt_bl','Sensible moisture flux','kg/m2/s',               &
             t_avg,d_bl,default_streams,'',RoutineName)

!       Stash 3,304
! DEPENDS ON: scmoutput
        Call SCMoutput(bl_top,                                          &
             'bl_top','Boundary layer top','m',                         &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 3,305
! DEPENDS ON: scmoutput
        Call SCMoutput(bl_type_1,                                       &
             'bl_type_1','Boundary layer type: stable','Indicator',     &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 3,306
! DEPENDS ON: scmoutput
        Call SCMoutput(bl_type_2,                                       &
             'bl_type_2','Boundary layer type: Sc over stable',         &
             'Indicator',                                               &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 3,307
! DEPENDS ON: scmoutput
        Call SCMoutput(bl_type_3,                                       &
             'bl_type_3','Boundary layer type: well mixed','Indicator', &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 3,308
! DEPENDS ON: scmoutput
        Call SCMoutput(bl_type_4,                                       &
             'bl_type_4','Boundary layer type: decoup Sc not over Cu',  &
             'Indicator',                                               &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 3,309
! DEPENDS ON: scmoutput
        Call SCMoutput(bl_type_5,                                       &
             'bl_type_5','Boundary layer type: decoup Sc over Cu',      &
             'Indicator',                                               &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 3,310
! DEPENDS ON: scmoutput
        Call SCMoutput(bl_type_6,                                       &
             'bl_type_6','Boundary layer type: cumulus capped',         &
             'Indicator',                                               &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 3,340
! DEPENDS ON: scmoutput
        Call SCMoutput(bl_type_7,                                       &
             'bl_type_7','Boundary layer type: shear driven',           &
             'Indicator',                                               &
             t_avg,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(bl_alltypes,                                     &
             'bl_alltypes','Boundary layer types','',                   &
             t_inst,d_sl,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_bl)
!
!-----------------------------------------------------------------------
!     SCM Boundary Layer OR Increments Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_bl) .OR. L_SCMDiags(SCMDiag_incs)) Then

!       Stash 3,181 and 9,181
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              work_3d(i,j,k) = T(i,j,k) - T_earliest(i,j,k)
            End Do ! i
          End Do ! j
        End Do ! k

! DEPENDS ON: scmoutput
        Call SCMoutput(work_3d,                                         &
             'dt_pc2blls','Temperature inc PC2+bdy layer+ls cld','K',   &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 3,182 and 9,182
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              work_3d_w(i,j,k) = q(i,j,k) - q_earliest(i,j,k)
            End Do ! i
          End Do ! j
        End Do ! k

! DEPENDS ON: scmoutput
        Call SCMoutput(work_3d_w,                                       &
             'dq_pc2blls',                                              &
             'Specific humidity inc PC2+bdy layer+ls cld','kg/kg',      &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 3,183 and 9,183
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              work_3d_w(i,j,k) = qcl(i,j,k) - qcl_earliest(i,j,k)
            End Do ! i
          End Do ! j
        End Do ! k

! DEPENDS ON: scmoutput
        Call SCMoutput(work_3d_w,                                       &
             'dqcl_pc2blls',                                            &
             'QCL increment PC2+bdy layer+ls cld','kg/kg',              &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 3,184
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              work_3d_w(i,j,k) = qcf(i,j,k) - qcf_earliest(i,j,k)
            End Do ! i
          End Do ! j
        End Do ! k

! DEPENDS ON: scmoutput
        Call SCMoutput(work_3d_w,                                       &
             'dqcf_pc2blls',                                            &
             'QCF increment PC2+bdy layer+ls cld','kg/kg',              &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 3,185
! DEPENDS ON: scmoutput
        Call SCMoutput(u_incr_diagnostic,                               &
           'du_bl','U wind increment bdy layer','m/s',                  &
           t_avg,d_all,default_streams,'',RoutineName)

!       Stash 3,186
! DEPENDS ON: scmoutput
        Call SCMoutput(v_incr_diagnostic,                               &
           'dv_bl','V wind increment bdy layer','m/s',                  &
           t_avg,d_all,default_streams,'',RoutineName)

!       Stash 3,192
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              work_3d_w(i,j,k) = cf(i,j,k) -   cf_earliest(i,j,k)
            End Do ! i
          End Do ! j
        End Do ! k

! DEPENDS ON: scmoutput
        Call SCMoutput(work_3d_w,                                       &
             'dbcf_bl','Bulk cloud fraction increment bdy layer',       &
             'Fraction',                                                &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 3,193
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              work_3d_w(i,j,k) = cfl(i,j,k) - cfl_earliest(i,j,k)
            End Do ! i
          End Do ! j
        End Do ! k

! DEPENDS ON: scmoutput
        Call SCMoutput(work_3d_w,                                       &
             'dcfl_bl','Liquid cloud fraction increment bdy layer',     &
             'Fraction',                                                &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 3,194
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              work_3d_w(i,j,k) = cff(i,j,k) - cff_earliest(i,j,k)
            End Do ! i
          End Do ! j
        End Do ! k

! DEPENDS ON: scmoutput
        Call SCMoutput(work_3d_w,                                       &
             'dcff_bl','Frozen cloud fraction increment bdy layer',     &
             'Fraction',                                                &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 3,189
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              work_3d(i,j,k) = T(i,j,k) -   T_earliest(i,j,k)           &
     &             - LCRCP * ( qcl(i,j,k) - qcl_earliest(i,j,k) )
            End Do ! i
          End Do ! j
        End Do ! k

        If (model_levels > wet_model_levels)  Then
          Do k = (wet_model_levels + 1), model_levels
            Do j = 1, rows
              Do i = 1, row_length
                work_3d(i,j,k) = T(i,j,k) - T_earliest(i,j,k)
              End Do ! i
            End Do ! j
          End Do ! k
        End If ! model_levels > wet_model_levels

! DEPENDS ON: scmoutput
        Call SCMoutput(work_3d,                                         &
             'dlwt_bl','Liquid water temp increment bdy layer','K',     &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 3,190
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              work_3d_w(i,j,k)  =   q(i,j,k) -   q_earliest(i,j,k)      &
     &                          + qcl(i,j,k) - qcl_earliest(i,j,k)
            End Do ! i
          End Do ! j
        End Do ! k

! DEPENDS ON: scmoutput
        Call SCMoutput(work_3d_w,                                       &
             'tqinc_bl','Total (liquid) water increment bdy layer',     &
             'kg/kg',                                                   &
             t_avg,d_all,default_streams,'',RoutineName)
!
      End If ! L_SCMDiags(SCMDiag_bl) .OR. L_SCMDiags(SCMDiag_incs)
!

! End of stuff that would be done in diagnostics_bl
!---------------------------------------------------------------------

      return
      END SUBROUTINE dgnstcs_imp_ctl2
#endif
