#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Calculate and output some SCM diagnostics

      Subroutine dgnstcs_scm_main(                                            &
        row_length, rows, land_points, model_levels,                          &
        wet_model_levels, sm_levels, st_levels, ntype,                        &
        r_theta_levels, r_rho_levels, rho, timestep, u, v, T, theta, q,       &
        qcl,qcf, layer_cloud, cca, ccw, t_deep_soil, p_star, tstar, smc,      &
        canopy_gb, snodep, zh, z0msea, smcl, sthu, sthf,                      &
        gs, lw_incs, surf_radflux, photosynth_act_rad, tstar_tile,            &
        aerosol, p_theta_levels, p, iccb, icct,                               &
        w, w_adv, area_cloud_fraction, bulk_cloud_fraction,                   &
        cloud_fraction_liquid, cloud_fraction_frozen, cclwp,                  &
!       EAK
        ftl_tile_cab, ftl_cab, le_tile_cab, le_cab,                           &
        tstar_tile_cab, tstar_cab, smcl_cab, tsoil_cab,                       &
        ustar_cab, surf_htf_cab, snow_rho1l, tot_alb,                         &
        u_s_cab, ch_cab, cd_cab,                                              &
        cd,ch, sw_down, radnet_tile, snow_tile, snow_tmp3l,                   &
        snow_rho3l, snow_depth3l, lw_down,                                    &
        nSCMDpkgs, L_SCMDiags)

      Use cv_cntl_mod, Only:                                                  &
          lcv_3d_ccw, lcv_ccrad

      Implicit None

! Description:
      ! Essentially just a lot of calls to SCMoutput
! Owner: R Wong

! Code Description:
      ! Language: Fortran77 with some bits of Fortran90

      ! All parameters are Intent IN

      ! Variables passed in through the parameter list...
      integer :: row_length, rows, land_points, model_levels,           &
                 wet_model_levels, sm_levels, st_levels, ntype
      real                                                              &
     &      r_theta_levels(row_length, rows, 0:model_levels)            &
     &     ,r_rho_levels(row_length,rows, model_levels)                 &
     &     ,rho(row_length, rows, model_levels)                         &
                                                ! density*r^2 (kg/m)
     &     ,u(row_length, rows,model_levels)                            &
                                             ! Zonal wind (m s^-1)
     &     ,v(row_length, rows,model_levels)                            &
                                             ! Meridional wind (m s^-1)
     &     ,w(row_length, rows, 0:model_levels)                         &
                                ! Vertical velocity (m/s)
     &     ,w_adv(row_length, rows, 0:model_levels)                     &
                                ! Advective vertical velocity (m/s)
     &     ,t(row_length, rows,model_levels)                            &
                                             ! Temperature(K)
     &     ,theta(row_length, rows,model_levels)                        &
                                ! Potential temperature (K)
     &     ,q(row_length, rows,wet_model_levels)                        &
                                ! Specific humidity (Kg Kg^-1)
     &     ,qcf(row_length, rows,wet_model_levels)                      &
                                ! Cloud ice content (Kg Kg^-1)
     &     ,qcl(row_length, rows,wet_model_levels)                      &
                                ! Cloud water content(Kg Kg^-1)
     &     ,layer_cloud(row_length, rows,wet_model_levels)              &
                                ! Layer cloud amount
                                ! (decimal fraction)
     &     ,cca(row_length, rows,model_levels)                          &
                                ! Convective cloud amount
     &     ,ccw(row_length, rows, wet_model_levels)                     &
                                ! Convective Cloud Water passed to
                                ! radiation scheme (kg/kg)

     &     ,area_cloud_fraction(row_length, rows, wet_model_levels)     &
                                ! Area cloud amount (decimal fraction)
     &     ,bulk_cloud_fraction(row_length, rows, wet_model_levels)     &
                                ! Cloud amount (decimal fraction)
     &     ,cloud_fraction_liquid(row_length, rows, wet_model_levels)   &
                                ! Liquid cloud amount (decimal fraction)
     &     ,cloud_fraction_frozen(row_length, rows, wet_model_levels)   &
                                ! Frozen cloud amount (decimal fraction)
     &     ,cclwp(row_length,rows)                                      &
                                   ! condensed water path (kg/m2)
     &     ,t_deep_soil(land_points, st_levels)                         &
                                ! Deep soil temperatures (K)
     &     ,p_star(row_length, rows)                                    &
                                ! Pressure at earth's surface
     &     ,tstar(row_length, rows)                                     &
                                    ! Surface temperature (K)
     &     ,smc(land_points)                                            &
                                  ! Soil moisture content(Kg m^-2)
     &     ,smcl(land_points, sm_levels)                                &
                                ! Soil moisture content in layers
                                ! (Kg m^-2)
     &     ,sthf(land_points, sm_levels)                                &
                                ! Frozen soil moisture content
                                ! of each layer as a fraction of
                                ! saturation.
     &     ,canopy_gb(land_points)                                      &
                                ! Canopy water content (Kg/m2)
     &     ,snodep(row_length, rows)                                    &
                                     ! Snow depth (Kg m^-2)
     &     ,p(row_length, rows, model_levels+1)                         &
                                ! Pressure on rho levels
     &     ,p_theta_levels(row_length, rows, model_levels)              &
                                ! Pressure on theta levels
     &     ,zh(row_length, rows)                                        &
                                ! Height above surface of top
                                ! of boundary layer (m)
     &     ,z0msea(row_length, rows)                                    &
                                     ! Sea surface roughness length
     &     ,sthu(land_points, sm_levels)                                &
                                ! Unfrozen soil moisture content
                                !  of each layer as a fraction of
                                !  saturation.
                                !  (Kg m^-2)
     &     ,gs(row_length*rows)                                         &
                                ! Stomatal conductance
     &     ,LW_incs(row_length, rows, 0:model_levels)                   &
     &     ,surf_radflux(row_length, rows)                              &
     &     ,photosynth_act_rad(row_length, rows)                        &
                                                 ! Net downward
                                ! shortwave radiation in band 1 (w/m2).
     &     ,tstar_tile(row_length*rows,ntype)                           &
                                ! Surface tile temperature
     &     ,aerosol(row_length, rows,wet_model_levels)                  &
                                ! Aerosol values ; only used if
                                ! l_murk=.true. ; default .false.
     &     ,timestep

!     EAK    diagnostic variables for CABLE output
      Real                                                              &
     & FTL_TILE_CAB(row_length*rows,ntype)                              &
     &,FTL_CAB(row_length*rows)                                         &
     &,LE_TILE_CAB(row_length*rows,ntype)                               &
     &,LE_CAB(row_length*rows)                                          &
     &,TSTAR_TILE_CAB(row_length*rows,ntype)                            &
     &,TSTAR_CAB(row_length*rows)                                       &
!     &,SMCL_CAB(row_length*rows,ntype,SM_LEVELS)                        &
!     &,TSOIL_CAB(row_length*rows,ntype,SM_LEVELS)                       &
     &,SMCL_CAB(row_length*rows,SM_LEVELS)                              &
     &,TSOIL_CAB(row_length*rows,SM_LEVELS)                             &
     &,USTAR_CAB(row_length*rows)                                       &
     &,SURF_HTF_CAB(row_length*rows)                                    &
!sxy
!     &,alb_tile(row_length*rows,ntype,4)
     &, SNOW_RHO1L(row_length*rows,ntype)                               &
     &, TOT_ALB(row_length*rows,ntype)                                  &
     &, CH_CAB(row_length*rows)                                         &
     &, CD_CAB(row_length*rows)                                         &
     &, U_S_CAB(row_length*rows)                                        &
     &, CD(row_length*rows)                                             &
     &, CH(row_length*rows)                                             &
     &, SW_DOWN(row_length*rows)                                        &
     &, LW_DOWN(row_length*rows)                                        & 
     &, RADNET_TILE(row_length*rows,ntype)                              &
     &, SNOW_TILE(row_length*rows,ntype)                                &
     &, SNOW_DEPTH3L(row_length*rows,ntype,3)                           &
     &, SNOW_TMP3L(row_length*rows,ntype,3)                             &
     &, SNOW_RHO3L(row_length*rows,ntype,3)                             
!   
!    local arguments    !sxy
     Real                                                               &
    &  snowdepth1(row_length*rows,ntype)                                &  
    &, snowdepth2(row_length*rows,ntype)                                & 
    &, snowdepth3(row_length*rows,ntype)                                & 
    &, snowtmp1(row_length*rows,ntype)                                  &
    &, snowtmp2(row_length*rows,ntype)                                  &
    &, snowtmp3(row_length*rows,ntype)                                  &
    &, snowrho1(row_length*rows,ntype)                                  &
    &, snowrho2(row_length*rows,ntype)                                  &
    &, snowrho3(row_length*rows,ntype)                                  
  
      integer iccb(row_length, rows)                                    &
                                     ! Convective cloud base and top
     &     ,icct(row_length, rows)   ! at levels 1 to model_levels

      Integer nSCMDpkgs              ! No of SCM diagnostics packages
      Logical L_SCMDiags(nSCMDpkgs)  ! Diagnostics packages logicals

! Include parameters necessary for calls to SCMoutput...
#include "s_scmop.h"

      ! Local variables...

      Character*(*), Parameter ::  RoutineName = 'dgnstcs_scm_main'

      real                                                              &
     &      sum_p_col(row_length,rows)                                  &
                                       ! Sums of pressures on theta
                                       ! levels
     &     ,coltemp(row_length,rows)                                    &
                                       ! Average pressure-weighted
                                       ! column temperature
     &     ,colq(row_length,rows)                                       &
                                       ! Average pressure-weighted
                                       ! column humidity
     &     ,accum_pptn(row_length,rows)                                 &
                                       ! Accumulated precipitation
     &     ,pptn_rate(row_length,rows)                                  &
                                       ! Precipitation rate
     &     ,ccb_m(row_length,rows)                                      &
                                       ! Convective cloud base (m)
     &     ,ccb_Pa(row_length,rows)                                     &
                                       ! Convective cloud base (Pa)
     &     ,cct_m(row_length,rows)                                      &
                                       ! Convective cloud top  (m)
     &     ,cct_Pa(row_length,rows)                                     &
                                       ! Convective cloud top  (Pa)
     &     ,rh(row_length,rows,wet_model_levels,2)                      &
                                                   ! Relative humidity
                                       ! (over land and water)
     &     ,qst                                                         &
                                       ! A saturation mixing ratio
! Heights above model surface in metres
     &     ,z_theta(row_length,rows,model_levels)                       &
                                                  ! Theta levels
     &     ,z_rho(row_length,rows,model_levels)                         &
                                                  ! Rho levels
! Actual density (kg/m3)
     &     ,rho_only(row_length, rows, model_levels)

      integer i,j,k

      ! Store diagnostics...
!
!-----------------------------------------------------------------------
!     SCM General Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_gen)) Then

       Do i = 1 , row_length*rows
         Do j=1,  ntype   
            snowdepth1(i,j)=snow_depth3l(i,j,1)
            snowdepth2(i,j)=snow_depth3l(i,j,2)
            snowdepth3(i,j)=snow_depth3l(i,j,3)
            snowtmp1  (i,j)=snow_tmp3l(i,j,1)
            snowtmp2  (i,j)=snow_tmp3l(i,j,2)
            snowtmp3  (i,j)=snow_tmp3l(i,j,3)
            snowrho1  (i,j)=snow_rho3l(i,j,1)
            snowrho2  (i,j)=snow_rho3l(i,j,2)
            snowrho3  (i,j)=snow_rho3l(i,j,3)  
         Enddo
       Enddo

!       Stash 0,002
! DEPENDS ON: scmoutput
        Call SCMoutput(u,                                               &
             'u','Zonal wind','m/s',                                    &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 0,003
! DEPENDS ON: scmoutput
        Call SCMoutput(v,                                               &
             'v','Meridional wind','m/s',                               &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 30,111
! DEPENDS ON: scmoutput
        Call SCMoutput(T,                                               &
             'T','Temperature','K',                                     &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 0,004
! DEPENDS ON: scmoutput
        Call SCMoutput(theta,                                           &
             'theta','Potential temperature','K',                       &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 0,010
! DEPENDS ON: scmoutput
        Call SCMoutput(q,                                               &
             'q','Specific humidity','kg/kg',                           &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 0,150
! DEPENDS ON: scmoutput
        Call SCMoutput(w,                                               &
             'w','Vertical velocity','m/s',                             &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 0,407
! DEPENDS ON: scmoutput
        Call SCMoutput(p,                                               &
             'p_rho','Pressure on rho levels','Pa',                     &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 0,408
! DEPENDS ON: scmoutput
        Call SCMoutput(p_theta_levels,                                  &
             'p_theta','Pressure on theta levels','Pa',                 &
             t_avg,d_all,default_streams,'',RoutineName)

        Do k=1,model_levels
          Do j=1,rows
            Do i=1,row_length
!             Height above model surface as calculation in ni_conv_ctl
              z_theta(i,j,k) = r_theta_levels(i,j,k)                    &
     &                         - r_theta_levels(i,j,0)
              z_rho(i,j,k)   = r_rho_levels(i,j,k)                      &
     &                         - r_theta_levels(i,j,0)
              rho_only(i,j,k) = rho(i,j,k) /                            &
     &                (r_rho_levels(i,j,k) * r_rho_levels(i,j,k))
            End Do ! i
          End Do ! j
        End Do ! k

!       Stash 15,101
! DEPENDS ON: scmoutput
        Call SCMoutput(z_theta,                                         &
             'h_theta','Height of model theta levels','m',              &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 15,102
! DEPENDS ON: scmoutput
        Call SCMoutput(z_rho,                                           &
             'h_rho','Height of model rho levels','m',                  &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 0,253
! DEPENDS ON: scmoutput
        Call SCMoutput(rho,                                             &
             'rho_r2','Density *r*r after timestep','kg/m',             &
             t_avg,d_all,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(rho_only,                                        &
             'rho_only','Density after timestep','kg/m3',               &
             t_avg,d_all,default_streams,'',RoutineName)

        ! Calculate relative humidities on theta levels
        Do j=1,rows
          Do i=1,row_length
            Do k=1,wet_model_levels
! DEPENDS ON: qsat
              Call QSAT(qst,T(i,j,k),p_theta_levels(i,j,k),1)
              rh(i,j,k,1) = q(i,j,k) / qst
! DEPENDS ON: qsat_wat
              Call QSAT_WAT(qst,T(i,j,k),p_theta_levels(i,j,k),1)
              rh(i,j,k,2) = q(i,j,k) / qst
            End Do ! k
          End Do ! i
        End Do ! j

!       Stash 30,113
! DEPENDS ON: scmoutput
        Call SCMoutput(rh(1,1,1,1),                                     &
             'rh','Relative humidity','%',                              &
             t_avg,d_wet,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(rh(1,1,1,2),                                     &
             'rh2','Relative humidity over liquid water','%',           &
             t_avg,d_wet,default_streams,'',RoutineName)

        ! Calculate pressure-weighted average temperature and humidity
        i=1
        j=1
        colq(i,j)=0.0
        coltemp(i,j)=0.0
        sum_p_col(i,j)=0.0
        Do k=1,model_levels
          sum_p_col(i,j) = sum_p_col(i,j) + p_theta_levels(i,j,k)
          coltemp(i,j)   = coltemp(i,j)                                 &
     &                   + T(i,j,k) * p_theta_levels(i,j,k)
          colq(i,j)      = colq(i,j)                                    &
     &                   + Q(i,j,k) * p_theta_levels(i,j,k)
        End Do ! k

! DEPENDS ON: scmoutput
        Call SCMoutput(sum_p_col,                                       &
             'sum_p_col','Sum of theta-level pressures','Pa',           &
             t_acc,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(coltemp,                                         &
             'tatmos','Mean column temperature','K',                    &
             t_acc_div,d_sl,default_streams,'sum_p_col',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(colq,                                            &
             'qatmos','Mean column specific humidity','kg/kg',          &
             t_acc_div,d_sl,default_streams,'sum_p_col',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(aerosol,                                         &
             'aerosol','Aerosol concentration','micro-g/kg',            &
             t_avg,d_wet,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_gen)
!
!-----------------------------------------------------------------------
!     SCM General OR Large Scale Cloud Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_lscld)) Then

!       Stash 0,254
! DEPENDS ON: scmoutput
        Call SCMoutput(qcl,                                             &
             'qcl','QCL cloud water content','kg/kg',                   &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 0,012
! DEPENDS ON: scmoutput
        Call SCMoutput(qcf,                                             &
             'qcf','QCF cloud ice content','kg/kg',                     &
             t_avg,d_wet,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(layer_cloud,                                     &
             'layer_cloud','Layer cloud amount','Fraction',             &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 0,265
! DEPENDS ON: scmoutput
        Call SCMoutput(area_cloud_fraction,                             &
             'acf','Area cloud fraction','Fraction',                    &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 0,266
! DEPENDS ON: scmoutput
        Call SCMoutput(bulk_cloud_fraction,                             &
             'bcf','Bulk cloud fraction','Fraction',                    &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 0,267
! DEPENDS ON: scmoutput
        Call SCMoutput(cloud_fraction_liquid,                           &
             'cfl','Liquid cloud fraction','Fraction',                  &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 0,268
! DEPENDS ON: scmoutput
        Call SCMoutput(cloud_fraction_frozen,                           &
             'cff','Frozen cloud fraction','Fraction',                  &
             t_avg,d_wet,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_lscld)
!
!-----------------------------------------------------------------------
!     SCM General OR Convection Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_conv)) Then

        ! Stash 0,016
        ! Condensed water path convective cloud passed to radiation
! DEPENDS ON: scmoutput
        Call SCMoutput(cclwp,                                           &
             'cclwp2rad','cclwp passed to Radiation','kg/m2',           &
             t_avg,d_sl,default_streams,'',RoutineName)

        ! Stash 0,211
        ! Convective Cloud Amount passed to Radiation
! DEPENDS ON: scmoutput
        Call SCMoutput(cca,                                             &
             'cca2rad','CCA passed to radiation', 'Fraction',           &
             t_avg,d_all,default_streams,'',RoutineName)

        ! Stash 0,212
        ! Convective Cloud Water passed to Radiation
        If (lcv_ccrad .AND. lcv_3d_ccw) Then

!DEPENDS ON: scmoutput
          Call SCMoutput(ccw,                                           &
             'ccw2rad','CCW passed to radiation', 'kg/kg',              &
             t_avg,d_wet,default_streams,'',RoutineName)
        End If

!       Stash 0,258
! DEPENDS ON: scmoutput
        Call SCMoutput(w_adv,                                           &
             'w_adv','Advective vertical velocity','m/s',               &
             t_avg,d_all,default_streams,'',RoutineName)


      End If ! L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_conv)
!
!-----------------------------------------------------------------------
!     SCM General OR Surface Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_surf)) Then

!       Stash 0,409
! DEPENDS ON: scmoutput
        Call SCMoutput(p_star,                                          &
             'pstar','Surface pressure','Pa',                           &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 0,024
! DEPENDS ON: scmoutput
        Call SCMoutput(tstar,                                           &
             'tstar','Surface temperature','K',                         &
             t_avg,d_sl,default_streams,'',RoutineName)

        Call SCMoutput(tstar_cab,                                       &
             'tstar_cab','Surface temperature for CABLE','K',           &
             t_avg,d_sl,default_streams,'',RoutineName)

        Call SCMoutput(ftl_cab,                                         &
             'ftl_cab','Sensible heat flux for CABLE','W/m2',           &
             t_avg,d_sl,default_streams,'',RoutineName)

        Call SCMoutput(le_cab,                                          &
             'le_cab','Latent heat flux for CABLE','W/m2',              &
             t_avg,d_sl,default_streams,'',RoutineName)

        Call SCMoutput(surf_htf_cab,                                    &
             'surf_htf_cab','Ground heat flux for CABLE','W/m2',        &
             t_avg,d_sl,default_streams,'',RoutineName)

        Call SCMoutput(ustar_cab,                                       &
             'ustar_cab','Surface friction velocity for CABLE','m/s',   &
             t_avg,d_sl,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_surf)
!
!-----------------------------------------------------------------------
!     SCM General OR Radiation Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_rad)) Then

!       Stash 2,232
! DEPENDS ON: scmoutput
        Call SCMoutput(lw_incs,                                         &
             'lwrate_day','LW heating rate','K/day',                    &
             t_mult,d_allxtra,default_streams,'ntspday',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(lw_incs,                                         &
             'lwrate_step','LW heating rate','K/timestep',              &
             t_avg,d_allxtra,default_streams,'ntspday',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(lw_incs,                                         &
             'lwrate_acc',                                              &
             'Accumulated LW heating rate over dumping period',         &
             'K/ts*dumps',                                              &
             t_acc,d_allxtra,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_rad)
!
!-----------------------------------------------------------------------
!     SCM Surface OR Radiation Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_surf) .OR. L_SCMDiags(SCMDiag_rad)) Then

! DEPENDS ON: scmoutput
        Call SCMoutput(surf_radflux,                                    &
             'net_surf_rad','Net radiation at surface','W/m2',          &
             t_avg,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(photosynth_act_rad,                              &
             'down_surf_sw_b1','Downward SW radn in band 1','W/m2',     &
             t_avg,d_sl,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_surf) .OR. L_SCMDiags(SCMDiag_rad)
!
!-----------------------------------------------------------------------
!     SCM Convection Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_conv)) Then

      ! Convert ccb as a model level to ccb in meters and Pascals
        Do j=1,rows
          Do i=1,row_length
            If (iccb(i,j) /= 0) Then
              ccb_m(i,j)  = r_theta_levels(i,j,iccb(i,j))               &
     &                    - r_theta_levels(i,j,0)
              cct_m(i,j)  = r_theta_levels(i,j,icct(i,j))               &
     &                    - r_theta_levels(i,j,0)
              ccb_Pa(i,j) = p_theta_levels(i,j,iccb(i,j))
              cct_Pa(i,j) = p_theta_levels(i,j,icct(i,j))
            Else
              ccb_m(i,j)  = 0.0
              cct_m(i,j)  = 0.0
              ccb_Pa(i,j) = 0.0
              cct_Pa(i,j) = 0.0
            End If
          End Do ! i
        End Do ! j

! DEPENDS ON: scmoutput
        Call SCMoutput(ccb_m*cca(1:row_length,1:rows,1),                &
             'ccb','Convective cloud base height (weighted average)',   &
             'm',                                                       &
             t_div,d_sl,default_streams,'cca2rad',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(ccb_Pa*cca(1:row_length,1:rows,1),               &
             'ccbpa',                                                   &
             'Convective cloud base pressure (weighted average)',       &
             'Pa',                                                      &
             t_div,d_sl,default_streams,'cca2rad',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(cct_m*cca(1:row_length,1:rows,1),                &
             'cct','Convective cloud top height (weighted average)',    &
             'm',                                                       &
             t_div,d_sl,default_streams,'cca2rad',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(cct_Pa*cca(1:row_length,1:rows,1),               &
             'cctpa',                                                   &
             'Convective cloud top pressure (weighted average)',        &
             'Pa',                                                      &
             t_div,d_sl,default_streams,'cca2rad',RoutineName)

      End If ! L_SCMDiags(SCMDiag_conv)
!
!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_bl)) Then

! DEPENDS ON: scmoutput
        Call SCMoutput(zh,                                              &
             'bl_height','Boundary layer height','m',                   &
             t_avg,d_sl,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_bl)
!
!-----------------------------------------------------------------------
!     SCM Sea Points Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_sea)) Then

! DEPENDS ON: scmoutput
        Call SCMoutput(z0msea,                                          &
             'sea_roughness','Sea surface roughness length','m',        &
             t_avg,d_sl,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_sea)
!
!-----------------------------------------------------------------------
!     SCM Land Points Diagnostics Package
!-----------------------------------------------------------------------

      If (L_SCMDiags(SCMDiag_land) .AND. land_points > 0) Then

!       Stash 0,020
! DEPENDS ON: scmoutput
        Call SCMoutput(t_deep_soil,                                     &
             'tsoildeep','Deep soil temperatures','K',                  &
             t_avg,d_soilt,default_streams,'',RoutineName)

        Call SCMoutput(tsoil_cab,                                       &
             'tsoil_cab','Deep soil temperatures for CABLE','K',        &
             t_avg,d_soilt,default_streams,'',RoutineName)

!       Stash 8,208
! DEPENDS ON: scmoutput
        Call SCMoutput(smc,                                             &
             'soilmoist','Soil moisture content','kg/m2',               &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 8,223
! DEPENDS ON: scmoutput
        Call SCMoutput(smcl,                                            &
             'soilmoistlay','Layer soil moisture content','kg/m2',      &
             t_avg,d_soilm,default_streams,'',RoutineName)

        Call SCMoutput(smcl_cab,                                        &
            'soilmoislay_cab','Layer soil moisture for CABLE','kg/m2',  &
             t_avg,d_soilm,default_streams,'',RoutineName)

!       Stash 0,214
! DEPENDS ON: scmoutput
        Call SCMoutput(sthu,                                            &
             'soilmoistunfroz','Unfrozen soil moisture content',        &
             'kg/m2',                                                   &
             t_avg,d_soilm,default_streams,'',RoutineName)

!       Stash 0,215
! DEPENDS ON: scmoutput
        Call SCMoutput(sthf,                                            &
             'soilmoistfroz','Frozen soil moisture content','kg/m2',    &
             t_avg,d_soilm,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(gs,                                              &
             'stoma_cond','Stomatal conductance','m/s',                 &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 0,022
! DEPENDS ON: scmoutput
        Call SCMoutput(canopy_gb,                                       &
             'canopy_gb','Canopy water content','kg/m2',                &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 0,023
! DEPENDS ON: scmoutput
        Call SCMoutput(snodep,                                          &
             'snowdepth','Snow depth','kg/m2',                          &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 3,316
! DEPENDS ON: scmoutput
        Call SCMoutput(tstar_tile,                                      &
             'tstartl','Tile surface temp','K',                         &
             t_avg,d_tile,default_streams,'',RoutineName)

! EAK
        Call SCMoutput(tstar_tile_cab,                                  &
             'tstar_tl_cab','Tile surface temp for CABLE','K',          &
             t_avg,d_tile,default_streams,'',RoutineName)

        Call SCMoutput(ftl_tile_cab,                                    &
             'ftl_tile_cab','Sensible heat flux for CABLE','W/m2',      &
             t_avg,d_tile,default_streams,'',RoutineName)

        Call SCMoutput(le_tile_cab,                                     &
             'le_tile_cab','Latent heat flux for CABLE','W/m2',         &
             t_avg,d_tile,default_streams,'',RoutineName)


!        Call SCMoutput(alb_tile,                                        &
!     &        'alb_tile','albedo for tiles','-',                        &
!     &        t_avg,d_alb,default_streams,'',RoutineName)                           

        Call SCMoutput(snow_rho1l,                                       &
             'snow_rho1l','Mean snow density','kg/m3',                   &
              t_avg,d_tile,default_streams,'',RoutineName)                           

        Call SCMoutput(tot_alb,                                          &
              'tot_alb','total albedo','-',                              &
              t_avg,d_tile,default_streams,'',RoutineName)   

        Call  SCMoutput(ch_cab,                                          &
              'ch_cab','Turbulent surface exchange_heat for cable','-',  &
               t_avg,d_sl,default_streams,'',RoutineName)   

        Call  SCMoutput(ch,                                              &
              'ch','Turbulent surface exchange_heat ','-',               &
               t_avg,d_sl,default_streams,'',RoutineName) 

        Call  SCMoutput(cd_cab,                                          &
              'cd_cab','Turbulent surface exchange_mom. for cable','-',  &
               t_avg,d_sl,default_streams,'',RoutineName) 

        Call  SCMoutput(cd,                                              &
              'cd','Turbulent surface exchange_heat ','-',               &
               t_avg,d_sl,default_streams,'',RoutineName)                  

        Call  SCMoutput(u_s_cab,                                         &
              'u_s_cab','surface friction velocity for CABLE','m/s',     &
               t_avg,d_sl,default_streams,'',RoutineName)

        Call  SCMoutput(sw_down,                                         &
              'sw_down','Surface downward SW radiation for CABLE','W/m2',          &
               t_avg,d_sl,default_streams,'',RoutineName)

        Call  SCMoutput(lw_down,                                         &
              'lw_down','Surface downward LW radiation E','W/m2',        &
               t_avg,d_sl,default_streams,'',RoutineName)                             

        Call SCMoutput(snow_tile,                                        &
             'snow_tile','Lying snow on tiles ','Kg/m2',                 &
              t_avg,d_tile,default_streams,'',RoutineName)

        Call SCMoutput(radnet_tile,                                      &
             'radnet_tile','Surface net radiation on land tiles ','W/m2',&
              t_avg,d_tile,default_streams,'',RoutineName)
        
        Call SCMoutput(snowrho1,                                         &
             'snowrho1','the first layer snow density','kg/m3',          &
              t_avg,d_tile,default_streams,'',RoutineName)    

        Call SCMoutput(snowrho2,                                         &
             'snowrho2','the second layer snow density','kg/m3',         &
              t_avg,d_tile,default_streams,'',RoutineName)

        Call SCMoutput(snowrho3,                                         &
             'snowrho3','the third layer snow density','kg/m3',          &
              t_avg,d_tile,default_streams,'',RoutineName)  

        Call SCMoutput(snowtmp1,                                         &
             'snowtmp1','the first layer snow temperature','K',          &
              t_avg,d_tile,default_streams,'',RoutineName)

        Call SCMoutput(snowtmp2,                                         &
             'snowtmp2','the second layer snow temperature','K',         &
              t_avg,d_tile,default_streams,'',RoutineName)

        Call SCMoutput(snowtmp3,                                         &
             'snowtmp3','the third layer snow temperature','K',          &
              t_avg,d_tile,default_streams,'',RoutineName)

        Call SCMoutput(snowdepth1,                                       &
             'snowdepth1','the first layer snow depth','m',              &
              t_avg,d_tile,default_streams,'',RoutineName)

        Call SCMoutput(snowdepth2,                                       &
             'snowdepth2','the second layer snow depth','m',             &
              t_avg,d_tile,default_streams,'',RoutineName)

        Call SCMoutput(snowdepth3,                                       &
             'snowdepth3','the third layer snow depth','m',              &
              t_avg,d_tile,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_land)/land_points 
!
      return
      END SUBROUTINE dgnstcs_scm_main
#endif
