#if defined(A06_0A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calls dummy version 0A of gravity wave drag scheme.
!
      SUBROUTINE g_wave(                                                &
        theta,u,v,row_length,rows,nrows,off_x,off_y,halo_i,halo_j,      &
        model_domain,at_extremity,                                      &
        levels,rho,r_rho_levels,r_theta_levels,sd_orog_land,            &
        orog_grad_xx_land,orog_grad_xy_land,orog_grad_yy_land,          &
        land_index,land_points,timestep,kay,frc,r_u,r_v,                &
        l_taus_scale,l_fix_gwsatn,l_gwd_40km,                           &
        sat_scheme,fsat,                                                &
! diagnostics
        stress_ud     ,stress_ud_on     , points_stress_ud    ,         &
        stress_vd     ,stress_vd_on     , points_stress_vd    ,         &
        stress_ud_satn,stress_ud_satn_on,points_stress_ud_satn,         &
        stress_vd_satn,stress_vd_satn_on,points_stress_vd_satn,         &
        stress_ud_wake,stress_ud_wake_on,points_stress_ud_wake,         &
        stress_vd_wake,stress_vd_wake_on,points_stress_vd_wake,         &
        du_dt_satn    ,du_dt_satn_on    ,points_du_dt_satn    ,         &
        dv_dt_satn    ,dv_dt_satn_on    ,points_dv_dt_satn    ,         &
        du_dt_wake    ,du_dt_wake_on    ,points_du_dt_wake    ,         &
        dv_dt_wake    ,dv_dt_wake_on    ,points_dv_dt_wake    ,         &
        u_s_d         ,u_s_d_on         ,points_u_s_d         ,         &
        v_s_d         ,v_s_d_on         ,points_v_s_d         ,         &
        nsq_s_d       ,nsq_s_d_on       ,points_nsq_s_d       ,         &
        fr_d          ,fr_d_on          ,points_fr_d          ,         &
        bld_d         ,bld_d_on         ,points_bld_d         ,         &
        bldt_d        ,bldt_d_on        ,points_bldt_d        ,         &
        num_lim_d     ,num_lim_d_on     ,points_num_lim_d     ,         &
        num_fac_d     ,num_fac_d_on     ,points_num_fac_d     ,         &
        tausx_d       ,tausx_d_on       ,points_tausx_d       ,         &
        tausy_d       ,tausy_d_on       ,points_tausy_d       ,         &
        taus_scale_d  ,taus_scale_d_on  ,points_taus_scale_d  ,         &
        iret)

      IMPLICIT NONE

!
! SUBROUTINE ARGUMENTS
!
      INTEGER                                                           &
                           !,intent(in):
       row_length          ! number of points per row

      INTEGER                                                           &
                           !,intent(in):
       rows                ! number of rows on theta and u grids

      INTEGER                                                           &
                           !,intent(in):
       nrows               ! number of rows on v grid

      INTEGER                                                           &
                           !,intent(in):
       model_domain        ! MPP switch - global model or not.


      INTEGER                                                           &
                           !,intent(in):
       off_x                                                            &
                           ! small x-halo
      ,off_y                                                            &
                           ! small y-halo
      ,halo_i                                                           &
                           ! large x-halo (for dynamics)
      ,halo_j              ! large y-halo (for dynamics)

      INTEGER                                                           &
                           !,intent(in):
       levels              ! number of model levels

      INTEGER                                                           &
                           !,intent(in):
       land_points         ! number of land points

      INTEGER                                                           &
                           !,intent(in):
       land_index(rows*row_length)
!                          ! index for land points

      INTEGER                                                           &
                           !,intent(in):
        iret               ! return code : iret=0 normal exit

!
! The integers below are set to size land_points if the corresponding
! diagnostic is called or to 1 if it is not. These are set in GWD_CTL2
!
      INTEGER                                                           &
                           !,intent(in):
        points_stress_ud                                                &
      , points_stress_vd                                                &
      , points_stress_ud_satn                                           &
      , points_stress_vd_satn                                           &
      , points_stress_ud_wake                                           &
      , points_stress_vd_wake                                           &
      , points_du_dt_satn                                               &
      , points_dv_dt_satn                                               &
      , points_du_dt_wake                                               &
      , points_dv_dt_wake                                               &
      , points_u_s_d                                                    &
      , points_v_s_d                                                    &
      , points_nsq_s_d                                                  &
      , points_fr_d                                                     &
      , points_bld_d                                                    &
      , points_bldt_d                                                   &
      , points_num_lim_d                                                &
      , points_num_fac_d                                                &
      , points_tausx_d                                                  &
      , points_tausy_d                                                  &
      , points_taus_scale_d
!

      REAL                                                              &
                           !,intent(in):
       theta(row_length,rows,levels)
                           ! primary theta field (K)

      REAL                                                              &
                           !,intent(in):
       rho(1-off_x:row_length+off_x,1-off_y:rows+off_y,levels)
!                          ! density*(radius of earth)^2

      REAL                                                              &
                           !,intent(inout):
       u(1-off_x:row_length+off_x,1-off_y:rows+off_y,levels)
!                          ! primary u field (ms**-1)

      REAL                                                              &
                           !,intent(inout):
       v(1-off_x:row_length+off_x,1-off_y:nrows+off_y,levels)
!                          ! primary v field (ms**-1)

      REAL                                                              &
                           !,intent(in):
       r_rho_levels(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j     &
                  ,levels) ! height of rho level above earth's
!                          ! centre        (m)

      REAL                                                              &
                           !,intent(in):
       r_theta_levels(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j   &
               ,0:levels)  ! height of theta level above earth's
!                          ! centre        (m)


      REAL                                                              &
                           !,intent(in):
       sd_orog_land(land_points)
!                          ! standard deviation of orography (m)

      REAL                                                              &
                           !,intent(in):
       orog_grad_xx_land(land_points)
!                          ! dh/dx squared gradient orography

      REAL                                                              &
                           !,intent(in):
       orog_grad_xy_land(land_points)
!                          ! (dh/dx)(dh/dy) gradient orography

      REAL                                                              &
                           !,intent(in):
       orog_grad_yy_land(land_points)
                           ! dh/dy squared gradient orography

      REAL                                                              &
                           !,intent(in):
       timestep            ! timestep (s)

      REAL                                                              &
                           !,intent(in):
       kay                 ! surface stress constant ( m**-1)

      REAL                                                              &
                           !,intent(in):
       frc                 ! critical Froude number

      REAL                                                              &
                            !,intent(in):
       fsat                 ! Froude number used to scale critical
!                           ! wave amplitude for breaking

!
!  Start of full field diagnostic arrays. This space is allocated
! in GWD_CTL2 only if a diagnostic is called.
!
      REAL                                                              &
                           !,intent(inout):
       r_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
              levels)      ! u wind increment diagnostic


      REAL                                                              &
                           !,intent(inout):
       r_v(1-off_x:row_length+off_x, 1-off_y:nrows+off_y,               &
              levels)      ! v wind increment diagnostic

      REAL                                                              &
                           !,intent(out):
       stress_ud(row_length,  rows,0:levels)                            &
                                                  !u   total stress
      ,stress_vd(row_length, nrows,0:levels)      !v   total stress

      REAL                                                              &
                           !,intent(out):
       stress_ud_satn(row_length,  rows,0:levels)                       &
                                                  !u   satn  stress
      ,stress_vd_satn(row_length, nrows,0:levels) !v   satn  stress

      REAL                                                              &
                           !,intent(out):
       stress_ud_wake(row_length,  rows,0:levels)                       &
                                                  !u   wake  stress
      ,stress_vd_wake(row_length, nrows,0:levels) !v   wake  stress

      REAL                                                              &
                           !,intent(out):
       du_dt_satn(row_length,  rows,levels)                             &
                                            !u acceln (saturation)
      ,dv_dt_satn(row_length, nrows,levels) !v acceln (saturation)

      REAL                                                              &
                           !,intent(out):
       du_dt_wake(row_length,  rows,levels)                             &
                                            !u acceln (blocked flow)
      ,dv_dt_wake(row_length, nrows,levels) !v acceln (blocked flow)

      REAL                                                              &
                           !,intent(out):
       u_s_d  (row_length,  rows)                                       &
                                     ! u_s  diag at theta pts
      ,v_s_d  (row_length,  rows)    ! v_s  diag at theta pts

      REAL                                                              &
                           !,intent(out):
       nsq_s_d(row_length,  rows)    ! nsq_s diag at theta pts

      REAL                                                              &
                           !,intent(out):
       fr_d   (row_length,  rows)    ! Fr diag at theta pts

      REAL                                                              &
                           !,intent(out):
       bld_d  (row_length,  rows)    ! blocked layer depth at theta pts

      REAL                                                              &
                           !,intent(out):
       bldt_d (row_length,  rows)    ! % of time blocked layer diagnosed

      REAL                                                              &
                           !,intent(out):
       num_lim_d (row_length,  rows) ! % of time numerical
                                     ! limiter invoked

      REAL                                                              &
                           !,intent(out):
       num_fac_d (row_length,  rows) ! % redn. of flow-blocking stress
                                     ! after numerical limiter invoked

      REAL                                                              &
                           !,intent(out):
       tausx_d(row_length,  rows)                                       &
                                   !x-component of surface stress
      ,tausy_d(row_length, nrows)  !y-component of surface stress

      REAL                                                              &
                           !,intent(out):
       taus_scale_d(row_length,rows) ! Factor surface stress scaled by
!                                    ! if Froude no. dependence is on.

      LOGICAL                                                           &
                           !,intent(in):
        at_extremity(4)    ! indicates if this processor is at north,
                           ! south, east or west of the processor grid

      LOGICAL                                                           &
                           !,intent(in):
        l_taus_scale                                                    &
                           ! if true then surface stress is made to
!                          ! depend on the low level Froude number
      , l_fix_gwsatn                                                    &
!                          ! if true then invoke minor bug fixes in     
!                          ! gwsatn
      , l_gwd_40km         ! if true then don't apply GWD above 40km 

      INTEGER                                                           &
                            !,intent(in):
       sat_scheme           ! Switch to determine whether to use
!                           ! amplitude or stress based saturation test


!
! END of full field diagnostic arrays
! Below are the stash flags for calculating diagnostics:
!
      LOGICAL                                                           &
                           !,intent(in):
       stress_ud_on                                                     &
                           !u      stress
      ,stress_vd_on                                                     &
                           !v      stress
      ,stress_ud_satn_on                                                &
                           !u satn stress
      ,stress_vd_satn_on                                                &
                           !v satn stress
      ,stress_ud_wake_on                                                &
                           !u wake stress
      ,stress_vd_wake_on                                                &
                           !v wake stress
      ,du_dt_satn_on                                                    &
                           !u accel (saturation)
      ,dv_dt_satn_on                                                    &
                           !v accel (saturation)
      ,du_dt_wake_on                                                    &
                           !u accel blocked flow
      ,dv_dt_wake_on                                                    &
                           !v accel blocked flow
      ,u_s_d_on                                                         &
                           !u_s_d   diag switch
      ,v_s_d_on                                                         &
                           !v_s_d   diag switch
      ,nsq_s_d_on                                                       &
                           !nsq_s_d diag switch
      ,fr_d_on                                                          &
                           !fr_d    switch
      ,bld_d_on                                                         &
                           !bld_d   switch
      ,bldt_d_on                                                        &
                           !bldt_d   switch
      ,num_lim_d_on                                                     &
                           !num_lim_d switch
      ,num_fac_d_on                                                     &
                           !num_fac_d switch
      ,tausx_d_on                                                       &
                           !tausx_d switch
      ,tausy_d_on                                                       &
                           !tausy_d switch
      ,taus_scale_d_on     !taus_scale_d switch

      INTEGER                       ::  ICODE
      CHARACTER (Len=80)            ::  CMESSAGE
      CHARACTER (Len=* ), Parameter ::  RoutineName='G_WAVE0A'

      CMESSAGE = 'Routine should not be callable'
      ICODE = 1
! DEPENDS ON: ereport
      CALL EReport(RoutineName,ICODE,CMESSAGE)

      Return
      END SUBROUTINE g_wave

#endif
