
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calls components of version 4A of gravity wave drag scheme.
!
      SUBROUTINE g_wave(                                                &
     &  theta,u,v,row_length,rows,nrows,off_x,off_y,halo_i,halo_j,      &
     &  model_domain,at_extremity,                                      &
     &  levels,rho,r_rho_levels,r_theta_levels,sd_orog_land,            &
     &  orog_grad_xx_land,orog_grad_xy_land,orog_grad_yy_land,          &
     &  land_index,land_points,timestep,kay,frc,r_u,r_v,                &
     &  l_taus_scale,l_fix_gwsatn,l_gwd_40km,                           &
     &  sat_scheme,fsat,                                                &
! diagnostics
     &  stress_ud     ,stress_ud_on     , points_stress_ud    ,         &
     &  stress_vd     ,stress_vd_on     , points_stress_vd    ,         &
     &  stress_ud_satn,stress_ud_satn_on,points_stress_ud_satn,         &
     &  stress_vd_satn,stress_vd_satn_on,points_stress_vd_satn,         &
     &  stress_ud_wake,stress_ud_wake_on,points_stress_ud_wake,         &
     &  stress_vd_wake,stress_vd_wake_on,points_stress_vd_wake,         &
     &  du_dt_satn    ,du_dt_satn_on    ,points_du_dt_satn    ,         &
     &  dv_dt_satn    ,dv_dt_satn_on    ,points_dv_dt_satn    ,         &
     &  du_dt_wake    ,du_dt_wake_on    ,points_du_dt_wake    ,         &
     &  dv_dt_wake    ,dv_dt_wake_on    ,points_dv_dt_wake    ,         &
     &  u_s_d         ,u_s_d_on         ,points_u_s_d         ,         &
     &  v_s_d         ,v_s_d_on         ,points_v_s_d         ,         &
     &  nsq_s_d       ,nsq_s_d_on       ,points_nsq_s_d       ,         &
     &  fr_d          ,fr_d_on          ,points_fr_d          ,         &
     &  bld_d         ,bld_d_on         ,points_bld_d         ,         &
     &  bldt_d        ,bldt_d_on        ,points_bldt_d        ,         &
     &  num_lim_d     ,num_lim_d_on     ,points_num_lim_d     ,         &
     &  num_fac_d     ,num_fac_d_on     ,points_num_fac_d     ,         &
     &  tausx_d       ,tausx_d_on       ,points_tausx_d       ,         &
     &  tausy_d       ,tausy_d_on       ,points_tausy_d       ,         &
     &  taus_scale_d  ,taus_scale_d_on  ,points_taus_scale_d  ,         &
     &  iret)

      IMPLICIT NONE
!
! Description:
! 1) Interpolate winds to theta points
!    gather data for land points only
! 2) Call surface stress routine
! 3) Calculate stress profiles due to different components of the
!    scheme. Calculate associated wind increments.
! 4) Interpolate acceleration to wind points and update winds
! 5) Gather diagnostics from land points and interpolate from p to
!    u,v staggering on 'c' grid
!
! current code owner: S.Webster
!
! history:
! version   date     comment
! -------   ----     -------
!  5.2   15/11/00   original deck.          Stuart Webster
!  5.3   16/10/01   Remove code no longer required because of
!                   simplifications to scheme. Stuart Webster
!  5.3   27/07/01   Permit a LAM model to run with GWD scheme on.
!                                                      S. Cusack
!  5.4   28/08/02   Introduce numerical limiter for flow blocking
!                   scheme.                            S. Webster
!  6.2   17/05/06   Remove t_to_p reference. P.Selwood.
!  6.2   21/02/06   Introduce surface stress diagnostics. S.Webster
!  6.2   21/02/06   Pass l_taus_scale switch thru' to GWSURF4A
!                   and l_fix_gwsatn thru to GWSATN4A. S. Webster
!
!
!   language: fortran 77 + common extensions.
!   this code is written to umdp3 v6 programming standards.
!
! suitable for single column use, with calls to: uv_to_p removed
!                                                p_to_uv removed
! suitable for rotated grids
!
! global variables (*called comdecks etc...):
!*L------------------COMDECK C_A----------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Replace variable A by more meaningful name for
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Convert to Fixed/Free format. P. Selwood

      ! Mean radius of Earth in m.
      Real, Parameter  :: Earth_Radius = 6371000.

!*----------------------------------------------------------------------
!
! SUBROUTINE ARGUMENTS
!
      INTEGER                                                           &
                           !,intent(in):
     & row_length          ! number of points per row

      INTEGER                                                           &
                           !,intent(in):
     & rows                ! number of rows on theta and u grids

      INTEGER                                                           &
                           !,intent(in):
     & nrows               ! number of rows on v grid

      INTEGER                                                           &
                           !,intent(in):
     & model_domain        ! mpp switch - global model or not.


      INTEGER                                                           &
                           !,intent(in):
     & off_x                                                            &
                           ! small x-halo
     &,off_y                                                            &
                           ! small y-halo
     &,halo_i                                                           &
                           ! large x-halo (for dynamics)
     &,halo_j              ! large y-halo (for dynamics)

      INTEGER                                                           &
                           !,intent(in):
     & levels              ! number of model levels

      INTEGER                                                           &
                           !,intent(in):
     & land_points         ! number of land points

      INTEGER                                                           &
                           !,intent(in):
     & land_index(rows*row_length)
!                          ! index for land points

      INTEGER                                                           &
                           !,intent(in):
     &  iret               ! return code : iret=0 normal exit

!
! The integers below are set to size land_points if the corresponding
! diagnostic is called or to 1 if it is not. These are set in GWD_CTL2
!
      INTEGER                                                           &
                           !,intent(in):
     &  points_stress_ud                                                &
     &, points_stress_vd                                                &
     &, points_stress_ud_satn                                           &
     &, points_stress_vd_satn                                           &
     &, points_stress_ud_wake                                           &
     &, points_stress_vd_wake                                           &
     &, points_du_dt_satn                                               &
     &, points_dv_dt_satn                                               &
     &, points_du_dt_wake                                               &
     &, points_dv_dt_wake                                               &
     &, points_u_s_d                                                    &
     &, points_v_s_d                                                    &
     &, points_nsq_s_d                                                  &
     &, points_fr_d                                                     &
     &, points_bld_d                                                    &
     &, points_bldt_d                                                   &
     &, points_num_lim_d                                                &
     &, points_num_fac_d                                                &
     &, points_tausx_d                                                  &
     &, points_tausy_d                                                  &
     &, points_taus_scale_d
!

      REAL                                                              &
                           !,intent(in):
     & theta(row_length,rows,levels)
                           ! primary theta field (K)

      REAL                                                              &
                           !,intent(in):
     & rho(1-off_x:row_length+off_x,1-off_y:rows+off_y,levels)
!                          ! density*(radius of earth)^2

      REAL                                                              &
                           !,intent(inout):
     & u(1-off_x:row_length+off_x,1-off_y:rows+off_y,levels)
!                          ! primary u field (ms**-1)

      REAL                                                              &
                           !,intent(inout):
     & v(1-off_x:row_length+off_x,1-off_y:nrows+off_y,levels)
!                          ! primary v field (ms**-1)

      REAL                                                              &
                           !,intent(in):
     & r_rho_levels(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j     &
     &            ,levels) ! height of rho level above earth's
!                          ! centre        (m)

      REAL                                                              &
                           !,intent(in):
     & r_theta_levels(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j   &
     &         ,0:levels)  ! height of theta level above earth's
!                          ! centre        (m)


      REAL                                                              &
                           !,intent(in):
     & sd_orog_land(land_points)
!                          ! standard deviation of orography (m)

      REAL                                                              &
                           !,intent(in):
     & orog_grad_xx_land(land_points)
!                          ! dh/dx squared gradient orography

      REAL                                                              &
                           !,intent(in):
     & orog_grad_xy_land(land_points)
!                          ! (dh/dx)(dh/dy) gradient orography

      REAL                                                              &
                           !,intent(in):
     & orog_grad_yy_land(land_points)
                           ! dh/dy squared gradient orography

      REAL                                                              &
                           !,intent(in):
     & timestep            ! timestep (s)

      REAL                                                              &
                           !,intent(in):
     & kay                 ! surface stress constant ( m**-1)

      REAL                                                              &
                           !,intent(in):
     & frc                 ! critical Froude number

      REAL                                                              &
                            !,intent(in):
     & fsat                 ! Froude number used to scale critical
!                           ! wave amplitude for breaking

!
!  Start of full field diagnostic arrays. This space is allocated
! in GWD_CTL2 only if a diagnostic is called.
!
      REAL                                                              &
                           !,intent(inout):
     & r_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &        levels)      ! u wind increment diagnostic


      REAL                                                              &
                           !,intent(inout):
     & r_v(1-off_x:row_length+off_x, 1-off_y:nrows+off_y,               &
     &        levels)      ! v wind increment diagnostic

      REAL                                                              &
                           !,intent(out):
     & stress_ud(row_length,  rows,0:levels)                            &
                                                  !u   total stress
     &,stress_vd(row_length, nrows,0:levels)      !v   total stress

      REAL                                                              &
                           !,intent(out):
     & stress_ud_satn(row_length,  rows,0:levels)                       &
                                                  !u   satn  stress
     &,stress_vd_satn(row_length, nrows,0:levels) !v   satn  stress

      REAL                                                              &
                           !,intent(out):
     & stress_ud_wake(row_length,  rows,0:levels)                       &
                                                  !u   wake  stress
     &,stress_vd_wake(row_length, nrows,0:levels) !v   wake  stress

      REAL                                                              &
                           !,intent(out):
     & du_dt_satn(row_length,  rows,levels)                             &
                                            !u acceln (saturation)
     &,dv_dt_satn(row_length, nrows,levels) !v acceln (saturation)

      REAL                                                              &
                           !,intent(out):
     & du_dt_wake(row_length,  rows,levels)                             &
                                            !u acceln (blocked flow)
     &,dv_dt_wake(row_length, nrows,levels) !v acceln (blocked flow)

      REAL                                                              &
                           !,intent(out):
     & u_s_d  (row_length,  rows)                                       &
                                     ! u_s  diag at theta pts
     &,v_s_d  (row_length,  rows)    ! v_s  diag at theta pts

      REAL                                                              &
                           !,intent(out):
     & nsq_s_d(row_length,  rows)    ! nsq_s diag at theta pts

      REAL                                                              &
                           !,intent(out):
     & fr_d   (row_length,  rows)    ! Fr diag at theta pts

      REAL                                                              &
                           !,intent(out):
     & bld_d  (row_length,  rows)    ! blocked layer depth at theta pts

      REAL                                                              &
                           !,intent(out):
     & bldt_d (row_length,  rows)    ! % of time blocked layer diagnosed

      REAL                                                              &
                           !,intent(out):
     & num_lim_d (row_length,  rows) ! % of time numerical
                                     ! limiter invoked

      REAL                                                              &
                           !,intent(out):
     & num_fac_d (row_length,  rows) ! % redn. of flow-blocking stress
                                     ! after numerical limiter invoked

      REAL                                                              &
                           !,intent(out):
     & tausx_d(row_length,  rows)                                       &
                                   !x-component of surface stress
     &,tausy_d(row_length, nrows)  !y-component of surface stress

      REAL                                                              &
                           !,intent(out):
     & taus_scale_d(row_length,rows) ! Factor surface stress scaled by
!                                    ! if Froude no. dependence is on.

      LOGICAL                                                           &
                           !,intent(in):
     &  at_extremity(4)    ! indicates if this processor is at north,
                           ! south, east or west of the processor grid

      LOGICAL                                                           &
                           !,intent(in):
     &  l_taus_scale                                                    &
                           ! if true then surface stress is made to
!                          ! depend on the low level Froude number
     &, l_fix_gwsatn                                                    &
!                          ! if true then invoke minor bug fixes in     
!                          ! gwsatn
     &, l_gwd_40km         ! if true then don't apply GWD above 40km 

      INTEGER                                                           &
                            !,intent(in):
     & sat_scheme           ! Switch to determine whether to use
!                           ! amplitude or stress based saturation test


!
! END of full field diagnostic arrays
! Below are the stash flags for calculating diagnostics:
!
      LOGICAL                                                           &
                           !,intent(in):
     & stress_ud_on                                                     &
                           !u      stress
     &,stress_vd_on                                                     &
                           !v      stress
     &,stress_ud_satn_on                                                &
                           !u satn stress
     &,stress_vd_satn_on                                                &
                           !v satn stress
     &,stress_ud_wake_on                                                &
                           !u wake stress
     &,stress_vd_wake_on                                                &
                           !v wake stress
     &,du_dt_satn_on                                                    &
                           !u accel (saturation)
     &,dv_dt_satn_on                                                    &
                           !v accel (saturation)
     &,du_dt_wake_on                                                    &
                           !u accel blocked flow
     &,dv_dt_wake_on                                                    &
                           !v accel blocked flow
     &,u_s_d_on                                                         &
                           !u_s_d   diag switch
     &,v_s_d_on                                                         &
                           !v_s_d   diag switch
     &,nsq_s_d_on                                                       &
                           !nsq_s_d diag switch
     &,fr_d_on                                                          &
                           !fr_d    switch
     &,bld_d_on                                                         &
                           !bld_d   switch
     &,bldt_d_on                                                        &
                           !bldt_d   switch
     &,num_lim_d_on                                                     &
                           !num_lim_d switch
     &,num_fac_d_on                                                     &
                           !num_fac_d switch
     &,tausx_d_on                                                       &
                           !tausx_d switch
     &,tausy_d_on                                                       &
                           !tausy_d switch
     &,taus_scale_d_on     !taus_scale_d switch

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
! DOMTYP contains different model domain types
!
! Author : P.Burton
! History:
! Version  Date      Comment.
! 5.0      15/04/99  New comdeck
! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

      INTEGER,PARAMETER:: mt_global        = 1
      INTEGER,PARAMETER:: mt_lam           = 2
      INTEGER,PARAMETER:: mt_cyclic_lam    = 3
      INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
      INTEGER,PARAMETER:: mt_single_column = 5
! DOMTYP end

!--------------------------------------------------------------------
! LOCAL DYNAMIC ARRAYS:
!--------------------------------------------------------------------

      INTEGER                                                           &
     & k_top(land_points)                                               &
                                     ! model level at mountain tops -
!                                    ! exact definition given in gwsurf
     &,k_top_max                     ! max(k_top)

! parameters for mpp
      INTEGER                                                           &
     &   pnorth,                                                        &
                      ! north processor address in the neighbor array
     &   peast,                                                         &
                      ! east processor address in the neighbor array
     &   psouth,                                                        &
                      ! south processor address in the neighbor array
     &   pwest,                                                         &
                      ! west processor address in the neighbor array
     &   nodomain     ! value in neighbor array if the domain has
                      !  no neighbor in this direction. otherwise
                      !  the value will be the tid of the neighbor
      PARAMETER (                                                       &
     &   pnorth   = 1,                                                  &
     &   peast    = 2,                                                  &
     &   psouth   = 3,                                                  &
     &   pwest    = 4,                                                  &
     &   nodomain = -1)

      INTEGER                                                           &
     & i,j,k,l                 ! loop counters in routine

! Work arrays
      REAL                                                              &
     & work_u(row_length,rows,levels)                                   &
     &,work_v(row_length,rows,levels)                                   &
     &,work_halo(1-off_x:row_length+off_x,1-off_y:rows+off_y,0:levels)  &
     &,work_on_v_grid(row_length,nrows,levels)

      REAL                                                              &
     & up_land(land_points,levels)                                      &
                                     ! interpolated u on theta grid
     &,vp_land(land_points,levels)   ! interpolated V on theta grid

      REAL                                                              &
     & theta_land(land_points,levels)! land theta field on theta levels

      REAL                                                              &
     & r_rho_levels_land(land_points,levels)
!                                    !  land field of heights of
!                                    !  rho levels above z=0

      REAL                                                              &
     & r_theta_levels_land(land_points,levels)!)
!                                    !  land field of heights of
!                                    !  theta levels above z=0

      REAL                                                              &
     & rho_land(land_points,levels)  ! density at land points

      REAL                                                              &
     & s_x_lin_stress(land_points)                                      & 
                                     ! 'surface'  x_lin_stress land pnts
     &,s_y_lin_stress(land_points)   ! 'surface'  y_lin_stress land pnts

      REAL                                                              &
     & s_x_wake_stress(land_points)                                     & 
                                     ! 'surface' x_wake_stress land pts
     &,s_y_wake_stress(land_points)  ! 'surface' y_wake_stress land pts

      REAL                                                              &
     & s_x_orog(land_points)                                            & 
                                     ! 'surface' x_orog on land points
     &,s_y_orog(land_points)         ! 'surface' y_orog on land points

      REAL                                                              &
     & du_dt(land_points,levels)                                        &
                                     ! total GWD du/dt on land/theta
     &,dv_dt(land_points,levels)     ! total GWD dv/dt on land/theta

      REAL                                                              &
     & lift(land_points)             ! depth of blocked layer

      REAL                                                              &
     & fr(land_points)               ! low level froude number

      REAL                                                              &
     & rho_s(land_points)            ! low level density

!
! Land points arrays below are for the total GWD stress and
! its 4 individual components. (x and y components for each)
!
      REAL                                                              &
     & stress_ud_land      ( points_stress_ud      , 0:levels )         &
     &,stress_vd_land      ( points_stress_vd      , 0:levels )         &
     &,stress_ud_satn_land ( points_stress_ud_satn , 0:levels )         &
     &,stress_vd_satn_land ( points_stress_vd_satn , 0:levels )         &
     &,stress_ud_wake_land ( points_stress_ud_wake , 0:levels )         &
     &,stress_vd_wake_land ( points_stress_vd_wake , 0:levels )
!
!  Land point arrays below are for the 4 individual components of
!  the GWD wind increment.  (x and y components for each)
!
      REAL                                                              &
     & du_dt_satn_land ( points_du_dt_satn , levels )                   &
     &,dv_dt_satn_land ( points_dv_dt_satn , levels )                   &
     &,du_dt_wake_land ( points_du_dt_wake , levels )                   &
     &,dv_dt_wake_land ( points_dv_dt_wake , levels )

!
!  Land point arrays below are for the 9 GWD 'surface' diagnostics.
!
      REAL                                                              &
     & u_s_d_land     ( points_u_s_d     )                              &
     &,v_s_d_land     ( points_v_s_d     )                              &
     &,nsq_s_d_land   ( points_nsq_s_d   )                              &
     &,fr_d_land      ( points_fr_d      )                              &
     &,bld_d_land     ( points_bld_d     )                              &
     &,bldt_d_land    ( points_bldt_d    )                              &
     &,num_lim_d_land ( points_num_lim_d )                              &
     &,num_fac_d_land ( points_num_fac_d )                              &
     &,tausx_d_land   ( points_tausx_d   )                              &
     &,tausy_d_land   ( points_tausy_d   )                              &
     &,taus_scale_d_land ( points_taus_scale_d )

      REAL                                                              &
     & recip_a2       ! 1/(radius of earth)^2

      PARAMETER                                                         &
     & (recip_a2=1./(earth_radius*earth_radius) )

      LOGICAL                                                           &
     & l_drag(land_points)           ! whether point has a non-zero
!                                    ! stress or not
!
! FUNCTION AND SUBROUTINE CALLS:
      EXTERNAL gw_surf,gw_vert, polar_vector_wind_n
      EXTERNAL p_to_u,p_to_v,u_to_p,v_to_p


!------------------------------------------------------------------
!l    1.1 interpolate winds to p/theta-grid
!------------------------------------------------------------------

! DEPENDS ON: u_to_p
      CALL u_to_p(u,row_length,rows,levels,                             &
     &            off_x, off_y, model_domain,                           &
     &            at_extremity, work_u)
!
! DEPENDS ON: v_to_p
      CALL v_to_p(v,row_length,rows,nrows,                              &
     &            levels, off_x, off_y, model_domain,                   &
     &            at_extremity, work_v)

! set polar winds to zero
      If(model_domain  ==  mt_global) Then

        If (at_extremity(psouth) ) Then
          Do k=1,levels
            Do i=1,row_length
              work_v(i,1,k) = 0.0
              work_u(i,1,k) = 0.0
            End do
          End do
        End If
        If (at_extremity(pnorth) ) Then
          Do k=1,levels
            Do i=1,row_length
              work_v(i,rows,k) = 0.0
              work_u(i,rows,k) = 0.0
            End do
          End do
        End If

      EndIf

!------------------------------------------------------------------
!l    1.2  gather winds at land points
!------------------------------------------------------------------

      Do k=1,levels
        Do l=1,land_points
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          up_land(l,k) =work_u(i,j,k)
          vp_land(l,k) =work_v(i,j,k)
        End do
      End do


!------------------------------------------------------------------
!l    1.3  gather theta, rho and heights at land points
!------------------------------------------------------------------

      If (land_points  >   0) Then


       Do k=1,levels
        Do l=1,land_points
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          r_rho_levels_land(l,k)  = r_rho_levels(i,j,k) -               &
     &                              r_theta_levels(i,j,0)
          r_theta_levels_land(l,k)= r_theta_levels(i,j,k) -             &
     &                              r_theta_levels(i,j,0)
          rho_land(l,k)           = rho(i,j,k)*recip_a2
          theta_land(l,k)         = theta(i,j,k)
        End do
      End do

!------------------------------------------------------------------
!l    2. calculate anisotropic 'surface' stress,CALL gw_surf
!------------------------------------------------------------------

! DEPENDS ON: gw_surf
      CALL gw_surf(                                                     &
     &            r_theta_levels_land,                                  &
     &            rho_land,theta_land,up_land,vp_land,timestep,         &
     &            sd_orog_land,orog_grad_xx_land,orog_grad_xy_land,     &
     &            orog_grad_yy_land,s_x_lin_stress,s_y_lin_stress,      &
     &            s_x_wake_stress,s_y_wake_stress,                      &
     &            s_x_orog,s_y_orog,levels,land_points,kay,rho_s,       &
     &            l_taus_scale, k_top,k_top_max,lift,l_drag,fr,frc,     &
     &            u_s_d_land     , u_s_d_on     , points_u_s_d    ,     &
     &            v_s_d_land     , v_s_d_on     , points_v_s_d    ,     &
     &            nsq_s_d_land   , nsq_s_d_on   , points_nsq_s_d  ,     &
     &            fr_d_land      , fr_d_on      , points_fr_d     ,     &
     &            bld_d_land     , bld_d_on     , points_bld_d    ,     &
     &            bldt_d_land    , bldt_d_on    , points_bldt_d   ,     &
     &            num_lim_d_land , num_lim_d_on , points_num_lim_d,     &
     &            num_fac_d_land , num_fac_d_on , points_num_fac_d,     &
     &            tausx_d_land   , tausx_d_on   , points_tausx_d  ,     &
     &            tausy_d_land   , tausy_d_on   , points_tausy_d  ,     &
     &            taus_scale_d_land    , taus_scale_d_on          ,     &
     &                                    points_taus_scale_d      )


!------------------------------------------------------------------
!l    3. calculate stress profile and accelerations,
!l       CALL gw_vert
!------------------------------------------------------------------

! DEPENDS ON: gw_vert
      CALL gw_vert(                                                     &
     &   rho_land,r_rho_levels_land,r_theta_levels_land,                &
     &   theta_land,up_land,vp_land,levels,land_points,                 &
     &   kay,sd_orog_land,s_x_lin_stress,s_y_lin_stress,                &
     &   s_x_wake_stress,s_y_wake_stress,                               &
     &   s_x_orog,s_y_orog,du_dt,dv_dt,                                 &
     &   k_top,k_top_max,lift,l_drag,fr,rho_s,l_fix_gwsatn,l_gwd_40km,  &
     &   sat_scheme,fsat,                                               &
     &   stress_ud_land ,points_stress_ud ,stress_ud_on,                &
     &   stress_vd_land ,points_stress_vd ,stress_vd_on,                &
     &   stress_ud_satn_land,points_stress_ud_satn,stress_ud_satn_on,   &
     &   stress_vd_satn_land,points_stress_vd_satn,stress_vd_satn_on,   &
     &   stress_ud_wake_land,points_stress_ud_wake,stress_ud_wake_on,   &
     &   stress_vd_wake_land,points_stress_vd_wake,stress_vd_wake_on,   &
     &   du_dt_satn_land,points_du_dt_satn,du_dt_satn_on,               &
     &   dv_dt_satn_land,points_dv_dt_satn,dv_dt_satn_on,               &
     &   du_dt_wake_land,points_du_dt_wake,du_dt_wake_on,               &
     &   dv_dt_wake_land,points_dv_dt_wake,dv_dt_wake_on )

      End If ! on land_points > 0

!------------------------------------------------------------------
!l    4. scatter accelerations to full area, interpolate to uv-grid
!l       and update winds
!------------------------------------------------------------------

! initialise work array: required to ensure zero values generated by
!   non-land points. [note that it is not necessary to re-initialise
!   the work array for each subsequent diagnostic since only land mask
!   points will differ, and these will be overwritten.]

      Do k=0,levels
       Do j=1,rows
        Do i=1,row_length
          work_halo(i,j,k)=0.
        End do !i
       End do  !j
      End do   !k


! expand from land points and interpolate to 'c' u,v grid
      Do k=1,levels

        Do l=1,land_points
           j=(land_index(l)-1)/row_length + 1
           i=land_index(l) - (j-1)*row_length

           work_halo(i,j,k) = du_dt(l,k)

        End do ! l

      End do ! k

! DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                 &
     &     work_halo(1-off_x,1-off_y,1), row_length, rows,              &
     &     levels, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
      Call FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,1),            &
     &                         row_length,rows,levels,off_x,off_y)

! DEPENDS ON: p_to_u
      CALL p_to_u(                                                      &
     &     work_halo(1-off_x,1-off_y,1), row_length, rows,              &
     &     levels, off_x, off_y,                                        &
     &                  work_u)

      Do k=1,levels

        Do l=1,land_points
           j=(land_index(l)-1)/row_length + 1
           i=land_index(l) - (j-1)*row_length

           work_halo(i,j,k) = dv_dt(l,k)

        End do ! l

      End do ! k

! DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                 &
     &     work_halo(1-off_x,1-off_y,1), row_length, rows,              &
     &     levels, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
      Call FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,1),            &
     &                         row_length,rows,levels,off_x,off_y)

! DEPENDS ON: p_to_v
      CALL p_to_v(                                                      &
     &     work_halo(1-off_x,1-off_y,1), row_length, rows,              &
     &     nrows,                                                       &
     &     levels, off_x, off_y,                                        &
     &     work_on_v_grid )

      Do k=1,levels
        Do j=1,rows
          Do i=1,row_length
            r_u(i,j,k)=r_u(i,j,k)+ timestep*work_u(i,j,k)
          End do
        End do
      End do

      Do k=1,levels
        Do j=1,nrows
          Do i=1,row_length
            r_v(i,j,k)=r_v(i,j,k)+ timestep*work_on_v_grid(i,j,k)
          End do
        End do
      End do

!------------------------------------------------------------------
!l    5. gather diagnostics from land points and interpolate from p to
!l       u,v staggering on 'c' grid.
!l       note that stress_ud_on,stress_vd_on are on theta levels
!l       whereas remaining diagnostics are on rho levels.
!------------------------------------------------------------------


      If (stress_ud_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=0,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = stress_ud_land(l,k)

          End do ! l

        End do ! k

! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
     &     work_halo(1-off_x,1-off_y,0), row_length, rows,              &
     &     levels+1, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        Call FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,0),          &
     &                         row_length,rows,levels+1,off_x,off_y)

! DEPENDS ON: p_to_u
        CALL p_to_u(                                                    &
     &     work_halo(1-off_x,1-off_y,0), row_length, rows,              &
     &     levels+1, off_x, off_y,                                      &
     &     stress_ud(1,1,0) )

      EndIf ! on stashflag

      If (stress_vd_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=0,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = stress_vd_land(l,k)

          End do ! l

        End do ! k

! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
     &     work_halo(1-off_x,1-off_y,0), row_length, rows,              &
     &     levels+1, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        Call FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,0),          &
     &                         row_length,rows,levels+1,off_x,off_y)

! DEPENDS ON: p_to_v
        CALL p_to_v(                                                    &
     &     work_halo(1-off_x,1-off_y,0), row_length, rows,              &
     &     nrows,                                                       &
     &     levels+1, off_x, off_y,                                      &
     &     stress_vd(1,1,0) )

      EndIf ! on stashflag


      If (stress_ud_satn_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=0,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = stress_ud_satn_land(l,k)

          End do ! l

        End do ! k

! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
     &     work_halo(1-off_x,1-off_y,0), row_length, rows,              &
     &     levels+1, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        Call FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,0),          &
     &                         row_length,rows,levels+1,off_x,off_y)

! DEPENDS ON: p_to_u
        CALL p_to_u(                                                    &
     &     work_halo(1-off_x,1-off_y,0), row_length, rows,              &
     &     levels+1, off_x, off_y,                                      &
     &     stress_ud_satn(1,1,0) )

      EndIf ! on stashflag

      If (stress_vd_satn_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=0,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = stress_vd_satn_land(l,k)

          End do ! l

        End do ! k

! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
     &     work_halo(1-off_x,1-off_y,0), row_length, rows,              &
     &     levels+1, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        Call FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,0),          &
     &                         row_length,rows,levels+1,off_x,off_y)

! DEPENDS ON: p_to_v
        CALL p_to_v(                                                    &
     &     work_halo(1-off_x,1-off_y,0), row_length, rows,              &
     &     nrows,                                                       &
     &     levels+1, off_x, off_y,                                      &
     &     stress_vd_satn(1,1,0) )

      EndIf ! on stashflag


      If (stress_ud_wake_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=0,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = stress_ud_wake_land(l,k)

          End do ! l

        End do ! k

! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
     &     work_halo(1-off_x,1-off_y,0), row_length, rows,              &
     &     levels+1, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        Call FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,0),          &
     &                         row_length,rows,levels+1,off_x,off_y)

! DEPENDS ON: p_to_u
        CALL p_to_u(                                                    &
     &     work_halo(1-off_x,1-off_y,0), row_length, rows,              &
     &     levels+1, off_x, off_y,                                      &
     &     stress_ud_wake(1,1,0) )

      EndIf ! on stashflag

      If (stress_vd_wake_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=0,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = stress_vd_wake_land(l,k)

          End do ! l

        End do ! k

! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
     &     work_halo(1-off_x,1-off_y,0), row_length, rows,              &
     &     levels+1, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        Call FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,0),          &
     &                         row_length,rows,levels+1,off_x,off_y)

! DEPENDS ON: p_to_v
        CALL p_to_v(                                                    &
     &     work_halo(1-off_x,1-off_y,0), row_length, rows,              &
     &     nrows,                                                       &
     &     levels+1, off_x, off_y,                                      &
     &     stress_vd_wake(1,1,0) )

      EndIf ! on stashflag


      If (du_dt_satn_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=1,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = du_dt_satn_land(l,k)

          End do ! l

        End do ! k

! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
     &     work_halo(1-off_x,1-off_y,1), row_length, rows,              &
     &     levels, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        Call FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,1),          &
     &                         row_length,rows,levels,off_x,off_y)

! DEPENDS ON: p_to_u
        CALL p_to_u(                                                    &
     &     work_halo(1-off_x,1-off_y,1), row_length, rows,              &
     &     levels, off_x, off_y,                                        &
     &     du_dt_satn )

      EndIf ! on stashflag

      If (dv_dt_satn_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=1,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = dv_dt_satn_land(l,k)

          End do ! l

        End do ! k

! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
     &     work_halo(1-off_x,1-off_y,1), row_length, rows,              &
     &     levels, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        Call FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,1),          &
     &                         row_length,rows,levels,off_x,off_y)

! DEPENDS ON: p_to_v
        CALL p_to_v(                                                    &
     &     work_halo(1-off_x,1-off_y,1), row_length, rows,              &
     &     nrows,                                                       &
     &     levels, off_x, off_y,                                        &
     &     dv_dt_satn )

      EndIf ! on stashflag


      If (du_dt_wake_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=1,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = du_dt_wake_land(l,k)

          End do ! l

        End do ! k

! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
     &     work_halo(1-off_x,1-off_y,1), row_length, rows,              &
     &     levels, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        Call FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,1),          &
     &                         row_length,rows,levels,off_x,off_y)

! DEPENDS ON: p_to_u
        CALL p_to_u(                                                    &
     &     work_halo(1-off_x,1-off_y,1), row_length, rows,              &
     &     levels, off_x, off_y,                                        &
     &     du_dt_wake )

      EndIf ! on stashflag


      If (dv_dt_wake_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=1,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = dv_dt_wake_land(l,k)

          End do ! l

        End do ! k

! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
     &     work_halo(1-off_x,1-off_y,1), row_length, rows,              &
     &     levels, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        Call FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,1),          &
     &                         row_length,rows,levels,off_x,off_y)

! DEPENDS ON: p_to_v
        CALL p_to_v(                                                    &
     &     work_halo(1-off_x,1-off_y,1), row_length, rows,              &
     &     nrows,                                                       &
     &     levels, off_x, off_y,                                        &
     &     dv_dt_wake )

      EndIf ! on stashflag


      If (u_s_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             u_s_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            u_s_d(i,j) = u_s_d_land(l)
         End Do ! l

      EndIf ! on stashflag

      If (v_s_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             v_s_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            v_s_d(i,j) = v_s_d_land(l)
         End Do ! l

      EndIf ! on stashflag

      If (nsq_s_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             nsq_s_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            nsq_s_d(i,j) = nsq_s_d_land(l)
         End Do ! l

      EndIf ! on stashflag


      If (fr_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             fr_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            fr_d(i,j) = fr_d_land(l)
         End Do ! l

      EndIf ! on stashflag

      If (bld_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             bld_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            bld_d(i,j) = bld_d_land(l)
         End Do ! l

      EndIf ! on stashflag


      If (bldt_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             bldt_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            bldt_d(i,j) = bldt_d_land(l)
         End Do ! l

      EndIf ! on stashflag


      If (num_lim_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             num_lim_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            num_lim_d(i,j) = num_lim_d_land(l)
         End Do ! l

      EndIf ! on stashflag


      If (num_fac_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             num_fac_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            num_fac_d(i,j) = num_fac_d_land(l)
         End Do ! l

      EndIf ! on stashflag

      If (tausx_d_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid

        Do l=1,land_points
           j=(land_index(l)-1)/row_length + 1
           i=land_index(l) - (j-1)*row_length

           work_halo(i,j,0) = tausx_d_land(l)

        End do ! l

! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
     &     work_halo(1-off_x,1-off_y,0), row_length, rows,              &
     &     1, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        Call FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,0),          &
     &                         row_length,rows,1,off_x,off_y)

! DEPENDS ON: p_to_u
        CALL p_to_u(                                                    &
     &     work_halo(1-off_x,1-off_y,0), row_length, rows,              &
     &     1, off_x, off_y,                                             &
     &     tausx_d(1,1) )

      EndIf ! on stashflag

      If (tausy_d_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,0) = tausy_d_land(l)

          End do ! l

! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
     &     work_halo(1-off_x,1-off_y,0), row_length, rows,              &
     &     1, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        Call FILL_EXTERNAL_HALOS(work_halo(1-off_x,1-off_y,0),          &
     &                         row_length,rows,1,off_x,off_y)

! DEPENDS ON: p_to_v
        CALL p_to_v(                                                    &
     &     work_halo(1-off_x,1-off_y,0), row_length, rows,              &
     &     nrows,                                                       &
     &     1, off_x, off_y,                                             &
     &     tausy_d(1,1) )

      EndIf ! on stashflag

      If (taus_scale_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             taus_scale_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            taus_scale_d(i,j) = taus_scale_d_land(l)
         End Do ! l

      EndIf ! on stashflag


      iret=0

      Return
      END SUBROUTINE g_wave

