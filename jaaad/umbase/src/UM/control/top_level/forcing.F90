#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!-----Subroutine FORCING
!
!     Purpose: Called by SCMMAIN (Single Column Model main routine)
!              to apply the appropriate Forcing (code was previously
!              in the main Calling routine SCMMAIN ).
!
!     Code Description:
!     Language - FORTRAN 77
!
!     Author: C. Bunton
!
!     Documentation: Single Column Model Guide - J. Lean
!=====================================================================
! Options to set initial profiles
!=====================================================================
! (i)   Observational large scale forcing (OBS=TRUE of namelist LOGIC)
!         Initial data is then from namelist INPROF
! (ii)  Statistical large scale forcing (STATS=TRUE of namelist LOGIC)
!         Initial data can either be derived from climate datasets
!         using subroutine INITSTAT or set from namelist
!         INPROF (set ALTDAT=TRUE in namelist LOGIC)
! (iii) No large-scale forcing initial data is set fron namelist
!         INPROF
! (iv)  Continuation from previous run stored on tape
!         (Set TAPEIN=TRUE in namelist LOGIC).  All other initial data
!         is overwritten
!=====================================================================
!---------------------------------------------------------------------
      Subroutine FORCING(                                               &
!     ! IN leading dimensions of arrays
     &  row_length, rows, nlevs, nwet,                                  &
     &  nfor, nbl_levs, nsoilt_levs, nsoilm_levs, ntrop,                &
!     ! IN
     &  sec_day,stats, obs, prindump_obs, prinstat,                     &
     &  dayno_wint, daycount, daysteps, stepcount,                      &
     &  timestep, ichgf,                                                &
     &  ad, at, avn, aw, cdbar, cdsd, ctbar, ctsd,                      &
     &  cvnbar, cvnsd, cwbar, cwsd, dbar, dsd, ddash,                   &
     &  deltan, p, rp,                                                  &
     &  px, py, tbar, tdash, tsd, vnbar, vnsd, vpbar,                   &
     &  wbar, wsd, t, q, u, v, w,                                       &
     &  nSCMDpkgs,L_SCMDiags,relaxT_eurocs,relaxq_eurocs,               &
     &  relaxuv_eurocs,tau_T_eurocs, tau_q_eurocs, tau_uv_eurocs,       &
!
!     ! INOUT
     &  qr, tr, vnr, vpr, wr,                                           &
     &  qcl, qcf,                                                       &
     &  flux_h, flux_e, t_inc, tstar_forcing, tstar, q_inc, u_inc,      &
! Add variables for qcf and qcl
     &  v_inc, w_inc, r_theta_levs, exner_theta_levels,qcl_inc,         &
     &   qcf_inc, ch_flux_h, ch_flux_e,                                 &
     &  ch_tls, ch_qls, ch_uls, ch_vls, ch_wls,                         &

!     ! OUT
     &  dap1, dab1, u_init, v_init, t_init, q_init, ilscnt, rhokh,      &
     &  factor_rhokh, iv, ntab, iy, idum, qls, uls,                     &
     &  vls, wls, tls, l_windrlx, l_vertadv, tau_rlx                    &

     &  )
!
      Implicit none
!
!---------------------------------------------------------------------
!     Arguments
!---------------------------------------------------------------------
!
      Integer                                                           &
     &  row_length, rows                                                &
                                          ! IN leading dimension of SCM
     &  ,nlevs                                                          &
                                ! IN no of levels.
     &  ,nwet                                                           &
                                ! IN no of model levels in which Q is
     &  ,nfor                                                           &
                                ! IN Number terms for observational
                                !    forcing
     &  ,nbl_levs                                                       &
                                ! IN Number of Boundary layer levels
     &  ,nsoilt_levs                                                    &
                                ! IN Number of soil temp levels
     &  ,nsoilm_levs                                                    &
                                ! IN Number of soil moisture levels
     &  ,ntrop                  ! IN Max number of levels in the
                                !     troposphere
      Integer                                                           &
     &  sec_day

      Integer                                                           &
     &  daycount                                                        &
                                ! IN Daynumber (1 represents
                                !  1st january)
     &  ,dayno_wint                                                     &
                                ! IN Daynumber relative to winter
                                !  solstice
     &  ,daysteps                                                       &
                                ! IN No. of timesteps in 1 day
     &  ,ichgf                                                          &
                                ! IN No. of timesteps between change
                                !    in observational forcing
     &  ,ilscnt                                                         &
                                ! INOUT Counts for observational
                                !    forcing
     &  ,ntab                                                           &
                                ! IN Dimension of array used in
                                !    random generator
     &  ,iv(ntab),iy,idum                                               &
                                ! IN state of number generator saved
                                !    on tape for continuation run
     &  ,stepcount              ! IN Timestep counter
      Logical                                                           &
     &  obs                                                             &
     &  ,prindump_obs                                                   &
                                ! T if printout of obs. results
                                !    required
     &  ,prinstat                                                       &
                                ! T if printout of increments
                                !    due to statistical forcing
                                !    required each timestep
     &  ,stats
      Real                                                              &
     &  ad(row_length, rows,nwet-1)                                     &
                                    ! IN term a of eqn 2.22 for dew
     &  ,at(row_length, rows,nlevs-1)                                   &
                                      !    pt depression and temp.
     &  ,avn(row_length, rows,nlevs-1)                                  &
                                      ! IN term a of equation 2.22 for
     &  ,aw(row_length, rows,ntrop-1)                                   &
                                      ! horiz. and vertical velocities
     &  ,cdbar(row_length, rows,nwet)                                   &
                                      ! Mean and SD of random variable
     &  ,cdsd(row_length, rows,nwet)                                    &
                                      !  for dew point depression
     &  ,ctbar(row_length, rows,nlevs)                                  &
                                      ! IN Mean and SD of random
     &  ,ctsd(row_length, rows,nlevs)                                   &
                                      !    variables for temp.
     &  ,cvnbar(row_length, rows,nlevs)                                 &
                                       ! IN Mean and SD of random
     &  ,cvnsd(row_length, rows,nlevs)                                  &
                                       !    variables for velocity VN
     &  ,cwbar(row_length, rows,ntrop)                                  &
                                       ! IN Mean and SD of random
     &  ,cwsd(row_length, rows,ntrop)                                   &
                                      ! variables for vertical velocity
     &  ,dab1(row_length, rows,44)                                      &
                                      ! OUT Observational diagnostics
     &  ,dap1(row_length, rows,36,nlevs)                                &
                                        ! OUT Observational diagnostics
     &  ,dbar(row_length, rows,nwet)                                    &
                                     ! IN Mean and SD dewpoint
     &  ,dsd(row_length, rows,nwet)                                     &
                                     !    depression at daycount days
                                !    from winter solstice (K)
     &  ,ddash(row_length, rows,nwet)                                   &
                                      ! IN Dew pt. corrections
     &  ,deltan(row_length, rows)                                       &
                                    ! IN Radius of area (m)
     &  ,factor_rhokh(row_length, rows)                                 &
     &  ,flux_h(row_length, rows,nfor)                                  &
                                       ! INOUT
     &  ,ch_flux_h(row_length, rows,nfor-1)                             &
                                   ! IN Change per sec in FLUX_H
     &  ,flux_e(row_length, rows,nfor)                                  &
                                          ! INOUT
     &  ,ch_flux_e(row_length, rows,nfor-1)                             &
                                   ! IN Change per sec in FLUX_E
     &  ,p(row_length, rows,nlevs)                                      &
                                   ! IN Pressure coordinates (Pa)
     &  ,r_theta_levs(row_length,rows,0:nlevs)                          &
                                               ! IN
                                   ! Distance of theta levels from
                                   ! centre of the earth (m)
     &  ,exner_theta_levels(row_length,rows,nlevs)                      &

     &  ,px(row_length, rows,ntrop)                                     &
                                   ! IN Reciprocal log functions for
     &  ,py(row_length, rows,ntrop-1)                                   &
                                      !    calc. of vert. advection
                                !    used in eqns 2.12 and 2.13
     &  ,q(row_length, rows,nwet)                                       &
                                   ! INOUT Specific humidity (kg kg^-1)
     &  ,u_init(row_length,rows,nlevs)                                  &
                                   ! INOUT u compt. of wind (m s^-1)
     &  ,v_init(row_length,rows,nlevs)                                  &
                                   ! INOUT v compt. of wind (m s^-1)
     &  ,qcl(row_length, rows, nwet)                                    &
                                     ! Cloud water content(Kg Kg^-1)
     &  ,qcf(row_length, rows, nwet)                                    &
                                     ! Cloud ice content(Kg Kg^-1)
     &  ,q_init(row_length, rows,nwet)                                  &
                                       ! OUT Initial specific humidity
                                !    (kg kg^-1)
     &  ,qr(row_length, rows,nwet,2)                                    &
                                     ! INOUT Randomly sampled humidity
                                ! (kg kg^-1)
     &  ,rhokh(row_length, rows,nbl_levs)                               &
     &  ,rp(row_length, rows,nlevs)                                     &
                                    ! IN Reciprocal pressure for rho
                                !    levels (HPa or mb ^-1)
     &  ,t(row_length, rows,nlevs)                                      &
                                   ! IN Temp (K)
     &  ,t_init(row_length, rows,nlevs)                                 &
                                        ! OUT Initial temp (K)
     &  ,tbar(row_length, rows,nlevs)                                   &
                                       ! IN Mean and SD temperature at
                                !    daycount days from
                                !    winter solstice (K)
     &  ,tdash(row_length, rows,nlevs)                                  &
                                       ! IN Temp. corrections (K)
     &  ,timestep                                                       &
                                ! IN model timestep (s)
     &  ,t_inc(row_length, rows,nfor,nlevs)                             &
                                ! INOUT Temp increment due to
     &  ,tls(row_length, rows,nfor,nlevs)                               &
     &  ,tstar_forcing(row_length, rows,nfor)                           &
                                ! INOUT Sea Surface temperature
                                ! increment due to forcing
     &  ,tstar(row_length, rows)                                        &
                                !    large-scale horizontal and
                                !    vertical advection
                                !    (K s^-1 day^-1)
     &  ,ch_tls(row_length, rows,nfor-1,nlevs)                          &
                                ! IN Change per sec in Temp incs
     &  ,q_inc(row_length, rows,nfor,nwet)                              &
                                ! Specific humidity increment
     &  ,qls(row_length, rows,nfor,nwet)                                &
                                !  due to large-scale horizontal
                                !  and vertical advection
                                !  (kg kg^-1 s^-1 day^-1)
     &  ,ch_qls(row_length, rows,nfor-1,nwet)                           &
                                ! IN change per sec in Spec humid
                                !    increment
     &  ,u_inc(row_length, rows,nfor,nlevs)                             &
                                ! Zonal and meridional wind
     &  ,uls(row_length, rows,nfor,nlevs)                               &
     &  ,v_inc(row_length, rows,nfor,nlevs)                             &
                                !  increment due to large-scale
     &  ,vls(row_length, rows,nfor,nlevs)                               &
                                !  horizontal and vertical
     &  ,w_inc(row_length, rows,nfor,nlevs)                             &
     &  ,wls(row_length, rows,nfor,nlevs)                               &
     &  ,ch_uls(row_length, rows,nfor-1,nlevs)                          &
                                ! IN change per sec in Zonal and
     &  ,ch_vls(row_length, rows,nfor-1,nlevs)                          &
                                               !    merid wind
     &  ,ch_wls(row_length, rows,nfor-1,nlevs)                          &
     &  ,tsd(row_length, rows,nlevs)                                    &
                                     ! IN SD of temp. at daycount days
                                !    from winter solstice (K)
     &  ,tr(row_length, rows,nlevs,2)                                   &
                                     ! INOUT Randomly sampled temp. (K)
     &  ,u(row_length, rows,nlevs)                                      &
                                     ! IN Zonal and meridional winds
     &  ,v(row_length, rows,nlevs)                                      &
                                     !    (m s^-1)
     &  ,w(row_length, rows,0:nlevs)                                    &
                                       !IN vert. vel.  (m s^-1)
     &  ,vnbar(row_length, rows,nlevs)                                  &
                                      ! IN Mean and SD velocity VN at
                                !    daycount days from
                                !    winter solstice (m s^-1)
     &  ,vnr(row_length, rows,nlevs,2)                                  &
                                      ! INOUT Randomly sampled horiz
                                !    velocity (m s^-1)
     &  ,vnsd(row_length, rows,nlevs)                                   &
                                     ! IN Mean and SD velocity VN at
                                !    daycount days from
                                !    winter solstice (m s^-1)
     &  ,vpbar(row_length, rows,nlevs)                                  &
                                      ! IN Mean  velocity VP at
                                !    daycount days from
                                !    winter solstice (m s^-1)
     &  ,vpr(row_length, rows,nlevs,2)                                  &
                                      ! INOUT Randomly sampled horiz
                                !    velocity (m s^-1)
     &  ,wbar(row_length, rows,ntrop)                                   &
                                      ! IN Mean and SD vertical
     &  ,wsd(row_length, rows,ntrop)                                    &
                                      !    velocity at daycount days
                                !    from winter solstice (mb s^-1)
     &  ,wr(row_length, rows,ntrop,2)                                   &
                                      ! INOUT Randomly sampled vertical
                                !    velocity (mb s^-1)

! Additional variables for qcl and qcf forcing
     &  ,qcl_inc(row_length, rows,nfor,nwet)                            &
     &  ,qcf_inc(row_length, rows,nfor,nwet)                            &
     &  ,qclls(row_length, rows,nfor,nwet)                              &
     &  ,qcfls(row_length, rows,nfor,nwet)

! Include parameters necessary for calls to SCMoutput...
#include "s_scmop.h"

! The increments to T,u,v,q,qcl,qcf on input
      Real                                                              &
     &   t_inc_in(row_length,rows,nlevs)                                &
     &  ,u_inc_in(row_length,rows,nlevs)                                &
     &  ,v_inc_in(row_length,rows,nlevs)                                &
     &  ,q_inc_in(row_length,rows,nlevs)                                &
     &  ,qcl_inc_in(row_length,rows,nlevs)                              &
     &  ,qcf_inc_in(row_length,rows,nlevs)                              &
!
! The increments to T,u,v,q due to observational forcing
     &  ,t_inc_obs(row_length,rows,nlevs)                               &
     &  ,u_inc_obs(row_length,rows,nlevs)                               &
     &  ,v_inc_obs(row_length,rows,nlevs)                               &
     &  ,q_inc_obs(row_length,rows,nlevs)

! Additional variables to update winds with stats forcing

      Real                                                              &
     &   rdaysteps                                                      &
     &  ,rstepcount
      Real                                                              &
     &   u_old(row_length, rows,nlevs)                                  &
     &  ,v_old(row_length, rows,nlevs)                                  &
     &  ,w_old(row_length, rows,nlevs)

      Logical                                                           &
     &  L_windrlx                                                       &
     &, L_vertadv

      Real                                                              &
     &  tau_rlx                                                         &
     & , th_p, th_m

      Real :: T_vertadv(row_length,rows,nlevs)
      Real :: q_vertadv(row_length,rows,nlevs)
      Real :: qcl_vertadv(row_length,rows,nlevs)
      Real :: qcf_vertadv(row_length,rows,nlevs)

      Integer                                                           &
     & nSCMDpkgs             ! No of SCM diagnostics packages

      Logical                                                           &
     & L_SCMDiags(nSCMDpkgs)                                             
                             ! Logicals for SCM diagnostics packages

      Logical :: relaxT_eurocs   !switches controlling EUROCS forcing 
      Logical :: relaxq_eurocs                                                
      Logical :: relaxuv_eurocs   

      Real :: tau_T_eurocs          !timescale for relaxation of T
      Real :: tau_q_eurocs          !timescale for relaxation of q
      Real :: tau_uv_eurocs         !timescale for relaxation of u and v

!---------------------------------------------------------------------
!     Local variables
!---------------------------------------------------------------------
      Integer                                                           &
     &  i, j, k, l
      Real                                                              &
     &  tstpfd
      Real                                                              &
     &  alpha

      Real ::                                                           &
        a2out_all (nlevs)                                               &
      , a2out_wet (nwet)

      CHARACTER*(80) CMESSAGE  ! Error message if ErrorStatus > 0
      CHARACTER (Len=*), Parameter :: RoutineName = 'forcing'

      Integer ErrorStatus

!---------------------------------------------------------------------
!     Control variable
!---------------------------------------------------------------------
      tstpfd = timestep / sec_day
      alpha  = timestep / 3600.

!---------------------------------------------------------------------
!     Make copies of the increments as they are on input
!---------------------------------------------------------------------

      t_inc_in=t_inc(1:row_length,1:rows,nfor,1:nlevs)
      u_inc_in  =u_inc(1:row_length,1:rows,nfor,1:nlevs)
      v_inc_in  =v_inc(1:row_length,1:rows,nfor,1:nlevs)
      q_inc_in  =q_inc(1:row_length,1:rows,nfor,1:nlevs)
      qcl_inc_in=qcl_inc(1:row_length,1:rows,nfor,1:nlevs)
      qcf_inc_in=qcf_inc(1:row_length,1:rows,nfor,1:nlevs)

!
!---------------------------------------------------------------------
!     Set instantaneous profiles and budgets to zero for OBS forcing
!---------------------------------------------------------------------
!
      If (obs) then
        If(.not. relaxT_eurocs) Then
          Do i = 1, nlevs
            Do j = 1, rows
              Do l = 1, row_length
                t_init(l,j,i) = t(l,j,i)
              End Do
            End Do
          End Do                   ! i
        End If                     !not relaxT_eurocs
        If(.not. relaxq_eurocs) Then
          Do i = 1, nwet
            Do j = 1, rows
              Do l = 1, row_length
                q_init(l,j,i) = q(l,j,i)
              End Do
            End Do
          End Do                   ! i
        End If                     !not relaxq_eurocs
        If (prindump_obs) Then
          Do i = 1, nlevs
            Do j = 1, 36
              Do k = 1, rows
                Do l = 1, row_length
                  dap1(l,k,j,i) = 0.0
                End Do
              End Do
            End Do               ! j
          End Do                 ! i
          Do j = 1, 44
            Do k = 1, rows
              Do l = 1, row_length
                dab1(l,k,j) = 0.0
              End Do
            End Do
          End Do                 ! j
        End If                   ! prindump_obs
      End If                     ! OBS
!
!---------------------------------------------------------------------
!     If statistical forcing required:-
!     Set up 2 profiles. 1 for start of day plus 1 for start of
!     following day and linearly interpolate between 2 values for
!     all forcing variables.  Increments to T and Q added and U
!     and V calculated
!---------------------------------------------------------------------
!

! Initialise qcl and qcf keep
        Do k = 1, nlevs
          Do j = 1, rows
            Do l = 1, nfor
              Do i = 1, row_length
                qclls(i,j,l,k) = 0.
                qcfls(i,j,l,k) = 0.
              End Do
            End Do
          End Do
        End Do

      If (stats) then
! _ls incs = Atmos_physics1 ; _inc = forced incs
!  We need to keep hold of the increments calculated from physics1 and
!  then add them on to the ones calculated from the stats forcing

        Do k = 1, nlevs
          Do j = 1, rows
            Do i = 1, row_length
              tls(i,j,nfor,k) = t_inc(i,j,nfor,k)
              uls(i,j,nfor,k) = u_inc(i,j,nfor,k)
              vls(i,j,nfor,k) = v_inc(i,j,nfor,k)
              wls(i,j,nfor,k) = w_inc(i,j,nfor,k)
            End Do
          End Do
        End Do
        Do k = 1, nwet
          Do j = 1, rows
            Do i = 1, row_length
              qls(i,j,nfor,k) = q_inc(i,j,nfor,k)
! qcl and qcf added
              qclls(i,j,nfor,k) = qcl_inc(i,j,nfor,k)
              qcfls(i,j,nfor,k) = qcf_inc(i,j,nfor,k)

            End Do
          End Do
        End Do

! DEPENDS ON: statstep
        Call STATSTEP(                                                  &
!       ! IN
     &    row_length, rows, nlevs, nwet, ntrop,                         &
!       !
     &    deltan, px, py, daysteps, stepcount, dayno_wint,              &
     &    tr, vnr, vpr, qr, wr, tbar, tsd, tdash,                       &
     &    dbar, dsd, ddash, vnbar, vpbar,                               &
     &    vnsd, wbar, wsd, ctbar,                                       &
     &    ctsd, at, cdbar, cdsd, ad, cvnbar, cvnsd, avn,                &
     &    cwbar, cwsd, aw, p, rp, u, v, w, t, q, prinstat,              &
     &    u_inc(1,1:rows,nfor:nfor,1:nlevs),                            &
     &    v_inc(1,1:rows,nfor:nfor,1:nlevs),                            &
     &    w_inc(1,1:rows,nfor:nfor,1:nlevs),                            &
     &    t_inc(1,1:rows,nfor:nfor,1:nlevs),                            &
     &    q_inc(1,1:rows,nfor:nfor,1:nwet), daycount, timestep,         &
     &    iv, ntab, iy, idum)

! qcl and qcf temporarily forced using q/100.

        Do k = 1, nlevs
          Do j = 1, rows
            Do i = 1, row_length
              qcl_inc(i,j,nfor,k) = q_inc(i,j,nfor,k)/100.
              qcf_inc(i,j,nfor,k) = q_inc(i,j,nfor,k)/100.
            End Do
          End Do
        End Do


!  Need to add on increments calculated in physics1

        Do k = 1, nlevs
          Do j = 1, rows
            Do i = 1, row_length
              t_inc(i,j,nfor,k) = t_inc(i,j,nfor,k) +                   &
     &                               tls(i,j,nfor,k)
              u_inc(i,j,nfor,k) = u_inc(i,j,nfor,k) +                   &
     &                               uls(i,j,nfor,k)
              v_inc(i,j,nfor,k) = v_inc(i,j,nfor,k) +                   &
     &                               vls(i,j,nfor,k)
            End Do
          End Do
        End Do

        Do k = 1, nwet
          Do j = 1, rows
            Do i = 1, row_length
              q_inc(i,j,nfor,k) = q_inc(i,j,nfor,k) +                   &
     &                               qls(i,j,nfor,k)
! Add qcl and qcf
              qcl_inc(i,j,nfor,k) = qcl_inc(i,j,nfor,k) +               &
     &                               qclls(i,j,nfor,k)
              qcf_inc(i,j,nfor,k) = qcf_inc(i,j,nfor,k) +               &
     &                               qcfls(i,j,nfor,k)

            End Do
          End Do
        End Do

! W not used in 4.5, added here, but will need someone to use
! stats forcing in anger to modify.
        Do k = 1, ntrop
          Do j = 1, rows
            Do i = 1, row_length
              w_inc(i,j,nfor,k) = w_inc(i,j,nfor,k) +                   &
     &                               wls(i,j,nfor,k)
            End Do
          End Do
        End Do


      elseif (obs) then
! _ls = forced increments ; _inc = Atmos_Physics1 incs

!
!---------------------------------------------------------------------
!       Select forcing value for time of day
!---------------------------------------------------------------------
!
! NOTE: ilscnt is used to count the number of observational profiles
!       as the run proceeds. Among other things, it is used to select
!       the appropriate change/gradient to the forcing tendencies when
!       observational intervals greater than the model timestep.
!       THIS FUNCTIONALITY IS FLAWED! It is recommended that users
!       have ilscnt set to 0 in their SCM namelists and have obs
!       forcings commence at the same time they wish to begin their
!       SCM run.                                    Wong/Kerr-Munslow
        If (mod((daycount-1) * int(sec_day)                             &
     &    + (stepcount-1) * int(timestep)                               &
     &    , ichgf * int(timestep))  ==  0)                              &
     &    ilscnt = ilscnt + 1
        If (ilscnt  ==  0) ilscnt = 1
        If (ilscnt  >=  nfor) then
          write(Cmessage,*) 'time exceeds forcing period'
          ErrorStatus = 1
! DEPENDS ON: ereport
          Call Ereport(RoutineName,ErrorStatus,Cmessage)
        End If
!
!---------------------------------------------------------------------
!       t_inc(nfor,k) etc contains present value of forcing
!---------------------------------------------------------------------
!
        Do l = 1, rows
          Do k = 1, row_length
            flux_h(k,l,nfor) = flux_h(k,l,nfor)                         &
     &      +              ch_flux_h(k,l,ilscnt) * timestep
            flux_e(k,l,nfor) = flux_e(k,l,nfor)                         &
     &      +              ch_flux_e(k,l,ilscnt) * timestep
            tstar_forcing(k,l,nfor) = tstar_forcing(k,l,nfor)           &
     &      +              tstar_forcing(k,l,ilscnt) * timestep
            rhokh(k,l,1) = flux_h(k,l,nfor)
            factor_rhokh(k,l) = flux_e(k,l,nfor)
            tstar(k,l) = tstar_forcing(k,l,nfor)
          Do i = 1, nlevs
              tls(k,l,nfor,i) = tls(k,l,nfor,i)                         &
     &          +             ch_tls(k,l,ilscnt,i) * timestep
              uls(k,l,nfor,i) = uls(k,l,nfor,i)                         &
     &          +             ch_uls(k,l,ilscnt,i) * timestep
              vls(k,l,nfor,i) = vls(k,l,nfor,i)                         &
     &          +             ch_vls(k,l,ilscnt,i) * timestep
              wls(k,l,nfor,i) = wls(k,l,nfor,i)                         &
     &          +             ch_wls(k,l,ilscnt,i) * timestep

! Observed value is K/day, so mult by timestep_day_fraction to
! put into Kelvin for addition to T

              t_inc(k,l,nfor,i) = t_inc(k,l,nfor,i) +                   &
     &                    tls(k,l,nfor,i) * tstpfd

                   !  If using forcing which relaxes back to
                   !  initial profiles
                   !  Damp back to initial profiles above 910 hPa
                   !
              If(relaxT_eurocs) Then
                If (p(k,l,i) .lt. 91000.) Then 
                  t_inc(k,l,nfor,i) = t_inc(k,l,nfor,i) -               &
     &               (t(k,l,i)-t_init(k,l,i))*timestep/tau_T_eurocs
                End If
              End If
              If(relaxuv_eurocs) Then
                If (p(k,l,i) .lt. 91000.) Then 
                  u_inc(k,l,nfor,i) = u_inc(k,l,nfor,i) -               &
     &               (u(k,l,i)-u_init(k,l,i))*timestep/tau_uv_eurocs
                  v_inc(k,l,nfor,i) = v_inc(k,l,nfor,i) -               &
     &               (v(k,l,i)-v_init(k,l,i))*timestep/tau_uv_eurocs
                End If 
              End If     !relaxuv_eurocs
!
! t_inc now contains the effects due to the observed forcing only
!
!-----------------------------------------------------------------------
!         SCM Forcing OR Increments Diagnostics Package
!-----------------------------------------------------------------------
          If (L_SCMDiags(SCMDiag_forc)                                  &
     &        .OR. L_SCMDiags(SCMDiag_incs)) Then
!
! Store increment due to observational forcing for output later
!
            t_inc_obs(k,l,i) = t_inc(k,l,nfor,i) - t_inc_in(k,l,i)
!
          End If


          T_vertadv(k,l,i) = 0.0

          If (L_vertadv) then
!
! large scale forcing does not include vertical advection
! approximate vertical advection using upstream approximation
! note that tls holds horizontal advection plus other fixed forcing.
!
            If(wls(k,l,nfor,i) <  0.0) Then
!
! subsidence
!
              If(i <  nlevs) Then
                th_p = t(k,l,i+1)/exner_theta_levels(k,l,i+1)
                th_m = t(k,l,i)  /exner_theta_levels(k,l,i)
                T_vertadv(k,l,i) = - timestep*wls(k,l,nfor,i)*          &
     &                 exner_theta_levels(k,l,i)*( th_p - th_m )/       &
     &                 (r_theta_levs(k,l,i+1)-r_theta_levs(k,l,i))
              Endif

            Else
!
! ascent
!
              If(i >  1) then
                th_p = t(k,l,i)  /exner_theta_levels(k,l,i)
                th_m = t(k,l,i-1)/exner_theta_levels(k,l,i-1)
                T_vertadv(k,l,i) = - timestep*wls(k,l,nfor,i)*          &
     &                 exner_theta_levels(k,l,i)*( th_p - th_m )/       &
     &                 (r_theta_levs(k,l,i)-r_theta_levs(k,l,i-1))
              Endif

            Endif ! wls<0

            t_inc(k,l,nfor,i) = t_inc(k,l,nfor,i) + T_vertadv(k,l,i)

          Endif ! vertadv

          IF (.not.L_windrlx) Then
!
! wind forcing specified as an increment to the winds
!
                u_inc(k,l,nfor,i)=u_inc(k,l,nfor,i) +                   &
     &                       uls(k,l,nfor,i) * tstpfd
                v_inc(k,l,nfor,i)=v_inc(k,l,nfor,i) +                   &
     &                       vls(k,l,nfor,i) * tstpfd
          Else

             If(tau_rlx  >   0.0) Then
!
! relax wind components to observations over given tinescale
!
                u_inc(k,l,nfor,i) = u_inc(k,l,nfor,i) -                 &
     &                        (u(k,l,i)-uls(k,l,nfor,i))*               &
     &                        timestep/tau_rlx
                v_inc(k,l,nfor,i) = v_inc(k,l,nfor,i) -                 &
     &                        (v(k,l,i)-vls(k,l,nfor,i))*               &
     &                        timestep/tau_rlx
             Else
!
! Set winds to values specified by ULS and VLS
!
                u_inc(k,l,nfor,i)=u_inc(k,l,nfor,i)+                    &
     &                            (uls(k,l,nfor,i)-u(k,l,i))
                v_inc(k,l,nfor,i)=v_inc(k,l,nfor,i)+                    &
     &                            (vls(k,l,nfor,i)-v(k,l,i))
             Endif ! tau_rlx > 0

          Endif ! L_windrlx

!
!-----------------------------------------------------------------------
!         SCM Forcing OR Increments Diagnostics Package
!-----------------------------------------------------------------------
          If (L_SCMDiags(SCMDiag_forc)                                  &
     &        .OR. L_SCMDiags(SCMDiag_incs)) Then
!
! Store increment due to observational forcing for output later
!
            u_inc_obs(k,l,i) = u_inc(k,l,nfor,i) - u_inc_in(k,l,i)
            v_inc_obs(k,l,i) = v_inc(k,l,nfor,i) - v_inc_in(k,l,i)
!
          End If

          enddo                 ! i
          Do i = 1, nwet
              qls(k,l,nfor,i) = qls(k,l,nfor,i)                         &
     &          +             ch_qls(k,l,ilscnt,i) * timestep
! Observed value is Kg/kg day^-1, so mult by timestep_day_fraction to
! make unitless for addition to q

              q_inc(k,l,nfor,i) = q_inc(k,l,nfor,i) +                   &
     &                  qls(k,l,nfor,i) * tstpfd
!
!-----------------------------------------------------------------------
!         SCM Forcing OR Increments Diagnostics Package
!-----------------------------------------------------------------------
          If (L_SCMDiags(SCMDiag_forc)                                  &
     &        .OR. L_SCMDiags(SCMDiag_incs)) Then
!
! Store increment due to observational forcing for output later
!
            q_inc_obs(k,l,i) = q_inc(k,l,nfor,i) - q_inc_in(k,l,i)
!
          End If

!             !  If relaxing back to input q profile
!             !  Damp back to initial profiles above 910 hPa
!             !
              If(relaxq_eurocs) then
                If (p(k,l,i) .lt. 91000.) Then 
                  q_inc(k,l,nfor,i) = q_inc(k,l,nfor,i) -               &
     &               (q(k,l,i)-q_init(k,l,i))*timestep/tau_q_eurocs
                End If
              End If           !relaxq_eurocs




          q_vertadv(k,l,i)   = 0.0
          qcl_vertadv(k,l,i) = 0.0
          qcf_vertadv(k,l,i) = 0.0

          If (L_vertadv) Then
!
! large scale forcing does not include vertical advection
! approximate vertical advection using upstream approximation
! note that tls holds horizontal advection plus other fixed forcing.
!
            If(wls(k,l,nfor,i) <  0.0) Then
!
! subsidence
!
              If(i <  nlevs) Then
                q_vertadv(k,l,i)   = - timestep*wls(k,l,nfor,i)*        &
     &                     (q(k,l,i+1)-q(k,l,i))/                       &
     &                     (r_theta_levs(k,l,i+1)-r_theta_levs(k,l,i))
                qcl_vertadv(k,l,i) = - timestep*wls(k,l,nfor,i)*        &
     &                     (qcl(k,l,i+1)-qcl(k,l,i))/                   &
     &                     (r_theta_levs(k,l,i+1)-r_theta_levs(k,l,i))
                qcf_vertadv(k,l,i) = - timestep*wls(k,l,nfor,i)*        &
     &                     (qcf(k,l,i+1)-qcf(k,l,i))/                   &
     &                     (r_theta_levs(k,l,i+1)-r_theta_levs(k,l,i))
              Endif

            Else
!
! ascent
!
              If(i >  1) Then
                q_vertadv(k,l,i) = - timestep*wls(k,l,nfor,i)*          &
     &                     (q(k,l,i)-q(k,l,i-1))/                       &
     &                     (r_theta_levs(k,l,i)-r_theta_levs(k,l,i-1))
                qcl_vertadv(k,l,i) = - timestep*wls(k,l,nfor,i)*        &
     &                     (qcl(k,l,i)-qcl(k,l,i-1))/                   &
     &                     (r_theta_levs(k,l,i)-r_theta_levs(k,l,i-1))
                qcf_vertadv(k,l,i) = - timestep*wls(k,l,nfor,i)*        &
     &                     (qcf(k,l,i)-qcf(k,l,i-1))/                   &
     &                     (r_theta_levs(k,l,i)-r_theta_levs(k,l,i-1))
              Endif

            Endif ! w<0

            q_inc(k,l,nfor,i)   = q_inc(k,l,nfor,i)                     &
     &                                + q_vertadv(k,l,i)
            qcl_inc(k,l,nfor,i) = qcl_inc(k,l,nfor,i)                   &
     &                                + qcl_vertadv(k,l,i)
            qcf_inc(k,l,nfor,i) = qcf_inc(k,l,nfor,i)                   &
     &                                + qcf_vertadv(k,l,i)


          Endif  ! vertadv
!




          enddo                 ! i


          If (prindump_obs) then
            Do i = 1, nlevs
                dap1(k,l,10,i) = t_inc(k,l,nfor,i) / sec_day
            enddo               ! i
            Do i = 1, nwet
                dap1(k,l,20,i) = q_inc(k,l,nfor,i) *                    &
     &                                  1000.0 / sec_day
            enddo               ! i
          endif
        enddo
        enddo
      endif                     ! stats or obs

!-----------------------------------------------------------------------
!     SCM Forcing OR Increments Diagnostics Packages
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_forc)                                      &
     &    .OR. L_SCMDiags(SCMDiag_incs)) Then

        Do k=1, nlevs
          a2out_all(k) = t_inc(1,1,nfor,k) - t_inc_in(1,1,k)
        End Do
! DEPENDS ON: scmoutput
        Call SCMoutput(                                                 &
             a2out_all,                                                 &
             'dt_totforc',                                              &
             'Total temperature increment from s_forcng','K',           &
             t_avg,d_all,default_streams,'',RoutineName)

        Do k=1, nlevs
          a2out_all(k) = u_inc(1,1,nfor,k) - u_inc_in(1,1,k)
        End Do
! DEPENDS ON: scmoutput
        Call SCMoutput(                                                 &
             a2out_all,                                                 &
             'du_totforc',                                              &
             'Total u increment from s_forcng','m/s',                   &
             t_avg,d_all,default_streams,'',RoutineName)

        Do k=1, nlevs
          a2out_all(k) = v_inc(1,1,nfor,k) - v_inc_in(1,1,k)
        End Do
! DEPENDS ON: scmoutput
        Call SCMoutput(                                                 &
             a2out_all,                                                 &
             'dv_totforc',                                              &
             'Total v increment from s_forcng','m/s',                   &
             t_avg,d_all,default_streams,'',RoutineName)

        Do k=1, nwet
          a2out_wet(k) = q_inc(1,1,nfor,k) - q_inc_in(1,1,k)
        End Do
! DEPENDS ON: scmoutput
        Call SCMoutput(                                                 &
             a2out_wet,                                                 &
             'dq_totforc',                                              &
             'Total humidity increment from s_forcng','kg/kg',          &
             t_avg,d_wet,default_streams,'',RoutineName)

        Do k=1, nwet
          a2out_wet(k) = qcl_inc(1,1,nfor,k) - qcl_inc_in(1,1,k)
        End Do
! DEPENDS ON: scmoutput
        Call SCMoutput(                                                 &
             a2out_wet,                                                 &
             'dqcl_totforc',                                            &
             'Total QCL increment from s_forcng','kg/kg',               &
             t_avg,d_wet,default_streams,'',RoutineName)

        Do k=1, nwet
          a2out_wet(k) = qcf_inc(1,1,nfor,k) - qcf_inc_in(1,1,k)
        End Do
! DEPENDS ON: scmoutput
        Call SCMoutput(                                                 &
             a2out_wet,                                                 &
             'dqcf_totforc',                                            &
             'Total QCF increment from s_forcng','kg/kg',               &
             t_avg,d_wet,default_streams,'',RoutineName)
!
!-----------------------------------------------------------------------
!       SCM Forcing OR Increments Diagnostics Packages
!       when observational based forcing specified
!-----------------------------------------------------------------------
        If (obs) Then

! DEPENDS ON: scmoutput
          Call SCMoutput(t_inc_obs,                                     &
               'dt_obsforc',                                            &
               'Temperature increment from observational forcing','K',  &
               t_avg,d_all,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
          Call SCMoutput(u_inc_obs,                                     &
               'du_obsforc',                                            &
               'U increment from observational forcing','m/s',          &
               t_avg,d_all,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
          Call SCMoutput(v_inc_obs,                                     &
               'dv_obsforc',                                            &
               'V increment from observational forcing','m/s',          &
               t_avg,d_all,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
          Call SCMoutput(q_inc_obs,                                     &
               'dq_obsforc',                                            &
               'Humidity increment from observational forcing','kg/kg', &
               t_avg,d_wet,default_streams,'',RoutineName)

        End If ! obs

!-----------------------------------------------------------------------
!       SCM Forcing OR Increments Diagnostics Packages
!       when vertical advection specified
!-----------------------------------------------------------------------
        If (L_vertadv) Then

! DEPENDS ON: scmoutput
          Call SCMoutput(t_vertadv,                                     &
               'dt_vertadv',                                            &
               'Temperature increment from vertical advection','K',     &
               t_avg,d_all,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
          Call SCMoutput(q_vertadv,                                     &
               'dq_vertadv',                                            &
               'Humidity increment from vertical advection','kg/kg',    &
               t_avg,d_wet,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
          Call SCMoutput(qcl_vertadv,                                   &
               'dqcl_vertadv',                                          &
               'QCL increment from vertical advection','kg/kg',         &
               t_avg,d_wet,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
          Call SCMoutput(qcf_vertadv,                                   &
               'dqcf_vertadv',                                          &
               'QCF increment from vertical advection','kg/kg',         &
               t_avg,d_wet,default_streams,'',RoutineName)

        End If ! L_vertadv

      End If ! L_SCMDiags(SCMDiag_forc).OR L_SCMDiags(SCMDiag_incs)

      Return
      END SUBROUTINE FORCING

#endif
