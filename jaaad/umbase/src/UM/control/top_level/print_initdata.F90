#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!-----SUBROUTINE PRINT_INITDATA---------------------------------------
!
! Version 1 written by: J.LEAN, 10/6/91.
!
! Code reviewed by: ??
!
!
!
!*********************************************************************
!     *                                                            *
!     *  Single Column Unified Model version.                      *
!     *                                                            *
!     **************************************************************
!
!
!     Purpose: To print out initial run data for SCM integration
!
!     Documentation: ??
!
!---------------------------------------------------------------------
      SUBROUTINE PRINT_INITDATA(                                        &
        row_length,rows, land_points, nlevs, nwet, nozone,              &
        nfor, nbl_levs, nsoilt_levs, nsoilm_levs, ntrop,                &
        year_init, dayno_init, ndayin, nminin, nsecin,                  &
        timestep, ntrad, a_sw_radstep_prog,a_sw_radstep_diag,           &
        lat, long, ancyc, local_time, change_clim, exname_in,           &
        runno_in, exname_out, runno_out, soil_type, land_mask, obs,     &
        geoforce, geoinit, stats, noforce, tapein, tapeout, smi_opt,    &
        smcli, fsmc, sth, ug,vg, conv_mode, altdat, tls, tstar_forcing, &
        qls, uls, vls, wls, ichgf, ilscnt, flux_h, flux_e,              &
        ui, vi, wi, ti, qi, ccai, iccbi, iccti, p_in, canopy_gbi, smci, &
        snodepi, t_deep_soili, tstari, z0mseai, ntrad1,                 &
        radcloud_fixed, cca_rad, iccb_rad, icct_rad, ccwpin_rad,        &
        layer_cloud_rad, qcl_rad, time_init, o3, nout, numnout)

        Implicit None

!---------------------------------------------------------------------
!     Arguments
!---------------------------------------------------------------------


      Integer                                                           &
     &   row_length                                                     &
                                ! IN leading dimension of SCM arrays.
     &  ,rows                                                           &
     &  ,nlevs                                                          &
                                ! IN no of levels.
     &  ,nwet                                                           &
                                ! IN no of model levels in which Q
                                !   is set.
     &  ,nfor                                                           &
                                ! IN Number terms for observational
                                !   forcing
     &  ,nbl_levs                                                       &
                                ! IN Number of Boundary layer levels
     &  ,nsoilt_levs                                                    &
                                ! IN Number of soil temperature
                                !   levels
     &  ,nsoilm_levs                                                    &
                                ! IN Number of soil moisture levels
     &  ,ntrop                                                          &
                                ! IN Max number of levels in the
                                !   troposphere
     &  ,nozone                 ! IN no of model levels in which
                                !   ozone is set.

      Integer :: land_points    ! IN no of land points

      Logical                                                           &
     &   altdat                                                         &
                                ! T if alternative initial profiles
                                !  are to be input for STATS forcing
     &  ,ancyc                                                          & 
                                ! T if annual cycle req'd
                                !  (ie. radiation input then
                                !  varies throughout year)
     &  ,geoforce                                                       &
                                ! T if geostrophic forcing.
     &  ,geoinit                                                        &
                                ! T if initialising dump to
                                !  geostrophic.
     &  ,land_mask(row_length, rows)                                    &
                                ! T for a land point
     &  ,local_time                                                     &
                                ! T if diagnostics required
                                !  for local time rather than GMT
     &  ,noforce                                                        &
                                ! T if no large-scale
                                !  forcing used
     &  ,obs                                                            &
                                ! T if observational
                                !  large-scale forcing used
     &  ,radcloud_fixed                                                 &
                                ! T if clouds to fixed to input values
                                !  for radiation
     &  ,stats                                                          &
                                ! T if statistical
                                !  large-scale forcing used
     &  ,tapein                                                         &
                                ! T if initial data is to be read
                                !  from previous run stored on
                                !  tape
     &  ,tapeout                ! T if restart information
                                !  plus diagnostic output to be
                                !  stored on tape

      Integer :: smi_opt        ! Selects method used to initialise soil
                                ! Moisture content
      Integer                                                           &
     &  change_clim                                                     &
                                ! No. of days between
                                !  changes of climatological data
                                !  (default=10)
     &  ,dayno_init                                                     &
                                ! Initial day (assume 360
                                !  day year)
     &  ,year_init                                                      &
                                ! Initial year
     &  ,iccbi(row_length, rows)                                        &
                                    ! Initial convective cloud base.
     &  ,iccb_rad(row_length, rows)                                     &
                                    ! Fixed convective cloud base.
     &  ,iccti(row_length, rows)                                        &
                                    ! Initial convective cloud top.
     &  ,icct_rad(row_length, rows)                                     &
                                    ! Fixed convective cloud base.
     &  ,ichgf                                                          &
                                ! No. of timesteps between change in
                                !  observational forcing
     &  ,ilscnt                                                         &
                                ! Counts for change in
                                !  observational forcing
     &  ,runno_in                                                       &
                                ! Number of run to be read
                                !  from previous run stored on tape
     &  ,runno_out                                                      &
                                ! Number of run to be written
                                !  to tape
     &  ,soil_type(land_points)                                         &
                                ! Soil type code
                                !  1 Ice
                                !  2 Fine
                                !  3 Medium
                                !  4 Coarse
     &  ,ndayin                                                         &
                                ! No. of days in run
     &  ,nminin                                                         &
                                ! No. minutes in run
     &  ,nsecin                                                         &
                                ! No. seconds in run
     &  ,ntrad                                                          &
                                ! No. of timestep
     &  ,a_sw_radstep_diag                                              &
                                ! No. of diagnostic timesteps
     &  ,a_sw_radstep_prog                                              &
                                ! No. of prognostic timesteps
                                !  between calls to radiation
     &  ,ntrad1                                                         &
                                ! st timestep radiation called.
     &  ,numnout                                                        &
                                ! No. units for output
     &  ,nout(numnout)          ! Units for output

      Integer                                                           &
     &  conv_mode               ! Mode to run convection
     
      Real                                                              &
     &  canopy_gbi(land_points)                                         &
                                ! Initial canopy water content (kg/m2)
     &, ccai(row_length, rows)                                          &
                                ! Initial convective cloud amount
     &, cca_rad(row_length, rows)                                       &
                                ! Fixed convective cloud amount
     &, ccwpin_rad(row_length, rows)                                    &
                                ! Fixed convective cloud water path
     &  ,flux_h(row_length, rows,nfor)                                  &
                                ! OBS calculation used for sensible
                                !  heat flux calculation
     &  ,flux_e(row_length, rows,nfor)                                  &
                                ! OBS calculation used for
                                !  evaporative flux calculation
     &  ,fsmc(land_points)                                              &
                                ! MOSES soil moisture stress factor
     &  ,lat(row_length, rows)                                          &
                                ! Lat. of gridpoiny chosen; read
                                !  automatically from climate
                                !  dataset if STATS forcing chosen
     &  ,layer_cloud_rad(row_length, rows,nwet)                         &
                                ! Layer cloud amount (fraction)
     &  ,long(row_length, rows)                                         &
                                ! Long. of gridpoint chosen;
                                !  read automatically from climate
                                !  dataset if STATS forcing chosen
     &  ,o3(row_length, rows,nozone)                                    &
                                ! Ozone values
     &  ,p_in(row_length, rows, nlevs+1)                                &
                                        ! Inital pressure.
     &  ,qcl_rad(row_length, rows,nwet)                                 &
                                ! Total cloud water and ice content in
                                !  cloud
     &  ,snodepi(row_length, rows)                                      &
                                ! Initial snow depth.
     &  ,smci(land_points)                                              &
                                ! Initial soil moisture content.
     &  ,smcli(land_points, nsoilm_levs)                                &
                                ! MOSES soil moisture in layers
                                !  initial value
     &  ,sth(land_points, nsoilm_levs)                                  &
                                ! MOSES soil moisture in layers as
                                !  a fraction of saturation.
     &  ,timestep                                                       &
                                ! Model timestep for all physics
                                !  subroutines except radiation
     &  ,ti(row_length, rows,nlevs)                                     &
                                ! Initial temperature profile
     &  ,time_init                                                      &
                                ! Start time of run(secs)
     &  ,tstari(row_length, rows)                                       &
                                 ! Initial surface temperature
     &  ,t_deep_soili(land_points, nsoilt_levs)                         &
                                ! Inital deep soil temperature
     &  ,qi(row_length, rows,nwet)                                      &
                                ! Initial specific humidity profile
     &  ,ug(row_length, rows)                                           &
                                ! Geostrophic U velocity (m/s)
     &  ,vg(row_length, rows)                                           &
                                ! Geostrophic V velocity (m/s)
     &  ,ui(row_length, rows,nlevs)                                     &
                                    ! Initial zonal wind profile
     &  ,vi(row_length, rows,nlevs)                                     &
                                    ! Initial meridional wind profile
     &  ,wi(row_length, rows,nlevs+1)                                   &
                                      ! Initial meridional wind profile
     &  ,tls(row_length, rows,nfor,nlevs)                               &
                                ! Temperature increment for OBS
                                !  forcing
     &  ,tstar_forcing(row_length, rows,nfor)                           &
                                ! Forcing for sea surface temperature
     &  ,qls(row_length, rows,nfor,nwet)                                &
                                ! Specific humidity increment for
                                !  OBS forcing
     &  ,uls(row_length, rows,nfor,nlevs)                               &
                                ! Zonal wind increment for OBS forcing
     &  ,vls(row_length, rows,nfor,nlevs)                               &
                                ! Meridional wind increment for OBS
                                ! forcing.
     &  ,wls(row_length, rows,nfor,nlevs)                               &
                                ! Vertical wind increment for OBS
                                ! forcing.
     &  ,z0mseai(row_length, rows)! Initial sea surface roughness
                                !  length.
      Character*8 exname_in,exname_out

!     Local Variables

      Integer                                                           &
     &  i, j, k, l, m, land_cnt                                         &
                                ! Counters
     &  ,multnl,resnl                                                   &
                                ! Counters for printing nlevs
                                !  variables
     &  ,multnf,resnf                                                   &
                                ! Counter for printing nfor variables
     &  ,multsl,ressl                                                   &
                                ! Counter for printing nsoil variables
     &  ,multoz,resoz                                                   &
                                ! Counter for printing nozone
                                !  variables
     &  ,multcl,rescl

      Character (Len=100) :: c1fmt, c2fmt, c3fmt, c4fmt, c5fmt
      Character (Len=100) :: c6fmt, c7fmt, c8fmt, c10fmt
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      Character (Len=100) :: c11fmt1, c11fmt2
      Character (Len=100) :: c12fmt, c13fmt, c14fmt, c15fmt
#else
      Character (Len=100) :: c11fmt, c12fmt, c13fmt, c14fmt, c15fmt
#endif
      Character (Len=100) :: c16fmt, c17fmt, c18fmt, c19fmt, c20fmt
      Character (Len=100) :: c31fmt, c32fmt, c33fmt, c34fmt, c35fmt
      Character (Len=100) :: c36fmt, c37fmt, c38fmt, c39fmt, c40fmt
      Character (Len=100) :: c41fmt, c42fmt, c43fmt, c44fmt, c45fmt
      Character (Len=100) :: c46fmt, c47fmt, c48fmt, c49fmt, c50fmt
      Character (Len=100) :: c51fmt, c52fmt, c53fmt, c54fmt, c55fmt
      Character (Len=100) :: c56fmt, c57fmt, c58fmt, c59fmt, c60fmt
      Character (Len=100) :: c61fmt, c62fmt, c63fmt, c64fmt, c65fmt
      Character (Len=100) :: c66fmt, c67fmt, c68fmt, c69fmt, c70fmt
      Character (Len=100) :: c71fmt, c72fmt, c73fmt, c74fmt, c75fmt
      Character (Len=100) :: c76fmt, c77fmt, c78fmt, c79fmt, c80fmt
      Character (Len=100) :: c81fmt, c82fmt, c83fmt, c84fmt, c85fmt
      Character (Len=100) :: c86fmt, c87fmt, c88fmt, c89fmt, c90fmt
      Character (Len=100) :: c91fmt, c92fmt, c93fmt, c94fmt, c95fmt
      Character (Len=100) :: c96fmt, c97fmt, c98fmt, c99fmt, c100fmt
      Character (Len=100) :: c101fmt, c102fmt, c103fmt

!---------------------------------------------------------------------
!     Set up format statements
!---------------------------------------------------------------------

      c1fmt = '('                                                       &
     &  //'" * Initial data used for scm integration      ",27x,"*")'
      c2fmt = '('                                                       &
     & //'"******************************************************"'     &
     & //',"*******************")'
      c3fmt = '('                                                       &
     & //'" * Model run is for a land point              ",26x," *")'
      c4fmt = '('                                                       &
     & //'" * Model run is for a sea point               ",26x," *")'
      c5fmt = '('                                                       &
     & //'" * Run started on year                        ",i26," *")'
      c6fmt = '('                                                       &
     & //'" * Run started on day                         ",i26," *")'
      c7fmt = '('                                                       &
     & //'" * Run started at time (secs)               ",f28.2," *")'
      c8fmt = '('                                                       &
     & //'" * length of run (days:mins:secs)",13x,i8,":",i8,":",'       &
     & //'i8, " *")'
      c10fmt = '('                                                      &
     & //'" * Timestep (seconds)                       ",f28.1," *")'
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      c11fmt1 = '('                                                     &
     & //'" * Prognostic Timestep for radiation (secs) ",f28.1," *")'
      c11fmt2 = '('                                                     &
     & //'" * Diagnostic Timestep for radiation (secs) ",f28.1," *")'
#else
      c11fmt = '('                                                      &
     & //'" * Timestep for radiation (seconds)         ",f28.1," *")'
#endif
      c12fmt = '('                                                      &
     & //'" * 1st timestep for radiation (seconds)     ",f28.1," *")'
      c13fmt = '('                                                      &
     & //'" * Latitude of gridpoint                    ",f28.1," *")'
      c14fmt = '('                                                      &
     & //'" * Longitude of gridpoint                   ",f28.1," *")'
      c15fmt = '('                                                      &
     & //'" * Number of model levels                     ",i26," *")'
      c16fmt = '('                                                      &
     & //'" * Number of atmospheric moist levels         ",i26," *")'
      c17fmt = '('                                                      &
     & //'" * Number of levels in boundary layer         ",i26," *")'
      c18fmt = '('                                                      &
     & //'" * Number of deep soil temperatures           ",i26," *")'
      c19fmt = '('                                                      &
     & //'" * Convection switched off                    ",26x," *")'
      c20fmt = '('                                                      &
     & //'" * Convection switched on for diagnostics only",26x," *")'
      c31fmt = '('                                                      &
     & //'" * Annual cycle included                      ",26x," *")'
      c32fmt = '('                                                      &
     & //'" * Annual cycle not included                  ",26x," *")'
      c33fmt = '('                                                      &
     & //'" * All diagnostics refer to local time        ",26x," *")'
      c34fmt = '('                                                      &
     & //'" * All diagnostics refer to gmt               ",26x," *")'
      c35fmt = '('                                                      &
     & //'" * Large-scale observational forcing selected ",26x," *")'
      c36fmt = '('                                                      &
     & //'" * Number of timesteps between change in forcing",i25,"*")'
      c37fmt = '('                                                      &
     & //'" * Counter for forcing                        ",i26," *")'
      c38fmt = '('                                                      &
     & //'" * Number of terms for observational forcing  ",i26," *")'
      c39fmt = '('                                                      &
     & //'" * Temperature increment tls (k)              ",26x," *")'
      c40fmt = '('                                                      &
     & //'" * Zonal wind increment uls     (m/s)         ",26x," *")'
      c41fmt = '('                                                      &
     & //'" * Meridional wind increment vls  (m/s)       ",26x," *")'
      c42fmt = '('                                                      &
     & //'" * Specific humidity increment qls (kg/kg)    ",26x," *")'
      c43fmt = '('                                                      &
     & //'" * Flux_h - observational forcing             ",26x," *")'
      c44fmt = '('                                                      &
     & //'" * Flux_e - observational forcing             ",26x," *")'
      c45fmt = '('                                                      &
     & //'" * Large-scale statistical forcing selected   ",26x," *")'
      c46fmt = '('                                                      &
     & //'" * Number of days between change in forcing   ",i26," *")'
      c47fmt = '('                                                      &
     & //'" * No large-scale forcing selected            ",26x," *")'
      c48fmt = '('                                                      &
     & //'" * Geostrophic forcing selected               ",26x," *")'
      c49fmt = '('                                                      &
     & //'" * Geostrophic initialisation selected        ",26x," *")'
      c50fmt = '('                                                      &
     & //'" * Initial zonal wind profile ui  (m/s)       ",26x," *")'
      c51fmt = '('                                                      &
     & //'" * Initial meridional wind profile vi  (m/s)  ",26x," *")'
      c52fmt = '('                                                      &
     & //'" * Initial geostrophic zonal wind ug (m/s)"'                 &
     & //',20x,f10.3," *")'
      c53fmt = '('                                                      &
     & //'" * Initial geostrophic meridional wind vg (m/s)",15x,'       &
     & //' f10.3," *")'
      c54fmt = '('                                                      &
     & //'" * Initial temperature profile ti  (k)        ",26x," *")'
      c55fmt = '('                                                      &
     & //'" * Initial specific humidity profile qi (kg/kg)",26x,"*")'
      c56fmt = '('                                                      &
     & //'" * Ozone                                      ",26x," *")'
      c57fmt = '('                                                      &
     & //'" * Initial convective cloud amount ccai     ",f28.4," *")'
      c58fmt = '('                                                      &
     & //'" * Initial convective cloud base iccbi        ",i26," *")'
      c59fmt = '('                                                      &
     & //'" * Initial convective cloud top iccti         ",i26," *")'
      c60fmt = '('                                                      &
     & //'" * Initial pressure p_in     ",10f7.5," *")'
      c61fmt = '('                                                      &
     & //'" * Fixed radiation conv. cloud amount cca_rad",f27.4," *")'
      c62fmt = '('                                                      &
     & //'" * Fixed radiation conv. cloud base iccb_rad  ",i26," *")'
      c63fmt = '('                                                      &
     & //'" * Fixed radiation conv. cloud top icct_rad   ",i26," *")'
      c64fmt = '('                                                      &
     & //'" * Fixed radiation conv. cloud ccwpin_rad    ",f27.4," *")'
      c65fmt = '('                                                      &
     & //'" * Fixed radiation layer cloud amount          ",25x," *")'
      c66fmt = '('                                                      &
     & //'" * Fixed radiation total water and ice in cloud",25x," *")'
      c67fmt = '('                                                      &
     & //'" * Initial canopy water content kg m-2      ",f28.3," *")'
      c68fmt = '('                                                      &
     & //'" * Initial soil moisture content kg m-2     ",f28.3," *")'
      c69fmt = '('                                                      &
     & //'" * Initial snow depth kg m-2                ",f28.3," *")'
      c70fmt = '('                                                      &
     & //'" * Initial deep soil temperatures (k)",6x,"*")'
      c71fmt = '('                                                      &
     & //'" * Initial surface temperature (K)          ",f28.3," *")'
      c72fmt = '('                                                      &
     & //'" * Moses initialised by input of actual soil "'              &
     & //'"moisture in ,"layers        *")'
      c73fmt = '('                                                      &
     & //'" * Initial smcl  (Kg m-2)                      ",26X,"*")'
      c74fmt = '('                                                      &
     & //'"* Moses initialisation by input of soil moisture ",'         &
     & //'"stress factor    *")'
      c75fmt = '('                                                      &
     & //'" * Initial soil moisture stress factor      ",f28.3," *")'
      c76fmt = '('                                                      &
     & //'" * Moses initialised by input of smcl as a fraction of ",'   &
     & //'"saturation        *")'
      c77fmt = '('                                                      &
     & //'" * Initial sth                                 ",26x,"*")'
      c78fmt = '('                                                      &
     & //'" * Initial sea surface roughness length     ",f28.3," *")'
      c79fmt = '('                                                      &
     & //'" * Initial data read from tape                ",26x," *")'
      c80fmt = '('                                                      &
     & //'" * With experment name                        ",a26," *")'
      c81fmt = '('                                                      &
     & //'" * With run number                            ",i26," *")'
      c82fmt = '('                                                      &
     & //'" * No data read from tape                     ",26x," *")'
      c83fmt = '('                                                      &
     & //'" * Data written to tape                       ",26x," *")'
      c84fmt = '('                                                      &
     & //'" * With experment name                        ",a26," *")'
      c85fmt = '('                                                      &
     & //'" * With run number                            ",i26," *")'
      c86fmt = '('                                                      &
     & //'" * No data written to tape                    ",26x," *")'
      c96fmt = '('                                                      &
     & //'" * Soil type ice                              ",26x," *")'
      c97fmt = '('                                                      &
     & //'" * Soil type fine                             ",26x," *")'
      c98fmt = '('                                                      &
     & //'" * Soil type medium                           ",26x," *")'
      c99fmt = '('                                                      &
     & //'" * Soil type coarse                           ",26x," *")'
      c100fmt = '('                                                     &
     & //'" * Initial soil temperature profile derived   ",26x," *")'
      c101fmt = '('                                                     &
     & //'" * From soil temp. cycle with mean annual   ",f28.1," *")'
      c102fmt = '('                                                     &
     & //'" * Daily amplitude of soil temperature cycle",f28.1," *")'
      c103fmt = '('                                                     &
     & //'" * Annual amplitude of soil temperature cycle",f27.1," *")'

!---------------------------------------------------------------------
!     Output the initial data to all the unit nos. where diagnostics
!     are required
!---------------------------------------------------------------------

      land_cnt = 0 

      Do l = 1, rows          ! loop on points
        Do k = 1, row_length
          
          If (land_mask(k,l)) Then
            land_cnt = land_cnt + 1
          End If 

          If (row_length*rows  >  1 )                                   &
     &      Write (nout(numnout),*) " COLUMN NUMBER : (",k,l, ")"
          If (nout(numnout)  /=  0) then

            Write (nout(numnout),c1fmt)
            Write (nout(numnout),c2fmt)
            If  (land_mask(k,l)) then
              Write (nout(numnout),c3fmt)
            else
              Write (nout(numnout),c4fmt)
            endif
            Write (nout(numnout),c5fmt) year_init
            Write (nout(numnout),c6fmt) dayno_init
            Write (nout(numnout),c7fmt) time_init
            Write (nout(numnout),c8fmt) ndayin, nminin, nsecin
            Write (nout(numnout),c10fmt) timestep
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
            Write (nout(numnout),c11fmt1) timestep*a_sw_radstep_prog
            Write (nout(numnout),c11fmt2) timestep*a_sw_radstep_diag
#else
            Write (nout(numnout),c11fmt) timestep*ntrad
#endif
            Write (nout(numnout),c12fmt) timestep*ntrad1
            Write (nout(numnout),c13fmt) lat(k,l)
            Write (nout(numnout),c14fmt) long(k,l)
            Write (nout(numnout),c15fmt) nlevs
            Write (nout(numnout),c16fmt) nwet
            Write (nout(numnout),c17fmt) nbl_levs
            If (land_mask(k,l)) write (nout(numnout),c18fmt) nsoilt_levs
!           If any physics routines are switched off or switched to
!           diagnostics only print out here
            If (conv_mode  ==  2) then
              Write (nout(numnout),c19fmt)
            elseif (conv_mode  ==  1) then
              Write (nout(numnout),c20fmt)
            endif

!           ancyc
            If (ancyc) then
              Write (nout(numnout),c31fmt)
            else
              Write (nout(numnout),c32fmt)
            endif
!           local time
            If (local_time) then
              Write (nout(numnout),c33fmt)
            else
              Write (nout(numnout),c34fmt)
            endif
!           obs
            If (obs) then
              Write (nout(numnout),c35fmt)
              Write (nout(numnout),c36fmt) ichgf
              Write (nout(numnout),c37fmt) ilscnt
              Write (nout(numnout),c38fmt) nfor

              multnf = int(nfor/10)
              resnf = mod(nfor,10)
!             tls
              Write (nout(numnout),c39fmt)
              If (multnf  /=  0) then
                Do j = 1, nlevs
                  Do m = 1, multnf
                    Write (nout(numnout),'(10f7.2)')                    &
     &                (tls(k,l,i,j), i = (m-1)*10+1, 10*m)
                  enddo
                  Write (nout(numnout),'(10f7.2)')                      &
     &             (tls(k,l,i,j), i = (m-1)*10+1, (m-1)*10+resnf)
                enddo
              else
                Do j = 1, nlevs
                  Write (nout(numnout),'(10f7.2)')                      &
     &              (tls(k,l,i,j), i = 1, resnf)
                enddo
              endif               ! multnf
!             tstar_forcing
              Write (nout(numnout),c39fmt)
              If (multnf  /=  0) then
                Do m = 1, multnf
                  Write (nout(numnout),'(10f7.2)')                      &
     &              (tstar_forcing(k,l,i), i = (m-1)*10+1, 10*m)
                enddo
                Write (nout(numnout),'(10f7.2)')                        &
     &           (tstar_forcing(k,l,i), i = (m-1)*10+1,                 &
     &                                             (m-1)*10+resnf)
              else
                Write (nout(numnout),'(10f7.2)')                        &
     &            (tstar_forcing(k,l,i), i = 1, resnf)
              endif               ! multnf
!             uls
              Write (nout(numnout),c40fmt)
              If (multnf  /=  0) then
                Do j = 1, nlevs
                  Do m = 1, multnf
                    Write (nout(numnout),'(10f7.2)')                    &
     &              (uls(k,l,i,j), i = (m-1)*10+1, 10*m)
                  enddo
                  Write (nout(numnout),'(10f7.2)')                      &
     &            (uls(k,l,i,j), i = (m-1)*10+1, (m-1)*10+resnf)
                enddo
              else
                Do j = 1, nlevs
                  Write (nout(numnout),'(10f7.2)')                      &
     &            (uls(k,l,i,j), i = 1, resnf)
                enddo
              endif               ! multnf
!             vls
              Write (nout(numnout),c41fmt)
              If (multnf  /=  0) then
                Do j = 1, nlevs
                  Do m = 1, multnf
                    Write (nout(numnout),'(10f7.2)')                    &
     &              (vls(k,l,i,j),i=(m-1)*10+1,10*m)
                  enddo
                  Write (nout(numnout),'(10f7.2)')                      &
     &            (vls(k,l,i,j), i = (m-1)*10+1, (m-1)*10+resnf)
                enddo
              else
                Do j = 1, nlevs
                  Write (nout(numnout),'(10f7.2)')                      &
     &             (vls(k,l,i,j),i=1,resnf)
                enddo
              endif               !multnf
!             wls
              Write (nout(numnout),c41fmt)
              If (multnf  /=  0) then
                Do j = 1, nlevs
                  Do m = 1, multnf
                    Write (nout(numnout),'(10f7.2)')                    &
     &              (wls(k,l,i,j),i=(m-1)*10+1,10*m)
                  enddo
                  Write (nout(numnout),'(10f7.2)')                      &
     &            (wls(k,l,i,j), i = (m-1)*10+1, (m-1)*10+resnf)
                enddo
              else
                Do j = 1, nlevs
                  Write (nout(numnout),'(10f7.2)')                      &
     &             (wls(k,l,i,j),i=1,resnf)
                enddo
              endif               !multnf
!             qls
              Write (nout(numnout),c42fmt)
              multnf = int(nfor/5)
              resnf = mod(nfor,5)
              If (multnf  /=  0) then
                Do j = 1, nwet
                  Do m = 1, multnf
                    Write (nout(numnout),'(5e10.2)')                    &
     &              (qls(k,l,i,j),i=(m-1)*5+1,5*m)
                  enddo
                  Write (nout(numnout),'(5e10.2)')                      &
     &            (qls(k,l,i,j), i = (m-1)*5+1, (m-1)*5+resnf)
                enddo
              else
                Do j = 1, nlevs
                  Write (nout(numnout),'(5e10.2)')                      &
     &              (qls(k,l,i,j),i=1,resnf)
                enddo
              endif               ! multnf
!             flux_h
              Write (nout(numnout),c43fmt)
              multnf = int(nfor/5)
              resnf = mod(nfor,5)
              If (multnf  /=  0) then
                Do m = 1, multnf
                  Write (nout(numnout),'(5f10.5)')                      &
     &            (flux_h(k,l,i),i=(m-1)*5+1,5*m)
                enddo
                Write (nout(numnout),'(5f10.5)')                        &
     &          (flux_h(k,l,i), i=(m-1)*5+1,(m-1)*5+resnf)
              else
                Write (nout(numnout),'(5f10.5)')                        &
     &          (flux_h(k,l,i),i=1,resnf)
              endif               !multnf
!             flux_e
              Write (nout(numnout),c44fmt)
              If (multnf  /=  0) then
                Do m = 1, multnf
                  Write (nout(numnout),'(5f10.5)')                      &
     &            (flux_e(k,l,i),i=(m-1)*5+1,5*m)
                enddo
                Write (nout(numnout),'(5f10.5)')                        &
     &            (flux_e(k,l,i), i = (m-1)*5+1, (m-1)*5+resnf)
              else
                Write (nout(numnout),'(5f10.5)')                        &
     &          (flux_e(k,l,i),i=1,resnf)
              endif               ! multnf
            endif                 ! obs
!           STATS
            If (stats) then
              Write (nout(numnout),c45fmt)
              Write (nout(numnout),c46fmt) change_clim

            endif                 ! stats

!           noforce

            If (noforce) then
              Write (nout(numnout),c47fmt)
            endif                 ! noforce

!           geoforce

            If (geoforce) then
              Write (nout(numnout),c48fmt)
              If (geoinit) then
                Write (nout(numnout),c49fmt)
              endif
            endif                 ! goforce

!           obs, noforce, altdat

            If (noforce .or.  obs .or. altdat .or. geoforce) then
              multnl = int((nlevs+1)/10)
              resnl = mod((nlevs+1),10)
              If (noforce .or. obs                                      &
     &         .or. (geoforce .and. .not. geoinit)) then
!               wi
                Write (nout(numnout),c50fmt)
                If (multnl /= 0) then
                  Do j=1,multnl
                    Write (nout(numnout),'("  ",10f7.2)')               &
     &              (wi(k,l,i), i = 10*(j-1)+1, 10*j)
                  enddo
                  Write (nout(numnout),'("  ",10f7.2)')                 &
     &            (wi(k,l,i), i = 10*(j-1)+1, 10*(j-1)+resnl)
                else
                  Write (nout(numnout),'("  ",10f7.2)')                 &
     &              (wi(k,l,i),i=1,resnl)
                endif             ! multnl
              endif                 !

              multnl = int(nlevs/10)
              resnl = mod(nlevs,10)

              If (noforce .or. obs                                      &
     &         .or. (geoforce .and. .not. geoinit)) then
!               ui
                Write (nout(numnout),c50fmt)
                If (multnl /= 0) then
                  Do j=1,multnl
                    Write (nout(numnout),'("  ",10f7.2)')               &
     &              (ui(k,l,i), i = 10*(j-1)+1, 10*j)
                  enddo
                  Write (nout(numnout),'("  ",10f7.2)')                 &
     &            (ui(k,l,i), i = 10*(j-1)+1, 10*(j-1)+resnl)
                else
                  Write (nout(numnout),'("  ",10f7.2)')                 &
     &              (ui(k,l,i),i=1,resnl)
                endif             ! multnl
!               vi
                Write (nout(numnout),c51fmt)
                If (multnl  /=  0) then
                  Do j = 1, multnl
                    Write (nout(numnout),'("  ",10f7.2)')               &
     &                (vi(k,l,i), i = 10*(j-1)+1, 10*j)
                  enddo
                  Write (nout(numnout),'("  ",10f7.2)')                 &
     &              (vi(k,l,i), i = 10*(j-1)+1, 10*(j-1)+resnl)
                else
                  Write (nout(numnout),'("  ",10f7.2)')                 &
     &              (vi(k,l,i), i = 1, resnl)
                endif             ! multnl

              elseif (geoforce .and. geoinit) then
!               ug,vg
                Write (nout(numnout),c52fmt) ug(k,l)
                Write (nout(numnout),c53fmt) vg(k,l)

              endif               ! obs or noforce or geoforce

!             ti

              Write (nout(numnout),c54fmt)
              If (multnl  /=  0) then
                Do j=1,multnl
                  Write (nout(numnout),'("  ",10f7.2)')                 &
     &              (ti(k,l,i),i=10*(j-1)+1,10*j)
                enddo
                Write (nout(numnout),'("  ",10f7.2)')                   &
     &            (ti(k,l,i),i=10*(j-1)+1,10*(j-1)+resnl)
              else
                Write (nout(numnout),'("  ",10f7.2)')                   &
     &          (ti(k,l,i),i=1,resnl)
              endif               !multnl
!             qi
              Write (nout(numnout),c55fmt)
              multnl = int(nwet/5)
              resnl = mod(nwet,5)
              If (multnl  /=  0) then
                Do j = 1, multnl
                  Write (nout(numnout),'("  ",5e11.5)')                 &
     &             (qi(k,l,i),i=5*(j-1)+1,5*j)
                enddo
                Write (nout(numnout),'("  ",5e11.5)')                   &
     &           (qi(k,l,i),i=5*(j-1)+1,5*(j-1)+resnl)
              else
                Write (nout(numnout),'("  ",5e11.5)')                   &
     &            (qi(k,l,i),i=1,resnl)
              endif               ! multnl
            endif                 ! end if altdat,obs,noforce

!           Ozone
            Write (nout(numnout),c56fmt)
            multoz = int(nozone/5)
            resoz = mod(nozone,5)
            If (multoz  /=  0) then
              Do j = 1, multoz
                Write (nout(numnout),'("  ",5e11.5)')                   &
     &            (o3(k,l,i),i=5*(j-1)+1,5*j)
              enddo
              Write (nout(numnout),'("  ",5e11.5)')                     &
     &          (o3(k,l,i),i=5*(j-1)+1,5*(j-1)+resoz)
            else
              Write (nout(numnout),'("  ",5e11.5)')                     &
     &              (o3(k,l,i),i=1,resoz)
            endif                 ! multoz
            Write (nout(numnout),c57fmt) ccai(k,l)
            Write (nout(numnout),c58fmt) iccbi(k,l)
            Write (nout(numnout),c59fmt) iccti(k,l)
            multnl = int(nlevs+1/10)
            resnl = mod(nlevs+1,10)
            Do j = 1, multnl
              Write (nout(numnout),c60fmt)                              &
     &         (p_in(k,l,i), i = 10*(j-1)+1, 10*j)
            End Do
            multnl = int(nlevs/10)
            resnl = mod(nlevs,10)
!           radcloud_fixed
            If (radcloud_fixed) then
              Write (nout(numnout),c61fmt) cca_rad(k,l)
              Write (nout(numnout),c62fmt) iccb_rad(k,l)
              Write (nout(numnout),c63fmt) icct_rad(k,l)
              Write (nout(numnout),c64fmt) ccwpin_rad(k,l)

              multcl = int(nwet/5)
              rescl = mod(nwet,5)
!             layer_cloud_rad
              Write (nout(numnout),c65fmt)
              If (multcl  /=  0) then
                Do j = 1, multcl
                  Write (nout(numnout),'("  ",5e11.5)')                 &
     &            (layer_cloud_rad(k,l,i), i=5*(j-1)+1, 5*j)
                enddo
                Write (nout(numnout),'("  ",5e11.5)')                   &
     &          (layer_cloud_rad(k,l,i), i = 5*(j-1)+1, 5*(j-1)+rescl)
              else
                Write (nout(numnout),'("  ",5e11.5)')                   &
     &          (layer_cloud_rad(k,l,i), i = 1, rescl)
              endif
!             qcl_rad
              Write (nout(numnout),c66fmt)
              If (multcl  /=  0) then
                Do j = 1, multcl
                  Write (nout(numnout),'("  ",5e11.5)')                 &
     &            (qcl_rad(k,l,i), i = 5*(j-1)+1, 5*j)
                enddo
                Write (nout(numnout),'("  ",5e11.5)')                   &
     &          (qcl_rad(k,l,i),i=5*(j-1)+1,5*(j-1)+rescl)
              else
                Write (nout(numnout),'("  ",5e11.5)')                   &
     &            (qcl_rad(k,l,i), i = 1, rescl)
              endif
            endif                 ! radcloud fixed

            ! landmask
            If (land_mask(k,l)) Then

              Write (nout(numnout),c67fmt) canopy_gbi(land_cnt)
              Write (nout(numnout),c68fmt) smci(land_cnt)
              Write (nout(numnout),c69fmt) snodepi(k,l)

              multsl = int((nsoilt_levs)/10)
              ressl  = mod((nsoilt_levs),10)

              ! t_deep_soil
              Write (nout(numnout),c70fmt)
              If (multsl  /=  0) Then

                Do j = 1, multsl
                  Write (nout(numnout),'("  ",10f7.2)')                       &
                        (t_deep_soili(land_cnt,i),i = 10*(j-1)+1, 10*j)
                End Do

                Write (nout(numnout),'("  ",10f7.2)')                         &
                      (t_deep_soili(land_cnt,i), i = 10*(j-1)+1,10*(j-1)+ressl)
              Else
                Write (nout(numnout),'("  ",10f7.2)')                         &
                      (t_deep_soili(land_cnt,i), i = 1, ressl)
              End If

              Write (nout(numnout),c71fmt) tstari(k,l)

            Else

              Write (nout(numnout),c78fmt) z0mseai

            End If      ! land point

            If (tapein) then
              Write (nout(numnout),c79fmt)
              Write (nout(numnout),c80fmt) exname_in
              Write (nout(numnout),c81fmt) runno_in
            else
              Write (nout(numnout),c82fmt)
            endif
            If (tapeout) then
              Write (nout(numnout),c83fmt)
              Write (nout(numnout),c84fmt) exname_out
              Write (nout(numnout),c85fmt) runno_out
            else
              Write (nout(numnout),c86fmt)
            endif

            If (land_mask(k,l)) Then
              If (soil_type(land_cnt) ==  1) Then
                Write (nout(numnout),c96fmt)
              Else If (soil_type(land_cnt)  ==  2) Then
                Write (nout(numnout),c97fmt)
              Else If (soil_type(land_cnt)  ==  3) Then
                Write (nout(numnout),c98fmt)
              Else If (soil_type(land_cnt)  ==  4) Then
                Write (nout(numnout),c99fmt)
              End If
            End If

            Write (nout(numnout),c2fmt)
          endif                   ! nout
        enddo                     ! k
      End do                      ! l

      Return
      END SUBROUTINE PRINT_INITDATA
#endif
