#if defined(OASIS3) || defined(OASIS4)
#define ACCESS_MASK1 1
#define ACCESS_SOLAR 1
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
      Subroutine oasis_inita2o(                                        &
#include "argd1.h"
#include "argsts.h"
#include "argduma.h"
#include "argptra.h"
     & icode,                                                           &
     & cmessage)

      USE oasis3_atm_data_mod
#if defined(ACCESS)
      USE auscom_cpl_data_mod
#endif

#if defined(ACCESS_SOLAR)
      USE atm_fields_mod
#endif

      Implicit none

!
! Description: Initialise pointers to STASH diagnostic arrays
!              which are used when coupling with external ocean 
!              and ice models.
!
! Author: R. Hill
! Current Code Owner : R. Hill 
!
!--------------------------------------------------------------------
#include "csubmodl.h"
#include "cntlatm.h"
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
#include "typsts.h"
#include "typduma.h"
#include "typptra.h"
#include "atm_lsm.h"
!
      Integer       icode       ! OUT - Error return code
      Character*(*) cmessage    ! OUT - Error return message
!*----------------------------------------------------------------------
!
#include "caoptr.h"
#include "c_mdi.h"
#include "stparam.h"
!
#if defined(ACCESS)
#include "c_0_dg_c.h"
#endif
!  Local variables
!
      Integer :: process_code   ! Processing code
      Integer :: freq_code      ! Frequency code
      Integer :: start,end,period ! Start, end and period step
      Integer :: gridpt_code,weight_code ! Gridpt and weighting codes
      Integer :: bottom_level,top_level ! Bottom and top input level
      Integer :: grid_n,grid_s,grid_w,grid_e ! Grid corner definitions
      Integer :: stashmacro_tag ! STASHmacro tag number
      Integer :: FIELD_CODE     ! Temporary STASH field code holder

      Integer :: nlandpt        ! Dummy argument for FROM_LAND_POINTS
                                ! output field.
      Integer :: I, J           ! Local loop indices
      INTEGER :: im_ident,im_index
!----------------------------------------------------------------------
!     Set grid definition information (undefined as search is on
!     STASHmacro tag number)
!
      process_code=imdi
      freq_code   =imdi
      start       =imdi
      end         =imdi
      period      =imdi
      gridpt_code =imdi
      weight_code =imdi
      bottom_level=imdi
      top_level   =imdi
      grid_n      =imdi
      grid_s      =imdi
      grid_e      =imdi
      grid_w      =imdi


!----------------------------------------------------------------------
! Get address for each field from its STASH section/item code
! and STASHmacro tag if a diagnostic, or from its primary pointer
! if prognostic or ancillary field
! Atmosphere -> Ocean (tag=10)
!
      stashmacro_tag=10


      ! Get the surface X-component of windstress (U grid)
      ! The actual field we need depends on whether
      ! coastal tiling is active or not.
      IF (L_CTILE) THEN
         FIELD_CODE = 392
      ELSE
         FIELD_CODE = 219
      END IF

!!$#if defined(ACCESS)
!!$! gol124: auscom coupling
!!$! enforce 219 because 392 is masked
!!$      FIELD_CODE = 219
!!$#endif

! DEPENDS ON: findptr
     Call findptr(atmos_im,3,FIELD_CODE,                                &
     &  process_code,freq_code,start,end,period,                        &
     &  gridpt_code,weight_code,                                        &
     &  bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,             &
     &  stashmacro_tag,imdi,ja_taux,                                    &
#include "argsts.h"
     &  icode,cmessage)


      If (ja_taux==0) Then
        icode=3000 + FIELD_CODE
        cmessage="oasis_inita2o: Coupling field not enabled - taux"
      END IF

      ! Get the surface Y-component of windstress (V grid)
      IF (L_CTILE) THEN
         FIELD_CODE = 394
      ELSE
         FIELD_CODE = 220
      END IF
 
!!$#if defined(ACCESS)
!!$! gol124: auscom coupling
!!$! enforce 220 because 394 is masked
!!$      FIELD_CODE = 220
!!$#endif

! DEPENDS ON: findptr
       Call findptr(                                                    &
     &  atmos_im,3,FIELD_CODE,                                          &
     &  process_code,freq_code,start,end,period,                        &
     &  gridpt_code,weight_code,                                        &
     &  bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,             &
     &  stashmacro_tag,imdi,ja_tauy,                                    &
#include "argsts.h"
     &  icode,cmessage)


      If (ja_tauy==0) Then
        icode=3000 + FIELD_CODE
        cmessage="oasis_inita2o: Coupling field not enabled - tauy"
      END IF


! Net integrated downward solar on atmos grid
! gol124: auscom_coupling
! this is mean over sea
! replace 203 by 201 to get unmasked field
! WARNING: 201 is accounting for albedo
#if defined(ACCESS_MASK1)
! DEPENDS ON: findptr
      Call findptr(                                                     &
     &  atmos_im,1,203,                                                 &
     &  process_code,freq_code,start,end,period,                        &
     &  gridpt_code,weight_code,                                        &
     &  bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,             &
     &  stashmacro_tag,imdi,ja_solar,                                   &
#include "argsts.h"
     &  icode,cmessage)
#else
! DEPENDS ON: findptr
      Call findptr(                                                     &
     &  atmos_im,1,201,                                                 &
     &  process_code,freq_code,start,end,period,                        &
     &  gridpt_code,weight_code,                                        &
     &  bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,             &
     &  stashmacro_tag,imdi,ja_solar,                                   &
#include "argsts.h"
     &  icode,cmessage)
#endif

      If (ja_solar==0) Then
! gol124: auscom_coupling
#if defined(ACCESS_MASK1)
        icode=1203
#else
        icode=1201
#endif
        cmessage="oasis_inita2o: Cpling field not enabled - Net solar"
      END IF

      ! Actual blue-band field used depends on coastal tiling switch.
      IF (L_CTILE) THEN
         FIELD_CODE = 260
      ELSE
         FIELD_CODE = 204
      END IF

      If (icode==0) Then
! Net downward blueband solar on atmos grid
! DEPENDS ON: findptr
          Call findptr(                                                 &
     &      atmos_im,1,FIELD_CODE,                                      &
     &      process_code,freq_code,start,end,period,                    &
     &      gridpt_code,weight_code,                                    &
     &      bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,         &
     &      stashmacro_tag,imdi,ja_blue,                                &
#include "argsts.h"
     &      icode,cmessage)


        If (ja_blue==0) Then
          icode=1000 + FIELD_CODE
          cmessage=                                                     &
     &   "oasis_inita2o: Coupling field not enabled - Blue solar atm"
        END IF
      END IF

      If (icode==0) Then
! Surface evaporation over sea weighted by fractional leads
! DEPENDS ON: findptr
        Call findptr(                                                   &
     &    atmos_im,3,232,                                               &
     &    process_code,freq_code,start,end,period,                      &
     &    gridpt_code,weight_code,                                      &
     &    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,           &
     &    stashmacro_tag,imdi,ja_evap,                                  &
#include "argsts.h"
     &    icode,cmessage)

        If (ja_evap==0) Then
          icode=3232
          cmessage=                                                     &
     &    "oasis_inita2o: Coupling field not enabled - Evap over sea"
        END IF
      END IF

      If (icode==0) Then
! Net downward longwave on atmos grid
! gol124: auscom_coupling
! this is mean over sea
! replace 203 by 201
#if defined(ACCESS_MASK1)
! DEPENDS ON: findptr
        Call findptr(                                                   &
     &    atmos_im,2,203,                                               &
     &    process_code,freq_code,start,end,period,                      &
     &    gridpt_code,weight_code,                                      &
     &    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,           &
     &    stashmacro_tag,imdi,ja_longwave,                              &
#include "argsts.h"
     &    icode,cmessage)
#else
! DEPENDS ON: findptr
        Call findptr(                                                   &
     &    atmos_im,2,201,                                               &
     &    process_code,freq_code,start,end,period,                      &
     &    gridpt_code,weight_code,                                      &
     &    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,           &
     &    stashmacro_tag,imdi,ja_longwave,                              &
#include "argsts.h"
     &    icode,cmessage)
#endif

        If (ja_longwave==0) Then
! gol124: auscom_coupling
#if defined(ACCESS_MASK1)
          icode=2203
#else
          icode=2201
#endif
        cmessage="oasis_inita2o: Coupling field not enabled - Longwave"
        END IF
      END IF

      If (icode==0) then
! Sensible heat on atmos grid, area mean over open sea
! gol124: auscom_coupling
! this is mean over sea
! replace 228 by 217 to get unmasked field
#if defined(ACCESS_MASK1)
! DEPENDS ON: findptr
        Call findptr(                                                   &
     &    atmos_im,3,228,                                               &
     &    process_code,freq_code,start,end,period,                      &
     &    gridpt_code,weight_code,                                      &
     &    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,           &
     &    stashmacro_tag,imdi,ja_sensible,                              &
#include "argsts.h"
     &    icode,cmessage)
#else
! DEPENDS ON: findptr
        Call findptr(                                                   &
     &    atmos_im,3,217,                                               &
     &    process_code,freq_code,start,end,period,                      &
     &    gridpt_code,weight_code,                                      &
     &    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,           &
     &    stashmacro_tag,imdi,ja_sensible,                              &
#include "argsts.h"
     &    icode,cmessage)
#endif

        If (ja_sensible==0) Then
! gol124: auscom_coupling
#if defined(ACCESS_MASK1)
          icode=3228
#else
          icode=3217
#endif
          cmessage=                                                     &
     &   "oasis_inita2o: Coupling field not enabled - Sensible heat"
        END IF
      END IF

      IF (ICODE==0) THEN

! Large-scale snowfall rate on atmos grid
! DEPENDS ON: findptr
      CALL FINDPTR(ATMOS_IM, 4,204,                                     &
     &             PROCESS_CODE,FREQ_CODE,START,END,PERIOD,             &
     &             GRIDPT_CODE,WEIGHT_CODE,                             &
     &             BOTTOM_LEVEL,TOP_LEVEL,GRID_N,GRID_S,GRID_W,GRID_E,  &
     &             STASHMACRO_TAG,IMDI,JA_LSSNOW,                       &
#include "argsts.h"
     &                       ICODE,CMESSAGE)
      IF (JA_LSSNOW == 0) THEN
        ICODE=4204
        CMESSAGE="oasis_inita2o: Coupling field not enabled - LS Snow"
      END IF
      END IF

      IF (ICODE==0) THEN

! Convective snowfall rate on atmos grid
! DEPENDS ON: findptr
      CALL FINDPTR(ATMOS_IM, 5,206,                                     &
     &             PROCESS_CODE,FREQ_CODE,START,END,PERIOD,             &
     &             GRIDPT_CODE,WEIGHT_CODE,                             &
     &             BOTTOM_LEVEL,TOP_LEVEL,GRID_N,GRID_S,GRID_W,GRID_E,  &
     &             STASHMACRO_TAG,IMDI,JA_CVSNOW,                       &
#include "argsts.h"
     &                       ICODE,CMESSAGE)
      IF (JA_CVSNOW == 0) THEN
        ICODE=5206
        CMESSAGE="oasis_inita2o: Coupling field not enabled - Conv Snow"
      END IF
      END IF

      IF (ICODE==0) THEN

! Large-scale rainfall rate on atmos grid
! DEPENDS ON: findptr
      CALL FINDPTR(ATMOS_IM, 4,203,                                     &
     &             PROCESS_CODE,FREQ_CODE,START,END,PERIOD,             &
     &             GRIDPT_CODE,WEIGHT_CODE,                             &
     &             BOTTOM_LEVEL,TOP_LEVEL,GRID_N,GRID_S,GRID_W,GRID_E,  &
     &             STASHMACRO_TAG,IMDI,JA_LSRAIN,                       &
#include "argsts.h"
     &                       ICODE,CMESSAGE)
      IF (JA_LSRAIN == 0) THEN
        ICODE=4203
        CMESSAGE="oasis_inita2o: Coupling field not enabled - LS Rain"
      END IF
      END IF


      IF (ICODE==0) THEN

! Convective rainfall rate on atmos grid
! DEPENDS ON: findptr
      CALL FINDPTR(ATMOS_IM, 5,205,                                     &
     &             PROCESS_CODE,FREQ_CODE,START,END,PERIOD,             &
     &             GRIDPT_CODE,WEIGHT_CODE,                             &
     &             BOTTOM_LEVEL,TOP_LEVEL,GRID_N,GRID_S,GRID_W,GRID_E,  &
     &             STASHMACRO_TAG,IMDI,JA_CVRAIN,                       &
#include "argsts.h"
     &                       ICODE,CMESSAGE)
      IF (JA_CVRAIN == 0) THEN
        ICODE=5205
        CMESSAGE="oasis_inita2o: Coupling field not enabled - Conv Rain"
      END IF
      END IF


! River routing (if required for Ocean
      IF (icode==0) THEN
! RIVER OUTFLOW on atmos grid
! DEPENDS ON: findptr
      CALL FINDPTR(ATMOS_IM, 26,004,                                    &
     &             PROCESS_CODE,FREQ_CODE,START,END,PERIOD,             &
     &             GRIDPT_CODE,WEIGHT_CODE,                             &
     &             BOTTOM_LEVEL,TOP_LEVEL,GRID_N,GRID_S,GRID_W,GRID_E,  &
     &             STASHMACRO_TAG,IMDI,JA_RIVEROUT,                     &
#include "argsts.h"
     &                           ICODE,CMESSAGE)

      IF (JA_RIVEROUT == 0) THEN
        ICODE=26004
        CMESSAGE=                                                       &
     &"oasis_inita2o: Coupling field not enabled - River Outflow ATMOS"
      END IF
      END IF

! WME for ocean KT scheme
      IF (icode==0) THEN
! Wind mixing energy on atmos grid
! DEPENDS ON: findptr
      CALL FINDPTR(ATMOS_IM, 3,224,                                     &
     &             PROCESS_CODE,FREQ_CODE,START,END,PERIOD,             &
     &             GRIDPT_CODE,WEIGHT_CODE,                             &
     &             BOTTOM_LEVEL,TOP_LEVEL,GRID_N,GRID_S,GRID_W,GRID_E,  &
     &             STASHMACRO_TAG,IMDI,JA_WINDMIX,                      &
#include "argsts.h"
     &                           ICODE,CMESSAGE)

      IF (JA_WINDMIX == 0) THEN
        ICODE=03224
        CMESSAGE=                                                       &
     &"oasis_inita2o: Coupling field not enabled - WME ATMOS"
      END IF
      END IF

! Sublimation rate
! Note that field 3353 is a rate(kg/m2/s) and 3231 is a total (kg/ts). 
! Sublimation usually needs to be passed to ocean models (eg. NEMO)  
! as rate. So if coastal tiling is not used, the sublimation total 
! field is converted to a rate in oasis_updatecpl.
! (See UM7.0 or earlier swap_a2o for historical perspective).

      IF (icode==0) THEN
! Sublimation
         IF (L_CTILE) THEN
            FIELD_CODE = 353
         ELSE
            FIELD_CODE = 231
         END IF

! DEPENDS ON: findptr
         CALL FINDPTR(ATMOS_IM, 3,FIELD_CODE,                           &
     &             PROCESS_CODE,FREQ_CODE,START,END,PERIOD,             &
     &             GRIDPT_CODE,WEIGHT_CODE,                             &
     &             BOTTOM_LEVEL,TOP_LEVEL,GRID_N,GRID_S,GRID_W,GRID_E,  &
     &             STASHMACRO_TAG,IMDI,JA_SUBLIM,                       &
#include "argsts.h"
     &                       ICODE,CMESSAGE)

         IF (JA_SUBLIM == 0) THEN
            ICODE=3000 + FIELD_CODE
            CMESSAGE="oasis_inita2o: Cpling field not enabled - Sublim"
         END IF
      END IF

      IF (icode==0) THEN
! Multi category BOTMELT
! NOTE: We're assuming that only multi category fields are being used!
! DEPENDS ON: findptr
        CALL FINDPTR(ATMOS_IM, 3,256,                                   &
     &             PROCESS_CODE,FREQ_CODE,START,END,PERIOD,             &
     &             GRIDPT_CODE,WEIGHT_CODE,                             &
     &             BOTTOM_LEVEL,TOP_LEVEL,GRID_N,GRID_S,GRID_W,GRID_E,  &
     &             STASHMACRO_TAG,IMDI,JA_BOTMELTN,                     &
#include "argsts.h"
     &                        ICODE,CMESSAGE)
        IF (JA_BOTMELTN == 0) THEN
          ICODE=3256
          CMESSAGE="oasis_inita2o: Cpling field not enabled - BOTMELTN"
        END IF
      END IF

      IF (icode==0) THEN
! Multi category TOPMELT
! NOTE: We're assuming that only multi category fields are being used!

! DEPENDS ON: findptr
        CALL FINDPTR(ATMOS_IM, 3,257,                                   &
     &             PROCESS_CODE,FREQ_CODE,START,END,PERIOD,             &
     &             GRIDPT_CODE,WEIGHT_CODE,                             &
     &             BOTTOM_LEVEL,TOP_LEVEL,GRID_N,GRID_S,GRID_W,GRID_E,  &
     &             STASHMACRO_TAG,IMDI,JA_TOPMELTN,                     &
#include "argsts.h"
     &                        ICODE,CMESSAGE)
        IF (JA_TOPMELTN == 0) THEN
          ICODE=3257
          CMESSAGE="oasis_inita2o: Cpling field not enabled - TOPMELTN"
        END IF
      END IF

#if defined(ACCESS)
! gol124: auscom coupling - additional variables
!
! (21) solar radiation                                  um_swflx
! (use ja_solar - Net solar) (mean over sea)

! (22) long wave radiation                              um_lwflx
! (use ja_longwave) (mean over sea)

! (23) sensible heat flux                               um_shflx
! (use ja_sensible) (mean over sea)

      IF (icode==0) THEN
! (24) Surface pressure um_press

! DEPENDS ON: findptr
        CALL FINDPTR(ATMOS_IM, 0,409,                                   &
     &             PROCESS_CODE,FREQ_CODE,START,END,PERIOD,             &
     &             GRIDPT_CODE,WEIGHT_CODE,                             &
     &             BOTTOM_LEVEL,TOP_LEVEL,GRID_N,GRID_S,GRID_W,GRID_E,  &
     &             STASHMACRO_TAG,IMDI,JA_PRESS,                        &
#include "argsts.h"
     &                        ICODE,CMESSAGE)
        IF (JA_PRESS == 0) THEN
          ICODE=409
          CMESSAGE="oasis_inita2o: Cpling field not enabled - PRESSURE"
        END IF
      END IF

! (25) CO2
      IF (icode==0) THEN
! DEPENDS ON: findptr
        CALL FINDPTR(ATMOS_IM, 0,252,                                   &
     &             PROCESS_CODE,FREQ_CODE,START,END,PERIOD,             &
     &             GRIDPT_CODE,WEIGHT_CODE,                             &
     &             BOTTOM_LEVEL,TOP_LEVEL,GRID_N,GRID_S,GRID_W,GRID_E,  &
     &             STASHMACRO_TAG,IMDI,jc_co2,                        &
#include "argsts.h"
     &                        ICODE,CMESSAGE)
        IF (jc_co2 == 0) THEN
          ICODE=252
          CMESSAGE="oasis_inita2o: Cpling field not enabled - 3D CO2"
        END IF
      END IF
! (26) Wind speed
      IF (icode==0) THEN
! DEPENDS ON: findptr
        CALL FINDPTR(ATMOS_IM, 3,225,                                   &
     &             PROCESS_CODE,FREQ_CODE,START,END,PERIOD,             &
     &             GRIDPT_CODE,WEIGHT_CODE,                             &
     &             BOTTOM_LEVEL,TOP_LEVEL,GRID_N,GRID_S,GRID_W,GRID_E,  &
     &             STASHMACRO_TAG,IMDI,JC_U10,                        &
#include "argsts.h"
     &                        ICODE,CMESSAGE)
        IF (JC_U10 == 0) THEN
          !ICODE=3209
          ICODE=3225
          CMESSAGE="oasis_inita2o: Cpling field not enabled - Wnd10m U"
        END IF
      END IF
      !v wind on C grid
      IF (icode==0) THEN
! DEPENDS ON: findptr
        CALL FINDPTR(ATMOS_IM, 3,210,                                   &
     &             PROCESS_CODE,FREQ_CODE,START,END,PERIOD,             &
     &             GRIDPT_CODE,WEIGHT_CODE,                             &
     &             BOTTOM_LEVEL,TOP_LEVEL,GRID_N,GRID_S,GRID_W,GRID_E,  &
     &             STASHMACRO_TAG,IMDI,JC_V10,                        &
#include "argsts.h"
     &                        ICODE,CMESSAGE)
        IF (JC_V10 == 0) THEN
          ICODE=3210
          CMESSAGE="oasis_inita2o: Cpling field not enabled - Wnd10m V"
        END IF
      END IF

#endif

      ! Set up dimensions for use with incoming and outgoing 
      ! transient arrays.
      oasis_imt=lasize(1,fld_type_p,halo_type_no_halo)
      oasis_jmt=lasize(2,fld_type_p,halo_type_no_halo)
      oasis_jmt_u=lasize(2,fld_type_u,halo_type_no_halo)
      oasis_jmt_v=lasize(2,fld_type_v,halo_type_no_halo)
      
      write (6,*) "oasis_imt, oasis_jmt, oasis_jmt_u, oasis_jmt_v=",  & 
         oasis_imt, oasis_jmt, oasis_jmt_u, oasis_jmt_v

      ! Allocate space for incoming arrays and associated
      ! fields required in processing them.  
      allocate(ocn_sst(oasis_imt,oasis_jmt))
      allocate(ocn_sst_orig(oasis_imt,oasis_jmt))
      allocate(tstar_local(oasis_imt,oasis_jmt))
      allocate(tstar_ssi(oasis_imt,oasis_jmt))
      allocate(tstar_sice_local(oasis_imt,oasis_jmt))
      allocate(tstar_land_local(oasis_imt,oasis_jmt))
      allocate(ocn_freeze(oasis_imt,oasis_jmt))
      allocate(ocn_freezen(oasis_imt,oasis_jmt,nice))
      allocate(ocn_hice(oasis_imt,oasis_jmt))
      allocate(ocn_hicen(oasis_imt,oasis_jmt,nice))
      allocate(ocn_snowthick(oasis_imt,oasis_jmt))
      allocate(ocn_snowthickn(oasis_imt,oasis_jmt,nice))
      allocate(ocn_u(oasis_imt,oasis_jmt_u))
      allocate(ocn_v(oasis_imt,oasis_jmt_v))   
      allocate(ocn_co2(oasis_imt,oasis_jmt))
      allocate(ocn_co2fx(oasis_imt,oasis_jmt))

#if defined(ACCESS)
! gol124: auscom coupling
      allocate(auscom_salinity(oasis_imt,oasis_jmt))
      auscom_salinity(:,:) = TFS/-0.054
#endif

      ! Allocate space for outgoing arrays 
      allocate(taux(oasis_imt,oasis_jmt_u))
      allocate(tauy(oasis_imt,oasis_jmt_v))
      allocate(solar2d(oasis_imt,oasis_jmt))
      allocate(blue2d(oasis_imt,oasis_jmt))
      allocate(evap2d(oasis_imt,oasis_jmt))
      allocate(sublim(oasis_imt,oasis_jmt))
      allocate(longwave2d(oasis_imt,oasis_jmt))
      allocate(sensible2d(oasis_imt,oasis_jmt))

      allocate(rainls(oasis_imt,oasis_jmt))
      allocate(snowls(oasis_imt,oasis_jmt))
      allocate(rainconv(oasis_imt,oasis_jmt))
      allocate(snowconv(oasis_imt,oasis_jmt))

      allocate(riverout(oasis_imt,oasis_jmt))
      allocate(wme(oasis_imt,oasis_jmt))
      allocate(um_co2(oasis_imt,oasis_jmt))
      allocate(c_u10(oasis_imt,oasis_jmt_u))
      allocate(c_v10(oasis_imt,oasis_jmt_v))
      allocate(um_wnd10(oasis_imt,oasis_jmt_v))

      ! Allocate arrays which actually gets used for putting data.
      allocate(heatflux(oasis_imt,oasis_jmt))
      allocate(latentflux(oasis_imt,oasis_jmt))
      allocate(totalrain(oasis_imt,oasis_jmt))
      allocate(totalsnow(oasis_imt,oasis_jmt))

      allocate(topmeltn(oasis_imt,oasis_jmt,nice))
      allocate(botmeltn(oasis_imt,oasis_jmt,nice))

#if defined(ACCESS)
! gol124: auscom coupling
      allocate(auscom_swflx(oasis_imt,oasis_jmt))
      allocate(auscom_lwflx(oasis_imt,oasis_jmt))
      allocate(auscom_shflx(oasis_imt,oasis_jmt))
      allocate(auscom_pressure(oasis_imt,oasis_jmt))
#endif

#if defined(ACCESS_SOLAR)
      allocate(c_solar(oasis_imt,oasis_jmt))
      allocate(c_blue(oasis_imt,oasis_jmt))
      allocate(c_longwave(oasis_imt,oasis_jmt))
      allocate(c_taux(oasis_imt,oasis_jmt))
      allocate(c_tauy(oasis_imt,oasis_jmt_v))
      allocate(c_windmix(oasis_imt,oasis_jmt))
      allocate(c_sensible(oasis_imt,oasis_jmt))
      allocate(c_sublim(oasis_imt,oasis_jmt))
      allocate(c_evap(oasis_imt,oasis_jmt))
      allocate(c_botmeltn(oasis_imt,oasis_jmt,5))
      allocate(c_topmeltn(oasis_imt,oasis_jmt,5))
      allocate(c_lsrain(oasis_imt,oasis_jmt))
      allocate(c_lssnow(oasis_imt,oasis_jmt))
      allocate(c_cvrain(oasis_imt,oasis_jmt))
      allocate(c_cvsnow(oasis_imt,oasis_jmt))
      allocate(c_riverout(oasis_imt,oasis_jmt))
      allocate(c_press(oasis_imt,oasis_jmt))
#endif

      ! Set up fractional land points
      allocate(fland_loc(oasis_imt,oasis_jmt))

! DEPENDS ON: from_land_points
      CALL FROM_LAND_POINTS(fland_loc,D1(JFRAC_LAND),                   &
     & atmos_landmask_local,                                            &
     & lasize(1,fld_type_p,halo_type_no_halo)*                          &
     & lasize(2,fld_type_p,halo_type_no_halo), nlandpt)

      ! Ensure land fraction is zero on any missing data points. 
      Do j=1,oasis_jmt
         Do i=1,oasis_imt
            If (fland_loc(I,J).eq.RMDI) fland_loc(I,J)=0.0
         Enddo
      Enddo

      End Subroutine oasis_inita2o
#endif
