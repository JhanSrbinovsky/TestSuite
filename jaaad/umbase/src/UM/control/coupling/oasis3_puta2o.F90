#if defined(OASIS3)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis3_puta2o(                                         &
#include "argd1.h"
#include "arg_atm_fields.h"
  PUT_STEP,cmessage)

  USE oasis3_atm_data_mod

#if defined(ACCESS)
  USE auscom_cpl_data_mod
  USE dump_sent
#endif

  IMPLICIT NONE

  !
  ! Description:
  ! Put data from atmosphere to be got by the Nemo or other Ocean model.
  !  
  ! Author: R. Hill
  ! Current Code Owner : R. Hill
  ! 
  !==================================================================

#include "parvars.h"
#include "decomptp.h"
#include "decompdb.h"

#include "atm_lsm.h"
#include "cmaxsize.h"
#include "typsize.h"
#include "typd1.h"
#include "typptra.h"
#include "typcona.h"
#include "typlndm.h"

#include "caoptr.h"
#include "c_lheat.h"
#include "c_mdi.h"
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"
#include "cntlall.h"
#include "cntlatm.h"
#include "ctracera.h"
#include "typ_atm_fields.h"
#include "nstypes.h"
#include "cruntimc.h"
#include "ccarbon.h"

  !     Subroutine arguments
  LOGICAL :: PUT_STEP       ! Proper data is to be put
  CHARACTER*(*) :: cmessage ! OUT - Error return message

  !     Local variables
  ! Some of the following constants really should
  ! come from somewhere central rather than being defined here 
  REAL :: latentHeatOfCond=lc
  REAL :: latentHeatOfFus=lf

  INTEGER :: i              ! Loop counter
  INTEGER :: j              ! Loop counter
  INTEGER :: k              ! Ice category counter
  INTEGER :: oasis_error, ft
  INTEGER :: icode          ! Error return code (=0 is OK)
  INTEGER :: info           ! Return code from MPP
  INTEGER :: gather_pe      ! Processor for gathering
  INTEGER :: s_cols, s_rows ! Send buffer sizes

  ! Perform necessary processing as described, for example,
  ! in 'Operations 1' on  HadGEM3 coupling diagram. (See HadGEM3 
  ! coupling documentation)

  IF (PUT_STEP) THEN
    ! We set up coupling data based on prognostic contents,
    ! Prepare fields for putting to the coupler
    ! by copying to our temporary arrays.


    ! Set up fields from the T grid
    DO j=1,oasis_jmt
      DO i=1,oasis_imt

!it is done in oasis_updatecpl.F90
!        if (l_co2_interactive) then
!          um_co2(i,j) = co2(i,j,1)
!        else
!          um_co2(i,j) = co2_mmr
!        endif
        ! Copy various heat flux components
        solar2d(i,j)=c_solar(i,j)
        blue2d(i,j)=c_blue(i,j)
        evap2d(i,j)=c_evap(i,j)
        sublim(i,j)=c_sublim(i,j)
        longwave2d(i,j)=c_longwave(i,j)
        sensible2d(i,j)=c_sensible(i,j)

        ! PME components
        rainls(i,j) = c_lsrain(i,j)
        snowls(i,j) = c_lssnow(i,j)
        rainconv(i,j) = c_cvrain(i,j)
        snowconv(i,j) = c_cvsnow(i,j)

        ! River runoff
        riverout(i,j) = c_riverout(i,j)

        ! As in old transo2a code, riverout needs scaling 
        ! to take take account of coastal tiling points
        IF (fland_loc(i,j).NE.1.0) THEN
          riverout(i,j)=riverout(i,j)/(1.0-fland_loc(i,j))
        END IF

        ! Wind mixing energy WME
        wme(i,j) = c_windmix(i,j)

#if defined(ACCESS)
! gol124: auscom coupling
! perhaps this copy is redundant
        ! Sw flux
        auscom_swflx(i,j) = c_solar(i,j)
        ! Lw flux
        auscom_lwflx(i,j) = c_longwave(i,j)
        ! Sensible heat flux
        auscom_shflx(i,j) = c_sensible(i,j)
        ! Surface pressure
        auscom_pressure(i,j) = c_press(i,j)
#endif

      END DO
    END DO

    ! Topmelt and Botmelt (category)
    DO k=1,nice
      DO j=1,oasis_jmt
        DO i=1,oasis_imt
          topmeltn(i,j,k) = c_topmeltn(i,j,k)
          botmeltn(i,j,k) = c_botmeltn(i,j,k)
        END DO
      END DO
    END DO

    !       SURFACE HEAT FLUXES; TOTAL RAIN AND SNOW
    !
    DO j=1,oasis_jmt
      DO i=1,oasis_imt
        ! We need to mask out land points
        ! otherwise we end up using MDIs
        IF (fland_loc(I,J) == 1.0) THEN
          latentflux(i,j)=0.0
          heatflux(i,j)=0.0
          totalrain(i,j)=0.0
          totalsnow(i,j)=0.0
#define ACCESS_MASK1 1
#if defined(ACCESS_MASK1)
! gol124: test using masked stash fields plus zero out land areas
          ! Sw flux
          auscom_swflx(i,j) = 0.0
          ! Lw flux
          auscom_lwflx(i,j) = 0.0
          ! Sensible heat flux
          auscom_shflx(i,j) = 0.0
          ! Surface pressure
          auscom_pressure(i,j) = 0.0
#endif
        ELSE
          latentflux(i,j)=-((latentHeatOfFus+                 &
             latentHeatOfCond)*                               &
             sublim(i,j))

          heatflux(i,j)=solar2d(i,j)+longwave2d(i,j)          &
             -(sensible2d(i,j)+latentHeatOfCond*evap2d(i,j))
          totalrain(i,j)=rainls(i,j)+rainconv(i,j)
          totalsnow(i,j)=snowls(i,j)+snowconv(i,j)
        END IF
      END DO
    END DO

    ! Set up fields from the U grid
    DO j=1,oasis_jmt
      DO i=1,oasis_imt
        taux(i,j)=c_taux(i,j)
      END DO
    END DO

    ! Set up fields from the V grid
    DO j=1,oasis_jmt_v
      DO i=1,oasis_imt
        tauy(i,j)=c_tauy(i,j)
      END DO
    END DO
    
    !wind at 10 metre on v point on C-grid
    um_wnd10 = sqrt(c_v10*c_v10+c_u10*c_u10)
    !change unit of co2 for ocean
    um_co2 = um_co2 * CO2CONV_A2O 

!    write(6,*) "um_co2=", um_co2

#if defined(ACCESS)
    ! dump coupling fields to netcdf file BEFORE sending them to OASIS3
    ! so that if something goes wrong at least we have a copy of data
    ! saved for post-mortem inspection
    IF (sdump_enable) THEN
        CALL write_sdata (dump_hflx, fld_type_p, heatflux,             &
                          oasis_imt, oasis_jmt, .false.)
        CALL write_sdata (dump_solflx, fld_type_p, blue2d,             &
                          oasis_imt, oasis_jmt, .false.)
        CALL write_sdata (dump_runoff, fld_type_p, riverout,           &
                          oasis_imt, oasis_jmt, .false. )
        CALL write_sdata (dump_wme, fld_type_p, wme,                   &
                          oasis_imt, oasis_jmt, .false.)
        CALL write_sdata (dump_train, fld_type_p, totalrain,           &
                          oasis_imt, oasis_jmt, .false.)
        CALL write_sdata (dump_tsnow, fld_type_p, totalsnow,           &
                          oasis_imt, oasis_jmt, .false.)
        CALL write_sdata (dump_evap, fld_type_p, evap2d,               &
                          oasis_imt, oasis_jmt, .false.)
        CALL write_sdata (dump_lhflx, fld_type_p, latentflux,          &
                          oasis_imt, oasis_jmt, .false.)
        DO K = 1, nice  ! Note: NICE MUST=5
           CALL write_sdata (dump_top(K), fld_type_p, topmeltn(1,1,k), &
                             oasis_imt, oasis_jmt, .false.)
           CALL write_sdata (dump_bot(K), fld_type_p, botmeltn(1,1,k), &
                             oasis_imt, oasis_jmt, .false.)
        END DO
        CALL write_sdata (dump_taux, fld_type_u, taux,                 &
                          oasis_imt, oasis_jmt, .false.)
        CALL write_sdata (dump_tauy, fld_type_v, tauy,                 &
                          oasis_imt, oasis_jmt, .true.)

        ! auscom coupling fields
        CALL write_sdata (dump_swflx, fld_type_p, auscom_swflx,        &
                          oasis_imt, oasis_jmt, .false.)
        CALL write_sdata (dump_lwflx, fld_type_p, auscom_lwflx,        &
                          oasis_imt, oasis_jmt, .false.)
        CALL write_sdata (dump_shflx, fld_type_p, auscom_shflx,        &
                          oasis_imt, oasis_jmt, .false.)
        CALL write_sdata (dump_press, fld_type_p, auscom_pressure,     &
                          oasis_imt, oasis_jmt, .false.)
        CALL write_sdata (dump_co2, fld_type_p, um_co2,     &
                          oasis_imt, oasis_jmt, .false.)
        CALL write_sdata (dump_wnd10, fld_type_v, um_wnd10,     &
                          oasis_imt, oasis_jmt_v, .true.)

        ! inc record counter, to be called after last write only!
        CALL inc_srec
    END IF
#endif


    ! Currently we have to observe very strict call sequences - i.e.
    ! the sequence of puts and fields referred to must be in the same
    ! order as the sequence of gets in the receiving component.
    ! If it isn't, OASIS3 will just hang without issuing any warnings!
    !
    ! Clearly this is not a flexible long term FLUME option since it needs
    ! each component to "know" what sequence of puts/gets to use.  
    ! In the long term we'll use the non SEQMODE option for better
    ! flexibility and FLUMEification. But for now a definite sequence
    ! helps with debugging while the system is still in relative infancy.
    !===============================================================

    IF (LTIMER) THEN
! DEPENDS ON: timer
       CALL TIMER('PUTO2A_COMM',3)
    END IF

    ft = fld_type_p
    IF (L_COUPLE_MASTER) THEN
      s_cols = glsize(1,ft)
      s_rows = glsize(2,ft)
    ELSE
      s_cols = oasis_imt
      s_rows = oasis_jmt
    END IF

    ! Put total heat flux
    write(6,*) "oasis3_puta2o: Put heatflux"
    var_ind = vind_heatflux
    ! DEPENDS ON: oasis3_put64
    CALL oasis3_put64(heatflux,oasis_imt,oasis_jmt,oasis_error,ft    &
       ,s_cols,s_rows,PUT_STEP                                       &
       ,halo_type_no_halo,gc_all_proc_group,mype)

    write(6,*) "oasis3_puta2o: Put pen_solar"
    var_ind = vind_pen_solar
    ! DEPENDS ON: oasis3_put64
    CALL oasis3_put64(blue2d,oasis_imt,oasis_jmt,oasis_error,ft      &
       ,s_cols,s_rows,PUT_STEP                                       &
       ,halo_type_no_halo,gc_all_proc_group,mype)

    ! Put riverrunoff
    write(6,*) "oasis3_puta2o: Put runoff"
    var_ind = vind_runoff
    ! DEPENDS ON: oasis3_put64
    CALL oasis3_put64(riverout,oasis_imt,oasis_jmt,oasis_error,ft    &
       ,s_cols,s_rows,PUT_STEP                                       &
       ,halo_type_no_halo,gc_all_proc_group,mype)

    ! Put wme
    write(6,*) "oasis3_puta2o: Put wme"
    var_ind = vind_wme
    ! DEPENDS ON: oasis3_put64
    CALL oasis3_put64(wme,oasis_imt,oasis_jmt,oasis_error,ft         &
       ,s_cols,s_rows,PUT_STEP                                       &
       ,halo_type_no_halo,gc_all_proc_group,mype)

    ! Put total rain
    write(6,*) "oasis3_puta2o: Put train"
    var_ind = vind_train
    ! DEPENDS ON: oasis3_put64
    CALL oasis3_put64(totalrain,oasis_imt,oasis_jmt,oasis_error,ft   &
       ,s_cols,s_rows,PUT_STEP                                       &
       ,halo_type_no_halo,gc_all_proc_group,mype)

    ! Put total snow
    write(6,*) "oasis3_puta2o: Put tsnow"
    var_ind = vind_tsnow
    ! DEPENDS ON: oasis3_put64
    CALL oasis3_put64(totalsnow,oasis_imt,oasis_jmt,oasis_error,ft   &
       ,s_cols,s_rows,PUT_STEP                                       &
       ,halo_type_no_halo,gc_all_proc_group,mype)

    ! Put evap
    write(6,*) "oasis3_puta2o: Put evap2d"
    var_ind = vind_evap2d
    ! DEPENDS ON: oasis3_put64
    CALL oasis3_put64(evap2d,oasis_imt,oasis_jmt,oasis_error,ft      &
       ,s_cols,s_rows,PUT_STEP                                       &
       ,halo_type_no_halo,gc_all_proc_group,mype)

    ! Put latent heat flux
    write(6,*) "oasis3_puta2o: Put lhflx"
    var_ind = vind_lhflx
    ! DEPENDS ON: oasis3_put64
    CALL oasis3_put64(latentflux,oasis_imt,oasis_jmt,oasis_error,ft  &
       ,s_cols,s_rows,PUT_STEP                                       &
       ,halo_type_no_halo,gc_all_proc_group,mype)

    ! Need to deal with each ice category separately.

    ! Put total topmeltn
    write(6,*) "oasis3_puta2o: Put topmelt 1-5"
    DO K = 1, nice  ! Note: NICE MUST=5
      var_ind = vind_topmeltn(k)
      ! DEPENDS ON: oasis3_put64
      CALL oasis3_put64(topmeltn(1,1,k),oasis_imt,oasis_jmt         &
         ,oasis_error,ft                                            &
         ,s_cols,s_rows,PUT_STEP                                    &
         ,halo_type_no_halo,gc_all_proc_group,mype)
    END DO

    ! Put total botmeltn
    write(6,*) "oasis3_puta2o: Put botmelt 1-5"
    DO K = 1, nice  ! Note: NICE MUST=5
      var_ind = vind_botmeltn(k)
      ! DEPENDS ON: oasis3_put64
      CALL oasis3_put64(botmeltn(1,1,k),oasis_imt,oasis_jmt         &
         ,oasis_error,ft                                            &
         ,s_cols,s_rows,PUT_STEP                                    &
         ,halo_type_no_halo,gc_all_proc_group,mype)
    END DO

    ft = fld_type_u
    IF (L_COUPLE_MASTER) THEN
      s_cols = glsize(1,ft)
      s_rows = glsize(2,ft)
    ELSE
      s_cols = oasis_imt
      s_rows = oasis_jmt_u
    END IF

    ! Put x windstress
    write(6,*) "oasis3_puta2o: Put taux"
    var_ind = vind_taux
    ! DEPENDS ON: oasis3_put64
    CALL oasis3_put64(taux,oasis_imt,oasis_jmt,oasis_error,ft        &
       ,s_cols,s_rows,PUT_STEP                                       &
       ,halo_type_no_halo,gc_all_proc_group,mype)

    ft = fld_type_v
    IF (L_COUPLE_MASTER) THEN
      s_cols = glsize(1,ft)
      s_rows = glsize(2,ft)
    ELSE
      s_cols = oasis_imt
      s_rows = oasis_jmt_v
    END IF

    ! Put y windstress
    write(6,*) "oasis3_puta2o: Put tauy"
    var_ind = vind_tauy
    ! DEPENDS ON: oasis3_put64
    CALL oasis3_put64(tauy,oasis_imt,oasis_jmt_v,oasis_error,ft      &
       ,s_cols,s_rows,PUT_STEP                                       &
       ,halo_type_no_halo,gc_all_proc_group,mype)

#if defined(ACCESS)
! gol124: auscom coupling
! additional fields are on t grid
    ft = fld_type_p
    IF (L_COUPLE_MASTER) THEN
      s_cols = glsize(1,ft)
      s_rows = glsize(2,ft)
    ELSE
      s_cols = oasis_imt
      s_rows = oasis_jmt
    END IF

    ! Put sw flux
    write(6,*) "oasis3_puta2o: Put swflx"
    var_ind = vind_swflx
    ! DEPENDS ON: oasis3_put64
    CALL oasis3_put64(auscom_swflx,oasis_imt,oasis_jmt               &
       ,oasis_error,ft,s_cols,s_rows,PUT_STEP                        &
       ,halo_type_no_halo,gc_all_proc_group,mype)

    ! Put lw flux
    write(6,*) "oasis3_puta2o: Put lwflx"
    var_ind = vind_lwflx
    ! DEPENDS ON: oasis3_put64
    CALL oasis3_put64(auscom_lwflx,oasis_imt,oasis_jmt               &
       ,oasis_error,ft,s_cols,s_rows,PUT_STEP                        &
       ,halo_type_no_halo,gc_all_proc_group,mype)

    ! Put sensible heat flux
    write(6,*) "oasis3_puta2o: Put shflx"
    var_ind = vind_shflx
    ! DEPENDS ON: oasis3_put64
    CALL oasis3_put64(auscom_shflx,oasis_imt,oasis_jmt               &
       ,oasis_error,ft,s_cols,s_rows,PUT_STEP                        &
       ,halo_type_no_halo,gc_all_proc_group,mype)

    ! Put surface pressure
    write(6,*) "oasis3_puta2o: Put press"
    var_ind = vind_press
    ! DEPENDS ON: oasis3_put64
    CALL oasis3_put64(auscom_pressure,oasis_imt,oasis_jmt            &
       ,oasis_error,ft,s_cols,s_rows,PUT_STEP                        &
       ,halo_type_no_halo,gc_all_proc_group,mype)

    write(6,*) "oasis3_puta2o: Put co2"
    var_ind = vind_co2
    ! DEPENDS ON: oasis3_put64
    CALL oasis3_put64(um_co2,oasis_imt,oasis_jmt            &
       ,oasis_error,ft,s_cols,s_rows,PUT_STEP                        &
       ,halo_type_no_halo,gc_all_proc_group,mype)

    write(6,*) "oasis3_puta2o: Put wnd10"

    ft = fld_type_v
    IF (L_COUPLE_MASTER) THEN
      s_cols = glsize(1,ft)
      s_rows = glsize(2,ft)
    ELSE
      s_cols = oasis_imt
      s_rows = oasis_jmt_v
    END IF
    var_ind = vind_wnd10
    ! DEPENDS ON: oasis3_put64
    CALL oasis3_put64(um_wnd10,oasis_imt,oasis_jmt_v            &
       ,oasis_error,ft,s_cols,s_rows,PUT_STEP                        &
       ,halo_type_no_halo,gc_all_proc_group,mype)

#endif

    IF (LTIMER) THEN
! DEPENDS ON: timer
       CALL TIMER('PUTO2A_COMM',4)
    END IF

  END IF  ! PUT_STEP=true



END SUBROUTINE oasis3_puta2o
#endif
