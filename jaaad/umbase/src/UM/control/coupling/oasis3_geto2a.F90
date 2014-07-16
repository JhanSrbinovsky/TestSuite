#if defined(OASIS3)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis3_geto2a(                                         &
#include "argd1.h"
#include "arg_atm_fields.h"
  &  GET_STEP,cmessage)


  USE mod_prism
  USE oasis3_atm_data_mod

#if defined(ACCESS)
  USE auscom_cpl_data_mod
  USE dump_received
#endif

  IMPLICIT NONE

  !
  ! Description:
  ! Receive incoming coupling data from the ocean (NEMO) to be used
  ! to drive the atmosphere. 
  !
  ! Author: R. Hill
  !
  ! Current Code Owner : R. Hill
  !
  !=======================================================================

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
#include "c_kappai.h"
#include "c_mdi.h"
! cntlall.h required for ltimer
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"
#include "cntlall.h"
#include "cntlatm.h"
#include "c_0_dg_c.h"
#include "ctracera.h"
#include "typ_atm_fields.h"

  !     Subroutine arguments
  LOGICAL :: GET_STEP       ! Whether to transfer incoming
  ! data to UM atmos D1 or just ignore it.
  CHARACTER*(*) :: cmessage ! OUT - Error return message

  !     Local variables
  INTEGER :: i, j, k, n           ! Loop counters
  INTEGER :: r_cols, r_rows ! Receive buffer dimensions

  INTEGER :: ft
  INTEGER :: icode          ! Error return code (=0 is OK)
  INTEGER :: info           ! Dummy arg for GCOM calls
  INTEGER :: oasis_info     ! Return code from prism calls
  INTEGER :: infosst, infou, infov, infofreeze
  INTEGER :: infofreezen, infosnowthickn, infohicen
  INTEGER :: infomean
  INTEGER :: infoco2, infoco2fx

! gol124: auscom coupling
  REAL :: ltfs
  INTEGER :: infosal

! measure time spent in different stages of geto2a:
! PREP - before communication
! COMM - receiving data from oasis
! MEAN_POLAR_ROW - meaning
! POST - post processing received data
  IF (LTIMER) THEN
! DEPENDS ON: timer
      CALL TIMER('GETO2A_PREP',3)
  END IF

#if defined(ACCESS)
  ltfs = access_tfs
#else
  ltfs = TFS
#endif

  ! Initialise potential incoming fields
  ocn_sst(:,:) = 0.0
  ocn_u(:,:) = 0.0
  ocn_v(:,:) = 0.0
  ocn_freeze(:,:) = 0.0
  ocn_freezen(:,:,:) = 0.0
  ocn_hice(:,:) = 0.0
  ocn_hicen(:,:,:) = 0.0
  ocn_snowthick(:,:) = 0.0
  ocn_snowthickn(:,:,:) = 0.0
  ocn_co2(:,:) = 0.0
  ocn_co2fx(:,:) = 0.0

  tstar_local(:,:) = 0.0
  tstar_ssi(:,:) = 0.0
  tstar_sice_local(:,:) = 0.0
  tstar_land_local(:,:) = 0.0

#if defined(ACCESS)
! gol124: auscom coupling
!  auscom_salinity = 0.0
  IF (.not.ocn_sss) THEN
      auscom_salinity(:,:) = (access_tfs-ZeroDegC)/-0.054
  END IF
#endif

  ! This saves us having to do an explicit mask because
  ! OASIS4 will just overlay the sea points
  IF (L_CTILE) THEN
    DO J = 1, oasis_jmt
      DO I = 1, oasis_imt
        ocn_sst(I,J)=TSTAR_SEA(I,J)
        ocn_sst_orig(I,J)=TSTAR_SEA(I,J)
      END DO
    END DO
  ELSE
    DO J = 1, oasis_jmt
      DO I = 1, oasis_imt
        ocn_sst(I,J)=TSTAR(I,J)
        ocn_sst_orig(I,J)=TSTAR(I,J)
      END DO
    END DO
  END IF


  DO J = 1, oasis_jmt
    DO I = 1, oasis_imt
      tstar_sice_local(I,J)=TSTAR_SICE(I,J)
      tstar_land_local(I,J)=TSTAR_LAND(I,J)
    END DO
  END DO

  DO K=1,nice
    DO J = 1, oasis_jmt
      DO I = 1, oasis_imt
        ocn_freezen(I,J,K)=ICE_FRACT_CAT(I,J,K)
        ocn_hicen(I,J,K)=ICE_THICK_CAT(I,J,K)
        ocn_snowthickn(I,J,K)=SNODEP_SEA_CAT(I,J,K)
      END DO
    END DO
  END DO

  IF (LTIMER) THEN
! DEPENDS ON: timer
      CALL TIMER('GETO2A_PREP',4)
  END IF

  IF (GET_STEP) THEN

    IF (LTIMER) THEN
! DEPENDS ON: timer
       CALL TIMER('GETO2A_COMM',3)
    END IF

    ft = fld_type_p

    IF (L_COUPLE_MASTER) THEN
      r_cols = glsize(1,ft)
      r_rows = glsize(2,ft)
    ELSE
      r_cols = oasis_imt
      r_rows = oasis_jmt
    END IF

    write(6,*) "oasis3_geto2a: Get ocn_sst"
    var_ind = vind_ocn_sst
    ! DEPENDS ON: oasis3_get64
    CALL oasis3_get64(ocn_sst,oasis_imt,oasis_jmt,infosst,ft          &
       ,r_cols,r_rows,get_step                          &
       ,halo_type_no_halo,gc_all_proc_group,mype)
#if defined (ACCESS)
    IF (rdump_enable) THEN
        CALL write_rdata (dump_sst,ft,ocn_sst,oasis_imt,oasis_jmt,.false.)
    END IF
#endif

    write(6,*) "oasis3_geto2a: Get ocn_freeze 1-5"
    infomean = 0
    DO K = 1, nice  ! Note: NICE must=5
      var_ind = vind_ocn_freezen(K)

      ! DEPENDS ON: oasis3_get64
      CALL oasis3_get64(ocn_freezen(1,1,k),oasis_imt,oasis_jmt       &
         ,infofreezen,ft                                  &
         ,r_cols,r_rows,get_step              &
         ,halo_type_no_halo,gc_all_proc_group,mype)
#if defined (ACCESS)
      IF (rdump_enable) THEN
          CALL write_rdata (dump_frzn(k), ft, ocn_freezen(1,1,k),    &
                          oasis_imt, oasis_jmt, .false.)
    END IF
#endif

      IF (K >  1.AND.infofreezen /= infomean) THEN
        WRITE(6,*) "INFORFREEZEN inconsistent"
        STOP
      END IF
      infomean=infofreezen
    END DO ! Over K

    write(6,*) "oasis3_geto2a: Get ocn snowthick 1-5"
    infomean = 0
    DO K = 1, nice  ! Note: NICE must=5
      var_ind = vind_ocn_snowthickn(K)

      ! DEPENDS ON: oasis3_get64
      CALL oasis3_get64(ocn_snowthickn(1,1,k),oasis_imt,oasis_jmt    &
         ,infosnowthickn,ft                               &
         ,r_cols,r_rows,get_step              &
         ,halo_type_no_halo,gc_all_proc_group,mype)
#if defined (ACCESS)
      IF (rdump_enable) THEN
          CALL write_rdata (dump_snwtn(k), ft, ocn_snowthickn(1,1,k),&
                            oasis_imt, oasis_jmt, .false.)
      END IF
#endif
      IF (K >  1.AND.infosnowthickn /= infomean) THEN
        WRITE(6,*) "INFOSNOWTHICKN inconsistent"
        STOP
      END IF
      infomean=infosnowthickn
    END DO

    write(6,*) "oasis3_geto2a: Get ocn_hice 1-5"
    infomean = 0
    DO K = 1, nice  ! Note: NICE must=5
      var_ind = vind_ocn_hicen(K)

      ! DEPENDS ON: oasis3_get64
      CALL oasis3_get64(ocn_hicen(1,1,k),oasis_imt,oasis_jmt         &
         ,infohicen,ft                                    &
         ,r_cols,r_rows,get_step              &
         ,halo_type_no_halo,gc_all_proc_group,mype)
#if defined (ACCESS)
      IF (rdump_enable) THEN
          CALL write_rdata (dump_hicn(k), ft, ocn_hicen(1,1,k),      &
                            oasis_imt, oasis_jmt, .false.)
      END IF
#endif
      IF (K >  1.AND.infohicen /= infomean) THEN
        WRITE(6,*) "INFOHICEN inconsistent"
        STOP
      END IF
      infomean=infohicen
    END DO

    ft = fld_type_u
    IF (L_COUPLE_MASTER) THEN
      r_cols = glsize(1,ft)
      r_rows = glsize(2,ft)
    ELSE
      r_cols = oasis_imt
      r_rows = oasis_jmt_u
    END IF

    write(6,*) "oasis3_geto2a: Get ocn_u"
    var_ind = vind_ocn_u

    ! DEPENDS ON: oasis3_get64
    CALL oasis3_get64(ocn_u,oasis_imt,oasis_jmt_u,infou,ft            &
       ,r_cols,r_rows,get_step                          &
       ,halo_type_no_halo,gc_all_proc_group,mype)
#if defined (ACCESS)
    IF (rdump_enable) THEN
        CALL write_rdata (dump_suno, ft, ocn_u,              &
                          oasis_imt, oasis_jmt, .false.)
    END IF
#endif

    ft = fld_type_v
    IF (L_COUPLE_MASTER) THEN
      r_cols = glsize(1,ft)
      r_rows = glsize(2,ft)
    ELSE
      r_cols = oasis_imt
      r_rows = oasis_jmt_v
    END IF

    write(6,*) "oasis3_geto2a: Get ocn_v"
    var_ind = vind_ocn_v

    ! DEPENDS ON: oasis3_get64
    CALL oasis3_get64(ocn_v,oasis_imt,oasis_jmt_v,infov,ft            &
       ,r_cols,r_rows,get_step                          &
       ,halo_type_no_halo,gc_all_proc_group,mype)
#if defined (ACCESS)
    IF (rdump_enable) THEN
        CALL write_rdata (dump_svno, ft, ocn_v,        &
                          oasis_imt, oasis_jmt, .true.)
    END IF
#endif

    ft = fld_type_p

    IF (L_COUPLE_MASTER) THEN
      r_cols = glsize(1,ft)
      r_rows = glsize(2,ft)
    ELSE
      r_cols = oasis_imt
      r_rows = oasis_jmt
    END IF

    write(6,*) "oasis3_geto2a: Get ocn_co2"
    var_ind = vind_ocn_co2
    ! DEPENDS ON: oasis3_get64
    CALL oasis3_get64(ocn_co2,oasis_imt,oasis_jmt,infoco2,ft          &
       ,r_cols,r_rows,get_step                          &
       ,halo_type_no_halo,gc_all_proc_group,mype)
#if defined (ACCESS)
    IF (rdump_enable) THEN
        CALL write_rdata (dump_ocn_co2,ft,ocn_co2,oasis_imt,oasis_jmt,.false.)
    END IF
#endif
    
    !the oasis3_get ocn_co2 interface is here, but received ocn_co2 is not used in atmosphere


    ft = fld_type_p

    IF (L_COUPLE_MASTER) THEN
      r_cols = glsize(1,ft)
      r_rows = glsize(2,ft)
    ELSE
      r_cols = oasis_imt
      r_rows = oasis_jmt
    END IF

    write(6,*) "oasis3_geto2a: Get ocn_co2fx"
    var_ind = vind_ocn_co2fx
    ! DEPENDS ON: oasis3_get64
    CALL oasis3_get64(ocn_co2fx,oasis_imt,oasis_jmt,infoco2fx,ft          &
       ,r_cols,r_rows,get_step                          &
       ,halo_type_no_halo,gc_all_proc_group,mype)
#if defined (ACCESS)
    IF (rdump_enable) THEN
        CALL write_rdata (dump_ocn_co2fx,ft,ocn_co2fx,oasis_imt,oasis_jmt,.false.)
    END IF
#endif
    !ocn_co2flux is assigned to atmosphere
    co2flux=ocn_co2fx

#if defined(ACCESS)
! gol124: auscom coupling
    IF (ocn_sss) THEN
        ft = fld_type_p

        IF (L_COUPLE_MASTER) THEN
            r_cols = glsize(1,ft)
            r_rows = glsize(2,ft)
        ELSE
            r_cols = oasis_imt
            r_rows = oasis_jmt
        END IF

        write(6,*) "oasis3_geto2a: Get ocn_sal"
        var_ind = vind_sal
! DEPENDS ON: oasis3_get64
        CALL oasis3_get64(auscom_salinity,oasis_imt,oasis_jmt,infosal,ft  &
                         ,r_cols,r_rows,get_step                          &
                         ,halo_type_no_halo,gc_all_proc_group,mype)
        IF (rdump_enable) THEN
            CALL write_rdata (dump_sss, ft, auscom_salinity,         &
                              oasis_imt, oasis_jmt, .false.)
        END IF
    ELSE
        infosal = infosst
    END IF
#else
    infosal = infosst
#endif

    IF (LTIMER) THEN
! DEPENDS ON: timer
       CALL TIMER('GETO2A_COMM',4)
    END IF

#if defined (ACCESS)
    IF (rdump_enable) THEN
        ! gol124: inc record counter, to be called after last write only!
        CALL inc_rrec
    END IF
#endif

    ! If this is a genuine get step.
    ! Should really be able to use the PRISM info
    ! flags to control this entirely, but it all adds to
    ! the run time.

    IF (L_COUPLE_MASTER) THEN
      ! Ensure all PEs know the return codes under all
      ! conditions - PE 0 should definitely have them
      CALL GC_IBCAST(805,1,0,nproc,info,infosst)
      CALL GC_IBCAST(810,1,0,nproc,info,infofreezen)
      CALL GC_IBCAST(825,1,0,nproc,info,infosnowthickn)
      CALL GC_IBCAST(830,1,0,nproc,info,infohicen)
      CALL GC_IBCAST(815,1,0,nproc,info,infou)
      CALL GC_IBCAST(820,1,0,nproc,info,infov)
#if defined(ACCESS)
! gol124: auscom coupling
      IF (ocn_sss) THEN
          CALL GC_IBCAST(835,1,0,nproc,info,infosal)
      ELSE
          infosal = infosst
      END IF
#endif
    END IF

    ! Check for successful receiving of fields from all possible
    ! sources.
    IF ((infosst == infofreezen).AND.                                 &
       (infosst == infohicen).AND.                                   &
       (infosst == infosnowthickn).AND.                              &
       (infosst == infosal).AND.                                     &
       (infosst == PRISM_Recvd.OR.                                   &
       infosst == PRISM_FromRest.OR.                                &
       infosst == PRISM_RecvOut.OR.                                 &
       infosst == PRISM_FromRestOut)) THEN

      ! If this PE deals with the North polar row we must
      ! set a uniform temperature in all grid points along that
      ! row to prevent crashes in the atmos dynamics.
      ! That's all jolly well in a 1xn decomposition, but a
      ! nx1 or nxm composition means we have to do some gathering
      ! and scattering.
      ! There shouldnt be any implications for conservation.
      ! Note: MEAN_POLAR_ROW must be called for all T grid point
      ! fields which we receive.

      ! This is a special adjustment to get sensible
      ! numbers in the polar row. One might regard this as a fix
      ! to make up for shortcomings in OASIS3/4 remapping where
      ! it fails to deal with the N pole as a singularity.
      IF (at_extremity(Pnorth)) THEN
        DO I = 1, oasis_imt
          ocn_sst(i,oasis_jmt)=ocn_sst(i,oasis_jmt-1)
          ocn_freezen(i,oasis_jmt,:)=ocn_freezen(i,oasis_jmt-1,:)
          ocn_hicen(i,oasis_jmt,:)=ocn_hicen(i,oasis_jmt-1,:)
          ocn_snowthickn(i,oasis_jmt,:)=                           &
             ocn_snowthickn(i,oasis_jmt-1,:)
        END DO
#if defined(ACCESS)
! gol124: auscom coupling
        IF (ocn_sss) THEN
            DO I = 1, oasis_imt
               auscom_salinity(i,oasis_jmt)=auscom_salinity(i,oasis_jmt-1)
            END DO
        END IF
#endif
      END IF

! measure time spent in mean_polar_row as whole group
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('MEAN_POLAR_ROW',3)
      END IF

      ! DEPENDS ON: mean_polar_row
      CALL MEAN_POLAR_ROW(OCN_SST)
#if defined(ACCESS)
! gol124: auscom coupling
      IF (ocn_sss) THEN
      ! DEPENDS ON: mean_polar_row
          CALL MEAN_POLAR_ROW(auscom_salinity)
      END IF
#endif

      DO K=1,nice
        ! DEPENDS ON: mean_polar_row
        CALL MEAN_POLAR_ROW(OCN_FREEZEN(1,1,K))
        ! DEPENDS ON: mean_polar_row
        CALL MEAN_POLAR_ROW(OCN_HICEN(1,1,K))
        ! DEPENDS ON: mean_polar_row
        CALL MEAN_POLAR_ROW(OCN_SNOWTHICKN(1,1,K))
      END DO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('MEAN_POLAR_ROW',4)
      END IF

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('GETO2A_POST',3)
      END IF

      ! Perform necessary processing as described, for example,
      ! in 'Operations 5' on  HadGEM3 coupling diagram. (See HadGEM3 
      ! coupling documentation)

      DO K=1,nice
        DO J = 1, oasis_jmt
          DO I = 1, oasis_imt

            ! Ensure ice-depth and snow thickness are not negative (which 
            ! they can be if 2nd order conservative coupling is being used).
            ! No need to do this here for ice fraction because any
            ! negative values will get picked up by the minimum 
            ! ice fraction check below. 

            ocn_hicen(I,J,K)=MAX(0.0,ocn_hicen(I,J,K))  
            ocn_snowthickn(I,J,K)=MAX(0.0,ocn_snowthickn(I,J,K)) 

            ! Here we ensure that regardless of what has happened to the ice 
            ! fraction and GBM ice depth during the coupling,
            ! the minimum local ice depth is 1cm (as in CICE). This is done by
            ! squashing the ice to a smaller area whilst conserving volume (so
            ! the GBM ice depth is unchanged).
            !
            ! Note that freezen is unchanged unless the local ice depth (i.e.
            ! hicen/freezen is less than 1.0e-02, in which case freezen is
            ! decreased so that hicen/freezen=1.0e-02

            ocn_freezen(I,J,K)=MIN(ocn_freezen(I,J,K),            &
               (ocn_hicen(I,J,K)/1.0e-02))

            ! Also impose a minimum ice fraction of 2.0e-04

            IF (ocn_freezen(I,J,K).LT.2.0e-04) THEN
              ocn_freezen(I,J,K)=0.0
              ocn_hicen(I,J,K)=0.0
              ocn_snowthickn(I,J,K)=0.0
            END IF

            ocn_freeze(I,J)=ocn_freeze(I,J)+ocn_freezen(I,J,K)
            ocn_hicen(I,J,K)=ocn_hicen(I,J,K)+                    &
               ocn_snowthickn(I,J,K)*(kappai/kappas)
            ocn_hice(I,J)=ocn_hice(I,J)+ocn_hicen(I,J,K)
            ocn_snowthick(I,J)=ocn_snowthick(I,J)+                &
               ocn_snowthickn(I,J,K)

            IF (ocn_freezen(I,J,K) >  0.0) THEN
              ocn_hicen(I,J,K)=ocn_hicen(I,J,K)/                 &
                 ocn_freezen(I,J,K)
              ocn_snowthickn(I,J,K)=ocn_snowthickn(I,J,K)*       &
                 rhosnow/ocn_freezen(I,J,K)

            END IF
          END DO
        END DO
      END DO

      DO J = 1, oasis_jmt
        DO I = 1, oasis_imt
          IF (ocn_freeze(I,J) >  0.0) THEN
            ocn_freeze(I,J)=MIN(ocn_freeze(I,J),1.0)
            ocn_snowthick(I,J)=ocn_snowthick(I,J)*                &
               rhosnow/ocn_freeze(I,J)
            ocn_hice(I,J)=ocn_hice(I,J)/ocn_freeze(I,J)
          END IF
        END DO
      END DO

      ! The incoming SST is expected in K NOT degrees C
      ! Any necessary transforms will be done in the PSMILE
      ! controlled via the SMIOC, namcouple or in the sending
      ! model component.
      DO J = 1, oasis_jmt
        DO I = 1, oasis_imt

          ! Special form of masking for T points - if
          ! value returned from coupling is approx abs zero
          ! we assume this is an unchanged point masked during
          ! the coupling process and we reinstate
          ! our original value.
          IF (ocn_sst(I,J) <  1.0) THEN
            ocn_sst(I,J)=ocn_sst_orig(I,J)
#if defined(ACCESS)
! gol124: auscom coupling
! for salinity enforce the same masking as for sea surface temp
            auscom_salinity(i,j) = 0.0
#endif
          END IF
        END DO
      END DO

      ! Process surface temperature fields here (as used to be done
      ! in SWAPA20 pre UM version 7.0).

      DO J = 1, oasis_jmt
        DO I = 1, oasis_imt
          tstar_ssi(I,J)=ocn_sst(I,J)

          IF(fland_loc(I,J) <  1.0) THEN
            IF(LTLEADS.EQV..FALSE.) THEN
#if defined(ACCESS)
               IF (ocn_sss) THEN
                   LTFS = ZeroDegC - 0.054 * auscom_salinity(I,J)
               END IF
#endif
              IF(ocn_freeze(I,J) >  0.0) ocn_sst(I,J)=LTFS
            END IF
            tstar_sice_local(I,J)=AMIN1(tstar_sice_local(I,J),TM)
            tstar_ssi(I,J)=ocn_freeze(I,J)*tstar_sice_local(I,J)+  &
               (1.0-ocn_freeze(I,J))*ocn_sst(I,J)
          END IF
          tstar_local(I,J)=fland_loc(I,J)*tstar_land_local(I,J)+    &
             (1.0-fland_loc(I,J))*tstar_ssi(I,J)
        END DO
      END DO

      ! Having got our new values we need to move them
      ! to the appropriate arrays where they'll be employed next 
      ! time the appropriate calculation is called.
      DO J = 1, oasis_jmt
        DO I = 1, oasis_imt

          IF (L_CTILE) THEN
            TSTAR_SEA(I,J)=ocn_sst(I,J)
          END IF

          TSTAR(I,J)=tstar_local(I,J)

          TSTAR_SICE(I,J)=tstar_sice_local(I,J)
          ICE_FRACTION(I,J,1)=ocn_freeze(I,J)
          ICE_THICKNESS(I,J,1)=ocn_hice(I,J)
          SNODEP_SEA(I,J)=ocn_snowthick(I,J)
        END DO
      END DO

      DO K = 1,nice
        DO J = 1, oasis_jmt
          DO I = 1, oasis_imt
            ICE_FRACT_CAT(I,J,K)=ocn_freezen(I,J,K)
            ICE_THICK_CAT(I,J,K)=ocn_hicen(I,J,K)
            SNODEP_SEA_CAT(I,J,K)=ocn_snowthickn(I,J,K)
          END DO
        END DO
      END DO

      IF (LTIMER) THEN
! DEPENDS ON: timer
          CALL TIMER('GETO2A_POST',4)
      END IF

    END IF

    ! Update our surface currents (U)
    IF ( infou == PRISM_Recvd        .OR.                             &
       infou == PRISM_FromRest     .OR.                             &
       infou == PRISM_RecvOut      .OR.                             &
       infou == PRISM_FromRestOut ) THEN

      ! Incoming currents expected in m/s

      ! This is a work around to get sensible
      ! numbers in the polar row -  it is NOT necessarily
      ! sensible for all configurations (e.g. LAMs). 
      ! It is a fix to make up for shortcomings in 
      ! remapping 
      IF (at_extremity(Pnorth)) THEN
        DO I = 1, oasis_imt
          ocn_u(i,oasis_jmt_u)=0.0
        END DO
      END IF

      DO J = 1, oasis_jmt_u
        DO I = 1, oasis_imt
          U_SEA(I,J)=ocn_u(I,J)
        END DO
      END DO
    END IF

    ! Update our surface currents (V)
    IF ( infov == PRISM_Recvd        .OR.                             &
       infov == PRISM_FromRest     .OR.                             &
       infov == PRISM_RecvOut      .OR.                             &
       infov == PRISM_FromRestOut ) THEN

      ! Incoming currents expected in m/s
      DO J = 1, oasis_jmt_v
        DO I = 1, oasis_imt
          V_SEA(I,J)=ocn_v(I,J)
        END DO
      END DO
    END IF

  END IF ! If GET_STEP = true

END SUBROUTINE oasis3_geto2a
#endif
