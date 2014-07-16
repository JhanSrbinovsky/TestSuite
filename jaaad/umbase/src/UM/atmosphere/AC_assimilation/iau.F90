#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Main routine for IAU scheme.

      SUBROUTINE IAU (                                                  &
#include "argppx.h"
#include "argcona.h"
#include "arglndm.h"
     &                 Weight, Timestep,                                &
                                                         ! in
     &                 L_LastCall,                                      &
                                                         ! in
     &                 l_mixing_ratio,                                  &
                                                         ! in
! IAU work arrays:
     &                 IAU_lookup,                                      &
                                                         ! inout
     &                 D1_IAU_k4,                                       &
                                                         ! inout
     &                 D1_IAU,                                          &
                                                         ! inout
! Model fields:
     &                 u,     v,      w,                                &
                                                         ! inout
     &                 u_adv, v_adv,  w_adv,                            &
                                                         ! inout
     &                 theta, exner,  rho,                              &
                                                         ! inout
     &                 q,     qCL,    qCF,                              &
                                                         ! inout
     &                 murk,  TSoil,  TStar, TStar_tile,                &
                                                         ! inout
     &                 p,     p_star, p_theta_levels,                   &
                                                         ! inout
     &                 exner_theta_levels,                              &
                                                         ! inout
     &                 snow_depth,                                      &
                                                         ! inout
     &                 area_cloud_fraction,                             &
                                                         ! inout
     &                 bulk_cloud_fraction,                             &
                                                         ! inout
     &                 cloud_fraction_liquid,                           &
                                                         ! inout
     &                 cloud_fraction_frozen,                           &
                                                         ! inout
     &                 ozone_tracer)                    ! inout 
     

      Use level_heights_mod
      Use trignometric_mod, Only : sec_v_latitude, tan_v_latitude,      &
     &                             sec_theta_latitude
      Use dyn_coriolis_mod, Only : f3_at_v
! Description:
!
!   Main routine for Incremental Analysis Update (IAU) scheme.
!
!   exner, theta, and rho may each be updated using one of two methods:
!
!     (a) direct update using equivalent fields in the IAU file, or
!     (b) indirectly via other fields in the IAU file using the
!         following relationships between atmospheric quantities:
!
!           exner - expression in terms of p
!           theta - hydrostatic equation
!           rho   - equation of state
!
!         The latter two are approximate relationships and are
!         therefore applied in incremental form.
!
!   The following fields may be also be updated directly using
!   equivalent fields in the IAU file:
!
!     u, v, w, q, qCL, qCF, murk (aerosol), ozone
!
!   Alternatively, the following fields may be updated indirectly
!   if the IAU file contains increments to total water, or to q and
!   this option is specifically switched on, and updates to p and
!   theta are available:
!
!   q, qCL, Cl
!
!   Increments to the wind fields u, v and w are copied to their
!   advected counterparts u_adv, v_adv and w_adv.
!
!   If required, the level-one temperature increments are also added to
!   the surface temperature fields (TStar and TStar_tile) at land
!   points, and also to the top-level soil temperature (TSoil(1)) field.
!
!   Unless there is just one update time, the IAU increment fields are
!   by default read in once and held in memory indefinitely. Since this
!   requires a significant amount of memory, there is also a low memory
!   option in which the IAU increments are read in one at a time as
!   needed, though obviously this adds a significant IO overhead.
!
! Method:
!
!   On each update:
!
!     1. Loop over fields in the IAU increment file:
!        a) If field is to be used, and the field is not already in
!           memory, read in field from the IAU increment file.
!        b) Add weighted IAU increment field to the relevant field in
!           the main D1 array.
!
!     2. Do indirect field updates as required.
!
!   On last call, reset stratospheric humidities if required and then
!   close IAU increment file.
!
! Current Code Owner: Adam Clayton.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! Declarations:

      IMPLICIT NONE

! Common blocks:

#include "varcld.h"
#include "cppxref.h"
#include "csubmodl.h"
#include "ppxlook.h"
#include "parvars.h"
#include "typsize.h"
#include "nstypes.h"
#include "cmaxsize.h"
#include "typcona.h"
#include "typlndm.h"
#include "clookadd.h"
#include "cprintst.h"
#include "ctfilt.h"
#include "cntlatm.h"
#include "cruntimc.h"
#include "c_g.h"
#include "c_lheat.h"
#include "c_r_cp.h"
#include "c_epslon.h"
#include "c_pi.h"
#include "c_visbty.h"
#include "c_kinds.h"

! Subroutine arguments:

      REAL,         INTENT(IN)    ::                                    &

     &  Weight                                                          &
                                     ! IAU weight.
     &, Timestep                     ! Model timestep (s)

      LOGICAL,      INTENT(IN)    ::                                    &

     &  L_LastCall                                                      &
                                     ! .TRUE. for last call.
     &, l_mixing_ratio               ! Use mixing ratio code


      INTEGER,      INTENT(INOUT) ::                                    &
                                     ! Lookup tables of IAU inc file

     &  IAU_lookup(IAU_Len1Lookup, IAU_Len2Lookup)

      REAL(KIND=real32), INTENT(INOUT) ::                               &
                                          ! Array for packed IAU fields

     &  D1_IAU_k4(D1_IAU_k4_len)

      REAL,         INTENT(INOUT) ::                                    &
                                     ! Array for unpacked IAU fields

     &  D1_IAU(D1_TFilt_len)

      REAL,         INTENT(INOUT) ::                                    &
                                     ! Model fields.

     &  u      ( 1 - offx   : row_length + offx,                        &
     &           1 - offy   : rows       + offy,                        &
     &           1          : model_levels ),                           &

     &  v      ( 1 - offx   : row_length + offx,                        &
     &           1 - offy   : n_rows     + offy,                        &
     &           1          : model_levels ),                           &

     &  w      ( 1 - offx   : row_length + offx,                        &
     &           1 - offy   : rows       + offy,                        &
     &           0          : model_levels ),                           &

     &  u_adv  ( 1 - halo_i : row_length + halo_i,                      &
     &           1 - halo_j : rows       + halo_j,                      &
     &           1          : model_levels ),                           &

     &  v_adv  ( 1 - halo_i : row_length + halo_i,                      &
     &           1 - halo_j : n_rows     + halo_j,                      &
     &           1          : model_levels ),                           &

     &  w_adv  ( 1 - halo_i : row_length + halo_i,                      &
     &           1 - halo_j : rows       + halo_j,                      &
     &           0          : model_levels ),                           &

     &  theta  ( 1 - offx   : row_length + offx,                        &
     &           1 - offy   : rows       + offy,                        &
     &           1          : model_levels ),                           &

     &  exner  ( 1 - offx   : row_length + offx,                        &
                                                  ! Exner on
     &           1 - offy   : rows       + offy,                        &
                                                  ! rho levels.
     &           1          : model_levels + 1),                        &

     &  rho    ( 1 - offx   : row_length + offx,                        &
     &           1 - offy   : rows       + offy,                        &
     &           1          : model_levels ),                           &

     &  q      ( 1 - halo_i : row_length + halo_i,                      &
     &           1 - halo_j : rows       + halo_j,                      &
     &           1          : wet_levels ),                             &

     &  qCL    ( 1 - halo_i : row_length + halo_i,                      &
     &           1 - halo_j : rows       + halo_j,                      &
     &           1          : wet_levels ),                             &

     &  qCF    ( 1 - halo_i : row_length + halo_i,                      &
     &           1 - halo_j : rows       + halo_j,                      &
     &           1          : wet_levels ),                             &

     &  murk   ( 1 - offx   : row_length + offx,                        &
     &           1 - offy   : rows       + offy,                        &
     &           1          : model_levels ),                           &

     &  TSoil      ( land_field ),                                      &
                                   ! Only need top level.

     &  TStar      ( theta_field_size ),                                &

     &  TStar_tile ( land_field, ntiles ),                              &

     &  p      ( 1 - offx   : row_length + offx,                        &
                                                  ! Pressure on
     &           1 - offy   : rows       + offy,                        &
                                                  ! rho levels.
     &           1          : model_levels + 1),                        &

     &  p_star ( row_length, rows ),                                    &

     &  p_theta_levels     ( 1 - offx : row_length + offx,              &
     &                       1 - offy : rows       + offy,              &
     &                       1        : model_levels ),                 &

     &  exner_theta_levels ( 1 - offx : row_length + offx,              &
     &                       1 - offy : rows       + offy,              &
     &                       1        : model_levels ),                 &

     &  snow_depth         ( theta_field_size ),                        &

     &  area_cloud_fraction   ( row_length, rows, wet_levels ),         &

     &  bulk_cloud_fraction ( 1 - halo_i : row_length + halo_i,         &
     &                        1 - halo_j : rows       + halo_j,         &
     &                        1          : wet_levels ),                &

     &  cloud_fraction_liquid ( 1 - halo_i : row_length + halo_i,       &
     &                          1 - halo_j : rows       + halo_j,       &
     &                          1          : wet_levels ),              &

     &  cloud_fraction_frozen ( 1 - halo_i : row_length + halo_i,       &
     &                          1 - halo_j : rows       + halo_j,       &
     &                          1          : wet_levels ),              &

     &  ozone_tracer (1- offx   : row_length + offx,                   &
     &                          1 - offy   : rows       + offy,         &
     &                          1          : model_levels )

! Local variables:

      INTEGER :: i, j, k,                                               &
     &           tile_num,                                              &
     &           field_size,                                            &
     &           FieldNum,                                              &
     &           D1_IAU_addr,                                           &
     &           s_addr,                                                &
     &           e_addr,                                                &
     &           Addr,                                                  &
     &           Code,                                                  &
     &           LocFldLen,                                             &
     &           ICode

      REAL    :: delta_r,                                               &
     &           thetaV_old_int,                                        &
     &           thetaV_new_int,                                        &
     &           Weight_1,                                              &
     &           Weight_2,                                              &
     &           Term_1,                                                &
     &           Term_2,                                                &
     &           kappa_recip,                                           &
     &           p_zero_recip,                                          &
     &           exner_exp

      ! These variables are for PC2 calculations
      REAL, DIMENSION(:,:,:), ALLOCATABLE::                             &
     &           q_work                                                 &
                                 ! Vapour work array (kg kg-1)
     &,          qcl_work                                               &
                                 ! Liquid work array (kg kg-1)
     &,          qcf_work                                               &
                                 ! Ice work array    (kg kg-1)
     &,          t_work                                                 &
                                 ! Temperature work array  (K)
     &,          p_work                                                 &
                                 ! Pressure (on theta levels)
!                                !             work array (Pa)
     &,          delta_q                                                &
                                 ! Change in q to force condensation
     &,          delta_qcl                                              &
                                 ! Change in qcl to force condensation
     &,          delta_qcf                                              &
                                 ! Change in qcf
     &,          delta_t                                                &
                                 ! Change in t to force condensation
     &,          delta_p         ! Change in p to force condensation

      REAL    :: P_UPPER_LIMIT,                                         &
     &           THETA_UPPER_LIMIT


      PARAMETER ( kappa_recip  = 1.0/kappa       )
      PARAMETER ( p_zero_recip = 1.0/p_zero      )
      PARAMETER ( exner_exp    = 1.0/kappa - 1.0 )

      ! Array for holding a single increment field:
      REAL :: IAUIncFld(theta_field_size)

      ! These large work arrays will not necessarily be needed, so
      ! allocate them on demand:
      REAL, ALLOCATABLE :: exner_old                  (:,:,:)
      REAL, ALLOCATABLE :: thetaV_old                 (:,:,:)
      REAL, ALLOCATABLE :: thetaV_new                 (:,:,:)
      REAL, ALLOCATABLE :: pv_at_theta                (:,:,:)
      REAL, ALLOCATABLE :: AreaDivBulk                (:,:,:)
      REAL, ALLOCATABLE :: cloud_fraction_liquid_0    (:)
      REAL, ALLOCATABLE :: cloud_fraction_liquid_plus (:)
      REAL, ALLOCATABLE :: cloud_fraction_frozen_0    (:)
      REAL, ALLOCATABLE :: cloud_fraction_frozen_plus (:)
      REAL, ALLOCATABLE :: Cl_Inc                     (:,:,:)
      REAL, ALLOCATABLE :: Cf_Inc                     (:,:,:)
      REAL, ALLOCATABLE :: Cb_Inc                     (:,:,:)
      REAL, ALLOCATABLE :: initialScaling             (:,:,:)
      REAL, ALLOCATABLE :: scalingCoeff               (:,:,:)
      REAL, ALLOCATABLE :: p_theta_levels_old         (:,:,:)
      REAL, ALLOCATABLE :: qcl_Inc                    (:,:,:)
      REAL, ALLOCATABLE :: qCL_0                      (:)
      REAL, ALLOCATABLE :: qCL_plus                   (:)
      REAL, ALLOCATABLE :: qCF_0                      (:)
      REAL, ALLOCATABLE :: qCF_max                    (:)
      REAL, ALLOCATABLE :: qCF_plus                   (:)
      REAL, ALLOCATABLE :: qSatW                      (:)
      REAL, ALLOCATABLE :: qT                         (:,:,:)
      REAL, ALLOCATABLE :: qT_plus                    (:,:,:)
      REAL, ALLOCATABLE :: VarTemp                    (:)
      REAL, ALLOCATABLE :: Var_qT                     (:)
      REAL, ALLOCATABLE :: Var_BGT                    (:)

      REAL :: t1_old    ( row_length, rows ),                           &
     &        t1_inc    ( theta_field_size ),                           &
     &        p_tmp     ( row_length, rows ),                           &
     &        q_sat_wat ( row_length, rows )

      CHARACTER(80) :: CMessage

      CHARACTER(*)  :: RoutineName
      PARAMETER      ( RoutineName='IAU' )

      LOGICAL, SAVE :: L_FirstUpdate = .TRUE.

! PC2 options are only important if L_pc2 = .true.
! Seek advice from the PC2 team before altering these parameters
! from .true., you may need to have put in place large amounts
! of extra code first.
      LOGICAL,PARAMETER:: L_pc2_cond=.true.  ! Do condensation and
! liquid cloud fraction changes.
      LOGICAL,PARAMETER:: L_pc2_cfl =.true.  ! Do liquid cloud
! fraction changes. This requires that the condensation as a
! result of assimilation is calculated either directly or via the
! estimate selected by setting l_pc2_cond to true.
! Note: One must also not run with condensation changes on but the
! liquid cloud fraction changes off.
      LOGICAL,PARAMETER:: L_pc2_cff =.true.  ! Do ice cloud fraction
! changes. This requires that the ice increment from assimilation is
! calculated directly.

! Artificially imposed upper limits on theta increments for pressures
! less than the P_UPPER_LIMIT. Above this level the maximum theta increment 
! is determined by THETA_UPPER_LIMIT. This restriction is only applied where
! L_IAU_UPPER_THETA is set to true.

      PARAMETER (P_UPPER_LIMIT = 200.0)
      PARAMETER (THETA_UPPER_LIMIT = 100.0)
! External functions called:

! External subroutines called:

!- End of header ------------------------------------------------------
      
!----------------------------------------------------------------------
! [1]: Add weighted increments to fields in the main D1 array.
!----------------------------------------------------------------------

      IF (ABS(Weight) < Weight_LL) THEN

        IF (PrintStatus >= PrStatus_Normal) THEN
          WRITE (6,*) ''
          WRITE (6,'(A,F13.10,A)') ' IAU: Weight (', Weight,            &
     &                             ') too small to bother with.'
          WRITE (6,*) ''
        END IF

      ELSE

!----------------------------------------------------------------------
! [1.1]: If necessary, make sure p is consistent with exner.
!----------------------------------------------------------------------

        IF (L_IAU_CalcExnerIncs) THEN ! exner will be calculated from p

          DO k = 1, model_levels+1
            DO j = 1, rows
              DO i = 1, row_length

                p(i,j,k) = exner(i,j,k)**kappa_recip * p_zero

              END DO
            END DO
          END DO

        END IF

!----------------------------------------------------------------------
! [1.2]: If required, save fields required for calculating increments
!        to theta, rho, level 1 temperature, q, qCl & Cl.
!----------------------------------------------------------------------

        ! Exner and thetaV:
        IF (L_IAU_CalcThetaIncs .OR. L_IAU_CalcRhoIncs) THEN

          ALLOCATE ( exner_old (row_length,rows,model_levels+1) )
          ALLOCATE ( thetaV_old(row_length,rows,model_levels)   )

          DO k = 1, model_levels+1
            DO j = 1, rows
              DO i = 1, row_length

                exner_old(i,j,k) = exner(i,j,k)

              END DO
            END DO
          END DO

          DO k = 1, wet_levels
            DO j = 1, rows
              DO i = 1, row_length

                thetaV_old(i,j,k) = theta(i,j,k)                        &
     &                            * ( 1.0 + C_Virtual * q(i,j,k) )

              END DO
            END DO
          END DO

          DO k = wet_levels+1, model_levels
            DO j = 1, rows
              DO i = 1, row_length

                thetaV_old(i,j,k) = theta(i,j,k)

              END DO
            END DO
          END DO

        END IF

        ! Level 1 temperature:
        IF (L_IAU_IncTStar) THEN

          DO j = 1, rows
            DO i = 1, row_length

              t1_old(i,j) = exner_theta_levels(i,j,1)                   &
     &                    * theta             (i,j,1)

            END DO
          END DO

        END IF
        ! q, qCL, qCF, Cl& Cf:

        ! Loop over fields in the IAU increment file:
        IF(.NOT. L_IAU_CalcCloudIncs) THEN

          DO FieldNum = 1, IAU_Len2Lookup

            Code = IAU_lookup(ITEM_CODE, FieldNum)
            ! If qT increments are present:
            IF (Code == 18001) THEN

              L_IAU_CalcCloudIncs = .TRUE.
              EXIT

            END IF

          END DO ! FieldNum

        END IF

        IF(L_IAU_CalcCloudIncs) THEN

          field_size = row_length * rows * wet_levels

          ALLOCATE ( AreaDivBulk        (row_length, rows, wet_levels) )
          ALLOCATE ( p_theta_levels_old (row_length, rows, wet_levels) )
          ALLOCATE ( qT                 (row_length, rows, wet_levels) )
          ALLOCATE ( qT_plus            (row_length, rows, wet_levels) )
          ALLOCATE ( cloud_fraction_liquid_0    (field_size          ) )
          ALLOCATE ( cloud_fraction_liquid_plus (field_size          ) )
          ALLOCATE ( qCL_0                      (field_size          ) )
          ALLOCATE ( qCL_plus                   (field_size          ) )
          ALLOCATE ( qSatW                      (field_size          ) )
          ALLOCATE ( Var_qT                     (field_size          ) )
          ALLOCATE ( VarTemp                    (field_size          ) )

          p_theta_levels_old = p_theta_levels                           &
     &                              (1:row_length, 1:rows, 1:wet_levels)

          qT                 = q    (1:row_length, 1:rows, 1:wet_levels)&
     &                         +                                        &
     &                         qCL  (1:row_length, 1:rows, 1:wet_levels)

          VarTemp            = RESHAPE(exner_theta_levels               &
     &                              (1:row_length, 1:rows, 1:wet_levels)&
     &                         *                                        &
     &                         theta(1:row_length, 1:rows, 1:wet_levels)&
     &                         -                                        &
     &                         qCL  (1:row_length, 1:rows, 1:wet_levels)&
     &                         *Lc/Cp, (/field_size/))
          IF (L_IAU_IncrementIce) THEN

            ALLOCATE ( cloud_fraction_frozen_0    (field_size        ) )
            ALLOCATE ( cloud_fraction_frozen_plus (field_size        ) )
            ALLOCATE ( qCF_0                      (field_size        ) )
            ALLOCATE ( qCF_max                    (field_size        ) )
            ALLOCATE ( qCF_plus                   (field_size        ) )
            ALLOCATE ( Var_BGT                    (field_size        ) )
            Var_BGT            = RESHAPE(exner_theta_levels             &
     &                                (1:row_length,1:rows,1:wet_levels)&
     &                           *                                      &
     &                           theta(1:row_length,1:rows,1:wet_levels)&
     &                           , (/field_size/))

          END IF

! DEPENDS ON: qsat_wat
          CALL QSAT_WAT (qSatW,                                         &
     &        VarTemp,                                                  &
     &        RESHAPE(p_theta_levels                                    &
     &                (1:row_length, 1:rows, 1:wet_levels),             &
     &                (/field_size/)), field_size)

        ELSE

          ! Report error if qcf increments requested but either
          ! L_IAU_CalcCloudIncs is set false or increment file does not
          ! contain qT'
          IF (L_IAU_IncrementIce) THEN
            ICode    = 1
            CMessage =                                                  &
     &        'qcf incrementing: set L_IAU_CalcCloudIncs or provide qT'
! DEPENDS ON: ereport
            CALL EReport (RoutineName, ICode, CMessage)

          END IF

        END IF


! Input values needed for the PC2 scheme
        If (L_pc2 .and.                                                 &
     &     (L_pc2_cond .or. L_pc2_cfl .or. L_pc2_cff)   ) then

          Allocate ( q_work  (row_length,rows,wet_levels) )
          Allocate ( qcl_work(row_length,rows,wet_levels) )
          Allocate ( qcf_work(row_length,rows,wet_levels) )
          Allocate ( t_work  (row_length,rows,wet_levels) )
          Allocate ( p_work  (row_length,rows,wet_levels) )
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i= 1, row_length
                ! Copy input values into the work arrays
                q_work  (i,j,k) = q  (i,j,k)
                qcl_work(i,j,k) = qcl(i,j,k)
                qcf_work(i,j,k) = qcf(i,j,k)
                p_work  (i,j,k) = p_theta_levels    (i,j,k)
                t_work  (i,j,k) = theta(i,j,k)                          &
     &                          * exner_theta_levels(i,j,k)
              End Do
            End Do
          End Do
        End If  ! L_pc2 .and. (L_pc2_cond.or.L_pc2_cfl.or.L_pc2_cff)

!----------------------------------------------------------------------
! [1.3]: DIRECT field updates using data in the IAU increment file.
!----------------------------------------------------------------------

        D1_IAU_addr = 1 ! Field address in IAU D1 array

        ! Loop over fields in the IAU increment file:
        DO FieldNum = 1, IAU_Len2Lookup

          Code = IAU_lookup(ITEM_CODE, FieldNum)
          k    = IAU_lookup(LBLEV,     FieldNum)

          ! Get local field length:
          LocFldLen = 0
          DO i = 1, IAU_NumFldCodes
            IF (IAU_FldCodes(i) == Code)                                &
     &        LocFldLen = IAU_LocFldLens(i)
          END DO

          IF (LocFldLen == 0) CYCLE ! Unused field

          Addr = 1 ! Address in IAUIncFld

!----------------------------------------------------------------------
! [1.3.1]: Obtain IAU increment field.
!----------------------------------------------------------------------

          IF (L_FirstUpdate .OR. L_IAU_LowMem) THEN

            ! Field not in memory, so read it from IAU increment file:
! DEPENDS ON: readiaufield
            CALL ReadIAUField (                                         &
#include "argppx.h"
#include "argcona.h"
#include "arglndm.h"
     &                          L_FirstUpdate,                          &
                                                         ! in
     &                          IAU_lookup,                             &
                                                         ! in
     &                          FieldNum,                               &
                                                         ! in
     &                          LocFldLen,                              &
                                                         ! in
     &                          IAUIncFld(1:LocFldLen) ) ! out

            IF (.NOT.L_IAU_LowMem) THEN

              s_addr = D1_IAU_addr
              e_addr = D1_IAU_addr+LocFldLen-1

              ! Save field in IAU D1 array:
              IF (L_Pack_D1_IAU) THEN

                ! Save in packed form:
! DEPENDS ON: pack21
                CALL pack21 ( LocFldLen,                                &
                                                         ! in
                              IAUIncFld(1:LocFldLen),                   &
                                                         ! in
                              D1_IAU_k4(s_addr:e_addr) ) ! out

              ELSE

                D1_IAU(s_addr:e_addr) = IAUIncFld(1:LocFldLen)

              END IF

            END IF

          ELSE ! Field already in memory, so retrieve it

            s_addr = D1_IAU_addr
            e_addr = D1_IAU_addr+LocFldLen-1

            IF (L_Pack_D1_IAU) THEN

              ! Unpack field:
! DEPENDS ON: expand21
              CALL expand21 ( LocFldLen,                                &
                                                        ! in
     &                        D1_IAU_k4(s_addr:e_addr),                 &
                                                        ! in
     &                        IAUIncFld(1:LocFldLen) )  ! out

            ELSE

              IAUIncFld(1:LocFldLen) = D1_IAU(s_addr:e_addr)

            END IF

          END IF

!----------------------------------------------------------------------
! [1.3.2]: u and u_adv.
!----------------------------------------------------------------------

          IF (Code == 2 .AND. k >= 1 .AND. k <= model_levels) THEN

            DO j = 1, rows
              DO i = 1, row_length

                u    (i,j,k) = u    (i,j,k) + Weight * IAUIncFld(Addr)
                u_adv(i,j,k) = u_adv(i,j,k) + Weight * IAUIncFld(Addr)
                Addr = Addr + 1

              END DO
            END DO

          END IF

!----------------------------------------------------------------------
! [1.3.3]: v and v_adv.
!----------------------------------------------------------------------

          IF (Code == 3 .AND. k >= 1 .AND. k <= model_levels) THEN

            DO j = 1, n_rows
              DO i = 1, row_length

                v    (i,j,k) = v    (i,j,k) + Weight * IAUIncFld(Addr)
                v_adv(i,j,k) = v_adv(i,j,k) + Weight * IAUIncFld(Addr)
                Addr = Addr + 1

              END DO
            END DO

          END IF

!----------------------------------------------------------------------
! [1.3.4]: w and w_adv.
!----------------------------------------------------------------------

          IF ( Code == 150 .AND.                                        &
     &         ( k == 9999 .OR.                                         &
     &           (k >= 1 .AND. k <= model_levels) ) ) THEN

            IF (k == 9999) k = 0

            DO j = 1, rows
              DO i = 1, row_length

                w    (i,j,k) = w    (i,j,k) + Weight * IAUIncFld(Addr)
                w_adv(i,j,k) = w_adv(i,j,k) + Weight * IAUIncFld(Addr)
                Addr = Addr + 1

              END DO
            END DO

          END IF

!----------------------------------------------------------------------
! [1.3.5]: theta.
!----------------------------------------------------------------------
          
          IF (Code == 4 .AND. k >= 1 .AND. k < model_levels) THEN
           
           IF(L_IAU_UPPER_THETA) THEN 
             DO j = 1, rows
               DO i = 1, row_length
                
                ! This logical has been added to impose a limit of 100 on the  
                ! theta increments at levels above 200.

                   IF (abs(p(i,j,k)) < P_UPPER_LIMIT .AND.             &
     &                (abs(Weight * IAUIncFld(Addr)) >                  &
     &                                          THETA_UPPER_LIMIT)) THEN

                      WRITE(6,*) 'IAUIncFld(Addr) is  ', IAUIncFld(Addr),&
     &                           'theta inc restricted for lev ', k,     &
     &                           ' pressure ', p(i,j,k)
                
 
                        IF (Weight * IAUIncFld(Addr) < 0) THEN
                           theta(i,j,k) = theta(i,j,k) - THETA_UPPER_LIMIT
                        ELSE
                           theta(i,j,k) = theta(i,j,k) + THETA_UPPER_LIMIT
                        END IF 
                    ELSE
                      theta(i,j,k) = theta(i,j,k) + Weight * IAUIncFld(Addr)
                    END IF

                 Addr = Addr + 1

               END DO
             END DO

           ELSE
             DO j = 1, rows
               DO i = 1, row_length      
                     theta(i,j,k) = theta(i,j,k) + Weight * IAUIncFld(Addr)

                     Addr = Addr + 1
               END DO
             END DO
      
           END IF ! Check on L_IAU_UPPER_THETA Switch
          
          END IF ! Check on  Code == 4

!----------------------------------------------------------------------
! [1.3.6]: rho.
!----------------------------------------------------------------------

          IF (Code == 253 .AND. k >= 1 .AND. k <= model_levels) THEN

            DO j = 1, rows
              DO i = 1, row_length

                rho(i,j,k) = rho(i,j,k) + Weight * IAUIncFld(Addr)
                Addr = Addr + 1

              END DO
            END DO

          END IF

!----------------------------------------------------------------------
! [1.3.7]: exner.
!----------------------------------------------------------------------

          IF (Code == 255 .AND.                                         &
     &        k >= 1 .AND. k <= model_levels+1) THEN

            DO j = 1, rows
              DO i = 1, row_length

                exner(i,j,k) = exner(i,j,k) + Weight * IAUIncFld(Addr)
                Addr = Addr + 1

              END DO
            END DO

          END IF

!----------------------------------------------------------------------
! [1.3.8]: p.
!----------------------------------------------------------------------

          IF (Code == 407 .AND.                                         &
     &        k >= 1 .AND. k <= model_levels+1) THEN

            DO j = 1, rows
              DO i = 1, row_length

                p(i,j,k) = p(i,j,k) + Weight * IAUIncFld(Addr)
                Addr = Addr + 1

              END DO
            END DO

          END IF

!----------------------------------------------------------------------
! [1.3.9]: q.
!----------------------------------------------------------------------

          IF (Code == 10 .AND. k >= 1 .AND. k <= wet_levels             &
     &                                .AND. k <  model_levels) THEN

            IF (L_IAU_CalcCloudIncs) THEN

              DO j = 1, rows
                DO i = 1, row_length

                  qT_plus(i,j,k) = qT(i,j,k) + Weight * IAUIncFld(Addr)
                  Addr = Addr + 1

                END DO
              END DO

            ELSE

              DO j = 1, rows
                DO i = 1, row_length

                  q(i,j,k) = q(i,j,k) + Weight * IAUIncFld(Addr)
                  ! Apply lower limit to q:
                  IF (q(i,j,k) < q_min) q(i,j,k) = q_min
                  Addr = Addr + 1

                END DO
              END DO

            END IF

          END IF

!----------------------------------------------------------------------
! [1.3.10]: qCL.
!----------------------------------------------------------------------

          IF (Code == 254 .AND. k >= 1 .AND. k <= wet_levels            &
     &             .AND. .NOT. L_IAU_CalcCloudIncs) THEN

            DO j = 1, rows
              DO i = 1, row_length

                qCL(i,j,k) = qCL(i,j,k) + Weight * IAUIncFld(Addr)
                Addr = Addr + 1

                ! Check for non-positive qCLs:
                IF (qCL(i,j,k) <= 0.0) THEN
                  qCL(i,j,k) = 0.0
                  cloud_fraction_liquid(i,j,k) = 0.0
                  IF ( qCF(i,j,k) <= q_CC_tol ) THEN
                    area_cloud_fraction(i,j,k) = 0.0
                    bulk_cloud_fraction(i,j,k) = 0.0
                  END IF
                END IF

              END DO
            END DO

          END IF

!----------------------------------------------------------------------
! [1.3.11]: qCF.
!----------------------------------------------------------------------

          IF (Code == 12 .AND. k >= 1 .AND. k <= wet_levels             &
     &             .AND. .NOT. L_IAU_CalcCloudIncs) THEN

            DO j = 1, rows
              DO i = 1, row_length

                qCF(i,j,k) = qCF(i,j,k) + Weight * IAUIncFld(Addr)
                Addr = Addr + 1

                ! Check for non-positive qCFs:
                IF (qCF(i,j,k) <= 0.0) THEN
                  qCF(i,j,k) = 0.0
                  cloud_fraction_frozen(i,j,k) = 0.0
                  IF ( qCL(i,j,k) <= q_CC_tol ) THEN
                    area_cloud_fraction(i,j,k) = 0.0
                    bulk_cloud_fraction(i,j,k) = 0.0
                  END IF
                END IF

              END DO
            END DO

          END IF

!----------------------------------------------------------------------
! [1.3.12]: Aerosol.
!----------------------------------------------------------------------

          IF (Code == 90 .AND. k >= 1 .AND. k <= bl_levels) THEN

            DO j = 1, rows
              DO i = 1, row_length

                murk(i,j,k) = murk(i,j,k) + Weight * IAUIncFld(Addr)
                ! Apply upper (AEROMAX) and lower (AERO0) limits:
                murk(i,j,k) = MIN(MAX(AERO0,murk(i,j,k)),AEROMAX)
                Addr = Addr + 1

              END DO
            END DO

          END IF

!----------------------------------------------------------------------
! [1.3.13]: qT
!----------------------------------------------------------------------

          IF (Code == 18001 .AND. k >= 1 .AND. k <= wet_levels) THEN

            DO j = 1, rows
              DO i = 1, row_length

                qT_plus(i,j,k) = qT(i,j,k) + Weight * IAUIncFld(Addr)
                Addr = Addr + 1

              END DO
            END DO

          END IF

!-----------------------------------------------------------------
! [1.3.14]: Ozone tracer.                              
!-----------------------------------------------------------------
                                                                    
          IF (Code == 480 .AND. k >= 1 .AND. k <= model_levels) THEN       
                                                      
            DO j = 1, rows                                
              DO i = 1, row_length                           
                                                                 
                ozone_tracer(i,j,k) = ozone_tracer(i,j,k) +             &
     &    Weight * IAUIncFld(Addr)                
                Addr = Addr + 1                                   
                ! check for non-positive ozone tracer:               
                IF (ozone_tracer(i,j,k) <= 0.0) THEN    
                    ozone_tracer(i,j,k) = 0.0             
                   IF(L_IAU_SetOzoneMin) THEN             
                     ozone_tracer(i,j,k) = oz_min           
                   ENDIF                                     
                ENDIF                                         
                                                               
              END DO                                              
            END DO                                                    
                                                                        
          END IF                                                        
                                                                        
!----------------------------------------------------------------------   
! [1.3.15]: Update address in IAU D1 array.
!----------------------------------------------------------------------

          D1_IAU_addr = D1_IAU_addr + LocFldLen

        END DO ! FieldNum

!----------------------------------------------------------------------
! [1.4]: INDIRECT field updates by calculation from other fields.
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! [1.4.1]: exner/p.
!----------------------------------------------------------------------

        IF (L_IAU_CalcExnerIncs) THEN ! Calculate exner from p

          DO k = 1, model_levels+1
            DO j = 1, rows
              DO i = 1, row_length

                exner(i,j,k) = ( p(i,j,k) * p_zero_recip )**kappa

              END DO
            END DO
          END DO

        ELSE                          ! Calculate p from exner

          DO k = 1, model_levels+1
            DO j = 1, rows
              DO i = 1, row_length

                p(i,j,k) = exner(i,j,k)**kappa_recip * p_zero

              END DO
            END DO
          END DO

        END IF

        ! Update pstar:
! DEPENDS ON: calc_p_star
        CALL Calc_P_star ( r_theta_levels, r_rho_levels,                &
                                                            ! in
     &                     p, rho, g,                                   &
                                                            ! in
     &                     row_length, rows, model_levels,              &
                                                            ! in
     &                     offx, offy, halo_i, halo_j,                  &
                                                            ! in
     &                     p_star )                         ! out

        ! Update exner on theta levels:
! DEPENDS ON: calc_exner_at_theta
        CALL Calc_Exner_at_theta ( r_theta_levels, r_rho_levels,        &
                                                                   ! in
     &                             exner,                               &
                                                                   ! in
     &                             row_length, rows, model_levels,      &
                                                                   ! in
     &                             offx, offy, halo_i, halo_j,          &
                                                                   ! in
     &                             exner_theta_levels, .FALSE. )   ! out

        ! Get p on theta levels from exner on theta levels:
! DEPENDS ON: calc_p_from_exner
        CALL Calc_P_from_Exner ( p_theta_levels,                        &
                                                                 ! out
     &                           kappa, p_zero,                         &
                                                                 ! in
     &                           row_length, rows, model_levels,        &
                                                                 ! in
     &                           offx, offy,                            &
                                                                 ! in
     &                           exner_theta_levels,.FALSE. )    ! in

!----------------------------------------------------------------------
! [1.4.2]: thetaV.
!----------------------------------------------------------------------

        IF (L_IAU_CalcThetaIncs .OR. L_IAU_CalcRhoIncs) THEN

          ALLOCATE ( thetaV_new(row_length,rows,model_levels) )

          ! Add hydrostatic thetaV increments:
          DO k = 1, model_levels
            DO j = 1, rows
              DO i = 1, row_length

                IF (k /= model_levels) THEN
                  delta_r = r_rho_levels(i,j,k+1)                       &
     &                    - r_rho_levels(i,j,k  )
                ELSE
                  delta_r = 2.0                                         &
     &                    * ( r_theta_levels(i,j,k)                     &
     &                      - r_rho_levels  (i,j,k)                     &
     &                      )
                END IF

                thetaV_new(i,j,k) = thetaV_old(i,j,k)                   &
     &                            - g/cp * delta_r                      &
     &                            * ( 1.0                               &
     &                              / ( exner    (i,j,k+1)              &
     &                                - exner    (i,j,k  )              &
     &                              )                                   &
     &                              - 1.0                               &
     &                              / ( exner_old(i,j,k+1)              &
     &                                - exner_old(i,j,k  )              &
     &                              ) )

              END DO
            END DO
          END DO

        END IF

!----------------------------------------------------------------------
! [1.4.3]: rho.
!----------------------------------------------------------------------

        IF (L_IAU_CalcRhoIncs) THEN

          DO k = 1, model_levels
            DO j = 1, rows
              DO i = 1, row_length

                ! Interpolate old and new thetaVs to rho points:
                IF (k /= 1) THEN

                  Weight_1 = r_rho_levels  (i,j,k  )                    &
     &                     - r_theta_levels(i,j,k-1)
                  Weight_2 = r_theta_levels(i,j,k  )                    &
     &                     - r_rho_levels  (i,j,k  )

                  thetaV_old_int = ( Weight_1 * thetaV_old(i,j,k  )     &
     &                             + Weight_2 * thetaV_old(i,j,k-1)     &
     &                             )                                    &
     &                           / ( Weight_1 + Weight_2 )

                  thetaV_new_int = ( Weight_1 * thetaV_new(i,j,k  )     &
     &                             + Weight_2 * thetaV_new(i,j,k-1)     &
     &                             )                                    &
     &                           / ( Weight_1 + Weight_2 )

                ELSE

                  thetaV_old_int = thetaV_old(i,j,k)
                  thetaV_new_int = thetaV_new(i,j,k)

                END IF

                Term_1 = exner    (i,j,k)**exner_exp / thetaV_new_int
                Term_2 = exner_old(i,j,k)**exner_exp / thetaV_old_int

                rho(i,j,k) = rho(i,j,k)                                 &
     &                     + ( Term_1 - Term_2 )                        &
     &                     * r_rho_levels(i,j,k)**2                     &
     &                     * p_zero / (kappa * cp)

              END DO
            END DO
          END DO

        END IF

!----------------------------------------------------------------------
! [1.4.4]: theta.
!----------------------------------------------------------------------

        IF (L_IAU_CalcThetaIncs) THEN

          DO k = 1, wet_levels
            DO j = 1, rows
              DO i = 1, row_length

                theta(i,j,k) = thetaV_new(i,j,k)                        &
     &                       / ( 1.0 + C_Virtual * q(i,j,k) )

              END DO
            END DO
          END DO


          DO k = wet_levels+1, model_levels
            DO j = 1, rows
              DO i = 1, row_length

                theta(i,j,k) = thetaV_new(i,j,k)

              END DO
            END DO
          END DO

        END IF

!----------------------------------------------------------------------
! [1.4.5]: Deallocate work arrays.
!----------------------------------------------------------------------

        IF (L_IAU_CalcThetaIncs .OR. L_IAU_CalcRhoIncs) THEN
          DEALLOCATE (exner_old )
          DEALLOCATE (thetaV_old)
          DEALLOCATE (thetaV_new)
        END IF

!----------------------------------------------------------------------
! [1.4.6]: q, qCL, qCF, Cl & Cf
!----------------------------------------------------------------------

        IF(L_IAU_CalcCloudIncs) THEN

!----------------------------------------------------------------------
! [1.4.6.1]: obtain initial and incremented values from Var scheme
!----------------------------------------------------------------------

          ! Use Var cloud diagnostic to obtain  initial T, from qT & Tl
          Var_qT  = RESHAPE(qT(1:row_length, 1:rows, 1:wet_levels),     &
     &                      (/field_size/))
! DEPENDS ON: gen_vardiagcloud
          CALL Gen_VarDiagCloud(                                        &
     &        field_size,                                               &
     &        RESHAPE(p_theta_levels_old, (/field_size/)),              &
     &        RESHAPE(SPREAD(SPREAD(                                    &
     &                rhcrit(1:wet_levels),1,rows),1,row_length),       &
     &                (/field_size/)),                                  &
     &        L_IAU_IncrementIce,                                       &
     &        Var_qT,                                                   &
     &        CMessage,                                                 &
     &        ICode,                                                    &
     &        TL=VarTemp)
          ! VarTemp now holds initial T values
          IF (ICode /= 0) THEN

! DEPENDS ON: ereport
            CALL EReport (RoutineName, ICode, CMessage)

          END IF

          IF (L_IAU_IncrementIce) THEN

          ! Use Var cloud diagnostic to obtain  initial qCL,qCF,Cl & Cf
            qcl_0  = RESHAPE(qCL (1:row_length, 1:rows, 1:wet_levels),  &
     &                       (/field_size/))
            qcf_0  = RESHAPE(qCF (1:row_length, 1:rows, 1:wet_levels),  &
     &                       (/field_size/))
! DEPENDS ON: gen_vardiagcloud
            CALL Gen_VarDiagCloud(                                      &
     &          field_size,                                             &
     &          RESHAPE(p_theta_levels_old, (/field_size/)),            &
     &          RESHAPE(SPREAD(SPREAD(                                  &
     &                  rhcrit(1:wet_levels),1,rows),1,row_length),     &
     &                  (/field_size/)),                                &
     &          L_IAU_IncrementIce,                                     &
     &          Var_qT,                                                 &
     &          CMessage,                                               &
     &          ICode,                                                  &
     &          CL    = cloud_fraction_liquid_0,                        &
     &          qCL   = qCL_0,                                          &
     &          BGqCL = RESHAPE(qCL (1:row_length,1:rows,1:wet_levels), &
     &                          (/field_size/)),                        &
     &          qCF   = qCF_0,                                          &
     &          CF    = cloud_fraction_frozen_0,                        &
     &          BGqCF = RESHAPE(qCF (1:row_length,1:rows,1:wet_levels), &
     &                          (/field_size/)),                        &
     &          BGT   = Var_BGT,                                        &
     &          T     = VarTemp)

          ELSE

            ! Use Var cloud diagnostic to obtain  initial qCL & Cl
! DEPENDS ON: gen_vardiagcloud
            CALL Gen_VarDiagCloud(                                      &
     &          field_size,                                             &
     &          RESHAPE(p_theta_levels_old, (/field_size/)),            &
     &          RESHAPE(SPREAD(SPREAD(                                  &
     &                  rhcrit(1:wet_levels),1,rows),1,row_length),     &
     &                  (/field_size/)),                                &
     &          L_IAU_IncrementIce,                                     &
     &          Var_qT,                                                 &
     &          CMessage,                                               &
     &          ICode,                                                  &
     &          CL    = cloud_fraction_liquid_0,                        &
     &          qCL   = qCL_0,                                          &
     &          T     = VarTemp)

          END IF
          IF (ICode /= 0) THEN

! DEPENDS ON: ereport
            CALL EReport (RoutineName, ICode, CMessage)

          END IF

          ! Assign VarTemp to hold incremented T values:
          IF (L_IAU_IncrementIce) THEN

            ! Assign Var_BGT to store the T associated with qcl_0
            ! Use qCL_Plus as a temporary variable
            ! This should be done in the qcl / Cl version too but
            !   will go in later as a bug-fix so that bit-comparison
            !   can be maintained when the code for qcf incrementing is
            !   implemented but switched off
            qCL_Plus = VarTemp
            VarTemp  = VarTemp + (                                      &
     &           RESHAPE(                                               &
     &           exner_theta_levels(1:row_length, 1:rows, 1:wet_levels) &
     &           * theta(1:row_length, 1:rows, 1:wet_levels),           &
     &          (/field_size/)) - Var_BGT)
            Var_BGT  = qCL_Plus

          ELSE

            VarTemp = RESHAPE(                                          &
     &           exner_theta_levels(1:row_length, 1:rows, 1:wet_levels) &
     &           * theta(1:row_length, 1:rows, 1:wet_levels),           &
     &          (/field_size/))

          END IF
          IF (L_IAU_IncrementIce) THEN

            ! Calculate maximum qcf values (for use later)
            qCL_Plus = 0.0
            qCF_Plus = 0.0
            Var_qT = RESHAPE(qT_plus(1:row_length,1:rows,1:wet_levels)  &
     &                       +  qCF (1:row_length,1:rows,1:wet_levels), &
     &                       (/field_size/))
            qCF_max = Var_qT
! DEPENDS ON: gen_vardiagcloud
            CALL Gen_VarDiagCloud(                                      &
     &          field_size,                                             &
     &          RESHAPE(p_theta_levels                                  &
     &                  (1:row_length, 1:rows, 1:wet_levels),           &
     &                  (/field_size/)),                                &
     &          RESHAPE(SPREAD(SPREAD(                                  &
     &                  rhcrit(1:wet_levels),1,rows),1,row_length),     &
     &                  (/field_size/)),                                &
     &          L_IAU_IncrementIce,                                     &
     &          Var_qT,                                                 &
     &          CMessage,                                               &
     &          ICode,                                                  &
     &          qCL    = qCL_Plus,                                      &
     &          BGqCL  = RESHAPE(qCL(1:row_length,1:rows,1:wet_levels), &
     &                           (/field_size/)),                       &
     &          qCF    = qCF_Plus,                                      &
     &          qCFmax = qCF_max,                                       &
     &          BGqCF  = RESHAPE(qCF(1:row_length,1:rows,1:wet_levels), &
     &                           (/field_size/)),                       &
     &          BGT    = Var_BGT,                                       &
     &          T      = VarTemp)
            IF (ICode /= 0) THEN

! DEPENDS ON: ereport
              CALL EReport (RoutineName, ICode, CMessage)

            END IF

            qCL_Plus  = qCL_0
            qCF_Plus  = qCF_0
            Var_qT = RESHAPE(qT     (1:row_length,1:rows,1:wet_levels)  &
     &                       +  qCF (1:row_length,1:rows,1:wet_levels), &
     &                       (/field_size/))
            WHERE (qCF_max <= 0.0)

              qcf_max = Var_qT

            ENDWHERE

          END IF

          Var_qT = RESHAPE(qT_plus(1:row_length, 1:rows, 1:wet_levels), &
     &                     (/field_size/))
          IF (L_IAU_IncrementIce) THEN

            ! Use Var cloud diagnostic to obtain  final qCL,qCF,Cl & Cf
! DEPENDS ON: gen_vardiagcloud
            CALL Gen_VarDiagCloud(                                      &
     &          field_size,                                             &
     &          RESHAPE(p_theta_levels                                  &
     &                  (1:row_length, 1:rows, 1:wet_levels),           &
     &                  (/field_size/)),                                &
     &          RESHAPE(SPREAD(SPREAD(                                  &
     &                  rhcrit(1:wet_levels),1,rows),1,row_length),     &
     &                  (/field_size/)),                                &
     &          L_IAU_IncrementIce,                                     &
     &          Var_qT,                                                 &
     &          CMessage,                                               &
     &          ICode,                                                  &
     &          qCL    = qCL_Plus,                                      &
     &          BGqCL  = RESHAPE(qCL(1:row_length,1:rows,1:wet_levels), &
     &                           (/field_size/)),                       &
     &          qCF    = qCF_Plus,                                      &
     &          BGqCF  = RESHAPE(qCF(1:row_length,1:rows,1:wet_levels), &
     &                           (/field_size/)),                       &
     &          qCFmax = qCF_max,                                       &
     &          CL     = cloud_fraction_liquid_plus,                    &
     &          CF     = cloud_fraction_frozen_plus,                    &
     &          BGT    = Var_BGT,                                       &
     &          T      = VarTemp)

          ELSE

            ! Use Var cloud diagnostic to obtain  final qCL & Cl
! DEPENDS ON: gen_vardiagcloud
            CALL Gen_VarDiagCloud(                                      &
     &          field_size,                                             &
     &          RESHAPE(p_theta_levels                                  &
     &                  (1:row_length, 1:rows, 1:wet_levels),           &
     &                  (/field_size/)),                                &
     &          RESHAPE(SPREAD(SPREAD(                                  &
     &                  rhcrit(1:wet_levels),1,rows),1,row_length),     &
     &                  (/field_size/)),                                &
     &          L_IAU_IncrementIce,                                     &
     &          Var_qT,                                                 &
     &          CMessage,                                               &
     &          ICode,                                                  &
     &          CL     = cloud_fraction_liquid_plus,                    &
     &          qCL    = qCL_Plus,                                      &
     &          T      = VarTemp)

          END IF
          IF (ICode /= 0) THEN

! DEPENDS ON: ereport
            CALL EReport (RoutineName, ICode, CMessage)

          END IF

!----------------------------------------------------------------------
! [1.4.6.2]: apply increments
!----------------------------------------------------------------------

          ! qcl
          IF (L_IAU_ScaleCloud) THEN
          ! Scale qcl increments to be physical values
            ALLOCATE (initialScaling (row_length, rows, wet_levels))
            ALLOCATE (scalingCoeff   (row_length, rows, wet_levels))
            ALLOCATE (qcl_Inc        (row_length, rows, wet_levels))
            WHERE(RESHAPE(qCL_0,(/row_length,rows,wet_levels/)) >       &
     &        qCL(1:row_length, 1:rows, 1:wet_levels) .AND.             &
     &        RESHAPE(qCL_plus - qCL_0,(/row_length,rows,wet_levels/))  &
     &        < 0.0)
              qcl_Inc = RESHAPE(qCL_plus - qCL_0,                       &
     &                          (/row_length, rows, wet_levels/))
              WHERE (qCL(1:row_length, 1:rows, 1:wet_levels) > 0.0)
                initialScaling = -qCL(1:row_length,1:rows,1:wet_levels)*&
     &              (1.0 - exp(qcl_Inc /                                &
     &                         qCL(1:row_length, 1:rows, 1:wet_levels)))
                scalingCoeff = -qCL(1:row_length, 1:rows, 1:wet_levels)*&
     &          (1.0-exp(-RESHAPE(qCL_0,(/row_length,rows,wet_levels/))/&
     &          qCL(1:row_length, 1:rows, 1:wet_levels)))
              ! The following calculation should remain >0 so no error
              ! trap needed
                scalingCoeff =                                          &
     &              (scalingCoeff +                                     &
     &               qCL(1:row_length, 1:rows, 1:wet_levels)) /         &
     &              (scalingCoeff +                                     &
     &               RESHAPE(qCL_0,(/row_length,rows,wet_levels/)))
              ELSEWHERE
                initialScaling = 0.0
                scalingCoeff = 0.0
              END WHERE
              qcl_Inc = initialScaling -                                &
     &                  (initialScaling - qcl_Inc) * scalingCoeff
            ELSEWHERE
              qcl_Inc = RESHAPE(qCL_plus - qCL_0,                       &
     &                          (/row_length, rows, wet_levels/))
            END WHERE
            qCL(1:row_length, 1:rows, 1:wet_levels) =                   &
     &          qCL(1:row_length, 1:rows, 1:wet_levels) + qcl_Inc
          ELSE
            qCL(1:row_length, 1:rows, 1:wet_levels) =                   &
     &          qCL(1:row_length, 1:rows, 1:wet_levels) +               &
     &          RESHAPE(qCL_plus - qCL_0,                               &
     &                  (/row_length, rows, wet_levels/))
          END IF

          ! qcf
          IF (L_IAU_IncrementIce) THEN

            qCF(1:row_length, 1:rows, 1:wet_levels) =                   &
     &          qCF(1:row_length, 1:rows, 1:wet_levels) +               &
     &          RESHAPE(qCF_plus - qCF_0,                               &
     &                  (/row_length, rows, wet_levels/))

          END IF
          
          ! Store input values of area / bulk cloud fractions
          IF (L_CLD_AREA) THEN
            
            WHERE (bulk_cloud_fraction(1:row_length,1:rows,1:wet_levels)&
     &          <= 0.0)

              AreaDivBulk = 1.0

            ELSEWHERE

              AreaDivBulk =                                             &
     &          area_cloud_fraction(1:row_length,1:rows,1:wet_levels) / &
     &          bulk_cloud_fraction(1:row_length,1:rows,1:wet_levels)

            ENDWHERE

          END IF

          ! Cl
          IF (L_IAU_ScaleCloud) THEN
          ! Scale Cl increments to be physical values
            ALLOCATE (Cl_Inc (row_length, rows, wet_levels))
            ALLOCATE (Cf_Inc (row_length, rows, wet_levels))
            WHERE(RESHAPE(cloud_fraction_liquid_0,                      &
     &                    (/row_length,rows,wet_levels/)) >             &
     &        cloud_fraction_liquid(1:row_length, 1:rows, 1:wet_levels) &
     &        .AND. RESHAPE(cloud_fraction_liquid_plus -                &
     &                      cloud_fraction_liquid_0,                    &
     &                      (/row_length, rows, wet_levels/)) < 0.0)
              Cl_Inc = RESHAPE(cloud_fraction_liquid_plus -             &
     &                         cloud_fraction_liquid_0,                 &
     &                         (/row_length, rows, wet_levels/))
              WHERE (cloud_fraction_liquid                              &
     &               (1:row_length, 1:rows, 1:wet_levels) > 0.0)
                initialScaling = - cloud_fraction_liquid                &
     &                             (1:row_length, 1:rows, 1:wet_levels) &
     &                           * (1.0 - exp(Cl_Inc /                  &
     &                              cloud_fraction_liquid               &
     &                            (1:row_length, 1:rows, 1:wet_levels)))
                scalingCoeff =                                          &
     &              - cloud_fraction_liquid                             &
     &                (1:row_length, 1:rows, 1:wet_levels) *            &
     &                (1.0 - exp(-RESHAPE(cloud_fraction_liquid_0,      &
     &                            (/row_length,rows,wet_levels/)) /     &
     &                cloud_fraction_liquid                             &
     &                (1:row_length, 1:rows, 1:wet_levels)))
              ! The following calculation should remain >0 so
              ! no error trap needed
                scalingCoeff =                                          &
     &              (scalingCoeff +                                     &
     &               cloud_fraction_liquid                              &
     &               (1:row_length, 1:rows, 1:wet_levels)) /            &
     &              (scalingCoeff + RESHAPE(cloud_fraction_liquid_0,    &
     &                              (/row_length,rows,wet_levels/)))
              ELSEWHERE
                initialScaling = 0.0
                scalingCoeff = 0.0
              END WHERE
              Cl_Inc = initialScaling - (initialScaling - Cl_Inc) *     &
     &                 scalingCoeff
            ELSE WHERE(RESHAPE(cloud_fraction_liquid_0,                 &
     &                         (/row_length,rows,wet_levels/)) <        &
     &                 cloud_fraction_liquid                            &
     &                 (1:row_length, 1:rows, 1:wet_levels) .AND.       &
     &                 RESHAPE(cloud_fraction_liquid_plus -             &
     &                         cloud_fraction_liquid_0,                 &
     &                         (/row_length, rows, wet_levels/)) > 0.0)
              Cl_Inc = RESHAPE(cloud_fraction_liquid_plus -             &
     &                         cloud_fraction_liquid_0,                 &
     &                         (/row_length, rows, wet_levels/))
              WHERE ((1.0 - cloud_fraction_liquid                       &
     &                      (1:row_length, 1:rows, 1:wet_levels)) > 0.0)
                initialScaling = (1.0 - cloud_fraction_liquid           &
     &                            (1:row_length, 1:rows, 1:wet_levels))*&
     &                           (1.0 - exp(-Cl_Inc / (1.0 -            &
     &                             cloud_fraction_liquid                &
     &                           (1:row_length, 1:rows, 1:wet_levels))))
                scalingCoeff = (1.0 - cloud_fraction_liquid             &
     &              (1:row_length, 1:rows, 1:wet_levels)) *             &
     &              (1.0 - exp(-(1.0 - RESHAPE(cloud_fraction_liquid_0, &
     &                                 (/row_length,rows,wet_levels/)))/&
     &              (1.0 - cloud_fraction_liquid                        &
     &                     (1:row_length, 1:rows, 1:wet_levels))))
                scalingCoeff =                                          &
                    ((1.0 - cloud_fraction_liquid                       &
     &                      (1:row_length, 1:rows, 1:wet_levels)) -     &
     &               scalingCoeff) / ((1.0 -                            &
     &               RESHAPE(cloud_fraction_liquid_0,                   &
     &                       (/row_length,rows,wet_levels/)))           &
     &               - scalingCoeff)
              ELSEWHERE
                initialScaling = 0.0
                scalingCoeff = 0.0
              END WHERE
              Cl_Inc = initialScaling + (Cl_Inc - initialScaling) *     &
     &                 scalingCoeff
            ELSEWHERE
              Cl_Inc = RESHAPE(cloud_fraction_liquid_plus -             &
     &                         cloud_fraction_liquid_0,                 &
     &                         (/row_length, rows, wet_levels/))
            END WHERE
            ! Reuse cloud_fraction_liquid_0 to store initial values for
            ! use in Cb incrementing
            cloud_fraction_liquid_0 = RESHAPE(cloud_fraction_liquid     &
     &          (1:row_length, 1:rows, 1:wet_levels), (/field_size/))
            cloud_fraction_liquid(1:row_length, 1:rows, 1:wet_levels) = &
     &          cloud_fraction_liquid(1:row_length,1:rows,1:wet_levels) &
     &          + Cl_Inc
          ELSE
            cloud_fraction_liquid(1:row_length, 1:rows, 1:wet_levels) = &
     &          cloud_fraction_liquid(1:row_length,1:rows,1:wet_levels) &
     &          + RESHAPE(cloud_fraction_liquid_plus -                  &
     &                    cloud_fraction_liquid_0,                      &
     &                    (/row_length, rows, wet_levels/))
          END IF


          ! Set Cl <= 1.0
          WHERE (                                                       &
     &        cloud_fraction_liquid(1:row_length, 1:rows, 1:wet_levels) &
     &        > 1.0)

            cloud_fraction_liquid(1:row_length, 1:rows, 1:wet_levels)   &
     &          = 1.0

          ENDWHERE

          ! Set qcl & Cl >= 0.0
          WHERE (                                                       &
     &        cloud_fraction_liquid(1:row_length, 1:rows, 1:wet_levels) &
     &        <= 0.0                                                    &
     &        .OR. qCL(1:row_length, 1:rows, 1:wet_levels) <= q_CC_tol)

            qCL(1:row_length, 1:rows, 1:wet_levels) = 0.0
            cloud_fraction_liquid(1:row_length, 1:rows, 1:wet_levels)   &
     &          = 0.0

          ENDWHERE

          IF (L_IAU_IncrementIce) THEN

            ! Cf
            IF (L_IAU_ScaleCloud) THEN
            ! Scale Cf increments to be physical values
              WHERE(RESHAPE(cloud_fraction_frozen_0,                    &
     &                      (/row_length,rows,wet_levels/)) >           &
     &          cloud_fraction_frozen(1:row_length,1:rows,1:wet_levels) &
     &          .AND. RESHAPE(cloud_fraction_frozen_plus -              &
     &                        cloud_fraction_frozen_0,                  &
     &                        (/row_length,rows,wet_levels/)) < 0.0)
                Cf_Inc = RESHAPE(cloud_fraction_frozen_plus -           &
     &                           cloud_fraction_frozen_0,               &
     &                           (/row_length,rows,wet_levels/))
                WHERE (cloud_fraction_frozen                            &
     &                 (1:row_length, 1:rows, 1:wet_levels) > 0.0)
                  initialScaling = - cloud_fraction_frozen              &
     &                               (1:row_length,1:rows,1:wet_levels) &
     &                             * (1.0 - exp(Cf_Inc /                &
     &                                cloud_fraction_frozen             &
     &                              (1:row_length,1:rows,1:wet_levels)))
                  scalingCoeff =                                        &
     &                - cloud_fraction_frozen                           &
     &                  (1:row_length, 1:rows, 1:wet_levels) *          &
     &                  (1.0 - exp(-RESHAPE(cloud_fraction_frozen_0,    &
     &                              (/row_length,rows,wet_levels/)) /   &
     &                  cloud_fraction_frozen                           &
     &                  (1:row_length, 1:rows, 1:wet_levels)))
                ! The following calculation should remain >0 so
                ! no error trap needed
                  scalingCoeff =                                        &
     &                (scalingCoeff +                                   &
     &                 cloud_fraction_frozen                            &
     &                 (1:row_length, 1:rows, 1:wet_levels)) /          &
     &                (scalingCoeff + RESHAPE(cloud_fraction_frozen_0,  &
     &                                (/row_length,rows,wet_levels/)))
                ELSEWHERE
                  initialScaling = 0.0
                  scalingCoeff = 0.0
                END WHERE
                Cf_Inc = initialScaling - (initialScaling - Cf_Inc) *   &
     &                   scalingCoeff
              ELSE WHERE(RESHAPE(cloud_fraction_frozen_0,               &
     &                           (/row_length,rows,wet_levels/)) <      &
     &                   cloud_fraction_frozen                          &
     &                   (1:row_length, 1:rows, 1:wet_levels) .AND.     &
     &                   RESHAPE(cloud_fraction_frozen_plus -           &
     &                           cloud_fraction_frozen_0,               &
     &                           (/row_length,rows,wet_levels/)) > 0.0)
                Cf_Inc = RESHAPE(cloud_fraction_frozen_plus -           &
     &                           cloud_fraction_frozen_0,               &
     &                           (/row_length,rows,wet_levels/))
                WHERE ((1.0 - cloud_fraction_frozen                     &
     &                        (1:row_length,1:rows,1:wet_levels)) > 0.0)
                  initialScaling = (1.0 - cloud_fraction_frozen         &
     &                              (1:row_length,1:rows,1:wet_levels))*&
     &                             (1.0 - exp(-Cf_Inc / (1.0 -          &
     &                               cloud_fraction_frozen              &
     &                             (1:row_length,1:rows,1:wet_levels))))
                  scalingCoeff = (1.0 - cloud_fraction_frozen           &
     &                (1:row_length, 1:rows, 1:wet_levels)) *           &
     &                (1.0 - exp(-(1.0 -RESHAPE(cloud_fraction_frozen_0,&
     &                             (/row_length,rows,wet_levels/))) /   &
     &                (1.0 - cloud_fraction_frozen                      &
     &                       (1:row_length, 1:rows, 1:wet_levels))))
                  scalingCoeff =                                        &
     &                ((1.0 - cloud_fraction_frozen                     &
     &                        (1:row_length, 1:rows, 1:wet_levels)) -   &
     &                 scalingCoeff) / ((1.0 -                          &
     &                 RESHAPE(cloud_fraction_frozen_0,                 &
     &                 (/row_length,rows,wet_levels/))) - scalingCoeff)
                ELSEWHERE
                  initialScaling = 0.0
                  scalingCoeff = 0.0
                END WHERE
                Cf_Inc = initialScaling + (Cf_Inc - initialScaling) *   &
     &                   scalingCoeff
              ELSEWHERE
                Cf_Inc = RESHAPE(cloud_fraction_frozen_plus -           &
     &                           cloud_fraction_frozen_0,               &
     &                           (/row_length,rows,wet_levels/))
              END WHERE
            ! Reuse cloud_fraction_frozen_0 to store initial values for
            ! use in Cb incrementing
              cloud_fraction_frozen_0 = RESHAPE(cloud_fraction_frozen   &
       &          (1:row_length, 1:rows, 1:wet_levels), (/field_size/))
              cloud_fraction_frozen(1:row_length,1:rows,1:wet_levels) = &
     &          cloud_fraction_frozen(1:row_length,1:rows,1:wet_levels) &
     &          + Cf_Inc
            ELSE
            cloud_fraction_frozen(1:row_length, 1:rows, 1:wet_levels) = &
     &          cloud_fraction_frozen(1:row_length,1:rows,1:wet_levels) &
     &          + RESHAPE(cloud_fraction_frozen_plus -                  &
     &           cloud_fraction_frozen_0,                               &
     &           (/row_length, rows, wet_levels/))
            END IF

            ! Set Cf <= 1.0
            WHERE (                                                     &
     &          cloud_fraction_frozen(1:row_length,1:rows,1:wet_levels) &
     &          > 1.0)

              cloud_fraction_frozen(1:row_length, 1:rows, 1:wet_levels) &
     &            = 1.0

            ENDWHERE

            ! Set qcf & Cf >= 0.0
            WHERE (                                                     &
     &          cloud_fraction_frozen(1:row_length,1:rows,1:wet_levels) &
     &          <= 0.0                                                  &
     &          .OR. qCF(1:row_length,1:rows,1:wet_levels) <= q_CC_tol)

              qCF(1:row_length, 1:rows, 1:wet_levels) = 0.0
              cloud_fraction_frozen(1:row_length, 1:rows, 1:wet_levels) &
     &            = 0.0

            ENDWHERE

          END IF

          ! Set other cloud variables to zero if qCL & qCF <= q_CC_tol
          WHERE (   qCL(1:row_length, 1:rows, 1:wet_levels) <= q_CC_tol &
     &        .AND. qCF(1:row_length, 1:rows, 1:wet_levels) <= q_CC_tol)

            qCF(1:row_length, 1:rows, 1:wet_levels) = 0.0
            cloud_fraction_frozen(1:row_length, 1:rows, 1:wet_levels)   &
     &          = 0.0
            area_cloud_fraction(1:row_length, 1:rows, 1:wet_levels)     &
     &          = 0.0
            bulk_cloud_fraction(1:row_length, 1:rows, 1:wet_levels)     &
     &          = 0.0

          ENDWHERE
          IF (L_IAU_IncrementIce) THEN

            WHERE (   qCL(1:row_length,1:rows,1:wet_levels) <= q_CC_tol &
     &          .AND. qCF(1:row_length,1:rows,1:wet_levels) <= q_CC_tol)

              qCL(1:row_length, 1:rows, 1:wet_levels) = 0.0
              cloud_fraction_liquid(1:row_length, 1:rows, 1:wet_levels) &
     &            = 0.0

            ENDWHERE

          END IF

          ! Apply stored area / bulk cloud fractions to output
          IF (L_IAU_ScaleCloud) THEN
          ! Scale Cb increments to be physical values
          ! Use Cb_inc and Cf_Inc as temporary storage for incremented
          ! and unincremented MIN(Cl+Cf, 1.0)
            ALLOCATE (Cb_Inc (row_length, rows, wet_levels))
            Cb_Inc = MIN(cloud_fraction_liquid                          &
     &                   (1:row_length, 1:rows, 1:wet_levels) +         &
     &                   cloud_fraction_frozen                          &
     &                   (1:row_length, 1:rows, 1:wet_levels), 1.0)
            IF (L_IAU_IncrementIce) THEN
              Cf_Inc = RESHAPE(MIN(cloud_fraction_liquid_0 +            &
     &                             cloud_fraction_frozen_0,1.0),        &
     &                         (/row_length,rows,wet_levels/))
            ELSE
              Cf_Inc = MIN(RESHAPE(cloud_fraction_liquid_0,             &
     &                         (/row_length,rows,wet_levels/)) +        &
     &                         cloud_fraction_frozen                    &
     &                         (1:row_length, 1:rows, 1:wet_levels), 1.0)
            END IF
            WHERE(Cf_Inc >                                              &
     &            bulk_cloud_fraction(1:row_length,1:rows,1:wet_levels) &
     &            .AND. Cb_Inc - Cf_Inc < 0.0)
              Cb_Inc = Cb_Inc - Cf_Inc
              WHERE (bulk_cloud_fraction                                &
     &               (1:row_length, 1:rows, 1:wet_levels) > 0.0)
                initialScaling = - bulk_cloud_fraction                  &
     &                             (1:row_length,1:rows,1:wet_levels) * &
     &                           (1.0 - exp(Cb_Inc /                    &
     &                                      bulk_cloud_fraction         &
     &                            (1:row_length, 1:rows, 1:wet_levels)))
                scalingCoeff = - bulk_cloud_fraction                    &
     &                               (1:row_length,1:rows,1:wet_levels) &
     &                             * (1.0 - exp(-Cf_Inc /               &
     &                               bulk_cloud_fraction                &
     &                            (1:row_length, 1:rows, 1:wet_levels)))
              ! The following calculation should remain >0 so no error
              !trap needed
                scalingCoeff =                                          &
     &           (scalingCoeff +                                        &
     &            bulk_cloud_fraction(1:row_length,1:rows,1:wet_levels))&
     &            / (scalingCoeff + Cf_Inc  )
              ELSEWHERE
                initialScaling = 0.0
                scalingCoeff = 0.0
              END WHERE
              Cb_Inc = initialScaling - (initialScaling - Cb_Inc) *     &
     &                                  scalingCoeff
            ELSE WHERE(Cf_Inc <                                         &
     &                 bulk_cloud_fraction                              &
     &                 (1:row_length, 1:rows, 1:wet_levels) .AND.       &
     &                 Cb_Inc - Cf_Inc > 0.0)
              Cb_Inc = Cb_Inc - Cf_Inc
              WHERE ((1.0 - bulk_cloud_fraction                         &
     &                      (1:row_length, 1:rows, 1:wet_levels)) > 0.0)
                initialScaling = (1.0 - bulk_cloud_fraction             &
     &                            (1:row_length, 1:rows, 1:wet_levels)) &
     &                         * (1.0 - exp(-Cb_Inc / (1.0 -            &
     &                                              bulk_cloud_fraction &
     &                           (1:row_length, 1:rows, 1:wet_levels))))
                scalingCoeff = (1.0 - bulk_cloud_fraction               &
     &                            (1:row_length, 1:rows, 1:wet_levels)) &
     &                          * (1.0 - exp(-(1.0 - Cf_Inc) / (1.0 -   &
     &                            bulk_cloud_fraction                   &
     &                           (1:row_length, 1:rows, 1:wet_levels))))
                scalingCoeff =                                          &
     &              ((1.0 - bulk_cloud_fraction                         &
     &                      (1:row_length, 1:rows, 1:wet_levels)) -     &
     &               scalingCoeff) /                                    &
     &              ((1.0 - Cf_Inc  ) - scalingCoeff)
              ELSEWHERE
                initialScaling = 0.0
                scalingCoeff = 0.0
              END WHERE
              Cb_Inc = initialScaling + (Cb_Inc - initialScaling) *     &
     &                                  scalingCoeff
            ELSEWHERE
              Cb_Inc = Cb_Inc - Cf_Inc
            END WHERE
            bulk_cloud_fraction(1:row_length, 1:rows, 1:wet_levels) =   &
     &          bulk_cloud_fraction(1:row_length, 1:rows, 1:wet_levels) &
     &          + Cb_Inc
          ELSE
            bulk_cloud_fraction(1:row_length, 1:rows, 1:wet_levels) =   &
     &         MIN(1.0,MAX(0.0,                                         &
     &          cloud_fraction_liquid                                   &
     &                              (1:row_length,1:rows,1:wet_levels)  &
     &          +                                                       &
     &          cloud_fraction_frozen                                   &
     &                              (1:row_length,1:rows,1:wet_levels)))
          END IF

          IF (.NOT. L_CLD_AREA) THEN

            area_cloud_fraction  (1:row_length, 1:rows, 1:wet_levels) = &
     &        bulk_cloud_fraction(1:row_length, 1:rows, 1:wet_levels)

          ELSE

            area_cloud_fraction(1:row_length, 1:rows, 1:wet_levels) =   &
     &          MIN(1.0,                                                &
     &          AreaDivBulk *                                           &
     &          bulk_cloud_fraction(1:row_length, 1:rows, 1:wet_levels))

          END IF

          ! q
          IF (L_IAU_IncrementIce) THEN

            IF (L_IAU_ScaleCloud) THEN

              q(1:row_length, 1:rows, 1:wet_levels) =                   &
     &            q(1:row_length, 1:rows, 1:wet_levels)     +           &
     &            (qT_plus - qT)                            -           &
     &            qcl_Inc -                                             &
     &            RESHAPE(qCF_plus - qCF_0,                             &
     &                    (/row_length, rows, wet_levels/))

            ELSE

              q(1:row_length, 1:rows, 1:wet_levels) =                   &
     &            q(1:row_length, 1:rows, 1:wet_levels)     +           &
     &            (qT_plus - qT)                            -           &
     &            RESHAPE(qCL_plus - qCL_0,                             &
     &                    (/row_length, rows, wet_levels/)) -           &
     &            RESHAPE(qCF_plus - qCF_0,                             &
     &                    (/row_length, rows, wet_levels/))

            END IF

          ELSE

            IF (L_IAU_ScaleCloud) THEN

              q(1:row_length, 1:rows, 1:wet_levels) =                   &
     &            q(1:row_length, 1:rows, 1:wet_levels)     +           &
     &            (qT_plus - qT)                            -           &
     &            qcl_Inc

            ELSE

              q(1:row_length, 1:rows, 1:wet_levels) =                   &
     &            q(1:row_length, 1:rows, 1:wet_levels)     +           &
     &            (qT_plus - qT)                            -           &
     &            RESHAPE(qCL_plus - qCL_0,                             &
     &                    (/row_length, rows, wet_levels/))

            END IF

          END IF
          ! Apply lower limit
          WHERE (q(1:row_length, 1:rows, 1:wet_levels) < q_min)

            q(1:row_length, 1:rows, 1:wet_levels) = q_min

          ENDWHERE

!----------------------------------------------------------------------
! [1.4.6.4]: Deallocate work arrays.
!----------------------------------------------------------------------

          DEALLOCATE (AreaDivBulk)
          DEALLOCATE (cloud_fraction_liquid_0)
          DEALLOCATE (cloud_fraction_liquid_plus)
          DEALLOCATE (p_theta_levels_old)
          DEALLOCATE (qCL_0)
          DEALLOCATE (qCL_plus)
          DEALLOCATE (qSatW)
          DEALLOCATE (qT)
          DEALLOCATE (qT_plus)
          DEALLOCATE (Var_qT)
          DEALLOCATE (VarTemp)
          IF (L_IAU_IncrementIce) THEN

            DEALLOCATE (cloud_fraction_frozen_0)
            DEALLOCATE (cloud_fraction_frozen_plus)
            DEALLOCATE (qCF_0)
            DEALLOCATE (qCF_max)
            DEALLOCATE (qCF_plus)
            DEALLOCATE (Var_BGT)

          END IF
          IF (L_IAU_ScaleCloud) THEN

            DEALLOCATE (initialScaling)
            DEALLOCATE (scalingCoeff)
            DEALLOCATE (qcl_Inc)
            DEALLOCATE (Cl_Inc)
            DEALLOCATE (Cf_Inc)
            DEALLOCATE (Cb_Inc)

          END IF
        ELSE
!----------------------------------------------------------------------
! [1.4.6.5]: PC2 Calculations
!----------------------------------------------------------------------

          If (L_pc2 .and.                                               &
     &        (L_pc2_cond .or. L_pc2_cfl .or. L_pc2_cff)) then
            Allocate ( delta_q  (row_length,rows,wet_levels) )
            Allocate ( delta_qcl(row_length,rows,wet_levels) )
            Allocate ( delta_qcf(row_length,rows,wet_levels) )
            Allocate ( delta_t  (row_length,rows,wet_levels) )
            Allocate ( delta_p  (row_length,rows,wet_levels) )

            Do k = 1, wet_levels
              Do j = 1, rows
                Do i= 1, row_length
                  ! Form forcing increments of vapour, liquid,
                  ! temperature and pressure
                  delta_q  (i,j,k) = q  (i,j,k) - q_work  (i,j,k)
                  delta_qcl(i,j,k) = qcl(i,j,k) - qcl_work(i,j,k)
                  delta_qcf(i,j,k) = qcf(i,j,k) - qcf_work(i,j,k)
                  delta_t  (i,j,k) = exner_theta_levels(i,j,k)          &
     &                             * theta(i,j,k) - t_work(i,j,k)
                  delta_p  (i,j,k) = p_theta_levels(i,j,k)              &
     &                                            - p_work(i,j,k)
                End Do
              End Do
            End Do

            ! Now force condensation and cloud fraction updates.
            ! PC2 is not set up for prognostic area_cloud_fraction.

! DEPENDS ON: pc2_assim
            Call pc2_assim( wet_levels, row_length, rows, timestep      &
     &,                     l_pc2_cfl, l_pc2_cff                        &
     &,                     l_mixing_ratio, t_work                      &
     &,          bulk_cloud_fraction  (1:row_length,1:rows,1:wet_levels)&
     &,          cloud_fraction_liquid(1:row_length,1:rows,1:wet_levels)&
     &,          cloud_fraction_frozen(1:row_length,1:rows,1:wet_levels)&
     &,                     q_work, qcl_work, qcf_work, p_work          &
     &,                     delta_t, delta_q, delta_qcl, delta_qcf      &
     &,                     delta_p                                     &
     &                    )

            ! q_work, qcl_work and t_work have now been updated by the
            ! forcing and condensation terms together.
            ! cloud_fraction_liquid (and _frozen and _bulk) has been
            ! updated by the condensation. qcf_work and p_work are not
            ! updated.

            ! Now copy q_work, qcl_work and t_work variables back
            ! into the inout variables if we haven't a confident
            ! direct estimate of condensation from another source.
            If (L_pc2_cond) then
              Do k = 1, wet_levels
                Do j = 1, rows
                  Do i= 1, row_length
                    q  (i,j,k) = q_work  (i,j,k)
                    qcl(i,j,k) = qcl_work(i,j,k)
                    ! Remember the inout variable is theta, not temp.
                    theta(i,j,k)=t_work(i,j,k)/exner_theta_levels(i,j,k)
                  End Do
                End Do
              End Do
            End If  ! L_pc2_cond

            ! Now deallocate the arrays
            Deallocate(q_work  )
            Deallocate(qcl_work)
            Deallocate(qcf_work)
            Deallocate(t_work  )
            Deallocate(p_work  )
            Deallocate(delta_q  )
            Deallocate(delta_qcl)
            Deallocate(delta_qcf)
            Deallocate(delta_t  )
            Deallocate(delta_p  )

          End If  ! L_pc2 .and. (L_pc2_cond.or.L_pc2_cfl.or.L_pc2_cff)
        ENDIF

!----------------------------------------------------------------------
! [1.4.7]: If required, add the level 1 temperature increments to
!          TSoil(1), TStar and TStar_tile.
!----------------------------------------------------------------------

        IF (L_IAU_IncTStar) THEN

          IF (PrintStatus >= PrStatus_Normal) THEN
            WRITE (6,*) ''
            WRITE (6,*) 'IAU: Adding t1 incs to TStar and TSoil(1).'
          END IF

          ! Get increment to level one temperature:
          Addr = 1
          DO j = 1, rows
            DO i = 1, row_length
              t1_inc(Addr) = exner_theta_levels(i,j,1)                  &
     &                     * theta             (i,j,1)                  &
     &                     - t1_old(i,j)
              Addr = Addr + 1
            END DO
          END DO

          ! Add increments to TSoil(1) and TStar:
          ! If snow_depth <= 0.05 (kg/m**2)
          DO i = 1, land_field

            IF ( snow_depth(land_index(i)) <= 0.05 ) THEN

              TSoil(i) = TSoil(i) + t1_inc(land_index(i))
              TStar(land_index(i)) =  TStar(land_index(i))              &
                                   + t1_inc(land_index(i))

            END IF

          END DO

          IF (h_sect(3) == "08A" .OR. h_sect(3) == "08B") THEN

            ! MOSES II run, so add increments to TStar_tile:
            ! If snow_depth <= 0.05 (kg/m**2)
            DO tile_num = 1, ntiles
              DO i = 1, land_field

                IF ( snow_depth(land_index(i)) <= 0.05 ) THEN

                  TStar_tile(i,tile_num) = TStar_tile(i,tile_num)       &
                                         + t1_inc(land_index(i))

                ENDIF

              END DO
            END DO

          END IF

        END IF ! (L_IAU_IncTStar)

!----------------------------------------------------------------------
! [1.5]: If required, remove supersaturation wrt water.
!----------------------------------------------------------------------

        IF (L_IAU_RemoveSS) THEN

          IF (PrintStatus >= PrStatus_Normal) THEN
            WRITE (6,*) ''
            WRITE (6,*) 'IAU: Removing supersaturation wrt water.'
          END IF

          DO k = 1, wet_levels

            ! Get temperature and pressure:
            DO j = 1, rows
              DO i = 1, row_length
                ! Use t1_old work array to hold temp on this level:
                t1_old(i,j) = exner_theta_levels(i,j,k) * theta(i,j,k)
                p_tmp (i,j) = p_theta_levels    (i,j,k)
              END DO
            END DO

            ! Calculate saturated specific humidity wrt water:
! DEPENDS ON: qsat_wat
            CALL QSAT_WAT (q_sat_wat, t1_old, p_tmp, row_length*rows)

            ! Limit q to q_sat:
            DO j = 1, rows
              DO i = 1, row_length
                q(i,j,k) = MIN(q(i,j,k), q_sat_wat(i,j))
              END DO
            END DO

          END DO

        END IF ! L_IAU_RemoveSS

        ! Update completed:
        L_FirstUpdate = .FALSE.

      END IF ! (ABS(Weight) < Weight_LL)

!----------------------------------------------------------------------
! [2]: Things to be done at the end of the IAU insertion period.
!----------------------------------------------------------------------

      IF (L_LastCall) THEN

!----------------------------------------------------------------------
! [2.1]: If required, reset stratospheric humidities.
!----------------------------------------------------------------------

        IF (L_IAU_CallStratQ) THEN

          IF (PrintStatus >= PrStatus_Normal) THEN
            WRITE (6,*) ''
            WRITE (6,*) 'IAU: Calling StratQ.'
          END IF

          ! Need to update halos ready for calculation of pv_at_theta:
! DEPENDS ON: swap_bounds
          CALL Swap_Bounds ( u,                                         &
                                                               ! inout
     &                       row_length, rows,   model_levels,          &
                                                               ! in
     &                       offx, offy, fld_type_u, .TRUE. )  ! in

! DEPENDS ON: swap_bounds
          CALL Swap_Bounds ( v,                                         &
                                                               ! inout
     &                       row_length, n_rows, model_levels,          &
                                                               ! in
     &                       offx, offy, fld_type_v, .TRUE. )  ! in

! DEPENDS ON: swap_bounds
          CALL Swap_Bounds ( theta,                                     &
                                                               ! inout
     &                       row_length, rows,   model_levels,          &
                                                               ! in
     &                       offx, offy, fld_type_p, .FALSE. ) ! in

! DEPENDS ON: swap_bounds
          CALL Swap_Bounds ( rho,                                       &
                                                               ! inout
     &                       row_length, rows,   model_levels,          &
                                                               ! in
     &                       offx, offy, fld_type_p, .FALSE. ) ! in

          ALLOCATE ( pv_at_theta(row_length, rows, model_levels) )

! DEPENDS ON: calc_pv_at_theta
          CALL Calc_PV_at_theta ( u, v, theta, rho,                     &
                                                                ! in
     &                            r_theta_levels, r_rho_levels,         &
                                                                ! in
     &                            r_at_u, r_at_v,                       &
                                                                ! in
     &                            sec_v_latitude,                       &
                                                                ! in
     &                            tan_v_latitude,                       &
                                                                ! in
     &                            sec_theta_latitude,                   &
                                                                ! in
     &                            f3_at_v,                              &
                                                                ! in
     &                            delta_lambda, delta_phi,              &
                                                                ! in
     &                            Model_domain,                         &
                                                                ! in
     &                            pv_at_theta )                 ! out

! DEPENDS ON: stratq
          CALL StratQ (                                                 &
#include "argcona.h"
     &                  q,                                              &
                                               ! inout
     &                  qCL,                                            &
                                               ! inout
     &                  qCF,                                            &
                                               ! inout
     &                  area_cloud_fraction,                            &
                                               ! inout
     &                  bulk_cloud_fraction,                            &
                                               ! inout
     &                  cloud_fraction_liquid,                          &
                                               ! inout
     &                  cloud_fraction_frozen,                          &
                                               ! inout
     &                  pv_at_theta,                                    &
                                               ! in
     &                  theta,                                          &
                                               ! in
     &                  p,                                              &
                                               ! in
     &                  p_theta_levels,                                 &
                                               ! in
     &                  exner_theta_levels,                             &
                                               ! in
     &                  IAU_LL_strat_pv,                                &
                                               ! in
     &                  IAU_UL_strat_p,                                 &
                                               ! in
     &                  IAU_LL_trop_p,                                  &
                                               ! in
     &                  L_IAU_Diags )          ! in

           DEALLOCATE (pv_at_theta)

        END IF ! (L_IAU_CallStratQ)

!----------------------------------------------------------------------
! [2.2]: Close IAU increment file and print exit message.
!----------------------------------------------------------------------

! DEPENDS ON: file_close
        CALL FILE_CLOSE (IAU_unit, 'IAU_inc', 7, 0, 0, ICode)
        IF (ICode > 0) THEN
          CMessage = 'Error closing IAU increment file.'
! DEPENDS ON: ereport
          CALL EReport (RoutineName, ICode, CMessage)
        END IF

        IF (PrintStatus >= PrStatus_Normal) THEN
          WRITE (6,*) ''
          WRITE (6,*) 'IAU: IAU completed successfully.'
          WRITE (6,*) ''
        END IF

      END IF ! (L_LastCall)


      RETURN
      END SUBROUTINE IAU

!+ Read in and check a single IAU increment field.




#endif
