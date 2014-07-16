#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Reset stratospheric humidities.

      SUBROUTINE StratQ (                                               &
#include "argcona.h"
     &                    q,                                            &
                                                 ! inout
     &                    qCL,                                          &
                                                 ! inout
     &                    qCF,                                          &
                                                 ! inout
     &                    area_cloud_fraction,                          &
                                                 ! inout
     &                    bulk_cloud_fraction,                          &
                                                 ! inout
     &                    cloud_fraction_liquid,                        &
                                                 ! inout
     &                    cloud_fraction_frozen,                        &
                                                 ! inout
     &                    pv_at_theta,                                  &
                                                 ! in
     &                    theta,                                        &
                                                 ! in
     &                    p,                                            &
                                                 ! in
     &                    p_theta_levels,                               &
                                                 ! in
     &                    exner_theta_levels,                           &
                                                 ! in
     &                    LL_strat_pv,                                  &
                                                 ! in
     &                    UL_strat_p,                                   &
                                                 ! in
     &                    LL_trop_p,                                    &
                                                 ! in
     &                    L_Diags )              ! in

      Use level_heights_mod, Only  :  r_theta_levels
      Use trignometric_mod,  Only  :  cos_theta_latitude

! Description:
!
!   Apply climatalogical bounds to stratospheric specific humidities:
!
!     1.0E-06 <= q <= 3.0E-06, RH < 10%
!
!   For consistency, qCL and qCF are set to zero in the stratosphere.
!
!   Stratospheric points are diagnosed using potential vorticity (PV)
!   and pressure. q points with ABS(PV) values greater or equal to a
!   supplied lower-limit, and pressure values less than or equal to a
!   supplied upper-limit, are assumed to be in the stratosphere. Also
!   supplied is a lower-limit for tropospheric pressures. q points with
!   pressure values lower than this limit are assumed to be in the
!   stratosphere.
!
!   Note that q, qCL and qCF haloes are not updated.
!
!
! Current Code Owner: Adam Clayton.
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 5.2       14/12/00 Original code based on Stuart Bell's pre-5.0
!                    version. Adam Clayton
! 5.4       13/11/02 Re-dimension bulk, liquid & frozen cloud fraction
!                    fields with extended haloes. Dave Robinson
! 5.5       09/01/03 Initialisation of q_new and change in RH test.
!                    Jorge Bornemann
! 6.2       21/10/05 Replaced GSYNC with SSYNC. P.Selwood
! 6.2       23/02/06 Increase q upper limit for consistency with
!                    methane oxidation scheme. David Jackson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! Declarations:

! Common blocks:

      IMPLICIT NONE

#include "parvars.h"
#include "typsize.h"
#include "cmaxsize.h"
#include "typcona.h"
#include "c_g.h"
#include "cntlatm.h"

! Subroutine arguments:

      REAL, INTENT(IN) ::                                               &
      
     &  pv_at_theta        ( 1          : row_length,                   &
     &                       1          : rows,                         &
     &                       1          : model_levels ),               &
      
     &  theta              ( 1 - offx   : row_length + offx,            &
     &                       1 - offy   : rows       + offy,            &
     &                       1          : model_levels ),               &
      
     &  p                  ( 1 - offx   : row_length + offx,            &
     &                       1 - offy   : rows       + offy,            &
     &                       1          : model_levels + 1 ),           &
      
     &  p_theta_levels     ( 1 - offx   : row_length + offx,            &
     &                       1 - offy   : rows       + offy,            &
     &                       1          : model_levels ),               &
      
     &  exner_theta_levels ( 1 - offx   : row_length + offx,            &
     &                       1 - offy   : rows       + offy,            &
     &                       1          : model_levels ),               &
      
     &  LL_strat_pv,                                                    &
                     ! Lower-limit for stratospheric ABS(PV)
     &  UL_strat_p,                                                     &
                     ! Upper-limit for stratospheric pressure
     &  LL_trop_p    ! Lower-limit for tropospheric  pressure

      REAL, INTENT(INOUT) ::                                            &
      
     &  q   ( 1 - halo_i : row_length + halo_i,                         &
     &        1 - halo_j : rows       + halo_j,                         &
     &        1          : wet_levels ),                                &
      
     &  qCL ( 1 - halo_i : row_length + halo_i,                         &
     &        1 - halo_j : rows       + halo_j,                         &
     &        1          : wet_levels ),                                &
      
     &  qCF ( 1 - halo_i : row_length + halo_i,                         &
     &        1 - halo_j : rows       + halo_j,                         &
     &        1          : wet_levels ),                                &
      
     &  area_cloud_fraction   ( row_length, rows, wet_levels ),         &
      
     &  bulk_cloud_fraction ( 1 - halo_i : row_length + halo_i,         &
     &                        1 - halo_j : rows       + halo_j,         &
     &                        1          : wet_levels ),                &
      
     &  cloud_fraction_liquid ( 1 - halo_i : row_length + halo_i,       &
     &                          1 - halo_j : rows       + halo_j,       &
     &                          1          : wet_levels ),              &
      
     &  cloud_fraction_frozen ( 1 - halo_i : row_length + halo_i,       &
     &                          1 - halo_j : rows       + halo_j,       &
     &                          1          : wet_levels )

      LOGICAL, INTENT(IN) ::                                            &
      
     &  L_Diags      ! If set, calculate and write out diagnostics.

! Local variables:

      INTEGER :: i, j, k, ICode

      LOGICAL :: Strat(row_length,rows)

      REAL :: t    (row_length,rows),                                   &
     &        p_tmp(row_length,rows),                                   &
     &        q_sat(row_length,rows),                                   &
     &        q_new(row_length,rows)

      REAL :: gb_mass,                                                  &
                            ! Total mass in gridbox
     &        strat_mass,                                               &
                            ! Total mass in stratosphere
     &        total_mass,                                               &
                            ! Total mass
     &        vmass_before,                                             &
                            ! Vapour mass on entry to   StratQ
     &        vmass_added,                                              &
                            ! Vapour mass added by      StratQ
     &        vmass_subtr,                                              &
                            ! Vapour mass subtracted by StratQ
     &        strat_frac    ! Fraction of total mass in stratosphere

      REAL :: LocalStats(5,wet_levels)

      REAL :: q_LL, q_UL, RH_UL, q_UL_used

      PARAMETER ( q_LL  = 1.0E-06 )
      PARAMETER ( q_UL  = 3.0E-06 )
      PARAMETER ( RH_UL = 0.1     )

! External subroutines called:

      EXTERNAL QSAT, GCG_RSUMR

!- End of header ------------------------------------------------------

! IF L_UPPER_STRATQ=T, then  increase the q upper limit to 3.75mg/kg
! in order to match the maximum value allowed by the
! methane oxidation scheme
      IF (L_UPPER_STRATQ) THEN
        q_UL_used = q_UL + 0.75E-06
      ELSE
        q_UL_used = q_UL
      ENDIF
      IF (L_Diags) THEN
        WRITE (6,*) ''
        WRITE (6,'(A,F6.3,A)')                                          &
     &              ' StratQ: Lower-limit for strat PV:       ',        &
     &              LL_strat_pv * 1.0E06, ' PVU'
        WRITE (6,'(A,F6.1,A)')                                          &
     &              '         Upper-limit for strat pressure: ',        &
     &              UL_strat_p / 100.0, ' hPa'
        WRITE (6,'(A,F6.1,A)')                                          &
     &              '         Lower-limit for trop  pressure: ',        &
     &              LL_trop_p  / 100.0, ' hPa'
        WRITE (6,*) ''
        WRITE (6,*) '        Layer-by-layer diags (masses in kg):'
        WRITE (6,*) ''
        WRITE (6,*) '                 Fraction of    Stratospheric  '   &
     &                             //'Vapour mass  Vapour mass'
        WRITE (6,*) '                 total mass in  vapour mass    '   &
     &                             //'added by     subtracted'
        WRITE (6,*) '        q layer  stratosphere   before StratQ  '   &
     &                             //'StratQ       by StratQ'

        WRITE (6,*) '        -------  -------------  -------------  '   &
     &                             //'-----------  -----------'
      END IF

      DO k = 1, wet_levels

        ! Get temperature and pressure.
        DO j = 1, rows
          DO i = 1, row_length
            t    (i,j) = exner_theta_levels(i,j,k) * theta(i,j,k)
            p_tmp(i,j) = p_theta_levels    (i,j,k)
          END DO
        END DO

        ! Calculate saturated specific humidity.
! DEPENDS ON: qsat
        CALL QSAT ( q_sat, t, p_tmp, row_length*rows )

        ! Initialise q_new:
        q_new(:,:) = q(1:row_length,1:rows,k)

        ! Apply bounds in stratosphere.
        DO j = 1, rows
          DO i = 1, row_length

            Strat(i,j) =                                                &
     &        ( (ABS(pv_at_theta(i,j,k)) >= LL_strat_pv) .AND.          &
     &          (p_theta_levels (i,j,k)  <= UL_strat_p)                 &
     &        ) .OR.                                                    &
     &        (p_theta_levels(i,j,k) < LL_trop_p)

            IF (Strat(i,j)) THEN
              qCL(i,j,k) = 0.0
              qCF(i,j,k) = 0.0
              cloud_fraction_liquid(i,j,k) = 0.0
              cloud_fraction_frozen(i,j,k) = 0.0
              area_cloud_fraction  (i,j,k) = 0.0
              bulk_cloud_fraction  (i,j,k) = 0.0
              IF (q(i,j,k) < q_LL) q_new(i,j) = q_LL
              IF (q(i,j,k) > q_UL_used) q_new(i,j) = q_UL_used

              IF (q_new(i,j)/q_sat(i,j) > RH_UL)                        &
     &          q_new(i,j) = q_sat(i,j) * RH_UL
            END IF

          END DO ! i
        END DO ! j

        ! If required, calculate local statistics.
        IF (L_Diags) THEN

          ! Initialise stats.
          strat_mass   = 0.0
          total_mass   = 0.0
          vmass_before = 0.0
          vmass_added  = 0.0
          vmass_subtr  = 0.0

          DO j = 1, rows
            DO i = 1, row_length

              ! Approximate mass in gridbox.
              gb_mass = ( p(i,j,k) - p(i,j,k+1) )                       &
     &                  * cos_theta_latitude(i,j)                       &
     &                  * delta_lambda * delta_phi                      &
     &                  * r_theta_levels(i,j,k)**2                      &
     &                  / g

              total_mass = total_mass + gb_mass

              IF (Strat(i,j)) THEN
                strat_mass   = strat_mass   + gb_mass
                vmass_before = vmass_before + gb_mass * q(i,j,k)
                IF (q_new(i,j) < q(i,j,k)) THEN
                  vmass_subtr = vmass_subtr                             &
     &                        + gb_mass * (q(i,j,k) - q_new(i,j))
                ELSE IF (q_new(i,j) > q(i,j,k)) THEN
                  vmass_added = vmass_added                             &
     &                        + gb_mass * (q_new(i,j) - q(i,j,k))
                END IF
              END IF

            END DO
          END DO

          LocalStats(1,k) = strat_mass
          LocalStats(2,k) = total_mass
          LocalStats(3,k) = vmass_before
          LocalStats(4,k) = vmass_added
          LocalStats(5,k) = vmass_subtr

        END IF ! (L_Diags)

        ! Update q.
        q(1:row_length,1:rows,k) = q_new(:,:)

      END DO ! k

      ! If required, get global stats and write them out.
      IF (L_Diags) THEN

        CALL GC_SSYNC ( nproc, ICode )

        CALL GCG_RSUMR ( wet_levels*5,                                  &
                                            ! in
     &                   gc_all_proc_group,                             &
                                            ! in
     &                   ICode,                                         &
                                            ! out
     &                   LocalStats )       ! inout

        DO k = wet_levels, 1, -1
          strat_frac   = LocalStats(1,k) / LocalStats(2,k)
          vmass_before = LocalStats(3,k)
          vmass_added  = LocalStats(4,k)
          vmass_subtr  = LocalStats(5,k)
          WRITE (6,'(A,I3,4(A,E12.5))')                                 &
     &      '          ', k, '    ',                                    &
     &      strat_frac,  '   ', vmass_before, '   ',                    &
     &      vmass_added, ' ', vmass_subtr
        END DO

        WRITE (6,*) ''

      END IF ! (L_Diags)


      RETURN
      END SUBROUTINE StratQ
#endif
