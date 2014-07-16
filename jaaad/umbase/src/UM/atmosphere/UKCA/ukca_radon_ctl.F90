#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  To manage Rn-222/Pb tracer experiments. Includes surface emissions
!  and radioactive decay.
!
!  Part of the UKCA model. UKCA is a community model supported
!  by The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Current Code Owner:       Colin Johnson/Olaf Morgenstern
!                           Fiona O'Connor
!
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
! Subroutine Interface:
      SUBROUTINE RADON_CTL(                                            &
!
! Parallel variables
        halo_i, halo_j, off_x, off_y, global_row_length, global_rows   &
      , proc_row_group, proc_col_group, at_extremity, n_proc, n_procx  &
      , n_procy, neighbour, g_rows, g_row_length, g_datastart, me      &
! model dimensions
      , row_length, rows, n_rows, land_points                          &
      , model_levels, wet_model_levels                                 &
      , tr_model_levels, bl_levels, n_cca_levels                       &
! Model switches
      , model_domain, L_CAL360, L_SEC_VAR, L_EqT, Ltimer               &
! Model parameters
! Physical constants
      , lc, lf, cp, two_Omega, p_zero, kappa                           &
      , R, g, Lapse_Rate, earth_radius, Pi                             &
! Co-ordinate information
      , r_rho_levels, r_theta_levels                                   &
      , eta_theta_levels, eta_rho_levels                               &
      , delta_lambda, delta_phi                                        &
      , lat_rot_NP, long_rot_NP                                        &
! Time stepping information
      , timestep                                                       &
      , val_year, val_day_number, val_hour, val_minute                 &
      , val_second, timestep_number                                    &
      , PREVIOUS_TIME                                                  &
      , CALL_CHEM_FREQ                                                 &
! Trig arrays
! Grid-dependent arrays
      , f3_at_u, true_longitude                                        &
!
! Data fields IN
      , tr_radon_em                                                    &
      , theta, q, qcl, qcf                                             &
      , rho,                                                           &
!     &  land_mask,
        p_theta_levels, exner_rho_levels, exner_theta_levels           &
! Logicals IN
      , L_Radon_surfem                                                 &
!
! Data fields IN/OUT
      , tr_radon, tr_lead                                              &
!
       )
!
      IMPLICIT NONE
!
! Arguments with intent IN:
!
! Parallel setup variables
      Integer                                                          &
        halo_i                                                         &
                   ! Size of halo in i direction.
      , halo_j                                                         &
                   ! Size of halo in j direction.
      , off_x                                                          &
                   ! Size of small halo in i
      , off_y                                                          &
                   ! Size of small halo in j.
      , global_row_length                                              &
                           ! number of points on a row
      , proc_row_group                                                 &
                       ! Group id for processors on the same row
      , proc_col_group                                                 &
                       ! Group id for processors on the same col
      , global_rows                                                    &
                           ! NUMBER OF global rows
      , n_proc                                                         &
                   ! Total number of processors
      , n_procx                                                        &
                   ! Number of processors in longitude
      , n_procy                                                        &
                   ! Number of processors in latitude
      , neighbour(4)                                                   &
                             ! Array with the Ids of the four neighbour
                             ! in the horizontal plane
      , g_rows (0:n_proc-1)                                            &
      , g_row_length (0:n_proc-1)                                      &
      , g_datastart (3,0:n_proc-1)                                     &
      , me         ! My processor number
!
      Logical                                                          &
        at_extremity(4) ! Indicates if this processor is at north, sout
                        ! east or west of the processor grid
! Model dimensions
      Integer                                                          &
        row_length                                                     &
      , rows                                                           &
      , n_rows                                                         &
      , land_points                                                    &
      , model_levels                                                   &
      , wet_model_levels                                               &
      , tr_model_levels                                                &
      , bl_levels                                                      &
      , n_cca_levels    ! No. conv cloud levels (1 if 2D, nlevs if 3D)
! Model switches
      Integer                                                          &
        model_domain
      Logical                                                          &
        L_CAL360                                                       &
                        ! T if using 360 day calendar
      , L_SEC_VAR                                                      & 
                        ! if T include secular varn of earth's orbit,
      , L_EqT                                                          &
                        ! if T include equn of time         in SOLPOS
      , Ltimer                                                           
                        ! if T then output some timing information

! Physical constants
      Real                                                             &
        lc, lf, cp                                                     &
      , two_Omega                                                      & 
                        ! twice Earth's rotation rate
      , p_zero                                                         &
      , kappa                                                          &
      , R, g, Lapse_Rate, earth_radius, Pi
! Co-ordinate arrays
      Real                                                             &
        r_theta_levels(1-halo_i:row_length+halo_i,                     &
                         1-halo_j:rows+halo_j,0:model_levels)          &
      , r_rho_levels(1-halo_i:row_length+halo_i,                       &
                       1-halo_j:rows+halo_j, model_levels)             &
      , eta_theta_levels(0:model_levels)                               &
      , eta_rho_levels(model_levels)                                   &
      , delta_lambda                                                   &
      , delta_phi                                                      &
      , lat_rot_NP                                                     &
      , long_rot_NP
! Trig arrays
! Grid-dependent arrays
      Real                                                             &
        f3_at_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y)         &
      , true_longitude(row_length, rows)                               &
      , FV_cos_theta_latitude (1-off_x:row_length+off_x,               &
                               1-off_y:rows+off_y)
! Time stepping information
      Real                                                             &
        timestep                       !atmosphere model timetsep
      Integer                                                          &
        val_year                                                       &
      , val_day_number                                                 &
      , val_hour                                                       &
      , val_minute                                                     &
      , val_second                                                     &
      , timestep_number                                                &
      , PREVIOUS_TIME(7)                                               &
      , CALL_CHEM_FREQ                  !frequency of calling chemistry
!                                         per atmos phys timestep
!
#include "csubmodl.h"
#include "typsts.h"
!
! Diagnostics info
!      Real
!        STASHwork17(*)  ! STASH workspace for section 17 (Aero_Ctl)
!
      Real                                                             &
        theta(1-off_x:row_length+off_x,1-off_y:rows+off_y,             &
                                                       model_levels)   &
      , q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
                                                   wet_model_levels)   &
      , qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
                                                   wet_model_levels)   &
      , qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
                                                   wet_model_levels)   &
      , rho(1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
                                                       model_levels)   &
      , p_theta_levels(1-off_x:row_length+off_x,                       &
                       1-off_y:rows+off_y, model_levels)               &
      , exner_rho_levels(1-off_x:row_length+off_x,                     &
                                 1-off_y:rows+off_y, model_levels+1)   &
      , exner_theta_levels(1-off_x:row_length+off_x,                   &
                                   1-off_y:rows+off_y, model_levels)   &
      ,tr_radon_em(row_length, rows)
                                         !radon emissions
!
      LOGICAL                                                          &
       LAND_MASK(row_length,rows)                                      &
                                         !T IF LAND, F IF SEA
      ,L_Radon_surfem                    !T if RADON surface emissions
!
! Arguments with intent IN/OUT:
      REAL                                                             &
       tr_radon(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
                                      tr_model_levels)                 &
                                                         !mmr radon
      ,tr_lead(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
                                      tr_model_levels)   !mmr lead

! Local variables
!
      INTEGER I,J,K           ! loop variables
      INTEGER FIRST_POINT,                                             &
                              ! First and last pts on which calcns done
              LAST_POINT      !  (omits N and S polar rows)
      INTEGER                                                          &
              NPNTS          ! no. of pts in 3D array on P LEVS


!      Integer i,j,k,n
        INTEGER N
!
      REAL CHEMSTEP                         ! chemistry timestep
!
        real max_em,min_em

! Calculate length of chemistry timestep and use it to control input
! of emissions and Rn Chemistry
!
      CHEMSTEP=timestep/(CALL_CHEM_FREQ*1.)
!
      DO n=1,CALL_CHEM_FREQ
!
! Call TRSRCE to insert emissions
!
        IF (L_RadoN_SURFEM) THEN        ! Insert surface radon emiss

          WRITE(6,*) 'MAX/MIN ems = ',maxval(tr_radon_em),             &
                      minval(tr_radon_em)
          WRITE(6,*) 'MAX/MIN tr_radon = ',maxval(tr_radon(:,:,1)),    &
                      minval(tr_radon(:,:,1))

          WRITE(6,*) 'Calling TRSRCE, halo_i/j=',halo_i,halo_j

          IF (tr_model_levels /= wet_model_levels)                     &
            write(6,*)'RADON_CTL: tr/wet_model_levels=',               &
            tr_model_levels,wet_model_levels

! DEPENDS ON: trsrce
            CALL TRSRCE(                                               &
              rows, row_length, off_x, off_y, halo_i, halo_j,          &
              model_levels, wet_model_levels,                          &
              halo_i, halo_j,                                          &
              r_rho_levels, r_theta_levels,                            &
              theta, q, qcl, qcf, exner_rho_levels, rho,               &
              tr_radon(:,:,1), tr_radon_em,                            &
              1, CHEMSTEP, 1, 1, 0.0,                                  &
              )

        ENDIF        ! end of IF (L_radon_surfem) statement

! DEPENDS ON: radon_decay
        CALL RADON_DECAY(                                              &
          off_x, row_length, off_y, rows, tr_model_levels,             &
          tr_Radon,tr_lead,CHEMSTEP)

      ENDDO            ! End of CALL_CHEM_FREQ loop
!
      RETURN
      END SUBROUTINE RADON_CTL
#endif
