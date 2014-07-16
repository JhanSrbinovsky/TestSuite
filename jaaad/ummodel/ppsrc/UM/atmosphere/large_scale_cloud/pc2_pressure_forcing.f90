
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ PC2 Cloud Scheme: Forcing due to advection (adiabatic cooling)
! Subroutine Interface:
      SUBROUTINE PC2_PRESSURE_FORCING(                                  &
!      Parallel variables
     & halo_i, halo_j, off_x, off_y                                     &
!      Pressure related fields
     &,P,PSTAR,p_theta_levels                                           &
!      Array dimensions
     &,wet_levels, row_length, rows, rhc_row_length, rhc_rows           &
!      timestep and rhcrit values
     &,timestep, rhcpt                                                  &
!      Prognostic Fields
     &,THETA, CF, CFL, CFF, Q, QCL, QCF                                 &
!      Forcing quantities for driving the homogeneous forcing
     &,EXNER_DEPARTURE, EXNER_CURRENT                                   &
!      Convective cloud base and convection flag
     &,ccb, cumulus, rhts, tlts, qtts, ptts, CF_AREA                    &
!      Diagnostics
     &, t_inc, q_inc, qcl_inc, qcf_inc, cf_inc, cfl_inc, cff_inc        &
     &, t_dini,q_dini,qcl_dini,qcf_dini,cf_dini,cfl_dini,cff_dini       &
!      Logical switches
     &, l_mixing_ratio, L_ACF_CUSACK                                    &
     & )
!
      IMPLICIT NONE
!
! Purpose:
!   This subroutine interfaces to the homogeneous forcing routine and
!   calculates the temperature forcing associated with adiabatic
!   changes in pressure across a timestep. This is mainly due to
!   large-scale vertical advection. It also checks values of water
!   contents and cloud fractions to check that they are sensible.
!
! Method:
!   Uses the departure and end values of exner to calculate the
!   pressure and, from theta, the temperature and temperature
!   forcing. Calls the homogeneous forcing routine to update values
!   of temperature, moisture and cloud, and converts the temperature
!   change to a theta increment.
!
! Current Owner of Code: D. R. Wilson
!
! History:
! Version   Date     Comment
!  5.4    22-07-02   Original Code (Damian Wilson)
!  5.5    16-12-02   Bring up to date with current PC2 formulation
!                                                        D. Wilson
!  6.1    25-05-04   Pass convective information into
!                    pc2_initiation_ctl                  D. Wilson
!  6.2    22-12-05   Include dummy arguments for passing down for
!                    SCM diagnostics               A. Kerr-Munslow
!  6.2    31-03-05   Pass RH(TL) to pc2_initiation_ctl   D. Wilson
!  6.4    18-08-06   Use mixing ratio formulation.  Damian Wilson
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: PC2 Cloud Scheme Documentation
!
!  Global Variables:----------------------------------------------------
!
!  Subroutine Arguments:------------------------------------------------
      INTEGER                                                           &
                        !, INTENT(IN)
     & halo_i, halo_j                                                   &
!       Size of large halos
     &,off_x, off_y                                                     &
!       Size of small halos
     &,wet_levels                                                       &
!       Number of levels being processed.
     &,row_length, rows                                                 &
!       Row length and number of rows being processed.
     &,rhc_row_length, rhc_rows                                         &
!       Row length and number of rows in the rhcpt variable
     &,ccb(row_length,rows)
!       convective cloud base
!
       LOGICAL                                                          &
                        !, INTENT(IN)
     & cumulus(row_length,rows)                                         &
     &,L_ACF_CUSACK                                                     &
     &,l_mixing_ratio  ! Use mixing ratio formulation
!
      REAL                                                              &
                        !, INTENT(IN)
     & p_theta_levels(1-off_x:row_length+off_x, 1-off_y:rows+off_y,     &
     &                wet_levels)                                       &
!       pressure at all points (Pa)
     &,p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                  &
     &    wet_levels+1)                                                 &
     &,pstar(row_length, rows)                                          &
     &,timestep                                                         &
!       Model timestep (s)
     &,EXNER_DEPARTURE(1-off_x:row_length+off_x, 1-off_y:rows+off_y,    &
     &               wet_levels)                                        &
!       Departure value of exner (before the advection)
     &,EXNER_CURRENT(1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &               wet_levels)                                        &
!       Current value of exner
     &,rhcpt(rhc_row_length,rhc_rows,wet_levels)                        &
!       Values of critical relative humidity
     &,rhts(row_length,rows,wet_levels)                                 &
!       Relative total humidity at start of timestep (time level n)
     &,tlts(row_length,rows,wet_levels)                                 &
!       TL at start of timestep
     &,qtts(row_length,rows,wet_levels)                                 &
!       qT at start of timestep
     &,ptts(row_length,rows,wet_levels)
!       pressure at theta levels at start of timestep
!
      REAL                                                              &
                        !, INTENT(INOUT)
     & THETA(1-off_x:row_length+off_x, 1-off_y:rows+off_y, wet_levels)  &
!       Potential temperature (K)
     &,CF(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j, wet_levels) &
!       Total cloud fraction (no units)
     &,CFL(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j, wet_levels)&
!       Liquid cloud fraction (no units)
     &,CFF(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j, wet_levels)&
!       Ice cloud fraction (no units)
     &,Q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j, wet_levels)  &
!       Vapour content (kg water per kg air)
     &,QCL(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j, wet_levels)&
!       Liquid content (kg water per kg air)
     &,QCF(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j, wet_levels)&
!       Ice content (kg water per kg air)
     &,CF_AREA(row_length, rows, wet_levels) !area cloud fraction
!
      REAL                                                              &
                        !, INTENT(OUT)
     & T_INC(row_length,rows,wet_levels)                                &
!       Change in temperature due to CONDENSATION (K)
     &,Q_INC(row_length,rows,wet_levels)                                &
!       Change in vapour due to CONDENSATION (kg kg-1)
     &,QCL_INC(row_length,rows,wet_levels)                              &
!       Change in liquid due to CONDENSATION (kg kg-1)
     &,QCF_INC(row_length,rows,wet_levels)                              &
!       Change in ice due to CONDENSATION (kg kg-1)
     &,CF_INC(row_length,rows,wet_levels)                               &
!       Change in total cloud fraction due to CONDENSATION (no units)
     &,CFL_INC(row_length,rows,wet_levels)                              &
!       Change in liquid cloud fraction due to CONDENSATION (no units)
     &,CFF_INC(row_length,rows,wet_levels)
!       Change in ice cloud fraction due to CONDENSATION (no units)
!
      REAL                                                              &
                        !, INTENT(OUT)
     & T_dini(row_length,rows,wet_levels)                               &
!       Change in temperature due to initiation CONDENSATION (K)
     &,Q_dini(row_length,rows,wet_levels)                               &
!       Change in vapour due to initiation CONDENSATION (kg kg-1)
     &,QCL_dini(row_length,rows,wet_levels)                             &
!       Change in liquid due to initiation CONDENSATION (kg kg-1)
     &,QCF_dini(row_length,rows,wet_levels)                             &
!       Change in ice due to initiation CONDENSATION (kg kg-1)
     &,CF_dini(row_length,rows,wet_levels)                              &
!       Change in total cloud fraction due to initiation CONDENSATION
     &,CFL_dini(row_length,rows,wet_levels)                             &
!       Change in liquid cloud fraction due to initiation CONDENSATION
     &,CFF_dini(row_length,rows,wet_levels)
!       Change in ice cloud fraction due to initiation CONDENSATION

!
!  External functions:
!
!  Local parameters and other physical constants------------------------
!
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
      REAL                                                              &
     & one_over_kappa
      PARAMETER ( one_over_kappa = 1.0 / kappa )
!
!  Local scalars--------------------------------------------------------
!
      REAL                                                              &
     & T_DEPARTURE                                                      &
                          ! Temperature at the departure point (K)
     &,P_DEPARTURE        ! Pressure at the departure point (Pa)
!
      INTEGER K,I,J       ! Loop counters: K - vertical level index
!                           I,J - horizontal position index
!
!  Local dynamic arrays-------------------------------------------------
      REAL                                                              &
     & T(row_length,rows,wet_levels)                                    &
!       Temperature (K)
     &,DELTAT(row_length,rows,wet_levels)                               &
!       Change in temperature between departure and arrival points (K)
     &,DELTAP(row_length,rows,wet_levels)                               &
!       Change in pressure between departure and arrival points (Pa)
     &,ZEROS(row_length,rows,wet_levels)                                &
!       Array of zero values
!
!  Local arrays which do not contain halos
     &,p_theta_levels_no_halos(row_length,rows,wet_levels)              &
     &,cf_no_halos(row_length,rows,wet_levels)                          &
     &,cfl_no_halos(row_length,rows,wet_levels)                         &
     &,cff_no_halos(row_length,rows,wet_levels)                         &
     &,q_no_halos(row_length,rows,wet_levels)                           &
     &,qcl_no_halos(row_length,rows,wet_levels)                         &
     &,qcf_no_halos(row_length,rows,wet_levels)
!
!  Local dummy variables for SCM Diagnostics to be passed to
!  PC2_initiation_ctl
      INTEGER, PARAMETER :: nSCMDpkgs = 12

      LOGICAL                                                           &
     & L_SCMDiags(nSCMDpkgs)
!
!- End of Header
!
! ==Main Block==--------------------------------------------------------
!
! ----------------------------------------------------------------------
! 1. Calculate temperature at departure and current (arrival) points
! ----------------------------------------------------------------------
!
      DO k=1,wet_levels
        DO j=1,rows
          DO i=1,row_length
            T_DEPARTURE   = THETA(i,j,k) * EXNER_DEPARTURE(i,j,k)
            T(i,j,k)      = THETA(i,j,k) * EXNER_CURRENT(i,j,k)
            DELTAT(i,j,k) = T(i,j,k)     - T_DEPARTURE
            P_DEPARTURE   = EXNER_DEPARTURE(i,j,k)**one_over_kappa *pref
            DELTAP(i,j,k) = p_theta_levels(i,j,k) - P_DEPARTURE
            ZEROS(i,j,k)  = 0.0
! Copy primary variables into non-halo variables
            p_theta_levels_no_halos(i,j,k) = p_theta_levels(i,j,k)
            cf_no_halos(i,j,k)  = cf(i,j,k)
            cfl_no_halos(i,j,k) = cfl(i,j,k)
            cff_no_halos(i,j,k) = cff(i,j,k)
            q_no_halos(i,j,k)   = q(i,j,k)
            qcl_no_halos(i,j,k) = qcl(i,j,k)
            qcf_no_halos(i,j,k) = qcf(i,j,k)
! Copy variables into increment variables for diagnostics
            t_inc(i,j,k)        = t(i,j,k)
            q_inc(i,j,k)        = q(i,j,k)
            qcl_inc(i,j,k)      = qcl(i,j,k)
            qcf_inc(i,j,k)      = qcf(i,j,k)
            cf_inc(i,j,k)       = cf(i,j,k)
            cfl_inc(i,j,k)      = cfl(i,j,k)
            cff_inc(i,j,k)      = cff(i,j,k)
          END DO
        END DO
      END DO
!
! ----------------------------------------------------------------------
! 2. Call homogeneous forcing routine with just the temperature forcing
! ----------------------------------------------------------------------
!
! DEPENDS ON: pc2_homog_plus_turb
      CALL PC2_HOMOG_PLUS_TURB(p_theta_levels_no_halos, wet_levels      &
     &,                        row_length, rows, timestep, T            &
     &,                        cf_no_halos, cfl_no_halos, cff_no_halos  &
     &,                        q_no_halos, qcl_no_halos                 &
     &,                        deltat, zeros, zeros, deltap, 0.0, 0.0   &
     &,                        l_mixing_ratio)
!
! ----------------------------------------------------------------------
! 3. Remember that the homogeneous forcing routine will include the
!    forcing increments in its changes. Deltat has already been included
!    in the model, so we need to subtract it again.
! ----------------------------------------------------------------------
!
      DO k=1,wet_levels
        DO j=1,rows
          DO i=1,row_length
            t(i,j,k) = t(i,j,k) - deltat(i,j,k)
          END DO
        END DO
      END DO
!
! ----------------------------------------------------------------------
! 4. Call initiation routine
! ----------------------------------------------------------------------
!
! Set dummy variables for SCM Diagnostics
      L_SCMDiags(1:nSCMDpkgs) = .FALSE.
!
! NB if you are changing the argument list to PC2_initiation_ctl
! please do an equivalent change in routine scm_main to keep the
! Single Column Model consistent.
!
! DEPENDS ON: pc2_initiation_ctl
      Call pc2_initiation_ctl(                                          &
     &  halo_i, halo_j, off_x, off_y                                    &
     &, row_length, rows, rhc_row_length, rhc_rows                      &
     &, wet_levels, wet_levels, .false., l_mixing_ratio                 &
     &, L_ACF_CUSACK,timestep                                           &
     &, nSCMDpkgs, L_SCMDiags                                           &
     &, T, q_no_halos, qcl_no_halos, qcf_no_halos                       &
     &, cf_no_halos,cfl_no_halos,cff_no_halos,rhts                      &
     &, tlts,qtts,ptts,CF_AREA,p,pstar                                  &
     &, p_theta_levels_no_halos, ccb, cumulus                           &
     &, rhcpt,t_dini,q_dini,qcl_dini,qcf_dini,cf_dini,cfl_dini          &
     &, cff_dini)
!
! ----------------------------------------------------------------------
! 5. Update variables
! ----------------------------------------------------------------------
!
      DO k=1,wet_levels
        DO j=1,rows
          DO i=1,row_length
! Convert change of temperature to change of theta
            THETA(i,j,k) = T(i,j,k) / EXNER_CURRENT(i,j,k)
! Copy changed non-halo variables into primary variables
            cf(i,j,k)  = cf_no_halos(i,j,k)
            cfl(i,j,k) = cfl_no_halos(i,j,k)
            cff(i,j,k) = cff_no_halos(i,j,k)
            q(i,j,k)   = q_no_halos(i,j,k)
            qcl(i,j,k) = qcl_no_halos(i,j,k)
            qcf(i,j,k) = qcf_no_halos(i,j,k)
! Increment in temperature due to the CONDENSATION associated with the
! change in pressure. This is not the same as the temperature change
! due to vertical advection, which is in delta_t
            T_inc(i,j,k)   = T(i,j,k)   - t_inc(i,j,k)  -  t_dini(i,j,k)
! Other increment variables
            q_inc(i,j,k)   = q(i,j,k)   - q_inc(i,j,k)  -  q_dini(i,j,k)
            qcl_inc(i,j,k) = qcl(i,j,k) - qcl_inc(i,j,k)-qcl_dini(i,j,k)
            qcf_inc(i,j,k) = qcf(i,j,k) - qcf_inc(i,j,k)-qcf_dini(i,j,k)
            cf_inc(i,j,k)  = cf(i,j,k)  - cf_inc(i,j,k) - cf_dini(i,j,k)
            cfl_inc(i,j,k) = cfl(i,j,k) - cfl_inc(i,j,k)-cfl_dini(i,j,k)
            cff_inc(i,j,k) = cff(i,j,k) - cff_inc(i,j,k)-cff_dini(i,j,k)
          END DO
        END DO
      END DO
!
! End of the subroutine
!
      RETURN
      END SUBROUTINE PC2_PRESSURE_FORCING
