
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Cloud Scheme: Forcing due to advection (adiabatic cooling)
! Subroutine Interface:
      SUBROUTINE PC2_ASSIM(                                             &
     & levels, row_length, rows                                         &
                                                  ! Array dimensions
     &,timestep                                                         &
                                                  ! Timestep
     &,L_pc2_homog, l_pc2_cff, l_mixing_ratio                           &
                                                  ! Control logicals
     &,T, cf, cfl, cff, q, qcl, qcf, p                                  &
                                                  ! Prognostic Fields
     &,deltat, deltaq, deltaqcl, deltaqcf, deltap                       &
                                                  ! Forcing quantities
!                                                 ! No diagnostics yet
     & )
!
      IMPLICIT NONE
!
! Purpose:
!   This subroutine acts to provide the PC2 cloud fraction and
!   condensation increments interfacting to the data assimilation
!   AC and VAR schemes.
!
! Method:
!   Uses the homogenous forcing routine to calculate condensation
!   and associated cloud fraction changes and performs an estimate
!   of ice cloud fraction changes based on the ice water content
!   forcing. Parts of the code can be switched with logicals to allow
!   future data assimilation schemes only to use the parts that they
!   need.
!
! Current Owner of Code: D. R. Wilson, PC2 team
!
! History:
! Version   Date     Comment
!  6.1    25-05-03   Original Code for UM vn6.1 (Damian Wilson)
!  6.4    18-08-06   Use mixing ratio formulation.  Damian Wilson
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: New Cloud Scheme Documentation
!
!  Global Variables:----------------------------------------------------
!
!  Subroutine Arguments:------------------------------------------------
      LOGICAL                                                           &
                        !, INTENT(IN)
     & L_pc2_homog, L_pc2_cff                                           &
!      Control on which sections of code need using
     &,L_mixing_ratio  ! Use mixing ratio formulation
!
      INTEGER                                                           &
                        !, INTENT(IN)
     & levels                                                           &
!       Number of levels being processed.
     &,row_length, rows
!       Row length and number of rows being processed.
!
      REAL                                                              &
                        !, INTENT(IN)
     & timestep                                                         &
!       Model timestep (s)
     &,deltat(row_length,rows,levels)                                   &
                                         ! Forcing of temperature over
!                                        ! the timestep            (K)
     &,deltaq  (row_length,rows,levels)                                 &
                                         ! Forcing of vapour (kg kg-1)
     &,deltaqcl(row_length,rows,levels)                                 &
                                         ! Forcing of liquid (kg kg-1)
     &,deltaqcf(row_length,rows,levels)                                 &
                                         ! Forcing of ice    (kg kg-1)
     &,deltap  (row_length,rows,levels)  ! Forcing of pressure    (Pa)
!
      REAL                                                              &
                        !, INTENT(IN)
     & qcf(row_length,rows,levels)                                      &
                                       ! Ice content          (kg kg-1)
     &,p(row_length,rows,levels)       ! Pressure at theta levels  (Pa)
!
      REAL                                                              &
                        !, INTENT(INOUT)
     & T(row_length,rows,levels)                                        &
                                       ! Dry-bulb temperature (K)
     &,cf(row_length,rows,levels)                                       &
                                       ! Bulk cloud fraction  (no units)
     &,cfl(row_length,rows,levels)                                      &
                                       ! Liquid cloud frac    (no units)
     &,cff(row_length,rows,levels)                                      &
                                       ! Ice cloud fraction   (no units)
     &,q(row_length,rows,levels)                                        &
                                       ! Vapour content       (kg kg-1)
     &,qcl(row_length,rows,levels)     ! Liquid content       (kg kg-1)

!  External subroutine calls: ------------------------------------------
      EXTERNAL pc2_homog_plus_turb
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
!
!  Local parameters-----------------------------------------------------
!
      REAL   ,PARAMETER:: qcf0         =1.0E-4 ! Specified in-cloud
!                                                 ice content (kg kg-1)
!
!  Local scalars--------------------------------------------------------
!
      INTEGER points                                                    &
                                               ! Counter for qcf changes
     &,       i,j,k                            ! loop counters
!
!  Local dynamic arrays-------------------------------------------------
      REAL                                                              &
     &        deltacff_c(row_length*rows*levels)                        &
                                                 ! Change in ice cloud
!                                                   fraction
     &,       cf_c (row_length*rows*levels)                             &
                                                ! Condensed points
     &,       cfl_c(row_length*rows*levels)                             &
                                                !   copies of cloud
     &,       cff_c(row_length*rows*levels)                             &
                                                !   fraction variables
     &,       zeros(row_length*rows*levels)     ! Array of zeros
!
      REAL                                                              &
     &        T_copy(row_length,rows,levels)                            &
                                                ! Copy of T
     &,       q_copy(row_length,rows,levels)    ! Copy of q
!
!- End of Header
!
! ==Main Block==--------------------------------------------------------
!
! ----------------------------------------------------------------------
! 1. Call homogenous forcing routine if needed
! ----------------------------------------------------------------------
      If (L_pc2_homog) Then
!
! In order that temperature (t) and vapour (q) are not updated by
! the condensation response to assimilation, we take copies of these
! variables to use in the forcing subroutine. The subroutine will
! update T_copy etc., not T. The cloud fractions and qcl *are* updated.
!
        DO k=1,levels
          DO j=1,rows
            DO i=1,row_length
                T_copy(i,j,k) = t(i,j,k)
                q_copy(i,j,k) = q(i,j,k)
            End Do
          End Do
        End Do
!
! DEPENDS ON: pc2_homog_plus_turb
      CALL PC2_HOMOG_PLUS_TURB(p, levels, row_length, rows,timestep     &
     &,                        T_copy, cf, cfl, cff, q_copy, qcl        &
     &,                        deltat, deltaq, deltaqcl, deltap         &
     &,                        0.0, 0.0, l_mixing_ratio)
!
! Now update q and T for the forcings
!
        DO k=1,levels
          DO j=1,rows
            DO i=1,row_length
                T(i,j,k) = T(i,j,k) + deltaT(i,j,k)
                q(i,j,k) = q(i,j,k) + deltaq(i,j,k)
            End Do
          End Do
        End Do
!
      End If  ! L_pc2_homog
!
      If (L_pc2_cff) Then
!
! ----------------------------------------------------------------------
! 2. Calculate changes to ice cloud fraction
! ----------------------------------------------------------------------
!
! Estimate the in-cloud ice content of the part of the cloud that is
! being added or removed to be qcf_ic = cf * (qcf/cf) + (1-cf) * qcf0
! where qcf0 is a specified in-cloud content.
!
! Loop over each point
!
        points=0
        DO k=1,levels
          DO j=1,rows
            DO i=1,row_length
              If (deltaqcf(i,j,k)  /=  0.0) then
                points=points+1
                ! Store condensed points versions
                cf_c (points)=cf (i,j,k)
                cfl_c(points)=cfl(i,j,k)
                cff_c(points)=cff(i,j,k)
                ! Adjust cff for each point
                deltacff_c(points) = deltaqcf(i,j,k)                    &
     &                       / ( qcf(i,j,k) + (1.0-cff(i,j,k))*qcf0 )
                ! Check that cff is sensible
                deltacff_c(points) =                                    &
     &            max(  min(  deltacff_c(points), (1.0 - cff(i,j,k))  ) &
     &                ,-cff(i,j,k))
              ! Form zero array for forcing of liquid cloud fraction.
              ! Note that pc2_homog_plus_turb calls pc2_total_cf for the
              ! liquid changes
              zeros(points) = 0.0
              End If
            End Do
          End Do
        End Do
!
! ----------------------------------------------------------------------
! 3. Calculate changes to total cloud fraction to update cf_c
! ----------------------------------------------------------------------
!
! DEPENDS ON: pc2_total_cf
        CALL PC2_TOTAL_CF(                                              &
     &        Points,cfl_c,cff_c,zeros,deltacff_c,cf_c)
!
! ----------------------------------------------------------------------
! 4. Scatter back changes to ice and total cloud fractions
! ----------------------------------------------------------------------
!
        points=0
        DO k=1,levels
          DO j=1,rows
            DO i=1,row_length
              If (deltaqcf(i,j,k)  /=  0.0) then
                points=points+1
                cff(i,j,k) = cff(i,j,k) + deltacff_c(points)
                cf (i,j,k) = cf_c(points)
              End If
            End Do
          End Do
        End Do
!
      End If  ! L_pc2_cff
!
! ----------------------------------------------------------------------
! 5. Call pc2_checks to make sure values are self consistent
! ----------------------------------------------------------------------
!
      If (L_pc2_homog .or. L_pc2_cff) Then
!
! DEPENDS ON: pc2_checks
        CALL PC2_CHECKS(p, levels,                                      &
     &      row_length, rows, t, cf, cfl, cff, q, qcl, qcf,             &
     &      l_mixing_ratio)
!
      End If  ! L_pc2_homog .or. L_pc2_cff
!
! End of the subroutine
!
 9999 CONTINUE ! Error exit
      RETURN
      END SUBROUTINE PC2_ASSIM
