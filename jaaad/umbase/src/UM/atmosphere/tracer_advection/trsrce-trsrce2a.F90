#if defined (A11_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!*LL  SUBROUTINE TRSRCE ----------------------------------------------
!LL
!LL  Purpose: Adds source increments to a single level of the aerosol
!LL  field.
!LL
!LL  Suitable for single-column use.
!LL
!LL Pete Clark  <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.4:
!LL version  Date
!LL   4.0  05/09/95 Model hour and minute added to call to TRSRCE.
!LL                 Diurnal variation of emissions, hard wired to
!LL                 peak at midday, with 10% amplitude, assumed.
!LL                 Programmer Pete Clark.
!LL   4.1  02/05/96 Amplitude of diurnal variation added to input
!LL                 argument list.             M.Woodage
!LL   5.2  02/10/00 New dynamics changes. 2A version.  P.Selwood
!  5.3  05/10/01  Correct halo sizes for q.       M Woodage, A Malcolm
!     6.2  11/08/05 Replace quotes with angle-brackets in includes.
!                   P.Selwood
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3,
!LL                        Version 7, dated 11/3/93.
!LL
!LL
!*L  Arguments:---------------------------------------------------------
      SUBROUTINE TRSRCE(                                                &
     & rows, row_length, offx, offy, halo_i, halo_j,                    &
     & model_levels, wet_levels,                                        &
     & q_halo_i, q_halo_j,                                              &
     & r_rho_levels, r_theta_levels,                                    &
     & theta, q, qcl, qcf, exner, rho, tracer, srce,                    &
     & level, timestep, i_hour, i_minute, amp                           &
     &)

      IMPLICIT NONE

      Integer, Intent(In)    :: rows         ! number of P/U rows
      Integer, Intent(In)    :: row_length   !
      Integer, Intent(In)    :: model_levels ! number of levels for P
      Integer, Intent(In)    :: wet_levels   ! number of wet-levels.
      Integer, Intent(In)    :: offx         ! EW size of std. halo
      Integer, Intent(In)    :: offy         ! NS size of std. halo
      Integer, Intent(In)    :: halo_i       ! EW extended halo
      Integer, Intent(In)    :: halo_j       ! NS extended halo
      Integer, Intent(In)    :: q_halo_i     ! EW extended halo for q
      Integer, Intent(In)    :: q_halo_j     ! NS extended halo for q
      Integer, Intent(In)    :: i_hour       ! Local time hour
      Integer, Intent(In)    :: i_minute     ! Local time minute
      Integer, Intent(In)    :: level        ! level of the tracer

      Real, Intent(In)       :: timestep     ! Timestep in seconds
      Real, Intent(In)       :: amp          ! Amplitude of diurnal
                                             ! variation of emission

      Real, Intent(In)       ::                                         &
     &      r_rho_levels  ( 1 - halo_i : row_length + halo_i,           &
     &                      1 - halo_j : rows + halo_j,                 &
     &                      model_levels )                              &
     &,     r_theta_levels( 1 - halo_i : row_length + halo_i,           &
     &                      1 - halo_j : rows + halo_j,                 &
     &                      0 : model_levels )

      Real, Intent(In)       ::                                         &
     &      theta( 1 - offx : row_length + offx,                        &
                                                     ! pot. temperature
     &             1 - offy : rows + offy,                              &
                                                     !
     &             model_levels )                                       &
     &,     q    ( 1-q_halo_i :row_length+q_halo_i,                     &
                                                        ! Q on theta
     &             1-q_halo_j :rows+q_halo_j,                           &
                                                        ! levels
     &             wet_levels )                                         &
     &,     qcl  ( 1-q_halo_i :row_length+q_halo_i,                     &
                                                        ! Qcl on theta
     &             1-q_halo_j :rows+q_halo_j,                           &
                                                        ! levels
     &             wet_levels )                                         &
     &,     qcf  ( 1-q_halo_i :row_length+q_halo_i,                     &
                                                        ! Qcf on theta
     &             1-q_halo_j :rows+q_halo_j,                           &
                                                        ! levels
     &             wet_levels )                                         &
     &,     exner( 1 - offx : row_length + offx,                        &
                                                   ! exner on rho
     &             1 - offy : rows + offy,                              &
                                                   ! levels
     &             model_levels + 1)                                    &
     &,     rho  ( 1 - offx : row_length + offx,                        &
                                                   ! density * r * r
     &             1 - offy : rows + offy,                              &
                                                   ! on rho levels
     &             model_levels )

      Real, Intent(In)       ::                                         &
     &      srce( 1 : row_length, 1 : rows)           ! tracer source

      Real, Intent(InOut)    ::                                         &
     &      tracer( 1 - offx : row_length + offx,                       &
                                                  ! level of tracer
     &              1 - offy : rows + offy )      ! to be updated

#include "c_g.h"
#include "c_pi.h"
#include "c_epslon.h"
#include "c_r_cp.h"

!L  Workspace usage----------------------------------------------------
!L  External subroutine called ----------------------------------------
      External Ereport

! Local, including SAVE'd, storage------------------------------------

!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
      Real             :: DM
      Real             :: TS
      Real             :: thetav        ! virtual potential temperature
      Real             :: exner_ave     ! an averaged exner term
      Real             :: rho_theta     ! rho on theta level
      Real             :: rho1          ! } values of rho after the
      Real             :: rho2          ! } r-squared factor removed

      Real, Parameter  :: Factor = 1.0  ! Factor to multiply source term
      Real, Parameter  :: TZero  = 12.0 ! Time of maximum emissions

!  (b) Others.
      Integer          :: i   ! Loop counter
      Integer          :: j   ! Loop counter

! Error Reporting
      Integer                      :: ErrorStatus
      Character (Len=*), Parameter :: RoutineName='TRSRCE'
      Character (Len=80)           :: Cmessage

!-----------------------------------------------------------------------
!  Check level number is not too large.
!-----------------------------------------------------------------------
      If (level > model_levels ) Then
        ErrorStatus = 10
        Write (Cmessage,*) 'Level for tracer updating is larger ',      &
     &                     ' than model_levels'
! DEPENDS ON: ereport
        Call Ereport( RoutineName, ErrorStatus, Cmessage )
      End If

!-----------------------------------------------------------------------
!L Subroutine structure :
!L Loop over field adding source term*timestep/level mass/unit area.
!-----------------------------------------------------------------------

!     Allow for source varying with time of day
      TS = 1.0
      If (AMP > 0.0) Then
        TS = 1.0 +                                                      &
     &  amp * Cos( (Real(i_hour) + Real(i_minute)/60.0 - TZero)         &
     &      * PI/12.0)
      End If

      TS = TS * timestep * factor

      Do j=1, rows
        Do i = 1, row_length
          If (level < model_levels) Then
! Remove the r squared factor from rho before interpolation
           rho1 =  rho(i,j,level)/(r_rho_levels(i,j,level) *            &
     &                             r_rho_levels(i,j,level) )
           rho2 =  rho(i,j,level+1)/(r_rho_levels(i,j,level+1) *        &
     &                               r_rho_levels(i,j,level+1) )


! DM = density (interpolated on to theta levels) * delta r

           DM = rho2 * (r_theta_levels(i,j,level) -                     &
     &                  r_rho_levels(i,j,level) ) +                     &
     &          rho1 * (r_rho_levels(i,j,level+1) -                     &
     &                  r_theta_levels(i,j,level) )
!
! Special case for lowest layer to get correct mass
           If (level == 1) Then
             DM = DM * (r_rho_levels(i,j,2) - r_theta_levels(i,j,0))    &
     &               / (r_rho_levels(i,j,2) - r_rho_levels(i,j,1))
           End If
!
! Convert DM to DRY density if level is wet
           If (level <= wet_levels) Then
             DM = DM * (1.0 -q(i,j,level)-qcl(i,j,level)-qcf(i,j,level))
           End If
!
          Else    ! level = model_level
!--------------------------------------------------------------------
! Cannot average here to get rho_theta. Hence calculate by using
! the equation of state vis
!         ___r                     P0
!         rho   = -----------------------------------------
!                 kappa * Cp * ____________________r * theta_v
!                              [          kappa-1 ]
!                              [ exner ** ------- ]
!                              [          kappa   ]
!-------------------------------------------------------------------

           If (wet_levels == model_levels) Then
             thetav = theta(i,j,level) * (1.0 +                         &
     &                  ( q(i,j,level) * C_Virtual )                    &
     &                - qcl(i,j,level)-qcf(i,j,level)  )
           Else
             thetav = theta(i,j,level)
           End If

           exner_ave = (exner(i,j,level) ** ((kappa - 1.0)/kappa) +     &
     &                  exner(i,j,level+1) ** ((kappa - 1.0)/kappa))/2.0
           rho_theta = Pref/( kappa * Cp * exner_ave * thetav )

! rho_theta is at the top theta level. We also need the value
! at the top rho level. This will be rho1
           rho1 = rho(i,j,model_levels)/(r_rho_levels(i,j,model_levels)*&
     &                                   r_rho_levels(i,j,model_levels))

! rho2 will be the average of rho1 and rho_theta
           rho2 = ( rho1 + rho_theta ) * 0.5

           DM = rho2 * (r_theta_levels(i,j,level) -                     &
     &                  r_rho_levels(i,j,level) )
          End If

          tracer(i, j) = tracer(i, j) + srce(i, j) * TS/DM
        End Do
      End Do

      Return
      END SUBROUTINE TRSRCE
#endif
