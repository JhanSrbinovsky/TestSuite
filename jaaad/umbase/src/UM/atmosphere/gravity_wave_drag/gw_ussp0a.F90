#if defined(A06_0A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Stub/Dummy Gravity Wave (USSP) Scheme.
! Subroutine Interface:
      SUBROUTINE GW_USSP(LEVELS, MODEL_DOMAIN, ROWS, NROWS,             &
     &  OFF_X, OFF_Y, HALO_I, HALO_J, ROW_LENGTH,                       &
     &  R_RHO_LEVELS, R_THETA_LEVELS, P_LAYER_BOUNDARIES,               &
     &  R_U, R_V, L_USSP_OPAQUE,                                        &
     &  SIN_THETA_LONGITUDE, SIN_THETA_LATITUDE,                        &
     &  THETA, RHO, TIMESTEP, U, V, AT_EXTREMITY,                       &
     &  GWSPEC_EFLUX,GWSPEC_SFLUX,GWSPEC_WFLUX,GWSPEC_NFLUX,            &
     &  GWSPEC_EWACC,GWSPEC_NSACC,ETA_THETA_LEVELS,                     &
     &  GWSPEC_EFLUX_ON, GWSPEC_SFLUX_ON, GWSPEC_WFLUX_ON,              &
     &  GWSPEC_NFLUX_ON, GWSPEC_EWACC_ON, GWSPEC_NSACC_ON)
!
! purpose: Dummy routine called when 4A GWD scheme is switched off
!
      IMPLICIT NONE

! ----------------------------------------------------------------------+-------
!     Subroutine arguments of GW_USSP
! ----------------------------------------------------------------------+-------
      INTEGER                                                           &
        LEVELS                                                          &
                             !IN Number of model levels
      , ROWS                                                            &
                             !IN Number of rows for u field
      , NROWS                                                           &
                             !IN Number of rows for v field
      , OFF_X                                                           &
                             !IN offset longitude
      , OFF_Y                                                           &
                             !IN offset latitude
      , HALO_I                                                          &
                             !IN Halo in longitude
      , HALO_J                                                          &
                             !IN Halo in latitude
      , ROW_LENGTH                                                      &
                             !IN Number of grid points in row
      , MODEL_DOMAIN         !IN Model type (global, LAM etc)
      REAL                                                              &
        SIN_THETA_LONGITUDE(ROW_LENGTH,ROWS)                            &
                                              !Grid point longitudes
      , SIN_THETA_LATITUDE(ROW_LENGTH,ROWS)                             &
                                              !P-GRID Latitudes
      , GWSPEC_EFLUX(ROW_LENGTH,ROWS,LEVELS)                            &
                                              ! Fp in each of 4
      , GWSPEC_SFLUX(ROW_LENGTH,NROWS,LEVELS)                           &
                                              ! azimuths for diags.
      , GWSPEC_WFLUX(ROW_LENGTH,ROWS,LEVELS)                            &
                                              !
      , GWSPEC_NFLUX(ROW_LENGTH,NROWS,LEVELS)                           &
                                              !
      , GWSPEC_EWACC(ROW_LENGTH,ROWS,LEVELS)                            &
                                              !
      , GWSPEC_NSACC(ROW_LENGTH,NROWS,LEVELS)                           &
                                              !
      , ETA_THETA_LEVELS(0:LEVELS)                                      &
      , THETA(ROW_LENGTH,ROWS,LEVELS)                                   &
                             !IN Primary model array for theta
      , RHO(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y,LEVELS)         &
                             !IN Primary model array for Density
                             !x (radius earth)^2.
      , TIMESTEP                                                        &
                             !IN Timestep
      , U(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y,LEVELS)           &
                             !INOUT Primary model array for U field
      , V(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:NROWS+OFF_Y,LEVELS)          &
                             !INOUT Primary model array for V field
      , R_RHO_LEVELS(1-HALO_I:ROW_LENGTH+HALO_I,1-HALO_J:ROWS           &
                     +HALO_J,LEVELS)                                    &
                   !Distance of rho levels from Earth centre.
      , R_THETA_LEVELS(1-HALO_I:ROW_LENGTH+HALO_I,1-HALO_J:             &
                       ROWS+HALO_J,0:LEVELS)                            &
                   !Distance of theta levels from Earth centre.
      , P_LAYER_BOUNDARIES(ROW_LENGTH,ROWS,0:LEVELS)                    &
      , R_U(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:ROWS+OFF_Y,LEVELS)         &
      , R_V(1-OFF_X:ROW_LENGTH+OFF_X,1-OFF_Y:NROWS+OFF_Y,LEVELS)
      LOGICAL                                                           &
        L_USSP_OPAQUE   !IN  Switch for Opaque Upper boundary condition

      LOGICAL                                                           &
        GWSPEC_EFLUX_ON                                                 &
                            !Switches for diagnostics of
      , GWSPEC_SFLUX_ON                                                 &
                            !Fp in each of 4 azimuths
      , GWSPEC_WFLUX_ON                                                 &
                            !
      , GWSPEC_NFLUX_ON                                                 &
                            !
      , GWSPEC_EWACC_ON                                                 &
                            !and accelerations
      , GWSPEC_NSACC_ON                                                 &
                            !
      , AT_EXTREMITY(4)     !Edge of domain indicator                   

      INTEGER                       ::  ICODE
      CHARACTER (Len=80)            ::  CMESSAGE
      CHARACTER (Len=* ), Parameter ::  RoutineName='GW_USSP_0A'

      CMESSAGE = 'Routine should not be callable'
      ICODE = 1
! DEPENDS ON: ereport
      CALL EReport(RoutineName,ICODE,CMESSAGE)

      RETURN
!
      END SUBROUTINE GW_USSP
#endif
