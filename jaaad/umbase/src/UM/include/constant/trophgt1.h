#if defined(ATMOS)
#if defined(CONTROL) || defined(SCMA)
!+ ---------------------------------------------------------------------
!  Description: Limits on the position of the tropopause used in
!               radiative calculations are defined.
!
!  Current Code Owner: J. M. Edwards
!
!  History:
!
!  Version  Date      Comment.
!  5.3      27/09/01  Original code.
!                     (J. M. Edwards)
!  6.2      27/01/06  Allows limits to be set in SCM (R. Wong)
!
!- ----------------------------------------------------------------------
!
!     Set limits on the height of the tropopause. Generous limits are
!     required as the tropopause can be very low in the Antarctic
!     winter. These limits should be reviewed for simulations of
!     climates very different from the present one, such as runs
!     with very high concentrations of CO2.
!
!     These limits are in part chosen to accord with the pressure
!     levels used in earlier configurations of the Unified Model.
!     In setting the upper limit consideration has been given to the
!     paper entitled "The tropical tropopause over the Western
!     Pacifiic: Wave driving, convection and the annual cycle,"
!     (J. Geophys. Res., 1996, Vol. 101 (D16), p. 21223) by
!     G. C. Reid and K. S. Gage.
      Real, parameter :: z_min_trop = 2.0e3 ! Lowest permitted height
      Real, parameter :: z_max_trop = 2.0e4 ! Maximum permitted height
!
! -----------------------------------------------------------------------
#endif

#endif
