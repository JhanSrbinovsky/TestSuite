#if defined (A11_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!*LL  SUBROUTINE TRBDRY-------------------------------------------------
!LL
!LL  Purpose: Special routine to add psuedo source terms to boundary
!LL           data in limited area model.
!LL  Method:  Sets the boundary of aerosol using the
!LL           function BDRYV. This is then copied outward to fill
!LL           the halo. The function BDRYV
!LL           is specific to a model configuration: the current
!LL           version (5.2) is specific to UK MES.
!LL
!LL
!LL Pete Clark  <- programmer of original code
!LL
!LL  Model            Modification history from model version 3.4:
!LL version  Date
!LL  4.2  15/08/96  Add MPP code. Remove unused variables.  RTHBarnes.
!LL  5.2  27/09/00  New Dynamics Changes. 2A version. P.Selwood.
!    5.3  27/07/01  Fix a bug to boundary updating.   Ian Culverwell.
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3,
!LL                        Version 7, dated 11/3/93.
!LL
!LL
!*L  Arguments:---------------------------------------------------------
      Subroutine Trbdry(                                                &
     & row_length, rows, n_rows, model_levels,                          &
     & offx, offy, at_extremity,                                        &
     & pressure,                                                        &
     & u, v,                                                            &
     & murk, timestep                                                   &
     &)

      Implicit None

      Integer, Intent(In) ::                                            &
     &  row_length                                                      &
                          ! Length of a model row
     &, rows                                                            &
                          ! Number of rows for theta,u fields
     &, n_rows                                                          &
                          ! Number of rows for v fields
     &, model_levels                                                    &
                          ! Number of model levels
     &, offx                                                            &
                          ! Size of "single" halo (EW direction)
     &, offy                                                            &
                          ! Size of "single" halo (NS direction)
     &, timestep          ! Timestep in seconds

      Logical, Intent(In) ::                                            &
     &  at_extremity(4)   ! At an edge?

      Real, Intent(In) ::                                               &
     &  u( 1 - offx : row_length + offx,                                &
                                                  ! U wind
     &     1 - offy : rows + offy,                                      &
     &     model_levels )                                               &
     &, v( 1 - offx : row_length + offx,                                &
                                                  ! V wind
     &     1 - offy : n_rows + offy,                                    &
     &     model_levels )                                               &
     &, pressure( 1 - offx : row_length + offx,                         &
                                                  ! Pressure on
     &            1 - offy : rows + offy,                               &
                                                  ! Rho levels
     &            model_levels + 1)

      Real, Intent(InOut) ::                                            &
     &  murk( 1 - offx : row_length + offx,                             &
                                                ! Aerosol for which the
     &        1 - offy : rows + offy,                                   &
                                                ! boundary will be set
     &        model_levels )

#include "c_pi.h"
#include "parparm.h"

!* Local, including SAVE'd, storage------------------------------------

!  (a) Scalars

      Real ::                                                           &
     & dc, wdir, wspeed, press

      Real ::                                                           &
     & BdryV        ! FUNCTION giving boundary value.

!  (b) Others.
      Integer  ::                                                       &
     & i, j, jj, k  ! Loop counters

!-----------------------------------------------------------------------
!L Loop across Northern row.
!-----------------------------------------------------------------------
!
      If ( at_extremity( PNorth ) ) Then
        Do k = 1, model_levels
          j = rows
          jj = n_rows
          Do i = 1, row_length

            If ( v(i, jj, k) > 0. ) Then       ! Outflow
              murk(i, j, k) = murk(i, j-1, k)
            Else                               ! Inflow
              press = pressure(i, j, k)
              wdir = Atan2( v(i, jj, k), u(i, j, k) ) *                 &
     &                                                Recip_Pi_Over_180
              wspeed = Sqrt( u(i, j, k) * u(i, j, k) +                  &
     &                       v(i, jj, k) * v(i, jj, k) )

! DEPENDS ON: bdryv
              murk(i, j, k) = Bdryv( wdir, wspeed, press )
            End If
          End Do
        End Do
      End If



!-----------------------------------------------------------------------
!L Loop across Southern row.
!-----------------------------------------------------------------------
!
      If ( at_extremity( PSouth ) ) Then
        Do k = 1, model_levels
          j = 1
          Do i = 1, row_length

            If ( v(i, j, k) < 0. ) Then             ! Outflow
              murk(i, j, k) = murk(i, j+1, k)
            Else                                    ! Inflow
              press = pressure(i, j, k)
              wdir = Atan2( v(i, j, k), u(i, j, k) ) * Recip_Pi_Over_180
              wspeed = Sqrt( u(i, j, k) * u(i, j, k) +                  &
     &                       v(i, j, k) * v(i, j, k) )

! DEPENDS ON: bdryv
              murk(i, j, k) = Bdryv( wdir, wspeed, press )
            End If
          End Do
        End Do
      End If

!-----------------------------------------------------------------------
!L Loop across Western column
!-----------------------------------------------------------------------
!
      If ( at_extremity( PWest ) ) Then
        Do k = 1, model_levels
          Do j = 1, rows

            ! jj is used for v indexing
            If ( j > n_rows + offy ) Then
              jj = n_rows +offy
            Else
              jj = j
            End If

            i = 1

            If ( u(i,j,k) < 0. ) Then           ! Outflow
              murk(i, j, k) = murk(i+1, j, k)
            Else                                ! Inflow
              press = pressure(i, j, k)
              wdir = Atan2( v(i, jj, k), u(i, j, k) ) *                 &
     &                                                Recip_Pi_Over_180
              wspeed = Sqrt( u(i, j, k) * u(i, j, k) +                  &
     &                       v(i, jj, k) * v(i, jj, k) )

! DEPENDS ON: bdryv
              murk(i, j, k) = Bdryv( wdir, wspeed, press )
            End If
          End Do
        End Do
      End If

!-----------------------------------------------------------------------
!L Loop across Eastern column
!-----------------------------------------------------------------------
!
      If ( at_extremity( PEast ) ) Then
        Do k = 1, model_levels
          Do j = 1, rows

            ! jj is used for v indexing
            If ( j > n_rows + offy ) Then
              jj = n_rows +offy
            Else
              jj = j
            End If

            i = row_length

            If ( u(i,j,k) > 0. ) Then          ! Outflow
              murk(i, j, k) = murk(i-1, j, k)
            Else                               ! Inflow
              press = pressure(i, j, k)
              wdir = Atan2( v(i, jj, k), u(i, j, k) ) *                 &
     &                                                Recip_Pi_Over_180
              wspeed = Sqrt( u(i, j, k) * u(i, j, k) +                  &
     &                       v(i, jj, k) * v(i, jj, k) )

! DEPENDS ON: bdryv
              murk(i, j, k) = Bdryv( wdir, wspeed, press )
            End If
          End Do
        End Do
      End If

!-----------------------------------------------------------------------
! Swap bounds to ensure halos full correctly
!-----------------------------------------------------------------------

! DEPENDS ON: swap_bounds
      Call Swap_Bounds( murk, row_length, rows, model_levels,           &
     &                  offx, offy, fld_type_p, .false. )

!-----------------------------------------------------------------------
!L Now need to fill full extended halos with copies of the
!L calculated data. Start at the North
!-----------------------------------------------------------------------

      If ( at_extremity( PNorth ) ) Then
        Do k = 1, model_levels
          Do j = rows + 1, rows + offy
            Do i = 1 - offx, row_length + offx
              murk(i, j, k) = murk(i, rows, k)
            End Do
          End Do
        End Do
      End If

!-----------------------------------------------------------------------
! Fill extended halos for Southern rows
!-----------------------------------------------------------------------

      If ( at_extremity( PSouth ) ) Then
        Do k = 1, model_levels
          Do j = 1 - offy, 0
            Do i = 1 - offx, row_length + offx
              murk(i, j, k) = murk(i, 1, k)
            End Do
          End Do
        End Do
      End If

!-----------------------------------------------------------------------
! Fill extended halos for Western columns
!-----------------------------------------------------------------------

      If ( at_extremity( PWest ) ) Then
        Do k = 1, model_levels
          Do j = 1 - offy, rows + offy
            Do i = 1 - offx, 0
              murk(i, j, k) = murk(1, j, k)
            End Do
          End Do
        End Do
      End If

!-----------------------------------------------------------------------
! Fill extended halos for Eastern columns
!-----------------------------------------------------------------------

      If ( at_extremity( PEast ) ) Then
        Do k = 1, model_levels
          Do j = 1 - offy, rows + offy
            Do i = row_length + 1, row_length + offx
              murk(i, j, k) = murk (row_length, j, k)
            End Do
          End Do
        End Do
      End If

      Return
      End Subroutine Trbdry


!*LL  FUNCTION BDRYV----------------------------------------------------
!LL
!LL  Purpose: Special routine to add psuedo source terms to boundary
!LL           data in limited area.
!LL  PARAMETERS ARE SPECIFIC TO UK MESOSCALE MODEL
!LL  Method:  The boundary concentrations are computed using a
!LL           simple model of transport from sources outside the
!LL           model. Analysis of the source distribution outside
!LL           the UK MES shows that it can be well represented by
!LL           a line source at constant radius from the centre of
!LL           the model, with a source distribution given by the
!LL           sum of two Gaussians. Concentrations from these are
!LL           computed assuming transport using the local windspeed u
!LL           or 1 m/s, whichever is stronger, over a distance
!LL           determined from the centroid of the source distribution, x
!LL           with a linear transformation rate k from emission to
!LL           aerosol, dry deposition at a rate determined from the
!LL           dry deposition velocity vd and mean mixed layer depth h.
!LL           Thus the max concentration is given by
!LL               Q/(uh)*k/(k+vd/h)*(1-exp(-k*x/u))
!LL           The source term is assumed to decrease with level
!LL           pressure. See forthcoming documentation for details.
!LL
!LL Pete Clark  <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.4:
!LL version  Date
!    6.2  26/05/06  Add missing comma. P.Selwood
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3,
!LL                        Version 7, dated 11/3/93.
!LL
!LL
!*L  Arguments:---------------------------------------------------------
#endif
