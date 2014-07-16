#if defined(A04_3D)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Autoconversion of graupel.
! Subroutine Interface:
      SUBROUTINE LSP_GRAUP_AUTOC(                                       &
     &  points, timestep                                                &
                                          ! Number of points and tstep
     &, qcf_agg, qgraup, T, rho                                         &
                                          ! Water contents and temp
     &, auto_graup_qcf_thresh, auto_graup_T_thresh, auto_graup_coeff    &
                                          ! Parametrization information
     &, l_seq                                                           &
                                          ! Code options
     &, psacw, psdep, ptransfer, iterations                             &
                                          ! Mass transfer diagnostic
     &  )
!
      Implicit None
!
! Purpose:
!   Update cloud prognostics as a result of the autoconversion of
!   snow aggregates to graupel.
!
!  Method:
!   This is a source term for graupel when snow growth is dominated
!   by riming liquid water cloud and is the only term that can create
!   graupel where there was no graupel before.
!   The conversion rate is proportional to how much the riming rate
!   of aggregates (PSACW) exceeds the rate of growth due to vapour
!   deposition (PSDEP) and collection/accretion of ice crystals (PSACI).
!   The coefficient reduces the conversion rate as the riming
!   aggregates will not immediately increase their density to that of
!   graupel.
!
! Current Owner of Code: Jonathan Wilkinson
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP 26.
!
!
! Subroutine Arguments
!
#include "c_0_dg_c.h"
!
      Integer, Intent(In) ::                                            &
     &  points                                                          &
                          ! Number of points to process
     &, iterations
                          ! Number of microphysics iterations
!
      Real, Intent(In) ::                                               &
     &  Timestep                                                        &
                          ! Timestep / s
     &, auto_graup_qcf_thresh                                           &
                          ! Threshold ice content for the formation
                          ! of graupel
     &, auto_graup_T_thresh                                             &
                          ! Threshold temperature for graupel
                          ! formation / deg C  
     &, auto_graup_coeff                                                &
                          ! Relationship between riming rate and
                          ! graupel production / no units
     &, T(points)                                                       &
                          ! Temperature / K
     &, rho(points)                                                     &
                          ! Air density / kg kg s-1
     &, psacw(points)                                                   &
                          ! Riming rate of snow aggs. / kg kg-1 s-1
     &, psdep(points)
                          ! Deposition rate of snow aggs. / kg kg-1 s-1
!
      Real, Intent(InOut) ::                                            &
     &  qgraup(points)                                                  &
                          ! Graupel content / kg kg-1
     &, qcf_agg(points)                                                 &
                          ! Ice water content of snow aggs. / kg kg-1
     &, ptransfer(points)
                          ! Autoconversion rate / kg kg-1 s-1
!
      Logical, Intent(In) ::                                            &
     &  l_seq
                          ! Carry out sequential updating
!
! Local Variables
!
      Integer                                                           &
     &  i

!
      Real                                                              &
     &  dqi(points)
                          ! Transfer amount from ice to snow / kg kg-1

      Do i=1, points

        If (rho(i)*qcf_agg(i)  >   auto_graup_qcf_thresh                &
     &      .and. T(i)  <  (zerodegc+auto_graup_T_thresh)) Then

          !-----------------------------------------------
          ! Calculate conversion rate / kg kg-1 s-1
          !-----------------------------------------------
          dqi(i) = auto_graup_coeff *                                   &
       &           max(0.0, psacw(i)-psdep(i))

          dqi(i) = min(dqi(i)*timestep,qcf_agg(i))

          !-----------------------------------------------
          ! Store process rate / kg kg-1 s-1
          !-----------------------------------------------
          ptransfer(i) = ptransfer(i) + dqi(i)/(timestep*iterations)

          !-----------------------------------------------
          ! Recalculate water contents
          !-----------------------------------------------
          If (l_seq) then
            qgraup(i)  = qgraup(i)  + dqi(i)
            qcf_agg(i) = qcf_agg(i) - dqi(i)
          End if  ! l_seq

          ! No cloud fraction updates

        End if !  qcf_agg > threshold etc.

      End Do  ! Points

      Return  ! End of the subroutine
      END SUBROUTINE LSP_GRAUP_AUTOC
#endif
