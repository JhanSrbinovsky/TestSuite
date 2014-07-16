
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Autoconversion of snow.
! Subroutine Interface:
      SUBROUTINE LSP_SNOW_AUTOC(                                        &
     &  points, timestep                                                &
                                          ! Number of points and tstep
     &, qcf_cry, qcf_agg, T, CTTemp                                     &
                                          ! Water contents and temp
     &, m0, T_scaling, qcf0                                             &
                                          ! Parametrization information
     &, l_seq                                                           &
                                          ! Code options
     &, ptransfer, iterations                                           &
                                          ! Mass transfer diagnostic
     &  )
!
      Implicit None
!
! Purpose:
!   Update cloud prognostics as a result of the autoconversion of
!   cloud ice to snow aggregates
!
!  Method:
!   Transfer some mass from ice crystals to snow aggregates depending
!   on the temperature and cloud top temperature.
!   Simple explicit Kessler type param. of autoconversion
!   with the autoconversion limit and rate set to emulate the
!   split-ice scheme.
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
      Integer, Intent(In) ::                                            &
     &  points                                                          &
                          ! Number of points to process
     &, iterations
                          ! Number of microphysics iterations
!
      Real, Intent(In) ::                                               &
     &  Timestep                                                        &
                          ! Timestep / s
     &, m0                                                              &
                          ! Seed ice mass / kg kg-1
     &, T_scaling                                                       &
                          ! Scaling temperature / K
     &, qcf0                                                            &
                          ! Prescribed ice content / kg kg-1
     &, T(points)                                                       &
                          ! Temperature / K
     &, CTTemp(points)
                          ! Cloud top temperature / K
!
      Real, Intent(InOut) ::                                            &
     &  qcf_cry(points)                                                 &
                          ! Ice water content of ice crystals / kg kg-1
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
     &  dpr(points)                                                     &
                          ! Transfer amount from ice to snow / kg kg-1
     &, qcfautolim(points)                                              &
                          ! Autoconversion limit / kg kg-1
     &, qcfautorate(points)                                             &
                          ! Rate of transfer / s-1
     &, qc(points)
                          ! Ice remaining after autoconversion

      Do i=1, points

        If (qcf_cry(i) > m0) Then

          !-----------------------------------------------
          ! Set autoconversion limit to emulate split-ice scheme
          !-----------------------------------------------
          qcfautolim(i) = (qcf_agg(i)+qcf_cry(i))                       &
     &               *max(exp(-T_scaling * max((T(I)-CTTemp(i)),0.0)    &
     &               *max(qcf_agg(i)+qcf_cry(i),0.0)*qcf0) , 0.0)

          !-----------------------------------------------
          ! Set rate to emulate spilt-ice scheme, i.e. infinite
          !-----------------------------------------------
          qcfautorate(i) = 1.0/timestep

          qc(i)  = min(qcfautolim(i) , qcf_cry(i))
          dpr(i) = min(qcfautorate(i) * timestep * (qcf_cry(i)-qc(i))   &
     &                , qcf_cry(i) - qc(i))

          !-----------------------------------------------
          ! Store process rate (kg kg-1 s-1)
          !-----------------------------------------------
          ptransfer(i) = ptransfer(i) + dpr(i)/(timestep*iterations)

          !-----------------------------------------------
          ! Update ice/snow variables
          !-----------------------------------------------
          If (l_seq) then
            qcf_cry(i) = qcf_cry(i)-dpr(i)
            qcf_agg(i) = qcf_agg(i)+dpr(i)
          End if  ! l_seq

          ! No cloud fraction updating is needed

       End if  ! qcf_cry > 0

      End Do  ! Points

      Return  ! End of the subroutine
      END SUBROUTINE LSP_SNOW_AUTOC
