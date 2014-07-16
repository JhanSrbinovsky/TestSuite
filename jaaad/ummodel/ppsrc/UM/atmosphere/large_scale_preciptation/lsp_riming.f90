
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Riming of ice particles
! Subroutine Interface:
      SUBROUTINE LSP_RIMING(                                            &
     &  points, timestep                                                &
                                          ! Number of points and tstep
     &, qcl, qcf, T                                                     &
                                          ! Water contents and temp
     &, area_liq, area_mix, cfliq, cficei                               &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
!    &, cf, cfl, cff                      ! Current cloud fractions for
!                                         ! updating
     &, rho, m0, tcg, tcgi, corr                                        &
                                          ! Parametrization information
     &, cx, constp, ai, bi, aic, bic, lfrcp , ice_type                  &
                                          ! Microphysical information
     &, l_update_cf, l_seq, l_psd                                       &
                                          ! Code options
     &, ptransfer, iterations                                           &
                                          ! Mass transfer diagnostic
!    &, cftransfer,cfltransfer,cfftransfer! Cloud transfer diagnostics
     &  )
!
      Implicit None
!
! Purpose:
!   Update cloud prognostics as a result of ice particle riming
!
! Method:
!   Solve the microphysical transfer equation for a specified
!   distribution of ice particles sweeping out liquid water.
!
! Current Owner of Code: Jonathan Wilkinson
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP 26.
!
! Riming is a source term for ice crystals
! and a sink term for supercooled water.
! Occurs if there are ice crystals and liquid water drops present
! in a mixed phase overlap region. The Numerical solution is implicit
! in liquid water.
! Riming does not, in this formulation, update cloud fractions
! which is why these variables have been commented out.
!
! Subroutine Arguments
!
! Obtain the size for cx and constp and define ice_type_offset
! Start C_LSPSIZ
! Description: Include file containing idealised forcing options
! Author:      R. Forbes
!
! History:
! Version  Date      Comment
! -------  ----      -------
!   6.1    01/08/04  Increase dimension for rain/graupel.  R.Forbes
!   6.2    22/08/05  Include the step size between ice categories.
!                                                   Damian Wilson

! Sets up the size of arrays for CX and CONSTP
      REAL CX(100),CONSTP(100)
      INTEGER,PARAMETER:: ice_type_offset=20

! End C_LSPSIZ
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
!
      Integer, Intent(In) ::                                            &
     &  points                                                          &
                          ! Number of points to calculate
     &, iterations                                                      &
                          ! Number of microphysics iterations
                          ! (for rate diagnostic)
     &, ice_type          ! Type of ice (0 - crystals, 1 - aggregates
                          !              3 - graupel)
!
      Real, Intent(In) ::                                               &
     &  Timestep                                                        &
                          ! Timestep / s
     &, area_liq(points)                                                &
                          ! Fraction of gridbox with liquid-only cloud
     &, area_mix(points)                                                &
                          ! Fraction of gridbox with mixed phase cloud
     &, cfliq(points)                                                   &
                          ! Fraction of gridbox with liquid cloud
     &, cficei(points)                                                  &
                          ! 1/Fraction of gridbox with ice cloud
!    &, cf(points)        ! Current cloud fraction
!    &, cfl(points)       ! Current liquid cloud fraction
!    &, cff(points)       ! Current ice cloud fraction
     &, rho(points)                                                     &
                          ! Air density / kg m-3
     &, m0                                                              &
                          ! Seed ice water content / kg kg-1
     &, tcg(points)                                                     & 
                          ! T dependent function in ice size dist'n
     &, tcgi(points)                                                    &
                          ! 1/tcg (no units)
     &, corr(points)                                                    &
                          ! Fall velocity correction factor (no units)
     &, lfrcp                                                           &
                          ! Latent heat of fusion
                          ! / heat capacity of air / K
     &, ai, bi, aic, bic
                          ! Ice mass-size relationships m(D) = ai D^bi
!
      Real, Intent(InOut) ::                                            &
     &  qcl(points)                                                     &
                          ! Liquid water content / kg kg-1
     &, qcf(points)                                                     &
                          ! Ice water content    / kg kg-1
     &, T(points)                                                       &
                          ! Temperature / K
     &, ptransfer(points) ! Mass rimed in this timestep / kg kg-1
!
!      Real, Intent(InOut) ::
!    &, cf_transfer_rate(points) ! Cloud fraction increment this tstep
!    &, cfl_transfer_rate(points)! Liquid cloud fraction inc this tstep
!    &, cff_transfer_rate(points)! Ice cloud fraction inc this tstep
!
      Logical, Intent(In) ::                                            &
     &  l_update_cf                                                     &
                          ! Update cloud fractions
     &, l_seq                                                           &
                          ! Carry out sequential updating
     &, l_psd
                          ! Use generic ice particle size distribution
!
! Local Variables
!
      Integer                                                           &
     &  i                                                               &
                          ! Loop counter for points
     &, cry_offset        ! Index offset for ice crystals
!
      Real                                                              &
     &  qclnew(points)                                                  &
                          ! For value of qcl after riming  / kg kg-1
     &, dqi(points)                                                     &
                          ! Amount of ice rimed  / kg kg-1
     &, m_2_di(points)
                          ! 2+DI moment of particle size distribution

      !-----------------------------------------------
      ! Select appropriate ice parametrization (see c_lspsiz)
      !-----------------------------------------------
      cry_offset=ice_type*ice_type_offset

      If (l_psd) then
        ! Use the generic ice particle size distribution 
        ! Calculate the 2+di (cx(81)) moment of the
        ! ice particle size distribution. 

! DEPENDS ON: lsp_moments
        Call lsp_moments(points,rho,T,qcf,cficei,ai,bi,cx(81),m_2_di)

      End if

      Do i = 1, points

        If (qcf(i) >  m0 .and. qcl(i) >  0.0 .and. T(i)  <   zerodegc   &
     &      .and. area_mix(i) >  0.0 .and. cfliq(i) >  0.0) then

          !-----------------------------------------------
          ! Calculate water content of mixed phase region
          !-----------------------------------------------
          If (l_psd) then

            ! Calculate the riming rate using the generic PSD
            ! constp(81) = (pi/4) ci
            qclnew(i) = qcl(i) /                                        &
     &                (cfliq(i)+cfliq(i)*constp(81)                     &
     &                *corr(i)*timestep*m_2_di(i))

          Else
            ! Use the defined gamma distribution

            qclnew(i) = qcl(i) /                                        &
     &                (cfliq(i)+cfliq(i)*constp(9+cry_offset)*tcg(i)    &
     &                *corr(i)*timestep*(rho(i)*qcf(i)*cficei(i)        &
     &                *constp(5+cry_offset)*tcgi(i))**cx(6+cry_offset))
          End if

          !-----------------------------------------------
          ! Convert to new grid box total water content
          !-----------------------------------------------
          qclnew(i)=qcl(i)*area_liq(i)/cfliq(i)+qclnew(i)*area_mix(i)
          dqi(i)=(qcl(i)-qclnew(i))


          !-----------------------------------------------
          ! Store process rate / kg kg-1 s-1 and cloud fraction changes
          !-----------------------------------------------
          ptransfer(i) = ptransfer(i) + dqi(i)/(timestep*iterations)

!         There are no cloud fraction updates associated
!         with the riming term. They will have already been set to
!         zero on input to this subroutine, which is why these
!         are commented out.
!          cf_transfer_rate(i)  = 0.0 / (timestep*iterations)
!          cfl_transfer_rate(i) = 0.0 / (timestep*iterations)
!          cff_transfer_rate(i) = 0.0 / (timestep*iterations)

          !-----------------------------------------------
          ! Update water contents
          !-----------------------------------------------
          If (l_seq) then
            qcf(i)=qcf(i)+dqi(i)
            T(i)=T(i)+lfrcp*dqi(i)
            qcl(i)=qclnew(i)
          End if

          !-----------------------------------------------
          ! Update cloud fractions
          !-----------------------------------------------
!          If (l_update_cf) Then
!            These are commented out since there is currently no
!            cloud fraction update associated with the riming term.
!            If (l_seq) Then
!              cf(i)  = cf(i) + cf_transfer_rate(i) *timestep*iterations
!              cfl(i) = cfl(i)+ cfl_transfer_rate(i)*timestep*iterations
!              cff(i) = cff(i)+ cff_transfer_rate(i)*timestep*iterations
!            End If  ! l_seq
!            cftransfer(i)  = cftransfer(i)  + cf_transfer_rate(i)
!            cfltransfer(i) = cfltransfer(i) + cfl_transfer_rate(i)
!            cfftransfer(i) = cfftransfer(i) + cff_transfer_rate(i)
!
!          End If  !_update_cf

         End If ! qcf(i) >  m0 etc.

      End Do  ! Points

      Return  ! End of the subroutine
      END SUBROUTINE LSP_RIMING
