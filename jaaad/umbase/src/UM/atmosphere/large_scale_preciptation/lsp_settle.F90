#if defined(A04_3D)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Accretion of droplets by raindrops
! Subroutine Interface:
      SUBROUTINE LSP_SETTLE(                                            &
     &  points, timestep, iterations                                    &
                                          ! Number of points and tstep
     &, q, qcl, T, droplet_flux                                         &
                                          ! Water contents and temp.
     &, l_update_cf, l_seq, bland                                       &
                                          ! Control logicals
     &, cfliq                                                           &
                                          ! Liquid cloud fraction
     &, rho, rhor, corr2, lcrcp                                         &
                                          ! Parametrization information
     &, dhi, dhir                                                       &
                                          ! Layer thickness information
     &, Ntot_land, Ntot_sea                                             &
                                          ! Droplet sizes
     &, ptransfer_qcl, ptransfer_q                                      &
                                          ! Transfer rates
     &  )
!
      Implicit None
!
! Purpose:
!   Update cloud prognostics as a result of gravitational droplet
!   settling.
!
! Method:
!   Solve for the bulk settling velocity based on an order 2 gamma
!   droplet size distribution.
!   Stokes Law is used to derive the fall speeds. This is accurate to
!   around 30 microns radius. The relationship from Rogers and Yau is
!   used.  vt = 2/9 r^2 g rho_wat / mu and mu is the dynamic viscosity
!   of air (= 1.717 x 10-5 at 0C).
!             = k1 r^2 where k1 = 1.27E8 m-1 s-1 / corr2
!   Integrating over the distribution we see that:
!   <vt> = 42 k1 [ 2 * lwc rho / (160 N_drop rho_wat pi) ]^(2/3)
!    = 1.339E6 kg^2/3 (lwc rho / N_drop)^(2/3) / corr2
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
                          ! Number of points to calculate
     &, iterations        ! Number of microphysics iterations
                          ! (for rate diagnostic)
!
      Real, Intent(In) ::                                               &
     &  Timestep                                                        &
                          ! Timestep / s
     &, cfliq(points)                                                   &
                          ! Liquid cloud fraction at start of timestep
     &, rho(points)                                                     &
                          ! Air density / kg m-3
     &, rhor(points)                                                    &
                          ! 1 / air density / m3 kg-1
     &, corr2(points)                                                   &
                          ! Air diffusivitiy correction (no units)
     &, dhi(points)                                                     &
                          ! Timestep / layer thickness / s m-1
     &, dhir(points)                                                    &
                          ! 1 / dhi  / m s-1
     &, lcrcp                                                           &
                          ! Latent heat condensation / cp / K
     &, Ntot_land                                                       &
                          ! Droplet concentration over land / m-3
     &, Ntot_sea          ! Droplet concentration over sea / m-3
!
      Real, Intent(InOut) ::                                            &
     &  q(points)                                                       &
                          ! Vapour content / kg kg-1
     &, qcl(points)                                                     &
                          ! Liquid water content / kg kg-1
     &, T(points)                                                       &
                          ! Temperature / K
     &, droplet_flux(points)                                            &
                              ! On input: Flux of water into layer
                              ! On output: Flux of water out of layer
                              ! / kg m-2 s-1
     &, ptransfer_qcl(points)                                           &
                              ! Transfer rate of qcl / kg kg-1 s-1
     &, ptransfer_q(points)   ! Transfer rate of q / kg kg-1 s-1
!
      Logical, Intent(In) ::                                            &
     &  l_update_cf                                                     &
                          ! Update cloud fractions
     &, l_seq                                                           &
                          ! Carry out sequential updating
     &, bland(points)     ! Land sea mask
!
! Local Variables
!
      Integer                                                           &
     &  i                 ! Loop counter

      Real                                                              &
     &  ndrop(points)                                                   &
                          ! Droplet number / m-3
     &, vt_droplet(points)                                              &
                          ! Average droplet settling velocity / m s-1
     &, fqirqi(points)                                                  &
                          ! Flux of liquid water
                          ! out of model layer / kg m-2 s-1
     &, flux_into_cloud(points)                                         &
                                ! Flux of liquid water from layer above
                          ! falling into cloud / kg m-2 s-1
     &, flux_into_clear_sky(points)                                     &
                                    ! Flux of liquid from layer above
                          ! falling into clear sky / kg m-2 s-1
     &, dqcl(points)                                                    &
                          ! Change in qcl this timestep / kg kg-1
     &, dq(points)        ! Change in q this timestep / kg kg-1

      Real, Parameter:: two_thirds = 2.0/3.0

      Do i=1,points

        !-----------------------------------------------
        ! Calculate droplet number
        !-----------------------------------------------
        ! At the moment we are going to used the fixed values
        ! given by Ntot_land and Ntot_sea rather than use
        ! the ndrop function.
        If (bland(i)) then
          ndrop(i) = Ntot_land
        Else
          ndrop(i) = Ntot_sea
        End if

        !-----------------------------------------------
        ! Calculate settling velocity
        !-----------------------------------------------
        vt_droplet(i) = 1.339E6 * (qcl(i) * rho(i) / ndrop(i))          &
     &                           ** two_thirds / corr2(i)

        !-----------------------------------------------
        ! Calculate flux of water downwards
        !-----------------------------------------------
        ! Limit to the size of model layer and timestep
        fqirqi(i) = rho(i)*qcl(i) * min( vt_droplet(i) , dhir(i) )

        !-----------------------------------------------
        ! Calculate transfer
        !-----------------------------------------------
        ! If we are using PC2 then we will need to evaporate drops that
        ! fall into clear sky. We assume at the moment a uniform
        ! distribution falling into the gridbox rather than pass around
        ! information about the cloud fraction in the layer above.

        flux_into_cloud(i) = cfliq(i) * droplet_flux(i)
        flux_into_clear_sky(i) = (1.0 - cfliq(i)) * droplet_flux(i)

        dqcl(i) = (flux_into_cloud(i)-fqirqi(i))*dhi(i)*rhor(i)
        dq(i)   =  flux_into_clear_sky(i) * dhi(i) * rhor(i)

        If (.not. l_update_cf) then
          ! We will allow the large scale cloud scheme to decide
          ! upon the cloud changes hence sum the two increments
          dqcl(i) = dqcl(i) + dq(i)
          dq(i)  = 0.0
        End if  ! .not. l_update_cf

        ! Update droplet flux leaving the layer
        droplet_flux(i) = fqirqi(i)

        ! Store transfer rate diagnostic
        ptransfer_qcl(i) = ptransfer_qcl(i)                             &
     &                     + dqcl(i) / (timestep*iterations)
        ptransfer_q(i)   = ptransfer_q(i)                               &
     &                     + dq(i) / (timestep*iterations)

        !------------------------------------------------
        ! Adjust liquid content
        !------------------------------------------------
        If (L_seq) Then
          qcl(i) = qcl(i) + dqcl(i)
          q(i) = q(i) + dq(i)
          T(i) = T(i) - lcrcp * dq(i)
          ! There is no change in the cloud fractions as we
          ! assume drops falling into clear sky are evaporated.
        End If  ! L_seq

      End Do  ! Points

      Return  ! End of the subroutine
      END SUBROUTINE LSP_SETTLE
#endif
