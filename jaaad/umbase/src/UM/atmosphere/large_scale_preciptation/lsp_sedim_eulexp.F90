#if defined(A04_3D)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Microphysics hydrometeor Eulerian sedimentation scheme

      SUBROUTINE LSP_SEDIM_EULEXP(                                      &
     &  LSITER,Points,M0,DHI,DHIR,Rho,RhoR                              &
     &, Flux_FromAbove, Fallspeed_ThisLayer                             &
     &, MixRatio_ThisLayer, FallSpeed_FromAbove                         &
     &, Total_Flux_Out)

      IMPLICIT NONE
!
! Description:
!   First order Eulerian advection scheme for hydrometeor sedimentation
!   with an exponential-based limiter to ensure stability for CFL>1
!   used by routine LSP_ICE for ice crystals, snow, rain and graupel.
!
! Method:
!   Based on method described in Rotstayn (1997)(QJRMS, 123, 1227-1282)
!
! Current Code Owner: Jonathan Wilkinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!

! Subroutine arguments

      ! Intent (In)
      Integer                                                           &
     &  LSITER                                                          &
                     ! number of iterations of microphysics param
     &, Points       ! number of points to process

      Real                                                              &
     &  M0                                                              &
                           ! Small mass (kg/kg) defined in c_lspmic
     &, DHI(Points)                                                     &
                           ! CFL limit (s m-1)
     &, DHIR(Points)                                                    &
                           ! 1.0/DHI (m s-1)
     &, Rho(Points)                                                     &
                           ! Air density (kg m-3)
     &, RhoR(Points)                                                    &
                           ! 1.0/Rho
     &, Flux_FromAbove(Points)                                          &
     &, FallSpeed_ThisLayer(Points)

      ! Intent (InOut)
      Real                                                              &
     &  MixRatio_ThisLayer(Points)                                      &
     &, Fallspeed_FromAbove(Points)

      ! Intent (Out)
      Real                                                              &
     &  Total_Flux_Out(Points)

! Local variables

      Real                                                              &
     &  MixRatio_FromAbove(Points)                                      &
                                   ! Mixing Ratio from above
     &, Flux_out(Points)                                                &
                                   ! Temporary flux out of layer
     &, ExpFactor(Points)          ! Exponential Factor

      Integer i                    ! Loop counter

!-----------------------------------------------------------------------

      Do i=1,Points

        ! -------------------------------------------------------------
        ! Adjust fall speeds to make a linear combination of the fall
        ! speed from the layer above. This ensures that the fall speed
        ! in this layer will not be calculated as zero even though
        ! there is mass falling into it from above.
        ! -------------------------------------------------------------

        MixRatio_FromAbove(i) = Flux_FromAbove(i)*DHI(i)*RhoR(i)

        If (MixRatio_ThisLayer(i) + MixRatio_FromAbove(i)  >   M0) THEN

          FallSpeed_ThisLayer(i) =                                      &
     &          (FallSpeed_ThisLayer(i)*MixRatio_ThisLayer(i)           &
     &          +FallSpeed_FromAbove(i)*MixRatio_FromAbove(i))          &
     &          /(MixRatio_ThisLayer(i)+MixRatio_FromAbove(i))

        Else

          FallSpeed_ThisLayer(i) = 0.0

        End If

        ! -------------------------------------------------------------
        ! Eulerian solution with exponential limiter
        ! -------------------------------------------------------------

        If (FallSpeed_ThisLayer(i)  >   0.0) Then

          ExpFactor(i) = EXP(-1.0*FallSpeed_ThisLayer(i)*DHI(i))

          ! Calculate flux out of this layer

          Flux_out(i) = Flux_FromAbove(i)+DHIR(i)                       &
     &                 *(Rho(i)*MixRatio_ThisLayer(i)-Flux_FromAbove(i) &
     &                 /FallSpeed_ThisLayer(i))*(1.0-ExpFactor(i))

          ! Calculate mass (kg/kg) that remains in this layer

          MixRatio_ThisLayer(i) = Flux_FromAbove(i)*RhoR(i)             &
     &                     / FallSpeed_ThisLayer(i)                     &
     &                     * (1.0-ExpFactor(i)) + MixRatio_ThisLayer(i) &
     &                     * ExpFactor(i)

        Else

          ! No fall out of the layer.
          ! Set MixingRatio to be the amount of mass falling in.
          ! FallSpeed can only be zero if MixingRatio_ThisLayer LE M0
          ! and MixingRatio_FromAbove LE M0
          ! so this is slightly inconsistent.
          Flux_Out(i)       = 0.0
          MixRatio_ThisLayer(i) = Flux_FromAbove(i)*RhoR(i)*DHI(i)

        End If

        ! No need to compute fall speed out of the layer in this method
        ! -------------------------------------------------------------
        !  Total_Flux_Out is the flux out of this layer that falls
        !  into the next layer down
        ! -------------------------------------------------------------

        Total_Flux_Out(i) = Total_Flux_Out(i) + Flux_Out(i)/LSITER

        ! -------------------------------------------------------------
        ! Store fall speed in this layer to be the fallspeed from above
        ! for the next layer down
        ! -------------------------------------------------------------

        FallSpeed_FromAbove(i) = FallSpeed_ThisLayer(i)

      End Do ! on loop over points

      Return
      END SUBROUTINE LSP_SEDIM_EULEXP
#endif
