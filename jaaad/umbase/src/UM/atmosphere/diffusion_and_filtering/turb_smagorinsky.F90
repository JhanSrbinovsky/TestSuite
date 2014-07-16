#if defined(A13_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine TURB_Smagorinsky
!
          Subroutine TURB_Smagorinsky(                                  &
     &                     rows, row_length, n_rows                     &
     &,                    model_levels                                 &
     &,                    r_theta_levels, r_rho_levels                 &
     &,                    u, v, w, visc, shear                         &
     &,                    z0, RNEUTML, timestep                        &
     &,                    diff_factor, mix_factor, max_diff            &
     &,                    cos_theta_latitude                           &
     &,                    cos_v_latitude                               &
     &,                    delta_lambda, delta_phi)

! Description: Calculates coefficients for use in subgrid turbulence
!              scheme
!
! Method:
!
! Current code owner: Carol Halliwell
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 6.2      05/01/06  New code
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      Implicit None

#include "c_a.h"
! for earth radius
#include "parvars.h"
! halo information
#include "c_vkman.h"
!Von Karman constant

! Variables with Intent (In)

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of v rows.
     &, model_levels     ! number of model levels.

      Real                                                              &
     &  u(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)       &
     &, v(1-offx:row_length+offx, 1-offy:n_rows+offy, model_levels)     &
     &, w(1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels)     &
     &, z0(row_length,rows)                                             &
                             ! roughness length
     &, timestep                                                        &
     &, diff_factor                                                     &
                    ! factor between 0 and 1 multiplied by max diff
!                   ! coeff for numerical stability to obtain max
!                   ! diffusion coeff for run
     &, max_diff                                                        &
                    ! max diffusion coeff for run
     &, mix_factor                                                      &
                    ! lambda0 = mix_factor * gridlength
     &, r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                    1-halo_j:rows+halo_j, 0:model_levels)         &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                    1-halo_j:rows+halo_j, model_levels)           &
     &, delta_lambda                                                    &
     &, delta_phi                                                       &
     &, cos_theta_latitude (1-Offx:row_length+Offx,                     &
     &                                           1-Offy:rows+Offy)      &
     &, cos_v_latitude (1-Offx:row_length+Offx,1-Offy:n_rows+Offy)

!Local parameters

      Integer                                                           &
     & i,j,k

      Real                                                              &
     &  delta_x, delta_y                                                &
     &, delta_z(row_length, rows, model_levels)                         &
     &, delta_zn(row_length+offx, rows+offy, model_levels)              &
     &, delta_zn_u(row_length, rows, model_levels)                      &
     &, delta_zn_v(row_length, rows, model_levels)                      &
     &, rdz(row_length, rows, model_levels)                             &
     &, rdzn_u(1-offx:row_length+offx, rows, model_levels)              &
     &, rdzn_v(row_length, 1-offy:rows+offy, model_levels)              &
     &, smallp                                                          &
                   ! A small number
     &, RNEUTML(row_length, rows, model_levels)                         &
!               ! mixing length scale (m) (lambda)
     &, RNEUTML_SQ(row_length, rows, model_levels)                      &
!               ! square of RNEUTML
     &, RMLMAX                                                          &
!               ! basic mixing length (lambda0)
     &, z(row_length, rows, model_levels)                               &
     &, shear(row_length, rows, model_levels)

      Real                                                              &
     &  dx_rho(row_length, rows, model_levels)                          &
!               ! dx on rho points
     &, dx_theta_u(row_length, rows, model_levels)                      &
!               ! dx above u points on theta levels
     &, dx_rho_centre(row_length, rows, model_levels)                   &
!               ! dx in centre of grid box on rho levels
     &, dy_rho(row_length, rows, model_levels)                          &
!               ! dy on rho points
     &, dy_theta_v(row_length, rows, model_levels)                      &
!               ! dy above v points on theta levels
     &, dy_rho_centre(row_length, rows, model_levels)                   &
!               ! dy in centre of grid box on rho levels
     &, cx_rho(row_length, rows, model_levels)                          &
!               ! reciprocal of dx_rho
     &, cx_theta_u(1-offx:row_length+offx, rows, model_levels)          &
!               ! reciprocal of dx_theta_u
     &, cx_rho_centre(1-offx:row_length+offx, 1-offy:rows+offy          &
     &                                          , model_levels)         &
!               ! reciprocal of dx_rho_centre
     &, cy_rho(row_length, rows, model_levels)                          &
!               ! reciprocal of dy_rho
     &, cy_theta_v(row_length, 1-offy:rows+offy, model_levels)          &
!               ! reciprocal of dy_theta_v
     &, cy_rho_centre(1-offx:row_length+offx, 1-offy:rows+offy          &
     &                                          , model_levels)         &
!               ! reciprocal of dy_rho_centre
     &, cx2_rho(1-offx:row_length+offx, rows, model_levels)             &
!               ! square of cx_rho
     &, cy2_rho(row_length, 1-offy:rows+offy, model_levels)
!               ! square of cy_rho

      Parameter (smallp=1.e-14)

      Real                                                              &
     &  ssq11                                                           &
                      ! ii component of shear stress
     &, ssq22                                                           &
                      ! jj component of shear stress
     &, ssq33                                                           &
                      ! kk component of shear stress
     &, ssq13                                                           &
                      ! ik component of shear stress
     &, ssq23                                                           &
                      ! jk component of shear stress
     &, ssq12                                                           &
                      ! ij component of shear stress
     &, ssq(row_length, rows, model_levels)                             &
                                     ! Half squared strain rate
     &, sum(row_length, rows, model_levels)   ! sum of Sij (ssqij)

! Variables with Intent (Out)

      Real                                                              &
     &  visc(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j           &
     &                  , model_levels)   ! OUT: lambda^2*S

!----------------------------------------------------------------------
! Initialise variables
!----------------------------------------------------------------------
      Do k= 1, model_levels

        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
! viscosity is set to zero at the top of the domain.
            visc(i,j,k) = 0.0
          End Do
        End Do

        Do j = 1, rows
          Do i = 1, row_length
            sum(i,j,k) = 0.0
            ssq(i,j,k) = 0.0
          End Do
        End Do

      End Do

      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            z(i,j,k) = r_theta_levels(i,j,k)                            &
     &                        - r_theta_levels(i,j,0)
           End Do
         End Do
      End Do
!----------------------------------------------------------------------
! Calculate grid spacings
!----------------------------------------------------------------------
! Horizontal grid spacings used in calculation of rate of strain term
! delta_lambda and delta_phi are in radians
!
      delta_x = Earth_radius * delta_lambda
      delta_y = Earth_radius * delta_phi

      Do k = 1, model_levels

        Do j = 1, rows
          Do i = 1, row_length

            dx_rho(i,j,k) =                                             &
     &               r_rho_levels(i,j,k)                                &
     &             * delta_lambda * cos_theta_latitude(i,j)

            dx_theta_u(i,j,k) =                                         &
     &         0.5*(r_theta_levels(i,j,k) + r_theta_levels(i+1,j,k))    &
     &       * delta_lambda * cos_theta_latitude(i,j)

            dy_rho(i,j,k) =                                             &
     &               r_rho_levels(i,j,k) * delta_phi

            dy_theta_v(i,j,k) =                                         &
     &         0.5*(r_theta_levels(i,j,k) + r_theta_levels(i,j+1,k))    &
     &       * delta_phi

            dy_rho_centre(i,j,k) = 0.25*(                               &
     &         r_rho_levels(i,j+1,k) + r_rho_levels(i+1,j+1,k)          &
     &        +r_rho_levels(i,j,k) + r_rho_levels(i+1,j,k) )            &
     &        * delta_phi

            cx_rho(i,j,k) = 1./dx_rho(i,j,k)
            cx_theta_u(i,j,k) =  1./dx_theta_u(i,j,k)
            cy_rho(i,j,k) =  1./dy_rho(i,j,k)
            cy_theta_v(i,j,k) =  1./dy_theta_v(i,j,k)
            cy_rho_centre(i,j,k) = 1./dy_rho_centre(i,j,k)
            cx2_rho(i,j,k) = cx_rho(i,j,k)*cx_rho(i,j,k)
            cy2_rho(i,j,k) =  cy_rho(i,j,k)*cy_rho(i,j,k)
          End Do
        End Do

        Do j = 1, n_rows
          Do i = 1, row_length
            dx_rho_centre(i,j,k) = 0.25*(                               &
     &         r_rho_levels(i,j+1,k) + r_rho_levels(i+1,j+1,k)          &
     &        +r_rho_levels(i,j,k) + r_rho_levels(i+1,j,k) )            &
     &                * delta_lambda * cos_v_latitude(i,j)
            cx_rho_centre(i,j,k) =  1./dx_rho_centre(i,j,k)
          End Do
        End Do

      End Do   !k

      Do k = 1, model_levels
        Do i = 1, row_length
          dx_rho_centre(i,rows,k) = dx_rho_centre(i,n_rows,k)
          cx_rho_centre(i,rows,k) =  1./dx_rho_centre(i,n_rows,k)
        End Do
      End Do

! Vertical grid spacings used in calculation of rate of strain term

      Do k = 1, model_levels

        Do j = 1, rows
          Do i = 1, row_length
            delta_z(i,j,k) = r_theta_levels(i,j,k)                      &
     &                -r_theta_levels(i,j,k-1)
            rdz(i,j,k) = 1./delta_z(i,j,k)
          End Do
        End Do

        Do j = 1, rows+offy
          Do i = 1, row_length+offx
            If (k  ==  1) then
              delta_zn(i,j,k) = r_rho_levels(i,j,k)                     &
     &                         -r_theta_levels(i,j,0)
            Else
              delta_zn(i,j,k) = r_rho_levels(i,j,k)                     &
     &                        - r_rho_levels(i,j,k-1)
            End If
          End Do
        End Do

      End Do

      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            delta_zn_u(i,j,k) =                                         &
     &             0.5*(delta_zn(i,j,k)+delta_zn(i+1,j,k))
            delta_zn_v(i,j,k) =                                         &
     &             0.5*(delta_zn(i,j,k)+delta_zn(i,j+1,k))
            rdzn_u(i,j,k) = 1./delta_zn_u(i,j,k)
            rdzn_v(i,j,k) = 1./delta_zn_v(i,j,k)
          End Do
        End Do
      End Do

! DEPENDS ON: swap_bounds
      call Swap_Bounds(cx2_rho, row_length, rows, model_levels,         &
     &                 offx, 0, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
      call Swap_Bounds(cx_theta_u, row_length, rows, model_levels,      &
     &                 offx, 0, fld_type_u, .false.)

! DEPENDS ON: swap_bounds
      call Swap_Bounds(cx_rho_centre, row_length, rows, model_levels,   &
     &                 offx, offy, fld_type_v, .false.)

! DEPENDS ON: swap_bounds
      call Swap_Bounds(cy2_rho, row_length, rows, model_levels,         &
     &                 0, offy, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
      call Swap_Bounds(cy_theta_v, row_length, rows, model_levels,      &
     &                 0, offy, fld_type_v, .false.)

! DEPENDS ON: swap_bounds
      call Swap_Bounds(cy_rho_centre, row_length, rows, model_levels,   &
     &                 offx, offy, fld_type_v, .false.)

! DEPENDS ON: swap_bounds
      call Swap_Bounds(rdzn_u, row_length, rows, model_levels,          &
     &                 offx, 0, fld_type_u, .false.)

! DEPENDS ON: swap_bounds
      call Swap_Bounds(rdzn_v, row_length, rows, model_levels,          &
     &                 0, offy, fld_type_v, .false.)

!--------------------------------------------------------
! As in the LEM code
! _Now calculate half-squared strain rate SSQ on w-points
! _CX=1./DX, CY=1./DY, RDZ(K)=1./DZ(K),
!       RDZN(K) =1./DZN(K)
! _SSQ= 0.5*^^DU_I/DX_J+DU_J/DX_I^^**2
! _SSQIJ= (DU_I/DX_J+DU_J/DX_I)**2
! _Hence SSQ= SUM(I,J) {0.5*(SSQIJ)}
! _Note that for a simple shear S, SSQ=S**2
!   (as in the notation of Mason and Callen 1986)
!--------------------------------------------------------

      Do k = 1, model_levels - 1

        Do j = 1, rows

          Do i = 1, row_length

          ssq11 =                                                       &
     &     cx2_rho(i,j,k+1)*(u(i,j,k+1)-u(i-1,j,k+1))**2 +              &
     &     cx2_rho(i,j,k)*(u(i,j,k)-u(i-1,j,k))**2
          ssq22 =                                                       &
     &     cy2_rho(i,j,k+1)*(v(i,j,k+1)-v(i,j-1,k+1))**2+               &
     &     cy2_rho(i,j,k)*(v(i,j,k)-v(i,j-1,k))**2
          ssq33 =                                                       &
     &       ((w(i,j,k)-w(i,j,k-1))*rdz(i,j,k))**2 +                    &
     &       ((w(i,j,k+1)-w(i,j,k))*rdz(i,j,k+1))**2
          ssq13= (                                                      &
     &       ((u(i,j,k+1)-u(i,j,k))*rdzn_u(i,j,k+1)+                    &
     &       (w(i+1,j,k)-w(i,j,k))*cx_theta_u(i,j,k))**2+               &
     &       ((u(i-1,j,k+1)-u(i-1,j,k))*rdzn_u(i-1,j,k+1)               &
     &       +(w(i,j,k)-w(i-1,j,k))*cx_theta_u(i-1,j,k))**2             &
     &             )*0.5      ! _averaging ssq13 over 2 points
          ssq23 = (                                                     &
     &       ((w(i,j,k)-w(i,j-1,k))*cy_theta_v(i,j-1,k)                 &
     &      +(v(i,j-1,k+1)-v(i,j-1,k))*rdzn_v(i,j-1,k+1))**2+           &
     &       ((w(i,j+1,k)-w(i,j,k))*cy_theta_v(i,j,k)+                  &
     &       (v(i,j,k+1)-v(i,j,k))*rdzn_v(i,j,k+1))**2                  &
     &            )*0.5        ! _averaging ssq23 over 2 points
          ssq12 = (  (  (                                               &
     &       ((u(i-1,j,k)-u(i-1,j-1,k))*cy_rho_centre(i-1,j-1,k)        &
     &      +(v(i,j-1,k)-v(i-1,j-1,k))*cx_rho_centre(i-1,j-1,k))**2 +   &
     &       ((u(i-1,j+1,k)-u(i-1,j,k))*cy_rho_centre(i-1,j,k)          &
     &      +(v(i,j,k)-v(i-1,j,k))*cx_rho_centre(i-1,j,k))**2           &
     &              )  +  (                                             &
     &       ((u(i,j,k)-u(i,j-1,k))*cy_rho_centre(i,j-1,k)              &
     &      +(v(i+1,j-1,k)-v(i,j-1,k))*cx_rho_centre(i,j-1,k))**2 +     &
     &       ((u(i,j+1,k)-u(i,j,k))*cy_rho_centre(i,j,k)                &
     &      +(v(i+1,j,k)-v(i,j,k))*cx_rho_centre(i,j,k))**2             &
     &        ) )  +  ( (                                               &
     &       ((u(i-1,j,k+1)-u(i-1,j-1,k+1))*cy_rho_centre(i-1,j-1,k+1)  &
     &      +(v(i,j-1,k+1)-v(i-1,j-1,k+1))                              &
     &        *cx_rho_centre(i-1,j-1,k+1))**2+                          &
     &       ((u(i-1,j+1,k+1)-u(i-1,j,k+1))*cy_rho_centre(i-1,j,k+1)    &
     &      +(v(i,j,k+1)-v(i-1,j,k+1))*cx_rho_centre(i-1,j,k+1))**2     &
     &              )  +  (                                             &
     &       ((u(i,j,k+1)-u(i,j-1,k+1))*cy_rho_centre(i,j-1,k+1)        &
     &      +(v(i+1,j-1,k+1)-v(i,j-1,k+1))                              &
     &       *cx_rho_centre(i,j-1,k+1))**2+                             &
     &       ((u(i,j+1,k+1)-u(i,j,k+1))*cy_rho_centre(i,j,k+1)          &
     &      +(v(i+1,j,k+1)-v(i,j,k+1))*cx_rho_centre(i,j,k+1))**2       &
     &         ) ) )*0.125      ! _averaging ssq12 over 8 points

          ssq(i,j,k)=ssq11+ssq22+ssq33+ssq13+ssq23+ssq12+smallp
          sum(i,j,k) = sum(i,j,k) + ssq(i,j,k)

          End Do   ! on i

        End Do    ! on j

      End Do ! on k

      Do k = 1,model_levels-1
        Do j = 1, rows
          Do i = 1, row_length
           shear(i,j,k) = sqrt(sum(i,j,k)) ! already been *0.5 in ssq
          End Do
        End Do
      End Do

      RMLMAX = mix_factor * MAX(delta_x, delta_y)

      Do k = 1, model_levels-1
        Do j = 1, rows
          Do i = 1, row_length

            RNEUTML(i,j,K)=SQRT(1./                                     &
     &     (1./(VKMAN*(z(i,j,k)+z0(i,j)))**2+1./RMLMAX**2) )
            RNEUTML_SQ(i,j,K)=RNEUTML(i,j,K)*RNEUTML(i,j,K)
            visc(i,j,k) = RNEUTML_SQ(i,j,k)*shear(i,j,k)
!
! visc is now lambda^2*S and still needs to be multiplied by
! stability functions (done in ATMSTEP)
!
          End Do
        End Do
      End Do
!
! maximum diffusion coefficient allowed in this run
!
      max_diff = diff_factor*delta_x*delta_x/(8.0*timestep)

      return
      END SUBROUTINE TURB_Smagorinsky

#endif
