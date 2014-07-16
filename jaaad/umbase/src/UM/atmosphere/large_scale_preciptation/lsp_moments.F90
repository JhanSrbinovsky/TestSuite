#if defined(A04_3D) || defined(A04_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Conversion between moments of PSD
! Subroutine Interface:
      SUBROUTINE lsp_moments(                                           &
     &  points                                                          &
     &, rho, T, qcf, cficei                                             &
     &, ai, bi, n_out, moment_out                                       &
     &  )
      Implicit None
!
! Purpose:
!   Obtains a specified moment of the in-cloud  ice particle size
!   distribution given grid box mean input ice water content.
!
! Method:
!   Follows Field et al, 2005. Parametrization of ice-particle size
!   distributions for mid-latitude stratiform cloud. Quart. J. Royal
!   Meterol. Soc., 131, pp 1997-2017.
!
! Current Owner of Code: Jonathan Wilkinson
!
!
! Description of Code:
!   Fortran 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP 26.
!       
!   This subroutine calculates the quantity p_moment_out = M(n_out)
!   given the quantity (rho qcf / cfice)=ai*M(bi) where M(bi) is the
!   bi'th moment of the particle size distribution (PSD), corresponding
!   to the mass diameter relationship m(D) = ai D**bi.
!
!   It first calculates M(2) given the input quantities rho, qcf and ai
!   and then uses M(2) to calculate the output quantity p_out*M(out),
!   which can be used directly in the transfer calculations.
!   The conversion is from Field et al and gives
!   M(n)=a(n,Tc)*M(2)**b(n,Tc) where a and b are specified parameters.
!
!   This subroutine can be used to convert directly from M(x) to M(y)
!   by inputting ai=1, bi=x, n_out=y, rho=[array of 1's], qcf=M(x),
!   cficei=[array of 1's] and outputting moment_out=M(y).
!
!
! Subroutine Arguments
!
#include "c_0_dg_c.h"
      Integer, Intent(In) ::                                            &
     &  points
                        ! Number of points to calculate

      Real, Intent(In) ::                                               &
     &  ai                                                              &
                        ! Prefactor in determining iwc from bi'th
                        ! moment of the particle size distribution
     &, bi                                                              &
                        ! Moment of the particle size distribution
                        ! corresponding to the ice water content
     &, n_out                                                           &
                        ! Moment of the PSD to be output
     &, rho(points)                                                     &
                        ! Air density / kg m-3
     &, T(points)                                                       &
                        ! Temperature / K
     &, qcf(points)                                                     &
                        ! Ice water content / kg kg-1
     &, cficei(points)
                        ! 1 / ice cloud fraction

      Real, Intent(Out) ::                                              &
     &  moment_out(points)
                        ! n_out moment of the in-cloud
                        ! particle size distn
!
! Local Variables
!
      Integer                                                           &
     &  i
                        ! Loop counter for points

      Real                                                              &
     &  one_over_ai                                                     &
                        ! 1 / ai
     &, Tc(points)                                                      &
                        ! Temperature in degrees Celsius
     &, log10_abi(points)                                               &
                        ! log10 of conversion factor a(bi,Tc)
     &, abi(points)                                                     &
                        ! Conversion factor a(bi,Tc)
     &, bbi(points)                                                     &
                        ! Conversion factor b(bi,Tc)
     &, m_2(points)                                                     &
                        ! Second moment of the PSD
     &, m_bi(points)                                                    &
                        ! bi'th moment of the PSD
     &, log10_an_out(points)                                            &
                        ! log10 of conversion factor a(n_out,Tc)
     &, an_out(points)                                                  &
                        ! Conversion factor a(n_out,Tc)
     &, bn_out(points)
                        ! Conversion factor b(n_out,Tc)

      ! The following values are from Table 1 of Field et al
      ! and represent the conversion parameters between moments.
      Real, Parameter:: a(10) = (/5.065339,-0.062659,-3.032362,         &
     &  0.029469,-0.000285,0.312550,0.000204,0.003199,0.000000,         &
     &  -0.015952 /)        
      Real, Parameter:: b(10) = (/0.476221,-0.015896,0.165977,          &
     &  0.007468,-0.000141,0.060366,0.000079,0.000594,0.000000,         &
     &  -0.003577 /)
!
! Start the subroutine
!
      one_over_ai = 1.0 / ai

      Do i = 1, points

        If (qcf(i) > 0.0) then
          !-----------------------------------------------
          ! Form the bi'th moment of the ice particle size
          ! distribution from the ice water content.
          !-----------------------------------------------
          m_bi(i) = rho(i) * qcf(i) * cficei(i) * one_over_ai
          m_bi(i) = max(m_bi(i) , 0.0)

          !-----------------------------------------------
          ! Calculate the second moment of the PSD
          !-----------------------------------------------
          Tc(i) = T(i) - zerodegc

          log10_abi(i) = a(1) + a(2)*Tc(i) + a(3)*bi                    &
     &              + a(4)*Tc(i)*bi + a(5)*Tc(i)*Tc(i)                  &
     &              + a(6)*bi*bi + a(7)*Tc(i)*Tc(i)*bi                  &
     &              + a(8)*Tc(i)*bi*bi + a(9)*Tc(i)*Tc(i)*Tc(i)         &
     &              + a(10)*bi*bi*bi
          abi(i) = 10.0**(log10_abi(i))

          bbi(i) =    b(1) + b(2)*Tc(i) + b(3)*bi                       &
     &              + b(4)*Tc(i)*bi + b(5)*Tc(i)*Tc(i)                  &
     &              + b(6)*bi*bi + b(7)*Tc(i)*Tc(i)*bi                  &
     &              + b(8)*Tc(i)*bi*bi + b(9)*Tc(i)*Tc(i)*Tc(i)         &
     &              + b(10)*bi*bi*bi

          m_2(i) = (m_bi(i) / abi(i))**(1.0/bbi(i)) 

          !-----------------------------------------------
          ! Calculate the n_out moment of the PSD
          !-----------------------------------------------
          log10_an_out(i) = a(1) + a(2)*Tc(i) + a(3)*n_out              &
     &              + a(4)*Tc(i)*n_out + a(5)*Tc(i)*Tc(i)               &
     &              + a(6)*n_out*n_out + a(7)*Tc(i)*Tc(i)*n_out         &
     &              + a(8)*Tc(i)*n_out*n_out + a(9)*Tc(i)*Tc(i)*Tc(i)   &
     &              + a(10)*n_out*n_out*n_out
          an_out(i) = 10.0**(log10_an_out(i))

          bn_out(i) = b(1) + b(2)*Tc(i) + b(3)*n_out                    &
     &              + b(4)*Tc(i)*n_out + b(5)*Tc(i)*Tc(i)               &
     &              + b(6)*n_out*n_out + b(7)*Tc(i)*Tc(i)*n_out         &
     &              + b(8)*Tc(i)*n_out*n_out + b(9)*Tc(i)*Tc(i)*Tc(i)   &
     &              + b(10)*n_out*n_out*n_out

          moment_out(i) = an_out(i) * m_2(i) ** bn_out(i)

        Else  ! qcf > 0
          !-----------------------------------------------
          ! Set the output moment to zero
          !-----------------------------------------------
          moment_out(i) = 0.0

        End if  ! qcf > 0

      End do  ! i

      Return  ! End of the subroutine
      END SUBROUTINE lsp_moments
#endif
