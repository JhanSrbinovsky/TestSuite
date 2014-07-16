#if defined(A09_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ PC2 Cloud Scheme: Calculation of RH(TL) for later use by initiation
! Subroutine Interface:
      SUBROUTINE PC2_RHTL(                                              &
!      Parallel variables
     &  halo_i, halo_j, off_x, off_y                                    &
!      Array dimensions
     &,wet_levels, row_length, rows                                     &
!      Prognostic arrays
     &,theta, exner_theta_levels, q, qcl, p_theta_levels                &
!      Output value of RH(TL)
     &,rhts                                                             &
!      Logical control
     &,l_mixing_ratio                                                   &
     & )
!
      IMPLICIT NONE
!
! Purpose:
!   Calculate total relative humidity wrt T_L
!
! Method:
!   Straight calculation using RH_T(TL) = (q+qcl)/qsatwat(TL)
!
! Current Owner of Code: Damian Wilson
!
! History:
! Version   Date     Comment
!  6.2    02-03-05   New Deck (Damian Wilson)
!  6.4    18-08-06   Use mixing ratio formulation.  Damian Wilson
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: PC2 Cloud Scheme Documentation
!
!  Global Variables:----------------------------------------------------
#include "c_lheat.h"
#include "c_r_cp.h"
!
!  Subroutine Arguments:------------------------------------------------
      INTEGER           , intent(in) ::                                 &
!      Model dimensions
     & wet_levels, row_length, rows                                     &
!
!      Parallel setup variables
     &, halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j                                                          &
                   ! Size of halo in j direction.
     &, off_x                                                           &
                   ! Size of small halo in i.
     &, off_y      ! Size of small halo in j.
!
      REAL              , intent(in) ::                                 &
!      Model prognostics
     &  theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &          wet_levels)                                             &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, wet_levels)              &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &        wet_levels)                                               &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_levels)                                               &
     &, p_theta_levels(1-off_x:row_length+off_x,                        &
     &                   1-off_y:rows+off_y, wet_levels)
!
      Logical           , intent(in) ::                                 &
     & l_mixing_ratio     ! Use mixing ratio formulation
!
      REAL              , intent(out) ::                                &
!      RH_T(T_L)
     & rhts(row_length,rows,wet_levels)
!
!  External functions:
!
!  Local parameters and other physical constants------------------------
      REAL                                                              &
     & LCRCP
!       Latent heat of condensation divided by heat capacity of air.
!
      PARAMETER(                                                        &
     &          LCRCP=LC/CP                                             &
     &          )
!
!  Local arrays---------------------------------------------------------
!
      REAL                                                              &
     & TL(row_length,rows)                                              &
                                   ! Liquid water temperature (K)
     &,QSL_TL(row_length,rows)                                          &
                                   ! Qsat wrt liquid water at temp TL
     &,p_no_halos(row_length,rows) ! Pressure without halo values (Pa)
!
!  Local scalars
!
      INTEGER K,I,J     ! Loop counters: K - vertical level index
!                         I,J - horizontal position index
!
! ==Main Block==--------------------------------------------------------
!
      do k=1,wet_levels
        do j=1,rows
          do i=1,row_length
!           Calculate liquid temperature TL
            TL(i,j)=theta(i,j,k)*exner_theta_levels(i,j,k)              &
     &              -LCRCP*QCL(i,j,k)
            p_no_halos(i,j)=p_theta_levels(i,j,k)
          end do
        end do
!
!       Calculate qsat(TL) with respect to liquid water
! DEPENDS ON: qsat_wat_mix
        call qsat_wat_mix(qsl_tl,tl,p_no_halos,row_length*rows,         &
     &                l_mixing_ratio)
!
        do j=1,rows
          do i=1,row_length
!           Calculate RH_T(TL)
            RHTS(i,j,k)=(Q(i,j,k)+QCL(i,j,k))/QSL_TL(i,j)
          end do
        end do
      end do  ! k
!
! End of the subroutine
!
      RETURN
      END SUBROUTINE PC2_RHTL
#endif
