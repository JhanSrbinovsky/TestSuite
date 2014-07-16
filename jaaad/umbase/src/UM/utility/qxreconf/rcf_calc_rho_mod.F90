#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculates RHO for the output dump

Module Rcf_Calc_Rho_Mod

!  Subroutine Rcf_Calc_Rho - calculates tho
!
! Description:
!   Calculates RHO for the 5.0/5.1 dump (rho should not be interpolated
!   so as to maintain dynamical balance)
!
! Method:
!   Code derived from New Dynamics 2.7
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.2   15/11/00   Allow use of P rather than Exner. P.Selwood.
!   5.3   27/07/01   Interpolate exner to theta levels directly, without
!                    using the intermediary pressure variable. S.Cusack
!   5.3   25/10/01   Cater for extra rho level at top. D.Robinson.
!   5.3   06/12/01   Correct inconsistency in rho calculation A.Malcolm
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Calc_Rho( theta, q, exner, p, theta_heights,        &
                         rho_heights, rho )

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Field_Equals_Mod, Only : &
    Rcf_Field_Equals

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

Use Rcf_Parvars_Mod, Only : &
    mype

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Use Rcf_Calc_Exner_Theta_Mod, Only : &
    Rcf_Calc_Exner_Theta

Implicit None

! Arguments
Type( field_type ), Intent( In )     :: theta
Type( field_type ), Intent( In )     :: q
Type( field_type ), Intent( In )     :: p           ! on rho levels
Type( field_type ), Intent( In )     :: exner       ! on rho levels
Type( field_type ), Intent( InOut )  :: rho

Real, Intent(In)                   :: rho_heights( theta % level_size, &
                                                 0 : theta % levels + 1)
Real, Intent(In)                   :: theta_heights(theta % level_size,&
                                                 0 : theta % levels + 1)

! Comdecks
#include "c_r_cp.h"
#include "c_epslon.h"
#include "cppxref.h"

! Local variables
Type( field_type )                 :: exner_theta ! on theta levels
Integer                            :: i
Integer                            :: k
Real                               :: weight1
Real                               :: weight2
Real                               :: weight3
Real                               :: temp
Real                               :: work_real ( theta % level_size )
Real                               :: work_real2( theta % level_size )

! Nullify data
Nullify (exner_theta % data)
Nullify (exner_theta % data_int)
Nullify (exner_theta % data_log)

!-------------------------------------------------------------------
! Need to calculate exner on theta levels
!-------------------------------------------------------------------
Call Rcf_field_equals( exner_theta, exner )
exner_theta % levels = theta % levels
Call Rcf_Alloc_Field( exner_theta )

! find exner on theta levels (held in exner_theta)
Call Rcf_calc_exner_theta( exner % level_size, exner_theta % levels, &
                           theta_heights, rho_heights(1:,1:),        &
                           exner % Data, exner_theta % Data )

!--------------------------------------------------------------------
! Now do the calculation of rho
!--------------------------------------------------------------------

If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
  Write (6,*) 'Calculating Rho'
End If

k = 1
Do i = 1, theta % level_size
! calculate thetav
  work_real(i) = theta % Data(i,k) * (1. +                             &
                                     (1./epsilon -1.) * q % Data(i,k) )
! calculate rho
  rho % Data(i,k) = rho_heights(i,k) * rho_heights(i,k) *           &
                 p % Data(i,k) / (R * work_real(i) * exner % Data(i,k))
End Do

Do k = 2, theta % levels
  Do i = 1, theta % level_size
    work_real2(i) = work_real(i)
  End Do

  If (k .le. q % levels ) Then
    Do i = 1, theta % level_size
    work_real(i) = theta % Data(i,k) * (1. +                           &
                                     (1./epsilon -1.) * q % Data(i,k) )
    End Do

  Else
    Do i = 1, theta % level_size
      work_real(i) = theta % Data(i,k)
    End Do
  End If

  If (k .ne. theta % levels) Then
    Do i = 1, theta % level_size
      weight1 = rho_heights(i,k) - theta_heights(i,k-1)
      weight2 = theta_heights(i,k) - rho_heights(i,k)
      weight3 = theta_heights(i,k) - theta_heights(i,k-1)

!      temp = ( weight2 * work_real(i) +             &
!               weight1 * work_real2(i) ) /          &
!               weight3

      temp = ( weight1 * work_real(i) +             &
               weight2 * work_real2(i) ) /          &
               weight3

      rho % Data(i,k) = rho_heights(i,k) * rho_heights(i,k)     &
                      * p % Data(i,k) / (R * temp * exner % Data(i,k))
    End Do

  Else

    Do i= 1, theta % level_size

      temp = work_real2(i) *                                   &
             exner_theta % Data(i, exner_theta % levels - 1) / &
             exner % Data(i,k)

      rho % Data(i,k) = rho_heights(i,k) * rho_heights(i,k) *      &
                        p % Data(i,k) / (R * temp * exner % Data(i,k) )

    End Do

  End If
End Do

!------------------------------------------------------------------
! Tidy up
!------------------------------------------------------------------
Call Rcf_DeAlloc_Field( exner_theta )

Return
End Subroutine Rcf_Calc_Rho
End Module Rcf_Calc_Rho_Mod
#endif
