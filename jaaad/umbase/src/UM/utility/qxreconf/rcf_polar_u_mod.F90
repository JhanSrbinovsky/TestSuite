#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Averaging of Polar U rows

Module Rcf_Polar_U_Mod

!  Subroutine Rcf_Polar_U - Averaging of Polar U rows
!
! Description:
!   Performs a vector mean on Polar U rows.
!
! Method:
!   Magnitude and direction of V rows are calculated and used to
!   set adjecent U rows.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.4   11/03/02   Remove comment lines from after #include
!                                                  S. Carroll
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Polar_U( u, v, delta_lambda )

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Parvars_Mod, Only : &
    mype,               &
    atNorth,            &
    atSouth,            &
    datastart,          &
    gc_proc_row_group

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Implicit None

! Arguments
Type( field_type), Intent(InOut)  :: u
Type( field_type), Intent(In)     :: v
Real,              Intent(In)     :: delta_lambda

! Comdecks
! P_Over_180
#include "c_pi.h"

! Local Variables.

Integer      :: i
Integer      :: k
Integer      :: info

Real         :: a_np
Real         :: b_np
Real         :: a_sp
Real         :: b_sp
Real         :: longitude

Real         :: wrk( 2 * v % levels)
Real         :: rwrk(v % row_len, 2 * v % levels)
Real         :: mag_vector
Real         :: dir_vector
Real         :: delta_lambda_rad

External gcg_rvecsumr


If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Finding vector mean for polar u rows'
End If

delta_lambda_rad = pi_over_180 * delta_lambda
! ----------------------------------------------------------------------
! Section 1.   Calculate magnitude and direction of polar vector wind
!              from v component of wind on row around pole.
! ----------------------------------------------------------------------

If (atNorth) then

   Do k = 1, 2 * v % levels
      wrk(k) = 0.
   End Do


   Do k = 1, v % levels
      Do i = 1, v % row_len
         ! This strange way of writing things is to fox the Cray
         ! optimiser and thus retain bit comparison (compiler bug)
         longitude = datastart(1) + i - 2
         longitude = longitude * delta_lambda_rad
         rwrk(i,2*(k-1)+1) = v % Data( (v % rows -1) * v % row_len+i,k)&
                                 * cos( longitude )
         rwrk(i,2*k) = v % Data( (v % rows -1) * v % row_len +i,k) * &
                                   sin( longitude )
      End Do
   End Do

   call gcg_rvecsumr(v % row_len, v % row_len, 1, 2 * v % levels,   &
                     rwrk, gc_proc_row_group, info, wrk)

   Do k = 1, v % levels

      a_np = 2. * wrk(2*(k-1)+1) / v % glob_row_len
      b_np = 2. * wrk(2*k) / v % glob_row_len

      mag_vector = sqrt (a_np*a_np + b_np*b_np)
      If (a_np .eq. 0. .and. b_np .eq. 0.) Then
         dir_vector = 0.
      Else
         dir_vector = atan2 (b_np, a_np)
      End If

      Do i = 1, u % row_len
        u % Data( (u % rows - 1) * u % row_len + i,k) =           &
                 mag_vector * sin( ( datastart(1) + i - 1 - .5) * &
                 delta_lambda_rad - dir_vector )
      End Do
   End Do

End If


If (atSouth) then

   Do k = 1, 2 * v % levels
      wrk(k) = 0.
   End Do

   Do k = 1, v % levels
      Do i = 1, v % row_len
         ! This strange way of writing things is to fox the Cray
         ! optimiser and thus retain bit comparison (compiler bug)
         longitude = datastart(1) + i - 2
         longitude = longitude * delta_lambda_rad
         rwrk(i,2*(k-1)+1) = v % Data(i,k) * cos( longitude )
         rwrk(i,2*k) = v % Data(i,k) * sin( longitude )
      End Do
   End Do

   call gcg_rvecsumr(v % row_len, v % row_len, 1, 2 * v % levels, &
             rwrk, gc_proc_row_group, info, wrk)

   Do k = 1, v % levels

      a_sp = 2. * wrk(2*(k-1)+1) / v % glob_row_len
      b_sp = 2. * wrk(2*k) / v % glob_row_len

      mag_vector = sqrt (a_sp*a_sp + b_sp*b_sp)
      If (a_sp .eq. 0. .and. b_sp .eq. 0.) Then
         dir_vector = 0.
      Else
         dir_vector = atan2 (b_sp, a_sp)
      End If

      Do i = 1, u % row_len
        u % Data(i,k) = -mag_vector * sin( (datastart(1) + i -1 -.5) * &
                        delta_lambda_rad - dir_vector)
      End Do

   End Do

End If


Return
End Subroutine Rcf_Polar_U

End Module Rcf_Polar_U_Mod
#endif
