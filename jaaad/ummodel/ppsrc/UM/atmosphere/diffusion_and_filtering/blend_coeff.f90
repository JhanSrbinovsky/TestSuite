
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Blend_coeff
!
          Subroutine Blend_coeff(visc_BL_m, visc_BL_h, RHOKM, RHOKH     &
     &,                         row_length, rows, BL_LEVELS )


! Purpose: Blends the mixing coefficients from the BL scheme and the
!          subfilter turbulence scheme which will then be passed out
!          of the BL scheme into the horizontal diffusion subroutines.
!
!          The details of the blending can be changed if desired.
!
! Method:     Calculates "blended" diffusion coefficient visc_BL_m
!             by adding the reciprocals of the diffusion coefficients
!             from the BL scheme and the subgrid turbulence scheme, i.e.
!             visc_BL_m = (1.0/RHOKM) + (1.0/visc_BL_m)
!
! Original Programmer:   Carol Halliwell
! Current code owner: Carol Halliwell

!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
      Implicit None

      Integer                                                           &
     & row_length                                                       &
     &,rows                                                             &
     &,BL_LEVELS                                                        &
     &,i,j,k

      Real                                                              &
     & RHOKM(1:row_length,1:rows,BL_LEVELS)                             &
!                             ! Layer K-1 - to - layer K
!                               turbulent mixing coefficient
!                               for momentum from BL scheme
     &,RHOKH(row_length,rows,BL_LEVELS)                                 &
!                             ! Layer K-1 - to - layer K
!                               turbulent mixing coefficient
!                               for heat and moisture from BL scheme
     &,visc_BL_m(1:row_length, 1:rows,bl_levels)                        &
!                 ! INOUT mixing coefficient for momentum
!                   IN:  from subgrid turbulence scheme
!                   OUT: values blended from BL and turbulence schemes
     &,visc_BL_h(1:row_length, 1:rows,bl_levels)                        &
!                 ! INOUT mixing coefficient for heat and moisture
!                   IN:  from subgrid turbulence scheme
!                   OUT: values blended from BL and turbulence schemes
     &,recip_RHOKM                                                      &
     &,recip_RHOKH                                                      &
     &,recip_visc_BL_m                                                  &
     &,recip_visc_BL_h      

      REAL,PARAMETER:: ZERO=0.0
      REAL,PARAMETER:: SMALLP=TINY(ZERO)

      Do k = 1, BL_LEVELS
        Do j = 1, rows
          Do i = 1, row_length

            If (RHOKM(i,j,k) >= SMALLP .AND.                            &
&               visc_BL_m(i,j,k) >= SMALLP) Then

              recip_RHOKM = 1.0/RHOKM(i,j,k)
              recip_visc_BL_m = 1.0/visc_BL_m(i,j,k)
              visc_BL_m(i,j,k) = 1.0/(recip_RHOKM + recip_visc_BL_m)

            Else

              visc_BL_m(i,j,k) = 0.0

            End If

          End Do
        End Do
      End Do

! End of routine
      return
      End Subroutine BLEND_COEFF


