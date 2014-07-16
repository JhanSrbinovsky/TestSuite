
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! Subroutine Set_coeff_lagrange --------------------------------------
!
! Purpose : Set up denominator terms for Lagrange interpolation
!
!         variable
!

      SUBROUTINE Set_coeff_lagrange(                                    &
     &                        lamphi, rowcol, halo, L_dprint,           &
     &                        lamphi_rm,  lamphi_rp,                    &
     &                        recip_lamphi_m,  recip_lamphi,            &
     &                        recip_lamphi_p, recip_lamphi_p2 )

      IMPLICIT NONE

! Input arguments
      
       Logical                                                          &
     &  L_dprint  ! True if diagnostic prints needed
      
     INTEGER                                                            &
     &  rowcol                                                          &
                 !IN total number of point in a row or column
     &, halo       
                 !IN halo size of lamphi array

      REAL:: lamphi ( 1-halo : rowcol + halo )

! Output arguments
      REAL                                                              &
     &  lamphi_rm  (1 - halo : rowcol + halo)                           &
     &, lamphi_rp  (1 - halo : rowcol + halo)                           &
     &, recip_lamphi_m  (1 - halo : rowcol + halo)                      &
     &, recip_lamphi    (1 - halo : rowcol + halo)                      &
     &, recip_lamphi_p  (1 - halo : rowcol + halo)                      &
     &, recip_lamphi_p2 (1 - halo : rowcol + halo)

! loop counters
       Integer                                                          &
     &  i                                                               &
     &, info

      REAL                                                              &
     &  rm, r, rp                                                       &
     &, x1, x2, x3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Section 1:  Set arrays and pointers for searching on variable grid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (L_dprint) then
        write(*,*)' Set reciprocal arrays for Lagrange interpolation'
        write(*,*)' Input field has ', rowcol,' points plus halo size'  &
     &                               , halo
      endif !  L_dprint

      do i = 1 - halo, rowcol + halo
        lamphi_rm(i) = 0.0
        lamphi_rp(i) = 0.0 
        recip_lamphi_m(i)  = 0.0
        recip_lamphi(i)    = 0.0
        recip_lamphi_p(i)  = 0.0
        recip_lamphi_p2(i) = 0.0
      end do ! i = 1 - halo, rowcol + halo

      i = 2 - halo
      r    = lamphi(i+1) - lamphi(i)
      rm   = lamphi(i) - lamphi(i-1)
      do i = 2 - halo, rowcol + halo - 2
        rp  = lamphi(i+2) - lamphi(i+1)
        x1 = rm / r
        x2 = rp / r
        x3 = x1 + x2
        lamphi_rm(i) = x1
        lamphi_rp(i) = x2
        recip_lamphi_m(i)    = 1.0 / (x1 * (x1 + 1.0) * (x3 + 1.0))
        recip_lamphi(i)      = 1.0 / (x1 * (x2 + 1.0))
        recip_lamphi_p(i)    = 1.0 / ((x1 + 1.0) * x2)
        recip_lamphi_p2(i)   = 1.0 / ((x3 + 1.0) * (x2 + 1.0) * x2)
        rm = r
        r  = rp
      end do ! i = 2 - halo, rowcol + halo - 2

      Return
      End SUBROUTINE Set_coeff_lagrange

