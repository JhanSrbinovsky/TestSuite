#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!=====================================================================
! Subroutine ABNEW
! Purpose:-           To calculate amplitude and mean of sinusoidal
!                     distribution for stats. Eqns. 10 and 11
!                     in SCM documentation.
! Programmer:-        J. LEAN - modified code from original SCM to
!                     meet UM standards
!     Modification History:
! Version  Date
!  4.5     07/98     SCM integrated as a standard UM configuration
!                    Introduce multicolumn SCM
!                    JC Thil.
!  5.3     01/05/01    Update to 5.3 with an extra dimesion in each
!                      of the arrays.            Z. Gardner
!
!
!=====================================================================
!
!
!=====================================================================
! FUNCTION DAYNEW
! PURPOSE:-           To calculate SIN of argument (in eqn. 12
!                     in SCM doc.) required in calculation of
!                     mean or SD of variable at day relative to winter
!                     solstice
! PROGRAMMER:-        J. LEAN - modified code from original SCM to
!                     meet UM standards
!     Modification History:
! Version  Date
!  4.5     07/98     SCM integrated as a standard UM configuration
!                    JC Thil.
!=====================================================================
!
!
!=====================================================================
! Subroutine XNEW
! Purpose:-           To calculate mean or SD of random variable
!                     at daynumber relative to winter solstice
!                     (eqn. 12 in SCM doc.)
! Programmer:-        J. LEAN - modified code from original SCM to
!                     meet UM standards
!     Modification History:
! Version  Date
!  4.5     07/98     SCM integrated as a standard UM configuration
!                    Introduce multicolumn SCM
!                    JC Thil.
!  5.3     01/05/01    Update to 5.3 with an extra dimesion in each
!                      of the arrays.            Z. Gardner
!
!=====================================================================
!
!
!=====================================================================
! SUBROUTINE PNEW
! PURPOSE:-           To calculate pressure and reciprocal pressure
!                     coordinates AK and BK values and P*
! PROGRAMMER:-        J. LEAN - modified code from original SCM to
!                     meet UM standards
!
!     Modification History:
! Version  Date
!  4.5     07/98      SCM integrated as a standard UM configuration
!                     Introduce multicolumn SCM
!                     JC Thil.
!
!=====================================================================
!
!
!=====================================================================
! SUBROUTINE ACINIT
! PURPOSE:-           To calculate mean and SD of a random variable
!                     eqns 6 and 7 in SCM doc.
! PROGRAMMER:-        J. LEAN - modified code from original SCM to
!                     meet UM standards
!
!     Modification History:
! Version  Date
!  4.5     07/98      SCM integrated as a standard UM configuration
!                     Introduce multicolumn SCM
!                     JC Thil.
!
!=====================================================================
!
       Subroutine ACINIT(xbar, xsd, a, cbar, csd, cor, n,               &
     &                                         row_length, rows)
!---------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------
      Implicit none
      Integer                                                           &
     & n, row_length, rows  ! IN no of model_levels, rows and columns
      Real                                                              &
     &  a(row_length, rows,n-1)                                         &
                                ! OUT term a of eqn. 2.22
     &  ,cbar(row_length, rows,n-1)                                     &
                                    ! OUT Mean of random variable C
     &  ,cor(row_length, rows)                                          &
                                ! IN Vertical correlation coefficient
     &  ,csd(row_length, rows,n-1)                                      &
                                    ! OUT SD of random variable C
     &  ,xbar(row_length, rows,n)                                       &
                                    ! IN Mean of forcing variable
     &  ,xsd(row_length, rows,n)    ! IN SD of forcing variable
!---------------------------------------------------------------------
!     Local variables
!---------------------------------------------------------------------
      Integer                                                           &
     &  i, j, k                 ! Loop counters
!
      Do  k = 1 ,n-1
        Do j = 1, rows
          Do i =  1, row_length
            a(i,j,k) = cor(i,j) * xsd(i,j,k+1) / xsd(i,j,k)
            cbar(i,j,k) = xbar(i,j,k+1) - a(i,j,k) * xbar(i,j,k)
            csd(i,j,k) = sqrt(1.-cor(i,j)*cor(i,j)) * xsd(i,j,k+1)
          end do
        enddo
      enddo
      Return
      END SUBROUTINE ACINIT

#endif
