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
      Subroutine ABNEW(x1, x2, xa, xb, row_length, rows, n)
      Implicit none
!---------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------
      Integer                                                           &
     &  i, j, k, row_length, rows, n         ! Loop counters
      Real                                                              &
     &  x1(row_length, rows,n)                                          &
                                ! IN SD or mean of forcing variable
                                !    for max. of annual cycle (July)
     & ,x2(row_length, rows,n)                                          &
                                ! IN SD or mean of forcing variable
                                !    for min. of annual cycle (Jan)
     &  ,xa(row_length, rows,n)                                         &
                                ! OUT Amplitude of seasonal variation
                                !    of forcing variable
     &  ,xb(row_length, rows,n) ! OUT Mean of seasonal variation
                                !    of forcing variable
!
      Do i = 1, row_length
        Do j = 1, rows
          Do k = 1, n
            xa(i,j,k) = (x1(i,j,k)-x2(i,j,k))/2.
            xb(i,j,k) = (x1(i,j,k)+x2(i,j,k))/2.
          enddo
        enddo
      enddo
      Return
      END SUBROUTINE ABNEW
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

#endif
