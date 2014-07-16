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
      Subroutine XNEW(x, xa, xb, row_length, rows, nlevs, xt)
!---------------------------------------------------------------------
!     Arguments
!---------------------------------------------------------------------
      Implicit none
      Integer                                                           &
     &  row_length,rows , nlevs  ! IN dimension of arrays,
                                 ! & model levels
      Real                                                              &
     &  x(row_length, rows,nlevs)                                       &
                                 ! OUT Mean or SD of forcing variable
                                !     at day relative to winter
                                !     solstice
     &  ,xa(row_length, rows,nlevs)                                     &
                                    ! IN Amplitude of seasonal variation
                                !     of forcing variable.
     &  ,xb(row_length, rows,nlevs)                                     &
                                    ! IN Mean of seasonal variation
                                !     of forcing variable.
     &  ,xt                     ! IN Sin of argument
!---------------------------------------------------------------------
!     Local variables
!---------------------------------------------------------------------
      Integer                                                           &
     &  i,j,k                     ! Loop counters

      do  k = 1, nlevs
        do j = 1, rows
          do  i = 1, row_length
            x(i,j,k) = xa(i,j,k) * xt + xb(i,j,k)
          enddo
        enddo
      enddo
      Return
      END SUBROUTINE XNEW
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
