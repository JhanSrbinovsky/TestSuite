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
      Subroutine PNEW(nlevs, p, rp, points, n, pstar, ak, bk)
!---------------------------------------------------------------------
!     Arguments
!---------------------------------------------------------------------
      implicit none
      Integer                                                           &
     &  nlevs                                                           &
                                ! IN no. of levs of the scm.
     &  ,n                                                              &
                                ! IN no. of levels to be processed
     &  ,points                 ! IN no. of model columns.
      Real                                                              &
     &  ak(nlevs)                                                       &
     &  ,bk(nlevs)                                                      &
                                ! IN AK and BK values at levels
     &  ,p(points,n)                                                    &
                                ! OUT Pressure coordinates (Pa)
     &  ,pstar(points)                                                  &
                                ! IN Surface pressure (Pa)
     &  ,rp(points,n)           ! OUT Reciprocal pressure
                                ! coordinates (HPa)
!---------------------------------------------------------------------
!     Local variables
!---------------------------------------------------------------------
      Integer                                                           &
     &  i,k                     ! Loop counters

!
      Do  k = 1, n
        Do  i = 1, points
          p(i,k) = ak(k) + bk(k) * pstar(i)
          rp(i,k) = 100. / p(i,k)
        enddo
      enddo
      Return
      END SUBROUTINE PNEW
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
