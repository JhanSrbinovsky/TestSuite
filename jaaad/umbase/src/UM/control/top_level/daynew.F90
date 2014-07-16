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
      Function DAYNEW(at,bt,itd)
!---------------------------------------------------------------------
!     Arguments
!---------------------------------------------------------------------
      implicit none
      Integer itd               ! IN Dayno. relative to winter
                                !    solstice
      Real                                                              &
     & at,bt                    ! IN Constants for calculating annual
                                !    cycle
!---------------------------------------------------------------------
!     Local variables
!---------------------------------------------------------------------
      Real                                                              &
     &  arg                                                             &
                                ! Argument
     &  ,daynew                 ! SIN of argument
!
      arg = at * float(itd) + bt
      daynew = sin(arg)
      Return
      END FUNCTION DAYNEW
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
