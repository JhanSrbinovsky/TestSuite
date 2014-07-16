!
! Description:
!  Comdeck contains parameters for horizontal interpolation routines
!
! Current Code Owner: D.M.Goddard
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   4.0  18/04/95   Original code. D.M. Goddard
!
! Declarations:

! Global parameters:
      INTEGER      BGRID            !I_GRID_IN/OUT value for B_GRID
      INTEGER      CGRID            !I_GRID_IN/OUT value for C_GRID
      INTEGER      IDLAT            !Posn of lat increment in REALHD
      INTEGER      IDLON            !Posn of lon increment in REALHD
      INTEGER      IPLAT            !Posn of N pole latitide in REALHD
      INTEGER      IPLON            !Posn of N pole longitude in REALHD
      INTEGER      ISLAT            !Posn of start latitide in REALHD
      INTEGER      ISLON            !Posn of start longitude in REALHD
      INTEGER      ISTAG            !Posn of grid staggering in FIXHD
      INTEGER      ITYPE            !Posn of domain type in FIXHD
      REAL         SMALL            !Used in IF test to check that
                                    !  target is within source area.

      PARAMETER(BGRID=1)
      PARAMETER(CGRID=2)
      PARAMETER(IDLAT=2)
      PARAMETER(IDLON=1)
      PARAMETER(IPLAT=5)
      PARAMETER(IPLON=6)
      PARAMETER(ISLAT=3)
      PARAMETER(ISLON=4)
      PARAMETER(ISTAG=9)
      PARAMETER(ITYPE=4)
      PARAMETER(SMALL=0.001)


!- End of COMDECK declaration
