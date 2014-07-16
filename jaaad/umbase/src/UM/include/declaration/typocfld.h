! TYPOCFLD
!  4.5  14/08/97  Removed the boundary arrays (TBOUND_N,...)
!                 C.G. Jones
!  5.2  31/07/00  Included fkmz. M J Bell
!  5.4  29/08/02  Included fkmqx. D. Storkey
      REAL :: HR (IMT,JMT) ! Reciprocal of depth at u,v points

      ! HR for a particular row which is outside the scope of halos in
      ! mpp code.
      REAL :: HRJ(IMT)

      REAL :: HRJP(IMT)     ! HR for J+2 outside halo
      REAL :: PJP(IMT)      ! P for J+3 outside halo
      REAL :: PBJP(IMT)     ! PB for J+3 outside halo
      REAL :: FKMQJP(IMT)     ! FKMQ for J+2 outside halo
      REAL :: UBTBBCJP(IMT)
      REAL :: VBTBBCJP(IMT)
      REAL :: UBTJP(IMT)
      REAL :: VBTJP(IMT)
      REAL :: FKMP(IMT,JMT) ! No. of vertical levels at T points
      REAL :: FKMQ(IMT,JMT) ! No, of vertical levels at u,v points
      REAL :: FKMQX(IMT,JMT) ! Equal to FKMQ but not masking straits
      REAL :: FKMZ(IMT,JMT) ! No. of vertical levels at interior points

      ! No of vert levs at u,v pnts(full)
      REAL :: FKMQ_GLOBAL(IMT,JMT_GLOBAL)

      REAL :: CORIOLIS(IMT,JMT) ! 2*OMEGA*SIN(lat) on rotated grid
      REAL ::  EMU(IMT,JMT)
      REAL :: EM (IMT,JMT)!used in delplus/delcross smoothing
! TYPOCFLD end
