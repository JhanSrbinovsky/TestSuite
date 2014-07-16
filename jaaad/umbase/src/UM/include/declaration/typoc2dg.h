! TYPOC2DG Defines 2D variables which are to be passed through all
! subroutine argument lists from ROW_CTL to ROWCALC.
      ! Net heat (not penetrating) into ocean through leads.
      REAL :: oceanheatflux(imt,jmt)

      REAL :: oceansnowrate(imt,jmt) ! Snowfall melting in ocean.

      ! Heat/W m**-2 lost at the ocean floor by mixing downwards after
      ! ice physics
      REAL :: heatsink(imt,jmt)

      ! Rate of change of salinity/s**-1 caused by surface water flux,
      ! calculated as surface salinity * water flux / density
      REAL :: diagsw(imt,jmt)

      ! tmin1, athkdftu, athkdftv are used in the Visbeck scheme
      ! Timescale, isopyc. thickness diffusion
      REAL :: tmin1(imt_vis,jmt_vis)

      ! thickness diffusion coeff. (u* pts)
      REAL :: athkdftu(imt_vis,jmt_vis)

      ! thickness diffusion coeff. (v* pts)
      REAL :: athkdftv(imt_vis,jmt_vis)
! TYPOC2DG end
