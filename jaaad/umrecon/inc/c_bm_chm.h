! C_BM_CHM start
! Contains constants required for biomass smoke conversion and
! diffusional scavenging by cloud droplets.
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   5.5    05/02/03 Original code, based on parameters used in the soot
!                   scheme.                         Paul Davison.

!     air parcel lifetime in cloud
      REAL,PARAMETER:: CLOUDTAU = 1.08E4            ! secs (=3 hours)

!     timescale for suspended aerosol to evaporate
      REAL,PARAMETER:: EVAPTAU = 300.0              ! secs  (=5 mins)

!     timescale for accumulation mode particles
      REAL,PARAMETER:: NUCTAU = 30.0                ! secs

!     Cloud liquid water threshold for nucleation scavenging to occur.
      REAL,PARAMETER:: THOLD = 1.0E-8               ! kg/kg

! C_BM_CHM end
