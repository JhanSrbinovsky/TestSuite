#if defined(A05_4A) || defined(A05_5A)
! CRITDEP start

      ! critical depth of cloud for the formation of
      ! convective precipitation over sea (m)
      REAL,PARAMETER:: CRITDSEA = 1.5E3

      ! critical depth of cloud for the formation of convective
      ! precipitation over land (m)
      REAL,PARAMETER:: CRITDLND = 4.0E3

      ! critical depth of a glaciated cloud for the formation of
      ! convective precipitation (m)
      REAL,PARAMETER:: CRITDICE = 1.0E3

! CRITDEP end
#endif
