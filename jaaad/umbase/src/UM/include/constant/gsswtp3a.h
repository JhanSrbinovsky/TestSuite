! GSSWTP3A defines arrays of gaussian points and weights for two-stream
! radiation code.

      ! max. order of gaussian quadrature
      INTEGER,PARAMETER:: NPD_GAUSS_ORD=10

      ! points of gaussian integration
      REAL :: GAUSS_POINT(NPD_GAUSS_ORD, NPD_GAUSS_ORD)

      ! weights of gaussian integration
      REAL :: GAUSS_WEIGHT(NPD_GAUSS_ORD, NPD_GAUSS_ORD)

      ! NB values defined in GSSWTD3A

! GSSWTP3A end
