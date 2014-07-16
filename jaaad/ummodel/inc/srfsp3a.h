! SRFSP3A defines permitted methods of specifying the  surface albedo and
! emissivity for two-stream radiation code.

      ! properties specified by surface type
      INTEGER,PARAMETER:: IP_SURFACE_SPECIFIED=1

      ! properties passed into code
      INTEGER,PARAMETER:: IP_SURFACE_INTERNAL=2

      ! direct albedo fitted as polynomial
      INTEGER,PARAMETER:: IP_SURFACE_POLYNOMIAL=3

      ! fit in the functional form used by payne
      INTEGER,PARAMETER:: IP_SURFACE_PAYNE=4
! SRFSP3A end
