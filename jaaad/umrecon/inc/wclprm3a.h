! WCLPRM3A defines numbers for water cloud schemes in two-stream
! radiation code.
      ! number of cloud fitting schemes
      INTEGER,PARAMETER:: NPD_CLOUD_FIT=3

      ! parametrization of slingo-schrecker
      INTEGER,PARAMETER:: IP_SLINGO_SCHRECKER=1

      ! parametrization of ackerman & stephens
      INTEGER,PARAMETER:: IP_ACKERMAN_STEPHENS=2

      ! unparametrized droplet data
      INTEGER,PARAMETER:: IP_DROP_UNPARAMETRIZED=3

      ! pade approximation of the second order (third order for the
      ! extinction)
      INTEGER,PARAMETER:: IP_DROP_PADE_2=5
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      ! parametrization of slingo-schrecker with phase func.
      INTEGER,PARAMETER:: IP_SLINGO_SCHR_PHF=6
#endif
! WCLPRM3A end
