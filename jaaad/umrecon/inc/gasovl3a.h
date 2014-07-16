! GASOVL3A defines treatments of overlapping gasaeous absorption for
! two-stream radiation code.

      ! one species only
      INTEGER,PARAMETER:: IP_OVERLAP_SINGLE     = 1

      ! random overlap
      INTEGER,PARAMETER:: IP_OVERLAP_RANDOM     = 2

      ! fast esft
      INTEGER,PARAMETER:: IP_OVERLAP_FESFT      = 3

      ! clear-sky fast esft
      INTEGER,PARAMETER:: IP_OVERLAP_CLR_FESFT  = 4

      ! equivalent extinction
      INTEGER,PARAMETER:: IP_OVERLAP_K_EQV      = 5

      ! interpolated treatment for principal species
      INTEGER,PARAMETER:: IP_OVERLAP_SINGLE_INT = 6

! GASOVL3A end
