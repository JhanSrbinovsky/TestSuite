! SOLVER3A defines reference numbers for solvers for two-stream
! radiation code.

      ! pentadiagonal scheme
      INTEGER,PARAMETER:: IP_SOLVER_PENTADIAGONAL=1

      ! mixed column scheme using full endekadiagonal matrix
      INTEGER,PARAMETER:: IP_SOLVER_MIX_11=6

      ! mixed column scheme with approximate scattering
      INTEGER,PARAMETER:: IP_SOLVER_MIX_APP_SCAT=9

      ! direct mixed column scheme for full fluxes
      INTEGER,PARAMETER:: IP_SOLVER_MIX_DIRECT=11

      ! direct solver for a homogeneous column
      INTEGER,PARAMETER:: IP_SOLVER_HOMOGEN_DIRECT=13

      ! direct solver for triple column
      INTEGER,PARAMETER:: IP_SOLVER_TRIPLE=14

      ! direct solver for triple column approximating scattering
      INTEGER,PARAMETER:: IP_SOLVER_TRIPLE_APP_SCAT=15

      ! direct mixed column scheme for full fluxes (modified
      !   for correct treatment of shadowing by Robin Hogan)
      INTEGER,PARAMETER:: IP_SOLVER_MIX_DIRECT_HOGAN=16

      ! direct solver for triple column (modified for
      !   correct treatment of shadowing by Robin Hogan)
      INTEGER,PARAMETER:: IP_SOLVER_TRIPLE_HOGAN=17

! SOLVER3A end
