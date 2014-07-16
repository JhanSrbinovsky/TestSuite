! CPRINTST defines print status for standard output messages

      ! Minimum output, only essential messages
      INTEGER,PARAMETER :: PrStatus_Min    = 1

      ! Normal informative messages + warnings
      INTEGER,PARAMETER :: PrStatus_Normal = 2

      ! Operational status, all informative messages
      INTEGER,PARAMETER :: PrStatus_Oper   = 3

      ! All informative + extra diagnostic messages
      INTEGER,PARAMETER :: PrStatus_Diag   = 4

      INTEGER PrintStatus ! Control volume of standard output messages
      COMMON/PrintSt/PrintStatus

! CPRINTST end
