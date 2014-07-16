!     ------------------------------------------------------------------
!     Module to set unit numbers for standard I/O.

      INTEGER, PARAMETER :: IU_STDIN = 5
!           Unit number for standard input
      INTEGER, PARAMETER :: IU_STDOUT = 6
!           Unit number for standard output
      INTEGER, PARAMETER :: IU_ERR = 6
!           Unit number for error messages

#if defined(STANDARD)
      INTEGER, PARAMETER :: IU_USER = 6
!           Unit number for prompts to user
      INTEGER, PARAMETER :: IU_MONITOR = 80
!           Unit number for monitoring output
      INTEGER, PARAMETER :: IU_MONITOR_2 = 81
!           Unit number for monitoring output
      INTEGER, PARAMETER :: IU_FILE_IN = 3
!           Unit number for input from a file
      INTEGER, PARAMETER :: IU_FILE_OUTPUT = 4
!           Unit number for output to a file
      INTEGER, PARAMETER :: IP_UNIT_USER_LOW = 10
!           Lowest unit number for user's input
      INTEGER, PARAMETER :: IP_UNIT_USER_HIGH = 99
!           Highest unit number for user's input
#endif
!     ------------------------------------------------------------------
