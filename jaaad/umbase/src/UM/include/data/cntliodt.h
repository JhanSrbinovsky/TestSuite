!====================== COMDECK CNTLIODT ========================
! Description:
!
!     Defines the sector size for well-formed transfers on Cray
!     Research systems.  Disk addresses must start on a sector
!     boundary, and transfers must be a number of sectors.  Disk
!     word addresses start at 0.
!
!     On the T3E, well-formed transfers must also start on a
!     cache-line boundary in memory.
!
!   4.4    24/10/97  New deck       C.P. Jones
!   4.5    02/10/98  Increase from 512 to 2048. D. Robinson.
!
!
      DATA UM_SECTOR_SIZE/2048/
