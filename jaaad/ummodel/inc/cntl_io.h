!====================== COMDECK CNTL_IO ========================
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
!   4.3    30/04/97  New deck       B. Carruthers, Cray Research
!   4.4    27/10/97  Remove DATA statement. C.P. Jones
!   5.1    07/06/00  Upon VAR request, provide alternative to the
!                    common statement for um_sector_size.
!                    JC Thil
!
!
      INTEGER UM_SECTOR_SIZE    ! Sector size on disk for I/O
!
#if !defined(VAROPSVER)
      COMMON / CNTL_IO / UM_SECTOR_SIZE
#else
      Parameter(UM_SECTOR_SIZE = 2048)
#endif
