!
!  4.1      04/09/96 New comdeck for Data Assimilation on T3E
!                                   Deborah Salmond
!  4.3      18/4/97 Increase OBSNUMDIM  Stuart Bell
!  4.4      19/11/97 Increase OBSNUMDIM,OBSDIM,inxdim  Stuart Bell
!  5.2      12/12/00 Restore COUNTA,B,C dropped at 5.0   B Macpherson
!LL  5.0      24/06/99 Alter names for C-P C-grid dynamics. M.L.Gallani
!LL  5.3      08/06/01 Remove duplicate declarations.  A van der Wal
!    5.3      05/12/01 Re-size obsdim to smaller value (500,000 to
!                      250,000) to reduce memory usage.    S. Cusack
!    6.0     17/10/03 Increase obsdim and obsnumdim for SX6. Clive Jones
!    6.2     11/01/06 Remove hard-wired obsdim and obsnumdim.
!                       Camilla Mathison/R Barnes

      ! dimensions for obs allocated in subroutine AC
      ! dynamically allocated in AC and RDOBS

      ! dimension inxdim allocated in subroutine HORINF
      INTEGER,PARAMETER:: inxdim    = 15000

      ! common for Statistics Calcs in DIAGO ; Prints in RDOBS,GETOBS
      REAL :: R_STAT(MODEL_LEVELS_MAX,0:8)
      REAL :: S_STAT(MODEL_LEVELS_MAX,0:8)
      INTEGER :: COUNTA(NOBTYPMX)
      INTEGER :: COUNTB(NOBTYPMX)
      INTEGER :: COUNTC(NOBTYPMX)

      COMMON /mpp_ac/ R_STAT,S_STAT,COUNTA,COUNTB,COUNTC

      ! common to pass longitudes and latitudes for edges of local
      ! box from setcona to RDOBS and HINTCF

      REAL :: LONG_E
      REAL :: LONG_W
      REAL :: LAT_N,LAT_S
      REAL :: LONG_W_MODEL
      REAL :: LONG_E_MODEL
      COMMON/latlonmax/                                                 &
     &  LONG_E,LONG_W,LAT_N,LAT_S,LONG_W_MODEL,LONG_E_MODEL
! MPPAC end
