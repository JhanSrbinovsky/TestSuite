! TYPOCFLW start
      INTEGER:: KAR(KM)           ! KAR(K)=K for vectorization
      INTEGER:: ISZ(JMT,LSEG)     ! Start points for streamfuction calc
      INTEGER:: IEZ(JMT,LSEG)     ! End points for streamfunction calc
      INTEGER:: ISE(JMT,LSEGFS)   ! start indicies for eta
      INTEGER:: IEE(JMT,LSEGFS)   ! stop indicies for eta
      INTEGER:: ISU(JMTM1,LSEGFS) ! st indicies for barotropic velys
      INTEGER:: IEU(JMTM1,LSEGFS) ! stop indicies for batotropic velys
      INTEGER:: LSE(JMT)       ! no of st/stop indicies for eta
      INTEGER:: LSU(JMTM1)     ! no of st/stop indicies for btpic velys
      INTEGER:: ISEG(NISLE)       ! Number of segments in island
      INTEGER:: ISIS(NISLE,ISEGM) ! Start of segment (col)
      INTEGER:: IEIS(NISLE,ISEGM) ! End of segment   (col)
      INTEGER:: JSIS(NISLE,ISEGM) ! Start of segment (row)
      INTEGER:: JEIS(NISLE,ISEGM) ! End of segment   (row)
      INTEGER :: CISLBDY(NISLE)   ! Islands on bdy: control array
      INTEGER :: IISLBDY(NISLE)   ! Isl on bdy: i-coord of sampling point
      INTEGER :: JISLBDY(NISLE)   ! Isl on bdy: j-coord of sampling point
      INTEGER :: PE_ISLBDY(NISLE) ! Isl on bdy: PE owner of sampling pt
      REAL :: t_strait(km,nt,n_strait,4) ! tracers stored at Strait pts
      REAL :: flux_strait(km,n_strait)   ! flux (cm^3/s) across Strait
                                     ! points FROM point 1 TO point 2
      REAL :: t_strait_clm(km,nt,12,n_strait_clm) ! climate tracer
                                     ! values for Straits exchange

      INTEGER:: ISTF(NJTBFT,LSEGF,KM) ! Start index for filtering(T)
      INTEGER:: IETF(NJTBFT,LSEGF,KM) ! End index for filtering (T)
      INTEGER:: ISUF(NJTBFU,LSEGF,KM) ! Start index for filtering (U)
      INTEGER:: IEUF(NJTBFU,LSEGF,KM) ! End index for filtering (U)
      INTEGER:: ISZF(NJTBFU,LSEGF)    ! Start index for filtering (vort)
      INTEGER:: IEZF(NJTBFU,LSEGF)    ! End index for filtering (vort)

      REAL :: SPSIN(IMT)          ! Turning coefficients for U filtering
      REAL :: SPCOS(IMT)          !    --------"--------------------

      ! Roussenov convective adj scheme
      INTEGER :: ISROUS   ! start i index
      INTEGER :: IEROUS   ! end i index
      INTEGER :: JSROUS   ! start j index
      INTEGER :: JEROUS   ! end j index
#include "comocflw.h"
#include "comocisl.h"
! TYPOCFLW end
