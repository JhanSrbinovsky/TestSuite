!--------------------------------------------------------------------
!LCOMDECK COMMG
!L-------------
      REAL          DLAT,DLONG,XLATN,XLONGW
#if defined(GLOBAL)
      COMMON/COMMG/ DLAT,DLONG,XLATN,XLONGW
#else
      REAL          ELFPLAT,ELFPLON
      COMMON/COMMG/ DLAT,DLONG,XLATN,XLONGW,ELFPLAT,ELFPLON
#endif
!--------------------------------------------------------------------
