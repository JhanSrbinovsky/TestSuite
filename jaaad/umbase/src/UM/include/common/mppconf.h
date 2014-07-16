! MPPCONF defines MPP configuration
!
! Owner: P. Selwood
!

      INTEGER :: extended_halo_size_EW
      INTEGER :: extended_halo_size_NS
      INTEGER :: gcom_coll_limit

      NAMELIST / NLST_MPP / extended_halo_size_EW,            &
                            extended_halo_size_NS,            &
                            gcom_coll_limit

      COMMON / COM_MPPCONF / extended_halo_size_EW,           &
                             extended_halo_size_NS,           &
                             gcom_coll_limit

! MPPCONF end
