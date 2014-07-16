!*L------------------ COMDECK CINTF -----------------------------------

!*L------------------ COMDECK CINTF  ----------------------------------
!
!   Contains Variables, Headers and Index blocks for control of
!   generation of boundary information for limited area models
!   (corresponds to CINTFA and CINTFO)


      INTEGER                                                           &
     &  INTF_ROW_LENGTH(N_INTF)                                         &
                                 ! row length  (interface field)
     & ,INTF_P_ROWS(N_INTF)                                             &
                                 ! no of rows          "
     & ,INTF_P_LEVELS(N_INTF)                                           &
                                 ! no of levels        "
     & ,INTF_Q_LEVELS(N_INTF)                                           &
                                 ! no of wet levels    "
     & ,INTF_TR_LEVELS(N_INTF)                                          &
                                 ! no of tracer levels "
     & ,INTFWIDTH(N_INTF)                                               &
                                 ! Width of interface zone
     & ,INTF_START_HR(N_INTF)                                           &
                               ! ) Start, Frequency and End time
     & ,INTF_FREQ_HR(N_INTF)                                            &
                               ! ) in hours for which  interface
     & ,INTF_END_HR(N_INTF)                                             &
                               ! ) data is to be generated.
     & ,LEN_INTF_P(N_INTF)                                              &
                                ! Length of interface p field
     & ,LEN_INTF_U(N_INTF)                                              &
                                ! Length of interface u field
     & ,LEN_INTF_DATA(N_INTF)                                           &
                                ! Length of interface data
     & ,INTF_PACK(N_INTF)                                               &
                                ! Packing Indicator for data
     & ,INTF_NTIMES(N_INTF)                                             &
                                ! number of times file to be output
     & ,INTF_LOOKUPS_OUT(N_INTF) ! number of field types to be output

      REAL                                                              &
     &  INTF_EWSPACE(N_INTF)                                            &
                                 ! E-W grid spacing (degrees)
     & ,INTF_NSSPACE(N_INTF)                                            &
                                 ! N-S grid spacing (degrees)
     & ,INTF_FIRSTLAT(N_INTF)                                           &
                                 ! Latitude of first row (degrees)
     & ,INTF_FIRSTLONG(N_INTF)                                          &
                                 ! Longitude of first row (degrees)
     & ,INTF_POLELAT(N_INTF)                                            &
                                 !  latitude of coordinate pole
     & ,INTF_POLELONG(N_INTF)                                           &
                                 !  longitude of coordinate pole
     & ,INTF_AKH(MAX_INTF_P_LEVELS+1,N_INTF)                            &
                                             ! A and B hybrid co-ords
     & ,INTF_BKH(MAX_INTF_P_LEVELS+1,N_INTF)                            &
                                            ! at  half levels
     & ,INTF_AK(MAX_INTF_P_LEVELS  ,N_INTF)                             &
                                            ! A and B hybrid co-ords
     & ,INTF_BK(MAX_INTF_P_LEVELS  ,N_INTF) ! at full levels

      LOGICAL                                                           &
     &  INTF_VERT_INTERP(N_INTF)                                        &
                                 ! T => do vertical interpoln
     & ,LNEWBND(N_INTF)                                                 &
                           ! T => initialise new bdy data file
     & ,ll_intf_n(N_INTF),                                              &
                           ! T=> corresponding         north
     &  ll_intf_e(N_INTF),                                              &
                           !     boundary data         south
     &  ll_intf_s(N_INTF),                                              &
                           !       to be               east
     &  ll_intf_w(N_INTF)                                               &
                           !      generated            west
     & ,ll_intf_tracer(N_INTF),                                         &
                                ! T=> corresponding    tracers
     &  ll_intf_uv(N_INTF),                                             &
                                !     boundary data    baroclinic uv
     &  ll_intf_stream(N_INTF),                                         &
                                !       to be          streamfunction
     &  ll_intf_ice(N_INTF),                                            &
                                !      generated       ice fields
     &  ll_intf_ssh(N_INTF),                                            &
                                !                      rigid-lid p/ssh
     &  ll_intf_sfuv(N_INTF)    !                      barotropic uv

!---------------------------------------------------------------------
