
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      SUBROUTINE DD_ALL_CALL (npnts,npossdd,klev,nlev,trlev,ntra        &
     &,                      kterm, iccb, icct, index1, l_tracer        &
     &,                      flg_dwn_flx, flg_entr_dwn, flg_detr_dwn    &
     &,                      bwater                                     &
     &,                      exner_layer_centres,exner_layer_boundaries &
     &,                      p_layer_centres, p_layer_boundaries,pstar  &
     &,                      recip_pstar,timestep , cca                 &
     &,                      thp, qp, the, qe, trap,trae, flx, precip   &
     &,                      dthbydt, dqbydt, dtrabydt                  &
     &,                      rain, snow ,rain_3d, snow_3d, dd_flux      &
     &,                      entrain_dwn, detrain_dwn)

!  Purpose : Calculate initial downdraught massflux.
!            Reset en/detrainment rates for downdraught
!            Compress/expand variables
!            Initialise downdrought
!            Call downdraught routine
!
!  Suitable for single column model use
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!

      IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!
      integer, intent(in) :: npnts     ! Vector length of full arrays

      integer, intent(in) :: npossdd   ! Length of array holding index
                                       ! of points where a downdraught
                                       ! is possible.
      integer, intent(in) :: nlev      ! Number of model levels
      integer, intent(in) :: trlev     ! Number of tracer levels
      integer, intent(in) :: ntra      ! Number of tracers
      integer, intent(in) :: klev     ! Number of levels (may be model
                                      ! levels or a reduced set required
                                      ! here equal to kmax_term+1).

      integer, intent(in) :: kterm(npnts)
                                       ! Convective cloud top layer

      integer, intent(in) :: iccb(npnts) ! Convective cloud base
                                         ! level (m)
      integer, intent(in) :: icct(npnts) ! Convective cloud top
                                         ! level (m)
      integer, intent(in) :: index1(npnts) ! index of points where
                                         ! downdraught possible
      logical, intent(in) :: l_tracer    ! Switch for tracers
      logical, intent(in) :: flg_dwn_flx ! stash flag for
                                         ! downdraught mass flux
      logical, intent(in) :: flg_entr_dwn ! stash flag for
                                          ! downdraught entrainment
      logical, intent(in) :: flg_detr_dwn ! stash flag for
                                          ! downdraught detrainment
      logical, intent(in) :: bwater(npnts,2:klev+1)
                                      ! Mask for points at which
                                      ! condensate is liquid

      real , intent(in) :: exner_layer_centres(npnts,0:klev+1) !exner

      real , intent(in) :: exner_layer_boundaries(npnts,0:klev+1)
                                      ! exner at half level above
                                      ! exner_layer_centres

      real , intent(in) :: p_layer_centres(npnts,0:klev+2)
                                      ! Pressure (Pa)

      real , intent(in) :: p_layer_boundaries(npnts,0:klev+1)
                                      ! Pressure at half level above
                                      ! p_layer_centres (Pa)

      real , intent(in) :: pstar(npnts)       ! Surface pressure (Pa)

      real , intent(in) :: recip_pstar(npnts) ! 1/pstar (Pa)

      real , intent(in) :: timestep    ! timestep
      real , intent(in) :: CCA(npnts)  ! 2d convective cloud amount

      real , intent(in) :: thp(npnts,klev+1)
                                 ! Parcel potential temperature (K)

      real , intent(in) :: qp(npnts,klev+1)
                                 ! Parcel mixing ratio (KG/KG)

      real , intent(in) :: the(npnts,klev+1)
                                 ! Model enviromental potential
                                 ! temperature (K)
      real , intent(in) :: qe(npnts,klev+1)
                                 ! Model enviromental mixing ratio
                                 ! (KG/KG)
      real , intent(in) :: trap(npnts,nlev,ntra)
                                 ! Parcel tracer (kg/kg)
      real , intent(in) :: trae(npnts,trlev,ntra)
                                 ! Environment tracer (kg/kg)
      real , intent(in) :: flx(npnts,klev+1)
                                 ! updraught mass flux (Pa/s)

!
! Arguments with intent INOUT:
!
      real , intent(inout) :: precip(npnts,klev+1)
                                 ! precipitation added when descending
                                 ! from layer k to k-1 (KG/M**2/S)
      real , intent(inout) :: dthbydt(npnts,klev+1)
                                 ! increment to model potential
                                 ! temperature (K/S)
      real , intent(inout) :: dqbydt(npnts,klev+1)
                                 ! increment to model mixing ratio
                                 ! (KG/KG/S)
      real , intent(inout) :: dtrabydt(npnts,nlev,ntra)
                                 ! increment to model tracers (KG/KG/S)

!
! Arguments with intent OUT:
!
      real , intent(inout) :: rain(npnts)
                                 ! rainfall at surface (KG/M**2/S)

      real , intent(inout) :: snow(npnts)
                                 ! snowfall at surface (KG/M**2/S)

      real , intent(inout) :: rain_3d(npnts,klev+1)
                                 ! rainfall flux (KG/M**2/S)

      real , intent(inout) :: snow_3d(npnts,klev+1)
                                 ! snowfall flux (KG/M**2/S)

      real , intent(out) :: dd_flux(npnts,klev+1)
                                       ! downdraught mass flux
      real , intent(out) :: entrain_dwn(npnts,klev+1)
                                       ! fractional entrainment rate for
                                       ! downdraught mass flux
      real , intent(out) :: detrain_dwn(npnts,klev+1)
                                       ! fractional detrainment rate for
                                       ! downdraught mass flux
!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------


      integer :: i,k ,i2,ktra            ! Loop counters

      integer :: ndd           ! Compressed vector length for
                               ! downdraught calculation
      integer :: nddon_tmp     ! Number of points with active
                               ! downdraught

!-----------------------------------------------------------------------
! Compressed arrays
!-----------------------------------------------------------------------
                             ! Compressed exner function at layer
      real :: exner_km12_c(npossdd)   ! k
      real :: exner_kp12_C(npossdd)   ! k+1
      real :: exner_KM32_C(npossdd)   ! k-1

      real :: pk(npossdd)             ! Pressure of layer k (PA)
      real :: p_km1(npossdd)          ! Pressure of layer K-1 (PA)

      real :: exk(npossdd)            ! exner ratio for layer K
      real :: exkm1(npossdd)          ! exner ratio for layer K-1

      real :: delpk(npossdd)       ! Pressure difference across layer K
                                   ! (PA)
      real :: delpkm1(npossdd)     ! Pressure difference across layer
                                   ! K-1 (PA)

      real :: amdetk(npossdd)      ! Mixing detrainment at level k
                                   ! multiplied by appropriate layer
                                   ! thickness
      real :: ekm14(npossdd)       ! exner ratio at layer k-1/4
      real :: ekm34(npossdd)       ! exner ratio at layer k-3/4

      real :: precip_k_c(npossdd)  ! Compressed precipitation added
                                   ! when descending from layer
                                   ! K to K-1 (KG/M**2/S)

      real :: q_k_c(npossdd)       ! Compressed parcel mixing ratio of
                                   ! layer K (KG/KG)

      real :: th_k_c(npossdd)      ! Compressed parcel potential
                                   ! temperature of layer k (K)

      real :: tra_k_c(npossdd,ntra)  ! Compressed parcel tracer in
                                     ! layer K (KG/KG)

      real :: pstar_c(npossdd)       ! Compressed surface pressure (PA)

      real :: recip_pstar_c(npossdd) ! Reciprocal of comp. pstar array

      real :: P_layer_centres_C(npossdd,0:klev+2)
      real :: P_layer_boundaries_C(npossdd,0:klev+1)
      real :: exner_layer_centres_C(npossdd,0:klev+1)

      real :: dthbydt_K_C(npossdd)  ! Compressed increment to model
                                    ! potential temperature of layer k
                                    ! (K/S)
      real :: dthbydt_km1_c(npossdd)  ! Compressed increment to model
                                    ! potential temperature of layer k-1
                                    ! (K/S)

      real :: dqbydt_k_c(npossdd)   ! Compressed increment to model
                                    ! mixing ratio of layer k (KG/KG/S)
      real :: dqbydt_km1_c(npossdd) ! Compressed increment to model
                                    ! mixing ratio of layer k-1(KG/KG/S)

      real :: dtra_k_c(npossdd,ntra) ! Compressed increment to model
                                     ! tracer of layer k (KG/KG/S)
      real :: dtra_km1_c(npossdd,ntra) ! Compressed increment to model
                                       ! tracer of layer k-1 (KG/KG/S)

      real :: deltd(npossdd)       ! Cooling necessary to achieve
                                   ! saturation (K)

      real :: delqd(npossdd)       ! moistening necessary to achieve
                                   ! saturation (K)

      real :: deltrad(npossdd,ntra) ! Depletion of environment tracer
                                  ! due to downdraught formation (KG/KG)

      real :: qdd_k(npossdd)       ! mixing ratio of downdraught in
                                   ! layer K (KG/KG)

      real :: thdd_k(npossdd)      ! Model potential temperature of
                                   ! downdraught in layer k (K)

      real :: tradd_k(npossdd,ntra) ! Model tracer of downdraught in
                                    ! layer K (KG/KG)

      real :: flx_dd_k(npnts)     ! Downdraught initial mass flux (Pa/s)
      real :: flx_dd_k_c(npossdd) ! Compressed downdraught initial mass
                                  ! flux (Pa/s)

      integer :: ICCB_c(npossdd)   ! Compressed cloud base level

      logical :: bwater_k_c(npossdd) ! Compressed mask for those points
                                     ! at which condensate is water in
                                     ! layer k
      logical :: bddi(npnts)        ! Mask for points where downdraught
                                    ! might occur
      logical :: bddi_c(npossdd)    ! Compressed mask for points where
                                    ! downdraught may initiate

      real :: qe_k_c(npossdd)       ! Compressed environment mixing
                                    ! ratio of layer k (KG/KG)
      real :: qe_km1_c(npossdd)     ! Compressed environment mixing
                                    ! ratio of layer k-1 (KG/KG)

      real :: the_k_c(npossdd)      ! Compressed potential temperature
                                    ! of environment in layer k (K)
      real :: the_km1_c(npossdd)    ! Compressed potential temperature
                                    ! of environment in layer k-1 (K)

      real :: trae_k_c(npossdd,ntra)  ! Compressed tracer of environment
                                      ! in layer k (KG/KG)
      real :: trae_km1_c(npossdd,ntra)! Compressed tracer of environment
                                      ! in layer k-1 (KG/KG)

      real :: rain_c(npossdd)  ! Compressed surface rainfall (kg/m**2/s)

      real :: snow_c(npossdd)  ! Compressed surface snowfall (kg/m**2/s)

      real :: flx_ud_k_c(npossdd)   ! updraught mass flux at layer K

      real :: rain_env(npossdd)     ! Amount of rainfall passing through
                                    ! environment (kg/m**2/s)

      real :: snow_env(npossdd)     ! Amount of snowfall passing through
                                    ! environment (kg/m**2/s)

      real :: rain_dd(npossdd)      ! Amount of rainfall passing through
                                    ! downdraught (KG/M**2/S)

      real :: snow_dd(npossdd)      ! Amount of snowfall passing through
                                    ! downdraught (KG/M**2/S)

      logical :: dbb_start(npnts)   ! mask for those points where
                                    ! downdraught is able to start
                                    ! from level k

      logical :: dbb_start_C(npossdd) ! Compressed mask for those points
                                    ! where downdraught is able to start
                                    ! from level k

      logical :: bddwt_k(npnts)     ! mask for points in downdraught
                                    ! where ppt in layer k is liquid

      logical :: bddwt_k_C(npossdd) ! Compressed mask for points in DD
                                    ! where ppt in layer k is liquid

      logical :: bddwt_kM1(npnts)   ! mask for points in downdraught
                                    ! where ppt in layer k-1 is liquid

      logical :: bddwt_km1_c(npossdd) ! Compressed mask for points in DD
                                    ! where ppt in layer k-1 is liquid

      logical :: bdd_on(npnts)      ! mask for those points where DD
                                    ! continues from layer K+1

      logical :: bdd_on_c(npossdd)  ! Compressed mask for points where
                                    ! DD continues from layer K+1

      integer :: kmin(npossdd)      ! Freezing level where Entrainment
                                    ! rates are increased

      real :: flx_strt(npnts)       ! Mass flux at level where
                                    ! downdraught starts (Pa/s)

      real :: flx_strt_c(npossdd)   ! Compressed value of flx_strt

      real :: CCA_C(npossdd)        ! Compressed convective cloud amount

      real :: Lr_ud_ref(npossdd)    ! precipitation mixing ratio at
                                    ! lowest precipitationg level of UD

!-----------------------------------------------------------------------
! Model constants
!-----------------------------------------------------------------------

!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------

!-----------------------------------------------------------------------
! EXTERNAL ROUTINES CALLED
!-----------------------------------------------------------------------

      EXTERNAL FLX_INIT, LAYER_DD, DD_INIT, DOWND


!-----------------------------------------------------------------------
! Compression to downdraught points  (all levels even above term level)
!-----------------------------------------------------------------------

      NDD=npossdd

      Do K=0,klev+1
        Do i=1,NDD
          p_layer_centres_c(i,k)     = p_layer_centres(index1(i),k)
          exner_layer_centres_c(i,k) = exner_layer_centres(index1(i),k)
          p_layer_boundaries_c(i,k)  = p_layer_boundaries(index1(i),k)
        End do
      End do

!----------------------------------------------------------------------
! INITIALISE logical ARRAYS AS FALSE
!-----------------------------------------------------------------------

      Do i=1,npnts
        bddi(i)      = .false.
        dbb_start(i) = .false.
        bddwt_k(i)   = .false.
        bddwt_km1(i) = .false.
        bdd_on(i)    = .false.
      End do
      Do i2=1,ndd
        bddi(INDEX1(I2)) = .TRUE.    ! all points with downdraught
      End do

!----------------------------------------------------------------------
! Calculate initial downdraught mass flux
! Done on all npnts but uses test on bddi (klev not used).
!-----------------------------------------------------------------------

! DEPENDS ON: flx_init
      call flx_init(npnts, klev, ICCB, ICCT, flx, flx_dd_k, bddi,       &
                    flx_strt)

!-----------------------------------------------------------------------
! Compress all input arrays for the downdraught calculation down to
! those points where a downdraught is possible when the convection
! terminates in the column.
!-----------------------------------------------------------------------
! Initialise various arrays before level loop

       Do i=1,NDD
         flx_dd_k_C(i) = flx_dd_k(index1(i))
         flx_strt_C(i) = flx_strt(index1(i))
         pstar_c(i) = pstar(index1(i))
         recip_pstar_c(i)=recip_pstar(index1(i))
         ICCB_C(i) = ICCB(index1(i))
         dbb_start_C(i) = dbb_start(index1(i))

         rain_c(i) = rain(index1(i))
         snow_c(i) = snow(index1(i))

         bddwt_k_C(i) = bddwt_k(index1(i))
         bddwt_km1_c(i) = bddwt_km1(index1(i))
         bdd_on_C(i) = bdd_on(index1(i))

         CCA_C(i) = CCA(index1(i))
         LR_UD_REF(i) = 0.0

! Initialise compressed downdraught indicator array to false

         bddi_c(i) = .false.
       End do


!-----------------------------------------------------------------------
! Main level loop working from top downwards
!-----------------------------------------------------------------------

      Do K = Klev+1,2,-1

! Need at this stage to reset bddi_c to true if reached level where
! convection terminated

        Do I=1,NDD
          If (kterm(index1(i))+1 == k) then
            bddi_c(i) = .true.
          End if

! Compress arrays to those points with downdraughts possible in the
! column
          th_k_c(i) = thp(index1(i),K)
          q_k_c(i)  = qp(index1(i),K)
          the_k_c(i)   = the(index1(i),K)
          the_km1_c(i) = the(index1(i),K-1)
          qe_k_c(i)   = qe(index1(i),K)
          qe_km1_c(i) = qe(index1(i),K-1)
          dthbydt_K_C(i)   = dthbydt(index1(i),K)
          dthbydt_km1_c(i) = dthbydt(index1(i),K-1)
          dqbydt_K_C(i)   = dqbydt(index1(i),K)
          dqbydt_km1_c(i) = dqbydt(index1(i),K-1)
          exner_KM12_C(i) = exner_layer_boundaries(index1(i),K-1)
          exner_KP12_C(i) = exner_layer_boundaries(index1(i),K)
          exner_KM32_C(i) = exner_layer_boundaries(index1(i),K-2)
          precip_K_C(i) = precip(index1(i),K)
          flx_ud_k_c(i) = flx(index1(i),K)
          bwater_k_c(i) = bwater(index1(i),K)
        End do

        If (L_TRACER) then

          Do ktra=1,ntra
             Do i=1,NDD
               tra_k_c(I,ktra)    = trap(index1(i),K,ktra)
               trae_k_c(I,ktra)   = trae(index1(i),K,ktra)
               trae_km1_c(I,ktra) = trae(index1(i),K-1,ktra)
               dtra_k_c(I,ktra)   = dtrabydt(index1(i),K,ktra)
               dtra_km1_c(I,ktra) = dtrabydt(index1(i),K-1,ktra)
             End do
          End do

        End if


!----------------------------------------------------------------------
! If below convective cloud base downdraught note allowed for form
!----------------------------------------------------------------------

        Do i=1,ndd
         If (k <  ICCB_C(i)) bddi_C(i)=.false.
        End do

!-----------------------------------------------------------------------
! Reset EN/DETRAINMENT rates for downdraught
!-----------------------------------------------------------------------
! Possible problem with calculation of kmin in layer_DD when looping
! over levels starting well above the termination level.
! kmin only calculated in layer_DD if k=klev+1

        if ( k /= klev+1) then
          Do i=1,ndd

            If (kterm(index1(i))+1 == k) then

              kmin(i)=klev+2      ! required to ensure amdetk calculated
                                  ! correctly in layer_dd
            End if
          End do
        End if

! DEPENDS ON: layer_dd
        call LAYER_DD (ndd,k,Klev,the_k_c,the_km1_c,flx_strt_c,         &
     &               P_layer_centres_c,P_layer_boundaries_c,            &
     &               exner_layer_centres_c,                             &
     &               exner_km12_c,exner_kp12_c,                         &
     &               exner_km32_c,pstar_c,pk,p_km1,delpk,delpkm1,exk,   &
     &               exkm1,amdetk,ekm14,ekm34,kmin,bddi_C,              &
     &               recip_pstar_c)

!----------------------------------------------------------------------
! If level k within 150mb of surface then downdraught not allowed to
! form
!----------------------------------------------------------------------

        Do i=1,ndd
         If (PK(i) >  (pstar_c(i)-15000.0)) bddi_c(i)=.false.
        End do


!-----------------------------------------------------------------------
! Initialise downdraught
! Downdraught not allowed to form from cloud top layer (Klev+1)
! or from below cloud base
!-----------------------------------------------------------------------

! DEPENDS ON: dd_init
        call DD_INIT(ndd,npossdd,th_k_c,q_k_c,the_k_c,qe_k_c,pk,exk,    &
     &              thdd_k,qdd_k,deltd,delqd,dbb_start_c,k,bddi_c,      &
     &              bdd_on_c,l_tracer,ntra,tra_k_c,                     &
     &              trae_k_c,tradd_k,deltrad)

!-----------------------------------------------------------------------
! Update mask for where downdraught occurs
!-----------------------------------------------------------------------

        Do i=1,ndd
          IF (dbb_start_c(i).or.bdd_on_c(i)) bdd_on_c(i)=.true.
        End do

!
! If downdraught initiated set diagnostic array
!
        If(FLG_DWN_FLX) then
          Do i=1,ndd
           If(dbb_start_C(i)) dd_flux(index1(i),K)=flx_dd_k(index1(i))
          End do
        End if

        nddon_tmp = 0
        Do i=1,ndd
          IF (bdd_on_c(i)) then
            nddon_tmp = nddon_tmp+1
          End if
        End do


!-----------------------------------------------------------------------
! Call downdraught routine
!-----------------------------------------------------------------------
        
! DEPENDS ON: downd
        call DOWND(ndd,npossdd,K,Klev,thdd_k,qdd_k,the_k_c,the_km1_c,   &
     &           qe_k_c,qe_km1_c,dthbydt_k_c,dthbydt_km1_c,dqbydt_k_c,  &
     &           dqbydt_km1_c,flx_dd_k_c,p_km1,delpk,delpkm1,exk,       &
     &           exkM1,deltd,delqd,amdetk,ekm14,ekm34,precip_k_c,       &
     &           rain_c,snow_c,ICCB_c,bwater_k_c,dbb_start_C,           &
     &           bddwt_k_c,bddwt_km1_c,bdd_on_C,rain_env,snow_env,      &
     &           rain_dd,snow_dd,flx_ud_k_c,timestep,CCA_C,nddon_tmp,   &
     &           Lr_UD_ref,L_tracer,ntra,tradd_k,                       &
     &           trae_k_c,trae_km1_c,dtra_k_c,dtra_km1_c,deltrad)

!-----------------------------------------------------------------------
! Decompress/expand those variables which areto be output
!-----------------------------------------------------------------------

        Do i=1,ndd
          dthbydt(index1(i),K)   = dthbydt_k_c(i)
          dthbydt(index1(i),K-1) = dthbydt_km1_c(i)
          dqbydt(index1(i),K)    = dqbydt_k_c(i)
          dqbydt(index1(i),K-1)  = dqbydt_km1_c(i)
        End do
!
! Need to check that point would be selected in S.R DOWND or else
! not sensible to set entrainment and detrainment rates  in diagnostics
!
        If(FLG_DWN_FLX) then
          Do i=1,ndd
            If(bdd_on_c(i)) then
              dd_flux(index1(i),K-1) = flx_dd_k_c(i)
            End if
          End do
        End if
        If(FLG_ENTR_DWN) then
          Do i=1,ndd
            If(bdd_on_C(i)) then
              entrain_dwn(index1(i),K)=(1.0-amdetk(i))*                 &
     &                           (ekm14(i)+ekm34(i)*(1.0+ekm14(i)))*    &
     &                             dd_flux(index1(i),k)
            End if
          End do
        End if
        If(FLG_DETR_DWN) then
          Do i=1,ndd
            If(bdd_on_c(i)) then
             detrain_dwn(index1(i),K)=-amdetk(i)*dd_flux(index1(i),K)
            End if
          End do
        End if

        Do i=1,ndd
          If (k == 2) then
            rain(index1(i)) = rain_c(i)
            snow(index1(i)) = snow_c(i)
          End if
          precip(index1(i),K) = precip_k_c(i)
        End do

        Do i=1, ndd
           rain_3d(index1(i),k-1) = rain_3d(index1(i),k-1)                    &
                                  + rain_dd(i) + rain_env(i)

           snow_3d(index1(i),k-1) = snow_3d(index1(i),k-1)                    &
                                  + snow_dd(i) + snow_env(i)
        End Do

        If(L_TRACER)then

          Do ktra=1,ntra
            Do i=1,ndd
              dtrabydt(index1(i),K,ktra)   = dtra_k_c(i,ktra)
              dtrabydt(index1(i),K-1,ktra) = dtra_km1_c(i,ktra)
            End do
          End do

        End if

!----------------------------------------------------------------------
!   End of main K loop
!----------------------------------------------------------------------

      End do    ! End of main level loop


      Return
      END SUBROUTINE DD_ALL_CALL

