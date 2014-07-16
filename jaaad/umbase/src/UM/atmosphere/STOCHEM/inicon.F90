#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE INICON(xx,ipos,nfill,month,year)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : Initialize species concentrations
!-
!-   Inputs  : XX
!-   Outputs : XX
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  4.0    09/12/93  Created.  W.J. Collins
!  5.5    04/01/02  New Dynamics verion. W.J. Collins
!  5.5    13/02/04  Minor changes to filename construction. K. Ketelsen
!  5.5    30/07/04  Now initialises methane based on year using IPCC
!                   values, scales other species to these values.
!                   C.E. Johnson
!  6.2    01/03/06  Initialisation for extra scenarios added, and
!                   a default.  M.G. Sanderson
!-
!VVV  V2.3  INICON 22/X/99 - Species magic numbers removed.
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER, DIMENSION(3,nclprc), intent(in) :: ipos
      INTEGER, INTENT(IN) :: nfill, month, year
      REAL, DIMENSION(nc,nclprc), intent(out) :: xx

      INTEGER :: i, j, k, l
      INTEGER :: ryear
      INTEGER, PARAMETER :: ndecades = 14            ! 1970-2100

      REAL :: height
      REAL :: ch4_global
      REAL :: r

! Pre-industrial CH4 (ppbv)
      REAL :: ch4_pi = 709.0

      REAL, DIMENSION(nc) :: c_init
      REAL, DIMENSION(mnlat,nlev):: zo3, zo3s
      REAL, DIMENSION(mnlat) :: ch4_ppb, co_ppb, c2h6_ppt

! "Standard" Initial concentrations. Note that they are ordered
! S -> N.
      REAL, DIMENSION(mnlat) :: co_ppb0 = (/                            &
     &    46.0,   46.0,   46.0,   46.0,   46.0,   46.0,   46.0,   46.0, &
     &    46.0,   46.0,   46.0,   46.0,   46.0,   46.0,   46.0,   46.0, &
     &    46.0,   46.0,   46.0,   46.0,   46.0,   46.0,   46.0,   46.0, &
     &    46.0,   46.0,   46.0,   46.0,   49.0,   49.0,   49.0,   49.0, &
     &    59.0,   59.0,   59.0,   59.0,   64.0,   64.0,   64.0,   64.0, &
     &    75.0,   75.0,   75.0,   75.0,   92.0,   92.0,   92.0,   92.0, &
     &   113.0,  113.0,  113.0,  113.0,  126.0,  126.0,  126.0,  126.0, &
     &   126.0,  126.0,  126.0,  126.0,  126.0,  126.0,  126.0,  126.0, &
     &   126.0,  126.0,  126.0,  126.0,  126.0,  126.0,  126.0,  126.0  &
     &   /)

      REAL, DIMENSION(mnlat) :: ch4_ppb0 = (/                           &
     &  1667.0, 1667.0, 1667.0, 1667.0, 1667.0, 1667.0, 1667.0, 1667.0, &
     &  1667.0, 1667.0, 1667.0, 1667.0, 1667.0, 1667.0, 1667.0, 1667.0, &
     &  1667.0, 1667.0, 1667.0, 1667.0, 1667.0, 1667.0, 1667.0, 1667.0, &
     &  1670.0, 1670.0, 1670.0, 1670.0, 1672.0, 1672.0, 1672.0, 1672.0, &
     &  1688.0, 1688.0, 1688.0, 1688.0, 1700.0, 1700.0, 1700.0, 1700.0, &
     &  1727.0, 1727.0, 1727.0, 1727.0, 1750.0, 1750.0, 1750.0, 1750.0, &
     &  1769.0, 1769.0, 1769.0, 1769.0, 1784.0, 1784.0, 1784.0, 1784.0, &
     &  1793.0, 1793.0, 1793.0, 1793.0, 1795.0, 1795.0, 1795.0, 1795.0, &
     &  1807.0, 1807.0, 1807.0, 1807.0, 1803.0, 1803.0, 1803.0, 1803.0  &
     &  /)

      REAL, DIMENSION(mnlat) :: c2h6_ppt0 = (/                          &
     &   360.0,  360.0,  360.0,  360.0,  360.0,  360.0,  360.0,  360.0, &
     &   360.0,  360.0,  360.0,  360.0,  360.0,  360.0,  360.0,  360.0, &
     &   360.0,  360.0,  360.0,  360.0,  360.0,  360.0,  360.0,  360.0, &
     &   380.0,  380.0,  380.0,  380.0,  400.0,  400.0,  400.0,  400.0, &
     &   550.0,  550.0,  550.0,  550.0,  600.0,  600.0,  600.0,  600.0, &
     &   700.0,  700.0,  700.0,  700.0,  900.0,  900.0,  900.0,  900.0, &
     &  1100.0, 1100.0, 1100.0, 1100.0, 1300.0, 1300.0, 1300.0, 1300.0, &
     &  1500.0, 1500.0, 1500.0, 1500.0, 1650.0, 1650.0, 1650.0, 1650.0, &
     &  1600.0, 1600.0, 1600.0, 1600.0, 1500.0, 1500.0, 1500.0, 1500.0  &
     &  /)

      INTEGER, DIMENSION(ndecades) :: ipcc_years =                      &
     &  (/ (1960+i*10, i=1,ndecades) /)

! IPCC (2001) Methane abundances (ppbv) 1970-2100

      REAL, DIMENSION(ndecades) :: ch4_ipcc_a1b = (/                    &
     &  1420.0, 1570.0, 1700.0, 1760.0, 1871.0, 2026.0, 2202.0,         &
     &  2337.0, 2400.0, 2386.0, 2301.0, 2191.0, 2078.0, 1974.0 /)

      REAL, DIMENSION(ndecades) :: ch4_ipcc_a1fi = (/                   &
     &  1420.0, 1570.0, 1700.0, 1760.0, 1851.0, 1986.0, 2175.0,         &
     &  2413.0, 2668.0, 2875.0, 3030.0, 3175.0, 3307.0, 3413.0 /)

      REAL, DIMENSION(ndecades) :: ch4_ipcc_a2 = (/                     &
     &  1420.0, 1570.0, 1700.0, 1760.0, 1851.0, 1997.0, 2163.0,         &
     &  2357.0, 2562.0, 2779.0, 3011.0, 3252.0, 3493.0, 3731.0 /)

      REAL, DIMENSION(ndecades) :: ch4_ipcc_b1 = (/                     &
     &  1420.0, 1570.0, 1700.0, 1760.0, 1827.0, 1891.0, 1927.0,         &
     &  1919.0, 1881.0, 1836.0, 1797.0, 1741.0, 1663.0, 1574.0 /)

      REAL, DIMENSION(ndecades) :: ch4_ipcc_b2 = (/                     &
     &  1420.0, 1570.0, 1700.0, 1760.0, 1839.0, 1936.0, 2058.0,         &
     &  2201.0, 2363.0, 2510.0, 2639.0, 2765.0, 2872.0, 2973.0 /)

      CHARACTER(len=80) :: cmessage    ! Error message
      CHARACTER(len=72) :: file_name
      CHARACTER(len=14) :: name
      CHARACTER(len=8)  :: species

! Reversed latitudes for New Dynamics
! Present day conditions, 72x20 crude for now
!      DATA ((O3_PPB(I,J),I=MNLAT,1,-1),J=1,NLEV) /
!     &              72*35.0, 72*35.0, 72*40.0,
!     &              72*55.0, 72*62.0, 72*62.0,
!     &              72*68.0, 72*68.0, 72*80.0,
!     &              72*80.0, 72*90.0, 72*90.0,
!     &              72*110.0,72*110.0,72*150.0,
!     &              72*150.0,72*200.0,72*200.0,
!     &              72*1000.0,72*1000.0/

      IF (l_emiss_current) THEN
        ryear = year
      ELSE
        ryear = emiss_year
      END IF
      WRITE(6,*) 'INICON: Initialising methane for year: ', ryear
      IF (scenario /= 'pi' .AND. (ryear < 1970 .OR. ryear > 2100)) THEN
        cmessage =                                                      &
     &  'Specified year out of range for IPCC scenarios'
! DEPENDS ON: ereport
        CALL EREPORT('STOCH_INICON1',1,cmessage)
      END IF

! Interpolate IPCC methane, then scale latitudinal profile.
      i = ((ryear - 1960) / 10) + 1
      SELECT CASE(Scenario)
        CASE ('pi')
          ch4_global = ch4_pi
        CASE ('ab')
          ch4_global = ch4_ipcc_a1b(i-1)+(ryear-ipcc_years(i-1)) *      &
     &      ((ch4_ipcc_a1b(i)-ch4_ipcc_a1b(i-1))/10.0)
        CASE ('fi')
          ch4_global = ch4_ipcc_a1fi(i-1)+(ryear-ipcc_years(i-1)) *     &
     &      ((ch4_ipcc_a1fi(i)-ch4_ipcc_a1fi(i-1))/10.0)
        CASE ('a2','A2')
          ch4_global = ch4_ipcc_a2(i-1)+(ryear-ipcc_years(i-1)) *       &
     &      ((ch4_ipcc_a2(i)-ch4_ipcc_a2(i-1))/10.0)
        CASE ('b1')
          ch4_global = ch4_ipcc_b1(i-1)+(ryear-ipcc_years(i-1)) *       &
     &      ((ch4_ipcc_b1(i)-ch4_ipcc_b1(i-1))/10.0)
        CASE ('b2')
          ch4_global = ch4_ipcc_b2(i-1)+(ryear-ipcc_years(i-1)) *       &
     &      ((ch4_ipcc_b2(i)-ch4_ipcc_b2(i-1))/10.0)
        CASE ('bu','mf') ! Assume same as B2 for now
          ch4_global = ch4_ipcc_b2(i-1)+(ryear-ipcc_years(i-1)) *       &
     &      ((ch4_ipcc_b2(i)-ch4_ipcc_b2(i-1))/10.0)
        CASE DEFAULT     ! Use B2
          ch4_global = ch4_ipcc_b2(i-1)+(ryear-ipcc_years(i-1)) *       &
     &      ((ch4_ipcc_b2(i)-ch4_ipcc_b2(i-1))/10.0)
      END SELECT

! Scale initial concentrations to IPCC methane values
      r = ch4_global / (SUM(ch4_ppb0)/mnlat)
      ch4_ppb = ch4_ppb0  * r
      co_ppb  = co_ppb0   * r
      c2h6_ppt= c2h6_ppt0 * r

      WRITE(6,*) 'Global Mean Initial CH4 (ppb): ',                     &
     &  SUM(ch4_ppb)/SIZE(ch4_ppb)
      WRITE(6,*) 'Global Mean Initial CO (ppb):',                       &
     &  SUM(co_ppb)/SIZE(co_ppb)
      WRITE(6,*) 'Global Mean Initial C2H6: (ppt)',                     &
     &  SUM(c2h6_ppt)/SIZE(c2h6_ppt)

! Initialise concentrations of tropospheric gases
! (volumetric mixing ratio units)

      c_init = 4.0e-18           ! Everything else
      c_init(i_od)=4.0e-23       ! O1D
      c_init(i_no)=7.5e-12       ! NO
      c_init(i_no2)=22.5e-12     ! NO2
      c_init(i_hcho)=0.2e-9      ! HCHO
      c_init(i_h2)=5.62e-07      ! H2
      c_init(i_hno3)=300.0e-12   ! HNO3
      c_init(i_c2h6)=734.0e-12   ! C2H6
      c_init(i_ch3cho)=4.0e-13   ! CH3CHO
      c_init(i_nc4h10)=80.0e-12  ! NC4H10
      c_init(i_c2h4)=59.0e-12    ! C2H4
      c_init(i_c3h6)=12.0e-12    ! C3H6
      c_init( (/i_be7,i_be10,i_rn222,i_pb210/) ) = 0.0   ! Radionuclides

! Read in zonal mean ozone field
      file_name = TRIM(datdir)//'zonal_o3_72.dat'
      WRITE(0,*) 'Open: ',file_name
      OPEN(8,FILE=file_name,STATUS='OLD')
      READ(8,*) name
      READ(8,*) species
      WRITE(6,*) 'Reading zonal mean ',species,' from ',name
      DO k=1,20
        READ(8,*) height
        DO i = 0, 64, 8
          READ(8,'(8e9.2)') zo3(i+1:i+8,k)
        END DO
      END DO
      zo3(:,:) = zo3(72:1:-1,:)
      CLOSE(8)

! Read in zonal mean ozone (S) field
      OPEN(8,file=TRIM(datdir)//'zonal_o3s_72.dat',STATUS='OLD')
      READ(8,*) name
      READ(8,*) species
      WRITE(6,*) 'Reading zonal mean ',species,' from ',name
      DO k = 1, 20
        READ(8,*) height
        DO I = 0, 64, 8
          READ(8,'(8e9.2)') zo3s(i+1:i+8,k)
        END DO
      END DO
      zo3s(:,:) = zo3s(72:1:-1,:)
      CLOSE(8)

! Copy concentrations to Lagrangian cells
      DO k=1,nc
        xx(k,:) = c_init(k)
      END DO

      DO l=1,nfill
        j=ipos(2,l)                              ! lat band
        k=ipos(3,l)                              ! height index
        xx(i_o3,l)      = zo3(j,k)
        xx(i_o3_strat,l)= zo3s(j,k)
        xx(i_co,l)      = co_ppb(j)*1.0e-9       ! CO
        xx(i_ch4,l)     = ch4_ppb(j)*1.0e-9      ! CH4
        xx(i_c2h6,l)    = c2h6_ppt(j)*1.0e-12    ! C2H6
      END DO

      END SUBROUTINE INICON
#endif
