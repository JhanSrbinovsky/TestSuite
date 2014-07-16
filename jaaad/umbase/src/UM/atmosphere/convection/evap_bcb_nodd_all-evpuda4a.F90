#if defined(A05_4A) || defined(A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      SUBROUTINE EVAP_BCB_NODD_ALL (npnts,n_nodd,klev,kterm             &
     &,                      iccb, index1, bwater                       &
     &,                      exner_layer_centres,exner_layer_boundaries &
     &,                      p_layer_centres, p_layer_boundaries,pstar  &
     &,                      timestep , cca, the, qe, precip            &
     &,                      dthbydt, dqbydt, rain, snow                &
     &,                      rain_3d, snow_3d)
!
!     Purpose: To calculate the convective precipitation reaching the
!              surface.
!
!     Method : the evaporation below cloud base follows that done in
!              the downdraught routine for the environmental part of the
!              column. the points which are gathered here are those
!              points which have an updraught, but no downdraught.
!
!   Called by DEEP_CONV & SHALLOW_CONV.
!
!   Current owners of code: R A Stratton
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
      integer, intent(in) :: npnts     ! Vector length for some arrays

      integer, intent(in) :: n_nodd    ! Compressed vector length for
                                       ! calculation.
      integer, intent(in) :: klev      ! Number of levels (may be model
                                       ! levels-1 or a reduced set
                                       ! required here).

      integer, intent(in) :: kterm(npnts)
                                       ! Convective cloud top layer

      integer, intent(in) :: iccb(npnts) ! Convective cloud base
                                         ! level (m)
      integer, intent(in) :: index1(npnts) ! index of points where
                                         ! downdraught not possible
      logical, intent(in) :: bwater(npnts,2:klev+1)
                                      ! Mask for points at which
                                      ! condensate is liquid

      real , intent(in) :: exner_layer_centres(npnts,0:klev+1) !exner

      real , intent(in) :: exner_layer_boundaries(npnts,0:klev+1)
                                      ! exner at half level above
                                      ! exner_layer_centres

      real , intent(in) :: p_layer_centres(npnts,0:klev+1)
                                      ! Pressure (Pa)

      real , intent(in) :: p_layer_boundaries(npnts,0:klev+1)
                                      ! Pressure at half level above
                                      ! p_layer_centres (Pa)

      real , intent(in) :: pstar(npnts) ! Surface pressure (Pa)

      real , intent(in) :: timestep    ! timestep
      real , intent(in) :: CCA(npnts)  ! 2d convective cloud amount

      real , intent(in) :: the(npnts,klev+1)
                                 ! Model enviromental potential
                                 ! temperature (K)
      real , intent(in) :: qe(npnts,klev+1)
                                 ! Model enviromental mixing ratio
                                 ! (KG/KG)

!
! Arguments with intent INOUT:
!
      real , intent(inout) :: precip(npnts,klev+1)
                                 ! precipitation added when descending
                                 ! from layer k to k-1 (KG/M**2/S)
                                 ! (KG/KG)
      real , intent(inout) :: dthbydt(npnts,klev+1)
                                 ! increment to model potential
                                 ! temperature (K/S)
      real , intent(inout) :: dqbydt(npnts,klev+1)
                                 ! increment to model mixing ratio
                                 ! (KG/KG/S)

      real , intent(inout) :: rain(npnts)
                                 ! rainfall at surface (KG/M**2/S)

      real , intent(inout) :: snow(npnts)
                                 ! snowfall at surface (KG/M**2/S)

      real , intent(inout) :: rain_3d(npnts,klev+1)
                                 ! rainfall profile (KG/M**2/S)

      real , intent(inout) :: snow_3d(npnts,klev+1)
                                 ! snowfall profile (KG/M**2/S)
!
! Arguments with intent OUT: None
!

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------


      integer :: i,k               ! Loop counters

      integer :: index2(n_nodd)    ! compress index for each level
      integer :: nsofar            ! points in level

      integer :: ICCB_C(n_nodd)    ! Compressed cloud base level

      real :: exner_km12_c(n_nodd) ! Compressed exner function at
                                   ! layer K
      real :: exner_kp12_c(n_nodd) ! Compressed exner function at
                                   ! layer K+1
      real :: exner_km32_c(n_nodd) ! Compressed exner function at
                                   ! layer K-1
      real :: pkm1_c(n_nodd)       ! pressure of layer k-1 (Pa)

      real :: exk_c(n_nodd)        ! exner ratio for layer k
      real :: exkm1_c(n_nodd)      ! exner ratio for layer k-1

      real :: delpkm1_c(n_nodd)    ! Pressure difference across layer
                                   ! K-1 (Pa)

      real :: precip_K_C(n_nodd)   ! Compressed precipitation added when
                                   ! descending from layer k to k-1
                                   ! (kg/m**2/s)

      real :: qe_k_c(n_nodd)       ! Compressed parcel mixing ratio of
                                   ! layer K (KG/KG)
      real :: qe_km1_c(n_nodd)     ! Compressed parcel mixing ratio of
                                   ! layer k (KG/KG)

      real :: the_k_c(n_nodd)      ! Compressed parcel potential
                                   ! temperature of layer k (K)
      real :: the_km1_c(n_nodd)    ! Compressed parcel potential
                                   ! temperature of layer k-1 (K)

      real :: pstar_c(n_nodd)      ! Compressed surface pressure (Pa)

      real :: dthbydt_km1_c(n_nodd) ! Compressed increment to model
                                   ! potential temperature of layer k-1
                                   ! (K/S)

      real :: dqbydt_km1_c(n_nodd) ! Compressed increment to model
                                   ! mixing ratio of layer k-1 (KG/KG/S)

      real :: rain_c(n_nodd)       ! Amount of rainfall passing through
                                   ! environment (kg/m**2/s)

      real :: snow_c(n_nodd)       ! Amount of snowfall passing through
                                   ! environment (kg/m**2/s)

      real :: CCA_c(n_nodd)        ! Compressed convective cloud amount
                                                                       !


#include "c_r_cp.h"

!
! External routines called:
!

      External                                                          &
     &   CHG_PHSE, PEVP_BCB
!----------------------------------------------------------------------
!
! Loop over levels working downwards from maximum termination level + 1
!
      Do k = klev+1,2,-1

! How many points terminated at or above this layer ?
! Only work on these points. For the bottom levels nsofar should
! equal n_nodd

        nsofar = 0
        Do i=1,n_nodd
          If (kterm(index1(i))+1 >= k) then
            nsofar = nsofar + 1
            index2(nsofar) = index1(i)
          End if
        End do

! Compress to required points

        Do i = 1, nsofar
          the_k_c(i)   = the(index2(i),k)
          the_km1_c(i) = the(index2(i),k-1)
          qe_k_c(i)    = qe(index2(i),k)
          qe_km1_c(i)  = qe(index2(i),k-1)
          dthbydt_km1_c(i)  = dthbydt(index2(i),k-1)
          dqbydt_km1_c(i)   = dqbydt(index2(i),k-1)
          exner_km12_c(i)   = exner_layer_boundaries(index2(i),k-1)
          exner_kp12_c(i)   = exner_layer_boundaries(index2(i),k)
          exner_km32_c(i)   = exner_layer_boundaries(index2(i),k-2)
          precip_k_c(i)   = precip(index2(i),k)
          pstar_c(i) = pstar(index2(i))
          rain_c(i)  = rain(index2(i))
          snow_c(i)  = snow(index2(i))
          iccb_c(i)  = iccb(index2(i))
          cca_c(i)   = cca(index2(i))
          pkm1_c(i)    = p_layer_centres(index2(i),k-1)
          delpkm1_c(i) = p_layer_boundaries(index2(i),k-2) -            &
     &                          p_layer_boundaries(index2(i),k-1)
          exk_c(i)   = exner_layer_centres(index2(i),k)
          exkm1_c(i) = exner_layer_centres(index2(i),k-1)
        End do

        Do i = 1,nsofar
          If (bwater(index2(i),k)) then
            rain_c(i) = rain_c(i) + precip(index2(i),k)
          Else
            snow_c(i) = snow_c(i) + precip(index2(i),k)
          End if
        End do

!----------------------------------------------------------------------
! Carry out change of phase calculation for precipitation falling
! through environment
!----------------------------------------------------------------------

! DEPENDS ON: chg_phse
        call CHG_PHSE (nsofar,k,rain_c,snow_c,dthbydt_km1_c,            &
     &                 exk_c,exkm1_c,delpkm1_c,the_k_c,the_km1_c,       &
     &                 timestep,CCA_c)

!----------------------------------------------------------------------
! Reset precipitation falling through environment if downdraught
! terminates
!----------------------------------------------------------------------

! DEPENDS ON: pevp_bcb
        call PEVP_BCB (nsofar,k-1,ICCB_c,the_km1_c,pkm1_c,qe_km1_c,     &
     &                 delpkm1_c,rain_c,snow_c,dthbydt_km1_c,           &
     &                 dqbydt_km1_c,exkm1_c,timestep,CCA_c)
!
! Expand output values back into full arrays
!
        Do i=1,nsofar
          dthbydt(index2(i),k-1) = dthbydt_km1_c(i)
          dqbydt(index2(i),k-1)  = dqbydt_km1_c(i)

! Zero precipitation, as is (slyly) done in downd3c

          precip(index2(i),k) = 0.0
          rain(index2(i)) = rain_c(i)
          snow(index2(i)) = snow_c(i)
        End do

! Capture 3d rain/snow profiles
        Do i=1,nsofar
          rain_3d(index2(i),k-1) = rain_3d(index2(i),k-1) + rain_c(i)
          snow_3d(index2(i),k-1) = snow_3d(index2(i),k-1) + snow_c(i)
        End do

      End do      ! main loop over levels
!----------------------------------------------------------------------
!
      Return
      END SUBROUTINE EVAP_BCB_NODD_ALL
#endif
