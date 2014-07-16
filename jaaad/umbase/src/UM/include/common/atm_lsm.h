!====================== COMDECK ATM_LSM ========================
! Description:
!   This comdeck contains a COMMON block which contains the
!   atmosphere land sea mask - both the full field, and the
!   local subdomain on this processor.
!   This data is required for various compression/decompression
!   algorithms.
!
!   Requires AMAXSIZE comdeck to be called first for Max2DFieldSize
!
! History:
!   Model    Date     Modification history
!   version
!   4.5      12/01/98 New comdeck created.                P.Burton
!   5.2      15/01/01 Only define common block in atmos model. R. Hill
!   5.4      06/09/02 SX now uses this deck - add def UTILIO.
!                                                    E.Leung
!   5.5      17/02/03 Upgrade Wave model from 4.1 to 5.5 D.Holmes-Bell
!   6.2      23/11/05 Removed all references to the wavemodel.
!                     T.Edwards
!

      LOGICAL                                                           &
!  Full-grid land-sea mask:
     &  atmos_landmask(Max2DFieldSize)                                  &
! Local subdomain area land-sea mask:
     &, atmos_landmask_local(Max2DFieldSize)

      INTEGER atmos_number_of_landpts ! total number of land points

#if defined (ATMOS) || defined(UTILIO)
      COMMON /Atmos_LSM_Common/                                         &
     &  atmos_landmask                                                  &
     &, atmos_landmask_local                                            &
     &, atmos_number_of_landpts

#if defined(CRAY)
!DIR$ CACHE_ALIGN /Atmos_LSM_Common/
#endif

! End of comdeck ATM_LSM
#endif
