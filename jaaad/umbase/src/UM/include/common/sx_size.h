!========================== COMDECK SX_SIZE ====================
!
!   Description:
!
!   This COMDECK contains COMMON blocks for landpt only and LBC
!   field for SX
!
! -------------------- COMMON BLOCKS --------------------
!
! =======================================================
! Common block for small executables
! =======================================================


! Variables and sizes for land-points only field
! ----------------------------------------------
      INTEGER:: LAND_FIELD           ! IN: No of land points in field
      INTEGER:: TR_VARS
      INTEGER:: TOT_LEVELS
      INTEGER:: global_ROWS          ! IN: No of global (theta) rows
      INTEGER:: global_ROW_LENGTH    ! IN: Points per global row
      INTEGER:: global_land_field    !  Global no of land points
      INTEGER:: local_land_field     !  Local no of land points
      COMMON /mpp_landpts/ global_land_field,local_land_field


! Variables and sizes for LBC
! ---------------------------

#include "rimtypes.h"

      INTEGER :: RIMWIDTHA(Nrima_max)
      INTEGER :: NRIM_TIMESA      ! IN: Max no of timelevels in rim flds

      ! Data structure sizes for OCEAN BOUNDARY file control routines
      INTEGER :: RIMWIDTHO     ! IN: No of points width in rim fields
      INTEGER :: NRIM_TIMESO   ! IN: Max no of timelevels in rim flds

      ! Data structure sizes for WAVE BOUNDARY file control routines
      INTEGER :: RIMWIDTHW    ! IN: No of points width in rim fields
      INTEGER :: NRIM_TIMESW  ! IN: Max no of timelevels in rim flds

      ! Size of atmos LBC for given field type, halo type and rimwidth
      ! type
      INTEGER:: LENRIMA(Nfld_max,NHalo_max,Nrima_max)

      ! Size of given side (PNorth,PEast,PSouth and PWest), field type,
      ! halo type and rimwidth type
      INTEGER:: LBC_SIZEA(4,Nfld_max,NHalo_max,Nrima_max)

      ! Start of a given side within the LBC
      INTEGER:: LBC_STARTA(4,Nfld_max,NHalo_max,Nrima_max)

      ! Start of a given side within the LBC on a given processor
      INTEGER:: g_LBC_STARTA(4,Nfld_max,NHalo_max,Nrima_max,0:Maxproc-1)

      ! Size of atmos LBC on disk for given field type, halo type and
      ! rimwidth type
      INTEGER:: global_LENRIMA(Nfld_max,NHalo_max,Nrima_max)

      ! Size of given side, field type and halo type
      INTEGER:: global_LBC_SIZEA(4,Nfld_max,NHalo_max,Nrima_max)

      ! Start of a given side within the LBC
      INTEGER:: global_LBC_STARTA(4,Nfld_max,NHalo_max,Nrima_max)

      ! Variables describing the Ocean Lateral Boundary Conditions
      INTEGER:: LENRIMO                ! Size of ocean LBC (theta)
      INTEGER:: LENRIMO_U              ! Size of ocean LBC (velocity)

      ! Variables that may be needed for vn5.2 but have not yet been
      ! dealt with at vn5.1
      INTEGER:: RIMFLDSA
      INTEGER:: RIMFLDSO
      INTEGER:: global_LENRIMDATA_A
      INTEGER:: global_LENRIMDATA_W
      INTEGER:: LENRIMDATA_A
      INTEGER:: LENRIMDATA_O
      INTEGER:: LENRIMDATA_W
      INTEGER:: BOUNDFLDS
      INTEGER:: RIM_LOOKUPSA
      INTEGER:: RIM_LOOKUPSO
      INTEGER:: RIM_LOOKUPSW
      INTEGER:: BOUND_LOOKUPSA
      INTEGER:: BOUND_LOOKUPSO
      INTEGER:: BOUND_LOOKUPSW

      COMMON/SX_NLSIZES/                                                &
     &  RIMWIDTHA, NRIM_TIMESA,                                         &
     &  RIMWIDTHO, NRIM_TIMESO,                                         &
     &  RIMWIDTHW, NRIM_TIMESW,                                         &
     &  LAND_FIELD, TR_VARS, TOT_LEVELS,                                &
     &  GLOBAL_ROW_LENGTH, GLOBAL_ROWS

      COMMON/DRSIZ_BO/                                                  &
      ! Atmosphere variables
     &  LENRIMA, LBC_SIZEA, LBC_STARTA, g_LBC_STARTA,                   &
     &  global_LENRIMA,global_LBC_SIZEA,global_LBC_STARTA,              &
      ! Wave model variables
     &  RIM_LOOKUPSW, LENRIMDATA_W, global_LENRIMDATA_W,                &
      ! Ocean variables
     &  LENRIMO, LENRIMO_U,                                             &
      ! Variables still to be dealt with
     &  RIMFLDSA,RIMFLDSO,BOUNDFLDS,RIM_LOOKUPSA,RIM_LOOKUPSO,          &
     &  BOUND_LOOKUPSA,BOUND_LOOKUPSO,BOUND_LOOKUPSW,                   &
     &  global_LENRIMDATA_A,                                            &
     &  LENRIMDATA_A,LENRIMDATA_O
! SX_SIZE end
