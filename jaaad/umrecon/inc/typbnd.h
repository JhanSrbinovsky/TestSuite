! TYPBND - needs TYPSIZE included first
!LL
!LL  4.5  04/08/97 Add O_BDY_STEP_PREV for ocean boundary routines
!LL                Delete FLOOR_STEPSO.          C.G. Jones
!LL  5.0  28/04/99 Remove references to FLOOR variables        P.Burton
!LL  5.2  20/09/00 Removed old LBC stuff
!LL                Removed 2nd dimension (for lower boundary)
!LL                Added LOOKUP_COMP_BOUNDA             P.Burton
!LL  5.3  15/10/01 Include additional ocean boundary data    M J Bell
!LL  5.5  17/02/03 Include Wave model boundary data. D.Holmes-Bell
!LL  6.2  23/11/05  Removed all references to the wavemodel.
!LL                 T.Edwards
#include "cbound.h"

#if defined(ATMOS) && !defined(GLOBAL)
      ! Headers from atmosphere boundary data sets
      ! Second index of header arrays = 1 Lateral boundary data
      ! = 2 Lower boundary data
      INTEGER :: FIXHD_BOUNDA(LEN_FIXHD,1)      ! Fixed header
      INTEGER :: INTHD_BOUNDA(A_LEN_INTHD,1)    ! Integer header

        ! Lookups for the first set of LBCs in the LBC file
      INTEGER :: LOOKUP_BOUNDA(LEN1_LOOKUP,RIM_LOOKUPSA)

      ! Varying items from the LOOKUP table for the entire LBC file
      INTEGER :: LOOKUP_COMP_BOUNDA(LEN1_LBC_COMP_LOOKUP,BOUND_LOOKUPSA)

      REAL ::  REALHD_BOUNDA(A_LEN_REALHD,1)   ! Real header
#endif
#if defined(OCEAN) && defined(BOUNDSO)
      ! Headers from ocean boundary data sets
      ! Second index of header arrays = 1 Lateral boundary data
      ! = 2 Lower boundary data
      INTEGER :: FIXHD_BOUNDO(LEN_FIXHD,2)      ! Fixed header
      INTEGER :: INTHD_BOUNDO(O_LEN_INTHD,2)    ! Integer header
      INTEGER :: LOOKUP_BOUNDO(LEN1_LOOKUP,BOUND_LOOKUPSO)  ! Lookups
      REAL    :: REALHD_BOUNDO(O_LEN_REALHD,2)   ! Real header
      integer joc_bdy_tracer(nt)  ! ptrs to fields in bdy data
      integer o_bdy_item_codes(rim_lookupso) ! item code for each fld
      real    o_bdy_prev(lenrimdata_o)  ! old boundary data
      real    o_bdy_next(lenrimdata_o)  ! next boundary data
! Additional scalar boundary data pointers
      INTEGER :: joc_bdy_u          ! u velocity boundary data
      INTEGER :: joc_bdy_v          ! v velocity boundary data
      INTEGER :: joc_bdy_stream     ! stream function bdy data
      INTEGER :: joc_bdy_tend       ! strm ftn tendency bdy data
      INTEGER :: joc_bdy_ztd        ! ztd bdy data
      INTEGER :: joc_bdy_snow       ! snow depth bdy data
      INTEGER :: joc_bdy_aice       ! ice concentration bdy data
      INTEGER :: joc_bdy_hice       ! ice mean depth bdy data
      COMMON/CBNDO/                                                     &
     &  joc_bdy_u, joc_bdy_v, joc_bdy_stream, joc_bdy_tend,             &
     & joc_bdy_ztd, joc_bdy_snow, joc_bdy_aice, joc_bdy_hice
#endif
      !  Control data calculated from namelist
      INTEGER :: RIM_STEPSA      ! Set by IN_BOUND from BOUND_FIELDCODE
      INTEGER :: RIM_STEPSO      ! Set by IN_BOUND from BOUND_FIELDCODE
      INTEGER::RIM_STEPSW      ! Set by IN_BOUND from BOUND_FIELDCODE
      INTEGER :: NBOUND_LOOKUP(2)
      INTEGER :: O_BDY_STEP_NEXT ! timestep for which next boundary data
                          ! is valid. Calculated in INBOUND / UPBOUND
      INTEGER::W_BDY_STEP_NEXT
      COMMON/CBND/                                                      &
     & RIM_STEPSA,RIM_STEPSO,O_BDY_STEP_NEXT                            &
     &  ,W_BDY_STEP_NEXT,RIM_STEPSW
! TYPBND end
