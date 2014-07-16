#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine UPDATE_VEG
!
!   Model           Modification history from model version 5.5:
!  version  Date
!   5.5   20/01/03  New Deck.  M. Best.
!
!  Programming standard : UM documentation paper no3,
!                         version no.1, dated 15/01/90
!
!  Purpose: The routine is entered when any of the ancillary
!           fields have to be updated. It then checks to see if
!           leaf area index and/or canopy height have been updated.
!           If this is the case, then the subroutine SPARM is called
!           to ensure that all other vegetation parameters are
!           consistent.

      Subroutine UPDATE_VEG (                                           &
#include "argd1.h"
     &                   submodel                                       &
     &                       )

#include "parparm.h"
#include "typsize.h"
#include "typd1.h"
#include "typptra.h"
#include "cntlatm.h"
#include "nstypes.h"
#include "cancila.h"

      Integer            :: submodel
!                        !  dummy variable to close subroutine
!                        !  argument list

      Integer, Parameter :: Anc_Ref_No_Leaf_Index = 84
      Integer, Parameter :: Anc_Ref_No_Canopy_Ht  = 85


!     Local variables
      Integer            :: TILE_PTS(NTYPE)
!                        !  Number of land points which
!                        !  include the nth surface type
      Integer            :: TILE_INDEX(LAND_FIELD,NTYPE)
!                        !  Indices of land points which
!                        !  include the nth surface type


! Update vegetation parameters if required

      If (       Update (Anc_Ref_No_Leaf_Index)                         &
     &      .or. Update (Anc_Ref_No_Canopy_Ht)                          &
     &    ) Then


!-----------------------------------------------------------------------
! Call TILEPTS to initialise TILE_PTS and TILE_INDEX
!-----------------------------------------------------------------------

! DEPENDS ON: tilepts
        CALL TILEPTS(LAND_FIELD,D1(JFRAC_TYP),TILE_PTS,TILE_INDEX)

!-----------------------------------------------------------------------
! Initialise tiled and gridbox mean vegetation parameters
!-----------------------------------------------------------------------

! DEPENDS ON: sparm
        CALL SPARM (LAND_FIELD,NTILES,CAN_MODEL,TILE_PTS,TILE_INDEX,    &
     &              D1(JFRAC_TYP),D1(JCANHT_PFT),                       &
     &              D1(JLAI_PFT),D1(JSAT_SOIL_COND),D1(JCATCH_SNOW),    &
     &              D1(JCATCH_TILE),D1(JINFIL_TILE),D1(JZ0_TILE))

      End If

      Return
      END SUBROUTINE UPDATE_VEG
#endif
