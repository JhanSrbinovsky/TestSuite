#if defined(A32_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculates LBC Grid Sizes
!
! Subroutine Interface:

      Subroutine LBC_Grid_Sizes ( jintf )

      Implicit NONE

!
! Description:
!   Computes the grid sizes for one level of LBC data for each grid
!
! Method:
!   <Say how it does it: include references to external documentation>
!   <If this routine is very complex, then include a "pseudo code"
!    description of it to make its structure and method clear>
!
! Current Code Owner: Dave Robinson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2    13/11/00   Original code. Dave Robinson
!   5.3    22/10/01   Remove lbc_rim_size + redundant code. D.Robinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!

#include "cmaxsize.h"
#include "cintfa.h"
#include "parvars.h"
#include "typsize.h"
#include "parlbcs.h"

      Integer  ::  jintf           !  Index to interface area

! Local variables

      Integer :: lbc_row_len    ! Row length including haloes
      Integer :: lbc_rows       ! No of rows including haloes
      Integer :: Ifld, Ihalo, Iside   !  Loop indices

!*---------------------------------------------------------------------
      integer intf_halosize(2,Nhalo_max)

      Integer   ::  ErrorStatus
      Character (Len=80)           :: CMessage
      Character (Len=*), Parameter :: RoutineName= 'LBC_Grid_Sizes'
!      -------------------------

      ErrorStatus = 0
      CMESSAGE=' '

!     Set up intf_halosize for this area
      intf_halosize(1,1)=1
      intf_halosize(2,1)=1

      intf_halosize(1,2)=intf_exthalo_ew(jintf)
      intf_halosize(2,2)=intf_exthalo_ns(jintf)

      intf_halosize(1,3)=0
      intf_halosize(2,3)=0

! --------------------------------------------------
! Determine global length of lbc field for one level
! --------------------------------------------------

      do ifld=1,nfld_max
        do ihalo=1,nhalo_max

          lbc_global_lenrima(ifld,ihalo) = 0

          do iside=1,4

            if ( (iside == PNorth) .or. (iside == PSouth) ) then

!     Full row length for p and v grid
!     One point less for u grid

              if (ifld == fld_type_u) then
                lbc_row_len = intf_row_length(jintf)-1
              else
                lbc_row_len = intf_row_length(jintf)
              endif
              lbc_row_len = lbc_row_len + 2 * intf_halosize(1,ihalo)

!     No of rows is rimwidth plus halo size

              lbc_rows   = intfwidtha(jintf) + intf_halosize(2,ihalo)

            else ! East or West Boundaries

!     Row length is rimwidth plus halo size

              lbc_row_len = intfwidtha(jintf) + intf_halosize(1,ihalo)

!     All rows used for p and u grid
!     One row less for v grid

              if (ifld == fld_type_v) then
                lbc_rows = intf_p_rows(jintf) - 1
              else
                lbc_rows = intf_p_rows (jintf)
              endif
              lbc_rows   = lbc_rows - 2 * intfwidtha(jintf)

            endif

!     Accumulate grid size

            if (intfwidtha(jintf) > 0 ) then
              lbc_global_lenrima(ifld,ihalo) =                          &
     &        lbc_global_lenrima(ifld,ihalo) + lbc_row_len * lbc_rows
            endif

          enddo  ! iside
        enddo    ! ihalo
      enddo      ! ifld

!  ------------------------------------------------------------
!  Calculate the sizes of the fields for interpolation
!  u and v fields differ from the actual sizes.
!  An extra P col or row is required to interpolate u/v to
!  the p grid before it is interpolated to the u or v grid.
!  ------------------------------------------------------------

      do ifld=1,nfld_max
        do ihalo=1,nhalo_max

          lbc_interp_lenrima(ifld,ihalo) = 0

          do iside=1,4

            if ( (iside == PNorth) .or. (iside == PSouth) ) then

              lbc_row_len = intf_row_length(jintf)
              lbc_row_len = lbc_row_len + 2 * intf_halosize(1,ihalo)

              lbc_rows   = intfwidtha(jintf) + intf_halosize(2,ihalo)
              if (ifld == fld_type_u .or. ifld == fld_type_v) then
                lbc_rows = lbc_rows + 1
              endif

            else ! East or West Boundaries

              lbc_row_len = intfwidtha(jintf) + intf_halosize(1,ihalo)
              if (ifld == fld_type_u .or. ifld == fld_type_v) then
                lbc_row_len = lbc_row_len + 1
              endif

              lbc_rows = intf_p_rows(jintf)
              lbc_rows = lbc_rows - 2 * intfwidtha(jintf)

            endif

            if (intfwidtha(jintf) > 0 ) then
              lbc_interp_lenrima(ifld,ihalo) =                          &
     &        lbc_interp_lenrima(ifld,ihalo) + lbc_row_len * lbc_rows
            endif

          enddo  ! iside
        enddo    ! ihalo
      enddo      ! ifld

      return
      END SUBROUTINE LBC_Grid_Sizes
#endif
