#if defined(A32_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Gather information from STASHmaster records for LBC variables.
!
! Subroutine Interface:

      SUBROUTINE LBC_SETUP (                                            &
#include "argsts.h"
#include "argppx.h"
     &  lbc_item_code,                                                  &
     &  lbc_fld_type,                                                   &
     &  lbc_halo_type,                                                  &
     &  lbc_rim_type,                                                   &
     &  lbc_level_type,                                                 &
     &  lbc_first_level,                                                &
     &  lbc_last_level,                                                 &
     &  lbc_levels,                                                     &
     &  n_lbc_vars,                                                     &
     &  jintf,                                                          &
     &  im_ident                                                        &
     & )

      IMPLICIT NONE
!
! Description:
!    For all LBC variables, gather information from STASHmaster records
!    and derive field type and number of levels.
!
! Method:
!    Extract relevant data from STASHmaster record.
!
! Current Code Owner: Dave Robinson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2    13/11/00   Original code. Dave Robinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):

#include "csubmodl.h"
#include "cmaxsize.h"
#include "cintfa.h"
#include "parvars.h"
#include "typsize.h"
#include "typsts.h"
#include "ppxlook.h"

! Subroutine arguments

      Integer  :: n_lbc_vars   !  No of LBC variables
      Integer  :: jintf        !  Index to interface area
      Integer  :: im_ident     !  Internal model identifier

      Integer  :: lbc_item_code  (n_lbc_vars)
      Integer  :: lbc_halo_type  (n_lbc_vars)
      Integer  :: lbc_fld_type   (n_lbc_vars)
      Integer  :: lbc_rim_type   (n_lbc_vars)
      Integer  :: lbc_level_type (n_lbc_vars)
      Integer  :: lbc_levels     (n_lbc_vars)
      Integer  :: lbc_first_level(n_lbc_vars)
      Integer  :: lbc_last_level (n_lbc_vars)

! Local parameters:

      Integer,           Parameter :: LBC_Sect   = 32
      Character (Len=*), Parameter :: RoutineName= 'LBC_SetUp'

! Local scalars:

      Integer  :: field_type        ! Extracted field type
      Integer  :: halo_type         ! Extracted halo type
      Integer  :: level_type        ! Extracted level type
      Integer  :: first_level       ! Extracted first level code
      Integer  :: last_level        ! Extracted last level code
      Integer  :: grid_type         ! Grid type derived from field code
      Integer  :: lbc_bottom_level  ! Bottom level derived from code
      Integer  :: lbc_top_level     ! Top level derived from code
      Integer  :: var               ! Loop index for lbc variables
      Integer  :: lbc_item          ! Item Code
      Integer  :: ErrorStatus       ! Error Code

      Character(Len=80)  :: CMessage

! Function & Subroutine calls:

      Integer get_fld_type
      Integer Exppxi

!- End of header

      ErrorStatus = 0
      CMessage    = ' '

! ---------------------------------
! Set up rim type for each variable
! ---------------------------------

!     Needs updating to use new grid type ??
      do var=1,n_lbc_vars
        if (var == 1) then
          lbc_rim_type(var) = rima_type_orog
        else
          lbc_rim_type(var) = rima_type_norm
        endif
      enddo

      Do Var = 1, n_lbc_vars

        lbc_item = lbc_item_code(var) - 32000

! -------------------------------------------------------
! Extract required information from stashmaster record
! -------------------------------------------------------

! Level type
! Top Level
! Bottom Level
! Grid Type
! Halo Type

! May be needed later
!      intf_pack_code  = exppxi
!    & (im_ident, lbc_sect, lbc_item_intfa(var), ppx_dump_packing,
!include <argppx/argppx.h>
!    &                  errorstatus, cmessage)

        level_type = exppxi                                             &
     &  (im_ident, lbc_sect, lbc_item, ppx_lv_code,                     &
#include "argppx.h"
     &   errorstatus, cmessage)

        first_level = exppxi                                            &
     &  (im_ident, lbc_sect, lbc_item, ppx_lb_code,                     &
#include "argppx.h"
     &   errorstatus, cmessage)

        last_level = exppxi                                             &
     &  (im_ident, lbc_sect, lbc_item, ppx_lt_code,                     &
#include "argppx.h"
     &   errorstatus, cmessage)

        grid_type  = exppxi                                             &
     &  (im_ident, lbc_sect, lbc_item, ppx_grid_type,                   &
#include "argppx.h"
     &   errorstatus, cmessage)

        halo_type  = exppxi                                             &
     &  (im_ident, lbc_sect, lbc_item, ppx_halo_type,                   &
#include "argppx.h"
     &   errorstatus, cmessage)

! ------------------------------
! Determine first level for LBCs
! ------------------------------

        if (first_level == 1) then        ! First atmos level
          lbc_bottom_level = 1
        elseif (first_level == 38) then   ! Surface level
          lbc_bottom_level = 0
        elseif (first_level == -1) then   ! Single level
          lbc_bottom_level = 1
        else
          write (cmessage,*) ' LBC First Level ',                       &
     &    first_level,' not recognised.'
          errorstatus = 10
! DEPENDS ON: ereport
          call ereport ( RoutineName ,ErrorStatus, CMessage)
        endif

! -----------------------------
! Determine last level for LBCs
! -----------------------------

        if (last_level == 2) then          ! Model levels
          lbc_top_level = intf_p_levels(jintf)
        elseif (last_level == 3) then      ! Wet levels
          lbc_top_level = intf_q_levels(jintf)
        elseif (last_level == 19) then     ! Model levels + 1
          lbc_top_level = intf_p_levels(jintf) + 1
        elseif (last_level == -1) then     ! Single level
          lbc_top_level = 1
        else
          write (cmessage,*) ' LBC Last Level ',                        &
     &    last_level,' not recognised.'
          errorstatus = 20
! DEPENDS ON: ereport
          call ereport ( RoutineName ,ErrorStatus, CMessage)
        endif

! --------------------------------
! Derive field type from grid type
! --------------------------------

! DEPENDS ON: get_fld_type
        field_type = get_fld_type (grid_type)

! --------------------------
! Keep extracted information
! --------------------------

        lbc_levels(var)     = lbc_top_level - lbc_bottom_level + 1
        lbc_fld_type(var)   = field_type
        lbc_halo_type(var)  = halo_type
        lbc_level_type(var) = level_type
        lbc_first_level(var)= lbc_bottom_level
        lbc_last_level(var) = lbc_top_level

      enddo   !  Var

      return
      END SUBROUTINE LBC_SETUP
#endif
