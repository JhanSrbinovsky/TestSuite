#if defined(A32_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Gather information from STASHmaster records for LBC source variables
!
! Subroutine Interface:

      SUBROUTINE LBC_SRC_SETUP (                                        &
#include "argsts.h"
#include "argppx.h"
     &  item_prog,                                                      &
     &  lbc_src_fld_type,                                               &
     &  lbc_src_halo_type,                                              &
     &  lbc_src_level_type,                                             &
     &  lbc_src_first_level,                                            &
     &  lbc_src_last_level,                                             &
     &  lbc_src_levels,                                                 &
     &  n_lbc_vars,                                                     &
     &  im_ident                                                        &
     & )

      IMPLICIT NONE
!
! Description:
!   Gather information from STASHmaster records for all model
!   prognostics from which LBCs will be generated. Also derive no
!   of model levels for the prognostics.
!
! Method:
!   Extract relevant data from STASHmaster records for prognostics.
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
! Global variables

#include "csubmodl.h"
#include "cmaxsize.h"
#include "parvars.h"
#include "typsize.h"
#include "typsts.h"
#include "ppxlook.h"

! Subroutine arguments

      Integer :: n_lbc_vars     !  No of LBC variables
      Integer :: item_prog(n_lbc_vars)
      Integer :: lbc_src_fld_type(n_lbc_vars)
      Integer :: lbc_src_halo_type(n_lbc_vars)
      Integer :: lbc_src_level_type(n_lbc_vars)
      Integer :: lbc_src_first_level(n_lbc_vars)
      Integer :: lbc_src_last_level(n_lbc_vars)
      Integer :: lbc_src_levels(n_lbc_vars)

! Local parameters:

      Integer,           Parameter :: Sect0   = 0
      Character (Len=*), Parameter :: RoutineName= 'LBC_Src_SetUp'

! Local scalars:

      Integer :: grid_type        ! Grid type from STASHmaster
      Integer :: field_type       ! Field type decoded from grid type
      Integer :: halo_type        ! Halo type from STASHmaster
      Integer :: level_type       ! Level type from STASHmaster
      Integer :: first_level      ! First level code from STASHmaster
      Integer :: last_level       ! Last  level code from STASHmaster
      Integer :: src_bottom_level ! Bottom level decoded
      Integer :: src_top_level    ! Top level decoded
      Integer :: src_levels       ! No of levels for variable

      Integer :: var              ! Loop index
      Integer :: im_ident         ! Internal model identifier

      Integer            ::  ErrorStatus
      Character (Len=80) ::  CMessage ! Error message

! Function & Subroutine calls:
      Integer Get_fld_type
      Integer Exppxi
      External LEVCOD

!- End of header
! -----------------------------------------------------------

      ErrorStatus=0
      CMessage = ' '

      Do Var = 1, n_lbc_vars

! ----------------------------------------------------
! Extract required information from stashmaster record
! ----------------------------------------------------

! Level type
! Top level
! Bottom level
! Grid type
! Halo Type

        level_type = exppxi                                             &
     &  (im_ident, sect0, item_prog(var), ppx_lv_code,                  &
#include "argppx.h"
     &   ErrorStatus, cmessage)

        first_level = exppxi                                            &
     &  (im_ident, sect0, item_prog(var), ppx_lb_code,                  &
#include "argppx.h"
     &   ErrorStatus, cmessage)

        last_level = exppxi                                             &
     &  (im_ident, sect0, item_prog(var), ppx_lt_code,                  &
#include "argppx.h"
     &   ErrorStatus, cmessage)

        grid_type  = exppxi                                             &
     &  (im_ident, sect0, item_prog(var), ppx_grid_type,                &
#include "argppx.h"
     &   ErrorStatus, cmessage)

        halo_type  = exppxi                                             &
     &  (im_ident, sect0, item_prog(var), ppx_halo_type,                &
#include "argppx.h"
     &   ErrorStatus, cmessage)

! ----------------------------------
! Derive first and last source level
! ----------------------------------

        If (level_type /= 5) Then
! DEPENDS ON: levcod
          call LEVCOD(first_level,src_bottom_level,ErrorStatus,CMESSAGE)
! DEPENDS ON: levcod
          call LEVCOD(last_level, src_top_level,   ErrorStatus,CMESSAGE)
        Else
          src_bottom_level = 1
          src_top_level    = 1
        End If

! --------------------------------
! Derive field type from grid type
! --------------------------------

! DEPENDS ON: get_fld_type
        field_type = get_fld_type(grid_type)

!----------------------------
! Store extracted information
!----------------------------

        lbc_src_levels(var)     = src_top_level - src_bottom_level + 1
        lbc_src_fld_type(var)   = field_type
        lbc_src_halo_type(var)  = halo_type
        lbc_src_level_type(var) = level_type
        lbc_src_first_level(var)= src_bottom_level
        lbc_src_last_level(var) = src_top_level

      enddo

      return
      END SUBROUTINE LBC_SRC_SETUP
#endif
