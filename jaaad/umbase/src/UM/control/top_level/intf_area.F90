#if defined(CONTROL) || defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      subroutine intf_area ( internal_model,  intf_unit_no,             &
     &                       intf_area_no )
!L  Purpose: calculates output interface file number from
!L           unit number and sub model number
      implicit none
      integer internal_model   ! IN  Internal sub-model number
      integer intf_unit_no     ! IN  Interface unit number
      integer intf_area_no     ! OUT Interface area number
!
      integer jintf            ! Local : Loop index
#include "csmid.h"
#include "cmaxsize.h"
#include "cmaxsizo.h"
#include "cintfa.h"

#if defined(ATMOS)
      if ( internal_model  ==  a_im) then
        do jintf=1,max_n_intf_a
          if ( intf_unit_no  ==  lbc_unit_no_a(jintf) ) then
            intf_area_no = jintf
          endif
        enddo
      endif
#endif

      return
      END SUBROUTINE intf_area
! ---------------------------------------------------------------------
#endif
! ---------------------------------------------------------------------
