#if defined(CONTROL) || defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      subroutine intf_unit ( internal_model,  intf_area_no,             &
     &                       intf_unit_no )
!L  Purpose: calculates output interface area number from
!L           file number and sub model number
      implicit none
      integer internal_model   ! IN  Internal sub-model number
      integer intf_area_no     ! IN  Interface area number
      integer intf_unit_no     ! OUT Interface unit number

#include "csmid.h"
#include "cmaxsize.h"
#include "cmaxsizo.h"
#include "cintfa.h"

#if defined(ATMOS)
      if ( internal_model  ==  a_im) then
        intf_unit_no = lbc_unit_no_a(intf_area_no)
      endif
#endif

      return
      END SUBROUTINE intf_unit
#endif
! ---------------------------------------------------------------------
