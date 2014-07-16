#if defined(FLUXPROC) || defined(FLXPLPR) || defined(FLXPLIN)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
! History:
! version  date         change
! 4.5      03/09/98     New code
! 5.3      10/10/00     Split to enable flux selection and grid
!                       interpolation to be split between 2 programs.
!                       M. J. Bell and A. Hines.
!
! Author:     M. J. Bell
!----------------------------------------------------------------------
! contains routines: amend_lookup_flux
!                    amend_lookup_grid
!                    amend_lookup
!
! Purpose: Flux processing Routine.
!          Amends grid information in Int_Head and Real_Head
!          Also amends packing and level information
!----------------------------------------------------------------------
#if !defined(FLXPLPR) && !defined(FLXPLIN)
!----------------------------------------------------------------------
#endif
#if defined(FLXPLPR) || defined(FLXPLIN)
!----------------------------------------------------------------------
#endif
#if defined(FLXPLIN)
      subroutine amend_lookup_grid (                                    &
     &        LookuplsmO, Int_Head, Real_Head, output_land_value)

      implicit none

! Purpose:  Amends the grid information in the real and integer
!           headers for a field

! declaration of parameters
#include "plookups.h"
#include "clookadd.h"

! declaration of argument list
      integer LookuplsmO(Len1_Lookup) ! IN lookup table for lsm
      integer Int_Head(Len_IntHd) ! IN/OUT int part of lookup table
      real Real_Head(Len_RealHd)  ! IN/OUT real part of lookup table
      real output_land_value      ! IN value at land points (real MDI)

! no other variables  used

! declaration of externals
      external copy_to_real
!----------------------------------------------------------------------
! 1. Set integer part of Lookup table
      Int_Head(LBLREC) = LookuplsmO(LBLREC)
      Int_Head(LBCODE) = LookuplsmO(LBCODE)
      Int_Head(LBHEM) = LookuplsmO(LBHEM)
      Int_Head(LBROW) = LookuplsmO(LBROW)
      Int_Head(LBNPT) = LookuplsmO(LBNPT)

! 2. Set the real part of Lookup table
! DEPENDS ON: copy_to_real
      call copy_to_real( LookuplsmO(BPLAT),                             &
     &                             Real_Head(BPLAT - Len_IntHd) )
! DEPENDS ON: copy_to_real
      call copy_to_real( LookuplsmO(BPLON),                             &
     &                             Real_Head(BPLON - Len_IntHd) )
! DEPENDS ON: copy_to_real
      call copy_to_real( LookuplsmO(BZY),                               &
     &                             Real_Head(BZY - Len_IntHd) )
! DEPENDS ON: copy_to_real
      call copy_to_real( LookuplsmO(BDY),                               &
     &                             Real_Head(BDY - Len_IntHd) )
! DEPENDS ON: copy_to_real
      call copy_to_real( LookuplsmO(BZX),                               &
     &                             Real_Head(BZX - Len_IntHd) )
! DEPENDS ON: copy_to_real
      call copy_to_real( LookuplsmO(BDX),                               &
     &                             Real_Head(BDX - Len_IntHd) )

      Real_Head(BMDI - Len_IntHd) = output_land_value

      return
      END SUBROUTINE amend_lookup_grid
#endif

#endif
