#if defined(FLUXPROC) || defined(FLXPLPR)
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
! 5.3      10/10/00     Changes to enable flux selection and grid
!                       interpolation to be split between 2 programs.
!                       M. J. Bell and A. Hines.
! 6.1      29/07/04     Correction to setting of LCyclic. A. Hines.
!
! Author:     M. J. Bell
!----------------------------------------------------------------------
! contains routines: read_field_headers, read_field_sizes
!
! Purpose: Flux processing routine.
!          Reads fixed headers and lookup tables of  flux files input
!          to FOAM_Flux_Process (i.e. Preferred or Previous fluxes and
!          climate fluxes)
!          Also works out dimensions ncols, nrowsu, nrowst
!----------------------------------------------------------------------
!----------------------------------------------------------------------
#if defined(FLXPLPR)
      subroutine read_field_sizes(                                      &
     &      Len1_Lookup, Len2_Lookup, Lookup,                           &
     &      ncols, nrowst, nrowsu, icode)

! Purpose: Flux processing routine.
!          Sets or checks dimensions of atmosphere fields

      implicit none

! declaration of INPUT arguments
      integer Len1_Lookup  ! 1st dimension of lookup table
      integer Len2_Lookup  ! 2nd dimension of lookup table
      integer Lookup(Len1_Lookup, Len2_Lookup)

! declaration of IN/OUT arguments
      integer ncols  ! number of columns on grid
      integer nrowst ! number of rows on tracer grid
      integer nrowsu ! number of rows on velocity grid
      integer icode  ! error indicator (>0 => fatal error)

! declaration of parameters used
#include "clookadd.h"
#include "cfdcodes.h"

! Globals used (as input only)
#include "cmess.h"
! local scalars
      integer ilook
      integer IStC
!-----------------------------------------------------------------
! 0. Preliminaries
!------------------------------------------------------------------

      CSub = 'read_field_sizes'

!------------------------------------------------------------------
! 1. Set ncols or check for inconsistency
!------------------------------------------------------------------

      do ilook = 1, Len2_Lookup

        if (ncols  ==  0) then
          ncols =  Lookup(LBNPT,ilook)
        else
          if (Lookup(LBNPT,ilook)  /=  ncols) then
            icode=0
            go to 9999
          end if
        end if

!------------------------------------------------------------------
! 2. decide on field type and set or check nrowst or nrowsu
!------------------------------------------------------------------

        IStC=Lookup(ITEM_CODE, ilook)
        if ( IStC  ==  StCWindSpeedU .or. IStC  ==  StCWindStressU      &
     &  .or. IStC  ==  StCWindStressV .or. IStC  ==  StCWindSpeedV)     &
     &  then
          if (nrowsu  ==  0) then
            nrowsu =  Lookup(LBROW,ilook)
          else
            if ( abs(Lookup(LBROW,ilook) - nrowsu)  >   1 ) then
              write(UnErr,*)CErr,CSub,                                  &
     &        ' Step 2.1: number of rows in field ', ilook,             &
     &        ' with stash code', IStC,                                 &
     &        ' do not match those in previous velocity fields '
              icode=1
              go to 9999
            end if
          end if
        else
          if (nrowst  ==  0) then
            nrowst =  Lookup(LBROW,ilook)
          else
            if (Lookup(LBROW,ilook)  /=  nrowst) then
              write(UnErr,*)CErr,CSub,                                  &
     &        ' Step 2.2: number of rows in field ', ilook,             &
     &        ' with stash code', IStC,                                 &
     &        ' do not match those in previous tracer fields '
              icode=1
              go to 9999
            end if
          end if
        end if

      end do !  ilook loop ends

9999  continue
      return
      END SUBROUTINE read_field_sizes
#endif
!----------------------------------------------------------------------
#endif
