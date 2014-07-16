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
!
! Author:     S. A. Spall
! 5.3      14/08/01     Allow atmosphere fields on B or C grids MJ
!----------------------------------------------------------------------
! contains routines: read_select_cntl
!
! Purpose: Flux processing routine.
!          Reads files controlling which fluxes to process
!----------------------------------------------------------------------
      subroutine read_select_cntl ( icode )

      implicit none

! declaration of argument list
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! declaration of parameters

! declaration of globals used
#include "cunitnos.h"
#include "cmess.h"
#include "cselct.h"

! No local arrays

! No local scalars

! namelist declaration
      NAMELIST /NmLstSlt/                                               &
     &           l_B_grid,                                              &
     &           l_winds_slt, l_heat_slt, l_moisture_slt,               &
     &           l_sea_ice_slt, l_references_slt, l_pressure_slt,       &
     &           l_windspd_slt

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_select_cntl'  ! subroutine name for error messages

! 1. set default values for variables in NmLstSlt

      l_B_grid = .false.

      l_winds_slt = .false.
      l_heat_slt = .false.
      l_moisture_slt = .false.
      l_sea_ice_slt = .false.
      l_references_slt = .false.
      l_pressure_slt = .false.
      l_windspd_slt = .false.

! 2. read debug control namelist
      read (UnitSlt, NmLstSlt, iostat = icode)

      if ( icode  /=  0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &     ' step 2. unable to read select control namelist'
        icode = 4008
        go to 9999
      end if

9999  continue
      return
      END SUBROUTINE read_select_cntl
!----------------------------------------------------------------------
#endif
