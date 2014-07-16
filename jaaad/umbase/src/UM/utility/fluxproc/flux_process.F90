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
!
! Author:     M. J. Bell
!----------------------------------------------------------------------
! contains routines: Flux_Process
!
! Purpose: Controls processing for Flux_process_main.
!----------------------------------------------------------------------
      subroutine Flux_Process (                                         &
#include "aflddims.h"
     &     ppxRecs,icode)

      implicit none

! parameters used
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"

! declaration of argument list
#include "cflddims.h"

      integer icode  ! IN/OUT error code ; > 0 => fatal error detected


! no local parameters

! declaration of globals used
#include "cunitnos.h"
#include "cmess.h"

! declaration of local arrays (all arrays in COMDECK CFIELDS)
#include "clsms.h"
#include "ccoords.h"
#include "cinterp.h"
#include "cfillin.h"
#include "crotgrd.h"
#include "cselct.h"

! local scalars
      integer IROW_NUMBER
      character*80 cmessage

! declaration of externals
      external winds, heat, moisture, sea_ice, reference, pressure
#if !defined(FLXPLPR)
      external read_lsms
#endif
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'Flux_Process'  ! subroutine name for error messages

! 0.1 Read StashMaster files

      IROW_NUMBER=0
! DEPENDS ON: getppx
      CALL GETPPX(22,2,'STASHmaster_A',IROW_NUMBER,                     &
#include "argppx.h"
     &  ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(22,2,'STASHmaster_O',IROW_NUMBER,                     &
#include "argppx.h"
     &  ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(22,2,'STASHmaster_S',IROW_NUMBER,                     &
#include "argppx.h"
     &  ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(22,2,'STASHmaster_W',IROW_NUMBER,                     &
#include "argppx.h"
     &  ICODE,CMESSAGE)

#if !defined(FLXPLPR)
! 1. Read in land sea masks and calculate grid coordinates and
!    coefficients for interpolation from atmosphere to ocean grids.

! DEPENDS ON: read_lsms
      call read_lsms (                                                  &
#include "afields.h"
#include "argppx.h"
     &                icode)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 1. error in read_lsms '
        go to 9999
      end if

#endif

! 2. Produce the output flux files
! 2.1 produce wind flux file

      if (l_winds_slt) then

! DEPENDS ON: winds
      call winds(                                                       &
#include "afields.h"
#include "argppx.h"
     &                icode)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2.1 error in winds '
        go to 9999
      end if

      end if ! l_winds_slt

! 2.2 produce heat flux file

      if (l_heat_slt) then

! DEPENDS ON: heat
      call heat(                                                        &
#include "afields.h"
#include "argppx.h"
     &                icode)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2.2 error in heat '
        go to 9999
      end if

      end if ! l_heat_slt

! 2.3 produce moisture flux file

      if (l_moisture_slt) then

! DEPENDS ON: moisture
      call moisture(                                                    &
#include "afields.h"
#include "argppx.h"
     &                icode)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2.3 error in moisture '
        go to 9999
      end if

      end if ! l_moisture_slt

! 2.4 produce sea ice flux file

      if (l_sea_ice_slt) then

! DEPENDS ON: sea_ice
      call sea_ice(                                                     &
#include "afields.h"
#include "argppx.h"
     &                icode)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2.4 error in sea_ice '
        go to 9999
      end if

      end if ! l_sea_ice_slt

! 2.5 produce reference flux file

      if (l_references_slt) then

! DEPENDS ON: reference
      call reference(                                                   &
#include "afields.h"
#include "argppx.h"
     &                icode)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2.5 error in reference '
        go to 9999
      end if

      end if ! l_references_slt

! 2.6 produce pressure flux file

      if (l_pressure_slt) then

! DEPENDS ON: pressure
      call pressure(                                                    &
#include "afields.h"
#include "argppx.h"
     &                icode)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2.6 error in pressure '
        go to 9999
      end if

      end if ! l_pressure_slt

! 2.7 produce wind speed flux file

      if (l_windspd_slt) then

! DEPENDS ON: windspd
      call windspd(                                                     &
#include "afields.h"
#include "argppx.h"
     &                icode)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2.7 error in windspd '
        go to 9999
      end if

      end if ! l_windspd_slt

9999  continue
      return
      END SUBROUTINE Flux_Process
!----------------------------------------------------------------------
#endif
