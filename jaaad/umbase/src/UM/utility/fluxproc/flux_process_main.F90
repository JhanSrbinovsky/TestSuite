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
! 6.1      30/07/04     Include parvars. A. Hines.
!
! Author:     M. J. Bell
!----------------------------------------------------------------------
! contains progam: Flux_Process_Main
!
! Purpose: Flux processing routine.
!          Takes fields files of fluxes
!          (a) output from NWP  and
!          (b) NWP or coupled model or observed "climatologies"
!          combines them into fluxes as required by ocean models
!
!          The version of the code to build is determined by which
!          of FLUXPROC, FLXPLPR or FLXPLIN is defined.
!
!          If FLUXPROC is defined, the original executable is built.
!
!          If FLXPLPR (flux parallel processing), the executable to
!          process fluxes to the original NWP grid in parallel is
!          built; interpolation to ocean grid is performed by a
!          separate program.
!
!          If FLXPLIN (flux parallel interpolation), the executable to
!          interpolate fluxes to the ocean grid in parallel is
!          built.
!----------------------------------------------------------------------
      Program Flux_Process_Main

      implicit none

#include "csubmodl.h"
#include "parvars.h"

! declaration of parameters

! declaration of globals used
#include "cmess.h"

! declaration of local scalars

#include "cflddims.h"
      integer icode  ! error code ; > 0 => fatal error detected
      integer ppxRecs
      character*80 cmessage

! declaration of routines used
      external read_control_files,                                      &
     &   read_field_headers, Flux_Process
#if !defined(FLXPLPR) && !defined(FLXPLIN)
      external open_control_files, close_files
#endif
#if defined(FLXPLPR)
      external open_flux_control_files, close_flux_files
#endif
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'Flux_Process_Main'  ! subroutine name for error messages

#if defined(FLXPLPR) && defined(T3E)
! Executable may be run in MPP mode on T3E in order to
! force the use of application PEs. Only 1 PE does the
! processing.
      IF ( MY_PE()  ==  0) then
#endif

      icode = 0   ! initialise icode

! 0.1 Initialise N_INTERNAL_MODEL/INTERNAL_MODEL_INDEX
      N_INTERNAL_MODEL=4
      INTERNAL_MODEL_INDEX(1)=1    !  Atmos
      INTERNAL_MODEL_INDEX(2)=2    !  Ocean
      INTERNAL_MODEL_INDEX(3)=3    !  Slab
      INTERNAL_MODEL_INDEX(4)=4    !  Wave

! 0.2 Read STASHmaster files
      ppxRecs=1
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(22,'STASHmaster_A',ppxRecs,ICODE,CMESSAGE)
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(22,'STASHmaster_O',ppxRecs,ICODE,CMESSAGE)
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(22,'STASHmaster_S',ppxRecs,ICODE,CMESSAGE)
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(22,'STASHmaster_W',ppxRecs,ICODE,CMESSAGE)

 ! 1. Open all control and log files
#if defined(FLXPLPR)
! DEPENDS ON: open_flux_control_files
      call open_flux_control_files( icode )
#else
! DEPENDS ON: open_control_files
      call open_control_files( icode )
#endif

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 1. Failed to open control and log files'
        go to 9999
      end if

! 2. Read all control files and open output flux files
! DEPENDS ON: read_control_files
      call read_control_files( icode )

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &' step 2. Failed to read control files or open output flux files'
        go to 9999
      end if

! 3. Open and read headers (fixed header & lookups) of flux fields
! DEPENDS ON: read_field_headers
      call read_field_headers(                                          &
#include "aflddima.h"
     &                        ppxRecs,icode )

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 3. Failed to read headers of flux fields'
        go to 9999
      end if

#if !defined(FLXPLPR)
! 4. Open and read lookup tables of land-sea masks and find
!    field dimensions
! DEPENDS ON: read_lsm_headers
      call read_lsm_headers(                                            &
#include "aflddims.h"
     &    ppxRecs,icode)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 4. Failed to read lookups of lsms'
        go to 9999
      end if

#else
! Set dummy values for unused variables
      ncolsO=1
      nrowstO=1
      nrowsuO=1
#endif

! 5. Do main processing at a lower level
! DEPENDS ON: flux_process
      call Flux_Process(                                                &
#include "aflddims.h"
     &     ppxRecs,icode)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 4. Failed doing main flux processing'
        go to 9999
      end if

! 6. close files opened in steps 1. - 3.

#if defined(FLXPLPR)
! DEPENDS ON: close_flux_files
      call close_flux_files
#else
! DEPENDS ON: close_files
      call close_files
#endif

9999  continue
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,Csub,                                        &
     &   'Flux Processing failed with error code = ',icode
        close ( UnErr )
      endif

#if defined(FLXPLPR) && defined(T3E)
      endif
#endif
      stop
      END PROGRAM Flux_Process_Main

!----------------------------------------------------------------------
#endif
