#if defined(FLXPLIN)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! contains progam: Flux_Transform_Main
!
! Purpose: Flux processing routine.
!          Takes pp files of fluxes as input, interpolates them
!          from an atmosphere to an ocean grid and fills in
!          missing data values
!
!    Model            Modification history:
!   version  Date
!    5.3  15/10/01  New deck. A. Hines
! 6.1      29/07/04     Allow NoInFiles to be specified as
!                       an environment variable. A. Hines
!
!    Programming standard :
!
!    Logical components covered :
!
!    System task:
!
!    External documentation:
!----------------------------------------------------------------------
      Program Flux_Transform_Main

      implicit none

#include "csubmodl.h"
#include "parvars.h"

! declaration of parameters

! declaration of globals used
#include "cmess.h"
#include "cunitnos.h"

! declaration of local scalars
#include "cflddims.h"
      integer icode  ! error code ; > 0 => fatal error detected
      integer ppxRecs
      integer iunit  ! loop counter
      integer ifile  ! loop counter
      integer NoInFiles ! number of input files
      character*80 cmessage
      Character*4  C_NoInFiles  ! Char variable to read env var

! declaration of routines used
      external open_grid_control_files, read_debug_cntl,                &
     &   read_lsm_headers, Flux_Transform_Grid,                         &
     &   close_grid_files
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'Flux_Transform_Main'  ! subroutine name for err messages

#if defined(T3E)
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

! 0.3 Read NoInFiles

      CALL FORT_GET_ENV('NO_FLUX_FILES',13,C_NoInFiles,4,icode)
      IF (icode  /=  0) THEN
        WRITE(6,*) 'Warning : Environment variable NO_FLUX_FILES ',     &
     &             'has not been set.'
        WRITE(6,*) 'Setting NoInFiles to 6'
        NoInFiles=6
      ELSE
        READ(C_NoInFiles,'(I2)') NoInFiles
        write (6,*) ' '
        write (6,*) 'NO_FLUX_FILES is set to ',NoInFiles
      ENDIF
 ! 1. Open all control and log files
! DEPENDS ON: open_grid_control_files
      call open_grid_control_files( icode )

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 1. Failed to open control and log files'
        go to 9999
      end if

! 2. Read debug control file and open debug ouput file
! DEPENDS ON: read_debug_cntl
      call read_debug_cntl( icode )

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &' step 2. Failed to read debug control file'
        go to 9999
      end if

! 3. Open and read lookup tables of land-sea masks and find
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

! 4. Open  input and output flux files

! 4.0 Set all logicals showing which files are open to false
      do iunit = IUnOutLow, IUnOutHi
        LUnOutOpen(iunit) = .False.
      end do

! 4.1 Open input flux file  ! should this status be 'old' or 'share' ?
      do ifile = 0, NoInFiles-1
        LUnOutOpen(IUnOutLow+10+ifile)= .True.
! DEPENDS ON: open_file
        call open_file (IUnOutLow+10+ifile,                             &
     &                  'unformatted', 'unknown', icode )
      enddo

! 4.2 Open output flux file
      LUnOutOpen(IUnOutLow+10+NoInFiles)= .True.
! DEPENDS ON: open_file
      call open_file (IUnOutLow+10+NoInFiles,                           &
     &                  'unformatted', 'unknown', icode )

! 5. Do main processing at a lower level
! DEPENDS ON: flux_transform_grid
      call Flux_Transform_Grid(                                         &
#include "aflddims.h"
     &     ppxRecs,NoInFiles,icode)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 5. Failed doing main processing'
        go to 9999
      end if

! 6. close files opened in steps 1. - 4.

! DEPENDS ON: close_grid_files
      call close_grid_files

9999  continue
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,Csub,                                        &
     &   'Flux Processing failed with error code = ',icode
        close ( UnErr )
      endif

#if defined(T3E)
      endif
#endif

      stop
      END PROGRAM Flux_Transform_Main
!----------------------------------------------------------------------
#endif
