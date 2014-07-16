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
! 5.3      22/10/01     New if defined added. A. Hines
!
! Author:     M. J. Bell
!----------------------------------------------------------------------
! contains routines: heat
!
! Purpose: Flux processing routine.
!          To produce a pp file containing:
!          Net Penetraing Solar Radiation (SOL)
!          Net non Penetraing Heat        (HTN)
!          for the times required.
!----------------------------------------------------------------------
      subroutine heat(                                                  &
#include "afields.h"
#include "argppx.h"
     &                 icode )

      implicit none

! declaration of argument list

! array dimensions, lsms, interpolation coeffs etc. : all intent IN
#include "cfields.h"

      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! declaration of parameters
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "cfdcodes.h"
#include "plookups.h"

      real lhevap               ! latent heat of evaporation
      parameter ( lhevap    = 2.25E6)

! declaration of globals used
#include "cunitnos.h"
#include "c_mdi.h"
#include "cmess.h"
#include "cvaloff.h"
#include "cdebug.h"

! declaration of local arrays
      integer Int_Head_SW1(Len_IntHd)   ! integer part of lookup table
      integer Int_Head_SW(Len_IntHd)    ! integer part of lookup table
      integer Int_Head_LW(Len_IntHd)    ! integer part of lookup table
      integer Int_Head_EVAP(Len_IntHd)  ! integer part of lookup table
      integer Int_Head_SH(Len_IntHd)    ! integer part of lookup table
      real Real_Head_SW1(Len_RealHd)    ! real part of lookup table
      real Real_Head_SW(Len_RealHd)     ! real part of lookup table
      real Real_Head_LW(Len_RealHd)     ! real part of lookup table
      real Real_Head_EVAP(Len_RealHd)   ! real part of lookup table
      real Real_Head_SH(Len_RealHd)     ! real part of lookup table
      real SW_radiation_band1(ncols, nrowst)! short wave flux (band 1)
      real SW_radiation(ncols, nrowst)      ! short wave flux
      real LW_radiation(ncols, nrowst)      ! long_wave_radiation
      real evaporation(ncols, nrowst)       ! evaporation
      real sensible_heat(ncols, nrowst)     ! sensible heat
      real latent_heat(ncols,nrowst)        ! latent heat
      real non_pen_heat(ncols,nrowst)       ! net non-penetrating heat
      real fieldint(ncols,nrowst)           ! intermediate field
! declaration of local scalars

      integer ivt           ! loop index over validity times
      integer IVTOffHr      ! offset of validity time from reference
      integer IOutUnit      ! output unit

      logical ldebug        ! T => output debugging info (set in 0.)
      logical l_leads       ! T => using minleadsfrac
                            ! F => using minicefrac

      character * 256 cmessage   ! error message

! declaration of externals
      external read_leads_flds, write_one_field,                        &
     &         ScalarMult,FieldSub,FieldAdd

!----------------------------------------------------------------------
! 0. Preliminaries
!----------------------------------------------------------------------
      CSub = 'heat'  ! subroutine name for error messages

      ldebug = l_heat_dbg      ! set by debug input control file

!----------------------------------------------------------------------
! 1. start loop over validity times
!----------------------------------------------------------------------
      do ivt = 1, NoValidTimes

        IVTOffHr = IValidOffHr(ivt)
        IOutUnit = IOutUnitOff(ivt) + UnitHeatOut

!----------------------------------------------------------------------
! 2. Read in net down short wave flux over open sea (band 1)
!----------------------------------------------------------------------
        l_leads = .true.                ! set to use minleadsfrac
! DEPENDS ON: read_leads_flds
        call read_leads_flds(StCSW1,StCAICE,                            &
     &                    IVTOffHr, ldebug,                             &
     &                    l_leads,Int_Head_SW1,                         &
     &                    Real_Head_SW1, ncols, nrowst,                 &
     &                    SW_radiation_band1,                           &
#include "argppx.h"
     &                    icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2. unable to read SW Radiation Flux (band 1)'
          icode = 1001
          go to 9999
        end if

! 2.2 Write out solar radiation
! DEPENDS ON: write_one_field
        call write_one_field (                                          &
#include "afields.h"
     &       OutStCSOL, FFSOL, PPSOL, IVTOffHr,                         &
     &       ITGrid, nrowst,                                            &
     &       Int_Head_SW1, Real_Head_SW1, IOutUnit,ldebug,              &
     &       SW_radiation_band1, icode)
        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2. unable to write penetrating '                    &
     &       ,'solar radiation (SOL)'
          icode = 1101
          go to 9999
        end if

!----------------------------------------------------------------------
! 3. Read in fields to calculate net non penetrating heat
!----------------------------------------------------------------------
! 3.1 Read net down short wave readiation
! DEPENDS ON: read_leads_flds
        call read_leads_flds(StCSW,StCAICE,                             &
     &                    IVTOffHr, ldebug,                             &
     &                    l_leads,Int_Head_SW,                          &
     &                    Real_Head_SW, ncols, nrowst,                  &
     &                    SW_radiation,                                 &
#include "argppx.h"
     &                    icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3. unable to read SW Radiation Flux'
          icode = 1002
          go to 9999
        end if
! 3.2 Use Field Sub to work out first component of HTN
! DEPENDS ON: fieldsub
        call FieldSub (ncols, nrowst, rmdi,                             &
     &            SW_radiation, SW_radiation_band1,                     &
     &            fieldint,                                             &
     &            icode, cmessage)
! 3.3 Read net down long wave flux
! DEPENDS ON: read_leads_flds
        call read_leads_flds(StCLongWave,StCAICE,                       &
     &                    IVTOffHr, ldebug,                             &
     &                    l_leads,Int_Head_LW,                          &
     &                    Real_Head_LW, ncols, nrowst,                  &
     &                    LW_radiation,                                 &
#include "argppx.h"
     &                    icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3. unable to read LW Radiation Flux'
          icode = 1003
          go to 9999
        end if

! 3.4 Use FieldAdd to do HTN = fieldint + LW_radiation
! DEPENDS ON: fieldadd
        call FieldAdd (ncols, nrowst, rmdi,                             &
     &            fieldint, LW_radiation,                               &
     &            non_pen_heat,                                         &
     &            icode, cmessage)

! 3.5 Read evaporation from sea
! DEPENDS ON: read_leads_flds
        call read_leads_flds(StCEvaporation,StCAICE,                    &
     &                    IVTOffHr, ldebug,                             &
     &                    l_leads,Int_Head_EVAP,                        &
     &                    Real_Head_EVAP, ncols, nrowst,                &
     &                    evaporation,                                  &
#include "argppx.h"
     &                    icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3. unable to read evaporation from sea'
          icode = 1004
          go to 9999
        end if

! 3.6 Use ScalarMult and FieldSub to work out latent heat
!     and subtract it from HTN
! DEPENDS ON: scalarmult
        call ScalarMult (ncols, nrowst, rmdi,                           &
     &            lhevap, evaporation,                                  &
     &            latent_heat,                                          &
     &            icode, cmessage)
! DEPENDS ON: fieldsub
        call FieldSub (ncols, nrowst, rmdi,                             &
     &            non_pen_heat, latent_heat,                            &
     &            fieldint,                                             &
     &            icode, cmessage)

! 3.7 Read Sensible Heat Flux
! DEPENDS ON: read_leads_flds
        call read_leads_flds(StCSensibleHeat,StCAICE,                   &
     &                    IVTOffHr, ldebug,                             &
     &                    l_leads,Int_Head_SH,                          &
     &                    Real_Head_SH, ncols, nrowst,                  &
     &                    sensible_heat,                                &
#include "argppx.h"
     &                    icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3. unable to read sensible heat flux'
          icode = 1005
          go to 9999
        end if

! 3.8 Use FieldSub to calculate final HTN
! DEPENDS ON: fieldsub
        call FieldSub (ncols, nrowst, rmdi,                             &
     &            fieldint, sensible_heat,                              &
     &            non_pen_heat,                                         &
     &            icode, cmessage)

! 3.9 Write out net non penetrating heat
! DEPENDS ON: write_one_field
        call write_one_field (                                          &
#include "afields.h"
     &       OutStCHTN, FFHTN, PPHTN, IVTOffHr,                         &
     &       ITGrid, nrowst,                                            &
     &       Int_Head_SH, Real_Head_SH, IOutUnit,ldebug,                &

     &       non_pen_heat, icode)
        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3. unable to write net non penetrating heat'
          icode = 1102
          go to 9999
        end if

!----------------------------------------------------------------------
! 4. end loop over validity times
!----------------------------------------------------------------------
        enddo    !  ivt

9999  continue
      return
      END SUBROUTINE heat
!----------------------------------------------------------------------
#endif
