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
! contains routines: reference
!
! Purpose: Flux processing routine.
!          To produce a pp file containing:
!             Reference Sea Surface temperature
!             Reference Sea Surface Salinity
!             Reference Ice Depth
!          for the times required.
!          Change to not output data if climatology
!          doesn't exist (S. Spall)
!----------------------------------------------------------------------
      subroutine reference(                                             &
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

! declaration of globals used
#include "clookadd.h"
#include "cunitnos.h"
#include "cmess.h"
#include "c_mdi.h"
#include "cvaloff.h"
#include "cdebug.h"
#include "c_0_dg_c.h"

! declaration of local arrays
      integer Int_Head_SST(Len_IntHd)  ! integer part of lookup table
      integer Int_Head_SSS(Len_IntHd)  ! integer part of lookup table
      integer Int_Head_HICE(Len_IntHd)  ! integer part of lookup table
      integer Int_Head_ICEFRAC(Len_IntHd) ! integer part of lookup table
      real Real_Head_SST(Len_RealHd)   ! real part of lookup table
      real Real_Head_SSS(Len_RealHd)   ! real part of lookup table
      real Real_Head_HICE(Len_RealHd)   ! real part of lookup table
      real Real_Head_ICEFRAC(Len_RealHd)   ! real part of lookup table
      real ref_sea_surface_temp(ncols, nrowst) ! ref SST
      real ref_sea_surface_salin(ncols, nrowst)! ref SSS
      real ref_ice_depth(ncols,nrowst)         ! reference ice depth
      real ice_depth(ncols,nrowst)             ! ice depth
      real icefrac(ncols,nrowst)               ! ice fraction

! declaration of local scalars

      integer ivt           ! loop index over validity times
      integer IVTOffHr      ! offset of validity time from reference
      integer IOutUnit      ! output unit

      logical ldebug        ! T => output debugging info (set in 0.)

      real Real_Add_value   ! real value to add field (of SSTs)
      real salinity_factor  ! to convert SSS from g/kg to kg/kg
      real salinity_offset  ! to take 0.035 from all salinity values

      parameter ( salinity_factor = 0.001 )
      parameter ( salinity_offset = -0.035 )

      character * 256 cmessage   ! error message


! declaration of externals
      external read_fields, write_one_field

!----------------------------------------------------------------------
! 0. Preliminaries
!----------------------------------------------------------------------
      CSub = 'reference'  ! subroutine name for error messages

      ldebug = l_references_dbg     ! set by debug input control file

!----------------------------------------------------------------------
! 1. start loop over validity times
!----------------------------------------------------------------------
      do ivt = 1, NoValidTimes

        IVTOffHr = IValidOffHr(ivt)
        IOutUnit = IOutUnitOff(ivt) + UnitReferencesOut

!----------------------------------------------------------------------
! 2. read in reference sea surface temperature
!----------------------------------------------------------------------
! DEPENDS ON: read_fields
        call read_fields(StCSST, IVTOffHr,                              &
     &               ldebug, Int_Head_SST, Real_Head_SST,               &
     &               ncols, nrowst,                                     &
     &               ref_sea_surface_temp,                              &
#include "argppx.h"
     &               icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2. unable to read reference SST'
          icode = 1013
          go to 9999
        end if

! 2.1 Change field from Kelvin to Celsius
        Real_Head_SST(BDATUM - Len_IntHd) = - ZERODEGC
        if ( Real_Head_SST(BDATUM - Len_IntHd)  /=  rmdi .and.          &
     &       Real_Head_SST(BDATUM - Len_IntHd)  /=  0.0 ) then
          Real_Add_value = Real_Head_SST(BDATUM - Len_IntHd)
! DEPENDS ON: scalaradd
          call ScalarAdd(ncols, nrowst, rmdi,                           &
     &                 Real_Add_value,                                  &
     &                 ref_sea_surface_temp,                            &
     &                 ref_sea_surface_temp, icode, cmessage)
        endif

! 2.2 Reset SSTs less than -1.8 deg C to -1.8 deg C
! DEPENDS ON: set_sst
        call set_sst(ncols, nrowst,                                     &
     &            ref_sea_surface_temp, rmdi,                           &
     &            ref_sea_surface_temp)

!----------------------------------------------------------------------
! 3. Write out reference SST
!----------------------------------------------------------------------
! DEPENDS ON: write_one_field
        call write_one_field (                                          &
#include "afields.h"
     &       OutStCSST, FFSST, PPSST, IVTOffHr,                         &
     &       ITGrid, nrowst,                                            &
     &       Int_Head_SST, Real_Head_SST, IOutUnit, ldebug,             &
     &       ref_sea_surface_temp, icode)
        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3. unable to write reference SST'
          icode = 1107
          go to 9999
        end if

!----------------------------------------------------------------------
! 4. read in reference sea surface salinity from climatology
!----------------------------------------------------------------------
        if (LClimate) then

! 4.1 read_climate field by calling read_climate field
! DEPENDS ON: read_climate_field
           call read_climate_field(StCSSS, IVTOffHr,                    &
     &           ldebug, Int_Head_SSS, Real_Head_SSS,                   &
     &           ncols, nrowst,                                         &
     &           ref_sea_surface_salin,                                 &
#include "argppx.h"
     &           icode)

           if ( icode  <=  0) then
             write(UnStd,*)CStd//CSub//'4. climate field extracted  ',  &
     &        ' for stash code =', StCSSS, '; IVTOffHr = ', IVTOffHr
           else

             write(UnErr,*)CErr//CSub//                                 &
     &        '4. failed to retrieve climate field ',                   &
     &        ' for stash code =', StCSSS, '; IVTOffHr = ', IVTOffHr
             icode = 1014
             goto 9999
           end if

! 4.2 Convert salinity units from g/kg to kg/kg
! DEPENDS ON: scalarmult
        call ScalarMult (ncols, nrowst, rmdi,                           &
     &            salinity_factor,                                      &
     &            ref_sea_surface_salin,                                &
     &            ref_sea_surface_salin,                                &
     &            icode, cmessage)

! 4.3 Subtract 0.035 from each salinity element in field
! DEPENDS ON: scalaradd
        call ScalarAdd                                                  &
     &           (ncols, nrowst, rmdi,                                  &
     &            salinity_offset, ref_sea_surface_salin,               &
     &            ref_sea_surface_salin,                                &
     &            icode, cmessage)

! DEPENDS ON: write_one_field
        call write_one_field (                                          &
#include "afields.h"
     &       OutStCSSS, FFSSS, PPSSS, IVTOffHr,                         &
     &       ITGrid, nrowst,                                            &
     &       Int_Head_SSS, Real_Head_SSS, IOutUnit,ldebug,              &
     &       ref_sea_surface_salin, icode)
        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 4. unable to write reference SSS'
          icode = 1108
          go to 9999
        end if

        end if !  LClimate

!----------------------------------------------------------------------
! 5. read in reference ice depth from climatology
!----------------------------------------------------------------------
        if (LClimate) then

! 5.1 read_climate field by calling read_climate field
! DEPENDS ON: read_climate_field
           call read_climate_field(StCHICE, IVTOffHr,                   &
     &           ldebug, Int_Head_HICE, Real_Head_HICE,                 &
     &           ncols, nrowst,                                         &
     &           ref_ice_depth,                                         &
#include "argppx.h"
     &           icode)

           if ( icode  <=  0) then
             write(UnStd,*)CStd//CSub//'5. climate field extracted  ',  &
     &        ' for stash code =', StCHICE, '; IVTOffHr = ', IVTOffHr
           else

             write(UnErr,*)CErr//CSub//                                 &
     &        '5. failed to retrieve climate field ',                   &
     &        ' for stash code =', StCHICE, '; IVTOffHr = ', IVTOffHr
             icode = 1015
             goto 9999
           end if

!----------------------------------------------------------------------
! 6. Read in ice fraction
!----------------------------------------------------------------------
! DEPENDS ON: read_fields
        call read_fields(StCAICE, IVTOffHr,                             &
     &               ldebug, Int_Head_ICEFRAC, Real_Head_ICEFRAC,       &
     &               ncols, nrowst,                                     &
     &               icefrac,                                           &
#include "argppx.h"
     &               icode)

        if ( icode  <=  0) then
          write(UnStd,*)CStd//CSub//'6. ice fraction extracted  ',      &
     &     ' for stash code =', StCAICE, '; IVTOffHr = ', IVTOffHr
        else

          write(UnErr,*)CErr//CSub//                                    &
     &     '6. failed to retrieve ice fraction field ',                 &
     &     ' for stash code =', StCAICE, '; IVTOffHr = ', IVTOffHr
          icode = 1016
          goto 9999
        end if

!----------------------------------------------------------------------
! 7. Use FieldMult to calculate HICE
!----------------------------------------------------------------------
! DEPENDS ON: fieldmult
        call FieldMult (ncols, nrowst, rmdi,                            &
     &            ref_ice_depth, icefrac,                               &
     &            ice_depth,                                            &
     &            icode, cmessage)
!----------------------------------------------------------------------
! 8. Write out HICE
!----------------------------------------------------------------------
! DEPENDS ON: write_one_field
        call write_one_field (                                          &
#include "afields.h"
     &       OutStCHICE, FFHICE, PPHICE, IVTOffHr,                      &
     &       ITGrid, nrowst,                                            &
     &       Int_Head_HICE, Real_Head_HICE, IOutUnit, ldebug,           &
     &       ice_depth, icode)
        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 8. unable to write ice_depth'
          icode = 1109
          go to 9999
        end if

        end if !  LClimate

!----------------------------------------------------------------------
! 9. end loop over validity times
!----------------------------------------------------------------------
        enddo    !  ivt

9999  continue
      return
      END SUBROUTINE reference
!----------------------------------------------------------------------
#endif
