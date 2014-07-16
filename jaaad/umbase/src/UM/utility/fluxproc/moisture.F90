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
! contains routines: moisture
!
! Purpose: Flux processing routine.
!          To produce a pp file containing:
!          Precipitation less evapration :-
!           calculated from input rainfall and snowfall fields.
!----------------------------------------------------------------------
      subroutine moisture(                                              &
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
#include "cunitnos.h"
#include "cmess.h"
#include "c_mdi.h"
#include "cvaloff.h"
#include "cdebug.h"
#include "creftim.h"
#include "cvaltim.h"

! declaration of local arrays
      integer Int_Head_evap(Len_IntHd)   ! integer part of lookup
                                         ! (evap)
      integer Int_Head_drain(Len_IntHd)  ! integer part of lookup
                                         ! (drain)
      integer Int_Head_convrain(Len_IntHd)! integer part of lookup
                                         ! (crain)
      integer Int_Head_dsnow(Len_IntHd)  ! integer part of lookup
                                         ! (dsnow)
      integer Int_Head_convsnow(Len_IntHd)! integer part of lookup
                                          ! (csnow)
      real Real_Head_evap(Len_RealHd)    ! real part of lookup (evap)
      real Real_Head_drain(Len_RealHd)   ! real part of lookup (drain)
      real Real_Head_convrain(Len_RealHd)! real part of lookup (crain)
      real Real_Head_dsnow(Len_RealHd)   ! real part of lookup (dsnow)
      real Real_Head_convsnow(Len_RealHd)! real part of lookup (csnow)
      real evaporation(ncols, nrowst)   ! evaporation field
      real dynamic_rain(ncols, nrowst)   ! large scale rain field
      real conv_rain(ncols,nrowst)      ! convective rain field
      real dynamic_snow(ncols, nrowst)   ! large scale snow field
      real conv_snow(ncols,nrowst)      ! convective snow field
      real Precip_less_evap(ncols,nrowst)! PLE field
      real fieldint(ncols,nrowst)        ! intermediate field

! declaration of local scalars

      integer ivt           ! loop index over validity times
      integer iadd          ! loop index over additional times
      integer IVTOffHr      ! offset of validity time from reference
      integer IOutUnit      ! output unit

      logical ldebug        ! T => output debugging info (set in 0.)
      logical l_leads       ! T => using minleadsfrac
                            ! F => using minicefrac
      logical lcalcprev     ! T => field has already been found for
                            !      additional time

      character * 256 cmessage   ! error message

! declaration of externals
      external read_leads_flds, read_accum_flds, write_one_field

!----------------------------------------------------------------------
! 0. Preliminaries
!----------------------------------------------------------------------
      CSub = 'moisture'  ! subroutine name for error messages

      ldebug = l_moisture_dbg    ! set by debug input control file

!----------------------------------------------------------------------
! 1. start loop over validity times
!----------------------------------------------------------------------
      do ivt = 1, NoValidTimes

        IVTOffHr = IValidOffHr(ivt)
        IOutUnit = IOutUnitOff(ivt) + UnitMoistureOut

!----------------------------------------------------------------------
! 2. Read in evaporation field
!----------------------------------------------------------------------
        lcalcprev = .false.
        if ( ivt  >   1 ) then
          do iadd = 1,NoAddTimesPreferred
            if ( IVTOffHr  ==  INewOffHrPreferred(iadd) ) then
              lcalcprev = .true.
            endif
          enddo
        endif
        if ( .not. lcalcprev ) then
          l_leads = .true.      ! set to true to use minleadsfrac
! DEPENDS ON: read_leads_flds
          call read_leads_flds (StCEvaporation,StCAICE,                 &
     &                    IVTOffHr, ldebug,                             &
     &                    l_leads,Int_Head_evap,                        &
     &                    Real_Head_evap, ncols, nrowst,                &
     &                    evaporation,                                  &
#include "argppx.h"
     &                    icode)

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' step 2. unable to read evaporation field'
            icode = 1008
            go to 9999
          end if

!----------------------------------------------------------------------
! 3. Read in large scale rain amount
!----------------------------------------------------------------------
! DEPENDS ON: read_accum_flds
          call read_accum_flds(StCdrain, IVTOffHr,                      &
     &               ldebug, Int_Head_drain,                            &
     &               Real_Head_drain,                                   &
     &               ncols, nrowst,                                     &
     &               dynamic_rain,                                      &
#include "argppx.h"
     &               icode)

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' step 3. unable to read dynamic rain'
            icode = 1009
            go to 9999
          end if

!----------------------------------------------------------------------
! 4. Calculate first part of PLE (dynamic_rain - evaporation)
!----------------------------------------------------------------------
! DEPENDS ON: fieldsub
          call FieldSub (ncols, nrowst, rmdi,                           &
     &            dynamic_rain, evaporation,                            &
     &            fieldint,                                             &
     &            icode, cmessage)
!----------------------------------------------------------------------
! 5. Read in covective rain field
!----------------------------------------------------------------------
! DEPENDS ON: read_accum_flds
          call read_accum_flds(StCconvrain, IVTOffHr,                   &
     &               ldebug, Int_Head_convrain,                         &
     &               Real_Head_convrain,                                &
     &               ncols, nrowst,                                     &
     &               conv_rain,                                         &
#include "argppx.h"
     &               icode)

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' step 5. unable to read convective rain'
            icode = 1010
            go to 9999
          end if

!----------------------------------------------------------------------
! 6. Continue PLE calculation (PLE = PLE + Conv_Rain)
!----------------------------------------------------------------------
! DEPENDS ON: fieldadd
          call FieldAdd(ncols, nrowst, rmdi,                            &
     &            fieldint, conv_rain,                                  &
     &            Precip_less_evap,                                     &
     &            icode, cmessage)

!----------------------------------------------------------------------
! 7. Read in large scale snow amount
!----------------------------------------------------------------------
! DEPENDS ON: read_accum_flds
          call read_accum_flds(StCdsnow, IVTOffHr,                      &
     &               ldebug, Int_Head_dsnow,                            &
     &               Real_Head_dsnow,                                   &
     &               ncols, nrowst,                                     &
     &               dynamic_snow,                                      &
#include "argppx.h"
     &               icode)

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' step 7. unable to read large scale snow field'
            icode = 1011
            go to 9999
          end if

!----------------------------------------------------------------------
! 8. Continue PLE calculation (PLE = PLE + dynamic_snow)
!----------------------------------------------------------------------
! DEPENDS ON: fieldadd
          call FieldAdd(ncols, nrowst, rmdi,                            &
     &            Precip_less_evap, dynamic_snow,                       &
     &            fieldint,                                             &
     &            icode, cmessage)


!----------------------------------------------------------------------
! 9. Read in convective snow field
!----------------------------------------------------------------------
! DEPENDS ON: read_accum_flds
          call read_accum_flds(StCconvsnow, IVTOffHr,                   &
     &               ldebug, Int_Head_convsnow,                         &
     &               Real_Head_convsnow,                                &
     &               ncols, nrowst,                                     &
     &               conv_snow,                                         &
#include "argppx.h"
     &               icode)

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' step 9. unable to read convective snow field'
            icode = 1012
            go to 9999
          end if

!----------------------------------------------------------------------
! 10. Final PLE calculation (PLE = PLE + conv_snow)
!----------------------------------------------------------------------
! DEPENDS ON: fieldadd
          call FieldAdd(ncols, nrowst, rmdi,                            &
     &            fieldint, conv_snow,                                  &
     &            Precip_less_evap,                                     &
     &            icode, cmessage)
!----------------------------------------------------------------------
! 11. Write out Precipitation less Evaporation
!----------------------------------------------------------------------
! DEPENDS ON: write_one_field
          call write_one_field (                                        &
#include "afields.h"
     &       OutStCPLE, FFPLE, PPPLE, IVTOffHr,                         &
     &       ITGrid, nrowst,                                            &
     &       Int_Head_convsnow, Real_Head_convsnow, IOutUnit, ldebug,   &
     &       Precip_less_evap, icode)
          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' step 11. unable to write PLE field'
            icode = 1106
            go to 9999
          end if
        else
! DEPENDS ON: add_hours
          call add_hours(                                               &
#include "areftim.h"
#include "avaltim.h"
     &       IVTOffHr)
! DEPENDS ON: amend_times
          call amend_times (                                            &
#include "avaltim.h"
     &                   Int_Head_convsnow,Len_IntHd )
! DEPENDS ON: write_one_field
          call write_one_field (                                        &
#include "afields.h"
     &       OutStCPLE, FFPLE, PPPLE, IVTOffHr,                         &
     &       ITGrid, nrowst,                                            &
     &       Int_Head_convsnow, Real_Head_convsnow, IOutUnit, ldebug,   &

     &       Precip_less_evap, icode)
          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' step 11. unable to write PLE field'
            icode = 1106
            go to 9999
          end if
        endif   ! .not. lcalcprev
!----------------------------------------------------------------------
! 12. end loop over validity times
!----------------------------------------------------------------------
        enddo    !  ivt

9999  continue
      return
      END SUBROUTINE moisture
!----------------------------------------------------------------------
#endif
