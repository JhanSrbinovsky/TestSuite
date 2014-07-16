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
! contains routines: sea_ice
!
! Purpose: Flux processing routine.
!          To produce a pp field containing:
!            Snowfall rate
!            Sublimation rate
!            Topmelt
!            Bottom melt
!          for each of the fields required
!----------------------------------------------------------------------
      subroutine sea_ice(                                               &
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
      integer Int_Head_drain(Len_IntHd)  ! integer part of lookup
                                         ! (drain)
      integer Int_Head_convrain(Len_IntHd)! integer part of lookup
                                          ! (crain)
      integer Int_Head_dsnow(Len_IntHd)  ! integer part of lookup
                                         ! (dsnow)
      integer Int_Head_convsnow(Len_IntHd)! integer part of lookup
                                          ! (csnow)
      integer Int_Head_subrate(Len_IntHd)! integer part of lookup
                                         ! (subrate)
      integer Int_Head_topmelt(Len_IntHd)! integer part of lookup
                                         ! (topmelt)
      integer Int_Head_botmelt(Len_IntHd)! integer part of lookup
                                         ! (botmelt)
      real Real_Head_drain(Len_RealHd)   ! real part of lookup (drain)
      real Real_Head_convrain(Len_RealHd)! real part of lookup (crain)
      real Real_Head_dsnow(Len_RealHd)   ! real part of lookup (dsnow)
      real Real_Head_convsnow(Len_RealHd)! real part of lookup (csnow)
      real Real_Head_subrate(Len_RealHd)! real part of lookup (subrate)
      real Real_Head_topmelt(Len_RealHd)! real part of lookup (topmelt)
      real Real_Head_botmelt(Len_RealHd)! real part of lookup (botmelt)
      real dynamic_rain(ncols, nrowst)   ! large scale rain field
      real conv_rain(ncols,nrowst)      ! convective rain field
      real dynamic_snow(ncols, nrowst)   ! large scale snow field
      real conv_snow(ncols,nrowst)      ! convective snow field
      real fieldint(ncols,nrowst)       ! intermediate field
      real total_snow_rate(ncols,nrowst)! total snow rate field
      real sublimation_rate(ncols,nrowst) ! sublimation rate
      real topmelt(ncols,nrowst)          ! top melt
      real bottommelt(ncols,nrowst)       ! bottom melt

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
      CSub = 'sea_ice'  ! subroutine name for error messages

      ldebug = l_sea_ice_dbg     ! set by debug input control file

!----------------------------------------------------------------------
! 1. start loop over validity times
!----------------------------------------------------------------------
      do ivt = 1, NoValidTimes

        IVTOffHr = IValidOffHr(ivt)
        IOutUnit = IOutUnitOff(ivt) + UnitSeaIceOut

!----------------------------------------------------------------------
! 2. Read in large scale rain amount
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
     &       ' step 2. unable to read dynamic rain'
            icode = 1017
            go to 9999
          end if
!----------------------------------------------------------------------
! 3. Read in convective rain field
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
     &       ' step 3. unable to read convective rain'
            icode = 1018
            go to 9999
          end if
!----------------------------------------------------------------------
! 4. Start Snowfall Rate Calculation (SNO = DRAIN + CRAIN)
!----------------------------------------------------------------------
! DEPENDS ON: fieldadd
          call FieldAdd(ncols, nrowst, rmdi,                            &
     &            dynamic_rain, conv_rain,                              &
     &            total_snow_rate,                                      &
     &            icode, cmessage)
!----------------------------------------------------------------------
! 5. Read in large scale snow amount
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
     &       ' step 5. unable to read large scale snow field'
            icode = 1019
            go to 9999
          end if
!----------------------------------------------------------------------
! 6. Continue SNO calculation (SNO = SNO + dynamic_snow)
!----------------------------------------------------------------------
! DEPENDS ON: fieldadd
          call FieldAdd(ncols, nrowst, rmdi,                            &
     &            total_snow_rate, dynamic_snow,                        &
     &            fieldint,                                             &
     &            icode, cmessage)

!----------------------------------------------------------------------
! 7. Read in convective snow field
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
     &       ' step 7. unable to read convective snow field'
            icode = 1020
            go to 9999
          end if
!----------------------------------------------------------------------
! 8. Final SNO calculation (SNO = SNO + conv_snow)
!----------------------------------------------------------------------
! DEPENDS ON: fieldadd
          call FieldAdd(ncols, nrowst, rmdi,                            &
     &            fieldint, conv_snow,                                  &
     &            total_snow_rate,                                      &
     &            icode, cmessage)
!----------------------------------------------------------------------
! 9. Write out Total Snow Rate
!----------------------------------------------------------------------
! DEPENDS ON: write_one_field
          call write_one_field (                                        &
#include "afields.h"
     &       OutStCSNO, FFSNO, PPSNO, IVTOffHr,                         &
     &       ITGrid, nrowst,                                            &
     &       Int_Head_convsnow, Real_Head_convsnow, IOutUnit, ldebug,   &
     &       Total_snow_rate, icode)
          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' step 9. unable to write SNO field'
            icode = 1110
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
     &       OutStCSNO, FFSNO, PPSNO, IVTOffHr,                         &
     &       ITGrid, nrowst,                                            &
     &       Int_Head_convsnow, Real_Head_convsnow, IOutUnit, ldebug,   &
     &       Total_snow_rate, icode)
          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' step 9. unable to write SNO field'
            icode = 1110
            go to 9999
          end if
        endif   ! .not. lcalcprev
!----------------------------------------------------------------------
! 10. Read in Sublimation Rate
!----------------------------------------------------------------------
        l_leads = .false.                ! set to use minicefrac
! DEPENDS ON: read_leads_flds
        call read_leads_flds(StCSublim,StCAICE,                         &
     &                    IVTOffHr, ldebug,                             &
     &                    l_leads,Int_Head_subrate,                     &
     &                    Real_Head_subrate, ncols, nrowst,             &
     &                    sublimation_rate,                             &
#include "argppx.h"
     &                    icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 10. unable to read sublimation rate'
          icode = 1021
          go to 9999
        end if
!----------------------------------------------------------------------
! 11. Write out Sublimation rate
!----------------------------------------------------------------------
! DEPENDS ON: write_one_field
        call write_one_field (                                          &
#include "afields.h"
     &       OutStCSUB, FFSUB, PPSUB, IVTOffHr,                         &
     &       ITGrid, nrowst,                                            &
     &       Int_Head_subrate, Real_Head_subrate, IOutUnit, ldebug,     &
     &       sublimation_rate, icode)
        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 11. unable to write sublimation rate'
          icode = 1111
          go to 9999
        end if
!----------------------------------------------------------------------
! 12. Read in Topmelt Rate
!----------------------------------------------------------------------
        l_leads = .false.                ! set to use minicefrac
! DEPENDS ON: read_leads_flds
        call read_leads_flds(StCTopmelt,StCAICE,                        &
     &                    IVTOffHr, ldebug,                             &
     &                    l_leads,Int_Head_topmelt,                     &
     &                    Real_Head_topmelt, ncols, nrowst,             &
     &                    topmelt,                                      &
#include "argppx.h"
     &                    icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 12. unable to read topmelt rate'
          icode = 1022
          go to 9999
        end if

!----------------------------------------------------------------------
! 13. Write out Topmelt rate
!----------------------------------------------------------------------
! DEPENDS ON: write_one_field
        call write_one_field (                                          &
#include "afields.h"
     &       OutStCTOP, FFTOP, PPTOP, IVTOffHr,                         &
     &       ITGrid, nrowst,                                            &
     &       Int_Head_topmelt, Real_Head_topmelt, IOutUnit, ldebug,     &
     &       topmelt, icode)
        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 13. unable to write topmelt rate'
          icode = 1112
          go to 9999
        end if

!----------------------------------------------------------------------
! 14. Read in Botmelt Rate
!----------------------------------------------------------------------
        l_leads = .false.                ! set to use minicefrac
! DEPENDS ON: read_leads_flds
        call read_leads_flds(StCBotmelt,StCAICE,                        &
     &                    IVTOffHr, ldebug,                             &
     &                    l_leads,Int_Head_botmelt,                     &
     &                    Real_Head_botmelt, ncols, nrowst,             &
     &                    bottommelt,                                   &
#include "argppx.h"
     &                    icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 14. unable to read botmelt rate'
          icode = 1023
          go to 9999
        end if

!----------------------------------------------------------------------
! 15. Write out Botmelt rate
!----------------------------------------------------------------------
! DEPENDS ON: write_one_field
        call write_one_field (                                          &
#include "afields.h"
     &       OutStCBOT, FFBOT, PPBOT, IVTOffHr,                         &
     &       ITGrid, nrowst,                                            &
     &       Int_Head_botmelt, Real_Head_botmelt, IOutUnit, ldebug,     &
     &       bottommelt, icode)
        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 15. unable to write botmelt rate'
          icode = 1113
          go to 9999
        end if

!----------------------------------------------------------------------
! 16. end loop over validity times
!----------------------------------------------------------------------
        enddo    !  ivt

9999  continue
      return
      END SUBROUTINE sea_ice
!----------------------------------------------------------------------
#endif
