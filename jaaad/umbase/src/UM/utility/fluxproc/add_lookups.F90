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
! Author:     M. J. Bell
!----------------------------------------------------------------------
! contains routines: add_lookups
!
! Purpose: Flux processing routine.
!          This routine modifies the lookup table for a file (copy in
!          memory only) to add further headers that are appropriate to
!          additional times. This is needed because the input fields
!          may be for a single time (eg SST) but the fields are needed
!          throughout the period of a FOAM forecast. By adding further
!          headers to the lookup table in memory, it is possible to
!          allow direct access to the fields on disk using the same
!          method whether the required time is actually in the file
!          or is obtained by re-labelling the field.
!
! Important: if readff used element LBEGIN in the lookup table to
!            determine the start location of a data field, the array
!            LookFldNo would not be necessary. But as things stand the
!            original lookup table has to be used to read in the data.
!            As there is no spare space in the Lookup table I have had
!            to use an additional array LookFldNo to store this field
!            number. The validity date of copied data read in by
!            read_one_field will thus be taken from the original lookup
!            table and will not be the validity date searched for !!
!                                          M. J. Bell 24/09/96
!----------------------------------------------------------------------
      subroutine add_lookups (                                          &
     &  NoAddTimes, ISrchOffHr, INewOffHr, l_climate_field,             &
     &  Len1_Lookup, Len2_Lookup,                                       &
     &  Len2_Lookup_Actual, Lookup, LookFldNo, icode )


      implicit none

! declaration of argument list
      integer NoAddTimes          ! IN # of additional lookups to make
      integer ISrchOffHr(NoAddTimes+1) ! IN offset hours to look for
      integer INewOffHr (NoAddTimes+1) ! IN new offset hours
      logical l_climate_field    ! IN F => not a climate field
      integer Len1_Lookup        ! IN length of lookup table (64)
      integer Len2_Lookup        ! IN length allocated for lookup table
      integer Len2_Lookup_Actual ! IN/OUT  actual # of lookup tables
      integer Lookup(Len1_Lookup, Len2_Lookup) ! IN/OUT lookups
      integer LookFldNo(Len2_Lookup) ! OUT field nos for lookups
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! declaration of parameters
#include "clookadd.h"
#include "cfdcodes.h"

! declaration of globals used
#include "cunitnos.h"
#include "cmess.h"
#include "creftim.h"
#include "cclm1tim.h"
#include "cclm2tim.h"

! declaration of local arrays

! declaration of local scalars
      integer Len2_Original        !
      integer i                    ! loop index over lookup tables
      integer iadd                 ! loop index over additional times
      integer inew                 ! index of new lookup table
      integer j                    ! loop index over elements in lookup
      integer StCode               ! stash code of field

      logical l_data_time  ! T => use data time; F => use validity time
      logical Itemfound    ! T => item found in lookup table

! externals used
      external add_hours, time_to_use
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'add_lookups'  ! subroutine name for error messages

! 1.0 Check that all Len2_Lookup_Actual lookup tables are non-empty;
!     Reset Len2_lookup_Actual and Len2_Original if an empty table is
!     found (empty lookup tables have -99 in all entries; entry 1 has
!     the year of validity in a non-empty table.)

      Len2_Original =  Len2_Lookup_Actual

      do i = 1, Len2_Original
        if ( Lookup (1, i)  ==  -99 ) then
          Len2_Lookup_Actual = i - 1
          go to 100
        end if
      end do   ! i

100   continue

      Len2_Original = Len2_Lookup_Actual

! 1.1 Set field numbers for original lookup tables

      do i = 1, Len2_Lookup_Actual
        LookFldNo(i) = i
      end do

! 2. Start Loop over validity times to copy headers and calculate
!    new validity dates

      do iadd = 1, NoAddTimes

! 3. Calculate old and new validity dates

! DEPENDS ON: add_hours
        call add_hours (                                                &
#include "areftim.h"
#include "aclm1tim.h"
     & ISrchOffHr(iadd) )

! DEPENDS ON: add_hours
        call add_hours (                                                &
#include "areftim.h"
#include "aclm2tim.h"
     & INewOffHr(iadd) )

! 4. Check that there are no original lookup tables present
!    which have the new validity date

        ItemFound = .false.
        do i = 1, Len2_Original

          StCode = Lookup(ITEM_CODE, i)
! DEPENDS ON: time_to_use
          call time_to_use ( StCode, l_climate_field, l_data_time)

          if ( l_data_time ) then

            if (Lookup(LBYRD,i)   ==  Clim2Year   .and.                 &
     &          Lookup(LBMOND,i)  ==  Clim2Month  .and.                 &
     &          Lookup(LBDATD,i)  ==  Clim2Day    .and.                 &
     &          Lookup(LBHRD,i)   ==  Clim2Hour   .and.                 &
     &          Lookup(LBMIND,i)  ==  Clim2Min          ) then

              ItemFound = .true.
            end if

          else ! .not. l_data_time

            if (Lookup(LBYR,i)   ==  Clim2Year   .and.                  &
     &          Lookup(LBMON,i)  ==  Clim2Month  .and.                  &
     &          Lookup(LBDAT,i)  ==  Clim2Day    .and.                  &
     &          Lookup(LBHR,i)   ==  Clim2Hour   .and.                  &
     &          Lookup(LBMIN,i)  ==  Clim2Min          ) then

              ItemFound = .true.
            end if

          end if ! l_data_time

        end do ! i

        if ( ItemFound ) then
          icode = 13
          write(UnErr,*)CErr,CSub,                                      &
     &    ' step 4. data already exists with new validity time '
          go to 9999
        end if

! 5. Start loop over original lookups to look for fields with
!    old validity date

        do i = 1, Len2_Original
          ItemFound = .false.

! 6. Look for a validity date

          StCode = Lookup(ITEM_CODE, i)
! DEPENDS ON: time_to_use
          call time_to_use ( StCode, l_climate_field, l_data_time)

! 6.1 Do not add lookup for accumulations fields
          if ( StCode  /=  StCDrain    .and.                            &
     &         StCode  /=  StCConvrain .and.                            &
     &         StCode  /=  StCDsnow    .and.                            &
     &         StCode  /=  StCConvsnow ) then

            if ( l_data_time ) then

              if (Lookup(LBYRD,i)   ==  Clim1Year   .and.               &
     &          Lookup(LBMOND,i)  ==  Clim1Month  .and.                 &
     &          Lookup(LBDATD,i)  ==  Clim1Day    .and.                 &
     &          Lookup(LBHRD,i)   ==  Clim1Hour   .and.                 &
     &          Lookup(LBMIND,i)  ==  Clim1Min          ) then

                ItemFound = .true.
              end if

            else ! .not. l_data_time

              if (Lookup(LBYR,i)   ==  Clim1Year   .and.                &
     &          Lookup(LBMON,i)  ==  Clim1Month  .and.                  &
     &          Lookup(LBDAT,i)  ==  Clim1Day    .and.                  &
     &          Lookup(LBHR,i)   ==  Clim1Hour   .and.                  &
     &          Lookup(LBMIN,i)  ==  Clim1Min          ) then

                ItemFound = .true.
              end if

            end if ! l_data_time

          end if  ! StCode
! 7. If found

          if ( ItemFound ) then

! 7.1  check space available for new lookup

             if ( Len2_Lookup_Actual + 1  >   Len2_Lookup) then
               icode = 14
               write(UnErr,*)CErr,CSub,                                 &
     &         ' step 4. not enough space to add an extra',             &
     &         ' field in lookups. Len2_Lookup = ', Len2_Lookup
               go to 9999
             end if

! 7.2 Copy the lookup table and amend dates in new lookup table

            inew = Len2_Lookup_Actual + 1

            do j = 1, Len1_Lookup
              Lookup(j,inew) = Lookup(j,i)
            end do


            if ( l_data_time ) then

              Lookup(LBYRD,inew)  = Clim2Year
              Lookup(LBMOND,inew) = Clim2Month
              Lookup(LBDATD,inew) = Clim2Day
              Lookup(LBHRD,inew)  = Clim2Hour
              Lookup(LBMIND,inew) = Clim2Min

            else ! .not. l_data_time

              Lookup(LBYR,inew)  = Clim2Year
              Lookup(LBMON,inew) = Clim2Month
              Lookup(LBDAT,inew) = Clim2Day
              Lookup(LBHR,inew)  = Clim2Hour
              Lookup(LBMIN,inew) = Clim2Min

            end if ! l_data_time

! 7.3 Amend number of non-empty lookup tables and set
!    pointer to lookup table to actually use

            Len2_Lookup_Actual = Len2_Lookup_Actual + 1
            LookFldNo(inew) = i

          end if ! ItemFound

! 8. End loop over  lookups (started in 5.)

         end do ! i = 1, Len2_Original

! 9. End loop over validity times (started in 2.)

       end do ! iadd = 1, NoAddTimes


9999   continue
       return
       END SUBROUTINE add_lookups
!----------------------------------------------------------------------
#endif
