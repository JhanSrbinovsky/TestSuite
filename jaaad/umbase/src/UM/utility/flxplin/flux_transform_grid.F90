#if defined(FLXPLIN)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! contains routines: Flux_Transform_Grid
!
! Purpose: Flux processing routine.
!          Controls main processing for Flux_Transform_Main
!
!    Model            Modification history:
!   version  Date
!    5.3  15/10/01  New deck. A. Hines
!    6.1  30/07/04  Update to input file units. A. Hines.
!
!    Programming standard :
!
!    Logical components covered :
!
!    System task:
!
!    External documentation:
!----------------------------------------------------------------------
      subroutine Flux_Transform_Grid (                                  &
#include "aflddims.h"
     &     ppxRecs,NoInFiles,icode)

      implicit none

! parameters used
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "clookadd.h"
#include "cfdcodes.h"
#include "plookups.h"

! declaration of argument list
#include "cflddims.h"
      integer NoInFiles ! IN number of input files
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

! local fields
      integer Int_Head(Len_IntHd)  ! integer header
      real Real_Head(Len_RealHd)   ! real header

      integer Int_Head_y(Len_IntHd) ! headers for y component
      real Real_Head_y(Len_RealHd)   ! of wind type fields

      real scalar(ncols, nrowst) ! scalar field

      real windx(ncols, nrowsuv)  ! wind field - x component
      real windy(ncols, nrowsuv)  ! wind field - x component
      real windx_tmp(ncols, nrowsuv) ! partially rotated
      real windy_tmp(ncols, nrowsuv) ! wind fields

! local scalars
      integer IROW_NUMBER
      character*80 cmessage
      integer IStC   ! stash code
      integer i, j, input_field, ifile   ! loop counters
      logical ldebug
      logical L_more_fields

! declaration of externals
      external GETPPX, read_lsms, field_interpolate_write
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'Flux_Transform_Grid'  ! subroutine name for err messages

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

! 2. Loop over fields to be read and written

! need to avoid unit 22 which is assigned to STASHMaster
      do ifile=IUnOutLow+10, IUnOutLow+NoInFiles+9

       print*,'processing file number ',ifile
       input_field = 0
       L_more_fields = .true.

       do while ( L_more_fields )
         input_field = input_field + 1

! 2.1 read the next field header

        read (ifile, IOStat = icode) Int_Head, Real_Head


        if ( icode  <   0 ) then
          write(UnWarn,*)CWarn,CSub,                                    &
     &       ' step 2.1 End of file ',ifile,' reached after',           &
     &       input_field-1,' fields'
          icode=0
          L_more_fields = .false.
        end if


        if ( L_more_fields ) then
          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &         ' step 2.1 Error reading ', input_field,'th field'
            go to 9999
          end if


          IStC = Int_Head(ITEM_CODE)

          print*,'STASH code for field ',input_field, ' is ',IStC
          print*,'scalar dimensions are ', ncols, nrowst

! 2.2 IF it is not a wind type field THEN

        if ( IStC  /=  StCWindSpeedU  .and. IStC  /=  OutStCTAUX        &
     & .and. IStC  /=  OutStCTAUY .and. IStC  /=  StCWindSpeedV )       &
     &  then

! 2.2.1 read in the scalar field

          read (ifile, IOStat = icode) scalar

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &         ' step 2.2.1 Error reading in ',                         &
     &         input_field, 'th field '
            go to 9999
          end if

! 2.2.2 call routine to interpolate to new grid and write it out

! DEPENDS ON: field_interpolate_write
          call field_interpolate_write(                                 &
#include "afields.h"
     &        Int_Head, Real_Head, ldebug, ITGrid, nrowst,              &
     &        IUnOutLow+NoInFiles+10, scalar, icode)

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' Step 2.2.2 error outputting scalar field '
            go to 9999
          end if

! 2.3  ELSE

        else ! not a scalar field

! 2.3.1 read in x component wind field

          read (ifile, IOStat = icode) windx

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &         ' step 2.3.1 Error reading in ',                         &
     &         input_field, 'th field - which is a wind field '
            go to 9999
          end if

! 2.3.2 read in the next field header

          read (ifile, IOStat = icode) Int_Head_y, Real_Head_y

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' step 2.3.2 Error reading header of ',                    &
     &       input_field, 'th field header '
            go to 9999
          end if

! 2.3.3 check that it is y component of a velocity type field

          IStC=Int_Head_y(ITEM_CODE)
          if (      IStC  /=  OutStCTAUY                                &
     &        .and. IStC  /=  StCWindSpeedV ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' step 2.3.3 Header must be y component wind field ',      &
     &       'STASH code is ',IStC
            go to 9999
          end if

! 2.3.4 read in the y component wind field

          read (ifile, IOStat = icode) windy

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' step 2.3.4 Error reading in ',                           &
     &       input_field, 'th field - which is a wind field '
            go to 9999
          end if

! 2.3.5 rotate wind components if necessary
          if (rotg) then
! DEPENDS ON: w_eqtoll
           call w_eqtoll(coef_angle1, coef_angle2, windx,               &
     &           windy, windx_tmp, windy_tmp, ncols*nrowsuv,            &
     &           ncols*nrowsuv)
          else
           do j = 1, nrowsuv
             do i = 1, ncols
               windx_tmp(i,j)=windx(i,j)
               windy_tmp(i,j)=windy(i,j)
             enddo
           enddo
          endif

          if (rotgO) then
! DEPENDS ON: w_lltoeq
           call w_lltoeq(coef_angle3, coef_angle4, windx_tmp,           &
     &           windy_tmp, windx, windy, ncols*nrowsuv,                &
     &           ncols*nrowsuv)
          else
           do j = 1, nrowsuv
             do i = 1, ncols
               windx(i,j)=windx_tmp(i,j)
               windy(i,j)=windy_tmp(i,j)
             enddo
           enddo
          endif

! 2.3.6 Output both fields

! DEPENDS ON: field_interpolate_write
          call field_interpolate_write(                                 &
#include "afields.h"
     &        Int_Head, Real_Head, ldebug, IUGrid, nrowsuv,             &
     &        IUnOutLow+NoInFiles+10, windx, icode)

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' Step 2.3.6 error outputting wind-x field '
            go to 9999
          end if

! DEPENDS ON: field_interpolate_write
          call field_interpolate_write(                                 &
#include "afields.h"
     &        Int_Head_y, Real_Head_y, ldebug, IUGrid, nrowsuv,         &
     &        IUnOutLow+NoInFiles+10, windy, icode)

          if ( icode  >   0 ) then
            write(UnErr,*)CErr,CSub,                                    &
     &       ' Step 2.3.6 error outputting wind-y field '
            go to 9999
          end if

        end if !   IStC  /=  wind codes

        end if ! L_more_fields

       end do
       print*,'finished processing file ',ifile

      end do

9999  continue
      return
      END SUBROUTINE Flux_Transform_Grid
!----------------------------------------------------------------------
#endif
