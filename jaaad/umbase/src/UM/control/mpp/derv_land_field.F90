#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C)
#if defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Subroutine DERV_LAND_FIELD : Computes no of land points in MPP jobs
!
! Subroutine Interface :
!
      SUBROUTINE DERV_LAND_FIELD (unit_no,icode,cmessage)

      implicit none
!
! Description : Calculates the no of land points on each PE.
!
! Method : Call READ_LAND_SEA to read in Land-Sea Mask from
!          Atmosphere Dump and then calculate no of land points.
!
! Current Code Owner : Dave Robinson, NWP
!
! History :
! Version    Date    Comment
! -------    ----    -------
!   4.5    15/04/98  Original Code
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!LL  5.4  03/09/02  Remove DEFS MPP.              E.Leung
!    5.5  01/08/00  Modification for parallelisation of WAM - read
!                   the land/sea mask for a Wave Job
!                   Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes
!  6.2  23/11/05  Removed all references to the wavemodel.
!                 T.Edwards
!
! Code Description :
! Language : FORTRAN 77 + common extensions
!
! Declarations :

!     Arguments
      integer unit_no        ! IN  Unit number for atmos dump
      integer icode          ! OUT Error Code
      character*(*) cmessage ! OUT Error message

!     Local variables
      integer ilen1_lookup   ! First dimesion of look-up table
      integer ilen2_lookup   ! Second dimension of look-up table
      integer fixhd(256)     ! Fixed header

#include "cenvir.h"
#include "parparm.h"
#include "typsize.h"

#if defined(ATMOS)
!     land_field is the global no of land-points.

!     Initialise global_land_field
      global_land_field = land_field

      write (6,*) ' global_land_field set to ',land_field
#endif

!     Open atmos/wave input dump
! DEPENDS ON: file_open
      call file_open (unit_no,ft_environ(unit_no),                      &
     &                len_ft_envir(unit_no),0,0,icode)

!     Check error code from file_open
      if (icode >  0) then
        write (6,*) 'Error in FILE_OPEN called from DERV_LAND_FIELD.'
        write (6,*) 'Trying to open atmos dump.'
        write (cmessage,*) 'DRLANDF1 : Error in FILE_OPEN.'
        go to 9999   !  Return
      endif

       if ( land_field == 0 ) print*,'land_field1',land_field
!     Read fixed header
! DEPENDS ON: read_flh
      call read_flh(unit_no,fixhd,256,icode,cmessage)

       if ( land_field == 0 ) print*,'land_field2',land_field
!     Check error code from read_flh
      if (icode >  0) then
        write (6,*) 'Error in READ_FLH called from DERV_LAND_FIELD.'
        write (6,*) 'Trying to read fixed header from atmos dump.'
        go to 9999   !  Return
      endif

!     Get dimensions of look-up table
      ilen1_lookup=fixhd(151)
      ilen2_lookup=fixhd(152)

!     Proceed to calculate no of land points on each PE.
! DEPENDS ON: calc_land_field
      CALL CALC_LAND_FIELD (unit_no,fixhd,ilen1_lookup,ilen2_lookup,    &
     &                      icode,cmessage)

       if ( land_field == 0 ) print*,'land_field3',land_field
#if defined(ATMOS)
!     land_field now contains the no of land_points for this PE.

!     Initialise local_land_field
      local_land_field = land_field

      write (6,*) ' local_land_field set to ',land_field
#endif

 9999 continue

      RETURN
      END SUBROUTINE DERV_LAND_FIELD
#endif
#endif
