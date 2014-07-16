

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Subroutine DERV_LAND_FIELD : Computes no of land points in mpp jobs
!
! Subroutine Interface :
!
      SUBROUTINE CALC_LAND_FIELD (unit_no,fixhd,                        &
     &                            len1_lookup,len2_lookup,              &
     &                            icode,cmessage)

      implicit none

!     Arguments
      integer unit_no       ! IN Unit Number
      integer fixhd(256)    ! IN Fixed header
      integer len1_lookup   ! IN First dimension of lookup table
      integer len2_lookup   ! IN Seconf dimension of lookup table
      integer icode         ! OUT Return code

      character*(*) cmessage ! OUT Error message

!     Local variables
      integer len_io        ! length of data returned from buffin
      integer lookup(len1_lookup,len2_lookup)   !  Lookup table
      real rcode            ! Real return code
!
!     Position atmos dump to read in lookup-table
! DEPENDS ON: setpos
      call setpos (unit_no,fixhd(150)-1,icode)

!     Check error code from setpos
      if (icode >  0) then
        write (6,*) 'Error in SETPOS called from CALC_LAND_FIELD.'
        write (6,*) 'Trying to point to start of lookup table '//       &

     &              'in atmos dump.'

        write (cmessage,*) 'DRLANDF1 : Error in SETPOS.'
        go to 9999   !  Return
      endif

!     Read in the look-up table
! DEPENDS ON: buffin
      call buffin (unit_no,lookup,len1_lookup*len2_lookup,len_io,rcode)

!     Check error code from buffin
      if (rcode /= -1.0) then
        write (6,*) 'Error in BUFFIN called from CALC_LAND_FIELD.'

        write (6,*) 'Trying to read lookup table from atmos dump.'

        write (6,*) 'Return code from BUFFIN ',rcode
        ICODE = 100
        write (cmessage,*) 'DRLANDF1 : Error in BUFFIN.'
        go to 9999   !  Return
      endif

!     Read in land-sea mask and then
!     compute the number of land points for each PE
! DEPENDS ON: read_land_sea
      CALL READ_LAND_SEA (unit_no,rcode,lookup,len1_lookup,len2_lookup, &
     &                    fixhd,256)

!     Check error code from read_land_sea
      if (rcode /= -1.0) then
        write (6,*) 'Error in READ_LAND_SEA.'
        write (6,*) 'Return code from READ_LAND_SEA ',rcode
        ICODE = 200
        write (cmessage,*) 'DRLANDF1 : Error in READ_LAND_SEA.'
        go to 9999   !  Return
      endif

 9999 continue
      RETURN
      END SUBROUTINE CALC_LAND_FIELD
