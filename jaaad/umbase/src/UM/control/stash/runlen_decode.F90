#if defined(C84_1A) || defined(FLDIO) || defined(UTILIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! 5.2      07/09/99     New code
!
! Author:     I. Edmond
!----------------------------------------------------------------------
! contains routines:  runlen_encode, runlen_decode
!
! Purpose: Reduce storage requirements by removing land points from
!          fields within ocean fieldsfiles. This is done using run
!          length encoding to compress the sequences of missing data
!          values that represent the land points. Using this method,
!          a sequence of missing data values are mapped into a single
!          missing data value followed by a value indicating the
!          length of the run of missing data values.
!
!----------------------------------------------------------------------


      subroutine runlen_decode(unpacked,unpacked_size,packed,           &
     &                         packed_size,bmdi,icode,cmessage)

      integer i,j,k
      integer nmdi

      integer unpacked_size    ! IN  Size of input field used to hold
                               !     unpacked data
      integer packed_size      ! IN  Size of field used to hold packed
                               !     data.
      integer icode, ocode     ! IN/OUT Error codes >0 indicates error
                               !                    =0 success.

      real unpacked(unpacked_size) ! OUT Input field representing
                               !     unpacked data.
      real packed(packed_size) ! IN  output field representing run
                               !     length decoded fields.
      real bmdi                ! IN  Real missing data values

      character*(80) cmessage  ! OUT Returned error message if icode >0

      i    = 0
      k    = 1
      nmdi = 0

!    Check if an error has already been encountered, and get out
!    if it has.
      ocode = 0
      if (icode  >   0) then
         goto 999
      else if (icode  <   0) then
         ocode = icode
         icode = 0
      end if

      do while(k <= packed_size)

        if (packed(k)  ==  bmdi) then
          nmdi = packed(k+1)

          do j=1, nmdi
            i = i + 1
            unpacked(i) = bmdi

            if (i  >   unpacked_size) then
              write (6,*)                                               &
     &          'RLENCODE : ERROR : Ocean Run length decoding failed'
              cmessage = 'Ocean Run length decoding failed'
              icode = 1
              goto 999
            end if
          end do

          k = k + 2
        else
          i = i + 1
          unpacked(i) = packed(k)
          k = k + 1

          if (i  >   unpacked_size) then
            write (6,*)                                                 &
     &        'RLENCODE : ERROR : Ocean Run length decoding failed'
            cmessage = 'Ocean Run length decoding failed'
            icode = 1
            goto 999
          end if

        end if

      end do

 999  continue

!     If we have found an error, leave it in icode.  If no error
!     occurred then check if the original input value of icode was]
!     non-zero (a previous untrapped error/warning), and copy this
!     back into ICODE before eaving the routine.
      IF (icode  ==  0 .and. ocode  /=  0) then
         icode = ocode
      END IF

      return
      END SUBROUTINE runlen_decode


!----------------------------------------------------------------------
#endif
