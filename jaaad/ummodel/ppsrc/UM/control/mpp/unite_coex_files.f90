

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Unifies many packed fragments of a data field into one

      SUBROUTINE unite_coex_files(input_buffer,                         &
     &                           output_buffer, num_out,                &
     &                           row_start_pe)

      IMPLICIT NONE
!
! Description:
!   Combines blocks of WGDOS packed data into a unified packed
!   data field.
!
! Method:
!   If the row block is the first one, the input buffer is just
!   copied into the output buffer. If it is a later one, the
!   input buffer is copied into the output buffer from word 4
!   (32-bit) onwards. After this, the coex header is updated
!   with the new sizes (see UMDP F3 for coex header/data structure)
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.1   13/07/04   Original code. P.Selwood/B.Carutthers(Cray)
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:

      Integer, Intent(In)   ::                                          &
     &  input_buffer(*)      ! Buffer holding the current
                             ! block of compressed data, which
                             ! is to be added to the o/p buffer


      Integer, Intent(Out)   ::                                         &
     &  output_buffer(*)     ! The output buffer that holds all
                             ! the collected, compressed data
!
      Integer, Intent(In) ::                                            &
     &  row_start_pe         ! First PE in a block - may or may not
                             ! be the gather PE

      Integer, Intent(Out) ::                                           &
     &  num_out              ! Length of the o/p data buffer

! Fuctions called
      integer,external :: ibm2ieee,ieee2ibm
!
! Work variables:
!
      Integer :: num_in,                                                &
     &           iy_in,                                                 &
     &           iy_out,                                                &
     &           ierr

!- End of header

! Find the total length in the o/p buffer thus far
!
      ierr=ibm2ieee(2,1,output_buffer(1),0,num_out,1,64,32)
!
! Find the length of the additional entry
!
      ierr=ibm2ieee(2,1,input_buffer(1),0,num_in,1,64,32)

!
! If this is the global_gather_pe, we need all the
! control information from the header
!
      if(row_start_pe == 0) then

        num_out=0
! DEPENDS ON: copy_buffer_32
        call copy_buffer_32(input_buffer, output_buffer, num_in, 0, 0)
        num_out = num_in

      else
!
! New data to add starts in word 4 (32-bit) of input_buffer
!
! DEPENDS ON: copy_buffer_32
         call copy_buffer_32(input_buffer, output_buffer, num_in,       &
     &                       num_out, 3)

         num_out = num_out + (num_in - 3)
!
! Now update the number of rows
!
!
        ierr=ibm2ieee(2,1,output_buffer(2),16,iy_out,1,64,16)
        ierr=ibm2ieee(2,1,input_buffer(2),16,iy_in,1,64,16)
        iy_out=iy_out+iy_in
!
        ierr=ieee2ibm(2,1,output_buffer(2),16,iy_out,1,64,16)

      endif
!
! Now update the number of 32-bit words
!
      ierr=ieee2ibm(2,1,output_buffer(1),0,num_out,1,64,32)
!
      return
      END SUBROUTINE unite_coex_files
