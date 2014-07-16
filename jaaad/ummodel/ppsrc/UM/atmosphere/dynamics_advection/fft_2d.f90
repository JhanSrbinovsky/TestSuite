
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! subroutine fft_2d
      SUBROUTINE FFT_2D(global_rows,global_row_length,spectra2          &
     &              ,spectra_im2,INIT,TRIGM,TRIGN,WORK,IFAIL)

      IMPLICIT NONE
!
! Description:
!   Dummy subroutine.  To be replaced by a subroutine which calculates
!   2D FFTs.
!
! Method:
!   Passes the original field back without performing any calculations
!   on it.
!
!
! Original Programmer: Carol Roadnight
! Current Code Owner: Carol Roadnight
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.5    26/02/03  This deck introduced            Carol Roadnight
!
!+ Dates should have leading zeroes in dd/mm/yy format
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Global variables (#include statements etc):

      Integer                                                           &
     & global_rows                                                      &
     &,global_row_length                                                &
     &,ifail

      Real                                                              &
     & SPECTRA2(global_row_length*global_rows)                          &
     &,SPECTRA_IM2(global_row_length*global_rows)                       &
     &,TRIGN(2*GLOBAL_ROW_LENGTH)                                       &
                                     !  HOLDS TRIGONOMETRIC TERMS
     &,TRIGM(2*GLOBAL_ROWS)                                             & 
                                     !  USED IN FFT'S
     &,WORK(2*GLOBAL_ROW_LENGTH*GLOBAL_ROWS)

      CHARACTER*1  INIT

!---------------------------------------------------------------------
! Section 1.
!---------------------------------------------------------------------

      write(6,*)'**************WARNING*****************'
      write(6,*)'******FFT_2D is a dummy routine*******'
      write(6,*)'*****original field is passed back****'
      write(6,*)'**************************************'

!!    END OF ROUTINE FFT_2D
      RETURN
      END SUBROUTINE FFT_2D
