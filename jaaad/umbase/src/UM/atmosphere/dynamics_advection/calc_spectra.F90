#if defined(A12_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! subroutine calc_spectra
      SUBROUTINE CALC_SPECTRA(local_field, local_norm,                  &
     &                local_row_length, local_rows,                     &
     &                global_row_length, global_rows,                   &
     &                fld_type,halo_type,                               &
     &                gath_proc)

      IMPLICIT NONE
!
! Description:
!   Calculates the square of the norm of the complex Fourier
!   transform.
!
! Method:
!   N.B. Subroutine FFT_2D is currently a dummy routine which returns
!   input field. This will be replaced with a suitable FFT routine.
!
!   Gathers the field and uses routine FFT_2D to obtain the
!   Fourier transforms. The required square of the norm is then
!   scattered over all processors.  Post processing is required
!   to obtain the power spectra
!   i.e If |F(n)|^2 is the square of the norm (where n is the
!   frequency) then the discrete spectral intensity (or energy),
!   E(n), is
!   For odd number of data points
!   E(n)=2.|F(n)|^2  (n=1,Nyquist frequency)
!   For even number of data points
!   E(n)=2.|F(n)|^2  (n=1,Nyquist frequency-1)
!   E(n)=F(n)|^2      n=Nyquist frequency
!
!   Documentation available.
!
! Original Programmer: Carol Roadnight
! Current Code Owner: Carol Roadnight
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.5    26/02/03   This deck introduced            Carol Roadnight
!
!+ Dates should have leading zeroes in dd/mm/yy format
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Global variables (#include statements etc):
#include "stparam.h"
#include "parvars.h"

      Integer                                                           &
     & local_rows                                                       &
     &,local_row_length                                                 &
     &,global_rows                                                      &
     &,global_row_length                                                &
     &,fld_type                                                         &
     &,halo_type                                                        &
     &,gath_proc

      Real                                                              &
     & local_field(local_row_length*local_rows)                         &
                                                !IN - original 2D field
     &,local_norm(local_row_length*local_rows)  !OUT - square of norm

! Local variables
      Integer                                                           &
     & icount                                                           &
     &,minval                                                           &
     &,ifail                                                            &
     &,k2                                                               &
     &,ICODE                                                            &
                          !OUT   RETURN CODE FROM ROUTINE
     &,klev                                                             &
     &,i,j,k

      CHARACTER*256  CMESSAGE
      parameter(minval=1e-12)

      Real                                                              &
     & SPECTRA(global_row_length,global_rows)                           &
     &,SPECTRA_IM(global_row_length,global_rows)                        &
     &,SPECTRA2(global_row_length*global_rows)                          &
     &,SPECTRA_IM2(global_row_length*global_rows)                       &
     &,TRIGN(2*GLOBAL_ROW_LENGTH)                                       &
                                     !  HOLDS TRIGONOMETRIC TERMS
     &,TRIGM(2*GLOBAL_ROWS)                                             & 
                                     !  USED IN FFT'S
     &,WORK(2*GLOBAL_ROW_LENGTH*GLOBAL_ROWS)                            &
     &,spec_energy(GLOBAL_ROW_LENGTH/2+1,global_rows/2+1)               &
     &,sq_of_norm(global_row_length*global_rows)                        &
     &,field(global_row_length,global_rows)

!---------------------------------------------------------------------
! Section 1.  Gather 2D field and prepare field for C06FUE
!---------------------------------------------------------------------

! DEPENDS ON: gather_field
      CALL GATHER_FIELD(local_field,spectra,                            &
     &                local_row_length,local_rows,                      &
     &                global_row_length,global_rows,                    &
     &                fld_type,halo_type,                               &
     &                gath_proc,gc_all_proc_group,                      &
     &                ICODE,CMESSAGE)

      if (mype  ==  gath_proc) then

        icount=1
        do k2 = 1,global_rows
          do k=1,global_row_length
            field(k,k2)=spectra(k,k2)
            SPECTRA_IM(k,k2)=0.0
            if (abs(spectra(k,k2))  <   minval) then
              spectra(k,k2) = 0.0
            endif
            SPECTRA2(icount)=spectra(k,k2)
            SPECTRA_IM2(icount)=spectra_im(k,k2)
            icount=icount+1
          enddo
        enddo

!---------------------------------------------------------------------
! Section 2. Call routine to calculate 2D Fourier transform
!---------------------------------------------------------------------
        IFAIL=0
! DEPENDS ON: fft_2d
        CALL FFT_2D(global_rows,global_row_length,spectra2              &
     &              ,spectra_im2,'I',TRIGM,TRIGN,WORK,IFAIL)

        do k =1,global_row_length*global_rows
!          sq_of_norm(k)=spectra2(k)**2.0+
!     &               spectra_im2(k)**2.0
          sq_of_norm(k)=spectra2(k)
        enddo

      endif ! on processor gath_pe

!---------------------------------------------------------------------
! Section 3. Scatter the square of the norm
!---------------------------------------------------------------------
! DEPENDS ON: scatter_field
      CALL SCATTER_FIELD(local_norm,sq_of_norm,                         &
     &                local_row_length,local_rows,                      &
     &                global_row_length,global_rows,                    &
     &                fld_type,halo_type,                               &
     &                gath_proc,gc_all_proc_group,                      &
     &                ICODE,CMESSAGE)

!!    END OF ROUTINE CALC_SPECTRA
      RETURN
      END SUBROUTINE CALC_SPECTRA
#endif
