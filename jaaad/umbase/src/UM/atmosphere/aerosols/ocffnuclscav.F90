#if defined(A17_2A) || defined(A17_2B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Perform nucleation scavenging of aged ocff to cloud ocff
!
!
      SUBROUTINE OCFFNUCLSCAV(                                          &
! Arguments IN
       rows, row_length, off_x, off_y, halo_i, halo_j,                  &
       model_levels, wet_model_levels, timestep,                        &
       cloudf, qcl, qcf,                                                &
! Arguments INOUT
       ocff_agd, ocff_cld,                                              &
! Arguments OUT
       delta_ocffnuclscav                                               &
       )
!
! Purpose:
!   To perform nccleation scavenging of aged ocff to form
!   the third mode of ocff, ocff in cloud water.
!
!   Called by Aero_ctl
!
! Current owner of code: Nicolas Bellouin
!
! Code Description:
!  Language: Fortran 90
!  This code is written to UMDP3 v8 programming standards
!
! Documentation: UMDP20
!
!
      IMPLICIT NONE
!
#include "c_ocff_chm.h"
!  includes parameters for rate of ocff nucleation/evaporation
!
#include "c_r_cp.h"
!  includes gas constant, R
!
#include "c_pi.h"
!
#include "c_g.h"
!
#include "c_densty.h"
!  includes density of water
!
! Arguments with intent IN:
!
      integer ::                                                        &
        rows,                                                           &
                          !no. of rows
        row_length,                                                     &
                          !no. of pts along a row
        off_x,                                                          &
                          !size of small halo in i
        off_y,                                                          &
                          !size of small halo in j
        halo_i,                                                         &
                          !EW halo size
        halo_j,                                                         &
                          !NS halo size
        model_levels,                                                   &
                          !no. of model levels
        wet_model_levels  !no. of wet model levels
!
      real ::                                                           &
        cloudf(row_length,rows,wet_model_levels),                       &
                                                 !Decimal cloud fraction
        qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
                              wet_model_levels),                        &
                                                 !Cloud liquid water
        qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
                              wet_model_levels),                        &
                                                 !Cloud frozen water
        timestep                                 !Atmos model timestep
!
! Arguments with intent IN:
!
      real ::                                                           &
         ocff_agd(1-off_x:row_length+off_x,1-off_y:rows+off_y,          &
                 model_levels),                                         &
                                  !mmr of aged OCFF
         ocff_cld(1-off_x:row_length+off_x,1-off_y:rows+off_y,          &
                 model_levels)    !mmr of OCFF-in-cloud
!
! Arguments with intent OUT:
!
      real ::                                                           &
        delta_ocffnuclscav(row_length,rows,model_levels)
                               ! cloud ocff increment due to
                               ! nucleation scavenging
!
! Local variables:
!
      integer :: i,j,k  ! loop counters
!
      real ::                                                           &
       delta_nuc(row_length,rows,model_levels),                         &
!              Increment to cloud OCFF
       delta_evap(row_length,rows,model_levels),                        &
!              Increment to aged OCFF
       qctotal(row_length,rows,wet_model_levels),                       &
!              Total cloud water
       clear_frac
!              Clear fraction of grid box (1.0 - cloudf)
!
!
      real ::                                                           &
        evaptime,                                                       &
                         ! timescale for cloud droplets to evaporate
        nuctime          ! timescale for particles to enter a cloud
!                          and nucleate.
!-----------------------------------------------------------------------
! 1. Initialise increments to zero
!-----------------------------------------------------------------------
      Do k=1,model_levels
        Do j=1,rows
          Do i=1,row_length
            delta_nuc(i,j,k)=0.0
            delta_evap(i,j,k)=0.0
          End Do
        End Do
      End Do
!
!-----------------------------------------------------------------------
! 2. Release ocff from evaporating cloud droplets in partly cloudy grid
!    boxes. Also release any ocff-in-cloud in cloud-free grid boxes.
!-----------------------------------------------------------------------
      Do k=1,wet_model_levels
        Do j=1,rows
          Do i=1,row_length
            qctotal(i,j,k)=qcl(i,j,k) + qcf(i,j,k)
            If (qctotal(i,j,k) .lt. thold) then
              delta_evap(i,j,k) = ocff_cld(i,j,k)
!                      evaporate all the cloud ocff in this grid box
            else if (cloudf(i,j,k) .lt. 0.95) then
              evaptime=evap_tau + 0.5*cloud_tau
              delta_evap(i,j,k) = (1.0 - exp(-timestep/evaptime)) *     &
                                  ocff_cld(i,j,k)
            else
              delta_evap(i,j,k) = 0.0
            end if
          End Do
        End Do
      End Do
!
!  Also evaporate any ocff-in-cloud on non-wet model levels
!
      If (wet_model_levels .lt. model_levels) then
        Do k=wet_model_levels+1,model_levels
          Do j=1,rows
            Do i=1,row_length
              delta_evap(i,j,k) = ocff_cld(i,j,k)
            End Do
          End Do
        End Do
      End If
!
!-----------------------------------------------------------------------
! 3. In-cloud scavenging of aged ocff particles by cloud droplets.
!    It is assumed that the ocff particles can nucleate nucleat 
!    droplets. 
!
!-----------------------------------------------------------------------
!
!
      Do k=1,wet_model_levels
        Do j=1,rows
          Do i=1,row_length

            clear_frac = 1.0 - cloudf(i,j,k)

            If ((qctotal(i,j,k) .ge. thold) .and.                       &
                               (cloudf(i,j,k) .gt.0.0)) then

!
      NUCTIME=NUC_TAU+((CLEAR_FRAC*CLOUD_TAU)/(2.0*CLOUDF(I,J,K)))
      DELTA_NUC(I,J,K)=(1.0-EXP(-TIMESTEP/NUCTIME))*OCFF_AGD(I,J,K)
!
            End If ! Test on qctotal and thold
          Enddo
        Enddo
      Enddo
!
!-----------------------------------------------------------------------
! 4. Calculate total increment for output
!-----------------------------------------------------------------------
!
      Do k=1,model_levels
        Do j=1,rows
          Do i=1,row_length
            delta_ocffnuclscav(i,j,k) =                                 &
                          delta_nuc(i,j,k) - delta_evap(i,j,k)
          End Do
        End Do
      End Do
!
      Return
      End subroutine OCFFNUCLSCAV
#endif
