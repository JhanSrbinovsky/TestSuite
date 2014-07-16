#if defined(A17_2A) || defined(A17_2B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Convert a proportion of fresh soot to aged soot
!
!
!-----------------------------------------------------------------------
!
!+ Perform diffusional scavenging of aged soot to cloud soot
!
      SUBROUTINE SOOTDIFFSCAV(                                          &
! Arguments IN
     & rows, row_length, off_x, off_y, halo_i, halo_j,                  &
     & model_levels, wet_model_levels, timestep,                        &
     & cloudf, qcl, qcf, p, t,                                          &
     & n_droplet,                                                       &
! Arguments INOUT
     & soot_agd, soot_cld,                                              &
! Arguments OUT
     & delta_sootdiffscav                                               &
     & )
!
! Purpose:
!   To perform diffusional scavenging of aged soot to form
!   the third mode of soot, soot in cloud water.
!
!   Called by Aero_ctl
!
! Current owners of code: P. Davison
!
! History:
! Version     Date     Comment
! -------     ----     -------
!
!   5.4      10/07/02  Original code, based on a re-write of
!                      v4.5 deck SOOT1A. Diffusional scavenging
!                      now included.                  P. Davison
!   5.5      07/02/03  Number_droplet calculation moved up to aero_ctl.
!                                                      P. Davison
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
! Documentation: UMDP20
!
!
      IMPLICIT NONE
!
#include "c_st_chm.h"
!  includes parameters for rate of soot nucleation/evaporation
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
      INTEGER                                                           &
     &  rows                                                            &
                          !no. of rows
     &, row_length                                                      &
                          !no. of pts along a row
     &, off_x                                                           &
                          !size of small halo in i
     &, off_y                                                           &
                          !size of small halo in j
     &, halo_i                                                          &
                          !EW halo size
     &, halo_j                                                          &
                          !NS halo size
     &, model_levels                                                    &
                          !no. of model levels
     &, wet_model_levels  !no. of wet model levels
!
      REAL                                                              &
     &  cloudf(row_length,rows,wet_model_levels)                        &
                                                 !Decimal cloud fraction
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &                        wet_model_levels)                         &
                                                 !Cloud liquid water
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &                        wet_model_levels)                         &
                                                 !Cloud frozen water
     &, p(1-off_x:row_length+off_x,1-off_y:rows+off_y,                  &
     &          model_levels)                                           &
                                                 !pressure
     &, t(row_length,rows,model_levels)                                 &
                                                 !temperature
     &, timestep                                                        &
                                                 !Atmos model timestep
     &, n_droplet(row_length,rows,wet_model_levels)
                                                 !droplet concentration
!
! Arguments with intent IN:
!
      REAL                                                              &
     &  soot_agd(1-off_x:row_length+off_x,1-off_y:rows+off_y,           &
     &          model_levels)                                           &
                                 !mmr of aged soot
     &, soot_cld(1-off_x:row_length+off_x,1-off_y:rows+off_y,           &
     &          model_levels)    !mmr of soot-in-cloud
!
! Arguments with intent OUT:
!
      REAL                                                              &
     &  delta_sootdiffscav(row_length,rows,model_levels)
                               ! cloud soot increment due to
                               ! diffusional scavenging
!
! Local variables:
!
      INTEGER                                                           &
     & i,j,k  ! loop counters
!
      REAL                                                              &
     & delta_nuc(row_length,rows,model_levels),                         &
!              Increment to cloud soot
     & delta_evap(row_length,rows,model_levels),                        &
!              Increment to aged soot
     & qctotal(row_length,rows,wet_model_levels),                       &
!              Total cloud water
     & clear_frac
!              Clear fraction of grid box (1.0 - cloudf)
!
!    Variables for improved diffusional scavenging.
!
      REAL                                                              &
     &  diffusivity                                                     &
                         ! Mean diffusion coefficent of aged soot.
     &, viscosity_air                                                   &
                         ! Dynamic viscosity of air (kg m-1 s-1).
     &, mean_free_path                                                  &
                         ! Mean free path of air molecules.
     &, mfp_ref                                                         &
                         ! Reference value of mean free path.
     &, tref_mfp                                                        &
                         ! Reference temperature for mean free path.
     &, pref_mfp                                                        &
                         ! Reference pressure for mean free path.
     &, Knudsen_Weber                                                   &
                         ! Expression related to the Cunningham slip
!                          flow correction to the diffusion coefficient.
     &, rad_age                                                         &
                         ! Median radius of aged soot distribution, (m).
     &, sigma_age                                                       &
                         ! Geometric standard deviation of the aged soot
!                          mode distribution.
     &, log_sigma_age                                                   &
                         ! Natural log of the previous parameter.
     &, sq_log_sigma_age                                                &
                         ! Square of the previous parameter.
     &, diff_con1                                                       &
                         ! Term in diffusion coefficent formula.
     &, diff_con2                                                       &
                         ! Term in diffusion coefficent formula.
     &, Boltzmann                                                       &
                         ! Boltzmanns constant.
     &, Pec                                                             &
                         ! Quantity associated with (not equal to)
!                          Peclet number.
     &, work_radius                                                     &
                         ! Variable linked to average droplet radius.
     &, scavcd           ! Scavenging coefficient for in-cloud
!                          advective-diffusive removal of aged
!                          soot particles.
!
      PARAMETER (Boltzmann = 1.3804E-23)  ! J K-1
      PARAMETER (mfp_ref = 6.6E-8                                       &
                                          ! m
     &          ,tref_mfp = 293.15                                      &
                                          ! K
     &          ,pref_mfp = 1.01325E5)    ! Pa
      PARAMETER (sigma_age = 2.0, rad_age = 4.0E-8)
!
      REAL                                                              &
     &  evaptime                                                        &
                         ! timescale for cloud droplets to evaporate
     &, nuctime                                                         &
                         ! timescale for particles to enter a cloud
!                          and nucleate.
     &, diffuse_tau                                                     &
                         ! diffusive lifetime of soot particles
!                          once they enter a cloud
     &, rho_cuberoot                                                    &
                         ! cube root of density of water.
     &, diffuse_tauratio                                                &
                         ! Cloudtau/Diffuse_tau
     &, probdiff_inv                                                    &
                         ! inverse of probability of a particle being
!                          diffusionallly scavenged in a cloud.
     &, probdiff_fn1                                                    &
                         ! Probdiff_inv - 0.5
     &, probdiff_fn2                                                    &
                         ! Probdiff_inv*EXP(diffuse_tauratio*0.5)
     &, probdiff_cloud                                                  &
                         ! probability of a soot particle
!                          being in cloud at the start of a step.
     &, probdiff_clear                                                  &
                         ! probability of a soot particle
!                          being in clear air at the start of a step.
     &, lambda_soot                                                     &
                         ! ratio of concentrations of soot
!                          particles in cloudy to clear air.
     &, diffrate         ! rate of diffusive capture of soot particles
!
      REAL                                                              &
     &  Term3                                                           &
                         ! Local workspace
     &, Term4                                                           &
                         ! Local workspace
     &, Denom            ! Local workspace
!
!     Extra parameters for improved diffusional scavenging.
!
      log_sigma_age = LOG(sigma_age)
      sq_log_sigma_age = log_sigma_age * log_sigma_age
      diff_con1 = EXP(-2.5*sq_log_sigma_age)/rad_age
      diff_con2 = EXP(-4.0*sq_log_sigma_age)/(rad_age*rad_age)
!
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
! 2. Release soot from evaporating cloud droplets in partly cloudy grid
!    boxes. Also release any soot-in-cloud in cloud-free grid boxes.
!-----------------------------------------------------------------------
      Do k=1,wet_model_levels
        Do j=1,rows
          Do i=1,row_length
            qctotal(i,j,k)=qcl(i,j,k) + qcf(i,j,k)
            If (qctotal(i,j,k)  <   thold) then
              delta_evap(i,j,k) = soot_cld(i,j,k)
!                      evaporate all the cloud soot in this grid box
            else if (cloudf(i,j,k)  <   0.95) then
              evaptime=evaptau + 0.5*cloudtau
              delta_evap(i,j,k) = (1.0 - exp(-timestep/evaptime)) *     &
     &                            soot_cld(i,j,k)
            else
              delta_evap(i,j,k) = 0.0
            end if
          End Do
        End Do
      End Do
!
!  Also evaporate any soot-in-cloud on non-wet model levels
!
      If (wet_model_levels  <   model_levels) then
        Do k=wet_model_levels+1,model_levels
          Do j=1,rows
            Do i=1,row_length
              delta_evap(i,j,k) = soot_cld(i,j,k)
            End Do
          End Do
        End Do
      End If
!
!-----------------------------------------------------------------------
! 3. In-cloud scavenging of aged soot particles by cloud droplets.
!    It is assumed that the soot particles do not themselves nucleate
!    droplets. The parametrisation is the same as that used for Aitken
!    mode sulphate aerosol.
!-----------------------------------------------------------------------
!
!
      Do k=1,wet_model_levels
        Do j=1,rows
          Do i=1,row_length

            clear_frac = 1.0 - cloudf(i,j,k)

            If ((qctotal(i,j,k)  >=  thold) .and.                       &
     &                         (cloudf(i,j,k)  >  0.0)) then

!     First compute in-cloud timescale for diffusional capture,
!     using total condensed water within the cloudy portion of the
!     grid box. The difference between liquid water and ice is
!     neglected here. This should be improved on eventually.
!
!     Compute mean free path of air molecules.
!     (See P.417 of Pruppacher and Klett, 2nd edition.)
!
              MEAN_FREE_PATH = (MFP_REF*PREF_MFP*T(i,j,k))/             &
     &                           (TREF_MFP*P(i,j,k))
!
!     Compute the Knudsen and Weber term for a particle of the median
!     radius (note approximation here: we do not average over the size
!     distribution). See P.450 of Pruppacher and Klett, 2nd edition.
!
              KNUDSEN_WEBER = 1.257 + 0.4*                              &
     &                          EXP(- ((1.1*RAD_AGE)/MEAN_FREE_PATH) )
!
!     Temporarily use DIFFUSIVITY to store working value.
!
              DIFFUSIVITY = DIFF_CON1 +                                 &
     &                        DIFF_CON2*MEAN_FREE_PATH*KNUDSEN_WEBER
!
!     Compute dynamic viscosity of air, using an approximate version of
!     the formula on P.417 of Pruppacher and Klett, 2nd edition.
!
              VISCOSITY_AIR = ( 1.718 + (T(i,j,k) - 273.15)*0.0049 )*   &
     &                            1.0E-5
!
!     Now compute mean diffusion coefficient.
!
              DIFFUSIVITY = (BOLTZMANN*T(i,j,k)*DIFFUSIVITY)/           &
     &                        (6.0*PI*VISCOSITY_AIR)
!
!     Now compute the term PEC related to (but not equal to) the cube
!     root of the Peclet Number.
!
              PEC= ((4.0*G*RHO_WATER)/(9.0*DIFFUSIVITY*VISCOSITY_AIR))  &
     &                                                       **0.333333
!
              WORK_RADIUS = QCTOTAL(i,j,k)/                             &
     &               (CLOUDF(i,j,k)*10.0*PI*RHO_WATER*N_DROPLET(i,j,k))
              WORK_RADIUS = WORK_RADIUS**0.333333
!
!     We can finally compute the timescale for diffusive
!     scavenging once inside cloud, DIFFUSE_TAU.
!
              SCAVCD = 6.0*PI*DIFFUSIVITY*N_DROPLET(i,j,k)*WORK_RADIUS* &
     &                   (1.0 + PEC*WORK_RADIUS)
              DIFFUSE_TAU = 1.0/SCAVCD
!
              DIFFUSE_TAURATIO = CLOUDTAU/DIFFUSE_TAU
              PROBDIFF_INV = 1.0/( 1.0 - EXP(-DIFFUSE_TAURATIO) )
              PROBDIFF_FN1 = PROBDIFF_INV - 0.5
              PROBDIFF_FN2 = PROBDIFF_INV*EXP(-(0.5*DIFFUSE_TAURATIO) )
!   CALCULATE LAMBDA_SOOT.
!
              TERM3 = (CLEAR_FRAC*DIFFUSE_TAURATIO)**2
              TERM3 = TERM3 + ( 2.0*DIFFUSE_TAURATIO                    &
     &                 *CLEAR_FRAC*(CLEAR_FRAC-CLOUDF(i,j,k)))
              TERM3 = SQRT(1.0+TERM3)
              TERM4=2.0*CLOUDF(i,j,k)-                                  &
     &                 DIFFUSE_TAURATIO*CLEAR_FRAC-1.0
              TERM4 = TERM4 + TERM3
              LAMBDA_SOOT = TERM4/(2.0*CLOUDF(i,j,k))
!
!   CALCULATE PROBDIFF_CLEAR AND PROBDIFF_CLOUD
!
              DENOM = CLEAR_FRAC+CLOUDF(i,j,k)*LAMBDA_SOOT
              PROBDIFF_CLEAR = CLEAR_FRAC/DENOM
              PROBDIFF_CLOUD = (CLOUDF(i,j,k)*LAMBDA_SOOT)/DENOM
!
!   CALCULATE EXPECTED LIFETIME OF A SOOT PARTICLE WITH
!   RESPECT TO DIFFUSIVE CAPTURE BY CLOUD DROPLETS.
!
              TERM3 = PROBDIFF_FN1*PROBDIFF_CLEAR
              TERM4 = PROBDIFF_FN2*PROBDIFF_CLOUD
              TERM4 = (TERM3+TERM4)*CLEAR_FRAC/CLOUDF(i,j,k)
              DENOM = TERM4*CLOUDTAU + DIFFUSE_TAU
              DIFFRATE = 1.0/DENOM
!
!   Now compute the amount of aged soot converted to cloud soot
!   in this timestep.
              DELTA_NUC(i,j,k)=(1.0-EXP(-DIFFRATE*TIMESTEP))*           &
     &                                                 SOOT_AGD(i,j,k)
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
            delta_sootdiffscav(i,j,k) =                                 &
     &                    delta_nuc(i,j,k) - delta_evap(i,j,k)
          End Do
        End Do
      End Do
!
      Return
      END SUBROUTINE SOOTDIFFSCAV
#endif
