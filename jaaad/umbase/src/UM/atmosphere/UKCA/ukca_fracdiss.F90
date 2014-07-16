#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Calculates the fraction of each species in each
!   dissolved state.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Current Code Owner:       Jamie Rae
!                            Colin Johnson
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
      SUBROUTINE UKCA_FRACDISS(n_pnts,levels,wet_levels,temp,rh,qcl,    &
                               frac_diss,kp_nh)


      USE ASAD_MOD,            ONLY: ddhr, dhr, kd298, k298
      USE UKCA_CONSTANTS,      ONLY: avogadro, rmol
      IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

      INTEGER, INTENT(IN) :: n_pnts         ! No of spatial points
      INTEGER, INTENT(IN) :: levels         ! No of vertical levels
      INTEGER, INTENT(IN) :: wet_levels     ! No of wet levels

      REAL, INTENT(IN)  :: temp(n_pnts,levels)  ! Temperature
      REAL, INTENT(IN)  :: rh(n_pnts,levels)    ! Relative Humidity
      REAL, INTENT(IN)  :: qcl(n_pnts,wet_levels)  ! Cloud liquid water

      REAL, INTENT(OUT) :: frac_diss(n_pnts,levels,jpdw,jpeq+1)
      REAL, INTENT(OUT) :: kp_nh(n_pnts,levels)

! Local variables
      INTEGER :: i                    ! Loop variable
      INTEGER :: k                    ! Loop variable
      INTEGER :: ns                   ! Loop variable

      REAL :: Rgas  !     Ideal gas constant in  atm.l/(mole.K)
      REAL, PARAMETER :: H_plus = 1.0E-05     ! [H+]

      REAL :: tmp1  ! Temporary variable used in calculation
      REAL :: tmp2  ! Temporary variable used in calculation
      REAL :: tmp3  ! Temporary variable used in calculation
      REAL :: inv_298 ! 1/298
      REAL :: Henry(n_pnts,levels,jpdw)     ! Henry's Law consts
      REAL :: Henry_eff(n_pnts,levels,jpdw) ! Effective  ---"---

! Aqueous equilibrium constants:
      REAL :: Equilib_const(n_pnts,levels,jpdw,jpeq)
      REAL :: Frac_aq(n_pnts,levels,jpdw)   ! Aqueous fraction

! For ammonium nitrate dissociation
      REAL :: rhd(n_pnts,wet_levels)        ! rhd deliquesance RH
      REAL :: kp(n_pnts)                    ! kp
      REAL :: kp1(n_pnts)                   ! kp1
      REAL :: p1(n_pnts)                    ! p1
      REAL :: p2(n_pnts)                    ! p2
      REAL :: p3(n_pnts)                    ! p3
      LOGICAL :: todo(n_pnts,wet_levels)

      inv_298 = 1.0/298.0
      rgas = rmol/100.0

      frac_diss(:,:,:,:)=0.0

      DO ns=1,jpdw
        DO k=1,wet_levels
          DO i=1,n_pnts
            tmp1 = (1.0/temp(i,k)) - inv_298
            Henry(i,k,ns) = k298(ns) * EXP(dhr(ns)*tmp1)
            Equilib_const(i,k,ns,1) = kd298(ns)                         &
                                    * EXP(ddhr(ns)*tmp1)
            tmp2 = Equilib_const(i,k,ns,1) / H_plus
            Henry_eff(i,k,ns) = Henry(i,k,ns) * (1.0 + tmp2)
            tmp3 = qcl(i,k) * Henry_eff(i,k,ns) * Rgas * temp(i,k)
            frac_aq(i,k,ns) = tmp3 / (1.0 + tmp3)
            frac_diss(i,k,ns,1) = frac_aq(i,k,ns) * Henry(i,k,ns)       &
                                / Henry_eff(i,k,ns)
            frac_diss(i,k,ns,2) = frac_diss(i,k,ns,1)*tmp2
          END DO
        END DO
      END DO

! Calculate NH4NO3 dissociation constant
      DO k=1,wet_levels
        DO i=1,n_pnts
          rhd(i,k) = exp((618.3/temp(i,k))-2.551)
        ENDDO
      ENDDO

      kp_nh(:,:)=0.0
      todo(:,:)=.false.
      todo = (rh < rhd)
      DO k=1,wet_levels
        DO i=1,n_pnts
          IF (todo(i,k)) THEN
            kp(i)=(temp(i,k)**(-6.025))*EXP(118.87-24084./temp(i,k))
          ELSE
            p1(i)=exp(-135.94+(8763/temp(i,k)))*temp(i,k)**19.12
            p2(i)=exp(-122.65+(9969/temp(i,k)))*temp(i,k)**16.22
            p3(i)=exp(-182.61+(13875/temp(i,k)))*temp(i,k)**24.46
            kp1(i)=(temp(i,k)**(-6.025))*EXP(118.87-24084./temp(i,k))
            kp(i)=(p1(i)-p2(i)*(1.0-rh(i,k))+p3(i)*(1.0-rh(i,k))**2)*   &
                  ((1.0-rh(i,k))**1.75)*kp1(i)
          ENDIF
          kp_nh(i,k)=kp(i)*1.0E-8/((rmol*1.0E6*temp(i,k)/avogadro)**2)
        ENDDO
      ENDDO

      END SUBROUTINE UKCA_FRACDISS
#endif
