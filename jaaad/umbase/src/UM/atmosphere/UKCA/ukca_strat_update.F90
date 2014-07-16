#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Module to contain routines concerned with stratospheric chemistry
!   Contains subroutines: relax_ozone, conserve, ukca_strat_photol,
!   and ukca_calc_ozonecol.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Current Code Owner:       Olaf Morgenstern
!                            Colin Johnson
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
     MODULE ukca_strat_update

     IMPLICIT NONE

     CONTAINS

!  UKCA stratospheric photolysis main routine.
!  Test version
! ######################################################################
!
! Subroutine Interface:
!
!---------------------------------------------------------------------------
! Subroutine STRAT_PHOTOL
!------------------------------------------------------------------------
!
! This routine computes stratospheric photolysis rates and merges the
! rates, where necessary, with the tropospheric rates. This is done for
! one level at a time. The stratospheric photolysis routines are taken
! from SLIMCAT.
!
! Version history:
! ----------------
!
! v1.1 Original code          Olaf Morgenstern 22/1/2004
! v1.2 Introduce pointers for photolysis reactions instead of searching
!      matching ASAD names at every call of the subroutine
!                             Olaf Morgenstern 7/7/2005
! v1.2 Modified to include CO2 photolysis
!                             Olaf Morgenstern 5/1/2006
! v1.3 Modified to include branching from HO2NO2 and BrONO2 photolysis.
!                             Olaf Morgenstern 12/1/2006
! v1.4 Modified to include branching from O2 photolysis
! Current code owner: Olaf Morgenstern.
!
      SUBROUTINE UKCA_STRAT_PHOTOL(pressure,                            &
                          temp,                                         &
                          ozonecol,                                     &
                          cos_zenith_angle,                             &
                          nphot,                                        &
                          photrates)

      USE ukca_photolib  ! SLIMCAT stratospheric photolysis routines
      USE ukca_dissoc    ! Holds photolysis rates (replaces dissoc.h)
      USE ASAD_MOD,      ONLY: jpab, jpat, jpaj, jpah, jpspb, jpspt,    &
                               jpspj, jpsph, spj
      IMPLICIT NONE

#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "parvars.h"
#include "typsize.h"

! Subroutine interface

      INTEGER, INTENT(IN) :: nphot

      REAL, INTENT(IN) :: pressure(row_length, rows)
      REAL, INTENT(IN) :: temp(row_length, rows)
      REAL, INTENT(IN) :: ozonecol(row_length, rows)
      REAL, INTENT(IN) :: cos_zenith_angle(row_length, rows)

      REAL, INTENT(INOUT) :: photrates(row_length, rows, nphot)

! local variables
      INTEGER :: i
      INTEGER :: j
      INTEGER :: l
      INTEGER :: k

      REAL :: frac
      LOGICAL, SAVE :: firstcall = .TRUE.

      INTEGER, SAVE :: ih2o2=0
      INTEGER, SAVE :: ihchoa=0
      INTEGER, SAVE :: ihchob=0
      INTEGER, SAVE :: iho2no2a=0
      INTEGER, SAVE :: ihno3=0
      INTEGER, SAVE :: imeooh=0
      INTEGER, SAVE :: in2o5=0
      INTEGER, SAVE :: ino2=0
      INTEGER, SAVE :: ino3a=0
      INTEGER, SAVE :: ino3b=0
      INTEGER, SAVE :: io2a=0
      INTEGER, SAVE :: io2b=0
      INTEGER, SAVE :: io3a=0
      INTEGER, SAVE :: io3b=0
      INTEGER, SAVE :: io3sa=0
      INTEGER, SAVE :: io3sb=0
      INTEGER, SAVE :: ich4=0
      INTEGER, SAVE :: ih2o=0
      INTEGER, SAVE :: ino=0
      INTEGER, SAVE :: in2o=0
      INTEGER, SAVE :: if11=0
      INTEGER, SAVE :: if12=0
      INTEGER, SAVE :: iclono2a=0
      INTEGER, SAVE :: iclono2b=0
      INTEGER, SAVE :: ihcl=0
      INTEGER, SAVE :: ihocl=0
      INTEGER, SAVE :: ioclo=0
      INTEGER, SAVE :: icl2o2=0
      INTEGER, SAVE :: ibro=0
      INTEGER, SAVE :: ihobr=0
      INTEGER, SAVE :: ibrono2a=0
      INTEGER, SAVE :: ibrcl=0
      INTEGER, SAVE :: imebr=0
      INTEGER, SAVE :: iccl4=0
      INTEGER, SAVE :: if113=0
      INTEGER, SAVE :: imecl=0
      INTEGER, SAVE :: imcfm=0
      INTEGER, SAVE :: if22=0
      INTEGER, SAVE :: ih1211=0
      INTEGER, SAVE :: ih1301=0
      INTEGER, SAVE :: icof2=0
      INTEGER, SAVE :: icofcl=0
      INTEGER, SAVE :: ico2=0
      INTEGER, SAVE :: ibrono2b=0
      INTEGER, SAVE :: iho2no2b=0
      INTEGER, SAVE :: icos=0

! upon first entry initialize positions of photolysis reactions
      IF (firstcall) THEN
        DO k=1,jppj
! Merge tropospheric and stratospheric photolysis rates
          SELECT CASE (spj(k,1))

! H2O2   + h nu -> 2 OH
            CASE ('H2O2      ')
              ih2o2 = k

! consider branching in case of HCHO
            CASE ('HCHO      ')
! HCHO   + h nu -> H2 + CO
              IF ((spj(k,3) == 'H2        ') .OR.                       &
                  (spj(k,4) == 'H2        ')) THEN
                 ihchob = k
              ELSE
! HCHO   + h nu -> H + CHO -> CO + 2 HO2
                 ihchoa = K
              ENDIF

! HO2NO2 + h nu -> NO2 + HO2
            CASE ('HO2NO2    ')
              IF ((spj(k,3) == 'HO2       ')                            &
              .OR.(spj(k,4) == 'HO2       ')) THEN
! HO2NO2 + h nu -> NO2 + HO2
                iho2no2a = k
              ELSE
! HO2NO2 + h nu -> NO3 + OH
                iho2no2b = k
              ENDIF

! HNO3   + h nu -> NO2 + OH
            CASE ('HONO2     ')
              ihno3 = k

! MeOOH  + h nu -> MeO + OH -> HCHO + HO2 + OH
            CASE ('MeOOH     ')
              imeooh = k

! N2O5   + h nu -> NO2 + NO3
            CASE ('N2O5      ')
              in2o5 = k

! NO2    + h nu -> NO  + O
            CASE ('NO2       ')
              ino2 = k

! consider branching for NO3
            CASE ('NO3       ')
              IF ((spj(k,3) == 'O(3P)     ')                            &
              .OR.(spj(k,4) == 'O(3P)     ')) THEN
! NO3    + h nu -> NO2 + O
                ino3b = k
              ELSE
! NO3    + h nu -> NO  + O2
                ino3a = k
              ENDIF

            CASE ('O2        ')
! O2     + h nu -> 2 O / O + O(1D)
              IF ((spj(k,3) == 'O(1D)     ')                            &
                .OR.(spj(k,4) == 'O(1D)     ')) THEN
                io2b = k
              ELSE
                io2a = k
              ENDIF

! consider branching for O3
            CASE ('O3        ')
              IF ((spj(k,3) == 'O(1D)     ')                            &
                .OR.(spj(k,4) == 'O(1D)     ')) THEN
! O3     + h nu -> O2  + O(1D)
                io3a = k
              ELSE
! O3     + h nu -> O2  + O
                io3b = k
              ENDIF

! consider branching for O3S
            CASE ('O3S       ')
              IF ((spj(k,3) == 'O(1D)S    ')                            &
              .OR.(spj(k,4) == 'O(1D)S    ')) THEN
! O3S    + h nu -> O2  + O(1D)S
                io3sa = k
              ELSE
! O3S    + h nu -> O2  + O(3P)S
                io3sb = k
              ENDIF

! CH4    + h nu -> CH3   + H -> MeOO + HO2
            CASE ('CH4       ')
              ich4 = k

! H2O    + h nu -> OH    + H -> OH   + HO2
            CASE ('H2O       ','H2OS      ')
              ih2o = k

! NO     + h nu -> O     + N -> 2O    + NO
            CASE ('NO        ')
              ino = k

! N2O    + h nu -> O(1D) + N2
            CASE ('N2O       ')
              in2o = k

! F11    + h nu -> 3Cl
            CASE ('CFCl3     ')
              if11 = k

! F12    + h nu -> 2Cl
            CASE ('CF2Cl2    ')
              if12 = k

! ClONO2 + h nu -> Cl   + NO3 / ClO + NO2
            CASE ('ClONO2    ')
              IF ((spj(k,3) == 'Cl        ')                            &
              .OR.(spj(k,4) == 'Cl        ')) THEN
! ClONO2 + h nu -> Cl + NO3
                iclono2a = k
              ELSE
! ClONO2 + h nu -> ClO + NO2
                iclono2b = k
              ENDIF

! HCl    + h nu -> H    + Cl
            CASE ('HCl       ')
              ihcl = k

! HOCl   + h nu -> OH   + Cl
            CASE ('HOCl      ')
              ihocl = k

! OClO   + h nu -> O    + ClO
            CASE ('OClO      ')
              ioclo = k

! Cl2O2  + h nu -> 2Cl  + O2
            CASE ('Cl2O2     ')
              icl2o2 = k

! BrO    + h nu -> Br   + O
            CASE ('BrO       ')
              ibro = k

! HOBr   + h nu -> Br   + OH
            CASE ('HOBr      ')
              ihobr = k

! BrONO2 + h nu -> Br   + NO3 / BrO + NO2
            CASE ('BrONO2    ')
              IF ((spj(k,3) == 'Br        ')                            &
              .OR.(spj(k,4) == 'Br        ')) THEN
! BrONO2 + h nu -> Br + NO3
                ibrono2a = k
              ELSE
! BrONO2 + h nu -> BrO + NO2
                ibrono2b = k
              ENDIF

! BrCl   + h nu -> Br   + Cl
            CASE ('BrCl      ')
              ibrcl = k

! MeBr   + h nu -> Br
            CASE ('MeBr      ')
              imebr = k

! CCl4   + h nu -> 4 Cl
            CASE ('CCl4      ')
              iccl4 = k

! CF2ClCFCl2+h nu->3 Cl
            CASE ('CF2ClCFCl2')
              if113 = k

! CH3Cl +h nu   -> Cl
            CASE ('MeCl      ')
              imecl = k

! CH3CCl3 +h nu -> 3 Cl
            CASE ('MeCCl3    ')
              imcfm = k

! CHF2Cl +h nu -> Cl
            CASE ('CHF2Cl    ')
              if22 = k

! CBrClF2 +h nu -> Cl + Br
            CASE ('CF2ClBr   ')
              ih1211 = k

! CBrF3 +h nu   -> Br
            CASE ('CF3Br     ')
              ih1301 = k

! COF2 +h nu    -> 2F + CO
            CASE ('COF2      ')
              icof2 = k

! COFCl +h nu   -> F + Cl
            CASE ('COFCl     ')
              icofcl = k

! CO2 +h nu     -> CO + O(3P)
            CASE ('CO2       ')
              ico2 = k

! COS + h nu    -> CO + S
            CASE ('COS       ')
              icos = k

          END SELECT
        END DO

        CALL inijtab(mype)
        firstcall = .FALSE.

      ENDIF

! CALCJS fills the stratospheric photolysis arrays AJxyz with sensible
! values.

      CALL calcjs(1, rows*row_length,                                   &
                  cos_zenith_angle,                                     &
                  pressure,                                             &
                  temp,                                                 &
                  ozonecol, theta_field_size)

! here: use only existing photolysis reactions where pressure is less than
! 300 hPa, with a linear transition into stratospheric rates

      DO i=1,rows
        DO j=1,row_length
          l=(i-1) * row_length + j

          IF (pressure(j,i) < 30000.) THEN
            IF (pressure(j,i) < 20000.) THEN
              frac = 1.
            ELSe
              frac = (30000. - pressure(j,i))/10000.
            ENDIF

! Merge tropospheric and stratospheric photolysis rates
! H2O2   + h nu -> 2 OH
            photrates(j,i,ih2o2) = frac  * ajh2o2(l)                    &
                               + (1.-frac) * photrates(j,i,ih2o2)

! consider branching in case of HCHO
! HCHO   + h nu -> H2 + CO
            photrates(j,i,ihchob) = frac  * ajc2ob(l)                   &
                               + (1.-frac) * photrates(j,i,ihchob)

! HCHO   + h nu -> H + CHO -> CO + 2 HO2
            photrates(j,i,ihchoa) = frac  * ajc2oa(l)                   &
                               + (1.-frac) * photrates(j,i,ihchoa)

! HO2NO2 + h nu. Branching ratio from JPL (2002)
            IF (iho2no2a > 0) THEN
              IF (iho2no2b > 0) THEN
! HO2NO2 + h nu -> NO2 + HO2
                photrates(j,i,iho2no2a) = 0.667 * frac  * ajpna(l)      &
                         + (1.-frac) * photrates(j,i,iho2no2a)
! HO2NO2 + h nu -> NO3 + OH
                photrates(j,i,iho2no2b) = 0.333 * frac  * ajpna(l)      &
                         + (1.-frac) * photrates(j,i,iho2no2b)
              ELSE
                photrates(j,i,iho2no2a) = frac  * ajpna(l)              &
                         + (1.-frac) * photrates(j,i,iho2no2a)
              END IF
            END IF

! HNO3   + h nu -> NO2 + OH
            photrates(j,i,ihno3) = frac  * ajhno3(l)                    &
                               + (1.-frac) * photrates(j,i,ihno3)

! MeOOH  + h nu -> MeO + OH -> HCHO + HO2 + OH
            photrates(j,i,imeooh) = frac  * ajmhp(l)                    &
                               + (1.-frac) * photrates(j,i,imeooh)

! N2O5   + h nu -> NO2 + NO3
            photrates(j,i,in2o5) = frac  * ajn2o5(l)                    &
                               + (1.-frac) * photrates(j,i,in2o5)

! NO2    + h nu -> NO  + O
            photrates(j,i,ino2) = frac  * ajno2(l)                      &
                               + (1.-frac) * photrates(j,i,ino2)

! consider branching for NO3
! NO3    + h nu -> NO2 + O
            photrates(j,i,ino3b) = frac  * ajno32(l)                    &
                               + (1.-frac) * photrates(j,i,ino3b)

! NO3    + h nu -> NO  + O2
            photrates(j,i,ino3a) = frac  * ajno31(l)                    &
                               + (1.-frac) * photrates(j,i,ino3a)

! O2     + h nu -> 2 O
            photrates(j,i,io2a) = frac  * aj2a(l)                       &
                               + (1.-frac) * photrates(j,i,io2a)

! consider branching for O3
! O3     + h nu -> O2  + O(1D)
            photrates(j,i,io3a) = frac  * aj3a(l)                       &
                                 + (1.-frac) * photrates(j,i,io3a)
! O3     + h nu -> O2  + O
            photrates(j,i,io3b) = frac  * aj3(l)                        &
                                 + (1.-frac) * photrates(j,i,io3b)

! consider branching for O3S
! O3S    + h nu -> O2  + O(1D)S
            IF (io3sa > 0)                                              &
              photrates(j,i,io3sa) = frac  * aj3a(l)                    &
                                 + (1.-frac) * photrates(j,i,io3sa)
! O3S    + h nu -> O2  + O(3P)S
            IF (io3sb > 0)                                              &
              photrates(j,i,io3sb) = frac  * aj3(l)                     &
                                 + (1.-frac) * photrates(j,i,io3sb)

          ENDIF

! purely stratospheric photolysis rates
! SLIMCAT specific photolysis reactions

! O2     + h nu -> O + O(1D)
          IF (io2b > 0) THEN
            photrates(j,i,io2b) = aj2b(l)
          ELSE
            photrates(j,i,io2a) = photrates(j,i,io2a) + aj2b(l)
          ENDIF
! CH4    + h nu -> CH3   + H -> MeOO + HO2
          IF (ich4 > 0)    photrates(j,i,ich4) = ajch4(l)

! H2O    + h nu -> OH    + H -> OH   + HO2
          IF (ih2o > 0)    photrates(j,i,ih2o) = ajh2o(l)

! NO     + h nu -> O     + N -> 2O    + NO
          IF (ino > 0)     photrates(j,i,ino)  = ajno(l)

! N2O    + h nu -> O(1D) + N2
          IF (in2o > 0)    photrates(j,i,in2o) = ajn2o(l)

! F11    + h nu -> 3Cl
          IF (if11 > 0)    photrates(j,i,if11) = ajf11(l)

! F12    + h nu -> 2Cl
          IF (if12 > 0)    photrates(j,i,if12) = ajf12(l)

! ClONO2 + h nu -> Cl   + NO3
          IF (iclono2a > 0) THEN
            photrates(j,i,iclono2a) = ajcnita(l)
            IF (iclono2b == 0)                                          &
               photrates(j,i,iclono2a) = photrates(j,i,iclono2a)        &
                                       + ajcnitb(l)
          ENDIF

! ClONO2 + h nu -> ClO  + NO2
          IF (iclono2b > 0) THEN
            photrates(j,i,iclono2b) = ajcnitb(l)
            IF (iclono2a == 0)                                          &
               photrates(j,i,iclono2b) = photrates(j,i,iclono2b)        &
                                       + ajcnita(l)
          ENDIF

! HCl    + h nu -> H    + Cl
          IF (ihcl > 0)    photrates(j,i,ihcl) = ajhcl(l)

! HOCl   + h nu -> OH   + Cl
          IF (ihocl > 0)   photrates(j,i,ihocl) = ajhocl(l)

! OClO   + h nu -> O    + ClO
          IF (ioclo > 0)   photrates(j,i,ioclo) = ajoclo(l)

! Cl2O2  + h nu -> 2Cl  + O2
          IF (icl2o2 > 0)  photrates(j,i,icl2o2) = ajcl2o2(l)

! BrO    + h nu -> Br   + O
          IF (ibro > 0)    photrates(j,i,ibro) = ajbro(l)

! HOBr   + h nu -> Br   + OH
          IF (ihobr > 0)   photrates(j,i,ihobr) = ajhobr(l)

! BrONO2 + h nu -> Br   + NO3 / BrO + NO2
! Consider branching ratio (JPL, 2002)
          IF (ibrono2a*ibrono2b > 0) THEN
            photrates(j,i,ibrono2a) = 0.29*ajbrno3(l) ! Br + NO3 channel
            photrates(j,i,ibrono2b) = 0.71*ajbrno3(l) ! BrO + NO2 channel
          ELSEIF (ibrono2a > 0) THEN ! no branching
            photrates(j,i,ibrono2a) = ajbrno3(l)
          ELSEIF (ibrono2b > 0) THEN
            photrates(j,i,ibrono2b) = ajbrno3(l)
          END IF

! BrCl   + h nu -> Br   + Cl
          IF (ibrcl > 0)   photrates(j,i,ibrcl) = ajbrcl(l)

! MeBr   + h nu -> Br
          IF (imebr > 0)   photrates(j,i,imebr) = ajch3br(l)

! CCl4   + h nu -> 4 Cl
          IF (iccl4 > 0)   photrates(j,i,iccl4) = ajccl4(l)

! CF2ClCFCl2+h nu->3 Cl
          IF (if113 > 0)   photrates(j,i,if113) = ajf113(l)

! CH3Cl +h nu   -> Cl
          IF (imecl > 0)   photrates(j,i,imecl) = ajch3cl(l)

! CH3CCl3 +h nu -> 3 Cl
          IF (imcfm > 0)   photrates(j,i,imcfm) = ajmcfm(l)

! CHF2Cl +h nu -> Cl
          IF (if22 > 0)    photrates(j,i,if22) = ajf22(l)

! CBrClF2 +h nu -> Cl + Br
          IF (ih1211 > 0)  photrates(j,i,ih1211) = ajf12b1(l)

! CBrF3 +h nu   -> Br
          IF (ih1301 > 0)  photrates(j,i,ih1301) = ajf13b1(l)

! COF2 +h nu    -> 2F + CO
          IF (icof2 > 0)   photrates(j,i,icof2) = ajcof2(l)

! COFCl +h nu   -> F + Cl
          IF (icofcl > 0)  photrates(j,i,icofcl) = ajcofcl(l)

! CO2 +h nu     -> CO + O(3P)
          IF (ico2 > 0)    photrates(j,i,ico2) = ajco2(l)

! COS +h nu     -> CO + S
          IF (icos > 0)    photrates(j,i,icos) = ajcos(l)

        ENDDO
      ENDDO


      RETURN
      END SUBROUTINE UKCA_STRAT_PHOTOL

!------------------------------------------------------------------
! Subroutine CALC_OZONECOL
!------------------------------------------------------------------
!
! This routine calculates overhead ozone columns in molecules/m^2 needed
! for photolysis calculations using the Lary scheme of SLIMCAT. The
! contribution from the area above the model top is correct for a model top
! of 65 km but is not changed here for simplicity.
!
      SUBROUTINE UKCA_CALC_OZONECOL(model_levels, rows, row_length,     &
                          p_layer_boundaries,                           &
                          p_layer_centres,                              &
                          ozone_vmr,                                    &
                          ozonecol)

      IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: model_levels
      INTEGER, INTENT(IN) :: rows
      INTEGER, INTENT(IN) :: row_length

      REAL, INTENT(IN) :: p_layer_boundaries(row_length, rows,          &
                                             0:model_levels)
      REAL, INTENT(IN) :: p_layer_centres(row_length, rows, model_levels)
      REAL, INTENT(IN) :: ozone_vmr(row_length, rows, model_levels)

      REAL, INTENT(OUT) :: ozonecol(row_length, rows, model_levels)

! local variables

! Ozone column at model top. At 38 km, calculated from 60-level ozone clima
! tology. At 85 km, extrapolated (hence a wild guess but should not matter
! too much...

      REAL, PARAMETER :: ozcol_39km = 5.0E17
      REAL, PARAMETER :: ozcol_85km = 6.7E13

      REAL, PARAMETER :: colfac = 2.132e20

      INTEGER :: l

      ozonecol = 0.
      DO l=model_levels-1,1,-1

! compute the contributions from layers above l
        ozonecol(:,:,l) = ozonecol(:,:,l+1) +                           &
            ozone_vmr(:,:,l+1) *                                        &
             (p_layer_boundaries(:,:,l) - p_layer_boundaries(:,:,l+1))  &
             * colfac
      END DO

! add contribution within top of level L

      ozonecol = ozonecol +                                             &
            ozone_vmr *                                                 &
           (p_layer_centres - p_layer_boundaries(:,:,1:model_levels))   &
            * colfac

! add contribution above model top. For 38 levels assume model top at
! 39 km with a climatological ozone column there of 5E17 molecules/cm^2.
! Otherwise, the top is assumed at 85 km with an ozone column there of
! 6.7E13 moelcules/cm^2.

      IF (model_levels == 38) THEN
        ozonecol = ozonecol + ozcol_39km
      ELSE
        ozonecol = ozonecol + ozcol_85km
      ENDIF

      RETURN
      END SUBROUTINE UKCA_CALC_OZONECOL


      SUBROUTINE conserve(row_length, rows, model_levels,               &
                          tr_vars, advt, tracers,                       &
                          pres, climoz, direction)
! Description:
!
! This routine calculates and conserves total chlorine, bromine, and
! hydrogen. For these elements closed chemistry should be prescribed.
! Called before chemistry, with direction = .TRUE., it calculates
! total bromine, chlorine, and hydrogen as 3-D fields. Called afer
! chemistry, with direction = .FALSE., it rescales the chlorine, bromine
! and hydrogen containing compounds so that total chlorine, bromine
! and hydrogen are conserved under chemistry. Where a compound contains
! more than one of the 3 elements. e.g, BrCl, it is only scaled for the
! less abundant of the two constituents, Br. It is then subtracted from
! the total chlorine.
!
! Method: Rescaling of tracer variables.
!
! Current code owner: Olaf Morgenstern (olaf.morgenstern@atm.ch.cam.ac.uk)
!
! History:
!
!  v1.1 26/10/2005 Original code  Olaf Morgenstern
!  v1.2 23/12/2005 Also include upper bounday condidition. Ox and RO3
!                  in the top 2 levels (for L60 version)
!                  are set equal to climatological
!                  ozone.         Olaf Morgenstern
!  v1.3 20/06/2006 Minor changes for code merging O. Morgenstern
!  v1.4 29/12/2006 Remove special code for H2OS   O. Morgenstern
! Code Description:
! Language: FORTRAN 90 + common extensions.
!
! Declarations:
! These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Global variables (*CALLed COMDECKs etc...):


! 6/8/2004    Olaf Morgenstern


      IMPLICIT NONE

#include "c_v_m.h"
#include "nstypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"

! Subroutine interface
      INTEGER, INTENT(IN) :: row_length
      INTEGER, INTENT(IN) :: rows
      INTEGER, INTENT(IN) :: model_levels
      INTEGER, INTENT(IN) :: tr_vars

      LOGICAL, INTENT(IN) :: direction

      CHARACTER*10, INTENT(IN) :: advt(jpctr)

      REAL, INTENT(IN) :: pres(row_length, rows, model_levels)

      REAL, INTENT(IN) :: climoz(row_length, rows, model_levels)

      REAL, INTENT(INOUT) :: tracers(row_length, rows, model_levels,    &
                                     tr_vars)

! Local variables

! Maximum number of Br, Cl, and H compounds permitted
      INTEGER, PARAMETER :: max_comp = 50

      REAL, PARAMETER :: adjust_level = 500. ! pressure below which allow
! for changes of total hydrogen due to dehydration, and above which only
! advective changes are allowed.

      REAL, PARAMETER :: washout_limit = 10000. ! pressure limit below which
! hydrogen conservation is not enforced.

      INTEGER :: m,i,j,k

! Reservoir tracer for chlorine
      INTEGER, SAVE :: n_hcl=0
! Reservoir tracer for bromine
      INTEGER, SAVE :: n_brx=0
      INTEGER, SAVE :: n_toth=0      ! position of total hydrogen tracer
      INTEGER, SAVE :: n_h2o=0
      INTEGER, SAVE :: n_h2os=0
      INTEGER, SAVE :: n_ro3=0
      INTEGER, SAVE :: n_ox=0
      INTEGER, SAVE :: n_n2o=0

      LOGICAL, SAVE :: firstcall = .TRUE. ! flag for first call of subr.

      INTEGER, SAVE ::       ncl_tracers     ! number of Cl tracers
      INTEGER, SAVE ::       nbr_tracers     ! number of Br tracers
      INTEGER, SAVE ::       nh_tracers      ! number of H tracers

      REAL, ALLOCATABLE, SAVE :: total_cl(:,:,:) ! total chlorine VMR
      REAL, ALLOCATABLE, SAVE :: total_br(:,:,:) ! total bromine VMR
      REAL, ALLOCATABLE, SAVE :: total_h (:,:,:) ! total hydrogen VMR


      INTEGER, SAVE ::  cl_tracers(max_comp) ! positions of Cl tracers
      INTEGER, SAVE ::  br_tracers(max_comp) ! positions of Br tracers
      INTEGER, SAVE ::   h_tracers(max_comp) ! positions of hydrogen tracers

      REAL   , SAVE ::c_cl_tracers(max_comp) ! conversion factors VMR/MMR
      REAL   , SAVE ::c_br_tracers(max_comp) ! conversion factors
      REAL   , SAVE ::c_h_tracers(max_comp)  ! conversion factors

      REAL   , SAVE :: cl_validity(max_comp) ! number of Cl atoms per mol.
      REAL   , SAVE :: br_validity(max_comp) ! number of Br atoms per mol.

      REAL   , SAVE :: h_validity(max_comp)  ! number of H atoms per mol.

      LOGICAL, SAVE ::contains_bromine(max_comp)! flag for Cl compounds als
                                             ! cont. bromine
      LOGICAL, SAVE ::do_not_change(max_comp)! leave unchanged

! Correction factor to achieve Cl or Br conservation
      REAL, ALLOCATABLE :: corrfac(:,:,:)

      IF (firstcall) THEN
!
! Calculate the number of Cl and Br tracers, their positions, mmr/vmr
! conversion ratios, and validities (numbers of Cl/Br atoms per molecule).
!
        ncl_tracers = 0
        nbr_tracers = 0
        nh_tracers = 0
        contains_bromine = .FALSE.
        do_not_change = .FALSE.
        cl_tracers = 0
        br_tracers = 0
        h_tracers = 0
        cl_validity = 1.
        br_validity = 1.
        h_validity = 1.
        c_br_tracers = 0.
        c_cl_tracers = 0.
        c_h_tracers = 0.
        DO m=1,jpctr
          SELECT CASE (advt(m))

! chlorine tracers
            CASE ('CF2Cl2    ')
              ncl_tracers = ncl_tracers + 1
              cl_tracers(ncl_tracers) = m
              c_cl_tracers(ncl_tracers) = c_cf2cl2
              cl_validity(ncl_tracers) = 2.
            CASE ('CFCl3     ')
              ncl_tracers = ncl_tracers + 1
              cl_tracers(ncl_tracers) = m
              c_cl_tracers(ncl_tracers) = c_cfcl3
              cl_validity(ncl_tracers) = 3.
            CASE ('Clx       ','ClO       ')
              ncl_tracers = ncl_tracers + 1
              cl_tracers(ncl_tracers) = m
              c_cl_tracers(ncl_tracers) = c_clo
            CASE ('Cl        ')
              ncl_tracers = ncl_tracers + 1
              cl_tracers(ncl_tracers) = m
              c_cl_tracers(ncl_tracers) = c_cl
            CASE ('Cl2O2     ')
              ncl_tracers = ncl_tracers + 1
              cl_tracers(ncl_tracers) = m
              c_cl_tracers(ncl_tracers) = c_cl2o2
              cl_validity(ncl_tracers) = 2.
            CASE ('HCl       ')
              ncl_tracers = ncl_tracers + 1
              cl_tracers(ncl_tracers) = m
              c_cl_tracers(ncl_tracers) = c_hcl
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_hcl
              do_not_change(nh_tracers) = .TRUE.
              n_hcl = m
            CASE ('HOCl      ')
              ncl_tracers = ncl_tracers + 1
              cl_tracers(ncl_tracers) = m
              c_cl_tracers(ncl_tracers) = c_hocl
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_hocl
              do_not_change(nh_tracers) = .TRUE.
            CASE ('OClO      ')
              ncl_tracers = ncl_tracers + 1
              cl_tracers(ncl_tracers) = m
              c_cl_tracers(ncl_tracers) = c_oclo
            CASE ('BrCl      ')
              ncl_tracers = ncl_tracers + 1
              cl_tracers(ncl_tracers) = m
              c_cl_tracers(ncl_tracers) = c_brcl
              contains_bromine(ncl_tracers) = .TRUE.
              nbr_tracers = nbr_tracers + 1
              br_tracers(nbr_tracers) = m
              c_br_tracers(nbr_tracers) = c_brcl
            CASE ('ClONO2    ')
              ncl_tracers = ncl_tracers + 1
              cl_tracers(ncl_tracers) = m
              c_cl_tracers(ncl_tracers) = c_clono2
            CASE ('CF2ClCFCl2')
              ncl_tracers = ncl_tracers + 1
              cl_tracers(ncl_tracers) = m
              c_cl_tracers(ncl_tracers) = c_cf2clcfcl2
              cl_validity(ncl_tracers) = 3.
            CASE ('CHF2Cl    ')
              ncl_tracers = ncl_tracers + 1
              cl_tracers(ncl_tracers) = m
              c_cl_tracers(ncl_tracers) = c_chf2cl
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_chf2cl
              do_not_change(nh_tracers) = .TRUE.
            CASE ('MeCCl3    ')
              ncl_tracers = ncl_tracers + 1
              cl_tracers(ncl_tracers) = m
              c_cl_tracers(ncl_tracers) = c_meccl3
              cl_validity(ncl_tracers) = 3.
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_meccl3
              do_not_change(nh_tracers) = .TRUE.
            CASE ('CCl4      ')
              ncl_tracers = ncl_tracers + 1
              cl_tracers(ncl_tracers) = m
              c_cl_tracers(ncl_tracers) = c_ccl4
              cl_validity(ncl_tracers) = 4.
            CASE ('MeCl      ')
              ncl_tracers = ncl_tracers + 1
              cl_tracers(ncl_tracers) = m
              c_cl_tracers(ncl_tracers) = c_mecl
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_mecl
              do_not_change(nh_tracers) = .TRUE.
            CASE ('CF2ClBr   ')
              ncl_tracers = ncl_tracers + 1
              cl_tracers(ncl_tracers) = m
              c_cl_tracers(ncl_tracers) = c_cf2clbr
              contains_bromine(ncl_tracers) = .TRUE.
              nbr_tracers = nbr_tracers + 1
              br_tracers(nbr_tracers) = m
              c_br_tracers(nbr_tracers) = c_cf2clbr
! Bromine tracers
            CASE ('Brx       ','BrO       ')
              nbr_tracers = nbr_tracers + 1
              br_tracers(nbr_tracers) = m
              c_br_tracers(nbr_tracers) = c_bro
              n_brx = m
            CASE ('Br        ')
              nbr_tracers = nbr_tracers + 1
              br_tracers(nbr_tracers) = m
              c_br_tracers(nbr_tracers) = c_br
            CASE ('HOBr      ')
              nbr_tracers = nbr_tracers + 1
              br_tracers(nbr_tracers) = m
              c_br_tracers(nbr_tracers) = c_hobr
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_hobr
              do_not_change(nh_tracers) = .TRUE.
            CASE ('BrONO2    ')
              nbr_tracers = nbr_tracers + 1
              br_tracers(nbr_tracers) = m
              c_br_tracers(nbr_tracers) = c_brono2
            CASE ('MeBr      ')
              nbr_tracers = nbr_tracers + 1
              br_tracers(nbr_tracers) = m
              c_br_tracers(nbr_tracers) = c_mebr
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_mebr
              do_not_change(nh_tracers) = .TRUE.
            CASE ('HBr       ')
              nbr_tracers = nbr_tracers + 1
              br_tracers(nbr_tracers) = m
              c_br_tracers(nbr_tracers) = c_hbr
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_hbr
              do_not_change(nh_tracers) = .TRUE.
            CASE ('CF3Br     ')
              nbr_tracers = nbr_tracers + 1
              br_tracers(nbr_tracers) = m
              c_br_tracers(nbr_tracers) = c_cf3br
! hydrogen tracers
            CASE ('TOTH      ')
              n_toth = m
            CASE ('H2O       ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_h2o
              h_validity(nh_tracers) = 2.
              n_h2o = m
            CASE ('CH4       ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_ch4
              h_validity(nh_tracers) = 4.
              do_not_change(nh_tracers) = .TRUE.
            CASE ('H2        ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_h2
              h_validity(nh_tracers) = 2.
            CASE ('HO2NO2    ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_ho2no2
              do_not_change(nh_tracers) = .TRUE.
            CASE ('HONO2     ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_hono2
              do_not_change(nh_tracers) = .TRUE.
            CASE ('H2O2      ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_h2o2
              h_validity(nh_tracers) = 2.
            CASE ('HCHO      ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_hcho
              h_validity(nh_tracers) = 2.
            CASE ('MeOOH     ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_meooh
              h_validity(nh_tracers) = 4.
            CASE ('HONO      ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_hono
              do_not_change(nh_tracers) = .TRUE.
            CASE ('C2H6      ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_c2h6
              h_validity(nh_tracers) = 6.
            CASE ('EtOOH     ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_etooh
              h_validity(nh_tracers) = 6.
            CASE ('MeCHO     ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_mecho
              h_validity(nh_tracers) = 4.
            CASE ('PAN       ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_pan
              h_validity(nh_tracers) = 3.
              do_not_change(nh_tracers) = .TRUE.
            CASE ('C3H8      ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_c3h8
              h_validity(nh_tracers) = 8.
            CASE ('n-PrOOH   ','i-PrOOH   ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_prooh
              h_validity(nh_tracers) = 8.
            CASE ('EtCHO     ','Me2CO     ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_etcho
              h_validity(nh_tracers) = 6.
            CASE ('MeCOCH2OOH')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_mecoch2ooh
              h_validity(nh_tracers) = 6.
            CASE ('PPAN      ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_ppan
              h_validity(nh_tracers) = 5.
              do_not_change(nh_tracers) = .TRUE.
            CASE ('MeONO2    ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_meono2
              h_validity(nh_tracers) = 3.
              do_not_change(nh_tracers) = .TRUE.
            CASE ('H         ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_h
            CASE ('OH        ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_oh
            CASE ('HO2       ','HOx       ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_ho2
            CASE ('MeOO      ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_meoo
              h_validity(nh_tracers) = 3.
            CASE ('EtOO      ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_etoo
              h_validity(nh_tracers) = 5.
            CASE ('i-PrOO    ','n-PrOO    ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_proo
              h_validity(nh_tracers) = 7.
            CASE ('EtCO3     ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_etco3
              h_validity(nh_tracers) = 5.
            CASE ('MeCOCH2OO ')
              nh_tracers = nh_tracers + 1
              h_tracers(nh_tracers) = m
              c_h_tracers(nh_tracers) = c_mecoch2oo
              h_validity(nh_tracers) = 5.
! special tracers
            CASE('Ox        ','O3        ')
              n_ox = m
            CASE('RO3       ')
              n_ro3 = m
            CASE('N2O       ')
              n_n2o = m
         END SELECT
        END DO

        firstcall = .FALSE.

      END IF  ! if firstcall

      IF (direction) THEN
! Calculate total chlorine and bromine vmrs.
        ALLOCATE(total_br(row_length, rows, model_levels))
        ALLOCATE(total_cl(row_length, rows, model_levels))
        total_br = 0.
        total_cl = 0.

        DO m=1,nbr_tracers
          total_br = total_br + tracers(:,:,:,br_tracers(m))*           &
                                br_validity(m) / c_br_tracers(m)
        END DO
        DO m=1,ncl_tracers
          total_cl = total_cl + tracers(:,:,:,cl_tracers(m))*           &
                                cl_validity(m) / c_cl_tracers(m)
        END DO
! Do not do hydrogen conservation unless at least two hydrogen
! reservoirs (H2O, CH4) are defined.
        IF (L_ukca_h2o_feedback) THEN
          ALLOCATE(total_h(row_length, rows, model_levels))
          total_h = 0.
          DO m=1,nh_tracers
            total_h = total_h + tracers(:,:,:,h_tracers(m))*            &
                                h_validity(m) / c_h_tracers(m)
          END DO
          IF (n_toth > 0) THEN
            WHERE (pres > adjust_level)                                 &
              tracers(:,:,:,n_toth) = total_h
            WHERE (pres <= adjust_level)                                &
              total_h = tracers(:,:,:,n_toth)
          END IF
        END IF

      ELSE
        IF (model_levels > 38) THEN
          IF (n_ro3 > 0)                                                &
            tracers(:,:,model_levels-1:model_levels,n_ro3) =            &
              climoz(:,:,model_levels-1:model_levels)
          tracers(:,:,model_levels-1:model_levels,n_ox) =               &
            climoz(:,:,model_levels-1:model_levels)
          IF (n_n2o > 0)                                                &
            tracers(:,:,model_levels,n_n2o) = 0.
        END IF

! Adjust tracers to match
        ALLOCATE(corrfac(row_length, rows, model_levels))
        IF (nbr_tracers > 0) THEN
          corrfac = 0.

! Calculate new total bromine
          DO m=1,nbr_tracers
            corrfac = corrfac + tracers(:,:,:,br_tracers(m))*           &
                                  br_validity(m) / c_br_tracers(m)
          END DO

! Adjust bromine tracers to match total bromine computed before
! chemistry

          corrfac = total_br / corrfac
          DO m=1,nbr_tracers
          tracers(:,:,:,br_tracers(m)) = tracers(:,:,:,br_tracers(m))   &
            * corrfac
          END DO
          IF ((minval(corrfac) < 0.9) .OR. (maxval(corrfac) > 1.1)) THEN
#if defined(IBM)
            WRITE (0,*) 'Correct bromine ',minval(corrfac),             &
                        maxval(corrfac)
#else
            WRITE (6,*) 'Correct bromine ',minval(corrfac),             &
                        maxval(corrfac)
!DEPENDS ON: um_fort_flush
            CALL um_fort_flush(101,1)
#endif
          END IF
        END IF
! Calculate new total chlorine, excluding BrCl and CF2ClBr

        IF (ncl_tracers > 0) THEN
          corrfac = 0.

          DO m=1,ncl_tracers
            IF (.NOT.(contains_bromine(m))) THEN
              corrfac = corrfac + tracers(:,:,:,cl_tracers(m))*         &
                                cl_validity(m) / c_cl_tracers(m)
            ELSE
              total_cl = total_cl - tracers(:,:,:,cl_tracers(m))*       &
                                cl_validity(m) / c_cl_tracers(m)
            END IF
          END DO

! Adjust chlorine species to match total chlorine computed before.
! Leave BrCl and CF2ClBr alone.

          corrfac = total_cl / corrfac
          DO m=1,ncl_tracers
            IF (.NOT.(contains_bromine(m)))                             &
              tracers(:,:,:,cl_tracers(m)) =                            &
                   tracers(:,:,:,cl_tracers(m)) * corrfac
          END DO
        IF ((minval(corrfac) < 0.9) .OR. (maxval(corrfac) > 1.1)) THEN
#if defined(IBM)
          WRITE(0,*)'Correct chlorine ',minval(corrfac),maxval(corrfac)
#else
          WRITE(*,*)'Correct chlorine ',minval(corrfac),maxval(corrfac)
! DEPENDS ON: um_fort_flush
          CALL um_fort_flush(101,1)
#endif
        END IF
        END IF

        IF (L_ukca_h2o_feedback) THEN
          corrfac = 0.

! Calculate new total hydrogen
          DO m=1,nh_tracers
            IF (.NOT.(do_not_change(m))) THEN
              corrfac = corrfac + tracers(:,:,:,h_tracers(m))*          &
                                  h_validity(m) / c_h_tracers(m)
            ELSE
              total_h = total_h - tracers(:,:,:,h_tracers(m))*          &
                                h_validity(m) / c_h_tracers(m)
            ENDIF
          END DO

! adjust upper boundary for hydrogen tracers only. This is needed to
! prevent model instability if too much water is present at the model
! top.
! Adjust hydrogen tracers to match total hydrogen computed before
! chemistry

          corrfac = total_h / corrfac

! Do not enforce hydrogen conservation in the troposphere, due to
! washout and dry deposition.

          WHERE (pres > washout_limit) corrfac = 1.

          DO m=1,nh_tracers
            IF (.NOT.(do_not_change(m)))                                &
              tracers(:,:,:,h_tracers(m)) = tracers(:,:,:,h_tracers(m)) &
                * corrfac
          END DO
          DEALLOCATE(total_h)
        END IF
! Deallocate fields
        DEALLOCATE(total_br)
        DEALLOCATE(total_cl)
        DEALLOCATE(corrfac)
      END IF

      RETURN
      END SUBROUTINE conserve

      SUBROUTINE relax_ozone(rows, row_length, model_levels,            &
                             pressure, ozone_mmr, climoz)

!
!  This subroutine replaces chemical ozone with climatological ozone
!  above 0.2 hPa, with a transition region between 0.2 and 0.3 hPa.
!  This is meant to improve ozone in this region and allow for better
!  transmission of UV radiation through the mesosphere.
!
!  O. Morgenstern      21/9/2005
!
      IMPLICIT NONE

! Subroutine interface

      INTEGER, INTENT(IN) :: rows
      INTEGER, INTENT(IN) :: row_length
      INTEGER, INTENT(IN) :: model_levels
      REAL, INTENT(IN) :: pressure(row_length, rows, model_levels)
      REAL, INTENT(IN) :: climoz(row_length, rows, model_levels)

      REAL, INTENT(INOUT) :: ozone_mmr(row_length, rows, model_levels)

! Local vaiables:

      REAL, PARAMETER :: maxpres = 30.
      REAL, PARAMETER :: minpres = 20.

      INTEGER :: i
      INTEGER :: j
      INTEGER :: k

      REAL :: frac

      DO k=1,model_levels
        IF (MINVAL(pressure(:,:,k)) < maxpres) THEN
          DO i=1,rows
            DO j=1,row_length
              IF (pressure(j,i,k) < maxpres) THEN
                IF (pressure(j,i,k) < minpres) THEN
                  ozone_mmr(j,i,k) = climoz(j,i,k)
                ELSE
                  frac = (maxpres-pressure(j,i,k))/(maxpres - minpres)
                  ozone_mmr(j,i,k) = frac  *    climoz(j,i,k) +        &
                               (1. - frac) * ozone_mmr(j,i,k)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE relax_ozone

      END MODULE ukca_strat_update
#endif
