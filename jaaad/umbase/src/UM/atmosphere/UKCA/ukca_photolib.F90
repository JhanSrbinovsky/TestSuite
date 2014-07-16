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
!     Module to hold collection of routines for stratospheric photolysis
!     in UKCA, originally written for the SLIMCAT stratospheric CTM.
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
      MODULE ukca_photolib

      USE UKCA_CONSTANTS, ONLY: pi, pi_over_180, Earth_radius
      IMPLICIT NONE

      REAL, PARAMETER :: RE = Earth_radius/1.0E3

      CONTAINS
!======================================================================
      SUBROUTINE CALCJS(IHMIN,IHMAX,F,PPA,T,O3COL,LNT)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                               C
!     LOOK UP PHOTOLYSIS RATES FROM  J TABLES                   C
!                                                               C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! 03/01/2006  Modified to include CO2 photolysis
!                                  Olaf Morgenstern

      USE UKCA_DISSOC
      IMPLICIT NONE

#include "parpho.h"
#include "tabjs.h"
#include "parvars.h"

! Subroutine interface
      INTEGER, INTENT(IN) :: IHMIN        ! first point
      INTEGER, INTENT(IN) :: IHMAX        ! last point
      INTEGER, INTENT(IN) :: LNT          ! No of points

      REAL, INTENT(IN) :: F(LNT)          ! Cos of zenith angle
      REAL, INTENT(IN) :: T(LNT)          ! Temperature
      REAL, INTENT(IN) :: O3COL(LNT)      ! Ozone column
      REAL, INTENT(IN) :: PPA(LNT)        ! Pressure

! Local variables
      LOGICAL, save :: init = .FALSE.

      INTEGER :: IH

      REAL :: CSUP
      REAL :: DLNP
      REAL :: TDIF
      REAL :: O3DIF
      REAL :: RADINT
      REAL :: COSINT

      INTEGER :: JX(LNT)
      INTEGER :: JXP1(LNT)
      INTEGER :: JY(LNT)
      INTEGER :: JYP1(LNT)
      INTEGER :: JZ(LNT)
      INTEGER :: JZP1(LNT)
      INTEGER :: JO(LNT)
      INTEGER :: JOP1(LNT)

      REAL :: ANGLE(LNT)
      REAL :: O3FAC(LNT)
      REAL :: P(LNT)
      REAL :: TU(LNT)
      REAL :: ZO (LNT)
      REAL :: ZX (LNT)
      REAL :: ZY (LNT)
      REAL :: ZZ (LNT)
      REAL :: CO (LNT)
      REAL :: CX (LNT)
      REAL :: CY (LNT)
      REAL :: CZ (LNT)
      REAL :: W1 (LNT)
      REAL :: W2 (LNT)
      REAL :: W3 (LNT)
      REAL :: W4 (LNT)
      REAL :: W5 (LNT)
      REAL :: W6 (LNT)
      REAL :: W7 (LNT)
      REAL :: W8 (LNT)
      REAL :: W9 (LNT)
      REAL :: W10(LNT)
      REAL :: W11(LNT)
      REAL :: W12(LNT)
      REAL :: W13(LNT)
      REAL :: W14(LNT)
      REAL :: W15(LNT)
      REAL :: W16(LNT)
      REAL :: SAO3C(JPLEV)

! do not read JTABLE file. Instead, call inijtab routine.
!
      IF (init) THEN
!
!       Open file for the photolysis tables and and read them.
        OPEN(NTAB,FILE='jtable',FORM='UNFORMATTED',STATUS='OLD')
!
        REWIND(NTAB)
!
!     Read the uv/vis radiation field.
        READ(NTAB) TABPRES ,TABANG   ,TABT     ,TABO3    ,SAO3C    ,    &
                   TABJ2A  ,TABJ2B   ,TABJ3    ,TABJCH4  ,              &
                   TABJ3A  ,TABJCCL4 ,TABJCL2O2,TABJCNITA,TABJCNITB,    &
                   TABJF11 ,                                            &
                   TABJF113,TABJF12  ,TABJF22  ,TABJH2O  ,TABJH2O2 ,    &
                   TABJHCL ,TABJHNO3 ,TABJHOCL ,TABJCH3CL,TABJN2O  ,    &
                   TABJN2O5,TABJNO   ,TABJNO2  ,TABJNO31 ,TABJNO32 ,    &
                   TABJPNA ,TABJBRCL ,TABJBRNO3,TABJBRO  ,TABJMCFM ,    &
                   TABJHOBR,TABJOCLO ,TABJC2OA ,TABJC2OB ,TABJMHP  ,    &
                   TABJCOF2,TABJCH3BR,TABJF12B1,TABJF13B1,TABJCOFCL,    &
                   TABJCO2, TABJCOS
!
        CLOSE(NTAB)
        init = .FALSE.

        IF (.FALSE.) THEN

!       Write out photolysis tables in ASCII format for test purposes.
          OPEN(NTAB,FILE='jtable.asc',FORM='FORMATTED',STATUS='NEW')
!
          WRITE(NTAB,'(A)') 'TABPRES'
          WRITE(NTAB,'(10(E11.3))') TABPRES
          WRITE(NTAB,'(A)') 'TABANG'
          WRITE(NTAB,'(10(E11.3))') TABANG
          WRITE(NTAB,'(A)') 'TABT'
          WRITE(NTAB,'(10(E11.3))') TABT
          WRITE(NTAB,'(A)') 'TABO3'
          WRITE(NTAB,'(10(E11.3))') TABO3
          WRITE(NTAB,'(A)') 'SAO3C'
          WRITE(NTAB,'(10(E11.3))') SAO3C
          WRITE(NTAB,'(A)') 'TABJ2A'
          WRITE(NTAB,'(10(E11.3))') TABJ2A
          WRITE(NTAB,'(A)') 'TABJ2B'
          WRITE(NTAB,'(10(E11.3))') TABJ2B
          WRITE(NTAB,'(A)') 'TABJ3'
          WRITE(NTAB,'(10(E11.3))') TABJ3
          WRITE(NTAB,'(A)') 'TABJCH4'
          WRITE(NTAB,'(10(E11.3))') TABJCH4
          WRITE(NTAB,'(A)') 'TABJ3A'
          WRITE(NTAB,'(10(E11.3))') TABJ3A
          WRITE(NTAB,'(A)') 'TABCCL4'
          WRITE(NTAB,'(10(E11.3))') TABJCCL4
          WRITE(NTAB,'(A)') 'TABCl2O2'
          WRITE(NTAB,'(10(E11.3))') TABJCL2O2
          WRITE(NTAB,'(A)') 'TABCNITA'
          WRITE(NTAB,'(10(E11.3))') TABJCNITA
          WRITE(NTAB,'(A)') 'TABCNITB'
          WRITE(NTAB,'(10(E11.3))') TABJCNITB
          WRITE(NTAB,'(A)') 'TABJF11'
          WRITE(NTAB,'(10(E11.3))') TABJF11
          WRITE(NTAB,'(A)') 'TABJF113'
          WRITE(NTAB,'(10(E11.3))') TABJF113
          WRITE(NTAB,'(A)') 'TABJF12'
          WRITE(NTAB,'(10(E11.3))') TABJF12
          WRITE(NTAB,'(A)') 'TABJF22'
          WRITE(NTAB,'(10(E11.3))') TABJF22
          WRITE(NTAB,'(A)') 'TABJH2O'
          WRITE(NTAB,'(10(E11.3))') TABJH2O
          WRITE(NTAB,'(A)') 'TABJH2O2'
          WRITE(NTAB,'(10(E11.3))') TABJH2O2
          WRITE(NTAB,'(A)') 'TABHCL'
          WRITE(NTAB,'(10(E11.3))') TABJHCL
          WRITE(NTAB,'(A)') 'TABHNO3'
          WRITE(NTAB,'(10(E11.3))') TABJHNO3
          WRITE(NTAB,'(A)') 'TABJHOCL'
          WRITE(NTAB,'(10(E11.3))') TABJHOCL
          WRITE(NTAB,'(A)') 'TABJCH3CL'
          WRITE(NTAB,'(10(E11.3))') TABJCH3CL
          WRITE(NTAB,'(A)') 'TABJN2O'
          WRITE(NTAB,'(10(E11.3))') TABJN2O
          WRITE(NTAB,'(A)') 'TABJN2O5'
          WRITE(NTAB,'(10(E11.3))') TABJN2O5
          WRITE(NTAB,'(A)') 'TABJNO'
          WRITE(NTAB,'(10(E11.3))') TABJNO
          WRITE(NTAB,'(A)') 'TABJNO2'
          WRITE(NTAB,'(10(E11.3))') TABJNO2
          WRITE(NTAB,'(A)') 'TABJNO31'
          WRITE(NTAB,'(10(E11.3))') TABJNO31
          WRITE(NTAB,'(A)') 'TABJNO32'
          WRITE(NTAB,'(10(E11.3))') TABJNO32
          WRITE(NTAB,'(A)') 'TABJPNA'
          WRITE(NTAB,'(10(E11.3))') TABJPNA
          WRITE(NTAB,'(A)') 'TABJBRCL'
          WRITE(NTAB,'(10(E11.3))') TABJBRCL
          WRITE(NTAB,'(A)') 'TABJBRNO3'
          WRITE(NTAB,'(10(E11.3))') TABJBRNO3
          WRITE(NTAB,'(A)') 'TABJBRO'
          WRITE(NTAB,'(10(E11.3))') TABJBRO
          WRITE(NTAB,'(A)') 'TABJMCFM'
          WRITE(NTAB,'(10(E11.3))') TABJMCFM
          WRITE(NTAB,'(A)') 'TABJHOBR'
          WRITE(NTAB,'(10(E11.3))') TABJHOBR
          WRITE(NTAB,'(A)') 'TABJOCLO'
          WRITE(NTAB,'(10(E11.3))') TABJOCLO
          WRITE(NTAB,'(A)') 'TABJC2OA'
          WRITE(NTAB,'(10(E11.3))') TABJC2OA
          WRITE(NTAB,'(A)') 'TABJC2OB'
          WRITE(NTAB,'(10(E11.3))') TABJC2OB
          WRITE(NTAB,'(A)') 'TABJMHP'
          WRITE(NTAB,'(10(E11.3))') TABJMHP
          WRITE(NTAB,'(A)') 'TABCOF2'
          WRITE(NTAB,'(10(E11.3))') TABJCOF2
          WRITE(NTAB,'(A)') 'TABJCH3BR'
          WRITE(NTAB,'(10(E11.3))') TABJCH3BR
          WRITE(NTAB,'(A)') 'TABF12B1'
          WRITE(NTAB,'(10(E11.3))') TABJF12B1
          WRITE(NTAB,'(A)') 'TABF13B1'
          WRITE(NTAB,'(10(E11.3))') TABJF13B1
          WRITE(NTAB,'(A)') 'TABJCOFCL'
          WRITE(NTAB,'(10(E11.3))') TABJCOFCL
          WRITE(NTAB,'(A)') 'TABJCO2'
          WRITE(NTAB,'(10(E11.3))') TABJCO2
          WRITE(NTAB,'(A)') 'TABJCOS'
          WRITE(NTAB,'(10(E11.3))') TABJCOS
!
          CLOSE(NTAB)
        END IF
      END IF
!
!     Spacing of angles in lookup table
      COSINT=1.0/REAL(JPCHI-1-JPS90)
      RADINT=(SZAMAX - 90.0)*PI/(180.0*REAL(JPS90))
!
!     Model pressure levels in hPa. SZA in radians. T within range.
      DO IH=IHMIN,IHMAX
        P (IH)=PPA(IH)*0.01
        P (IH)=MIN(P (IH), TABPRES(    1))
        P (IH)=MAX(P (IH), TABPRES(JPLEV))
        TU(IH)=T(IH)
        TU(IH)=MIN(TU(IH), TMAX)
        TU(IH)=MAX(TU(IH), TMIN)
        ANGLE(IH)=ACOS(F(IH))
      END DO
!
!     Find the location of the pressure levels,
!     zenith angle, temperature and O3 in the photolysis look up tables.
!
!     i)  X. Pressure levels equally spaced in log(p)
!     ii) Y. Zenith angle
      DLNP=(LOG(TABPRES(1))-LOG(TABPRES(JPLEV)))/REAL(JPLEV-1)

      DO IH = IHMIN, IHMAX
        JX(IH)=INT((LOG(TABPRES(1))-LOG(P(IH)))/DLNP) + 1
!
        IF(ANGLE(IH) < TABANG(JPCHIN)) THEN
          JY(IH)=INT((1.0 - F(IH))/COSINT) + 1
        ELSE
          JY(IH)=JPCHIN + INT((ANGLE(IH)-0.5*PI)/RADINT)
        END IF
      END DO
!
!     iii) Z  temperature evenly spaced
      IF (JPTEM == 1) THEN
        DO IH = IHMIN, IHMAX
          JZ(IH)=1
        END DO
      ELSE
        TDIF=(TMAX-TMIN)/REAL(JPTEM-1)
        DO IH = IHMIN, IHMAX
          JZ(IH)=INT((TU(IH)-TMIN)/TDIF) + 1
        END DO
      END IF
!
!     iv) O3 factor
      IF (JPO3P == 1) THEN
        DO IH = IHMIN, IHMAX
          JO(IH)=1
        END DO
      ELSE
!       Find O3 factor using tabulated O3 column
        O3DIF=(O3MAX-O3MIN)/REAL(JPO3P-1)
        DO IH=IHMIN,IHMAX
          CSUP=(P(IH)            -TABPRES(JX(IH)))/                     &
               (TABPRES(JX(IH)+1)-TABPRES(JX(IH)))
!
!         Ratio O3 to standard atmosphere O3 above model pressure.
          O3FAC(IH)=O3COL(IH)/                                          &
                    ((1.0-CSUP)*SAO3C(JX(IH))+CSUP*SAO3C(JX(IH)+1))
!
!         Check that O3 column ratio is within range
          O3FAC(IH)=MIN(O3FAC(IH), O3MAX)
          O3FAC(IH)=MAX(O3FAC(IH), O3MIN)
!
          JO(IH)=INT((O3FAC(IH)-O3MIN)/O3DIF) + 1
        END DO
      END IF
!
!     If out of range set pointers accordingly.
      DO IH=IHMIN,IHMAX
        JX(IH) = MAX0(JX(IH),      1)
        JX(IH) = MIN0(JX(IH),JPLEV-1)
        JY(IH) = MAX0(JY(IH),      1)
        JY(IH) = MIN0(JY(IH),JPCHI-1)
        JZ(IH) = MAX0(JZ(IH),      1)
        JZ(IH) = MIN0(JZ(IH),JPTEM-1)
        JO(IH) = MAX0(JO(IH),      1)
        JO(IH) = MIN0(JO(IH),JPO3P-1)
!
        JZ(IH) = MAX0(JZ(IH),1)
        JO(IH) = MAX0(JO(IH),1)
!
        JXP1(IH) =      JX(IH) + 1
        JYP1(IH) =      JY(IH) + 1
        JZP1(IH) = MIN0(JZ(IH) + 1, JPTEM)
        JOP1(IH) = MIN0(JO(IH) + 1, JPO3P)
      END DO
!
!     Find points used in interpolation
      DO IH=IHMIN,IHMAX
!        i) X: pressure
         ZX(IH)=(LOG(P      (IH))       - LOG(TABPRES(JX(IH))))/       &
                (LOG(TABPRES(JXP1(IH))) - LOG(TABPRES(JX(IH))))
!
!        ii) Y: zenith angle
         ZY(IH)=(ANGLE(IH)         - TABANG (JY(IH)))/                  &
                (TABANG (JYP1(IH)) - TABANG (JY(IH)))
!
!        iii) Z: temperature
         IF (JPTEM == 1) THEN
           ZZ(IH)=1.0
         ELSE
           ZZ(IH)=(TU   (IH)       - TABT (JZ(IH)))/                    &
                  (TABT (JZP1(IH)) - TABT (JZ(IH)))
         END IF
!
!        iv) O: O3 profile
         IF (JPO3P == 1) THEN
           ZO(IH)=1.0
         ELSE
           ZO(IH)=(O3FAC(IH)       - TABO3(JO(IH)))/                    &
                  (TABO3(JOP1(IH)) - TABO3(JO(IH)))
         END IF
!
         CX(IH)=1.0 - ZX(IH)
         CY(IH)=1.0 - ZY(IH)
         CZ(IH)=1.0 - ZZ(IH)
         CO(IH)=1.0 - ZO(IH)
      END DO
!
!     Calculate weights for interpolation
      DO IH=IHMIN,IHMAX
        W1 (IH)=CX(IH)*CY(IH)*CZ(IH)*CO(IH)
        W2 (IH)=ZX(IH)*CY(IH)*CZ(IH)*CO(IH)
        W3 (IH)=ZX(IH)*ZY(IH)*CZ(IH)*CO(IH)
        W4 (IH)=CX(IH)*ZY(IH)*CZ(IH)*CO(IH)
        W5 (IH)=CX(IH)*CY(IH)*ZZ(IH)*CO(IH)
        W6 (IH)=ZX(IH)*CY(IH)*ZZ(IH)*CO(IH)
        W7 (IH)=ZX(IH)*ZY(IH)*ZZ(IH)*CO(IH)
        W8 (IH)=CX(IH)*ZY(IH)*ZZ(IH)*CO(IH)
        W9 (IH)=CX(IH)*CY(IH)*CZ(IH)*ZO(IH)
        W10(IH)=ZX(IH)*CY(IH)*CZ(IH)*ZO(IH)
        W11(IH)=ZX(IH)*ZY(IH)*CZ(IH)*ZO(IH)
        W12(IH)=CX(IH)*ZY(IH)*CZ(IH)*ZO(IH)
        W13(IH)=CX(IH)*CY(IH)*ZZ(IH)*ZO(IH)
        W14(IH)=ZX(IH)*CY(IH)*ZZ(IH)*ZO(IH)
        W15(IH)=ZX(IH)*ZY(IH)*ZZ(IH)*ZO(IH)
        W16(IH)=CX(IH)*ZY(IH)*ZZ(IH)*ZO(IH)
      END DO
!
!     Look up photolysis rates. 38 rates.
      DO IH=IHMIN,IHMAX
!
!     If angle in range
      IF(ANGLE(IH) <= TABANG(JPCHI)) THEN
!1
        AJ2A  (IH)=                                                     &
           W1 (IH)*TABJ2A   (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJ2A   (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJ2A   (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJ2A   (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJ2A   (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJ2A   (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJ2A   (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJ2A   (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJ2A   (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJ2A   (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJ2A   (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJ2A   (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJ2A   (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJ2A   (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJ2A   (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJ2A   (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!1B
        AJ2B  (IH)=                                                     &
           W1 (IH)*TABJ2B   (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJ2B   (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJ2B   (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJ2B   (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJ2B   (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJ2B   (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJ2B   (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJ2B   (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJ2B   (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJ2B   (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJ2B   (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJ2B   (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJ2B   (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJ2B   (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJ2B   (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJ2B   (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!2
        AJ3   (IH)=                                                     &
           W1 (IH)*TABJ3    (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJ3    (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJ3    (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJ3    (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJ3    (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJ3    (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJ3    (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJ3    (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJ3    (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJ3    (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJ3    (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJ3    (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJ3    (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJ3    (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJ3    (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJ3    (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!3
        AJ3A  (IH)=                                                     &
           W1 (IH)*TABJ3A   (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJ3A   (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJ3A   (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJ3A   (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJ3A   (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJ3A   (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJ3A   (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJ3A   (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJ3A   (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJ3A   (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJ3A   (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJ3A   (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJ3A   (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJ3A   (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJ3A   (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJ3A   (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!4
        AJNO  (IH)=                                                     &
           W1 (IH)*TABJNO   (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJNO   (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJNO   (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJNO   (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJNO   (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJNO   (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJNO   (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJNO   (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJNO   (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJNO   (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJNO   (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJNO   (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJNO   (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJNO   (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJNO   (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJNO   (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!5
        AJNO2 (IH)=                                                     &
           W1 (IH)*TABJNO2  (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJNO2  (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJNO2  (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJNO2  (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJNO2  (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJNO2  (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJNO2  (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJNO2  (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJNO2  (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJNO2  (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJNO2  (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJNO2  (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJNO2  (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJNO2  (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJNO2  (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJNO2  (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!6
        AJNO31(IH)=                                                     &
           W1 (IH)*TABJNO31 (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJNO31 (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJNO31 (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJNO31 (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJNO31 (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJNO31 (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJNO31 (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJNO31 (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJNO31 (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJNO31 (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJNO31 (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJNO31 (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJNO31 (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJNO31 (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJNO31 (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJNO31 (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!7
        AJNO32(IH)=                                                     &
           W1 (IH)*TABJNO32 (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJNO32 (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJNO32 (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJNO32 (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJNO32 (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJNO32 (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJNO32 (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJNO32 (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJNO32 (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJNO32 (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJNO32 (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJNO32 (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJNO32 (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJNO32 (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJNO32 (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJNO32 (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!8
        AJN2O (IH)=                                                     &
           W1 (IH)*TABJN2O  (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJN2O  (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJN2O  (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJN2O  (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJN2O  (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJN2O  (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJN2O  (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJN2O  (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJN2O  (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJN2O  (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJN2O  (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJN2O  (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJN2O  (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJN2O  (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJN2O  (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJN2O  (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!9
        AJN2O5(IH)=                                                     &
           W1 (IH)*TABJN2O5 (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJN2O5 (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJN2O5 (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJN2O5 (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJN2O5 (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJN2O5 (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJN2O5 (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJN2O5 (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJN2O5 (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJN2O5 (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJN2O5 (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJN2O5 (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJN2O5 (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJN2O5 (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJN2O5 (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJN2O5 (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!10
        AJHNO3(IH)=                                                     &
           W1 (IH)*TABJHNO3 (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJHNO3 (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJHNO3 (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJHNO3 (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJHNO3 (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJHNO3 (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJHNO3 (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJHNO3 (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJHNO3 (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJHNO3 (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJHNO3 (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJHNO3 (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJHNO3 (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJHNO3 (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJHNO3 (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJHNO3 (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!11
        AJCNITA(IH)=                                                    &
           W1 (IH)*TABJCNITA(JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJCNITA(JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJCNITA(JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJCNITA(JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJCNITA(JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJCNITA(JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJCNITA(JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJCNITA(JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJCNITA(JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJCNITA(JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJCNITA(JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJCNITA(JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJCNITA(JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJCNITA(JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJCNITA(JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJCNITA(JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
        AJCNITB(IH)=                                                    &
           W1 (IH)*TABJCNITB(JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJCNITB(JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJCNITB(JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJCNITB(JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJCNITB(JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJCNITB(JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJCNITB(JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJCNITB(JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJCNITB(JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJCNITB(JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJCNITB(JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJCNITB(JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJCNITB(JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJCNITB(JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJCNITB(JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJCNITB(JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!12
        AJPNA (IH)=                                                     &
           W1 (IH)*TABJPNA  (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJPNA  (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJPNA  (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJPNA  (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJPNA  (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJPNA  (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJPNA  (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJPNA  (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJPNA  (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJPNA  (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJPNA  (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJPNA  (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJPNA  (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJPNA  (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJPNA  (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJPNA  (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!13
        AJH2O2(IH)=                                                     &
           W1 (IH)*TABJH2O2 (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJH2O2 (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJH2O2 (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJH2O2 (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJH2O2 (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJH2O2 (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJH2O2 (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJH2O2 (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJH2O2 (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJH2O2 (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJH2O2 (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJH2O2 (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJH2O2 (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJH2O2 (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJH2O2 (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJH2O2 (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!14
        AJH2O (IH)=                                                     &
           W1 (IH)*TABJH2O  (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJH2O  (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJH2O  (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJH2O  (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJH2O  (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJH2O  (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJH2O  (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJH2O  (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJH2O  (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJH2O  (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJH2O  (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJH2O  (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJH2O  (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJH2O  (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJH2O  (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJH2O  (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!15
        AJHOCL(IH)=                                                     &
           W1 (IH)*TABJHOCL (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJHOCL (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJHOCL (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJHOCL (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJHOCL (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJHOCL (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJHOCL (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJHOCL (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJHOCL (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJHOCL (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJHOCL (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJHOCL (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJHOCL (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJHOCL (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJHOCL (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJHOCL (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!16
        AJHCL (IH)=                                                     &
           W1 (IH)*TABJHCL  (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJHCL  (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJHCL  (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJHCL  (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJHCL  (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJHCL  (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJHCL  (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJHCL  (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJHCL  (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJHCL  (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJHCL  (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJHCL  (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJHCL  (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJHCL  (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJHCL  (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJHCL  (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!17
        AJCL2O2(IH)=                                                    &
           W1 (IH)*TABJCL2O2(JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJCL2O2(JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJCL2O2(JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJCL2O2(JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJCL2O2(JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJCL2O2(JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJCL2O2(JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJCL2O2(JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJCL2O2(JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJCL2O2(JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJCL2O2(JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJCL2O2(JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJCL2O2(JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJCL2O2(JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJCL2O2(JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJCL2O2(JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!18
        AJBRO (IH)=                                                     &
           W1 (IH)*TABJBRO  (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJBRO  (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJBRO  (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJBRO  (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJBRO  (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJBRO  (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJBRO  (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJBRO  (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJBRO  (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJBRO  (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJBRO  (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJBRO  (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJBRO  (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJBRO  (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJBRO  (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJBRO  (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!19
        AJBRCL(IH)=                                                     &
           W1 (IH)*TABJBRCL (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJBRCL (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJBRCL (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJBRCL (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJBRCL (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJBRCL (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJBRCL (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJBRCL (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJBRCL (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJBRCL (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJBRCL (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJBRCL (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJBRCL (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJBRCL (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJBRCL (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJBRCL (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!20
        AJCO2(IH)=                                                      &
           W1 (IH)*TABJCO2  (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJCO2  (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJCO2  (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJCO2  (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJCO2  (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJCO2  (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJCO2  (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJCO2  (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJCO2  (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJCO2  (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJCO2  (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJCO2  (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJCO2  (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJCO2  (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJCO2  (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJCO2  (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
! COS
        AJCOS(IH)=                                                      &
           W1 (IH)*TABJCOS  (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJCOS  (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJCOS  (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJCOS  (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJCOS  (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJCOS  (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJCOS  (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJCOS  (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJCOS  (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJCOS  (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJCOS  (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJCOS  (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJCOS  (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJCOS  (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJCOS  (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJCOS  (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!
!     Else darkness
      ELSE
        AJ2A   (IH) = 0.0
        AJ2B   (IH) = 0.0
        AJ3    (IH) = 0.0
        AJ3A   (IH) = 0.0
        AJNO   (IH) = 0.0
        AJNO2  (IH) = 0.0
        AJNO31 (IH) = 0.0
        AJNO32 (IH) = 0.0
        AJN2O  (IH) = 0.0
        AJN2O5 (IH) = 0.0
        AJHNO3 (IH) = 0.0
        AJCNITA(IH) = 0.0
        AJCNITB(IH) = 0.0
        AJPNA  (IH) = 0.0
        AJH2O2 (IH) = 0.0
        AJH2O  (IH) = 0.0
        AJHOCL (IH) = 0.0
        AJHCL  (IH) = 0.0
        AJCL2O2(IH) = 0.0
        AJBRO  (IH) = 0.0
        AJBRCL (IH) = 0.0
        AJCO2  (IH) = 0.0
        AJCOS  (IH) = 0.0
      END IF
!
      END DO
!
      DO IH=IHMIN,IHMAX
!
!     If angle in range
      IF(ANGLE(IH) <= TABANG(JPCHI)) THEN
!20
        AJBRNO3(IH)=                                                    &
           W1 (IH)*TABJBRNO3(JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJBRNO3(JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJBRNO3(JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJBRNO3(JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJBRNO3(JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJBRNO3(JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJBRNO3(JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJBRNO3(JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJBRNO3(JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJBRNO3(JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJBRNO3(JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJBRNO3(JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJBRNO3(JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJBRNO3(JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJBRNO3(JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJBRNO3(JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!21
        AJHOBR(IH)=                                                     &
           W1 (IH)*TABJHOBR (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJHOBR (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJHOBR (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJHOBR (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJHOBR (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJHOBR (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJHOBR (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJHOBR (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJHOBR (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJHOBR (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJHOBR (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJHOBR (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJHOBR (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJHOBR (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJHOBR (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJHOBR (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!22
        AJOCLO(IH)=                                                     &
           W1 (IH)*TABJOCLO (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJOCLO (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJOCLO (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJOCLO (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJOCLO (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJOCLO (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJOCLO (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJOCLO (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJOCLO (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJOCLO (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJOCLO (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJOCLO (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJOCLO (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJOCLO (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJOCLO (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJOCLO (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!23
        AJC2OA(IH)=                                                     &
           W1 (IH)*TABJC2OA (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJC2OA (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJC2OA (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJC2OA (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJC2OA (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJC2OA (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJC2OA (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJC2OA (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJC2OA (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJC2OA (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJC2OA (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJC2OA (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJC2OA (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJC2OA (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJC2OA (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJC2OA (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!24
        AJC2OB(IH)=                                                     &
           W1 (IH)*TABJC2OB (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJC2OB (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJC2OB (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJC2OB (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJC2OB (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJC2OB (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJC2OB (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJC2OB (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJC2OB (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJC2OB (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJC2OB (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJC2OB (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJC2OB (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJC2OB (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJC2OB (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJC2OB (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!25
        AJMHP (IH)=                                                     &
           W1 (IH)*TABJMHP  (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJMHP  (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJMHP  (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJMHP  (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJMHP  (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJMHP  (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJMHP  (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJMHP  (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJMHP  (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJMHP  (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJMHP  (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJMHP  (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJMHP  (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJMHP  (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJMHP  (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJMHP  (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!26
        AJCH3CL(IH)=                                                    &
           W1 (IH)*TABJCH3CL(JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJCH3CL(JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJCH3CL(JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJCH3CL(JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJCH3CL(JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJCH3CL(JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJCH3CL(JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJCH3CL(JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJCH3CL(JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJCH3CL(JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJCH3CL(JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJCH3CL(JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJCH3CL(JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJCH3CL(JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJCH3CL(JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJCH3CL(JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!27
        AJMCFM(IH)=                                                     &
           W1 (IH)*TABJMCFM (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJMCFM (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJMCFM (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJMCFM (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJMCFM (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJMCFM (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJMCFM (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJMCFM (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJMCFM (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJMCFM (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJMCFM (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJMCFM (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJMCFM (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJMCFM (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJMCFM (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJMCFM (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!28
        AJF11 (IH)=                                                     &
           W1 (IH)*TABJF11  (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJF11  (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJF11  (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJF11  (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJF11  (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJF11  (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJF11  (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJF11  (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJF11  (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJF11  (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJF11  (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJF11  (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJF11  (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJF11  (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJF11  (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJF11  (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!29
        AJF12 (IH)=                                                     &
           W1 (IH)*TABJF12  (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJF12  (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJF12  (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJF12  (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJF12  (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJF12  (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJF12  (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJF12  (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJF12  (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJF12  (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJF12  (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJF12  (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJF12  (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJF12  (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJF12  (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJF12  (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!30
        AJF22 (IH)=                                                     &
           W1 (IH)*TABJF22  (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJF22  (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJF22  (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJF22  (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJF22  (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJF22  (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJF22  (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJF22  (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJF22  (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJF22  (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJF22  (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJF22  (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJF22  (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJF22  (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJF22  (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJF22  (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!31
        AJF113(IH)=                                                     &
           W1 (IH)*TABJF113 (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJF113 (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJF113 (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJF113 (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJF113 (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJF113 (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJF113 (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJF113 (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJF113 (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJF113 (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJF113 (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJF113 (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJF113 (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJF113 (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJF113 (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJF113 (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!32
        AJCCL4(IH)=                                                     &
           W1 (IH)*TABJCCL4 (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJCCL4 (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJCCL4 (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJCCL4 (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJCCL4 (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJCCL4 (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJCCL4 (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJCCL4 (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJCCL4 (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJCCL4 (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJCCL4 (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJCCL4 (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJCCL4 (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJCCL4 (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJCCL4 (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJCCL4 (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!33
        AJCH3BR(IH)=                                                    &
           W1 (IH)*TABJCH3BR(JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJCH3BR(JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJCH3BR(JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJCH3BR(JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJCH3BR(JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJCH3BR(JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJCH3BR(JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJCH3BR(JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJCH3BR(JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJCH3BR(JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJCH3BR(JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJCH3BR(JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJCH3BR(JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJCH3BR(JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJCH3BR(JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJCH3BR(JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!34
        AJF12B1(IH)=                                                    &
           W1 (IH)*TABJF12B1(JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJF12B1(JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJF12B1(JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJF12B1(JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJF12B1(JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJF12B1(JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJF12B1(JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJF12B1(JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJF12B1(JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJF12B1(JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJF12B1(JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJF12B1(JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJF12B1(JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJF12B1(JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJF12B1(JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJF12B1(JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!35
        AJF13B1(IH)=                                                    &
           W1 (IH)*TABJF13B1(JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJF13B1(JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJF13B1(JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJF13B1(JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJF13B1(JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJF13B1(JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJF13B1(JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJF13B1(JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJF13B1(JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJF13B1(JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJF13B1(JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJF13B1(JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJF13B1(JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJF13B1(JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJF13B1(JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJF13B1(JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!36
        AJCOF2(IH)=                                                     &
           W1 (IH)*TABJCOF2 (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJCOF2 (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJCOF2 (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJCOF2 (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJCOF2 (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJCOF2 (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJCOF2 (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJCOF2 (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJCOF2 (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJCOF2 (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJCOF2 (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJCOF2 (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJCOF2 (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJCOF2 (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJCOF2 (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJCOF2 (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!37
        AJCOFCL(IH)=                                                    &
           W1 (IH)*TABJCOFCL(JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJCOFCL(JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJCOFCL(JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJCOFCL(JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJCOFCL(JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJCOFCL(JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJCOFCL(JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJCOFCL(JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJCOFCL(JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJCOFCL(JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJCOFCL(JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJCOFCL(JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJCOFCL(JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJCOFCL(JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJCOFCL(JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJCOFCL(JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!38
        AJCH4 (IH)=                                                     &
           W1 (IH)*TABJCH4  (JX  (IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W2 (IH)*TABJCH4  (JXP1(IH),JY  (IH),JZ  (IH),JO  (IH))       &
         + W3 (IH)*TABJCH4  (JXP1(IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W4 (IH)*TABJCH4  (JX  (IH),JYP1(IH),JZ  (IH),JO  (IH))       &
         + W5 (IH)*TABJCH4  (JX  (IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W6 (IH)*TABJCH4  (JXP1(IH),JY  (IH),JZP1(IH),JO  (IH))       &
         + W7 (IH)*TABJCH4  (JXP1(IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W8 (IH)*TABJCH4  (JX  (IH),JYP1(IH),JZP1(IH),JO  (IH))       &
         + W9 (IH)*TABJCH4  (JX  (IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W10(IH)*TABJCH4  (JXP1(IH),JY  (IH),JZ  (IH),JOP1(IH))       &
         + W11(IH)*TABJCH4  (JXP1(IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W12(IH)*TABJCH4  (JX  (IH),JYP1(IH),JZ  (IH),JOP1(IH))       &
         + W13(IH)*TABJCH4  (JX  (IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W14(IH)*TABJCH4  (JXP1(IH),JY  (IH),JZP1(IH),JOP1(IH))       &
         + W15(IH)*TABJCH4  (JXP1(IH),JYP1(IH),JZP1(IH),JOP1(IH))       &
         + W16(IH)*TABJCH4  (JX  (IH),JYP1(IH),JZP1(IH),JOP1(IH))
!
!     Else darkness
      ELSE
        AJBRNO3(IH) = 0.0
        AJHOBR (IH) = 0.0
        AJOCLO (IH) = 0.0
        AJC2OA (IH) = 0.0
        AJC2OB (IH) = 0.0
        AJMHP  (IH) = 0.0
        AJCH3CL(IH) = 0.0
        AJMCFM (IH) = 0.0
        AJF11  (IH) = 0.0
        AJF12  (IH) = 0.0
        AJF22  (IH) = 0.0
        AJF113 (IH) = 0.0
        AJCCL4 (IH) = 0.0
        AJCH3BR(IH) = 0.0
        AJF12B1(IH) = 0.0
        AJF13B1(IH) = 0.0
        AJCOF2 (IH) = 0.0
        AJCOFCL(IH) = 0.0
        AJCH4  (IH) = 0.0
      END IF
!
      END DO
!
      END SUBROUTINE CALCJS

!  Initialization of stratospheric photolysis tables.
!
!  Test version
! ######################################################################
!
! Subroutine Interface:
!
      SUBROUTINE INIJTAB(mype)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                     C
!     SUBROUTINE TO INITIALISE PHOTOLYSIS TABLES                      C
!                                                                     C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! 3/1/2006 Modified to include CO2 photolysis
!                                     Olaf Morgenstern
!
!      AJ2A:     O2       +  hv             ->   O     +   O
!      AJ2B:     O2       +  hv             ->   O     +   O(1D)
!      AJ3:      O3       +  hv             ->   O2    +   O(3P)
!      AJ3A:     O3       +  hv(<310nm)     ->   O2    +   O(1D)
!      AJNO:     NO       +  hv             ->   N     +   O
!      AJNO2:    NO2      +  hv             ->   NO    +   O(3P)
!      AJNO31:   NO3      +  hv             ->   NO    +   O2
!      AJNO32:   NO3      +  hv             ->   NO2   +   O
!      AJHNO3:   HNO3     +  hv             ->   NO2   +   OH
!      AJPNA:    HO2NO2   +  hv             ->   NO2   +   HO2 / NO3 + OH
!      AJN2O:    N2O      +  hv             ->   N2    +   O(1D)
!      AJN2O5:   N2O5     +  hv             ->   NO2   +   NO3
!      AJH2O:    H2O      +  hv             ->   OH    +   H
!      AJH2O2:   H2O2     +  hv             ->   OH    +   OH
!      AJCNITA:  ClONO2   +  hv             ->   Cl    +   NO3
!      AJCNITB:  ClONO2   +  hv             ->   ClO   +   NO2
!      AJF11:    CFCl3    +  hv             ->   3Cl
!      AJF12:    CF2Cl2   +  hv             ->   2Cl
!      AJF22:    CHF2Cl   +  hv             ->   Cl
!      AJF113:   CF2CFCl3 +  hv             ->   3Cl
!      AJCH3CL:  CH3Cl    +  hv             ->   CH3   +   Cl
!      AJCCL4:   CCl4     +  hv             ->   4Cl
!      AJHCL:    HCl      +  hv             ->   H     +   Cl
!      AJHOCL:   HOCl     +  hv             ->   OH    +   Cl
!      AJOCLO:   OClO     +  hv             ->   O     +   ClO
!      AJCL2O2:  Cl2O2    +  hv             ->   O2    +   2Cl
!      AJBRO:    BrO      +  hv             ->   Br    +   O
!      AJHOBR:   HOBr     +  hv             ->   Br    +   OH
!      AJBRNO3:  BrNO3    +  hv             ->   Br    +   NO3 / Bro + NO2
!      AJBRCL:   BrCl     +  hv             ->   Br    +   Cl
!      AJC2OA:   CH2O     +  hv             ->   H     +   CHO
!      AJC2OB:   CH2O     +  hv             ->   H2    +   CO
!      AJMHP :   CH3OOH   +  hv             ->   CH3O  +   OH
!      AJCH3BR:  CH3Br    +  hv             ->   CH3   +   Br
!      AJMCFM:   CH3CCl3  +  hv             ->   CHx   +   3Cl
!      AJCH4:    CH4      +  hv             ->   CH3   +   H
!      AJF12B1:  CBrClF2  +  hv             ->   Cxx   +   Cl   + Br
!      AJF13B1:  CBrF3    +  hv             ->   Cxx   +   Br
!      AJCOF2:   COF2     +  hv             ->   Cx    +   2HF
!      AJCOFCL:  COFCl    +  hv             ->   Cx    +   Cl + HF
!      AJCO2:    CO2      +  hv             ->   CO    +   O(3P)
!      AJCOS:    COS      +  hv             ->   CO    +   S
!
!     Name      Type              Description.
!
!     ALT      Array of real      Altitude at the edge of each level
!                                 in km.
!
!     ALTC     Array of real      Altitude at the centre of each level
!                                 in km.
!
!     SAO2     Array of real      Vertical O2 column above the edge of
!                                 each level in molecules per cm^2.
!
!     SAO3C    Array of real      Vertical O3 column above the centre
!                                 of each level in molecules per cm^2.
!
!     DO2C     Array of real      Average number density of O2 within
!                                 each level in molecules per cm^3.
!
!     DO3C     Array of real      Average number density of O3 within
!                                 each level in molecules per cm^3.
!
!     DO3CU                       Scaled O3 profiles
!
!     DRSC                        [M] within level
!
!     DALT                        Thickness of level (km)
!
!     TEMPC    Array of real      Temperature (K) at the centre of each level.
!
!     PRESC    Array of real      Pressure (mb) at the centre of each level.
!
!     TABPRES  Array of real      Pressure (mb) of each level in the
!                                 look up table.
!
!     PRES     Array of real      Pressure (mb) at the edge of each level.
!
! Current code owner: Olaf Morgenstern.
! Version history:
! v1.1    11/2003  Original code adapted to UM. Introduce fill_spectra
!                  subroutine.     Olaf Morgenstern
! v1.2 20/12/2005  Introduce switch to read standard temperature and
!                  ozone profiles from string common instead of from
!                  external file.  Olaf Morgenstern
! v1.3 05/01/2006  Introduce CO2 photolysis  Olaf Morgenstern
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT NONE

#include "parpho.h"
#include "tabjs.h"
#include "crossec.h"
#include "stdto3.h"

!
! Subroutine interface
      INTEGER, INTENT(IN) :: mype  ! processor number

! Local variables
! do not read STDTO3 input file. Instead, use hard-wired values.
      LOGICAL, PARAMETER :: read_stdto3 = .FALSE.
!
      REAL, PARAMETER :: GG=9.81,RA=287.06,COLFAC=2.132E20
!     O2 mixing ratio
      REAL, PARAMETER :: R2O2=0.209
!     Ground albedo
      REAL, PARAMETER :: ALBEDO=0.30
!
!     Whether scattering is on or not.
      LOGICAL, PARAMETER :: LSCAT =.TRUE.
!
      INTEGER :: JC
      INTEGER :: JL
      INTEGER :: JLU
      INTEGER :: JLL
      INTEGER :: JO
      INTEGER :: JT
      INTEGER :: JW
      INTEGER :: L

      REAL :: O3ABOVE
!
      REAL :: CINF
      REAL :: CSUP
      REAL :: DPDIF
      REAL :: PDIF
      REAL :: QU

      REAL :: FACTOR(JPWAV)
      REAL :: DO2C  (JPLEV)
      REAL :: DO3C  (JPLEV)
      REAL :: DRSC  (JPLEV)
      REAL :: DO3CU (JPLEV)
      REAL :: PRES  (JPLEVP1)
      REAL :: PRESC (JPLEV)
      REAL :: O3VMRC(JPLEV)
      REAL :: SAO2  (JPLEVP1)
      REAL :: SAO3C (JPLEV)
!
      REAL :: TABS  (JPLEV,JPCHI,JPWAV)
      REAL :: TEMPC (JPLEV)
      REAL :: TSPO2 (JPLEV,JPCHI)
      REAL :: ZENMAX(JPLEV)
      REAL :: ALT   (JPLEVP1)
      REAL :: ALTC  (JPLEV)
      REAL :: DALT  (JPLEV)
!
      REAL :: RDTEMP(NIN)
      REAL :: RDPRES(NIN)
      REAL :: RDO3V (NIN)
!
!
!
      IF (read_stdto3) THEN
!     Read in standard atmosphere temperature and O3
        OPEN(NSTD,file='STDTO3')
        REWIND(NSTD)
        DO L=1,NIN
          READ(NSTD,*) RDPRES(L),RDTEMP(L),RDO3V(L)
        END DO
        CLOSE(NSTD)
      ELSE
! Read standard atmosphere data from string variable
        DO L=1,NIN
          READ(stdto3(l),*) RDPRES(L),RDTEMP(L),RDO3V(L)
        END DO
      END IF

!     Set up j table inter levels equally spaced in log(p)
      PRES(      1)=RDPRES(  1)
      PRES(JPLEVP1)=RDPRES(NIN)
      PDIF=LOG(RDPRES(1)) - LOG(RDPRES(NIN))
      DPDIF=PDIF/(REAL(JPLEVP1-1))
      DO L=2,JPLEVP1-1
        PRES (L)=EXP(LOG(PRES(L-1)) - DPDIF)
      END DO
!     Pressure at jtable mid levels
      DO L=1,JPLEV
        PRESC  (L)=EXP(0.5*(LOG(PRES(L)) + LOG(PRES(L+1))))
        TABPRES(L)=PRESC(L)
      END DO
!
!     Interpolate temperature and O3 to level centres
      DO L=1,JPLEV
        DO JL=1,NIN
          IF(RDPRES(JL) < PRESC(L)) EXIT
        END DO
        JLU=JL
        JLL=JL-1
        CSUP=(LOG(RDPRES(JLL)) - LOG(PRESC (L  )))/                    &
             (LOG(RDPRES(JLL)) - LOG(RDPRES(JLU)))
        CINF=1.0 - CSUP
        TEMPC (L)=CSUP*RDTEMP(JLU) + CINF*RDTEMP(JLL)
        O3VMRC(L)=CSUP*RDO3V (JLU) + CINF*RDO3V (JLL)
      END DO
!
!     Work out geopotential height (km). Bottom level is ground.
      PRES(1)=1013.0
      ALT (1)=0.0
      DO L=2,JPLEVP1
        ALT(L)=ALT(L-1) +0.001*RA*TEMPC(L-1)*LOG(PRES(L-1)/PRES(L))/GG
      END DO
!     Height at centre of level and thickness of level (km).
      DO L=1,JPLEV
        ALTC(L)=0.5*(ALT(L+1) + ALT(L))
        DALT(L)=     ALT(L+1) - ALT(L)
      END DO
!
!     Work out O2 column above interfaces
      DO L=1,JPLEVP1
        SAO2(L)=R2O2*100.0*PRES(L)*COLFAC
      END DO
!
!     Average [O2] and [M] within a level
      DO L=1,JPLEV
        DO2C(L)=(SAO2(L)-SAO2(L+1))/(1.0E5*DALT(L))
        DRSC(L)= DO2C(L)/R2O2
      END DO
!     Add O2 above top level to top level
      DO2C(JPLEV)=DO2C(JPLEV) + SAO2(JPLEVP1)/(1.0E5*DALT(JPLEV))
!
!     Average [O3] within a level
      DO L=1,JPLEV
        DO3C(L)=O3VMRC(L)*DO2C(L)/R2O2
      END DO
!
!     O3 above top level with scale height of 3.5 km
      O3ABOVE=0.5*O3VMRC(JPLEV)*SAO2(JPLEVP1)/R2O2
      DO3C (JPLEV)=DO3C(JPLEV) + O3ABOVE/(1.0E5*DALT(JPLEV))
!
!     Column O3 above centre:
      SAO3C(JPLEV)=0.5*(O3ABOVE + 1.E5*DO3C(JPLEV)*DALT(JPLEV))
      DO L=JPLEV-1,1,-1
         SAO3C(L)=SAO3C(L+1) + 0.5E5*(DALT(L+1)*DO3C(L+1)               &
                                     +DALT(L  )*DO3C(L  ))
      END DO
!
! Call fill_spectra routine here to define spectral information etc.
      CALL fill_spectra
!
!     Calculate the Rayleigh scattering cross section.
      CALL SCATCS(JPWAV,SCS,WAVENM)
!
!     Set zenith angle, max za, temperature and O3 grid.
      CALL SETZEN(TABANG,ZENMAX,ALTC,TABT,TABO3)
!
!     Loop over O3 profiles
      DO JO=1,JPO3P
!
!     Scale O3 profile
      DO JL=1,JPLEV
        DO3CU(JL)=DO3C(JL)*TABO3(JO)
      END DO
!
      CALL SETTAB(ALT,ALTC,DO2C,DO3CU,TEMPC,                            &
           PRESC,DRSC,TABS,TABANG,ALBEDO,ZENMAX,LSCAT,                  &
           TSPO2,TABPRES,QUANTA,WAVENM,SCS,AO2,AO2SR,AO3,DALT)
!
!  Set up the photolysis rate tabulations. 39 rates.
!
      TABJ2A   (:,:,:,JO)=0.0
      TABJ2B   (:,:,:,JO)=0.0
      TABJ3    (:,:,:,JO)=0.0
      TABJ3A   (:,:,:,JO)=0.0
      TABJNO2  (:,:,:,JO)=0.0
      TABJNO   (:,:,:,JO)=0.0
      TABJNO31 (:,:,:,JO)=0.0
      TABJNO32 (:,:,:,JO)=0.0
      TABJN2O  (:,:,:,JO)=0.0
      TABJN2O5 (:,:,:,JO)=0.0
      TABJHNO3 (:,:,:,JO)=0.0
      TABJPNA  (:,:,:,JO)=0.0
      TABJCNITA(:,:,:,JO)=0.0
      TABJCNITB(:,:,:,JO)=0.0
      TABJH2O2 (:,:,:,JO)=0.0
      TABJH2O  (:,:,:,JO)=0.0
      TABJHOCL (:,:,:,JO)=0.0
      TABJHCL  (:,:,:,JO)=0.0
      TABJCL2O2(:,:,:,JO)=0.0
      TABJCH3CL(:,:,:,JO)=0.0
      TABJCCL4 (:,:,:,JO)=0.0
      TABJF11  (:,:,:,JO)=0.0
      TABJF12  (:,:,:,JO)=0.0
      TABJF22  (:,:,:,JO)=0.0
      TABJF113 (:,:,:,JO)=0.0
      TABJBRO  (:,:,:,JO)=0.0
      TABJHOBR (:,:,:,JO)=0.0
      TABJBRNO3(:,:,:,JO)=0.0
      TABJBRCL (:,:,:,JO)=0.0
      TABJOCLO (:,:,:,JO)=0.0
      TABJC2OA (:,:,:,JO)=0.0
      TABJC2OB (:,:,:,JO)=0.0
      TABJMHP  (:,:,:,JO)=0.0
      TABJMCFM (:,:,:,JO)=0.0
      TABJCH3BR(:,:,:,JO)=0.0
      TABJF12B1(:,:,:,JO)=0.0
      TABJF13B1(:,:,:,JO)=0.0
      TABJCOF2 (:,:,:,JO)=0.0
      TABJCOFCL(:,:,:,JO)=0.0
      TABJCH4  (:,:,:,JO)=0.0
      TABJCO2  (:,:,:,JO)=0.0
      TABJCOS  (:,:,:,JO)=0.0
!
!     Lyman alpha photolysis
      CALL LYMANA(TSPO2,TABJ2B,TABJH2O,TABJCH4,QUANTA,JO)
!
!     Temperature loop.
      DO JT = 1,JPTEM
!
!        Set up the cross-sections for this pressure & temperature.
!        N2O5
         CALL ACSN2O5(TABT(JT),JPWAV,WAVENM,AN2O5)
!
!        N2O
         CALL ACSN2O (TABT(JT),JPWAV,WAVENM,AN2O)
!
!        CCl4
         CALL ACSCCL4(TABT(JT),JPWAV,WAVENM,ACCL4)
!
!        F11 (CFCl3)
         CALL ACSF11 (TABT(JT),WAVENM,AF11)
!
!        F12 (CF2Cl2)
         CALL ACSF12 (TABT(JT),WAVENM,AF12)
!
!        F22 (CHF2Cl)
         CALL ACSF22 (TABT(JT),JPWAV,WAVENM,AF22)
!
!        CH3Cl
         CALL ACSMC  (TABT(JT),JPWAV,WAVENM,ACH3CL)
!
!        O[1D] Quantum yield.
         CALL QUANTO12(TABT(JT),JPWAV,WAVENM,QEO1D)
!
!        O3
         CALL ACSO3  (TABT(JT),JPWAV,AO3)
!
!        HNO3
         CALL ACSHNO3(TABT(JT),WAVENM,AHNO3)
!
!        ClONO2
         CALL ACSCNIT(TABT(JT),WAVENM,ACNITA, ACNITB)
!
!        H2O2
         CALL ACSH2O2(TABT(JT),JPWAV,WAVENM,AH2O2)
!
!        BrCl
         CALL ACSBRCL(TABT(JT),JPWAV,WAVENM,ABRCL)
!
!        NO2
         CALL ACSNO2 (TABT(JT),WAVENM,ANO2 )
!
!        COS
         CALL ACSCOS (TABT(JT),WAVENM,AOCS )
!
!        Level loop.
         DO JL = 1,JPLEV
!
!        Zenith angle loop.
         DO JC = 1 , JPCHI
!
!          Calculate the NO pressure dependent cross section.
           CALL ACSNO(TABANG(JC),PRESC(JL),TSPO2(JL,1),ANO)
!
!          Schumann-Runge bands.
           IF ( TABANG(JC) <= ZENMAX(JL) ) THEN
             CALL ACSSR(TSPO2(JL,JC),JPWAV,AO2SR)
           ELSE
             DO JW = JPLO, JPHI
               AO2SR(JW) = 0.0
             END DO
           END IF
!
!          Set up enhancement factor array.
           DO JW = JPLO , JPHI
              FACTOR(JW) = MAX(0.0,TABS(JL,JC,JW))
           END DO
!
!          Wavelength loop.
           DO JW = JPLO , JPHI
!
!           Number of photons into the volume element.
!           (need to attenuate quanta above top level?)
             QU = QUANTA(JW)*FACTOR(JW)

      TABJ2A (JL,JC,JT,JO)=TABJ2A(JL,JC,JT,JO)+QU*(AO2A(JW)+AO2SR(JW))
      TABJ2B (JL,JC,JT,JO)=TABJ2B(JL,JC,JT,JO)+QU*AO2B(JW)
      TABJ3  (JL,JC,JT,JO)=TABJ3(JL,JC,JT,JO)+QU*AO3(JW)*(1.0-QEO1D(JW))
      TABJ3A (JL,JC,JT,JO)=TABJ3A (JL,JC,JT,JO)+QU*AO3 (JW)*QEO1D(JW)
      TABJNO2(JL,JC,JT,JO)=TABJNO2(JL,JC,JT,JO)+QU*ANO2(JW)*QENO2(JW)
      TABJNO   (JL,JC,JT,JO)=TABJNO   (JL,JC,JT,JO) + QU*ANO   (JW)
      TABJNO31 (JL,JC,JT,JO)=TABJNO31 (JL,JC,JT,JO) + QU*ANO31 (JW)
      TABJNO32 (JL,JC,JT,JO)=TABJNO32 (JL,JC,JT,JO) + QU*ANO32 (JW)
      TABJN2O  (JL,JC,JT,JO)=TABJN2O  (JL,JC,JT,JO) + QU*AN2O  (JW)
      TABJN2O5 (JL,JC,JT,JO)=TABJN2O5 (JL,JC,JT,JO) + QU*AN2O5 (JW)
      TABJHNO3 (JL,JC,JT,JO)=TABJHNO3 (JL,JC,JT,JO) + QU*AHNO3 (JW)
      TABJPNA  (JL,JC,JT,JO)=TABJPNA  (JL,JC,JT,JO) + QU*APNA  (JW)
      TABJCNITA(JL,JC,JT,JO)=TABJCNITA(JL,JC,JT,JO) + QU*ACNITA(JW)
      TABJCNITB(JL,JC,JT,JO)=TABJCNITB(JL,JC,JT,JO) + QU*ACNITB(JW)
      TABJH2O2 (JL,JC,JT,JO)=TABJH2O2 (JL,JC,JT,JO) + QU*AH2O2 (JW)
      TABJH2O  (JL,JC,JT,JO)=TABJH2O  (JL,JC,JT,JO) + QU*AH2O  (JW)
      TABJHOCL (JL,JC,JT,JO)=TABJHOCL (JL,JC,JT,JO) + QU*AHOCL (JW)
      TABJHCL  (JL,JC,JT,JO)=TABJHCL  (JL,JC,JT,JO) + QU*AHCL  (JW)
      TABJCL2O2(JL,JC,JT,JO)=TABJCL2O2(JL,JC,JT,JO) + QU*ACL2O2(JW)
      TABJCH3CL(JL,JC,JT,JO)=TABJCH3CL(JL,JC,JT,JO) + QU*ACH3CL(JW)
      TABJCCL4 (JL,JC,JT,JO)=TABJCCL4 (JL,JC,JT,JO) + QU*ACCL4 (JW)
      TABJF11  (JL,JC,JT,JO)=TABJF11  (JL,JC,JT,JO) + QU*AF11  (JW)
      TABJF12  (JL,JC,JT,JO)=TABJF12  (JL,JC,JT,JO) + QU*AF12  (JW)
      TABJF22  (JL,JC,JT,JO)=TABJF22  (JL,JC,JT,JO) + QU*AF22  (JW)
      TABJF113 (JL,JC,JT,JO)=TABJF113 (JL,JC,JT,JO) + QU*AF113 (JW)
      TABJBRO  (JL,JC,JT,JO)=TABJBRO  (JL,JC,JT,JO) + QU*ABRO  (JW)
      TABJHOBR (JL,JC,JT,JO)=TABJHOBR (JL,JC,JT,JO) + QU*AHOBR (JW)
      TABJBRNO3(JL,JC,JT,JO)=TABJBRNO3(JL,JC,JT,JO) + QU*ABRNO3(JW)
      TABJBRCL (JL,JC,JT,JO)=TABJBRCL (JL,JC,JT,JO) + QU*ABRCL (JW)
      TABJOCLO (JL,JC,JT,JO)=TABJOCLO (JL,JC,JT,JO) + QU*AOCLO (JW)
      TABJC2OA (JL,JC,JT,JO)=TABJC2OA (JL,JC,JT,JO) + QU*AC2OA (JW)
      TABJC2OB (JL,JC,JT,JO)=TABJC2OB (JL,JC,JT,JO) + QU*AC2OB (JW)
      TABJMHP  (JL,JC,JT,JO)=TABJMHP  (JL,JC,JT,JO) + QU*AMHP  (JW)
      TABJCH3BR(JL,JC,JT,JO)=TABJCH3BR(JL,JC,JT,JO) + QU*ACH3BR(JW)
      TABJF12B1(JL,JC,JT,JO)=TABJF12B1(JL,JC,JT,JO) + QU*AF12B1(JW)
      TABJF13B1(JL,JC,JT,JO)=TABJF13B1(JL,JC,JT,JO) + QU*AF13B1(JW)
      TABJCOF2 (JL,JC,JT,JO)=TABJCOF2 (JL,JC,JT,JO) + QU*ACOF2 (JW)
      TABJCOFCL(JL,JC,JT,JO)=TABJCOFCL(JL,JC,JT,JO) + QU*ACOFCL(JW)
      TABJMCFM (JL,JC,JT,JO)=TABJMCFM (JL,JC,JT,JO) + QU*AMCFM (JW)
      TABJCH4  (JL,JC,JT,JO)=TABJCH4  (JL,JC,JT,JO) + QU*ACH4  (JW)
      TABJCO2  (JL,JC,JT,JO)=TABJCO2  (JL,JC,JT,JO) + QU*ACO2  (JW)
      TABJCOS  (JL,JC,JT,JO)=TABJCOS  (JL,JC,JT,JO) + QU*AOCS  (JW)
!
           END DO
!
!
!        End of the zenith angle loop.
         END DO
!
!     End of level loop.
      END DO
!
!     End of temperature loop.
      END DO
!
!     End of O3 profile loop.
      END DO
!
!
!      IF (mype == 0) THEN
!     Open file for the j tables.
!        OPEN(NTAB,FILE='jtable',FORM='UNFORMATTED',STATUS='UNKNOWN')
!
!     Write out the j rates.
!        WRITE (NTAB) TABPRES ,TABANG   ,TABT     ,TABO3    ,SAO3C    ,
!     1               TABJ2A  ,TABJ2B   ,TABJ3    ,TABJCH4  ,
!     2               TABJ3A  ,TABJCCL4 ,TABJCL2O2,TABJCNITA,TABCNITB ,TABJF11  ,
!     3               TABJF113,TABJF12  ,TABJF22  ,TABJH2O  ,TABJH2O2 ,
!     4               TABJHCL ,TABJHNO3 ,TABJHOCL ,TABJCH3CL,TABJN2O  ,
!     5               TABJN2O5,TABJNO   ,TABJNO2  ,TABJNO31 ,TABJNO32 ,
!     6               TABJPNA ,TABJBRCL ,TABJBRNO3,TABJBRO  ,TABJMCFM ,
!     7               TABJHOBR,TABJOCLO ,TABJC2OA ,TABJC2OB ,TABJMHP  ,
!     8               TABJCOF2,TABJCH3BR,TABJF12B1,TABJF13B1,TABJCOFCL,
!     9               TABJCO2 ,TABJCOS
!
!        CLOSE(NTAB)
!      END IF
!
      END SUBROUTINE INIJTAB

!
      SUBROUTINE LYMANA(TSPO2,TABJ2,TABJH2O,TABJCH4,QUANTA,JO)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Photodissociation due to Lyman alpha for O2, H2O, CH4.
!
!     Parameterisation based on Nicolet (see Brasseur + Solomon).
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
#include "parpho.h"
!
! Subroutine interface
      INTEGER, INTENT(IN) :: JO
      REAL, INTENT(IN) :: TSPO2  (JPLEV,JPCHI)
      REAL, INTENT(IN) :: QUANTA(JPWAV)
      REAL, INTENT(INOUT) :: TABJ2  (JPLEV,JPCHI,JPTEM,JPO3P)
      REAL, INTENT(INOUT) :: TABJH2O(JPLEV,JPCHI,JPTEM,JPO3P)
      REAL, INTENT(INOUT) :: TABJCH4(JPLEV,JPCHI,JPTEM,JPO3P)

! Local variables
      REAL :: EXPFAC
      REAL :: QULY
      REAL :: SIGLYA
      REAL :: PLY
!
      INTEGER :: JC
      INTEGER :: JL
      INTEGER :: JT
!
!
!     Only for zenith angles less than 90 degrees
      DO JT=1,JPTEM
      DO JC=1,JPCHIN-1
      DO JL=1,JPLEV
!
!       slant column of O2
        PLY=TSPO2(JL,JC)

        EXPFAC= 4.17E-19*(PLY**0.917)
        SIGLYA=EXPFAC/PLY
        QULY=QUANTA(1)*EXP(-EXPFAC)
!
        TABJ2  (JL,JC,JT,JO) = SIGLYA*QULY
        TABJH2O(JL,JC,JT,JO) = 1.40E-17*0.85*QULY
        TABJCH4(JL,JC,JT,JO) = 1.37E-17*0.85*QULY
!
      END DO
      END DO
      END DO
!
      END SUBROUTINE LYMANA

      SUBROUTINE ACSCCL4(T,JPWAV,WAVENM,ACCL4)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name     Type               Description.
!     T        Real               Temperature in kelvin.
!
!     WAVENM   Array of real      Wavelength of each interval in nm.
!
!     ACCL4    Array of real      Absorption cross-section of CCl4 in
!                                 cm^2 for temperature T for each
!                                 interval wavelength.
!
!     A subroutine which calculates the carbon tetra chloride (CCl4)
!     absorption cross section based on P. C. Simon et al. (1988).
!
!     Journal of atmospheric chemistry, Vol. 7, pp. 107-135, 1988.
!
!     This is done via a polynomial expression of the form;
!     log10(sigma)=A(lamda) + T B(lamda)
!
!     The expression is valid for the wavelength range; 194-250 nm
!                            and the temperature range; 210-300 K.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: JPWAV
      REAL, INTENT(IN)    :: T
      REAL, INTENT(IN)    :: WAVENM(JPWAV)
      REAL, INTENT(INOUT) :: ACCL4(JPWAV)
!
! Local variables.
!     Polynomial coefficients.
      REAL, PARAMETER :: A0=-37.104
      REAL, PARAMETER :: A1=-5.8218E-1
      REAL, PARAMETER :: A2= 9.9974E-3
      REAL, PARAMETER :: A3=-4.6765E-5
      REAL, PARAMETER :: A4= 6.8501E-8
!
      REAL, PARAMETER :: B0= 1.0739
      REAL, PARAMETER :: B1=-1.6275E-2
      REAL, PARAMETER :: B2= 8.8141E-5
      REAL, PARAMETER :: B3=-1.9811E-7
      REAL, PARAMETER :: B4= 1.5022E-10

      INTEGER :: JW
      REAL :: LAM
      REAL :: LAM2
      REAL :: LAM3
      REAL :: LAM4
      REAL :: ARG
      REAL :: TC
!
!     Check that temperature is in range.
      TC = T
      IF ( TC > 300.0 ) TC = 300.0
      IF ( TC < 210.0 ) TC = 210.0
!
!     Wavelength loop.
      DO JW = 57 , 79
!
!        Wavelength in nm.
         LAM = WAVENM(JW)
!
         LAM2 = LAM*LAM
         LAM3 = LAM*LAM2
         LAM4 = LAM*LAM3
!
         ARG = A0 + A1*LAM + A2*LAM2 + A3*LAM3 + A4*LAM4 +              &
               TC*(B0+B1*LAM+B2*LAM2+B3*LAM3+B4*LAM4)
!
         ACCL4(JW) = 10.0**ARG
!
!     End of the wavelength loop.
      END DO
!
      END SUBROUTINE ACSCCL4

      SUBROUTINE ACSF11(T,WAVENM,AF11)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Calculate T-dependent F11 cross sections
!
!     Name     Type               Description.
!     T        Real               Temperature in kelvin.
!
!     WAVENM   Array of real      Wavelength of each interval in nm.
!
!     AF11     Array of real      Absorption cross-section of F11 in
!                                 cm^2 for temperature T for each
!                                 interval wavelength.
!
!     AF11T    Array of real      Absorption cross-section of F11 at 298K
!
!     The expression is valid for the wavelength range; 174-230 nm
!                            and the temperature range; 210-300 K.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!  v1.2 17/1/2006 Small bug corrected  Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
#include "parpho.h"
! Subroutine interface
      REAL, INTENT(IN) :: T
      REAL, INTENT(IN) :: WAVENM(JPWAV)
      REAL, INTENT(INOUT) :: AF11(JPWAV)

! Local variables
      INTEGER :: JW
      REAL :: TC
      REAL :: ARG

      REAL :: AF11T(JPWAV)
!
!     CFCl3: JPL 1992 T=298K values.
      DATA AF11T/                                                       &
        36*0.0,                                                         &
        0.000E+00,0.000E+00,0.000E+00,0.000E+00,3.139E-19,2.478E-18,    &
        3.182E-18,3.168E-18,3.141E-18,3.098E-18,3.042E-18,3.096E-18,    &
        2.968E-18,2.765E-18,2.558E-18,2.319E-18,2.107E-18,1.839E-18,    &
        1.574E-18,1.332E-18,1.092E-18,8.911E-19,7.221E-19,5.751E-19,    &
        4.389E-19,3.340E-19,2.377E-19,1.700E-19,1.171E-19,7.662E-20,    &
        5.082E-20,3.184E-20,1.970E-20,1.206E-20,8.000E-21,4.834E-21,    &
        2.831E-21,1.629E-21,9.327E-22,5.209E-22,3.013E-22,1.617E-22,    &
        9.035E-23,5.427E-23,3.474E-23,2.141E-23,9.102E-24,1.499E-25,    &
        119*0.0/
!
!     Check that temperature is in range.
      TC = MIN (T, 300.0)
      TC = MAX (TC,210.0)
!
!     Wavelength loop.
      DO JW = 45 , 72
         ARG = 1.0E-4*(WAVENM(JW)-184.9)*(TC-298.0)
         AF11(JW) = AF11T(JW)*EXP(ARG)
      END DO
!
      END SUBROUTINE ACSF11

      SUBROUTINE ACSF12(T,WAVENM,AF12)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Calculate T-dependent F12 cross sections
!
!     Name     Type               Description.
!     T        Real               Temperature in kelvin.
!
!     WAVENM   Array of real      Wavelength of each interval in nm.
!
!     AF12     Array of real      Absorption cross-section of F12 in
!                                 cm2 for temperature T for each
!                                 interval wavelength.
!
!     AF12T    Array of real      Absorption cross-section of F12 at 298K
!
!     The expression is valid for the wavelength range; 174-226 nm
!                            and the temperature range; 210-300 K.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
#include "parpho.h"
! Subroutine interface
      REAL, INTENT(IN) :: T
      REAL, INTENT(IN) :: WAVENM(JPWAV)
      REAL, INTENT(INOUT) :: AF12(JPWAV)

! Local variables
      INTEGER :: JW
      REAL :: TC
      REAL :: ARG
      REAL :: AF12T(JPWAV)
!
!     CF2Cl2: JPL 1992 T=298K values.
      DATA AF12T/                                                       &
        36*0.0,                                                         &
        0.000E+00,0.000E+00,0.000E+00,0.000E+00,9.716E-20,8.639E-19,    &
        1.370E-18,1.630E-18,1.752E-18,1.846E-18,1.894E-18,1.778E-18,    &
        1.655E-18,1.515E-18,1.312E-18,1.030E-18,8.615E-19,6.682E-19,    &
        4.953E-19,3.567E-19,2.494E-19,1.659E-19,1.088E-19,7.081E-20,    &
        4.327E-20,2.667E-20,1.753E-20,9.740E-21,5.336E-21,2.976E-21,    &
        2.572E-21,3.840E-21,5.644E-22,3.270E-22,1.769E-22,8.850E-23,    &
        4.328E-23,2.236E-23,1.040E-23,3.751E-24,1.146E-24,0.000E+00,    &
        125*0.0/
!
!
!     Check that temperature is in range.
      TC = MIN (T, 300.0)
      TC = MAX (T, 210.0)
!
!     Wavelength loop.
      DO JW = 45 , 71
         ARG = 4.1E-4*(WAVENM(JW)-184.9)*(TC-298.0)
         AF12(JW) = AF12T(JW)*EXP(ARG)
      END DO
!
      END SUBROUTINE ACSF12

      SUBROUTINE ACSF22(T,JPWAV,WAVENM,AF22)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name     Type               Description.
!     T        Real               Temperature in kelvin.
!
!     WAVENM   Array of real      Wavelength of each interval in nm.
!
!     AF22     Array of real      Absorption cross-section of F22 in
!                                 cm^2 for temperature T for each
!                                 interval wavelength.
!
!     A subroutine which calculates the F22 (CHF2Cl) absorption
!     cross section based on P. C. Simon et al. (1988).
!
!     Journal of atmospheric chemistry, Vol. 7, pp. 107-135, 1988.
!
!     This is done via a polynomial expression of the form;
!     log10(sigma)=A(lamda) + T B(lamda)
!
!     The expression is valid for the wavelength range; 174-204 nm
!                            and the temperature range; 210-300 K.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: JPWAV
      REAL, INTENT(IN)    :: T
      REAL, INTENT(IN)    :: WAVENM(JPWAV)
      REAL, INTENT(INOUT) :: AF22(JPWAV)

! Local variables.
!
!     Polynomial coefficients.
      REAL, PARAMETER :: A0=-106.029
      REAL, PARAMETER :: A1= 1.5038
      REAL, PARAMETER :: A2=-8.2476E-3
      REAL, PARAMETER :: A3= 1.4206E-5
!
      REAL, PARAMETER :: B0=-1.3399E-1
      REAL, PARAMETER :: B1= 2.7405E-3
      REAL, PARAMETER :: B2=-1.8028E-5
      REAL, PARAMETER :: B3= 3.8504E-8

      INTEGER :: JW
      REAL :: TC
!
      REAL :: LAM
      REAL :: LAM2
      REAL :: LAM3
      REAL :: ARG
!
!
!     Check that temperature is in range.
      TC = T
      IF ( TC > 300.0 ) TC = 300.0
      IF ( TC < 210.0 ) TC = 210.0
!
!     Wavelength.
      DO JW = 45 , 61
!
!        Wavelength in nm.
         LAM = WAVENM(JW)
!
         LAM2 = LAM*LAM
         LAM3 = LAM*LAM2
         ARG = A0 + A1*LAM + A2*LAM2 + A3*LAM3 +                        &
               TC*(B0+B1*LAM+B2*LAM2+B3*LAM3)
         AF22(JW) = 10.0**ARG
!
      END DO
!
      END SUBROUTINE ACSF22

      SUBROUTINE ACSHNO3(TEMP,WAVENM,AHNO3)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Calculate HNO3 temperature dependent cross sections,
!
!     AHNO3    Array of real   Contains the cross sections at model
!                              wavelengths.Contains the result on exit.
!
!     AHNO3T   Array of real   Contains the cross sections at 298K.
!
!     B        Array of real   Contains B coefficient.
!
!     WAVENM   Array of real   Contains the wavelength intervals.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
#include "parpho.h"
!
! Subroutine interface
      REAL, INTENT(IN)    :: TEMP
      REAL, INTENT(IN)    :: WAVENM(JPWAV)
      REAL, INTENT(INOUT) :: AHNO3(JPWAV)

! Local variables
      REAL :: TC
      REAL :: AHNO3T(JPWAV)
      REAL :: B     (JPWAV)
!
!     Intercept B from Burkholder
      DATA B/54*0.0,                                                    &
        0.000E+00,0.000E+00,1.152E-03,1.668E-03,1.653E-03,1.673E-03,    &
        1.720E-03,1.750E-03,1.817E-03,1.935E-03,2.060E-03,2.168E-03,    &
        2.178E-03,2.195E-03,2.106E-03,1.987E-03,1.840E-03,1.782E-03,    &
        1.838E-03,1.897E-03,1.970E-03,1.978E-03,1.855E-03,1.655E-03,    &
        1.416E-03,1.247E-03,1.162E-03,1.121E-03,1.136E-03,1.199E-03,    &
        1.315E-03,1.493E-03,1.637E-03,1.767E-03,1.928E-03,2.139E-03,    &
        2.380E-03,2.736E-03,3.139E-03,3.695E-03,4.230E-03,5.151E-03,    &
        6.450E-03,7.327E-03,9.750E-03,1.013E-02,1.180E-02,1.108E-02,    &
        9.300E-03,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,    &
        95*0.0/
!
!     AHNO3 AT 298K
      DATA AHNO3T/48*0.0,                                               &
        0.000E+00,8.305E-18,1.336E-17,1.575E-17,1.491E-17,1.385E-17,    &
        1.265E-17,1.150E-17,1.012E-17,8.565E-18,6.739E-18,5.147E-18,    &
        3.788E-18,2.719E-18,1.796E-18,1.180E-18,7.377E-19,4.487E-19,    &
        2.810E-19,1.826E-19,1.324E-19,1.010E-19,8.020E-20,6.479E-20,    &
        5.204E-20,4.178E-20,3.200E-20,2.657E-20,2.298E-20,2.086E-20,    &
        1.991E-20,1.962E-20,1.952E-20,1.929E-20,1.882E-20,1.804E-20,    &
        1.681E-20,1.526E-20,1.335E-20,1.136E-20,9.242E-21,7.186E-21,    &
        5.320E-21,3.705E-21,2.393E-21,1.442E-21,8.140E-22,4.131E-22,    &
        1.970E-22,9.434E-23,4.310E-23,2.204E-23,1.030E-23,5.841E-24,    &
        4.170E-24,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,    &
        95*0.0/
!
!
      TC=TEMP
      TC=MAX(TC, 200.0)
      TC=MIN(TC, 360.0)
!
      AHNO3(50:103) = AHNO3T(50:103)*EXP(B(50:103)*(TC-298.0))
!
      END SUBROUTINE ACSHNO3

      SUBROUTINE ACSCNIT(TEMP,WAVENM,ACNITA, ACNITB)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Calculate ClONO2 temperature dependent cross sections,
!
!     ACNITA, ACNITB Arrays of real
!                      Contain the cross sections times quantum yields
!                              at model
!                              wavelengths.Contains the result on exit.
!
!     ACNITT   Array of real   Contains the cross sections at 296K.
!
!     A1       Array of real   Contains A1 coefficient.
!     A2       Array of real   Contains A2 coefficient.
!
!     WAVENM   Array of real   Contains the wavelength intervals.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!  v1.2 Introduce branching ratio for channels ClO + NO2 / Cl + NO3
!                                      Olaf Morgenstern 12/1/2006
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
#include "parpho.h"
!
! Subroutine interface
      REAL, INTENT(IN)    :: TEMP
      REAL, INTENT(IN)    :: WAVENM(JPWAV)
      REAL, INTENT(INOUT) :: ACNITA(JPWAV)
      REAL, INTENT(INOUT) :: ACNITB(JPWAV)

! Local variables
      REAL :: TC
      LOGICAL, SAVE :: first = .TRUE.
      REAL :: ACNITT(JPWAV)
      REAL :: A1    (JPWAV)
      REAL :: A2    (JPWAV)
      REAL, SAVE :: QUANTY(JPWAV) ! Quantum yield for Cl + NO3 channel
!
!     Coeffs A1, A2 from Burkholder et al GRL 1994
      DATA A1/54*0.0,                                                   &
         1.73E-05, 7.50E-05, 1.00E-04, 8.82E-05, 3.61E-05,-5.88E-05,    &
        -1.95E-04,-3.44E-04,-5.11E-04,-6.59E-04,-7.85E-04,-8.71E-04,    &
        -9.03E-04,-8.73E-04,-7.83E-04,-6.38E-04,-4.53E-04,-2.35E-04,    &
        -6.97E-06, 2.19E-04, 4.16E-04, 5.64E-04, 6.78E-04, 7.81E-04,    &
         9.08E-04, 1.08E-03, 1.26E-03, 1.44E-03, 1.59E-03, 1.74E-03,    &
         1.88E-03, 2.03E-03, 2.21E-03, 2.37E-03, 2.55E-03, 2.74E-03,    &
         2.95E-03, 3.28E-03, 3.72E-03, 4.12E-03, 4.53E-03, 4.98E-03,    &
         5.40E-03, 5.75E-03, 5.92E-03, 5.85E-03, 5.51E-03, 4.92E-03,    &
         4.02E-03, 3.27E-03, 2.70E-03, 2.08E-03, 1.33E-03, 7.65E-04,    &
         3.53E-04, 2.39E-04, 4.10E-04, 7.77E-04, 1.38E-03, 2.15E-03,    &
         3.38E-03, 4.88E-03, 6.70E-03, 8.10E-03, 9.72E-03, 9.96E-03,    &
         83*0.0/
      DATA A2/54*0.0,                                                   &
        -1.16E-06,-5.32E-06,-7.85E-06,-8.21E-06,-7.81E-06,-7.52E-06,    &
        -7.46E-06,-7.62E-06,-7.93E-06,-8.30E-06,-8.69E-06,-9.03E-06,    &
        -9.25E-06,-9.37E-06,-9.37E-06,-9.27E-06,-9.06E-06,-8.65E-06,    &
        -7.98E-06,-7.13E-06,-6.36E-06,-6.00E-06,-6.09E-06,-6.45E-06,    &
        -6.57E-06,-6.03E-06,-5.08E-06,-4.20E-06,-3.50E-06,-2.74E-06,    &
        -1.87E-06,-9.85E-07, 6.15E-08, 1.13E-06, 2.14E-06, 3.05E-06,    &
         3.74E-06, 5.29E-06, 7.78E-06, 9.88E-06, 1.20E-05, 1.49E-05,    &
         1.84E-05, 2.27E-05, 2.70E-05, 3.01E-05, 3.11E-05, 2.86E-05,    &
         2.07E-05, 1.38E-05, 8.59E-06, 2.01E-06,-7.40E-06,-1.44E-05,    &
        -1.91E-05,-2.11E-05,-2.05E-05,-1.87E-05,-1.42E-05,-7.14E-06,    &
         4.47E-06, 1.93E-05, 3.87E-05, 5.57E-05, 7.52E-05, 7.81E-05,    &
         83*0.0/
!
!     ACNIT at 296K
      DATA ACNITT/54*0.0,                                               &
        4.320E-19,1.980E-18,2.911E-18,3.015E-18,2.871E-18,2.785E-18,    &
        2.780E-18,2.839E-18,2.956E-18,3.097E-18,3.264E-18,3.386E-18,    &
        3.448E-18,3.392E-18,3.236E-18,2.971E-18,2.640E-18,2.268E-18,    &
        1.922E-18,1.591E-18,1.314E-18,1.086E-18,8.967E-19,7.444E-19,    &
        6.074E-19,5.129E-19,4.352E-19,3.703E-19,3.152E-19,2.662E-19,    &
        2.213E-19,1.840E-19,1.498E-19,1.211E-19,9.519E-20,7.333E-20,    &
        5.500E-20,4.007E-20,2.969E-20,2.190E-20,1.600E-20,1.142E-20,    &
        8.310E-21,6.114E-21,4.660E-21,3.657E-21,3.020E-21,2.576E-21,    &
        2.290E-21,2.079E-21,2.000E-21,1.795E-21,1.590E-21,1.414E-21,    &
        1.210E-21,1.056E-21,9.090E-22,7.588E-22,6.380E-22,5.376E-22,    &
        4.440E-22,3.672E-22,3.160E-22,2.314E-22,1.890E-22,5.264E-23,    &
        83*0.0/
!
!
      IF (first) THEN
        QUANTY = MIN(MAX(7.143e-3 * WAVENM - 1.6, 0.6), 1.0)
        first = .FALSE.
      END IF

      TC=TEMP
      TC=MAX(TC, 220.0)
      TC=MIN(TC, 298.0)
!
! ClONO2 cross section, following JPL (2002)
      ACNITA(55:120) = ACNITT(55:120)*(1.0 + A1(55:120)*(TC-296.0)      &
                     + A2(55:120)*((TC-296.0)**2.0))
! ClONO2 cross section, times quantum yield for ClO + NO2 channel
      ACNITB(55:120) = (1.0 - QUANTY(55:120)) * ACNITA(55:120)
! ClONO2 cross section, times quantum yield for Cl  + NO3 channel
      ACNITA(55:120) = QUANTY(55:120) * ACNITA(55:120)
!
      END SUBROUTINE ACSCNIT

      SUBROUTINE ACSH2O2(TEMP,JPWAV,WAVENM,AH2O2)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Calculate H2O2 temperature dependent cross sections,
!
!     AH2O2    Array of real   Contains the cross sections at model
!                              wavelengths.Contains the result on exit.
!
!     WAVENM   Array of real   Contains the wavelength intervals.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: jpwav
      REAL, INTENT(IN)    :: temp
      REAL, INTENT(IN)    :: wavenm(jpwav)
      REAL, INTENT(INOUT) :: ah2o2(jpwav)

! Local variables
      REAL, PARAMETER :: A0= 6.4761E4
      REAL, PARAMETER :: A1=-9.2170972E2
      REAL, PARAMETER :: A2= 4.535649
      REAL, PARAMETER :: A3=-4.4589016E-3
      REAL, PARAMETER :: A4=-4.035101E-5
      REAL, PARAMETER :: A5= 1.6878206E-7
      REAL, PARAMETER :: A6=-2.652014E-10
      REAL, PARAMETER :: A7= 1.5534675E-13
      REAL, PARAMETER :: B0= 6.8123E3
      REAL, PARAMETER :: B1=-5.1351E1
      REAL, PARAMETER :: B2= 1.1522E-1
      REAL, PARAMETER :: B3=-3.0493E-5
      REAL, PARAMETER :: B4=-1.0924E-7
!
      REAL :: CHI
      REAL :: LAM
      REAL :: LAM2
      REAL :: LAM3
      REAL :: LAM4
      REAL :: LAM5
      REAL :: LAM6
      REAL :: LAM7
      REAL :: TC

      INTEGER JW
!
!
!     Check that temperature is in range.
      TC = TEMP
      IF ( TC > 400.0 ) TC = 400.0
      IF ( TC < 200.0 ) TC = 200.0
!
      CHI=1.0/(1.0 + EXP(-1265.0/TC))
!
!     Wavelength loop 260nm - 350nm
      DO JW = 83,103
         LAM =WAVENM(JW)
         LAM2=LAM*LAM
         LAM3=LAM*LAM2
         LAM4=LAM*LAM3
         LAM5=LAM*LAM4
         LAM6=LAM*LAM5
         LAM7=LAM*LAM6
         AH2O2(JW)= 1.0E-21*(CHI*(A0 + A1*LAM  + A2*LAM2 + A3*LAM3      &
           +                 A4*LAM4 + A5*LAM5 + A6*LAM6 + A7*LAM7)     &
           + (1.0-CHI)*(B0 + B1*LAM  + B2*LAM2 + B3*LAM3 + B4*LAM4))
      END DO
!
      END SUBROUTINE ACSH2O2

      SUBROUTINE ACSBRCL(TEMP,JPWAV,WAVENM,ABRCL)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Calculate BrCl temperature dependent cross sections,
!
!     ABRCL    Array of real   Contains the cross sections at model
!                              wavelengths.Contains the result on exit.
!
!     WAVENM   Array of real   Contains the wavelength intervals.
!
!     Parameterisation from Maric et al [1994]
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: jpwav
      REAL, INTENT(IN)    :: temp
      REAL, INTENT(IN)    :: wavenm(jpwav)
      REAL, INTENT(INOUT) :: abrcl(jpwav)

! Local variables
      REAL, PARAMETER :: A1=7.34E-20
      REAL, PARAMETER :: A2=4.35E-19
      REAL, PARAMETER :: A3=1.12E-19
      REAL, PARAMETER :: B1= 68.6
      REAL, PARAMETER :: B2=123.6
      REAL, PARAMETER :: B3= 84.8
      REAL, PARAMETER :: C1=227.6
      REAL, PARAMETER :: C2=372.5
      REAL, PARAMETER :: C3=442.4
      REAL, PARAMETER :: WE=443.1
      REAL, PARAMETER :: H=6.63E-34
      REAL, PARAMETER :: C=3.0E10
      REAL, PARAMETER :: BOLTZ=1.38E-23
!
      REAL :: LAM
      REAL :: TANT
      REAL :: TC

      INTEGER :: JW
!
!
!     Check that temperature is in range.
      TC = TEMP
      TC=MIN(TC, 300.0)
      TC=MAX(TC, 200.0)
!
      TANT=TANH(H*C*WE/(2.0*BOLTZ*TC))
!
!     Wavelength loop 200nm - 516nm (estimated upper limit for Br-Cl bond)
      DO JW = 60,136
        LAM =WAVENM(JW)
        ABRCL(JW)=                                                      &
          A1*SQRT(TANT)*EXP(-B1*TANT*(LOG(C1/LAM))**2.0)                &
        + A2*SQRT(TANT)*EXP(-B2*TANT*(LOG(C2/LAM))**2.0)                &
        + A3*SQRT(TANT)*EXP(-B3*TANT*(LOG(C3/LAM))**2.0)
      END DO
!
      END SUBROUTINE ACSBRCL

      SUBROUTINE ACSNO2(TEMP,WAVENM,ANO2)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Calculate NO2 temperature dependent cross sections,
!
!     ANO2      Array of real   Contains the cross sections at model
!                              wavelengths.Contains the result on exit.
!
!     ANO2T    Array of real   Contains the cross sections at 298K.
!
!     A        Array of real   Contains A coefficient.
!
!     WAVENM   Array of real   Contains the wavelength intervals.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
#include "parpho.h"
!
! Subroutine interface
      REAL, INTENT(IN)    :: temp
      REAL, INTENT(IN)    :: wavenm(jpwav)
      REAL, INTENT(INOUT) :: ano2(jpwav)

! Local variables
      REAL ::  ANO2T (JPWAV)
      REAL ::  A     (JPWAV)
!
!     slope A
      DATA A/                                                           &
        86*0.0,0.075, 0.082,-0.053,-0.043,-0.031,-0.162,-0.284,         &
       -0.357,-0.536,-0.686,-0.786,-1.105,-1.355,-1.277,-1.612,         &
       -1.890,-1.219,-1.921,-1.095,-1.322,-1.102,-0.806,-0.867,         &
       -0.945,-0.923,-0.738,-0.599,-0.545,-1.129, 0.001,-1.208,         &
        86*0.0/
!
!     ANO2 at 298K
      DATA ANO2T/                                                       &
        48*0.0,                                                         &
        0.000E+00,0.000E+00,0.000E+00,2.670E-19,2.780E-19,2.900E-19,    &
        2.790E-19,2.600E-19,2.420E-19,2.450E-19,2.480E-19,2.750E-19,    &
        4.145E-19,4.478E-19,4.454E-19,4.641E-19,4.866E-19,4.818E-19,    &
        5.022E-19,4.441E-19,4.713E-19,3.772E-19,3.929E-19,2.740E-19,    &
        2.778E-19,1.689E-19,1.618E-19,8.812E-20,7.472E-20,3.909E-20,    &
        2.753E-20,2.007E-20,1.973E-20,2.111E-20,2.357E-20,2.698E-20,    &
        3.247E-20,3.785E-20,5.030E-20,5.880E-20,7.000E-20,8.150E-20,    &
        9.720E-20,1.154E-19,1.344E-19,1.589E-19,1.867E-19,2.153E-19,    &
        2.477E-19,2.807E-19,3.133E-19,3.425E-19,3.798E-19,4.065E-19,    &
        4.313E-19,4.717E-19,4.833E-19,5.166E-19,5.315E-19,5.508E-19,    &
        5.644E-19,5.757E-19,5.927E-19,5.845E-19,6.021E-19,5.781E-19,    &
        5.999E-19,5.651E-19,5.812E-19,0.000E+00,0.000E+00,0.000E+00,    &
        83*0.0/
!
!
      ANO2(87:117) = ANO2T(87:117) + 1.0E-22*A(87:117)*(TEMP-273.0)
!
      END SUBROUTINE ACSNO2

      SUBROUTINE ACSMC(T,JPWAV,WAVENM,AMC)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name     Type               Description.
!     T        Real               Temperature in kelvin.
!
!     WAVENM   Array of real      Wavelength of each interval in nm.
!
!     AMC      Array of real      Absorption cross-section of CH3Cl in
!                                 cm^2 for temperature T for each
!                                 interval wavelength.
!
!     Subroutine which calculates the methyl chloride (CH3Cl) absorption
!     cross section based on P. C. Simon et al. (1988).
!
!     Journal of atmospheric chemistry, Vol. 7, pp. 107-135, 1988.
!
!     This is done via a polynomial expression of the form;
!     log10(sigma)=A(lamda) + T B(lamda)
!
!     The expression is valid for the wavelength range; 174-226 nm
!                            and the temperature range; 210-300 K.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Subroutine interface
      REAL, INTENT(IN)    :: T
      INTEGER, INTENT(IN) :: JPWAV
      REAL, INTENT(IN)    :: WAVENM(JPWAV)
      REAL, INTENT(INOUT) :: AMC(JPWAV)

! Local variables
!
!     Polynomial coefficients.
      REAL, PARAMETER :: A0=-299.80
      REAL, PARAMETER :: A1= 5.1047
      REAL, PARAMETER :: A2=-3.363E-2
      REAL, PARAMETER :: A3= 9.5805E-5
      REAL, PARAMETER :: A4=-1.0135E-7
!
      REAL, PARAMETER :: B0=-7.1727
      REAL, PARAMETER :: B1= 1.4837E-1
      REAL, PARAMETER :: B2=-1.1463E-3
      REAL, PARAMETER :: B3= 3.9188E-6
      REAL, PARAMETER :: B4=-4.9994E-9
!
      INTEGER :: JW

!     Temperature in kelvin.
      REAL :: TC
      REAL :: LAM
      REAL :: LAM2
      REAL :: LAM3
      REAL :: LAM4
      REAL :: ARG
!
!     Check that temperature is in range.
      TC = T
      IF ( TC > 300.0 ) TC = 300.0
      IF ( TC < 210.0 ) TC = 210.0
!
!     Wavelength.
      DO JW = 45 , 67
!
!        Wavelength in nm.
         LAM = WAVENM(JW)
!
         LAM2 = LAM*LAM
         LAM3 = LAM*LAM2
         LAM4 = LAM*LAM3
         ARG = A0 + A1*LAM + A2*LAM2 + A3*LAM3 + A4*LAM4 +              &
               TC*(B0+B1*LAM+B2*LAM2+B3*LAM3+B4*LAM4)
         AMC(JW) = 10.0**ARG
!
      END DO
!
      END SUBROUTINE ACSMC

      SUBROUTINE ACSN2O(T,JPWAV,WAVENM,AN2O)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name     Type               Description.
!     T        Real               Temperature in kelvin.
!
!     WAVENM   Array of real      Wavelength of each interval in nm.
!
!     AN2O     Array of real      Absorption cross-section of N2O in
!                                 cm^2 for temperature T for each
!                                 interval wavelength.
!
!      Calculate the temperature dependent N2O absorption cross
!      section between 173 nm and 240 nm.  This parameterisation
!      applies for temperatures between 194 and 320 K.  It is taken
!      from the work of Selwyn et al. (1977) presented in JPL EVAL 8
!      (1987) pp. 108. .
!
!      When the temperature is out of the parameterisation range;
!      < 194 K, the cross section for 194 K is used,
!      > 320 K, the cross section for 320 K is used.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: JPWAV
      REAL, INTENT(IN)    :: T
      REAL, INTENT(IN)    :: WAVENM(JPWAV)
      REAL, INTENT(INOUT) :: AN2O(JPWAV)

! Local variables
!     Polynomial coefficients used.
      REAL, PARAMETER :: A1N2O=68.210230
      REAL, PARAMETER :: A2N2O=-4.071805
      REAL, PARAMETER :: A3N2O= 4.301146E-2
      REAL, PARAMETER :: A4N2O=-1.777846E-4
      REAL, PARAMETER :: A5N2O= 2.520672E-7
!
      REAL, PARAMETER :: B1N2O=123.401400
      REAL, PARAMETER :: B2N2O=-2.116255
      REAL, PARAMETER :: B3N2O= 1.111572E-2
      REAL, PARAMETER :: B4N2O=-1.881058E-5
!
      INTEGER :: JW
!
!     Temperature in kelvin.
      REAL :: TC
      REAL :: TM300
      REAL :: lam1
      REAL :: lam2
      REAL :: lam3
      REAL :: lam4
      REAL :: arg
      REAL :: arga
      REAL :: argb
!
!     Check is not out of parameterisation range.
      TC = T
      IF ( TC < 194.0 ) TC = 194.0
      IF ( TC > 320.0 ) TC = 320.0
!
!     Temperature - 300.
      TM300 = TC - 300.0
!
!     Wavelengths between 173 nm & 240 nm.
      DO JW = 44 , 76
!
!        Various powers of wavelength in nm.
         LAM1 = WAVENM(JW)
         LAM2 = LAM1*LAM1
         LAM3 = LAM2*LAM1
         LAM4 = LAM3*LAM1
!
!        Evaluate the two polynomial expressions.
         ARGA = A1N2O + A2N2O*LAM1 + A3N2O*LAM2 + A4N2O*LAM3 +          &
                A5N2O*LAM4
         ARGB = B1N2O + B2N2O*LAM1 + B3N2O*LAM2 + B4N2O*LAM3
!
         ARG = ARGA + TM300*EXP(ARGB)
!
         AN2O(JW) = EXP(ARG)
!
      END DO
!
      END SUBROUTINE ACSN2O

      SUBROUTINE ACSN2O5(T,JPWAV,WAVENM,AN2O5)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name     Type               Description.
!     T        Real               Temperature in kelvin.
!
!     WAVENM   Array of real      Wavelength of each interval in nm.
!
!     AN2O5    Array of real      Absorption cross-section of N2O5 in
!                                 cm^2 for temperature T for each
!                                 interval wavelength.
!
!      Calculate the temperature dependent N2O5 cross section,
!      using the parameterisation of Yao et al. (1982), given in JPL
!      EVAL 8 (1987) pp. 110-111. .
!
!      This applies for wavelengths between 285 and 380 nm, and
!      temperatures 225 K to 300 K.
!
!      If the temperature is below 225 K, the value for 225 K is used.
!      If the temperature is above 300 K, the value for 300 K is used.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: JPWAV
      REAL, INTENT(IN)    :: T
      REAL, INTENT(IN)    :: WAVENM(JPWAV)
      REAL, INTENT(INOUT) :: AN2O5(JPWAV)

! Local variables
      INTEGER :: JW
      INTEGER :: LAM
!
!     Temperature in kelvin.
      REAL :: TC
      REAL :: ARG
!
!
!     Check is not out of parameterisation range.
      TC = T
      IF ( TC < 225.0 ) TC = 225.0
      IF ( TC > 300.0 ) TC = 300.0
!
!     Wavelengths between 285 nm & 380 nm.
      DO JW = 89 , 109
!
!        Wavelength in nm.
         LAM = WAVENM(JW)
!
!        Evaluate the parameterisation expression.
         ARG = 2.735 + ((4728.5-17.127*LAM)/TC)
!
!        Calculate the cross section.
         AN2O5(JW) = 1.0E-20*EXP(ARG)
!
      END DO
!
      END SUBROUTINE ACSN2O5

      SUBROUTINE ACSNO(ANGLE,P,VC,ANO)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name     Type               Description.
!     ANGLE    Real               Solar zenith angle in radians.
!
!     P        Real               Pressure in mb.
!
!     VC       Real               Vertical O2 column above the point
!                                 in molecules per cm^2.
!
!     ANO      Array of real      Absorption cross-section of NO in
!                                 cm^2 for pressure P and solar zenith
!                                 angle ANGLE.
!
!     Calculate the NO absorption cross section accounting for
!     the NO absorption in the DEL(0-0) & the DEL(0-1) bands.
!     Calculated from the Mark Allen & John E Frederick
!     parameterisation, taken from;
!
!     The Journal of atmospheric sciences, Vol. 39, September 82.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
#include "parpho.h"
!
! Subroutine interface
      REAL, INTENT(IN)    :: ANGLE
      REAL, INTENT(IN)    :: P
      REAL, INTENT(IN)    :: VC
      REAL, INTENT(INOUT) :: ANO(JPWAV)

! Local variables
      REAL, PARAMETER :: TORAD=pi_over_180
      INTEGER, PARAMETER :: JPNOB=9
      INTEGER, PARAMETER :: JPNOZ=5
!
      REAL :: CZ0
      REAL :: CZ1
      REAL :: SEC
      REAL :: SIGMAE0
      REAL :: SIGMAE1
      REAL :: ZLOGN
      REAL :: ZLOGN2
      REAL :: ZLOGN3
      REAL :: ZLOGN4
      REAL :: ZLOGP
      REAL :: ZLOGP2
      REAL :: ZLOGP3
      REAL :: ZLOGP4
      REAL :: ZLOGP5
      REAL :: ZLOGP6
      REAL :: ZLOGP7
      REAL :: ZLOGP8
      REAL :: C00
      REAL :: C01
!
!
!     NO parameterisation Allen & Frederick 1983.
      REAL :: BAND00(JPNOB)
      REAL :: ZEN00(JPNOZ)
      REAL :: BAND01(JPNOB)
      REAL :: ZEN01(JPNOZ)
!
!     Band del(0-0)
      DATA BAND00/ - 1.790868E1 , -1.924701E-1 , -7.217717E-2 ,         &
           5.648282E-2 , 4.569175E-2 , 8.353572E-3 , 0.0 , 0.0 , 0.0/
!
!     Band del(0-1)
      DATA BAND01/ - 1.654245E1 , 5.836899E-1 , 3.449436E-1 ,           &
           1.700653E-1 , -3.324717E-2 , -4.952424E-2 , 1.579306E-2 ,    &
           1.835462E-2 , 3.368125E-3/
!
!     Polynomial coefficients for the NO effective absorption cross
!     section Zenith angle dependence.
!
!     Band del(0-0)
      DATA ZEN00/7836.832 , -1549.88 , 114.8342 , -3.777754 ,           &
           4.655696E-2/
!
!     Band del(0-1)
      DATA ZEN01/12975.81 , -2582.981 , 192.7709 , -6.393008 ,          &
           7.949835E-2/
!
!
!     Initialise cross sections
      ANO(49:50)=0.0
      ANO(54:55)=0.0
!
!     Set up the two logical conditions.
!     i.e. Only continue if the parameterisation applies.
!
!          > 20 KM           < 84 degrees.
      IF ( (P <= 50.0) .AND. (ANGLE <= (84.0*TORAD)) ) THEN
!
!       Set values of the NO crossection for the TWO wavelength
!       intervals where the parameterisation applies over the
!       altitudes which it applies, i.e. ABOVE 20 kilometres ONLY.
!       This parameterisation was performed by Allen & Frederick
!       using high spectral resolution calculations of the delta
!       band predissociation of nitric oxide.
!
!       N.B. It is very important that the clause (ANGLE <= (84.0 deg))
!           is included. This limits the calculation to angles less than
!           84 degrees. At angles greater than this the
!           parameterisation is NOT valid.
!
!        Set Log(P), and Log(N), P being Pressure, N being the VERTICAL
!        O2 Column above the given point.
         ZLOGP = LOG10(P)
         ZLOGN = LOG10(VC)
!
!        Initialise Zenith angle dependence variables for both bands.
         C00 = 0.
         C01 = 0.
!
         SEC = (1/COS(ANGLE))
!
         ZLOGP2 = ZLOGP*ZLOGP
         ZLOGP3 = ZLOGP2*ZLOGP
         ZLOGP4 = ZLOGP3*ZLOGP
         ZLOGP5 = ZLOGP4*ZLOGP
         ZLOGP6 = ZLOGP5*ZLOGP
         ZLOGP7 = ZLOGP6*ZLOGP
         ZLOGP8 = ZLOGP7*ZLOGP
!
         ZLOGN2 = ZLOGN*ZLOGN
         ZLOGN3 = ZLOGN2*ZLOGN
         ZLOGN4 = ZLOGN3*ZLOGN
!
!
!        For BAND DEL(1-0) calculate the NO effective cross section for
!        an OVERHEAD sun using a simple polynomial expression.
         SIGMAE1 = BAND01(1) + BAND01(2)*ZLOGP + BAND01(3)              &
                   *ZLOGP2 + BAND01(4)*ZLOGP3 + BAND01(5)               &
                   *ZLOGP4 + BAND01(6)*ZLOGP5 + BAND01(7)               &
                   *ZLOGP6 + BAND01(8)*ZLOGP7 + BAND01(9)*ZLOGP8
!
         SIGMAE1 = 10.**SIGMAE1
!
!        For BAND DEL(1-0) calculate the Zenith angle dependence of the
!        NO effective cross section using a simple polynomial expression.
         CZ1 = ZEN01(1) + ZEN01(2)*ZLOGN + ZEN01(3)*ZLOGN2 + ZEN01(4)   &
               *ZLOGN3 + ZEN01(5)*ZLOGN4
!
         C01 = (SEC)**CZ1
!
!        Effective NO absorption cross sections for BAND DEL(1-0)
         ANO(49) = SIGMAE1*C01
         ANO(50) = ANO(49)
!
!
!        For BAND DEL(0-0) calculate the NO effective cross section for
!        an OVERHEAD sun using a simple polynomial expression.
         SIGMAE0 = BAND00(1) + BAND00(2)*ZLOGP + BAND00(3)              &
                   *ZLOGP2 + BAND00(4)*ZLOGP3 + BAND00(5)               &
                   *ZLOGP4 + BAND00(6)*ZLOGP5
!
         SIGMAE0 = 10.**SIGMAE0
!
!        For BAND DEL(0-0) calculate the Zenith angle dependence of the
!        NO effective cross section using a simple polynomial expression.
         CZ0 = ZEN00(1) + ZEN00(2)*ZLOGN + ZEN00(3)*ZLOGN2 + ZEN00(4)   &
               *ZLOGN3 + ZEN00(5)*ZLOGN4
!
         C00 = (SEC)**CZ0
!
!        Effective NO absorption cross sections for BAND DEL(0-0)
         ANO(54) = SIGMAE0*C00
         ANO(55) = ANO(54)
!
!     End of in range if statement.
      END IF
!
      END SUBROUTINE ACSNO

      SUBROUTINE ACSO3(T,JPWAV,AO3)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name     Type               Description.
!     T        Real               Temperature in kelvin.
!
!     WAVENM   Array of real      Wavelength of each interval in nm.
!
!     AO3      Array of real      Absorption cross-section of O3 in
!                                 cm^2 for temperature T for each
!                                 interval wavelength.
!
!     A subroutine which calculates the ozone absorption cross section
!     based on John E. Frederick (1985). The temperature dependent cross
!     section data set is that of A. M. Bass of the National Bureau of
!     Standards provided by R. D. McPeters.
!
!     The temperature dependence is in the 3rd significant figure between
!     263.158 - 266.167 nm. (Intervals 84-102).
!
!     The temperature range covered is 203 to 298 K.
!
!     The fit used is a quadratic fit of the form;
!     sigma(O3,t)={C0(i)+C1(i)(T-230)+C2(i)(T-230)^2}10^-n
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: JPWAV
      REAL, INTENT(IN)    :: T
      REAL, INTENT(INOUT) :: AO3(JPWAV)

! Local varaiables
!     Temperature in kelvin.
      REAL :: TC
!
!     Local variables & Polynomial coefficients.
      REAL :: TM230
      REAL :: TM2302

      REAL :: C(3,19)
      REAL :: N(19)
!
      DATA C/9.6312E+0 , 1.1875E-3 , -1.7386E-5 , 8.3211E+0 ,           &
           3.6495E-4 , 2.4691E-6 , 6.8810E+0 , 2.4598E-4 , 1.1692E-5 ,  &
           5.3744E+0 , 1.0325E-3 , 1.2573E-6 , 3.9575E+0 , 1.6851E-3 ,  &
          -6.8648E-6 , 2.7095E+0 , 1.4502E-3 ,-2.8925E-6 , 1.7464E+0 ,  &
           8.9350E-4 , 3.5914E-6 , 1.0574E+0 , 7.8270E-4 , 2.0024E-6 ,  &
           5.9574E+0 , 4.9448E-3 , 3.6589E-5 , 3.2348E+0 , 3.5392E-3 ,  &
           2.4769E-5 , 1.7164E+0 , 2.4542E-3 , 1.6913E-5 , 8.9612E+0 ,  &
           1.4121E-2 , 1.2498E-4 , 4.5004E+0 , 8.4327E-3 , 7.8903E-5 ,  &
           2.1866E+0 , 4.8343E-3 , 5.1970E-5 , 1.0071E+1 , 3.3409E-2 ,  &
           2.6621E-4 , 5.0848E+0 , 1.8178E-2 , 1.6301E-4 , 2.1233E+0 ,  &
           8.8453E-3 , 1.2633E-4 , 8.2861E+0 , 4.2692E-2 , 8.7057E-4 ,  &
           2.9415E+0 , 5.3051E-2 , 3.4964E-4/
!
      DATA N/8*18 , 3*19 , 3*20 , 3*21 , 2*22/
!
!
!     Check that temperature is in range.
!
      TC = T
      IF ( TC > 298.0 ) TC = 298.0
      IF ( TC < 203.0 ) TC = 203.0
!
      TM230 = TC - 230.0
      TM2302 = TM230*TM230
!
      AO3(84:102) = (C(1,1:19)+C(2,1:19)*TM230+C(3,1:19)*TM2302)        &
                   *(10.**(-N(1:19)))
!
      END SUBROUTINE ACSO3

      SUBROUTINE ACSO3W(JW,T,JPWAV,AO3)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     As ACSO3 but for a single wavelength interval JW
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: JW
      INTEGER, INTENT(IN) :: JPWAV
      REAL, INTENT(IN)    :: T
      REAL, INTENT(INOUT) :: AO3(JPWAV)

! Local variables
      INTEGER :: JJ
!
!     Temperature in kelvin.
      REAL :: TC
!
!     Local variables & Polynomial coefficients.
      REAL :: TM230
      REAL :: TM2302
      REAL :: C(3,19)
      REAL :: N(19)
!
      DATA C/9.6312E+0 , 1.1875E-3 , -1.7386E-5 , 8.3211E+0 ,           &
           3.6495E-4 , 2.4691E-6 , 6.8810E+0 , 2.4598E-4 , 1.1692E-5 ,  &
           5.3744E+0 , 1.0325E-3 , 1.2573E-6 , 3.9575E+0 , 1.6851E-3 ,  &
          -6.8648E-6 , 2.7095E+0 , 1.4502E-3 ,-2.8925E-6 , 1.7464E+0 ,  &
           8.9350E-4 , 3.5914E-6 , 1.0574E+0 , 7.8270E-4 , 2.0024E-6 ,  &
           5.9574E+0 , 4.9448E-3 , 3.6589E-5 , 3.2348E+0 , 3.5392E-3 ,  &
           2.4769E-5 , 1.7164E+0 , 2.4542E-3 , 1.6913E-5 , 8.9612E+0 ,  &
           1.4121E-2 , 1.2498E-4 , 4.5004E+0 , 8.4327E-3 , 7.8903E-5 ,  &
           2.1866E+0 , 4.8343E-3 , 5.1970E-5 , 1.0071E+1 , 3.3409E-2 ,  &
           2.6621E-4 , 5.0848E+0 , 1.8178E-2 , 1.6301E-4 , 2.1233E+0 ,  &
           8.8453E-3 , 1.2633E-4 , 8.2861E+0 , 4.2692E-2 , 8.7057E-4 ,  &
           2.9415E+0 , 5.3051E-2 , 3.4964E-4/
!
      DATA N/8*18 , 3*19 , 3*20 , 3*21 , 2*22/
!
!
!     Check if the calculation is required.
      IF ( (JW >= 84) .AND. (JW <= 102) ) THEN
!
!        Check that temperature is in range.
         TC = T
         IF ( TC > 298.0 ) TC = 298.0
         IF ( TC < 203.0 ) TC = 203.0
!
         TM230 = TC - 230.0
         TM2302 = TM230*TM230
!
         JJ = JW - 83
         AO3(JW) = (C(1,JJ)+C(2,JJ)*TM230+C(3,JJ)*TM2302)               &
                   *(10.**(-N(JJ)))
!
      END IF
!
      END SUBROUTINE ACSO3W

      SUBROUTINE ACSSR_OLD(TC,JPWAV,AO2SR)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name     Type               Description.
!     TC       Real               Total slant path O2 column above the
!                                 point in molecules per cm^2.
!
!     AO2SR    Array of real      Absorption cross-section of O2 in
!                                 cm^2 for each wavelength interval
!                                 in the Schumann-Runge bands.
!
!     A subroutine which calculates the SR absorption cross section for
!     wavelength interval JW and slantpath column TC.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: JPWAV
      REAL, INTENT(IN)    :: TC      ! Total O2 column
!     O2 cross-section in the Schumann Runge bands:
      REAL, INTENT(INOUT) :: AO2SR(JPWAV)

! Local variables
      INTEGER, PARAMETER :: JPWAVESR=17
!
      REAL :: LOGC
      REAL :: LOGC2
!
      REAL :: SR0(JPWAVESR)
      REAL :: SR1(JPWAVESR)
      REAL :: SR2(JPWAVESR)
!
!     Absorption cross sections.
!
!     Schumann Runge Cross section parameterisation (Frederick 1985)
!     WMO intervals (46 to 62)
!
      DATA SR0/8.587E1 , 2.670E1 , 1.858E1 , 4.791E1 , 3.796E1 ,        &
           5.381E1 , 5.440E1 , 5.467E1 , 5.941E1 , 8.666E1 , 1.238E2 ,  &
           1.220E2 , 1.229E2 , 1.345E2 , 4.000E1 , 4.769E1 , 1.965E1/
!
      DATA SR1/ - 2.239 , 4.390E-1 , 8.151E-1 , -4.853E-1 , -5.609E-2 , &
           -6.840E-1 , -6.537E-1 , -6.389E-1 , -8.120E-1 , -1.838 ,     &
           -3.403 , -3.196 , -3.120 , -3.390 , 4.420E-1 , 2.376E-1 ,    &
           1.353/
!
      DATA SR2/2.922E-2 , -4.684E-4 , -4.577E-3 , 9.756E-3 , 5.366E-3 , &
           1.173E-2 , 1.097E-2 , 1.081E-2 , 1.258E-2 , 2.231E-2 ,       &
           3.868E-2 , 3.560E-2 , 3.418E-2 , 3.548E-2 , -2.869E-3 ,      &
           -1.604E-3 , -1.253E-2/
!
!
!     Natural LOG of O2 column used in the Frederick parameterisation
!     of O2 cross section.
      LOGC  = LOG(TC)
      LOGC2 = LOGC*LOGC
!
      AO2SR(46:62) =                                                    &
        EXP(-(SR0(1:17) + SR1(1:17)*LOGC + SR2(1:17)*LOGC2))
! Olaf check
!     &  EXP(-(SR0(1:17) + SR1(1:17)*LOGC + SR2(1:17)*LOGC2ARG))
!
      END SUBROUTINE ACSSR_OLD

      SUBROUTINE ACSSR(TC,JPWAV,AO2SR)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name     Type               Description.
!     TC       Real               Total slant path O2 column above the
!                                 point in molecules per cm^2.
!
!     AO2SR    Array of real      Absorption cross-section of O2 in
!                                 cm^2 for each wavelength interval
!                                 in the Schumann-Runge bands.
!
!     O2 cross section in Schumann-Runge region. Origin of data: Table
!     7.7 in WMO (1985) gives transmission in the Schumann-Runge window
!     (175-205 nm). d ln T/ d (column oxygen) is displayed here, which
!     is the O2 cross section.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!  v1.2 Parameterization replaced with numbers derived transmission
!       data from WMO (1985)           Olaf Morgenstern 10/1/2006
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE

! Subroutine interface
!     Total O2 column.
      REAL, INTENT(IN)    :: tc
!     Number of wavelength intervals
      INTEGER, INTENT(IN) :: jpwav
!     O2 cross section in Schumann-Runge interval
      REAL, INTENT(INOUT) :: AO2SR(JPWAV)

! Local parameters and variables
!     Number of Schumann-Runge wavelength intervals.
      INTEGER, PARAMETER :: JPWAVESR=17
! Number of O2 column intervals
      INTEGER, PARAMETER :: NO2COL = 20
!
!     Oxygen column intervol
      INTEGER :: jo2
      REAL :: frac
      LOGICAL, SAVE :: first = .TRUE.
!
! O2 cross section in Schumann-Runge bands
      REAL, SAVE :: sr(JPWAVESR, no2col)
      REAL, SAVE :: logsr(JPWAVESR, no2col)
      REAL, SAVE :: o2col(no2col)
      REAL, SAVE :: logo2col(no2col)

! Main block: Initialize data upon first entry
      IF (first) THEN
        sr(:,1) = (/                                                    &
      2.0596E-19,1.3273E-19,6.6113E-20,5.1048E-20,4.5025E-20,1.4983E-20,&
      8.9858E-21,2.9953E-21,2.9932E-21,1.0000E-40,1.0000E-40,1.0000E-40,&
      1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40 /)
        sr(:,2) = (/                                                    &
      1.9882E-19,1.1754E-19,5.8335E-20,5.0977E-20,4.0015E-20,1.4504E-20,&
      6.0382E-21,3.6212E-21,1.2071E-21,1.0000E-40,1.0000E-40,1.0000E-40,&
      1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40 /)
        sr(:,3) = (/                                                    &
      1.8704E-19,9.5822E-20,5.1962E-20,4.8400E-20,3.6629E-20,1.2836E-20,&
      6.4058E-21,2.9866E-21,1.2793E-21,8.5257E-22,4.2614E-22,1.0000E-40,&
      1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40 /)
        sr(:,4) = (/                                                    &
      1.6133E-19,6.4581E-20,4.2558E-20,4.2452E-20,2.9852E-20,1.1903E-20,&
      6.1489E-21,2.8427E-21,1.1954E-21,4.4801E-22,5.9735E-22,1.4923E-22,&
      1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40 /)
        sr(:,5) = (/                                                    &
      1.2039E-19,3.7382E-20,2.9637E-20,3.1606E-20,2.0773E-20,1.0706E-20,&
      5.7048E-21,2.7256E-21,1.1954E-21,5.4263E-22,5.4263E-22,1.0845E-22,&
      1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40 /)
        sr(:,6) = (/                                                    &
      8.1089E-20,2.2320E-20,1.7193E-20,2.0043E-20,1.2370E-20,8.3474E-21,&
      4.9413E-21,2.4020E-21,1.1291E-21,4.5832E-22,5.2091E-22,1.4555E-22,&
      4.1557E-23,1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40 /)
        sr(:,7) = (/                                                    &
      5.2569E-20,1.4185E-20,9.5606E-21,1.1560E-20,7.1082E-21,5.4135E-21,&
      3.7658E-21,1.8750E-21,1.0425E-21,4.5744E-22,5.1708E-22,1.2648E-22,&
      3.3682E-23,8.4137E-24,1.0000E-40,1.0000E-40,1.0000E-40 /)
        sr(:,8) = (/                                                    &
      3.3180E-20,9.4479E-21,5.6342E-21,6.7482E-21,4.2844E-21,3.1560E-21,&
      2.4997E-21,1.3218E-21,8.5520E-22,4.1620E-22,4.9829E-22,1.3120E-22,&
      3.6328E-23,1.0000E-40,3.6300E-24,3.6300E-24,3.6274E-24  /)
        sr(:,9) = (/                                                    &
      2.1783E-20,6.5782E-21,3.6155E-21,4.2548E-21,2.7253E-21,1.9454E-21,&
      1.5785E-21,9.1690E-22,6.2044E-22,3.6428E-22,4.5727E-22,1.2660E-22,&
      3.4726E-23,3.2997E-24,1.6499E-24,1.0000E-40,1.0000E-40  /)
        sr(:,10) = (/                                                   &
      1.6102E-20,4.7863E-21,2.5785E-21,2.8755E-21,1.8324E-21,1.2969E-21,&
      1.0247E-21,6.3919E-22,4.0578E-22,2.7130E-22,3.8879E-22,1.1894E-22,&
      3.2220E-23,4.6906E-24,1.5631E-24,1.0000E-40,1.0000E-40  /)
        sr(:,11) = (/                                                   &
      1.4089E-20,3.6881E-21,2.0761E-21,3.7937E-21,1.3054E-21,9.1906E-22,&
      7.2141E-22,4.5887E-22,2.7029E-22,2.1300E-22,3.0038E-22,1.0554E-22,&
      2.9401E-23,4.2685E-24,3.8789E-24,3.8765E-25,7.7525E-25  /)
        sr(:,12) = (/                                                   &
      1.4608E-20,3.1322E-21,1.8465E-21,8.0188E-22,8.9807E-22,6.8966E-22,&
      5.5373E-22,3.4768E-22,1.9326E-22,1.6128E-22,2.1394E-22,8.7938E-23,&
      2.4211E-23,5.0370E-24,1.4084E-24,6.0272E-25,2.0093E-25  /)
        sr(:,13) = (/                                                   &
      1.4608E-20,2.9988E-21,1.7814E-21,1.3742E-21,7.9193E-22,5.3877E-22,&
      4.6246E-22,2.8948E-22,1.4895E-22,1.2973E-22,1.5169E-22,6.8160E-23,&
      2.0371E-23,5.7585E-24,2.4894E-24,6.4781E-25,5.3968E-25  /)
        sr(:,14) = (/                                                   &
      1.4608E-20,3.1296E-21,1.7853E-21,1.1892E-21,6.3662E-22,4.2117E-22,&
      3.8823E-22,1.8296E-22,1.1938E-22,1.0509E-22,1.0961E-22,5.1871E-23,&
      1.7125E-23,6.7396E-24,2.6245E-24,7.7167E-25,4.1525E-25  /)
        sr(:,15) = (/                                                   &
      1.4608E-20,2.9777E-21,1.5883E-21,9.6037E-22,5.0665E-22,3.0382E-22,&
      2.7925E-22,1.3249E-22,8.6811E-23,7.6305E-23,7.3698E-23,3.7378E-23,&
      1.3751E-23,5.5964E-24,2.6709E-24,7.9697E-25,4.1381E-25  /)
        sr(:,16) = (/                                                   &
      1.4608E-20,2.9777E-21,1.2210E-21,7.0899E-22,3.6453E-22,2.0328E-22,&
      1.7512E-22,8.4507E-23,5.6116E-23,4.8982E-23,4.3889E-23,2.4642E-23,&
      9.6018E-24,4.3275E-24,2.4343E-24,6.6242E-25,4.1886E-25  /)
        sr(:,17) = (/                                                   &
      1.4608E-20,2.9777E-21,1.2210E-21,5.5897E-22,2.7411E-22,1.4980E-22,&
      1.1815E-22,5.8790E-23,3.6523E-23,3.2044E-23,2.5235E-23,1.5842E-23,&
      6.5220E-24,3.2660E-24,2.2322E-24,6.0985E-25,4.1236E-25  /)
        sr(:,18) = (/                                                   &
      1.4608E-20,2.9777E-21,1.2210E-21,5.5897E-22,2.3540E-22,1.2884E-22,&
      1.7942E-22,4.6834E-23,2.6335E-23,2.2983E-23,1.5504E-23,1.0513E-23,&
      4.5950E-24,2.6892E-24,2.1383E-24,5.7770E-25,3.9394E-25  /)
        sr(:,19) = (/                                                   &
      1.4608E-20,2.9777E-21,1.2210E-21,5.5897E-22,2.3540E-22,1.1905E-22,&
      5.0663E-23,4.0295E-23,2.1386E-23,1.7647E-23,1.0693E-23,7.3588E-24,&
      3.3738E-24,2.3488E-24,2.0769E-24,5.6322E-25,3.9207E-25  /)
        sr(:,20) = (/                                                   &
      1.4608E-20,2.9777E-21,1.2210E-21,5.5897E-22,2.3540E-22,1.1905E-22,&
      5.0663E-23,3.5399E-23,1.8521E-23,1.4160E-23,8.0985E-24,5.4163E-24,&
      2.5977E-24,2.1190E-24,2.0413E-24,5.4415E-25,3.8368E-25  /)

        logsr = LOG(sr)

        o2col = (/                                                      &
      5.3368E+16,1.0627E+17,2.4629E+17,6.4304E+17,1.7548E+18,4.7351E+18,&
      1.2299E+19,3.0403E+19,7.1301E+19,1.5943E+20,3.4126E+20,6.9993E+20,&
      1.3797E+21,2.6309E+21,4.9365E+21,9.3681E+21,1.8360E+22,3.7372E+22,&
      7.8501E+22,1.6851E+23 /)

        logo2col = LOG(o2col)
        first = .FALSE.
      END IF

! Perform linear interpolation in log-log space of cross section in
! Schumann-Runge window. In case the O2 column is outside the window
! covered, take

      IF (tc < o2col(1)) THEN
        jo2 = 1
        frac = 0.
      ELSEIF (tc  >=  o2col(no2col)) THEN
        jo2 = no2col-1
        frac = 1.
      ELSE
        jo2 = no2col-1
        DO WHILE (o2col(jo2) > tc)
          jo2 = jo2 - 1
        END DO
        frac = (LOG(tc)          - logo2col(jo2))/                      & 
               (logo2col(jo2 + 1) - logo2col(jo2))
      END IF

      AO2SR(46:45 + jpwavesr) =                                         &
        EXP((1. - frac) * logsr(:,jo2) + frac * logsr(:,jo2 + 1))
!
      END SUBROUTINE ACSSR

      SUBROUTINE ACSSRW_OLD(TC,JW,JPWAV,AO2SR)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     As ACSSR but for a single wavelength interval JW
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!  v1.2 Made redundant                 Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: JPWAV
!     Wavelength interval.
      INTEGER, INTENT(IN) :: JW
!     Total O2 column.
      REAL, INTENT(IN)    :: TC
!     O2 cross-section in the Schumann Runge bands.
      REAL, INTENT(INOUT) :: AO2SR(JPWAV)

! Local variables
!     Number of Schumann-Runge wavelength intervals.
      INTEGER, PARAMETER :: JPWAVESR=17

      INTEGER :: J

      REAL :: ARG
      REAL :: LOGC
      REAL :: LOGC2
!
      REAL :: SR0(JPWAVESR)
      REAL :: SR1(JPWAVESR)
      REAL :: SR2(JPWAVESR)
!
!     Absorption cross sections.
!
!     Schumann Runge Cross section parameterisation (Frederick 1985)
!     WMO intervals (46 to 62)
!
      DATA SR0/8.587E1 , 2.670E1 , 1.858E1 , 4.791E1 , 3.796E1 ,        &
           5.381E1 , 5.440E1 , 5.467E1 , 5.941E1 , 8.666E1 , 1.238E2 ,  &
           1.220E2 , 1.229E2 , 1.345E2 , 4.000E1 , 4.769E1 , 1.965E1/
!
      DATA SR1/ - 2.239 , 4.390E-1 , 8.151E-1 , -4.853E-1 , -5.609E-2 , &
           -6.840E-1 , -6.537E-1 , -6.389E-1 , -8.120E-1 , -1.838 ,     &
           -3.403 , -3.196 , -3.120 , -3.390 , 4.420E-1 , 2.376E-1 ,    &
           1.353/
!
      DATA SR2/2.922E-2 , -4.684E-4 , -4.577E-3 , 9.756E-3 , 5.366E-3 , &
           1.173E-2 , 1.097E-2 , 1.081E-2 , 1.258E-2 , 2.231E-2 ,       &
           3.868E-2 , 3.560E-2 , 3.418E-2 , 3.548E-2 , -2.869E-3 ,      &
           -1.604E-3 , -1.253E-2/
!
!
      IF ( (JW >= 46) .AND. (JW <= 62) ) THEN
!
!        Natural LOG of O2 column used in the Frederick parameterisation
!        of O2 cross section.
         LOGC = LOG(TC)
         LOGC2 = LOGC*LOGC
!
         J = JW - 45
         ARG = SR0(J) + SR1(J)*LOGC + SR2(J)*LOGC2
         AO2SR(JW) = EXP(-ARG)
!
      END IF
!
      END SUBROUTINE ACSSRW_OLD

      SUBROUTINE ACSSRW(TC,JW,JPWAV,AO2SR)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name     Type               Description.
!     TC       Real               Total slant path O2 column above the
!                                 point in molecules per cm^2.
!
!     AO2SR    Array of real      Absorption cross-section of O2 in
!                                 cm^2 for each wavelength interval
!                                 in the Schumann-Runge bands.
!
!     A subroutine which calculates the SR absorption cross section for
!     wavelength interval JW and slantpath column TC.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!  v1.2 Parameterization replaced with numbers derived transmission
!       data from WMO (1985)           Olaf Morgenstern 10/1/2006
!  v1.3 Changes made for single wavenumber use.
!                                      Olaf Morgenstern 10/3/2006
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE

! Subroutine interface
!     Total O2 column.
      REAL, INTENT(IN)    :: tc
!     Number of wavelength intervals
      INTEGER, INTENT(IN) :: jpwav
      INTEGER, INTENT(IN) :: jw
!     O2 cross section in Schumann-Runge interval
      REAL, INTENT(INOUT) :: AO2SR(JPWAV)

! Local parameters and variables
!     Number of Schumann-Runge wavelength intervals.
      INTEGER, PARAMETER :: JPWAVESR=17
! Number of O2 column intervals
      INTEGER, PARAMETER :: NO2COL = 20
!
!     Oxygen column intervol
      INTEGER :: jo2
      REAL :: frac

      LOGICAL, SAVE :: first = .TRUE.

      INTEGER :: j
!
! O2 cross section in Schumann-Runge bands
      REAL, SAVE :: sr(JPWAVESR, no2col)
      REAL, SAVE :: logsr(JPWAVESR, no2col)
      REAL, SAVE :: o2col(no2col)
      REAL, SAVE :: logo2col(no2col)

! Main block: Initialize data upon first entry
      IF (first) THEN
        sr(:,1) = (/                                                    &
      2.0596E-19,1.3273E-19,6.6113E-20,5.1048E-20,4.5025E-20,1.4983E-20,&
      8.9858E-21,2.9953E-21,2.9932E-21,1.0000E-40,1.0000E-40,1.0000E-40,&
      1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40 /)
        sr(:,2) = (/                                                    &
      1.9882E-19,1.1754E-19,5.8335E-20,5.0977E-20,4.0015E-20,1.4504E-20,&
      6.0382E-21,3.6212E-21,1.2071E-21,1.0000E-40,1.0000E-40,1.0000E-40,&
      1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40 /)
        sr(:,3) = (/                                                    &
      1.8704E-19,9.5822E-20,5.1962E-20,4.8400E-20,3.6629E-20,1.2836E-20,&
      6.4058E-21,2.9866E-21,1.2793E-21,8.5257E-22,4.2614E-22,1.0000E-40,&
      1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40 /)
        sr(:,4) = (/                                                    &
      1.6133E-19,6.4581E-20,4.2558E-20,4.2452E-20,2.9852E-20,1.1903E-20,&
      6.1489E-21,2.8427E-21,1.1954E-21,4.4801E-22,5.9735E-22,1.4923E-22,&
      1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40 /)
        sr(:,5) = (/                                                    &
      1.2039E-19,3.7382E-20,2.9637E-20,3.1606E-20,2.0773E-20,1.0706E-20,&
      5.7048E-21,2.7256E-21,1.1954E-21,5.4263E-22,5.4263E-22,1.0845E-22,&
      1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40 /)
        sr(:,6) = (/                                                    &
      8.1089E-20,2.2320E-20,1.7193E-20,2.0043E-20,1.2370E-20,8.3474E-21,&
      4.9413E-21,2.4020E-21,1.1291E-21,4.5832E-22,5.2091E-22,1.4555E-22,&
      4.1557E-23,1.0000E-40,1.0000E-40,1.0000E-40,1.0000E-40 /)
        sr(:,7) = (/                                                    &
      5.2569E-20,1.4185E-20,9.5606E-21,1.1560E-20,7.1082E-21,5.4135E-21,&
      3.7658E-21,1.8750E-21,1.0425E-21,4.5744E-22,5.1708E-22,1.2648E-22,&
      3.3682E-23,8.4137E-24,1.0000E-40,1.0000E-40,1.0000E-40 /)
        sr(:,8) = (/                                                    &
      3.3180E-20,9.4479E-21,5.6342E-21,6.7482E-21,4.2844E-21,3.1560E-21,&
      2.4997E-21,1.3218E-21,8.5520E-22,4.1620E-22,4.9829E-22,1.3120E-22,&
      3.6328E-23,1.0000E-40,3.6300E-24,3.6300E-24,3.6274E-24  /)
        sr(:,9) = (/                                                    &
      2.1783E-20,6.5782E-21,3.6155E-21,4.2548E-21,2.7253E-21,1.9454E-21,&
      1.5785E-21,9.1690E-22,6.2044E-22,3.6428E-22,4.5727E-22,1.2660E-22,&
      3.4726E-23,3.2997E-24,1.6499E-24,1.0000E-40,1.0000E-40  /)
        sr(:,10) = (/                                                   &
      1.6102E-20,4.7863E-21,2.5785E-21,2.8755E-21,1.8324E-21,1.2969E-21,&
      1.0247E-21,6.3919E-22,4.0578E-22,2.7130E-22,3.8879E-22,1.1894E-22,&
      3.2220E-23,4.6906E-24,1.5631E-24,1.0000E-40,1.0000E-40  /)
        sr(:,11) = (/                                                   &
      1.4089E-20,3.6881E-21,2.0761E-21,3.7937E-21,1.3054E-21,9.1906E-22,&
      7.2141E-22,4.5887E-22,2.7029E-22,2.1300E-22,3.0038E-22,1.0554E-22,&
      2.9401E-23,4.2685E-24,3.8789E-24,3.8765E-25,7.7525E-25  /)
        sr(:,12) = (/                                                   &
      1.4608E-20,3.1322E-21,1.8465E-21,8.0188E-22,8.9807E-22,6.8966E-22,&
      5.5373E-22,3.4768E-22,1.9326E-22,1.6128E-22,2.1394E-22,8.7938E-23,&
      2.4211E-23,5.0370E-24,1.4084E-24,6.0272E-25,2.0093E-25  /)
        sr(:,13) = (/                                                   &
      1.4608E-20,2.9988E-21,1.7814E-21,1.3742E-21,7.9193E-22,5.3877E-22,&
      4.6246E-22,2.8948E-22,1.4895E-22,1.2973E-22,1.5169E-22,6.8160E-23,&
      2.0371E-23,5.7585E-24,2.4894E-24,6.4781E-25,5.3968E-25  /)
        sr(:,14) = (/                                                   &
      1.4608E-20,3.1296E-21,1.7853E-21,1.1892E-21,6.3662E-22,4.2117E-22,&
      3.8823E-22,1.8296E-22,1.1938E-22,1.0509E-22,1.0961E-22,5.1871E-23,&
      1.7125E-23,6.7396E-24,2.6245E-24,7.7167E-25,4.1525E-25  /)
        sr(:,15) = (/                                                   &
      1.4608E-20,2.9777E-21,1.5883E-21,9.6037E-22,5.0665E-22,3.0382E-22,&
      2.7925E-22,1.3249E-22,8.6811E-23,7.6305E-23,7.3698E-23,3.7378E-23,&
      1.3751E-23,5.5964E-24,2.6709E-24,7.9697E-25,4.1381E-25  /)
        sr(:,16) = (/                                                   &
      1.4608E-20,2.9777E-21,1.2210E-21,7.0899E-22,3.6453E-22,2.0328E-22,&
      1.7512E-22,8.4507E-23,5.6116E-23,4.8982E-23,4.3889E-23,2.4642E-23,&
      9.6018E-24,4.3275E-24,2.4343E-24,6.6242E-25,4.1886E-25  /)
        sr(:,17) = (/                                                   &
      1.4608E-20,2.9777E-21,1.2210E-21,5.5897E-22,2.7411E-22,1.4980E-22,&
      1.1815E-22,5.8790E-23,3.6523E-23,3.2044E-23,2.5235E-23,1.5842E-23,&
      6.5220E-24,3.2660E-24,2.2322E-24,6.0985E-25,4.1236E-25  /)
        sr(:,18) = (/                                                   &
      1.4608E-20,2.9777E-21,1.2210E-21,5.5897E-22,2.3540E-22,1.2884E-22,&
      1.7942E-22,4.6834E-23,2.6335E-23,2.2983E-23,1.5504E-23,1.0513E-23,&
      4.5950E-24,2.6892E-24,2.1383E-24,5.7770E-25,3.9394E-25  /)
        sr(:,19) = (/                                                   &
      1.4608E-20,2.9777E-21,1.2210E-21,5.5897E-22,2.3540E-22,1.1905E-22,&
      5.0663E-23,4.0295E-23,2.1386E-23,1.7647E-23,1.0693E-23,7.3588E-24,&
      3.3738E-24,2.3488E-24,2.0769E-24,5.6322E-25,3.9207E-25  /)
        sr(:,20) = (/                                                   &
      1.4608E-20,2.9777E-21,1.2210E-21,5.5897E-22,2.3540E-22,1.1905E-22,&
      5.0663E-23,3.5399E-23,1.8521E-23,1.4160E-23,8.0985E-24,5.4163E-24,&
      2.5977E-24,2.1190E-24,2.0413E-24,5.4415E-25,3.8368E-25  /)

        logsr = LOG(sr)

        o2col = (/                                                      &
      5.3368E+16,1.0627E+17,2.4629E+17,6.4304E+17,1.7548E+18,4.7351E+18,&
      1.2299E+19,3.0403E+19,7.1301E+19,1.5943E+20,3.4126E+20,6.9993E+20,&
      1.3797E+21,2.6309E+21,4.9365E+21,9.3681E+21,1.8360E+22,3.7372E+22,&
      7.8501E+22,1.6851E+23 /)

        logo2col = LOG(o2col)
        first = .FALSE.
      END IF

      IF ( (JW >= 46) .AND. (JW <= 62) ) THEN
! Perform linear interpolation in log-log space of cross section in
! Schumann-Runge window. In case the O2 column is outside the window
! covered, take

      IF (tc < o2col(1)) THEN
        jo2 = 1
        frac = 0.
      ELSEIF (tc  >=  o2col(no2col)) THEN
        jo2 = no2col-1
        frac = 1.
      ELSE
        jo2 = no2col-1
        DO WHILE (o2col(jo2) > tc)
          jo2 = jo2 - 1
        END DO
        frac = (LOG(tc)          - logo2col(jo2))/                      &
               (logo2col(jo2 + 1) - logo2col(jo2))
      END IF

      j = jw - 45
      AO2SR(jw) =                                                       &
        EXP((1. - frac) * logsr(j,jo2) + frac * logsr(j,jo2 + 1))
!
      END IF
      END SUBROUTINE ACSSRW

      SUBROUTINE ACSCOS(T,WAVENM,AOCS)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Calculate T-dependent COS cross sections, following JPL (2002)
!
!     Name     Type               Description.
!     T        Real               Temperature in kelvin.
!
!     WAVENM   Array of real      Wavelength of each interval in nm.
!
!     AOCS     Array of real      Absorption cross-section of COS in
!                                 cm^2 for temperature T for each
!                                 interval wavelength.
!
!     The expression is valid for the wavelength range; 186.1-296.3 nm
!                            and the temperature range; 225-295 K.
!
!  v1.1 17/1/2006  Original code   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
#include "parpho.h"
!
! Subroutine interface
      REAL, INTENT(IN)  :: T
      REAL, INTENT(IN)  :: WAVENM(JPWAV)
      REAL, INTENT(OUT) :: AOCS(JPWAV)

! Local variables
      REAL :: ARG
      REAL :: TC
      REAL, SAVE :: ACOS295(JPWAV)
      REAL, SAVE :: ACOS225(JPWAV)
!
!     COS: JPL 2002 T=295K values.
      DATA ACOS295/                                                     &
        51*0.0,                                                         &
        18.9, 8.33, 3.75, 2.21, 1.79, 1.94, 2.48, 3.30, 4.48, 6.12,     &
        8.19, 10.8, 14.1, 17.6, 21.8, 25.5, 28.2, 30.5, 31.9, 30.2,     &
        26.8, 22.1, 17.1, 12.5, 8.54, 5.61, 3.51, 2.11, 1.21, .674,     &
        .361, .193,.0941,.0486,.0248,.0119,.0584,.0264,.0012,.0005,     &
       .0002,                                                           &
       111*0.0/
!     COS: JPL 2002 T=225K values.
      DATA ACOS225/                                                     &
        51*0.0,                                                         &
        13.0, 5.63, 2.50, 1.61, 1.53, 1.84, 2.44, 3.30, 4.50, 6.17,     &
        8.27, 10.9, 14.2, 17.6, 21.8, 25.3, 27.7, 29.4, 29.5, 27.4,     &
        23.7, 18.8, 14.0, 9.72, 6.24, 3.89, 2.29, 1.29, .679, .353,     &
        .178,.0900,.0419,.0199,.0101,.0048,.0021,.0009,.0005,.0002,     &
       112*0.0/
!
!     Check that temperature is in range.
      TC = MAX(MIN (T, 295.0), 225.0)
!
      ARG = (TC - 225.)/70.
! Linear interpolation between the cross sections given.
      AOCS = 1.0e-20 *                                                  &
            (ACOS295 * ARG + ACOS225 * (1. - ARG))
!
      END SUBROUTINE ACSCOS

      FUNCTION EI1(X)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name     Type               Description.
!     X        Real               Argument of the first exponential
!                                 integral.
!
!     This function calculates the first exponential integral using a
!     series expansion.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Function interface
      REAL, INTENT(IN) :: X
      REAL :: EI1

! Local variables
      REAL, PARAMETER :: THRESH=1.0E-35
      REAL, PARAMETER :: XMAX=50.0
!
      REAL, PARAMETER :: A011= -0.57721566
      REAL, PARAMETER :: A012=  0.99999193
      REAL, PARAMETER :: A013= -0.24991055
      REAL, PARAMETER :: A014=  0.05519968
      REAL, PARAMETER :: A015= -0.00976004
      REAL, PARAMETER :: A016=  0.00107857
!
      REAL, PARAMETER :: A1101= 8.5733285301
      REAL, PARAMETER :: A1102=18.0590169520
      REAL, PARAMETER :: A1103= 8.6347608925
      REAL, PARAMETER :: A1104= 0.2677737343
!
      REAL, PARAMETER :: B1101= 9.5733223454
      REAL, PARAMETER :: B1102=25.6329561486
      REAL, PARAMETER :: B1103=21.0996530827
      REAL, PARAMETER :: B1104= 3.9584969228
!
      REAL, PARAMETER :: A10INF1=4.03640
      REAL, PARAMETER :: A10INF2=1.15198
      REAL, PARAMETER :: B10INF1=5.03637
      REAL, PARAMETER :: B10INF2=4.19160

      REAL :: TOP
      REAL :: BOT
      REAL :: X2
      REAL :: X3
      REAL :: X4
!
!
      IF ( X < 0.0 ) WRITE (6,*)                                        &
              '** EI1 WARNING ** Negative argument >',X
!
!     Check for minimum value
      IF ( X <= THRESH ) THEN
         EI1 = 86.9210129
      ELSEIF ( (X > THRESH) .AND. (X <= 1.0) ) THEN
         EI1 = A011                                                     &
             + A012* X                                                  &
             + A013*(X**(2))                                            &
             + A014*(X**(3))                                            &
             + A015*(X**(4))                                            &
             + A016*(X**(5))                                            &
             - LOG(X)
      ELSEIF ( (X > 1.0) .AND. (X <= 10.0) ) THEN
         X2 = X*X
         X3 = X2*X
         X4 = X3*X
         TOP = X4 + A1101*X3 + A1102*X2 + A1103*X + A1104
         BOT = X4 + B1101*X3 + B1102*X2 + B1103*X + B1104
         EI1 = TOP/(BOT*X*EXP(X))
      ELSEIF ( (X > 10.0) .AND. (X < XMAX) ) THEN
         X2 = X*X
         TOP = X2 + A10INF1*X + A10INF2
         BOT = X2 + B10INF1*X + B10INF2
         EI1 = TOP/(BOT*X*EXP(X))
      ELSEIF ( X >= XMAX ) THEN
         EI1 = 0.0
      END IF
!
      END FUNCTION EI1

      FUNCTION EI2(X)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name     Type               Description.
!     X        Real               Argument of the second exponential
!                                 integral.
!
!     This function calculates the second exponential integral of
!     the real number x.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Function interface
      REAL, INTENT(IN) :: X

! Local variables
      REAL, PARAMETER :: THRESH=0.10E-15
      REAL, PARAMETER :: XMAX=50.0
      REAL :: EI2
!
!
!     Calculate integral
      IF ( X <= THRESH ) THEN
         EI2 = 1.00
      ELSEIF ( X >= XMAX ) THEN
         EI2 = 0.0
      ELSE
         EI2 = EXP(-X) - X*EI1(X)
      END IF
!
      END FUNCTION EI2

      FUNCTION EI3(X)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name     Type               Description.
!     X        Real               Arguement of the third exponential
!                                 integral.
!
!     This function calculates the third exponential integral of
!     the real number x.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Function interface
      REAL, INTENT(IN) :: X

! Local variables
      REAL, PARAMETER :: THRESH=0.10E-7
      REAL, PARAMETER :: XMAX=50.0
      REAL :: EI3
!
!
!     Calculate integral.
      IF ( X <= THRESH ) THEN
         EI3 = 0.50
      ELSEIF ( X >= XMAX ) THEN
         EI3 = 0.0
      ELSE
         EI3 = (EXP(-X)-X*EI2(X))*0.5
      END IF
!
      END FUNCTION EI3

!  Calculate quantum yield of O(1D) from O3 photolysis.
!  Test version
! ######################################################################
!
! Subroutine Interface:
!
      SUBROUTINE QUANTO1D(T,JPWAV,WAVENM,QEO1D)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name     Type               Description.
!     T        Real               Temperature in Kelvin.
!
!     WAVENM   Array of real      Wavelength of each interval in nm.
!
!     QEO1D    Array of real      The quantum yield of O(1D) from O3
!                                 photolysis at wavelengths less than
!                                 310 nm.
!
!     A subroutine which calculates the quantum yield of O(1D) from O3
!     photolysis at wavelengths less than 310 nm.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: JPWAV
      REAL, INTENT(IN)    :: T
      REAL, INTENT(IN)    :: WAVENM(JPWAV)
      REAL, INTENT(INOUT) :: QEO1D(JPWAV)

! Local variables
      INTEGER :: JW
      REAL :: TC
      REAL :: A
      REAL :: ARG
      REAL :: B
      REAL :: C
      REAL :: TAU
      REAL :: TAU2
      REAL :: TAU3
!
!
!     Set values of PHI at values not covered by the parameterisation.
      QEO1D(1:91) = 0.9
!
      QEO1D(98:JPWAV) = 0.0
!
!     Intervals 92 to 97.
!     Set up parameter TAU used in the parameterisation.
!     [A Temperature parameter.]
      TC = T
      TAU = TC - 230.0
      TAU2 = TAU*TAU
      TAU3 = TAU2*TAU
!
!     Set up the quantum yield of O(1D) from O3 photolysis at
!     wavelengths less than 310 nm using the JPL (evaluation 8)
!     parameterisation due to Moorgat & Kudzus (1978).
!
      DO JW = 92 , 97
!
         A =  0.332 + 2.5650E-4*TAU + 1.152E-5*TAU2 + 2.3130E-8*TAU3
         B = -0.575 + 5.5900E-3*TAU - 1.439E-5*TAU2 - 3.2700E-8*TAU3
         C =  0.466 + 8.8830E-4*TAU - 3.546E-5*TAU2 + 3.5190E-7*TAU3
         ARG= 308.2 + 4.4871E-2*TAU + 6.938E-5*TAU2 - 2.5452E-6*TAU3
!
         QEO1D(JW) = A*ATAN(B*(WAVENM(JW)-ARG)) + C
         QEO1D(JW) = MIN(QEO1D(JW),0.9)
         QEO1D(JW) = MAX(QEO1D(JW),0.0)
!
      END DO
!
      END SUBROUTINE QUANTO1D

      SUBROUTINE QUANTO12(T,JPWAV,WAVENM,QEO1D)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name     Type               Description.
!     T        Real               Temperature in Kelvin.
!
!     WAVENM   Array of real      Wavelength of each interval in nm.
!
!     QEO1D    Array of real      The quantum yield of O(1D) from O3
!                                 photolysis
!
!     A subroutine which calculates the quantum yield of O(1D) from O3
!     photolysis based on parameteristn of Michelson et al GRL Oct 1994.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: JPWAV
      REAL, INTENT(IN)    :: T
      REAL, INTENT(IN)    :: WAVENM(JPWAV)
      REAL, INTENT(INOUT) :: QEO1D(JPWAV)

! Local variables
      REAL :: phifac
      REAL :: PHI1DA(5)
      REAL :: PHI1DB(5)
!
!     Data from Michelsen et al. 1994
      DATA PHI1DA/1.01,  2.45, 16.55,  13.82,   11.8/
      DATA PHI1DB/3.93, 294.2, 846.3, 1061.9, 1435.0/
!
!
!     Factor = hc (planck (Js) x speed of light (cms-1))
      PHIFAC = 6.626E-34*2.99E10
!
      QEO1D(1:85) = 0.9
!
      QEO1D(86:93) = 1.98 - (301.0/WAVENM(86:93))
!
      QEO1D(94:98) = PHI1DA(1:5)*EXP((-PHI1DB(1:5)*PHIFAC)           &
                                           /(1.38E-23*T))
!
      QEO1D(99:JPWAV) = 0.0
!
      END SUBROUTINE QUANTO12

      SUBROUTINE SCATCS(JPWAV,SCS,WAVENM)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name     Type               Description.
!
!     SCS      Array of real      The Rayleigh scattering cross
!                                 section for air molecules in
!                                 molecules per cm^2.
!
!     WAVENM   Array of real      Wavelength of each interval in nm.
!
!     Subroutine which calculates the Rayleigh scattering cross section
!     based on Nicolet (1984).
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: JPWAV
      REAL, INTENT(IN)    :: WAVENM(JPWAV)
      REAL, INTENT(INOUT) :: SCS(JPWAV)

! Local variables
      INTEGER :: JW
      REAL :: XLAMDA
      REAL :: CHI
!
!
      SCS(1:45) = 0.0
      DO JW = 46 , 143
         XLAMDA  = WAVENM(JW)*1.0E-3
         CHI     = 0.389*XLAMDA + (0.09426/XLAMDA) - 0.3228
         SCS(JW) = 4.02E-28*(XLAMDA**(-4.0 - CHI))
      END DO
      DO JW = 144 , JPWAV
         XLAMDA  = WAVENM(JW)*1.0E-3
         SCS(JW) = 4.02E-28*(XLAMDA**(-4.04))
      END DO
!
      END SUBROUTINE SCATCS

      SUBROUTINE SETTAB(ALT,ALTC,DO2C,DO3C,TEMPC,                       &
                 PRESC,DRSC,TABS,TABANG,ALBEDO,ZENMAX,LSCAT,            &
                 TSPO2,TABPRES,QUANTA,WAVENM,SCS,AO2,AO2SR,AO3,DALT)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Calculate photolysis factors
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
#include "parpho.h"
!
! Subroutine interface
      LOGICAL, INTENT(IN) :: LSCAT
      REAL, INTENT(IN) :: ALBEDO             !     Ground albedo.
      REAL, INTENT(IN) :: ALT(JPLEVP1)
      REAL, INTENT(IN) :: ALTC(JPLEV)
      REAL, INTENT(IN) :: DO2C(JPLEV)
      REAL, INTENT(IN) :: DO3C (JPLEV)
      REAL, INTENT(IN) :: TEMPC(JPLEV)
      REAL, INTENT(IN) :: PRESC  (JPLEV)
      REAL, INTENT(IN) :: DRSC(JPLEV)
      REAL, INTENT(OUT):: TABS (JPLEV,JPCHI,JPWAV)
      REAL, INTENT(IN) :: TABANG(JPCHI)
      REAL, INTENT(IN) :: ZENMAX(JPLEV)!  Maximum zenith angle in radians.
      REAL, INTENT(OUT):: TSPO2(JPLEV,JPCHI)
      REAL, INTENT(IN) :: TABPRES(JPLEV)
      REAL, INTENT(IN) :: QUANTA(JPWAV)
      REAL, INTENT(IN) :: WAVENM(JPWAV)
      REAL, INTENT(IN) :: SCS(JPWAV)
      REAL, INTENT(IN) :: AO2(JPWAV)
! Olaf check, ao2sr and ao3 are intent(inout) in CSO2O3
      REAL, INTENT(INOUT) :: AO2SR(JPWAV)
      REAL, INTENT(INOUT) :: AO3(JPWAV)
      REAL, INTENT(IN) :: DALT (JPLEV)

! Local variables
      INTEGER :: I
      INTEGER :: J
      INTEGER :: JC
      INTEGER :: JCI
      INTEGER :: JK
      INTEGER :: JN
      INTEGER :: JT
      INTEGER :: JW
      REAL :: ALPHA
      REAL :: ALTA
      REAL :: ALTACM
      REAL :: ALTCUR
      REAL :: ALTB
      REAL :: ALTBCM
      REAL :: ALT90
      REAL :: ARG
      REAL :: BETA
      REAL :: DTN
      REAL :: SPL
      REAL :: SPL2
      REAL :: TAUK
      REAL :: TAUN
      REAL :: TAUV
      REAL :: TEAK
      REAL :: TEAN
      REAL :: TEAV
      REAL :: TWOAMU
      REAL :: TEANP
      REAL :: TAUNP

      INTEGER :: IJT(JPS90,JPLEV)
      INTEGER :: INDX(JPLEV)

!     Slant path between centre of levels.
      REAL :: TABLEN(JPLEV,JPLEV,JPCHIN)
      REAL :: DEPTHA(JPLEV,JPCHI)
      REAL :: DEPTHS(JPLEV,JPCHI)
      REAL :: TABLENOA(JPLEV,JPLEV,JPS90)
      REAL :: TABLENOB(JPLEV,JPLEV,JPS90)
!
!
!     Vertical optical depths etc.
      REAL :: TEAC(JPLEVP1)
      REAL :: TAUC(JPLEVP1)
      REAL :: JAC (JPLEV)
      REAL :: TEA (JPLEVP1)
      REAL :: TAU (JPLEVP1)
!
!     Enhancement factor table.
      REAL :: TABS0(JPLEV,JPCHI,JPWAV)
!
      REAL :: B(JPLEV,JPLEV)
      REAL :: BINV (JPLEV,JPLEV)
      REAL :: BCOPY(JPLEV,JPLEV)
      REAL :: A(JPLEV,JPLEV)
      REAL :: DELTA(JPLEV,JPLEV)
!
!
!     Set up path lengths (cm) up to and including 90 degrees.
      DO JC = 1, JPCHIN
         DO I = 1, JPLEV
!
!           Zenith angle at the starting point.
            ALPHA = TABANG(JC)
!
!           Follow the path back in spherical geometry.
            DO J = I, JPLEV
!
!              If at the starting level
               IF(J == I)THEN
                 ALTA = ALTC(J)
               ELSE
!                For all other levels.
                 ALTA = ALT (J)
               END IF
!
               ALTB = ALT(J+1)
!
!              Calculate zenith angle for next level up.
               BETA = ASIN((ALTA+RE)*SIN(ALPHA)/(ALTB+RE))
!
!              Calculate the slant path within this level in cm
               ALTACM = (ALTA+RE)*1.0E5
               ALTBCM = (ALTB+RE)*1.0E5
               SPL2 = ALTACM*ALTACM + ALTBCM*ALTBCM -                   &
                   2.0*ALTACM*ALTBCM*COS(ALPHA-BETA)
               SPL  = SQRT(MAX(0.0,SPL2))
!
!              Distance in cm within level.
               TABLEN(I,J,JC) = SPL
!
!              Reinitialise variables.
               ALPHA = BETA
!
            END DO
         END DO
      END DO
!
!
!     Find the level index for each zenith angle and each level for
!     which the zenith angle is 90 degrees.
      DO JC = JPCHIN + 1, JPCHI
         JCI = JC - JPCHIN
         DO J = 1, JPLEV
            ALT90 = (ALTC(J)+RE)*SIN(TABANG(JC)) - RE
!           Altitudes monotonically increasing
            JT=ISRCHFGT(JPLEVP1,ALT,1,ALT90) - 1
            IJT(JCI,J) = JT
         END DO
      END DO
!
!
!     When the zenith angle TABANG(JC) > 90.0 degrees.
      DO JC = JPCHIN + 1, JPCHI
!
!        A zenith angle index.
         JCI = JC - JPCHIN
!
         DO J = 1, JPLEV
!
!           Initialise to zero.
            SPL = 0.0
            DO I = 1, JPLEV
               TABLENOA(I,J,JCI) = 0.0
               TABLENOB(I,J,JCI) = 0.0
            END DO
!
!           Level, JT, which contains tangent point.
            JT = IJT(JCI,J)
!
!           Check if light can get there.
            IF (JT > 0 .AND. J >= JT) THEN
!
!              Start by calculating path lengths just to the left of
!              the tangent point.
!
!              Check if this part of the calculation is needed.
               IF (J >= JT+1) THEN
!
!              Altitude where zenith angle is 90 degrees.
               ALT90 = (ALTC(J)+RE)*SIN(TABANG(JC))
!
!              Altitude of interface JT+1.
               ALTCUR = ALT(JT+1) + RE
!
!              Set zenith angle at current point.
               ALPHA = ASIN(ALT90/ALTCUR)
!
!              Go from current level
!              to the level above JT summing the path lengths.
               DO I = JT + 1, J - 1
!
                  ALTA = ALT(I  )
                  ALTB = ALT(I+1)
!
!                 Calc. zenith angle at the centre of the next level down.
                  BETA = ASIN((ALTA+RE)*SIN(ALPHA)/(ALTB+RE))
!
!                 Calculate the slant path between the centre of
!                 this level and the centre of the next level down in cm.
                  ALTACM = (ALTA+RE)*1.0E5
                  ALTBCM = (ALTB+RE)*1.0E5
                  SPL2 = ALTACM*ALTACM + ALTBCM*ALTBCM -                &
                     2.0*ALTACM*ALTBCM*COS(ALPHA-BETA)
                  SPL  = SQRT(MAX(0.0,SPL2))
!
                  TABLENOA(I,J,JCI) = SPL
!
!                 Reintialise variables.
!                 adjust zenith angle for spherical geometry.
                  ALPHA = BETA
!
               END DO
!
!              End of check if this part of the calculation is needed.
               END IF
!
!              Path length within level containing tangent point
!              Case (A) J >= JT +1
               IF (J >= JT+1) THEN
!
                 ALT90 = (ALTC(J)+RE)*SIN(TABANG(JC))
                 ALTCUR = ALT(JT+1) + RE
                 SPL = 2.0E5*SQRT(ALTCUR*ALTCUR - ALT90*ALT90)
!
!              Case (B) J = JT. Tangent height within level J
               ELSE IF (J == JT) THEN
!
                 ALT90 = (ALTC(J)+RE)*SIN(TABANG(JC))
!
                 ALTCUR = ALT(JT+1) + RE
                 SPL =       1.0E5*SQRT(ALTCUR*ALTCUR - ALT90*ALT90)
!
                 ALTCUR = ALTC(J) + RE
                 SPL = SPL + 1.0E5*SQRT(ALTCUR*ALTCUR - ALT90*ALT90)
!
               END IF
!
               TABLENOA(JT,J,JCI) = SPL
!
!              Now calculate path lengths to the right of the tangent point.
!
!              Set zenith angle at the current point.
               ALPHA = ASIN(ALT90/ALTCUR)
!
!              Now go from JT to the top of the atmosphere.
               DO I = JT + 1, JPLEV
!
!                 Initialise variables.
                  ALTA = ALT(I  )
                  ALTB = ALT(I+1)
!
!                 Calculate zenith angle for the next level up.
                  BETA = ASIN((ALTA+RE)*SIN(ALPHA)/(ALTB+RE))
!
!                 Calculate the slant path between this level and the
!                 next level down in cm.
                  ALTACM = (ALTA+RE)*1.0E5
                  ALTBCM = (ALTB+RE)*1.0E5
                  SPL2 = ALTACM*ALTACM + ALTBCM*ALTBCM -                &
                     2.0*ALTACM*ALTBCM*COS(ALPHA-BETA)
                  SPL  = SQRT(MAX(0.0,SPL2))
!
                  TABLENOB(I,J,JCI) = SPL
!
!                 Reintialise variables.
!                 Adjust zenith angle for spherical geometry.
                  ALPHA = BETA
!
               END DO
            END IF
         END DO
      END DO
!
!
!     Set up slant path O2 column. Used in AO2SR and ANO parameterisations.
!
!     Set up path lengths up to and including 90 degrees.
      DO JC = 1, JPCHIN
         DO I = 1, JPLEV
!
            TSPO2(I,JC) = 0.0
!
!           Follow the path back in spherical geometry.
            DO J = I, JPLEV
               TSPO2(I,JC) = TSPO2(I,JC) + TABLEN(I,J,JC)*DO2C(J)
            END DO
!
         END DO
      END DO
!
!     When the zenith angle TABANG(JC) > 90.0 degrees.
      DO JC = JPCHIN + 1, JPCHI
!
!        A zenith angle index.
         JCI = JC - JPCHIN
!
         DO J = 1, JPLEV
!
!           Set O2 column to small value to prevent LOG problems in ACSSRW
            TSPO2(J,JC) = 1.0E10
!
!           Find level, JT, just below where ALPHA=90 degrees.
            JT = IJT(JCI,J)
!
!           Check if light can get there.
            IF (JT > 0 .AND. J >= JT) THEN
!
!              Having found this level first go from the current level
!              down to the level above JT summing the path lengths.
               DO I = JT + 1, J - 1
                 TSPO2(J,JC)=TSPO2(J,JC) + TABLENOA(I,J,JCI)*DO2C(I)
               END DO
!
!              Now calculate the path length between the level
!              JT + 1 as it stradles the zenith angle of 90 degrees.
!
               TSPO2(J,JC) = TSPO2(J,JC) + TABLENOA(JT,J,JCI)*DO2C(JT)
!
!              Now calculate path lengths to the right of the tangent point.
!
!              Now go from JT to the top of the atmosphere.
               DO I = JT + 1, JPLEV
                  TSPO2(J,JC) = TSPO2(J,JC) + TABLENOB(I,J,JCI)*DO2C(I)
               END DO
!
            END IF
         END DO
      END DO
!
!     Set up optical depths
!
!     Wavelength loop.
      DO JW = JPLO, JPHI
!
!        Up to and including 90 degrees.
         DO JC = 1, JPCHIN
            DO I = 1, JPLEV
!
               DEPTHA(I,JC) = 0.0
               DEPTHS(I,JC) = 0.0
               DO J = I, JPLEV
!
!                 Reset the O2 and O3 absorption cross sections.
                  CALL CSO2O3(AO2SR,AO3,TEMPC(J),TSPO2(J,JC),           &
                              JW,JPWAV)
                  DEPTHA(I,JC) = DEPTHA(I,JC) + TABLEN(I,J,JC)*         &
                    (DO3C(J)*AO3(JW) + DO2C(J)*(AO2(JW)+AO2SR(JW)))
                  DEPTHS(I,JC) = DEPTHS(I,JC) + TABLEN(I,J,JC)*         &
                     DRSC(J)*SCS(JW)
               END DO
            END DO
         END DO
!
!
!        Greater than 90 degrees.
         DO JC = JPCHIN + 1, JPCHI
!
!           A zenith angle index.
            JCI = JC - JPCHIN
!
            DO J = 1, JPLEV
!
!              Initialise DEPTHA
               DEPTHA(J,JC) = 0.0
               DEPTHS(J,JC) = 0.0
!
!              Recall what JT was
               JT = IJT(JCI,J)
!
!              Having found this level first go from the current level
!              down to the level above JT summing the path lengths.
!
!              If the calculation is required:
               IF (JT > 0 .AND. J >= JT) THEN
                  DO I = JT + 1, J - 1
!
!                    Reset the O2 and O3 absorption cross sections.
                     CALL CSO2O3(AO2SR,AO3,TEMPC(I),TSPO2(I,JC),        &
                                 JW,JPWAV)
!
                     DEPTHA(J,JC) = DEPTHA(J,JC) + TABLENOA(I,J,JCI)*   &
                      (DO3C(I)*AO3(JW) + DO2C(I)*(AO2(JW)+AO2SR(JW)))
                     DEPTHS(J,JC) = DEPTHS(J,JC) + TABLENOA(I,J,JCI)*   &
                       DRSC(I)*SCS(JW)
                  END DO
!
!                 Reset the O2 and O3 absorption cross sections.
                  CALL CSO2O3(AO2SR,AO3,TEMPC(JT),TSPO2(JT,JC),         &
                              JW,JPWAV)
!
                  DEPTHA(J,JC) = DEPTHA(J,JC) + TABLENOA(JT,J,JCI)*     &
                    (DO3C(JT)*AO3(JW) + DO2C(JT)*(AO2(JW)+AO2SR(JW)))
                  DEPTHS(J,JC) = DEPTHS(J,JC) + TABLENOA(JT,J,JCI)*     &
                     DRSC(JT)*SCS(JW)
!
!                 Now go from JT+1 all the way to the top of the atmosphere.
                  DO I = JT + 1, JPLEV
!
!                    Reset the O2 and O3 absorption cross sections.
                     CALL CSO2O3(AO2SR,AO3,TEMPC(I),TSPO2(I,JC),        &
                                 JW,JPWAV)
!
                     DEPTHA(J,JC) = DEPTHA(J,JC) + TABLENOB(I,J,JCI)*   &
                      (DO3C(I)*AO3(JW) + DO2C(I)*(AO2(JW)+AO2SR(JW)))
                     DEPTHS(J,JC) = DEPTHS(J,JC) + TABLENOB(I,J,JCI)*   &
                       DRSC(I)*SCS(JW)
                  END DO
               END IF
            END DO
         END DO
!
!
!        Set up the vertical optical depths etc.
         DO J = 1, JPLEV
!
            TEAC(J) = DEPTHA(J,1)
            TAUC(J) = DEPTHS(J,1)
!
         END DO
         TEAC(JPLEVP1) = 0.0
         TAUC(JPLEVP1) = 0.0
!
!
!        Assign the total vertical optical depth at the ground;
!
         TAUV = 0.0
         TEAV = 0.0
!
         DO J = 1, JPLEV
!
!           Reset O2 and O3 absorption cross sections.
            CALL CSO2O3(AO2SR,AO3,TEMPC(J),TSPO2(J,1),                  &
                        JW,JPWAV)
!
!           Scattering.
            TAUV = TAUV + DALT(J)*1.0E5* DRSC(J)*SCS(JW)
!
!           Absorption.
            TEAV = TEAV + DALT(J)*1.0E5*(DO3C(J)*AO3(JW) +              &
                                DO2C(J)*(AO2(JW)+AO2SR(JW)))
!
         END DO
!
!
         TAU(1) = TAUV
         TEA(1) = TEAV
         TAU(JPLEVP1) = 0.0
         TEA(JPLEVP1) = 0.0
!
         DO J = JPLEV, 2, -1
!           Reset O2 and O3 absorption cross sections.
            CALL CSO2O3(AO2SR,AO3,TEMPC(J),TSPO2(J,1),                  &
                        JW,JPWAV)
!
            TEA(J) = TEA(J+1) + DALT(J)*1.0E5*(DO3C(J)*AO3(JW) +        &
                                DO2C(J)*(AO2(JW)+AO2SR(JW)))
            TAU(J) = TAU(J+1) + DALT(J)*1.0E5* DRSC(J)*SCS(JW)
         END DO
!
!        Calculate d tau / d delta tau
!
         DO I = 1, JPLEV
!
!           Reset O2 and O3 absorption cross sections.
            CALL CSO2O3(AO2SR,AO3,TEMPC(I),TSPO2(I,1),                  &
                        JW,JPWAV)
!
            JAC(I) = 1.0/(1.0 + (DO3C(I)*AO3(JW) + DO2C(I)*             &
                          (AO2(JW)+AO2SR(JW)))/(DRSC(I)*SCS(JW)))
!
         END DO
!
!
!        Calculate the initial scattering rate.
!
!        Zenith angle loop.
         DO JC = 1, JPCHI
!
!          Ground reflected component.
           IF (JC  <=  JPCHIN) THEN
             ARG = DEPTHA(1,JC) + DEPTHS(1,JC)
             TWOAMU = 2.0*ALBEDO*COS(TABANG(JC))*EXP(-ARG)
           ELSE
             TWOAMU = 0.0
           END IF
!
!          Level loop.
           DO JK = 1, JPLEV
!
             ARG = DEPTHA(JK,JC) + DEPTHS(JK,JC)
             TABS0(JK,JC,JW) = EXP(-ARG)
!
             IF (ALBEDO > 0.0) THEN
               TEAK = TEAC(JK)
               TAUK = TAUC(JK)
               TABS0(JK,JC,JW) = TABS0(JK,JC,JW) + TWOAMU*EI2(          &
                   MAX(0.0,TAUV - TAUK + TEAV - TEAK))
             END IF
!
!          End of level loop.
           END DO
!
!        End of zenith angle loop.
         END DO
!
!
         IF (LSCAT) THEN
!
!           Set up identity matrix
            DO JN = 1, JPLEV
              DO JK = 1, JPLEV
                DELTA(JN,JK)=0.0
              END DO
            END DO
            DO JK = 1, JPLEV
                DELTA(JK,JK)=1.0
            END DO
!
!           Two nested level loops.
            DO JN = 1, JPLEV
!
              TEAN = TEA(JN)
              TAUN = TAU(JN)
              TEANP = TEA(JN+1)
              TAUNP = TAU(JN+1)
!
              DO JK = 1, JPLEV
!
                TEAK = TEAC(JK)
                TAUK = TAUC(JK)
!
                IF(JN.NE.JK) THEN
                  A(JN,JK) = 0.5*ABS(                                   &
                  EI2(ABS(TAUK - TAUN)  + ABS(TEAK - TEAN))             &
                - EI2(ABS(TAUK - TAUNP) + ABS(TEAK - TEANP)))*JAC(JN)
                ELSE
                  A(JN,JK) = 0.5*(2.0 -                                 &
                  EI2(ABS(TAUK - TAUN)  + ABS(TEAK - TEAN))             &
                - EI2(ABS(TAUK - TAUNP) + ABS(TEAK - TEANP)))*JAC(JN)
                END IF
!
                IF (ALBEDO > 0.0) A(JN,JK) = A(JN,JK) + ALBEDO*         &
                     EI2(MAX(0.0,TAUV - TAUK  + TEAV - TEAK))*          &
                 ABS(EI3(MAX(0.0,TAUV - TAUN  + TEAV - TEAN))           &
                   - EI3(MAX(0.0,TAUV - TAUNP + TEAV - TEANP)))*        &
                 JAC(JN)
!
!               Set up B matrix
                B(JN,JK) = DELTA(JN,JK) - A(JN,JK)
!
!             End of nested level loops.
              END DO
            END DO
!
!
!           Invert the matrix B.
            CALL INVERT(B,BINV,BCOPY,INDX,JPLEV)
!
!
!           For greater than 90 degrees make sure where there is no direct
!           sun that the initial scattering rate is zero.
            DO JC = JPCHIN + 1, JPCHI
!
!              A zenith angle index.
               JCI = JC - JPCHIN
!
               DO J = 1, JPLEV
!
!                 Recall what JT was.
                  JT = IJT(JCI,J)
                  IF (JT  <=  0) TABS0(J,JC,JW) = 0.0
               END DO
            END DO
!
!
!           Zenith angle loop.
            DO JC = 1, JPCHI
!
               DO JN = 1, JPLEV
!
!                 Initialise S.
                  TABS(JN,JC,JW) = 0.0
!
!                 Sum: S(JN)=S0(JK)*BINV(JK,JN)
                  DO JK = 1, JPLEV
                     TABS(JN,JC,JW) = TABS(JN,JC,JW) + TABS0(JK,JC,JW)* &
                        BINV(JK,JN)
                  END DO
               END DO
!
!           End of zenith angle loop.
            END DO
!
!
!           For greater than 90 degrees make sure where there is no direct
!           sun that the scattering rate is zero.
            DO JC = JPCHIN + 1, JPCHI
!
!              A zenith angle index.
               JCI = JC - JPCHIN
!
               DO J = 1, JPLEV
!
!                 Recall what JT was.
                  JT = IJT(JCI,J)
                  IF (JT  <=  0) TABS(J,JC,JW) = 0.0
               END DO
            END DO
!
!
!           Check values are not negative.
            DO JC = 1, JPCHI
               DO J = 1, JPLEV
                  TABS(J,JC,JW) = MAX(0.0,TABS(J,JC,JW))
               END DO
            END DO
!
!
!        Else if no scattering.
         ELSE
!
!
!           Zenith angle loop.
            DO JC = 1, JPCHI
!
               DO JN = 1, JPLEV
!
!                 Just initial scattering rate.
                  TABS(JN,JC,JW) = TABS0(JN,JC,JW)
!
               END DO
!
            END DO
!
         END IF
!
!
!     End of wavelength loop.
      END DO
!
      END SUBROUTINE SETTAB

      SUBROUTINE INVERT(A,AINV,ACOPY,INDX,NP)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Invert photolysis matrix
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN)    :: NP
      INTEGER, INTENT(INOUT) :: INDX(NP)
      REAL, INTENT(IN)       :: A(NP,NP)
      REAL, INTENT(OUT)      :: AINV(NP,NP)
      REAL, INTENT(OUT)      :: ACOPY(NP,NP)

! Local variables
      REAL :: D
      INTEGER :: I
      INTEGER :: J
!
!     Set up the identity matrix.
      DO I = 1 , NP
         DO J = 1 , NP
            AINV(I,J) = 0.0
         END DO
         AINV(I,I) = 1.0
      END DO
!
!     Take copy
      DO I = 1 , NP
         DO J = 1 , NP
            ACOPY(I,J) = A(I,J)
         END DO
      END DO
!
!     Decompose the matrix once.
      CALL LUDCMP(ACOPY,NP,NP,INDX,D)
!
!     Find the inverse by columns.
      DO J = 1 , NP
         CALL LUBKSB(ACOPY,NP,NP,INDX,AINV(1,J))
      END DO
!
      END SUBROUTINE INVERT

      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     For matrix inversion
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: NP
      INTEGER, INTENT(IN) :: INDX(N)
      REAL, INTENT(IN)    :: A(NP,NP)
      REAL, INTENT(INOUT) :: B(N)

! Local variables
      INTEGER :: I
      INTEGER :: II
      INTEGER :: J
      INTEGER :: LL

      REAL :: SUM
!
      II=0
      DO I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II /= 0)THEN
          DO J=II,I-1
            SUM=SUM-A(I,J)*B(J)
          END DO
        ELSE IF (SUM /= 0.) THEN
          II=I
        END IF
        B(I)=SUM
      END DO
      DO I=N,1,-1
        SUM=B(I)
        IF(I < N)THEN
          DO J=I+1,N
            SUM=SUM-A(I,J)*B(J)
          END DO
        END IF
        B(I)=SUM/A(I,I)
      END DO
      END SUBROUTINE LUBKSB

      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     For matrix inversion
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: np
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(INOUT) :: d
      INTEGER, INTENT(OUT):: indx(n)
      REAL, INTENT(INOUT) :: a(np,np)

! Local variables
      INTEGER, PARAMETER :: NMAX=100
      REAL, PARAMETER :: TINY=1.0E-20

      INTEGER :: I
      INTEGER :: IMAX
      INTEGER :: J
      INTEGER :: K

      REAL :: AAMAX
      REAL :: DUM
      REAL :: SUM
      REAL :: VV(NMAX)
!
      D=1.
      DO I=1,N
        AAMAX=0.
        DO J=1,N
          IF (ABS(A(I,J)) > AAMAX) AAMAX=ABS(A(I,J))
        END DO
        IF (AAMAX == 0.) PAUSE 'SINGULAR MATRIX.'
        VV(I)=1./AAMAX
      END DO
      DO J=1,N
        IF (J > 1) THEN
          DO I=1,J-1
            SUM=A(I,J)
            IF (I > 1)THEN
              DO K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
              END DO
              A(I,J)=SUM
            END IF
          END DO
        END IF
        AAMAX=0.
        DO I=J,N
          SUM=A(I,J)
          IF (J > 1)THEN
            DO K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
            END DO
            A(I,J)=SUM
          END IF
          DUM=VV(I)*ABS(SUM)
          IF (DUM >= AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          END IF
        END DO
        IF (J /= IMAX)THEN
          DO K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
          END DO
          D=-D
          VV(IMAX)=VV(J)
        END IF
        INDX(J)=IMAX
        IF(J /= N)THEN
          IF(A(J,J) == 0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO I=J+1,N
            A(I,J)=A(I,J)*DUM
          END DO
        END IF
      END DO
      IF(A(N,N) == 0.)A(N,N)=TINY
      END SUBROUTINE LUDCMP

      SUBROUTINE CSO2O3(AO2SR,AO3,TEMPC,TSPO2,JW,JPWAV)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Calculate current O2 and O3 cross sections
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: JPWAV
      INTEGER, INTENT(IN) :: JW
      REAL, INTENT(IN)    :: TEMPC
      REAL, INTENT(IN)    :: TSPO2
      REAL, INTENT(INOUT) :: AO2SR(JPWAV)
      REAL, INTENT(INOUT) :: AO3(JPWAV)
!
!     Calculate the T dependent O3 cross-section.
      IF ( JW >= 84 .AND. JW <= 102 ) THEN
        CALL ACSO3W(JW,TEMPC,JPWAV,AO3)
      END IF
!
!     Calculate O2 cross section using the Frederick parameterisation.
!     Note - this depends on the slant path column.
      IF ( JW >= 46 .AND. JW <= 62 ) THEN
        CALL ACSSRW(TSPO2,JW,JPWAV,AO2SR)
      END IF
!
      END SUBROUTINE CSO2O3

      SUBROUTINE SETZEN(TABANG,ZENMAX,ALTC,TABT,TABO3)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name     Type               Description.
!
!     TABANG    Array of real     Solar zenith angle in radians.
!
!     ZENMAX    Array of real     Mazimum zenith angle at each altitude
!                                 at which direct sunlight is
!                                 received.
!
!     ALTC     Array of real      Altitude of each level centre in km.
!
!     TABT     Array of real      Temperature
!
!     TABO3    Array of real      O3 factor
!
!     This routine calculates the zenith angle grid used, also max zenith
!     angle and sets up temperature grid and O3 factor grid.
!
!  v1.1 Original code, adapted to UM   Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
#include "parpho.h"
!
! Subroutine interface
      REAL, INTENT(IN)  :: ALTC(JPLEV)

      REAL, INTENT(OUT) :: TABANG(JPCHI)
      REAL, INTENT(OUT) :: TABT(JPTEM)
      REAL, INTENT(OUT) :: TABO3(JPO3P)
      REAL, INTENT(OUT) :: ZENMAX(JPLEV)

! Local variables
      REAL :: COSINT
      REAL :: RADINT
      REAL :: TDIF
      REAL :: O3DIF
!
      INTEGER :: J
      INTEGER :: JC
      INTEGER :: JT
!
!
!     Spacing of angles in lookup table
      COSINT=1.0/REAL(JPCHI-1-JPS90)
      RADINT=(SZAMAX - 90.0)*PI/(180.0*REAL(JPS90))
!
!     First angle = 0 degrees
      TABANG(1) = 0.0
!
!     From 0 to 90 degrees
      DO J=2,JPCHIN-1
        JC = ABS(J-JPCHIN)
        TABANG(J) = ACOS(COSINT*REAL(JC))
      END DO
!
!     Angle 90 degrees
      TABANG(JPCHIN) = 0.5*PI
!
      DO J = JPCHIN+1,JPCHI
        TABANG(J) = 0.5*PI + REAL(J-JPCHIN)*RADINT
      END DO
!
!     Maximum zenith angle
      DO J=1,JPLEV
        ZENMAX(J) = (PI - ASIN(RE/(RE + ALTC(J))))
      END DO
!
!     Temperature grid
      IF (JPTEM == 1) THEN
        TABT(    1)=0.5*(TMIN + TMAX)
      ELSE
        TABT(    1)=TMIN
        TABT(JPTEM)=TMAX
        TDIF=(TMAX-TMIN)/REAL(JPTEM-1)
        DO JT=2,JPTEM-1
          TABT(JT) = TMIN + TDIF*REAL(JT-1)
        END DO
      END IF
!
!     O3 grid
      IF (JPO3P == 1) THEN
        TABO3(   1)=1.0
      ELSE
        TABO3(    1)=O3MIN
        TABO3(JPO3P)=O3MAX
        O3DIF=(O3MAX-O3MIN)/REAL(JPO3P-1)
        DO JT=2,JPO3P-1
          TABO3(JT) = O3MIN + O3DIF*REAL(JT-1)
        END DO
      END IF
!
      END SUBROUTINE SETZEN

      SUBROUTINE fill_spectra
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                     C
!      INITIALISE T INDEPENDENT RATE CONSTANTS AND CROSS SECTIONS     C
!                    DATA FROM JPL 94-26                              C
!                                                                     C
!  v1.1 2/12/2003 Original code, adapted to UM.
!                 Remove BLOCK DATA structures.   Olaf Morgenstern
!  v1.2 9/3/2006  Adapted to include branching for O2 photolysis
!                                                 Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT NONE

#include "parpho.h"
#include "crossec.h"

! Local variables
      INTEGER :: i
!
      WAVENM(1:108) = (/                                                &
        1.215E+02,1.167E+02,1.173E+02,1.179E+02,1.187E+02,1.194E+02,    &
        1.202E+02,1.208E+02,1.216E+02,1.223E+02,1.231E+02,1.238E+02,    &
        1.246E+02,1.254E+02,1.262E+02,1.270E+02,1.278E+02,1.286E+02,    &
        1.294E+02,1.303E+02,1.311E+02,1.320E+02,1.329E+02,1.338E+02,    &
        1.346E+02,1.356E+02,1.365E+02,1.374E+02,1.384E+02,1.399E+02,    &
        1.418E+02,1.439E+02,1.459E+02,1.481E+02,1.504E+02,1.526E+02,    &
        1.550E+02,1.574E+02,1.600E+02,1.626E+02,1.653E+02,1.681E+02,    &
        1.709E+02,1.731E+02,1.746E+02,1.762E+02,1.778E+02,1.794E+02,    &
        1.810E+02,1.827E+02,1.843E+02,1.861E+02,1.878E+02,1.896E+02,    &
        1.914E+02,1.932E+02,1.951E+02,1.970E+02,1.990E+02,2.010E+02,    &
        2.031E+02,2.051E+02,2.073E+02,2.094E+02,2.116E+02,2.139E+02,    &
        2.162E+02,2.186E+02,2.210E+02,2.235E+02,2.260E+02,2.286E+02,    &
        2.312E+02,2.339E+02,2.367E+02,2.395E+02,2.424E+02,2.454E+02,    &
        2.485E+02,2.516E+02,2.548E+02,2.581E+02,2.614E+02,2.649E+02,    &
        2.685E+02,2.721E+02,2.759E+02,2.797E+02,2.837E+02,2.878E+02,    &
        2.920E+02,2.963E+02,3.008E+02,3.053E+02,3.100E+02,3.150E+02,    &
        3.200E+02,3.250E+02,3.300E+02,3.350E+02,3.400E+02,3.450E+02,    &
        3.500E+02,3.550E+02,3.600E+02,3.650E+02,3.700E+02,3.750E+02/)
      WAVENM(109:203) = (/                                              &
        3.800E+02,3.850E+02,3.900E+02,3.950E+02,4.000E+02,4.050E+02,    &
        4.100E+02,4.150E+02,4.200E+02,4.250E+02,4.300E+02,4.350E+02,    &
        4.400E+02,4.450E+02,4.500E+02,4.550E+02,4.600E+02,4.650E+02,    &
        4.700E+02,4.750E+02,4.800E+02,4.850E+02,4.900E+02,4.950E+02,    &
        5.000E+02,5.050E+02,5.100E+02,5.150E+02,5.200E+02,5.250E+02,    &
        5.300E+02,5.350E+02,5.400E+02,5.450E+02,5.500E+02,5.550E+02,    &
        5.600E+02,5.650E+02,5.700E+02,5.750E+02,5.800E+02,5.850E+02,    &
        5.900E+02,5.950E+02,6.000E+02,6.050E+02,6.100E+02,6.150E+02,    &
        6.200E+02,6.250E+02,6.300E+02,6.350E+02,6.400E+02,6.450E+02,    &
        6.500E+02,6.550E+02,6.600E+02,6.650E+02,6.700E+02,6.750E+02,    &
        6.800E+02,6.850E+02,6.900E+02,6.950E+02,7.000E+02,7.050E+02,    &
        7.100E+02,7.150E+02,7.200E+02,7.250E+02,7.300E+02,7.350E+02,    &
        7.400E+02,7.450E+02,7.500E+02,7.550E+02,7.600E+02,7.650E+02,    &
        7.700E+02,7.750E+02,7.800E+02,7.850E+02,7.900E+02,7.950E+02,    &
        8.000E+02,8.050E+02,8.100E+02,8.150E+02,8.200E+02,8.250E+02,    &
        8.300E+02,8.350E+02,8.400E+02,8.450E+02,8.500E+02/)
      WAVECM(1:108) = (/                                                &
        8.230E+04,8.573E+04,8.525E+04,8.478E+04,8.428E+04,8.375E+04,    &
        8.323E+04,8.275E+04,8.224E+04,8.173E+04,8.123E+04,8.074E+04,    &
        8.026E+04,7.974E+04,7.924E+04,7.874E+04,7.825E+04,7.776E+04,    &
        7.725E+04,7.675E+04,7.625E+04,7.576E+04,7.527E+04,7.477E+04,    &
        7.427E+04,7.377E+04,7.326E+04,7.275E+04,7.225E+04,7.151E+04,    &
        7.052E+04,6.952E+04,6.852E+04,6.752E+04,6.651E+04,6.551E+04,    &
        6.452E+04,6.351E+04,6.250E+04,6.150E+04,6.050E+04,5.949E+04,    &
        5.850E+04,5.775E+04,5.726E+04,5.675E+04,5.624E+04,5.574E+04,    &
        5.525E+04,5.473E+04,5.426E+04,5.373E+04,5.325E+04,5.274E+04,    &
        5.225E+04,5.176E+04,5.126E+04,5.076E+04,5.025E+04,4.975E+04,    &
        4.924E+04,4.876E+04,4.824E+04,4.776E+04,4.726E+04,4.675E+04,    &
        4.625E+04,4.575E+04,4.525E+04,4.474E+04,4.425E+04,4.374E+04,    &
        4.325E+04,4.275E+04,4.225E+04,4.175E+04,4.125E+04,4.075E+04,    &
        4.024E+04,3.975E+04,3.925E+04,3.874E+04,3.826E+04,3.775E+04,    &
        3.724E+04,3.675E+04,3.624E+04,3.575E+04,3.525E+04,3.475E+04,    &
        3.425E+04,3.375E+04,3.324E+04,3.275E+04,3.226E+04,3.175E+04,    &
        3.125E+04,3.077E+04,3.030E+04,2.985E+04,2.941E+04,2.899E+04,    &
        2.857E+04,2.817E+04,2.778E+04,2.740E+04,2.703E+04,2.667E+04/)
      WAVECM(109:203) = (/                                              &
        2.632E+04,2.597E+04,2.564E+04,2.532E+04,2.500E+04,2.469E+04,    &
        2.439E+04,2.410E+04,2.381E+04,2.353E+04,2.326E+04,2.299E+04,    &
        2.273E+04,2.247E+04,2.222E+04,2.198E+04,2.174E+04,2.151E+04,    &
        2.128E+04,2.105E+04,2.083E+04,2.062E+04,2.041E+04,2.020E+04,    &
        2.000E+04,1.980E+04,1.961E+04,1.942E+04,1.923E+04,1.905E+04,    &
        1.887E+04,1.869E+04,1.852E+04,1.835E+04,1.818E+04,1.802E+04,    &
        1.786E+04,1.770E+04,1.754E+04,1.739E+04,1.724E+04,1.709E+04,    &
        1.695E+04,1.681E+04,1.667E+04,1.653E+04,1.639E+04,1.626E+04,    &
        1.613E+04,1.600E+04,1.587E+04,1.575E+04,1.562E+04,1.550E+04,    &
        1.538E+04,1.527E+04,1.515E+04,1.504E+04,1.493E+04,1.481E+04,    &
        1.471E+04,1.460E+04,1.449E+04,1.439E+04,1.429E+04,1.418E+04,    &
        1.408E+04,1.399E+04,1.389E+04,1.379E+04,1.370E+04,1.361E+04,    &
        1.351E+04,1.342E+04,1.333E+04,1.324E+04,1.316E+04,1.307E+04,    &
        1.299E+04,1.290E+04,1.282E+04,1.274E+04,1.266E+04,1.258E+04,    &
        1.250E+04,1.242E+04,1.235E+04,1.227E+04,1.220E+04,1.212E+04,    &
        1.205E+04,1.198E+04,1.190E+04,1.183E+04,1.176E+04/)
!     Quanta: intervals 1-45 Ackerman. 46-203 WMO 1985.
      QUANTA(1:108) = (/                                                &
        4.000E+11,1.030E+08,2.660E+08,1.120E+08,1.240E+08,1.820E+08,    &
        1.900E+08,7.400E+08,2.280E+09,3.670E+09,1.360E+09,1.610E+09,    &
        1.320E+09,1.410E+09,3.110E+09,1.060E+09,1.370E+09,1.020E+09,    &
        1.140E+09,7.290E+09,2.200E+09,1.590E+09,2.210E+09,1.240E+10,    &
        1.990E+09,3.090E+09,2.570E+09,2.740E+09,3.100E+09,7.600E+09,    &
        1.010E+10,1.300E+10,1.820E+10,2.330E+10,2.660E+10,2.900E+10,    &
        3.600E+10,4.750E+10,6.400E+10,5.490E+10,1.190E+11,1.760E+11,    &
        2.320E+11,1.440E+11,1.830E+11,1.740E+11,2.100E+11,2.380E+11,    &
        3.040E+11,3.190E+11,2.930E+11,3.620E+11,4.730E+11,5.610E+11,    &
        6.630E+11,6.900E+11,9.560E+11,1.150E+12,1.270E+12,1.520E+12,    &
        1.780E+12,2.200E+12,2.690E+12,4.540E+12,7.140E+12,8.350E+12,    &
        8.390E+12,1.080E+13,1.180E+13,1.600E+13,1.340E+13,1.410E+13,    &
        1.570E+13,1.380E+13,1.600E+13,1.450E+13,2.200E+13,1.990E+13,    &
        1.970E+13,1.940E+13,2.910E+13,4.950E+13,4.530E+13,1.070E+14,    &
        1.200E+14,1.100E+14,1.040E+14,8.240E+13,1.520E+14,2.150E+14,    &
        3.480E+14,3.400E+14,3.220E+14,4.230E+14,4.950E+14,5.440E+14,    &
        5.930E+14,6.950E+14,8.150E+14,7.810E+14,8.350E+14,8.140E+14,    &
        8.530E+14,9.170E+14,8.380E+14,1.040E+15,1.100E+15,9.790E+14/)
      QUANTA(109:203) = (/                                              &
        1.130E+15,8.890E+14,1.140E+15,9.170E+14,1.690E+15,1.700E+15,    &
        1.840E+15,1.870E+15,1.950E+15,1.810E+15,1.670E+15,1.980E+15,    &
        2.020E+15,2.180E+15,2.360E+15,2.310E+15,2.390E+15,2.380E+15,    &
        2.390E+15,2.440E+15,2.510E+15,2.300E+15,2.390E+15,2.480E+15,    &
        2.400E+15,2.460E+15,2.490E+15,2.320E+15,2.390E+15,2.420E+15,    &
        2.550E+15,2.510E+15,2.490E+15,2.550E+15,2.530E+15,2.540E+15,    &
        2.500E+15,2.570E+15,2.580E+15,2.670E+15,2.670E+15,2.700E+15,    &
        2.620E+15,2.690E+15,2.630E+15,2.680E+15,2.660E+15,2.590E+15,    &
        2.690E+15,2.610E+15,2.620E+15,2.620E+15,2.630E+15,2.600E+15,    &
        2.550E+15,2.480E+15,2.570E+15,2.610E+15,2.610E+15,2.620E+15,    &
        2.620E+15,2.570E+15,2.520E+15,2.600E+15,2.580E+15,2.520E+15,    &
        2.510E+15,2.480E+15,2.450E+15,2.480E+15,2.450E+15,2.440E+15,    &
        2.390E+15,2.400E+15,2.410E+15,2.400E+15,2.380E+15,2.340E+15,    &
        2.320E+15,2.300E+15,2.330E+15,2.340E+15,2.290E+15,2.290E+15,    &
        2.270E+15,2.270E+15,2.200E+15,2.220E+15,2.180E+15,2.200E+15,    &
        2.140E+15,2.140E+15,2.130E+15,2.090E+15,2.050E+15/)
!     QENO2: JPL 1992.
      QENO2 = (/                                                        &
        (1.0,i=1,84),                                                   &
        1.000E+00,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.990E-01,    &
        9.990E-01,9.980E-01,9.970E-01,9.960E-01,9.950E-01,9.940E-01,    &
        9.930E-01,9.920E-01,9.910E-01,9.900E-01,9.890E-01,9.880E-01,    &
        9.870E-01,9.860E-01,9.840E-01,9.830E-01,9.810E-01,9.790E-01,    &
        9.740E-01,9.690E-01,9.600E-01,9.270E-01,6.940E-01,3.550E-01,    &
        1.340E-01,6.000E-02,1.800E-02,1.000E-03,0.000E+00,0.000E+00,    &
        (0.0,i=1,83)/)
!     NO2: intervals 61-203 JPL 1992 T=298K
      ANO2 = (/                                                         &
        (0.0,i=1,48),                                                   &
        0.000E+00,0.000E+00,0.000E+00,2.670E-19,2.780E-19,2.900E-19,    &
        2.790E-19,2.600E-19,2.420E-19,2.450E-19,2.480E-19,2.750E-19,    &
        4.145E-19,4.478E-19,4.454E-19,4.641E-19,4.866E-19,4.818E-19,    &
        5.022E-19,4.441E-19,4.713E-19,3.772E-19,3.929E-19,2.740E-19,    &
        2.778E-19,1.689E-19,1.618E-19,8.812E-20,7.472E-20,3.909E-20,    &
        2.753E-20,2.007E-20,1.973E-20,2.111E-20,2.357E-20,2.698E-20,    &
        3.247E-20,3.785E-20,5.030E-20,5.880E-20,7.000E-20,8.150E-20,    &
        9.720E-20,1.154E-19,1.344E-19,1.589E-19,1.867E-19,2.153E-19,    &
        2.477E-19,2.807E-19,3.133E-19,3.425E-19,3.798E-19,4.065E-19,    &
        4.313E-19,4.717E-19,4.833E-19,5.166E-19,5.315E-19,5.508E-19,    &
        5.644E-19,5.757E-19,5.927E-19,5.845E-19,6.021E-19,5.781E-19,    &
        5.999E-19,5.651E-19,5.812E-19,0.000E+00,0.000E+00,0.000E+00,    &
        (0.0,i=1,83)/)
!
      ANO31 = (/                                                        &
        (0.0,i=1,144),                                                  &
        0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,9.700E-19,    &
        2.200E-18,1.620E-18,1.160E-18,9.000E-19,4.600E-19,2.400E-19,    &
        3.800E-19,1.000E-19,0.000E+00,0.000E+00,0.000E+00,0.000E+00,    &
        (0.0,i=1,41)/)
!
      ANO32 = (/                                                        &
        (0.0,i=1,126),                                                  &
        5.900E-19,7.200E-19,7.000E-19,8.600E-19,1.000E-18,1.150E-18,    &
        1.060E-18,1.240E-18,1.610E-18,1.450E-18,1.740E-18,1.810E-18,    &
        2.070E-18,2.420E-18,1.830E-18,2.560E-18,2.650E-18,3.710E-18,    &
        3.430E-18,3.270E-18,3.330E-18,4.050E-18,4.130E-18,4.830E-18,    &
        4.980E-18,3.330E-18,2.300E-18,1.740E-18,9.900E-19,7.900E-19,    &
        1.950E-18,9.100E-19,2.000E-19,0.000E+00,0.000E+00,0.000E+00,    &
        (0.0,i=1,41)/)
!
      ANO = 0.0
!
      AN2O =(/                                                          &
        (0.0,i=1,42),                                                   &
        0.000E+00,0.000E+00,1.250E-19,1.330E-19,1.400E-19,1.450E-19,    &
        1.460E-19,1.450E-19,1.430E-19,1.380E-19,1.290E-19,1.180E-19,    &
        1.000E-19,9.000E-20,7.500E-20,5.900E-20,4.600E-20,3.500E-20,    &
        2.600E-20,1.950E-20,1.350E-20,9.000E-21,5.600E-21,3.500E-21,    &
        2.100E-21,1.250E-21,7.500E-22,4.400E-22,2.300E-22,1.000E-22,    &
        7.000E-23,3.650E-23,2.150E-23,1.150E-23,0.000E+00,0.000E+00,    &
        (0.0,i=1,125)/)
!
      AN2O5 = (/                                                        &
        (0.0,i=1,54),                                                   &
        0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,9.100E-18,    &
        8.800E-18,8.100E-18,6.950E-18,5.900E-18,5.050E-18,4.200E-18,    &
        3.450E-18,2.650E-18,2.100E-18,1.650E-18,1.350E-18,1.100E-18,    &
        9.000E-19,8.000E-19,7.150E-19,6.350E-19,5.650E-19,5.100E-19,    &
        4.300E-19,3.750E-19,3.250E-19,2.800E-19,2.400E-19,2.000E-19,    &
        1.700E-19,1.450E-19,1.250E-19,1.150E-19,1.000E-19,0.000E+00,    &
        (0.0,i=1,113)/)
!
      ACH4 = 0.0
!
      AHNO2 = (/                                                        &
        (0.0,i=1,90),                                                   &
        0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,4.200E-21,    &
        2.100E-20,4.040E-20,7.290E-20,6.450E-20,1.050E-19,8.540E-20,    &
        6.830E-20,2.460E-19,6.830E-20,1.330E-19,1.190E-19,2.700E-20,    &
        7.780E-20,1.900E-19,1.200E-20,0.000E+00,0.000E+00,0.000E+00,    &
        (0.0,i=1,89)/)
!      HNO3 at 298 K
       AHNO3 = (/(0.0,i=1,48),                                          &
        0.000E+00,8.305E-18,1.336E-17,1.575E-17,1.491E-17,1.385E-17,    &
        1.265E-17,1.150E-17,1.012E-17,8.565E-18,6.739E-18,5.147E-18,    &
        3.788E-18,2.719E-18,1.796E-18,1.180E-18,7.377E-19,4.487E-19,    &
        2.810E-19,1.826E-19,1.324E-19,1.010E-19,8.020E-20,6.479E-20,    &
        5.204E-20,4.178E-20,3.200E-20,2.657E-20,2.298E-20,2.086E-20,    &
        1.991E-20,1.962E-20,1.952E-20,1.929E-20,1.882E-20,1.804E-20,    &
        1.681E-20,1.526E-20,1.335E-20,1.136E-20,9.242E-21,7.186E-21,    &
        5.320E-21,3.705E-21,2.393E-21,1.442E-21,8.140E-22,4.131E-22,    &
        1.970E-22,9.434E-23,4.310E-23,2.204E-23,1.030E-23,5.841E-24,    &
        4.170E-24,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,    &
        (0.0,i=1,95)/)
!     HO2NO2: JPL 1992
      APNA = (/                                                         &
        (0.0,i=1,48),                                                   &
        0.000E+00,0.000E+00,0.000E+00,3.255E-18,7.139E-18,9.735E-18,    &
        9.616E-18,8.926E-18,8.104E-18,7.080E-18,6.090E-18,5.184E-18,    &
        4.335E-18,3.639E-18,3.019E-18,2.517E-18,2.102E-18,1.752E-18,    &
        1.490E-18,1.282E-18,1.123E-18,9.978E-19,9.038E-19,8.287E-19,    &
        7.648E-19,7.054E-19,6.467E-19,5.907E-19,5.350E-19,4.823E-19,    &
        4.335E-19,3.925E-19,3.525E-19,3.085E-19,2.690E-19,2.310E-19,    &
        1.955E-19,1.605E-19,1.259E-19,9.518E-20,6.929E-20,4.813E-20,    &
        3.240E-20,2.092E-20,1.313E-20,8.704E-21,5.000E-21,3.000E-21,    &
        2.000E-21,1.000E-21,0.000E+00,0.000E+00,0.000E+00,0.000E+00,    &
        (0.0,i=1,101)/)
!     ClONO2: JPL 1994
      ACNITA = (/                                                       &
         (0.0,i=1,54),                                                  &
         4.32E-19, 1.98E-18, 2.91E-18, 3.01E-18, 2.87E-18, 2.79E-18,    &
         2.78E-18, 2.84E-18, 2.96E-18, 3.10E-18, 3.26E-18, 3.39E-18,    &
         3.45E-18, 3.39E-18, 3.24E-18, 2.97E-18, 2.64E-18, 2.27E-18,    &
         1.92E-18, 1.59E-18, 1.31E-18, 1.09E-18, 8.97E-19, 7.44E-19,    &
         6.07E-19, 5.13E-19, 4.35E-19, 3.70E-19, 3.15E-19, 2.66E-19,    &
         2.21E-19, 1.84E-19, 1.50E-19, 1.21E-19, 9.52E-20, 7.33E-20,    &
         5.50E-20, 4.01E-20, 2.97E-20, 2.19E-20, 1.60E-20, 1.14E-20,    &
         8.31E-21, 6.11E-21, 4.66E-21, 3.66E-21, 3.02E-21, 2.58E-21,    &
         2.29E-21, 2.08E-21, 2.00E-21, 1.79E-21, 1.59E-21, 1.41E-21,    &
         1.21E-21, 1.06E-21, 9.09E-22, 7.59E-22, 6.38E-22, 5.38E-22,    &
         4.44E-22, 3.67E-22, 3.16E-22, 2.31E-22, 1.89E-22, 5.26E-23,    &
         (0.0,i=1,83)/)
!     O2 Herzberg continuum. WMO 1985.
      AO2 = (/                                                          &
        1.000E+00,2.000E-20,1.250E-18,2.550E-19,3.000E-20,3.750E-19,    &
        4.450E-18,8.350E-18,6.000E-19,2.350E-19,4.500E-19,3.350E-19,    &
        1.750E-17,8.950E-19,4.300E-19,1.100E-19,2.050E-19,4.430E-19,    &
        5.550E-19,4.200E-19,6.850E-19,1.450E-18,2.250E-18,2.300E-18,    &
        4.550E-18,7.230E-18,9.500E-18,1.230E-17,1.320E-17,1.360E-17,    &
        1.400E-17,1.480E-17,1.410E-17,1.290E-17,1.150E-17,9.910E-18,    &
        8.240E-18,6.580E-18,4.970E-18,3.450E-18,2.080E-18,1.230E-18,    &
        7.220E-19,4.580E-19,2.740E-19,4.610E-24,5.030E-24,5.460E-24,    &
        5.880E-24,6.290E-24,6.680E-24,7.040E-24,7.360E-24,7.640E-24,    &
        7.870E-24,8.040E-24,8.140E-24,8.170E-24,8.130E-24,8.010E-24,    &
        7.840E-24,7.630E-24,7.330E-24,6.990E-24,6.450E-24,5.810E-24,    &
        5.230E-24,4.710E-24,4.260E-24,3.800E-24,3.350E-24,2.900E-24,    &
        2.450E-24,2.050E-24,1.690E-24,1.300E-24,9.300E-25,0.000E+00,    &
        (0.0,i=1,125)/)
! Herzberg continuum, with cross section for lambda > 204 nm
! updated with JPL 2002 recommendation.
      AO2 = (/                                                          &
        1.000E+00,2.000E-20,1.250E-18,2.550E-19,3.000E-20,3.750E-19,    &
        4.450E-18,8.350E-18,6.000E-19,2.350E-19,4.500E-19,3.350E-19,    &
        1.750E-17,8.950E-19,4.300E-19,1.100E-19,2.050E-19,4.430E-19,    &
        5.550E-19,4.200E-19,6.850E-19,1.450E-18,2.250E-18,2.300E-18,    &
        4.550E-18,7.230E-18,9.500E-18,1.230E-17,1.320E-17,1.360E-17,    &
        1.400E-17,1.480E-17,1.410E-17,1.290E-17,1.150E-17,9.910E-18,    &
        8.240E-18,6.580E-18,4.970E-18,3.450E-18,2.080E-18,1.230E-18,    &
        7.220E-19,4.580E-19,2.740E-19,4.610E-24,5.030E-24,5.460E-24,    &
        5.880E-24,6.290E-24,6.680E-24,7.040E-24,7.360E-24,7.640E-24,    &
        7.870E-24,8.040E-24,8.140E-24,8.170E-24,8.130E-24,8.010E-24,    &
        7.840E-24,                                                      &
      7.32052e-24,7.00070e-24,6.60698e-24,6.11726e-24,5.73530e-24,      &
      5.30105e-24,                                                      &
      4.73948e-24,4.25983e-24,3.78649e-24,3.21096e-24,2.68894e-24,      &
      2.21678e-24,                                                      &
      1.79184e-24,1.38277e-24,1.05234e-24,0.0,0.0,                      &
        (0.0,i=1,125)/)
! AO2A: O2 -> 2 O(3P) channel
      AO2A(46:203)  = AO2(46:203)
      AO2A(1:45)    = 0.
! AO2B: O2 -> O(3P) + O(1D)
      AO2B(46:203)  = 0.
      AO2B(1:45)    = AO2(1:45)

!     O2 Schumann Runge bands
      AO2SR = 0.0
!     O3: intervals 46-105 JPL 1992. 106-203 WMO 1985. T=298K
      AO3(1:108) = (/                                                   &
        0.000E+00,7.800E-18,7.970E-18,8.660E-18,9.510E-18,1.250E-17,    &
        1.840E-17,2.190E-17,2.300E-17,2.260E-17,2.060E-17,1.300E-17,    &
        8.910E-18,7.240E-18,6.090E-18,5.660E-18,5.870E-18,6.470E-18,    &
        8.140E-18,1.240E-17,1.520E-17,1.470E-17,1.510E-17,1.510E-17,    &
        1.650E-17,1.540E-17,1.350E-17,1.050E-17,7.970E-18,7.170E-18,    &
        6.280E-18,5.660E-18,5.230E-18,4.470E-18,3.690E-18,2.930E-18,    &
        2.190E-18,1.630E-18,1.200E-18,9.770E-19,8.660E-19,8.140E-19,    &
        8.170E-19,8.570E-19,8.400E-19,8.110E-19,7.990E-19,7.860E-19,    &
        7.630E-19,7.290E-19,6.880E-19,6.220E-19,5.760E-19,5.260E-19,    &
        4.760E-19,4.280E-19,3.830E-19,3.470E-19,3.230E-19,3.140E-19,    &
        3.260E-19,3.640E-19,4.340E-19,5.420E-19,6.990E-19,9.210E-19,    &
        1.190E-18,1.550E-18,1.990E-18,2.560E-18,3.230E-18,4.000E-18,    &
        4.830E-18,5.790E-18,6.860E-18,7.970E-18,9.000E-18,1.000E-17,    &
        1.080E-17,1.130E-17,1.150E-17,1.120E-17,1.060E-17,9.650E-18,    &
        8.340E-18,6.920E-18,5.420E-18,4.020E-18,2.770E-18,1.790E-18,    &
        1.090E-18,6.240E-19,3.430E-19,1.850E-19,9.800E-20,5.010E-20,    &
        2.490E-20,1.200E-20,6.170E-21,2.740E-21,1.170E-21,5.880E-22,    &
        2.660E-22,1.090E-22,5.490E-23,0.000E+00,0.000E+00,0.000E+00/)
      AO3(109:203) = (/                                                 &
        0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,    &
        2.910E-23,3.140E-23,3.990E-23,6.540E-23,6.830E-23,8.660E-23,    &
        1.250E-22,1.490E-22,1.710E-22,2.120E-22,3.570E-22,3.680E-22,    &
        4.060E-22,4.890E-22,7.110E-22,8.430E-22,8.280E-22,9.090E-22,    &
        1.220E-21,1.620E-21,1.580E-21,1.600E-21,1.780E-21,2.070E-21,    &
        2.550E-21,2.740E-21,2.880E-21,3.070E-21,3.170E-21,3.360E-21,    &
        3.880E-21,4.310E-21,4.670E-21,4.750E-21,4.550E-21,4.350E-21,    &
        4.420E-21,4.610E-21,4.890E-21,4.840E-21,4.540E-21,4.240E-21,    &
        3.900E-21,3.600E-21,3.430E-21,3.170E-21,2.740E-21,2.610E-21,    &
        2.420E-21,2.200E-21,2.020E-21,1.850E-21,1.670E-21,1.540E-21,    &
        1.420E-21,1.250E-21,1.120E-21,1.020E-21,9.200E-22,8.400E-22,    &
        7.700E-22,6.900E-22,6.300E-22,5.700E-22,5.250E-22,4.750E-22,    &
        4.470E-22,4.200E-22,3.750E-22,3.250E-22,2.920E-22,2.760E-22,    &
        2.700E-22,2.800E-22,2.850E-22,2.520E-22,2.220E-22,1.820E-22,    &
        1.630E-22,1.750E-22,1.900E-22,1.850E-22,1.700E-22,1.520E-22,    &
        1.420E-22,1.400E-22,1.400E-22,1.420E-22,1.450E-22/)
!
      AH2O = (/                                                         &
        (0.0,i=1,42),                                                   &
        0.000E+00,0.000E+00,2.800E-18,2.300E-18,1.700E-18,1.020E-18,    &
        4.800E-19,2.050E-19,7.300E-20,3.000E-20,1.300E-20,5.400E-21,    &
        2.400E-21,1.000E-21,4.200E-22,1.800E-22,7.400E-23,0.000E+00,    &
        (0.0,i=1,143)/)
!     H2O2: JPL 1992 T=298K
      AH2O2 =(/                                                         &
        (0.0,i=1,48),                                                   &
        0.000E+00,0.000E+00,0.000E+00,2.148E-19,4.724E-19,6.469E-19,    &
        6.398E-19,6.007E-19,5.620E-19,5.258E-19,4.910E-19,4.603E-19,    &
        4.316E-19,4.070E-19,3.844E-19,3.631E-19,3.409E-19,3.179E-19,    &
        2.945E-19,2.709E-19,2.493E-19,2.287E-19,2.098E-19,1.915E-19,    &
        1.738E-19,1.565E-19,1.407E-19,1.264E-19,1.131E-19,1.004E-19,    &
        8.839E-20,7.766E-20,6.760E-20,5.797E-20,4.972E-20,4.220E-20,    &
        3.549E-20,2.994E-20,2.485E-20,2.033E-20,1.611E-20,1.332E-20,    &
        1.070E-20,8.380E-21,6.494E-21,5.022E-21,3.900E-21,2.900E-21,    &
        2.200E-21,1.600E-21,1.300E-21,1.000E-21,7.000E-22,5.000E-22,    &
        4.000E-22,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,    &
        (0.0,i=1,95)/)
!
      ACH2O = (/                                                        &
        (0.0,i=1,84),                                                   &
        0.000E+00,0.000E+00,0.000E+00,0.000E+00,2.400E-20,2.400E-20,    &
        3.200E-20,3.200E-20,3.300E-20,3.300E-20,3.100E-20,3.100E-20,    &
        2.400E-20,2.400E-20,2.400E-20,2.400E-20,2.000E-20,2.000E-20,    &
        8.000E-21,8.000E-21,2.000E-21,2.000E-21,0.000E+00,0.000E+00,    &
        (0.0,i=1,95)/)
!
      ACO2 = (/                                                         &
        (0.0,i=1,42),                                                   &
        0.000E+00,0.000E+00,0.000E+00,6.000E-21,6.000E-21,2.000E-21,    &
        2.000E-21,2.000E-21,5.800E-22,5.800E-22,5.800E-22,9.000E-23,    &
        9.000E-23,9.000E-23,3.200E-23,3.200E-23,3.000E-24,3.000E-24,    &
        3.000E-24,3.000E-24,3.000E-24,3.000E-24,0.000E+00,0.000E+00,    &
        (0.0,i=1,137)/)
!     CFCl3: JPL 1992
      AF11 = (/                                                         &
        (0.0,i=1,36),                                                   &
        0.000E+00,0.000E+00,0.000E+00,0.000E+00,3.139E-19,2.478E-18,    &
        3.182E-18,3.168E-18,3.141E-18,3.098E-18,3.042E-18,3.096E-18,    &
        2.968E-18,2.765E-18,2.558E-18,2.319E-18,2.107E-18,1.839E-18,    &
        1.574E-18,1.332E-18,1.092E-18,8.911E-19,7.221E-19,5.751E-19,    &
        4.389E-19,3.340E-19,2.377E-19,1.700E-19,1.171E-19,7.662E-20,    &
        5.082E-20,3.184E-20,1.970E-20,1.206E-20,8.000E-21,4.834E-21,    &
        2.831E-21,1.629E-21,9.327E-22,5.209E-22,3.013E-22,1.617E-22,    &
        9.035E-23,5.427E-23,3.474E-23,2.141E-23,9.102E-24,1.499E-25,    &
        (0.0,i=1,119)/)
!     CF2Cl2: JPL 1992
      AF12 = (/                                                         &
        (0.0,i=1,36),                                                   &
        0.000E+00,0.000E+00,0.000E+00,0.000E+00,9.716E-20,8.639E-19,    &
        1.370E-18,1.630E-18,1.752E-18,1.846E-18,1.894E-18,1.778E-18,    &
        1.655E-18,1.515E-18,1.312E-18,1.030E-18,8.615E-19,6.682E-19,    &
        4.953E-19,3.567E-19,2.494E-19,1.659E-19,1.088E-19,7.081E-20,    &
        4.327E-20,2.667E-20,1.753E-20,9.740E-21,5.336E-21,2.976E-21,    &
        2.572E-21,3.840E-21,5.644E-22,3.270E-22,1.769E-22,8.850E-23,    &
        4.328E-23,2.236E-23,1.040E-23,3.751E-24,1.146E-24,0.000E+00,    &
        (0.0,i=1,125)/)
!
      AF22 = (/                                                         &
        (0.0,i=1,42),                                                   &
        0.000E+00,0.000E+00,5.180E-20,3.890E-20,2.950E-20,2.250E-20,    &
        1.620E-20,1.130E-20,7.350E-21,4.710E-21,3.500E-21,2.520E-21,    &
        1.800E-21,1.190E-21,8.150E-22,5.220E-22,3.160E-22,2.400E-22,    &
        1.710E-22,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,    &
        (0.0,i=1,137)/)
!
      AF113 = (/                                                        &
        (0.0,i=1,42),                                                   &
        0.000E+00,0.000E+00,2.000E-18,1.860E-18,1.780E-18,1.660E-18,    &
        1.510E-18,1.350E-18,1.230E-18,1.020E-18,8.510E-19,6.760E-19,    &
        5.370E-19,4.070E-19,2.880E-19,2.090E-19,1.450E-19,1.000E-19,    &
        7.080E-20,4.790E-20,3.090E-20,2.000E-20,1.260E-20,7.760E-21,    &
        4.790E-21,3.020E-21,1.820E-21,1.150E-21,0.000E+00,0.000E+00,    &
        (0.0,i=1,131)/)
!
      ACH3CL = (/                                                       &
        (0.0,i=1,48),                                                   &
        0.000E+00,0.000E+00,0.000E+00,2.450E-19,1.810E-19,1.360E-19,    &
        9.930E-20,7.010E-20,4.800E-20,3.240E-20,2.160E-20,1.390E-20,    &
        8.540E-21,5.630E-21,3.760E-21,2.340E-21,1.430E-21,8.970E-22,    &
        5.500E-22,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,    &
        (0.0,i=1,131)/)
!
      ACCL4 = (/                                                        &
        (0.0,i=1,42),                                                   &
        0.000E+00,0.000E+00,1.000E-17,1.010E-17,9.890E-18,8.410E-18,    &
        6.750E-18,5.400E-18,4.290E-18,3.150E-18,2.270E-18,1.570E-18,    &
        1.110E-18,8.090E-19,6.980E-19,6.690E-19,6.550E-19,6.350E-19,    &
        6.130E-19,5.840E-19,5.370E-19,4.840E-19,4.090E-19,3.360E-19,    &
        2.660E-19,2.030E-19,1.490E-19,1.040E-19,7.200E-20,5.040E-20,    &
        3.390E-20,2.190E-20,1.360E-20,8.520E-21,5.390E-21,3.410E-21,    &
        2.260E-21,1.390E-21,6.910E-22,3.570E-22,2.050E-22,1.280E-22,    &
        7.670E-23,4.270E-23,2.400E-23,0.000E+00,0.000E+00,0.000E+00,    &
        (0.0,i=1,113)/)
!
      AF12B1 = (/                                                       &
               (0.0,i=1,53),                                            &
               4.20E-19,5.00E-19,6.08E-19,7.18E-19,8.24E-19,            &
               9.20E-19,9.99E-19,1.06E-18,1.09E-18,1.09E-18,1.06E-18,   &
               1.01E-18,9.36E-19,8.49E-19,7.56E-19,6.62E-19,5.72E-19,   &
               4.75E-19,3.84E-19,3.12E-19,2.56E-19,1.99E-19,1.51E-19,   &
               1.14E-19,8.53E-20,6.24E-20,4.45E-20,3.05E-20,2.01E-20,   &
               1.29E-20,8.02E-21,4.88E-21,2.90E-21,1.67E-21,9.23E-22,   &
               4.90E-22,2.58E-22,1.33E-22,6.65E-23,3.20E-23,1.55E-23,   &
               6.60E-24,2.60E-24,9.80E-25,3.70E-25,1.50E-25,            &
               (0.0,i=1,104)/)
!
      AF13B1 = (/(0.0,i=1,54),                                          &
               7.53E-20,8.63E-20,9.67E-20,1.06E-19,1.14E-19,            &
               1.21E-19,1.26E-19,1.29E-19,1.28E-19,1.25E-19,1.18E-19,   &
               1.08E-19,9.65E-20,8.31E-20,6.95E-20,5.65E-20,4.32E-20,   &
               3.20E-20,2.35E-20,1.72E-20,1.15E-20,7.48E-21,4.82E-21,   &
               3.05E-21,1.82E-21,1.06E-21,6.12E-22,3.37E-22,1.78E-22,   &
               9.22E-23,4.66E-23,2.35E-23,1.18E-23,5.79E-24,2.80E-24,   &
               1.27E-24,5.43E-25,2.29E-25,(0.0,i=1,111)/)
!
      ACH3BR = (/(0.0,i=1,54),                                          &
               5.11E-19,6.02E-19,6.79E-19,7.38E-19,7.76E-19,            &
               7.89E-19,7.86E-19,7.61E-19,7.18E-19,6.65E-19,6.03E-19,   &
               5.39E-19,4.69E-19,4.00E-19,3.33E-19,2.72E-19,2.17E-19,   &
               1.69E-19,1.27E-19,9.24E-20,6.44E-20,4.39E-20,2.94E-20,   &
               1.93E-20,1.22E-20,7.32E-21,4.21E-21,2.34E-21,1.24E-21,   &
               6.34E-22,3.18E-22,1.57E-22,7.61E-23,3.49E-23,1.49E-23,   &
               5.51E-24,(0.0,i=1,113)/)
!     HOCl: JPL 1992.
!     DATA AHOCL/
!    &  60*0.0,
!    &  0.000E+00,0.000E+00,0.000E+00,0.000E+00,3.240E-20,7.151E-20,
!    &  9.698E-20,1.187E-19,1.434E-19,1.717E-19,2.024E-19,2.357E-19,
!    &  2.691E-19,3.023E-19,3.347E-19,3.610E-19,3.774E-19,3.915E-19,
!    &  4.038E-19,4.014E-19,3.867E-19,3.595E-19,3.239E-19,2.845E-19,
!    &  2.485E-19,2.217E-19,2.033E-19,1.925E-19,1.850E-19,1.781E-19,
!    &  1.691E-19,1.578E-19,1.416E-19,1.219E-19,1.043E-19,8.600E-20,
!    &  6.950E-20,5.540E-20,4.350E-20,3.320E-20,2.480E-20,1.830E-20,
!    &  1.340E-20,9.200E-21,6.100E-21,4.200E-21,2.700E-21,1.500E-21,
!    &  1.500E-21,1.500E-21,1.500E-21,9.400E-22,5.000E-22,4.800E-22,    &
!    &  4.000E-22,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,    &
!    &  83*0.0/
!     HOCL: JPL 1994
      AHOCL = (/                                                        &
         (0.0,i=1,54),                                                  &
         0.00E+00, 0.00E+00, 2.76E-21, 4.49E-20, 6.78E-20, 6.54E-20,    &
         5.79E-20, 5.45E-20, 5.45E-20, 5.62E-20, 6.01E-20, 6.57E-20,    &
         7.57E-20, 8.80E-20, 1.03E-19, 1.19E-19, 1.35E-19, 1.54E-19,    &
         1.72E-19, 1.86E-19, 1.99E-19, 2.06E-19, 2.09E-19, 1.99E-19,    &
         1.83E-19, 1.62E-19, 1.40E-19, 1.17E-19, 9.64E-20, 7.88E-20,    &
         6.44E-20, 5.48E-20, 4.91E-20, 4.70E-20, 4.79E-20, 5.08E-20,    &
         5.40E-20, 5.81E-20, 5.99E-20, 6.02E-20, 5.90E-20, 5.51E-20,    &
         4.90E-20, 4.29E-20, 3.50E-20, 2.87E-20, 2.40E-20, 1.81E-20,    &
         1.50E-20, 1.26E-20, 8.00E-21, 9.50E-21, 8.00E-21, 8.25E-21,    &
         8.00E-21, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00,    &
         (0.0,i=1,89)/)
!
      AHCL = (/                                                         &
        (0.0,i=1,30),                                                   &
        0.000E+00,0.000E+00,2.110E-18,2.810E-18,2.810E-18,3.450E-18,    &
        3.450E-18,3.820E-18,3.820E-18,3.320E-18,3.320E-18,2.480E-18,    &
        2.480E-18,1.630E-18,1.630E-18,1.090E-18,1.090E-18,5.880E-19,    &
        5.880E-19,5.880E-19,3.130E-19,3.130E-19,3.130E-19,1.450E-19,    &
        1.450E-19,1.450E-19,6.180E-20,6.180E-20,2.560E-20,2.560E-20,    &
        2.560E-20,9.800E-21,9.800E-21,9.800E-21,3.900E-21,3.900E-21,    &
        1.400E-21,1.400E-21,5.000E-22,5.000E-22,0.000E+00,0.000E+00,    &
        (0.0,i=1,131)/)
!     Cl2O2: JPL 1992 with tail from Burkholder et al (1990)
!     DATA ACL2O2/
!    &  54*0.0,
!    &  0.000E+00,0.000E+00,1.411E-19,2.323E-18,3.594E-18,3.678E-18,
!    &  3.376E-18,3.103E-18,2.829E-18,2.579E-18,2.352E-18,2.176E-18,
!    &  2.071E-18,2.070E-18,2.187E-18,2.430E-18,2.802E-18,3.326E-18,
!    &  3.950E-18,4.649E-18,5.344E-18,5.921E-18,6.293E-18,6.431E-18,
!    &  6.266E-18,5.865E-18,5.292E-18,4.610E-18,3.932E-18,3.292E-18,
!    &  2.767E-18,2.363E-18,2.024E-18,1.745E-18,1.491E-18,1.263E-18,
!    &  1.046E-18,8.580E-19,6.958E-19,5.610E-19,4.330E-19,3.250E-19,
!    &  2.560E-19,2.020E-19,1.670E-19,1.374E-19,1.210E-19,1.085E-19,
!    &  1.060E-19,9.712E-20,7.800E-20,7.662E-20,5.400E-20,4.638E-20,
!    &  3.700E-20,2.812E-20,2.600E-20,2.075E-20,2.100E-20,1.487E-20,
!    &  0.900E-20,0.600E-20,0.300E-20,0.000E-00,0.000E+00,0.000E+00,
!    &  83*0.0/
!     Cl2O2: JPL 1994
      ACL2O2 = (/                                                       &
         (0.0,i=1,54),                                                  &
         0.00E+00, 0.00E+00, 1.41E-19, 2.32E-18, 3.59E-18, 3.68E-18,    &
         3.38E-18, 3.10E-18, 2.83E-18, 2.58E-18, 2.35E-18, 2.18E-18,    &
         2.07E-18, 2.07E-18, 2.19E-18, 2.43E-18, 2.80E-18, 3.33E-18,    &
         3.95E-18, 4.65E-18, 5.34E-18, 5.92E-18, 6.29E-18, 6.43E-18,    &
         6.27E-18, 5.86E-18, 5.29E-18, 4.61E-18, 3.93E-18, 3.29E-18,    &
         2.77E-18, 2.36E-18, 2.02E-18, 1.75E-18, 1.49E-18, 1.26E-18,    &
         1.05E-18, 8.58E-19, 6.96E-19, 5.61E-19, 4.33E-19, 3.25E-19,    &
         2.56E-19, 2.02E-19, 1.67E-19, 1.37E-19, 1.21E-19, 1.05E-19,    &
         8.20E-20, 6.40E-20, 5.50E-20, 3.96E-20, 3.18E-20, 2.55E-20,    &
         2.05E-20, 1.64E-20, 1.32E-20, 1.05E-20, 8.50E-21, 6.79E-21,    &
         5.40E-21, 4.40E-21, 3.50E-21, 2.80E-21, 2.30E-21, 1.79E-21,    &
         1.50E-21, 1.15E-21, 9.00E-22, 0.00E+00, 0.00E+00, 0.00E+00,    &
         (0.0,i=1,77)/)
!     BrCl: Cox. Ignores long wavelengths.
!     DATA ABRCL/
!    &  96*0.0,
!    &  0.000E+00,0.000E+00,4.500E-20,9.100E-20,1.470E-19,2.080E-19,
!    &  2.720E-19,3.310E-19,3.860E-19,4.240E-19,4.450E-19,4.450E-19,
!    &  4.350E-19,4.120E-19,3.740E-19,3.360E-19,2.950E-19,2.540E-19,
!    &  2.190E-19,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,
!    &  83*0.0/
!     BrCl: Seery and Britten. Overwritten in ACSBRCL.
      ABRCL = (/                                                        &
        (0.0,i=1,92),                                                   &
        5.800E-21,1.029E-20,1.490E-20,2.650E-20,3.820E-20,              &
        6.420E-20,9.020E-20,1.318E-19,1.735E-19,2.231E-19,              &
        2.728E-19,3.177E-19,3.626E-19,3.865E-19,4.104E-19,              &
        4.083E-19,4.062E-19,3.821E-19,3.580E-19,3.253E-19,              &
        2.927E-19,2.610E-19,2.293E-19,2.062E-19,1.830E-19,              &
        1.672E-19,1.513E-19,1.412E-19,1.311E-19,1.232E-19,              &
        1.154E-19,1.089E-19,1.024E-19,9.510E-20,8.790E-20,              &
        7.950E-20,7.110E-20,6.270E-20,5.430E-20,4.720E-20,              &
        4.010E-20,3.420E-20,2.830E-20,                                  &
        (0.0,i=1,68)/)
!     BrONO2: JPL 1992
!     DATA ABRNO3/
!    &  48*0.0,
!    &  0.000E+00,7.282E-18,1.208E-17,1.495E-17,1.414E-17,1.322E-17,
!    &  1.214E-17,1.106E-17,9.945E-18,8.892E-18,7.768E-18,6.476E-18,
!    &  5.190E-18,4.272E-18,3.719E-18,3.300E-18,3.018E-18,2.793E-18,
!    &  2.628E-18,2.484E-18,2.332E-18,2.179E-18,2.060E-18,1.956E-18,
!    &  1.870E-18,1.761E-18,1.553E-18,1.336E-18,1.146E-18,9.806E-19,
!    &  8.408E-19,7.212E-19,6.160E-19,5.247E-19,4.508E-19,3.914E-19,
!    &  3.529E-19,3.262E-19,3.064E-19,2.912E-19,2.762E-19,2.520E-19,
!    &  2.332E-19,2.103E-19,1.897E-19,1.776E-19,1.500E-19,1.400E-19,
!    &  1.200E-19,1.100E-19,1.000E-19,9.500E-20,8.700E-20,8.500E-20,
!    &  7.700E-20,6.925E-20,6.200E-20,5.500E-20,4.900E-20,4.475E-20,
!    &  4.000E-20,4.233E-20,2.900E-20,0.000E+00,0.000E+00,0.000E+00,
!    &  89*0.0/
!     BrONO2: Burkholder
      ABRNO3 = (/                                                       &
         (0.0,i=1,48),                                                  &
         8.42E-18, 1.37E-17, 1.58E-17, 1.50E-17, 1.41E-17, 1.32E-17,    &
         1.22E-17, 1.11E-17, 9.94E-18, 8.88E-18, 5.54E-18, 5.14E-18,    &
         4.73E-18, 4.34E-18, 3.90E-18, 3.49E-18, 3.22E-18, 3.01E-18,    &
         2.84E-18, 2.71E-18, 2.58E-18, 2.46E-18, 2.33E-18, 2.19E-18,    &
         2.03E-18, 1.85E-18, 1.65E-18, 1.44E-18, 1.22E-18, 1.01E-18,    &
         8.31E-19, 6.83E-19, 5.57E-19, 4.72E-19, 4.04E-19, 3.52E-19,    &
         3.20E-19, 2.96E-19, 2.76E-19, 2.60E-19, 2.45E-19, 2.28E-19,    &
         2.10E-19, 1.91E-19, 1.71E-19, 1.52E-19, 1.33E-19, 1.17E-19,    &
         1.03E-19, 9.27E-20, 8.44E-20, 7.78E-20, 7.24E-20, 6.77E-20,    &
         6.31E-20, 5.84E-20, 5.39E-20, 4.88E-20, 4.40E-20, 3.90E-20,    &
         3.42E-20, 2.98E-20, 2.58E-20, 2.24E-20, 1.98E-20, 1.74E-20,    &
         1.58E-20, 1.45E-20, 1.35E-20, 1.28E-20, 1.22E-20, 1.15E-20,    &
         1.08E-20, 1.00E-20, 9.07E-21, 8.00E-21, 7.04E-21, 5.95E-21,    &
         4.92E-21, 3.84E-21, 2.92E-21, 2.17E-21, 1.42E-21, 8.60E-22,    &
         (0.0,i=1,71)/)
!     BrO: JPL 1992 averages
      ABRO = (/                                                         &
        (0.0,i=1,90),                                                   &
        0.000E+00,0.000E+00,2.000E-18,2.590E-18,4.540E-18,3.910E-18,    &
        6.000E-18,7.530E-18,6.280E-18,5.890E-18,5.150E-18,3.990E-18,    &
        2.280E-18,1.720E-18,1.610E-18,9.200E-19,5.100E-19,0.000E+00,    &
        (0.0,i=1,95)/)
!     HOBr: Burkholder
!     DATA AHOBR/
!    &  72*0.0,
!    &  0.000E-20,0.000E-20,0.000E-20,0.000E-20,1.000E-20,2.000E-20,
!    &  3.800E-20,5.550E-20,8.000E-20,1.100E-19,1.420E-19,1.810E-19,
!    &  2.230E-19,2.590E-19,2.920E-19,3.080E-19,3.070E-19,2.920E-19,
!    &  2.550E-19,2.070E-19,1.618E-19,1.475E-19,1.090E-19,7.020E-20,
!    &  3.950E-20,4.190E-20,4.890E-20,5.330E-20,5.690E-20,5.980E-20,
!    &  5.980E-20,5.890E-20,5.700E-20,5.210E-20,4.530E-20,3.810E-20,
!    &  3.020E-20,2.250E-20,1.670E-20,1.080E-20,1.800E-21,2.900E-21,
!    &  89*0.0/
!     HOBr: Rattigan. No tail
!     DATA AHOBR/
!    &   72*0.0,
!    &   0.00E+00, 0.00E+00, 0.00E+00, 0.00E+00, 1.94E-20, 4.79E-20,
!    &   6.02E-20, 7.56E-20, 9.71E-20, 1.23E-19, 1.55E-19, 1.89E-19,
!    &   2.25E-19, 2.58E-19, 2.86E-19, 3.03E-19, 3.10E-19, 2.99E-19,
!    &   2.75E-19, 2.42E-19, 2.01E-19, 1.66E-19, 1.38E-19, 1.18E-19,
!    &   1.08E-19, 1.06E-19, 1.10E-19, 1.15E-19, 1.20E-19, 1.23E-19,
!    &   1.25E-19, 1.22E-19, 1.16E-19, 1.07E-19, 9.60E-20, 8.37E-20,
!    &   7.02E-20, 5.73E-20, 4.48E-20, 3.36E-20, 2.42E-20, 1.67E-20,
!    &   1.11E-20, 6.54E-21, 3.39E-21, 1.61E-21, 0.00E+00, 0.00E+00,
!    &   83*0.0/
!     HOBr: Rattigan. With tail
      AHOBR = (/                                                        &
        (0.0,i=1,72),                                                   &
        0.000E+00,0.000E+00,0.000E+00,0.000E+00,1.938E-20,4.794E-20,    &
        6.024E-20,7.559E-20,9.708E-20,1.235E-19,1.546E-19,1.887E-19,    &
        2.249E-19,2.580E-19,2.860E-19,3.031E-19,3.104E-19,2.991E-19,    &
        2.751E-19,2.418E-19,2.014E-19,1.657E-19,1.383E-19,1.176E-19,    &
        1.082E-19,1.062E-19,1.096E-19,1.150E-19,1.196E-19,1.233E-19,    &
        1.246E-19,1.220E-19,1.162E-19,1.068E-19,9.597E-20,8.370E-20,    &
        7.400E-20,6.220E-20,5.080E-20,4.130E-20,3.270E-20,2.560E-20,    &
        2.040E-20,1.590E-20,1.280E-20,1.060E-20,9.200E-21,8.400E-21,    &
        7.400E-21,7.100E-21,6.700E-21,6.500E-21,6.100E-21,5.300E-21,    &
        4.900E-21,4.000E-21,3.400E-21,2.800E-21,2.100E-21,1.400E-21,    &
        9.000E-22,5.000E-22,0.000E+00,0.000E+00,0.000E+00,0.000E+00,    &
        (0.0,i=1,65)/)
!
!     OClO
!     Maximum values from JPL 1994. To be scaled by 0.33
!     DATA AOCLO/
!    &  84*0.0,    &
!    &  0.000E+00,0.000E+00,5.940E-19,7.230E-19,8.970E-19,1.050E-18,    &
!    &  1.400E-18,1.890E-18,2.550E-18,3.740E-18,4.960E-18,6.750E-18,    &
!    &  8.930E-18,1.090E-17,1.270E-17,1.420E-17,1.490E-17,1.520E-17,    &
!    &  1.530E-17,1.510E-17,1.450E-17,1.400E-17,1.340E-17,1.260E-17,    &
!    &  1.170E-17,1.090E-17,9.940E-18,8.850E-18,7.690E-18,6.560E-18,    &
!    &  5.570E-18,4.750E-18,4.010E-18,3.280E-18,2.560E-18,1.940E-18,    &
!    &  1.440E-18,1.040E-18,7.170E-19,4.540E-19,2.220E-19,0.000E+00,    &
!    &  77*0.0/
!     Wahner et al. [1987]. Averaged directly onto model grid
      AOCLO = (/                                                        &
        (0.0,i=1,76),                                                   &
        3.76E-19, 3.66E-19, 3.42E-19, 3.44E-19, 3.47E-19, 3.50E-19,     &
        3.70E-19, 3.89E-19, 4.30E-19, 4.78E-19, 5.41E-19, 6.25E-19,     &
        7.38E-19, 8.82E-19, 1.06E-18, 1.31E-18, 1.67E-18, 2.05E-18,     &
        2.48E-18, 2.59E-18, 2.08E-18, 3.58E-18, 4.57E-18, 5.22E-18,     &
        2.51E-18, 5.48E-18, 6.51E-18, 2.56E-18, 6.54E-18, 3.54E-18,     &
        4.76E-18, 4.32E-18, 3.09E-18, 4.26E-18, 2.14E-18, 2.77E-18,     &
        2.21E-18, 1.60E-18, 1.92E-18, 7.75E-19, 1.48E-18, 2.98E-19,     &
        5.18E-19, 5.78E-19, 1.93E-19, 2.57E-19, 7.62E-20, 6.56E-20,     &
        9.04E-20, 5.90E-20, 8.38E-20, 1.31E-20, 0.00E+00, 0.00E+00,     &
        (0.0,i=1,73)/)
!
!     OBrO
      AOBRO = 0.0
!
!     HONO
      AHONO = 0.0
!
!     CH3OOH
      AMHP = (/(0.0,i=1,63),                                            &
            3.96E-19,3.44E-19,3.06E-19,2.70E-19,2.37E-19,               &
            2.10E-19,1.89E-19,1.68E-19,1.48E-19,1.32E-19,               &
            1.16E-19,1.03E-19,8.98E-20,7.97E-20,7.06E-20,               &
            6.20E-20,5.43E-20,4.75E-20,4.13E-20,3.59E-20,               &
            3.12E-20,2.68E-20,2.25E-20,1.86E-20,1.52E-20,               &
            1.25E-20,1.01E-20,8.30E-21,6.90E-21,5.60E-21,               &
            4.40E-21,3.40E-21,2.60E-21,1.90E-21,1.50E-21,               &
            1.10E-21,8.00E-22,6.00E-22,5.00E-22,                        &
            4.00E-22,3.90E-22,                                          &
            (0.0,i=1,99)/)
!
!     CH2O -> H + CHO
      AC2OA = (/(0.0,i=1,76),                                           &
            9.67E-23,1.02E-22,1.82E-22,4.00E-22,7.34E-22,               &
            1.13E-21,1.49E-21,1.96E-21,2.86E-21,4.99E-21,               &
            7.93E-21,1.10E-20,1.39E-20,1.65E-20,1.85E-20,1.98E-20,      &
            2.04E-20,2.03E-20,1.87E-20,1.46E-20,1.10E-20,               &
            8.02E-21,4.96E-21,1.63E-21,                                 &
            (0.0,i=1,103)/)
!
!     CH2O -> H2 + CO
      AC2OB = (/(0.0,i=1,72),                                           &
             6.57E-24,6.12E-23,1.48E-22,2.58E-22,1.94E-22,              &
             2.04E-22,3.55E-22,7.55E-22,1.30E-21,1.90E-21,              &
             2.34E-21,2.79E-21,3.54E-21,4.85E-21,5.98E-21,              &
             6.62E-21,6.75E-21,6.63E-21,6.29E-21,5.94E-21,              &
             5.62E-21,5.33E-21,5.58E-21,5.83E-21,6.76E-21,              &
             9.90E-21,1.10E-20,1.02E-20,7.87E-21,3.95E-21,              &
             1.62E-21,4.37E-22,4.80E-23,6.75E-24,                       &
             (0.0,i=1,97)/)
!
!     CH3CCl3
      AMCFM = (/(0.0,i=1,51),                                           &
        2.49E-18,2.23E-18,1.98E-18,1.73E-18,1.49E-18,                   &
        1.28E-18,1.08E-18,8.95E-19,7.29E-19,5.79E-19,                   &
        4.53E-19,3.50E-19,2.61E-19,1.85E-19,1.25E-19,                   &
        8.26E-20,5.36E-20,3.40E-20,2.30E-20,1.47E-20,                   &
        9.11E-21,5.61E-21,3.44E-21,2.04E-21,1.13E-21,                   &
        (0.0,i=1,127)/)
!
      ACOF2 = (/                                                        &
        (0.0,i=1,50),                                                   &
        0.000E+00,5.500E-20,4.800E-20,4.200E-20,3.700E-20,3.100E-20,    &
        2.600E-20,2.100E-20,1.600E-20,1.300E-20,0.950E-20,0.740E-20,    &
        0.520E-20,0.400E-20,0.280E-20,0.200E-20,0.120E-20,0.080E-20,    &
        0.049E-20,0.035E-20,0.024E-20,0.018E-20,(0.0,i=1,131)/)
!
      ACOFCL = (/                                                       &
        (0.0,i=1,50),                                                   &
        0.000E+00,1.560E-19,1.400E-19,1.340E-19,1.290E-19,1.270E-19,    &
        1.250E-19,1.240E-19,1.230E-19,1.250E-19,1.200E-19,1.150E-19,    &
        1.080E-19,9.900E-20,9.000E-20,7.900E-20,6.800E-20,5.800E-20,    &
        4.800E-20,3.800E-20,2.900E-20,2.200E-20,(0.0,i=1,131)/)
!
      END SUBROUTINE fill_spectra

      FUNCTION ISRCHFLT(N,XX,INC,X)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name      Type              Description.
!     XX        Array of real     Monotonic array of length N.
!
!     N         Integer           Length of array XX.
!
!     X         Real              Value whose position in XX is required
!
!     Given an array XX of length N, and given a value X, this returns a
!     value J such that X lies between XX(J) and XX(J+1). XX must be
!     monotonically  decreasing. J=0 or J=N is returned
!     to indicate that X is out of range.
!
!     The table entry J is found by bisection.
!
!     This routine is a substitute for the efficient Cray routine
!
!  v1.1 2/12/2003 Original code, adapted to UM  Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: INC
      REAL, INTENT(IN)    :: XX(N)
      REAL, INTENT(IN)    :: X

! Local variables
      INTEGER :: JL
      INTEGER :: JU
      INTEGER :: JM
      INTEGER :: ISRCHFLT
!
!     Initialize upper & lower limits.
      JL = 0
      JU = N + 1
!
      DO WHILE (JU - JL > 1)
        JM = (JU + JL) / 2
        IF ( (XX(N) > XX(1)) .EQV. (X > XX(JM)) ) THEN
          JL = JM
        ELSE
          JU = JM
        END IF
      END DO
      ISRCHFLT = JL + 1
      END FUNCTION ISRCHFLT

      FUNCTION ISRCHFGT(N,XX,INC,X)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     Name      Type              Description.
!     XX        Array of real     Monotonic array of length N.
!
!     N         Integer           Length of array XX.
!
!     X         Real              Value whose position in XX is required
!
!     Given an array XX of length N, and given a value X, this returns a
!     value J such that X lies between XX(J) and XX(J+1). XX must be
!     monotonically  increasing. J=0 or J=N is returned
!     to indicate that X is out of range.
!
!     The table entry J is found by bisection.
!
!     This routine is a substitute for the efficient Cray routine
!
!  v1.1 2/12/2003 Original code, adapted to UM    Olaf Morgenstern
!
!  Current code owner: Olaf Morgenstern.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: INC
      REAL, INTENT(IN)    :: XX(N)
      REAL, INTENT(IN)    :: X

! Local variables
      INTEGER :: JL
      INTEGER :: JU
      INTEGER :: JM
      INTEGER :: ISRCHFGT
!
!     Initialize upper & lower limits.
      JL = 0
      JU = N + 1

      DO WHILE ( JU-JL > 1 )
        JM = (JU + JL) / 2
        IF ( (XX(N) > XX(1)) .EQV. (X > XX(JM)) ) THEN
          JL = JM
        ELSE
          JU = JM
        END IF
      END DO
      ISRCHFGT = JL + 1
!
      END FUNCTION ISRCHFGT
!====================================================================
      END MODULE ukca_photolib
#endif
