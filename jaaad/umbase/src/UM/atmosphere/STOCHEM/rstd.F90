#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE RSTD(P,Z)
!
! To calculate R from standard atmosphere table of P & R
!
! C.E. Johnson    4/X/01

!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      REAL, INTENT(IN)  :: P        ! Pressure in Pa
      REAL, INTENT(OUT) :: Z        ! metres above surface

      INTEGER, PARAMETER :: NREF=59
      INTEGER :: I

      REAL, DIMENSION(NREF) :: PREF=(/101325., 95461., 89876.,          &
     &         84559., 79501., 74691., 70121., 65780., 61660.,          &
     &         57752., 54048., 50539., 47217., 44075., 41105.,          &
     &         38299., 35651., 33154., 30800., 28584., 26499.,          &
     &         22699., 19399., 16579., 14170., 12111., 10352.,          &
     &         8849.7, 7565.2, 6467.4, 5529.3, 4728.9, 4047.5,          &
     &         3466.8, 2971.7, 2549.2, 2188.3, 1879.9, 1616.1,          &
     &         1390.4, 1197.0, 1031.2,  889.1,  663.4,  498.5,          &
     &          282.9,  215.4,  165.0,  127.1,   98.5,   76.7,          &
     &           59.8,   31.9,   22.0,   10.9,   5.22,   2.38,          &
     &           1.05,  0.446/)

      REAL, DIMENSION(NREF) :: ZREF=(/0., 500., 1000., 1500., 2000.,    &
     &       2500., 3000., 3500., 4000., 4500., 5000., 5500., 6000.,    &
     &       6500., 7000., 7500., 8000., 8500., 9000., 9500.,10000.,    &
     &      11000.,12000.,13000.,14000.,15000.,16000.,17000.,18000.,    &
     &      19000.,20000.,21000.,22000.,23000.,24000.,25000.,26000.,    &
     &      27000.,28000.,29000.,30000.,31000.,32000.,34000.,36000.,    &
     &      38000.,40000.,42000.,44000.,46000.,48000.,50000.,55000.,    &
     &      60000.,65000.,70000.,75000.,80000.,85000./)


      IF (P < PREF(NREF)) THEN
        Write(6,*) 'P out of range: ',P
        Z=ZREF(NREF)
      ELSEIF (P > PREF(1)) THEN
        Write(6,*) 'P out of range: ',P
        Z=ZREF(1)
      ELSE
        I=1
        DO
          IF(P > PREF(I+1)) EXIT
          I=I+1
        ENDDO
!        write(6,*) 'I=',I,PREF(I),PREF(I+1)
        Z=ZREF(I)+(ZREF(I+1)-ZREF(I))*(PREF(I)-P)/(PREF(I)-PREF(I+1))
      ENDIF

      END SUBROUTINE RSTD
#endif
