#if defined(A70_1C) || defined(A70_1Z) \
 || defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   This module defines the elements required for radiative forcing
!
! Current Code Owner: J.-C. Thelen
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.2   20/08/04   Original code.
!                        J.-C. Thelen
!
! Code Description:
!   Language: FORTRAN 90
!
!- End of header
!
MODULE CORADOCA
!
      USE MAX_CALLS
      INTEGER, Parameter :: C2C_size=npd_swcall-1
!     Whether to change each quantity in diagnostic calls to radiation.
!     The arrays will normally be of size 1, but dimensioning them as
!     arrays will greatly simplify any changes to have more than 1 call.
      LOGICAL :: C2C_O2(C2C_size), C2C_O3(C2C_size), C2C_CO2(C2C_size),&
       C2C_N2O(C2C_size), C2C_CH4(C2C_size), C2C_CFC11(C2C_size),&
       C2C_CFC12(C2C_size), C2C_C113(C2C_size), C2C_HCFC22(C2C_size),&
       C2C_HFC125(C2C_size), C2C_HFC134(C2C_size),&
       C2C_AEROSOL(C2C_size), C2C_SULPC_D(C2C_size),&
       C2C_SEAS_D(C2C_size), C2C_SOOT_D(C2C_size), C2C_BMB_D(C2C_size),&
       C2C_OCFF_D(C2C_size), &
       C2C_SUN(C2C_size), C2C_VOL(C2C_size), C2C_LAND_S(C2C_size),&
       C2C_ALL(C2C_size), C2C_WMG(C2C_size)
END MODULE CORADOCA
#endif
