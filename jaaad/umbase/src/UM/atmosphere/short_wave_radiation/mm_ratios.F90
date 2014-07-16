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
!  Define the reference mixing ratios for the radiatively active
!  gases for the diagnostic call in order to calculate the
!  radiative forcings.
!
! Method:
!
! Current Code Owner: Jean-Claude Thelen
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 6.2       13/02/06  Original code.  Jean-Claude Thelen
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.
!
! Declarations:
!
! Global variables (#include statements etc):
!
MODULE MM_RATIOS
!     Define the reference mixing ratios for the radiatively active
!      gases in the diagnostic call.
!     Values taken from adutv by WJI 31/3/04.
!
      Real, Parameter :: CO2_MMR_D = 4.34800e-04
      Real, Parameter :: n2o_mix_ratio_d = 4.205e-07
      Real, Parameter :: ch4_mix_ratio_d = 4.461e-07
      Real, Parameter :: O2_MMR_D = 0.2314
      Real, Parameter :: cfc11_mix_ratio_d = 0.0
      Real, Parameter :: cfc12_mix_ratio_d = 0.0
      Real, Parameter :: C113MMR_D    = 0.0
      Real, Parameter :: HCFC22MMR_D  = 0.0
      Real, Parameter :: HFC125MMR_D  = 0.0
      Real, Parameter :: HFC134AMMR_D = 0.0
END MODULE MM_RATIOS
#endif
!
