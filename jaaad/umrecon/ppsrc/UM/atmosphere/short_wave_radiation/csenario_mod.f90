! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      MODULE CSENARIO_MOD
      IMPLICIT NONE

      Integer, Parameter :: NWMGHG   = 9                ! Number of well-mixed greenhouse gases
      Integer, Parameter :: NSULPAT  = 2                ! Number of sulphate loading patterns
      Integer, Parameter :: LENSCEN  = 800              ! Maximum length of scenarios
      Integer, Parameter :: NCLMFCGS = NWMGHG + NSULPAT ! Number of such scenarios, made up of:

      Integer, Parameter :: S_CO2     = 1               ! Carbon dioxide (CO2)
      Integer, Parameter :: S_CH4     = 2               ! methane (CH4)
      Integer, Parameter :: S_N2O     = 3               ! nitrous oxide (N2O)
      Integer, Parameter :: S_CFC11   = 4               ! Trichlorofluoromethane (CCl3F, "CFC-11")
      Integer, Parameter :: S_CFC12   = 5               ! dichlorodifluoromethane (CCl2F2, "CFC-12")
      Integer, Parameter :: S_SO4     = 6               ! HadCM2-style anthropogenic sulphate loading pattern
      Integer, Parameter :: S_CFC113  = 8               !
      Integer, Parameter :: S_HCFC22  = 9               !
      Integer, Parameter :: S_HFC125  = 10              !
      Integer, Parameter :: S_HFC134A = 11              !

      END MODULE
