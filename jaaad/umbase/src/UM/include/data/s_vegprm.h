!include file: s_vegprm.h
! Description:
!   defines vegetation parameters for SCM
!
! Current Code Owner: M Hughes
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 5.?      ??/??/?? Original code. Z Gardner
! 6.1      17/08/04 Add SSFM code   M. Hughes
! 6.2      02/02/06 Remove SSFM.   P.Selwood
! 6.2      13/02/06 Remove MOSES I variables.  A lock
!
! Declarations:

!**** Vegetation parameters ******************************************
      Real                                                               &
!     &   ALBSNC(ntypes)                                                  &
!                                ! Cold deep snow-covered albedo
!     &  ,ALBSNF(ntypes)                                                  &
!                                ! Snow free albedo
!     &  ,CATCH(ntypes)                                                   &
!                                ! Surface/canopy water capacity
!                                !  (kg m^-2)
!     &  ,RESIST(ntypes)                                                  &
!                                ! Stomatal resistance to evaporation
!                                !  (s m^-1)
!     &  ,ROOTDEP(ntypes)                                                 &
!                                ! Root depth (metres)
!     &  ,Z0(ntypes)                                                      &
!                                ! Vegetative roughness length (metres)
!     &  ,VEG_FRAC(ntypes)                                                &
!                                ! Vegetation fraction used in
!                                !  calculation of infiltration rate
     &  ,INFIL_FAC(ntypes)                                               &
                                ! Infiltration enhancement factor
                                !  used in calculation of
                                !  infiltration rate.
     &  ,CANHT(ntypes)                                                   &
                                ! Height of the vegetation canopy (m)
     &  ,LAI(ntypes)            ! Leaf area index of vegetation canopy

      Data ALBSNC                                                        &
     &  /0.23,  0.635, 0.255, 0.53,  0.65,  0.78,  0.80,  0.52,  0.80/

      Data ALBSNF                                                        &
     &  /0.12,  0.179, 0.144, 0.185, 0.194, 0.192, 0.162, 0.192, 0.35/

      Data CATCH                                                         &
     &  /0.74,  0.54,  1.13,  0.69,  0.63,  0.50,  0.58,  0.62,  0.50/

      Data RESIST                                                        &
     &  /128.5,69.0,  84.5,  91.0,  79.0,  68.0, 116.,  108.00,100.  /

      Data ROOTDEP                                                      &
     &  /1.43,  0.605,  .84,  0.83,  0.57,  0.54,  0.16,  0.78,  0.10/

      Data Z0                                                           &
     &  /1.0,   0.041, 0.76,  0.104, 0.02,  0.0164,0.00171,0.127,0.003  &
     &  /
      Data VEG_FRAC                                                     &
     &  /0.95,  0.95,  0.95,  0.95,  0.85,  0.80,  0.40,  0.60,  0.00/

      Data INFIL_FAC                                                    &
     &  /5.73,  2.20,  5.50,  2.73,  1.68,  1.60,  0.70,  3.00,  0.50/

!    Only used by MOSES code but need to be declared for input
!     to BL_INTCT
!---------------------------------------------------------------------
