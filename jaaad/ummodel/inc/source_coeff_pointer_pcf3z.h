!     ------------------------------------------------------------------
!     Module to set pointers to source coefficients.
!
      INTEGER                                                           &
     &    IP_SCF_SOLAR_UP                                               &
!           Pointer to source coeficient for upward solar beam
     &  , IP_SCF_SOLAR_DOWN                                             &
!           Pointer to source coeficient for downward solar beam
     &  , IP_SCF_IR_1D                                                  &
!           Pointer to source coeficient
!           For 1st difference of planckian
     &  , IP_SCF_IR_2D
!           Pointer to source coeficient
!           For 2nd difference of planckian
!
      PARAMETER(                                                        &
     &    IP_SCF_SOLAR_UP=1                                             &
     &  , IP_SCF_SOLAR_DOWN=2                                           &
     &  , IP_SCF_IR_1D=1                                                &
     &  , IP_SCF_IR_2D=2                                                &
     &  )
!

!    ------------------------------------------------------------------
