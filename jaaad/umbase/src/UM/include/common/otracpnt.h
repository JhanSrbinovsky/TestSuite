! Description:
! Common block for Ocean model to include pointers into TA
! for specific Cox tracers; TCO2, alkalinity etc, in accordance
! with the 'Conventionally TCO2' type descriptions in the ppxref
! records.
!
! INDICES FOR TRACERS... TEMPERATURE, SALINITY,
!                        INORGANIC CARBON, ALKALINITY
!                        NUTRIENT,PHYTOPLANKTON,
!                        ZOOPLANKTON, DETRITUS,
!                        TRITIUM, 3H + 3HE, CFC11 TO 13,
!                        CARBON14
!
! These indices are set in the SET_OCN_POINTERS routine and
! are used in the carbon model in TRACER and below, esp. BIOLOGY.
! History:
! Version  Date     Comment
! -------  ----     -------
!  UM4.4   9/97   Original code.   Jonathan Palmer
!    4.5   7/98   Deleted O_ADVECT_SCHEME and put into CNTLOCN.cdk
!                 D.Storkey
!    5.3   9/01   Add new indices for tracers in the different
!                 ocean carbon models.     S.Spall
!    5.5   2/03   Add new indices for new tracers (Helium, Oxygen etc)
!                                                          S.Liddicoat
!
      INTEGER T_TRACER,                                                 &
     &        S_TRACER,                                                 &
     &        TCO2_TRACER,                                              &
     &        ALK_TRACER,                                               &
     &        NUTRIENT_TRACER,                                          &
     &        PHYTO_TRACER,                                             &
     &        ZOO_TRACER,                                               &
     &        BOMBHE3_TRACER,                                           &
     &        TRITIUM_TRACER,                                           &
     &        BOMC14_TRACER,                                            &
     &        CFC11_TRACER,                                             &
     &        CFC12_TRACER,                                             &
     &        C14_TRACER,                                               &
     &        C13_TRACER,                                               &
     &        HE3_TRACER,                                               &
     &        HE4_TRACER,                                               &
     &        O2_TRACER,                                                &

     &        DN_TRACER,                                                &
     &        DC_TRACER,                                                &
     &        SN_TRACER,                                                &
     &        SC_TRACER,                                                &
     &        DIATOM_TRACER,                                            &
     &        DTMSI_TRACER,                                             &
     &        DSI_TRACER,                                               &
     &        SILICATE_TRACER,                                          &
     &        FET_TRACER,                                               &
     &        CHLPHY_TRACER,                                            &
     &        CHLDTM_TRACER,                                            &
     &        REGNUT_TRACER,                                            &
     &        CEX_TRACER

      COMMON /TRACPNT/ T_TRACER,                                        &
     &        S_TRACER,                                                 &
     &        TCO2_TRACER,                                              &
     &        ALK_TRACER,                                               &
     &        NUTRIENT_TRACER,                                          &
     &        PHYTO_TRACER,                                             &
     &        ZOO_TRACER,                                               &
     &        BOMBHE3_TRACER,                                           &
     &        TRITIUM_TRACER,                                           &
     &        BOMC14_TRACER,                                            &
     &        CFC11_TRACER,                                             &
     &        CFC12_TRACER,                                             &
     &        C14_TRACER,                                               &
     &        C13_TRACER,                                               &
     &        HE3_TRACER,                                               &
     &        HE4_TRACER,                                               &
     &        O2_TRACER,                                                &
     &        DN_TRACER,                                                &
     &        DC_TRACER,                                                &
     &        SN_TRACER,                                                &
     &        SC_TRACER,                                                &
     &        DIATOM_TRACER,                                            &
     &        DTMSI_TRACER,                                             &
     &        DSI_TRACER,                                               &
     &        SILICATE_TRACER,                                          &
     &        FET_TRACER,                                               &
     &        CHLPHY_TRACER,                                            &
     &        CHLDTM_TRACER,                                            &
     &        REGNUT_TRACER,                                            &
     &        CEX_TRACER
