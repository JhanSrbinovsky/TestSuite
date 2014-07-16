#if defined(ATMOS) || defined(MAKEBC)
!+ COMDECK CNAMINFA
!
!    Description:
!       This COMDECK contains the INTFCNSTA namelist which
!       contains all the variables required to define the grids
!       of the interface areas for Atmosphere Boundary data.
!
!       All variables are set up in the UMUI.
!       All variables are declared in comdeck CINTFA
!
!   History:
!
!   Model    Date     Modification history
!  version
!   4.5    03/08/98   New COMDECK created.
!   5.2    10/11/00   Add Intf_ExtHalo_NS, Intf_ExtHalo_EW,
!                     Intf_RimW_Orog, LBC_ND and LBC_Q_MIN.
!                     Add new namelist VERTLEVS. D.Robinson
!   5.3    22/10/01   Remove Intf_RimW_Orog and Boundary_Layer_Levels
!                     D.Robinson
!   5.3    18/09/01   Add A_INTF_FREQ_MN,A_INTF_FREQ_SC to INTFCNSTA.
!                     Peter Clark
!   6.1    18/08/04   Add OLDVERT namelist for 4.5 LBCs. D Robinson
!
      NAMELIST/INTFCNSTA/                                               &
     &         INTF_EWSPACE,INTF_NSSPACE,INTF_FIRSTLAT,INTF_FIRSTLONG,  &
     &         INTF_POLELAT,INTF_POLELONG,                              &
     &         INTF_ROW_LENGTH,INTF_P_ROWS,INTF_P_LEVELS,INTF_Q_LEVELS, &
     &         INTF_TR_LEVELS,                                          &
     &         INTFWIDTHA, Intf_ExtHalo_NS, Intf_ExtHalo_EW,            &
     &         Intf_Pack, LBC_ND, LBC_Q_MIN,                            &
     &         A_INTF_FREQ_HR,A_INTF_FREQ_MN,A_INTF_FREQ_SC,            &
     &         A_INTF_START_HR,A_INTF_END_HR                            &
     &        ,LBC_Stream_A, Intf_VertLevs                              &
     &        ,INTF_L_VAR_LBC, INTF_HorzGrid
 
      ! --------------------------
      ! HORZGRID Namelist for LBCs
      ! --------------------------
      
      Namelist/Lbcgrids/                                                &
     &          LAMBDA_LBC_P, LAMBDA_LBC_U                              &
     &,         PHI_LBC_P, PHI_LBC_V
      
      ! --------------------------
      ! VERTLEVS Namelist for LBCs
      ! --------------------------

      Integer                                                           &
     &   first_constant_r_rho_level

      Real                                                              &
     &   z_top_of_model                                                 &
     &,  eta_theta (max_model_levels+1)                                 &
     &,  eta_rho (max_model_levels)

      Namelist /vertlevs/                                               &
     &   first_constant_r_rho_level,                                    &
     &   z_top_of_model, eta_theta, eta_rho

      ! -----------------------------
      ! OLDVERT Namelist for 4.5 LBCs
      ! -----------------------------

      Logical Vert_Interp
      Integer Meth_Lev_Calc
      Integer Max_sig_hlev, Min_prs_hlev
      Real Etah (max_model_levels+1)

      Namelist /OLDVERT/ Vert_Interp, Meth_Lev_Calc,                    &
     &                   Max_sig_hlev, Min_prs_hlev, etah

!- End of comdeck CNAMINFA
#endif
