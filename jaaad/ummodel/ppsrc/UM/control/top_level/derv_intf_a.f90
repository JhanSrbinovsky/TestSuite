
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Subroutine DERV_INTF_A : Calculates Interface array dimensions.
!
! Subroutine Interface :
!
      SUBROUTINE DERV_INTF_A (TOT_LEN_INTFA_P,TOT_LEN_INTFA_U,          &
     &           MAX_INTF_MODEL_LEVELS,MAX_LBCROW_LENGTH,MAX_LBCROWS,   &
     &           N_INTF_A,U_FIELD,U_FIELD_INTFA)

      IMPLICIT NONE
!
! Description : Calculate array dimensions for boundary data output.
!
! Method : Reads in INTFCNSTA namelist to get grid dimensions of
!          interface area. Calculates array dimensions for boundary
!          data. Also sets dimensions to 1 if no interface areas
!          required to prevent zero dynamic allocation.
!
! Current Code Owner : Dave Robinson, NWP
!
! History :
! Version    Date    Comment
! -------    ----    -------
!   4.5    03/08/98  Original Code
!   5.2    10/11/00  Cater for 5.x LBCs. D.Robinson
!   6.0    05/09/03  Add new def for makebc. R.Sempers
!   6.0    04/11/03  Correct CPP directives.
!                       A. A. Dickinson
!
! Code Description :
! Language : FORTRAN 77 + common extensions
!
! Declarations :

!     Arguments
      Integer TOT_LEN_INTFA_P   ! OUT  Total length of interface data
                                !      on P grid
      Integer TOT_LEN_INTFA_U   ! OUT  Total length of interface data
                                !      on U grid
      Integer MAX_INTF_MODEL_LEVELS ! OUT  Max no of lbc levels
      Integer MAX_LBCROW_LENGTH     ! OUT  Max no of lbc row_length
      Integer MAX_LBCROWS           ! OUT  Max no of lbc rows
      Integer N_INTF_A          ! IN   No of interface areas
      Integer U_FIELD           ! IN   Dimension of U_FIELD
      Integer U_FIELD_INTFA     ! OUT  Dimension of U_FIELD for dynamic
                                !      allocation.

! CMAXSIZE maximum sizes for dimensioning arrays
! of model constants whose sizes are configuration dependent. This
! allows constants to be read in from a NAMELIST file and maintain
! the flexibility of dynamic allocation for primary variables. The
! maximum sizes should agree with the maximum sizes implicit in the
! front-end User Interface.

!
!  Model            Modification history:
! version  Date
! 3.2  26/03/93  New COMDECK. Author R.Rawlins
! 3.4  06/08/94: Parameter MAX_NO_OF_SEGS used to dimension addresses
!                in macro-tasked calls to SWRAD, LWRAD & CONVECT.
!                Authors: A.Dickinson, D.Salmond, Reviewer: R.Barnes
! 3.5  22/05/95  Add MAX_N_INTF. D. Robinson
! 4.5  29/07/98  Increase MAX_N_INTF/MAX_N_INTF_A to 8. D. Robinson.
! 5.0  20/04/99  Changes for conversion to C-P C dynamics grid.
!                R. Rawlins
!  6.1   04/08/04  Add diffusion variable max_power     Terry Davies
! 6.2  25/12/05  Add max_updiff_levels/max_sponge_width   Terry Davies

      INTEGER,PARAMETER::max_model_levels = 100 ! Maximum no. of levels

      ! Max levels in boundary layer
      INTEGER,PARAMETER:: max_bl_levels = max_model_levels

      ! Max size of alpha_Cd
      INTEGER,PARAMETER :: max_number_alpha_cds = max_bl_levels

      ! Max no. of levels for pvort output
      INTEGER,PARAMETER :: MAX_REQ_THPV_LEVS = max_model_levels

      ! Max no. 1-2-1 rows in polar filter
      INTEGER,PARAMETER ::  max_121_rows =  8
      ! 0 is used for horizontal diffusion pointer

      ! Max no. of levels (from top) to apply upper level diffusion
      INTEGER,PARAMETER ::  max_updiff_levels = 10

      ! Max size of any sponge zones
      INTEGER,PARAMETER ::  max_sponge_width = 10

      ! Max size of look-up tables for searches
      INTEGER,PARAMETER ::  max_look = 2048

      ! Max no. of atmos interface areas
      INTEGER,PARAMETER :: MAX_N_INTF_A =  8

      ! Max no. of points in LBC      
      INTEGER,PARAMETER :: max_intf_lbcrow_length = 1000
      INTEGER,PARAMETER :: max_intf_lbcrows = 1000
        
      ! Max no. of atmos interface levels
      INTEGER,PARAMETER :: MAX_INTF_LEVELS = max_model_levels

      ! Maximum number of physics segments
      INTEGER,PARAMETER :: MAX_NO_OF_SEGS = 200
      ! MAX_N_INTF/MAX_N_INTF_A to be sorted out in later version
      ! Max no. of interface areas
      INTEGER, PARAMETER :: MAX_N_INTF =  8
! CMAXSIZE end
!*L------------------ COMDECK CINTFA ----------------------------------
!L CMAXSIZE should be called first.
!
!   Contains Variables, Headers and Index blocks for control of
!   generation of boundary information for the limited area model.
!
!   Interfaces to all other models are handled by STASH, and there is
!   no explicit coding written for them in the model.
!
!Interface variables initialised through INTFCNSTA
!namelist read in the interface control routine INTF_CTL.
!L
!L 29/07/98  CINTF comdeck renamed to CINTFA. New arrays LBC_STREAM_A
!L           and LBC_UNIT_NO_A added. INTF_AK/BK/AKH/BKH removed - now
!L           in ARGINFA/TYPINFA. D. Robinson.
!L 10/11/00  5.2 Add Intf_ExtHalo_NS, Intf_ExtHalo_EW, Intf_RimW_Orog,
!L               LBC_ND, LBC_Z_TOP_MODEL, LBC_BL_LEVELS,
!L               LBC_FIRST_R_RHO and LBC_Q_MIN. D.Robinson
!  22/10/01  5.3 Remove Intf_RimW_Orog. D. Robinson
!L 18/09/01  5.3 Add A_INTF_FREQ_MN,A_INTF_FREQ_SC. Peter Clark
!L
      INTEGER                                                           &
     &  INTF_ROW_LENGTH                                                 &
                         ! Interface field row length
     & ,INTF_P_ROWS                                                     &
                         ! Interface field no of rows
     & ,INTF_P_LEVELS                                                   &
                         ! Interface field no of levels
     & ,INTF_Q_LEVELS                                                   &
                         ! Interface field no of wet levels
     & ,INTF_TR_LEVELS                                                  &
                         ! Interface field no of tracer levels
     & ,INTFWIDTHA                                                      &
                         ! Width of interface zone (atmosphere)
     & ,Intf_ExtHalo_NS                                                 &
                         ! Extended Halo in NS direction
     & ,Intf_ExtHalo_EW                                                 &
                         ! Extended Halo in EW direction
     & ,LBC_ND                                                          &
                         ! LBCs for old UM (0) or ND (1)
     & ,A_INTF_START_HR                                                 &
                         ! ) Start and End time in
     & ,A_INTF_FREQ_HR                                                  &
                         ! ) hours, Frequency in h,m,s for which
     & ,A_INTF_FREQ_MN                                                  &
                         ! ) atmosphere interface data
     & ,A_INTF_FREQ_SC                                                  &
                         ! ) is to be generated.
     & ,A_INTF_END_HR                                                   &
                         ! )
     & ,LEN_INTFA_P                                                     &
                         ! Length of interface p field
     & ,LEN_INTFA_U                                                     &
                         ! Length of interface u field
     & ,LEN_INTFA_DATA                                                  &
                         ! Length of interface data
     & ,INTF_PACK                                                       &
                         ! Packing Indicator for boundary data
     & ,LBC_STREAM_A                                                    &
                         ! Output streams in UMUI
     & ,LBC_UNIT_NO_A                                                   &
                         ! Unit Nos for Atmos Boundary Dataset
     & ,LBC_FIRST_R_RHO                                                 &
                         ! First rho level at which height is constant
     & ,LBC_BL_LEVELS                                                   &
                         ! No of Boundary Layer levels
!
! Following 3 variables not in common ; in namelist
     & ,INTF_METH_LEV_CALC(MAX_N_INTF_A)                                &
!                              !Method of calculating Eta level (ETAK)
!                              !from layers (ETAH)
     & ,INTF_MAX_SIG_HLEV(MAX_N_INTF_A)                                 &
!                              !level below which sigma coordinates used
     & ,INTF_MIN_PRS_HLEV(MAX_N_INTF_A)
!                              !level above which pressure coordinates

      REAL                                                              &
     &  INTF_EWSPACE                                                    &
                         ! E-W grid spacing (degrees)
     & ,INTF_NSSPACE                                                    &
                         ! N-S grid spacing (degrees)
     & ,INTF_FIRSTLAT                                                   &
                         ! Latitude of first row (degrees)
     & ,INTF_FIRSTLONG                                                  &
                         ! Longitude of first row (degrees)
     & ,INTF_POLELAT                                                    &
                         ! Real latitude of coordinate pole (degrees)
     & ,INTF_POLELONG                                                   &
                         ! Real longitude of coordinate pole (degrees)
     & ,LBC_Z_TOP_MODEL                                                 &
                         ! Height of top of model
     & ,LBC_Q_MIN                                                       &
                         ! Minimum value for q
!
! VarRes grid spacing
      , LAMBDA_INTF_P(MAX_INTF_LBCROW_LENGTH, MAX_N_INTF_A)             &
      , LAMBDA_INTF_U(MAX_INTF_LBCROW_LENGTH, MAX_N_INTF_A)             &    
      , PHI_INTF_P(MAX_INTF_LBCROWS, MAX_N_INTF_A)                      &
      , PHI_INTF_V(MAX_INTF_LBCROWS, MAX_N_INTF_A)                      &
      , LAMBDA_LBC_P(MAX_INTF_LBCROW_LENGTH)                            &
      , LAMBDA_LBC_U(MAX_INTF_LBCROW_LENGTH)                            &    
      , PHI_LBC_P(MAX_INTF_LBCROWS)                                     &
      , PHI_LBC_V(MAX_INTF_LBCROWS)

! Following variable not in common ; in namelist
      REAL INTF_ETAH(MAX_INTF_LEVELS+1,MAX_N_INTF_A)
!                          !Eta values at model layer boundaries ETAKH

      LOGICAL                                                           &
     &  INTF_VERT_INTERP                                                &
                         ! Switch to request vertical interpolation
     & ,LNEWBND          ! True for initialising new boundary data file

! Switch for variable resolution LBC output
      LOGICAL  INTF_L_VAR_LBC(MAX_N_INTF_A)

! Files for VERTLEVS namelist     
      Character(Len=80) :: INTF_VERTLEVS

! Files for HorzGrid namelist  
      Character(Len=80) :: INTF_HorzGrid(MAX_N_INTF_A)
!*----------------------------------------------------------------------
      COMMON /INTFCTL_ATMOS/                                            &
     &  INTF_EWSPACE(MAX_N_INTF_A)    ,INTF_NSSPACE(MAX_N_INTF_A)       &
     & ,INTF_FIRSTLAT(MAX_N_INTF_A)   ,INTF_FIRSTLONG(MAX_N_INTF_A)     &
     & ,INTF_POLELAT(MAX_N_INTF_A)    ,INTF_POLELONG(MAX_N_INTF_A)      &
     & ,INTF_ROW_LENGTH(MAX_N_INTF_A) ,INTF_P_ROWS(MAX_N_INTF_A)        &
     & ,INTF_P_LEVELS(MAX_N_INTF_A)   ,INTF_Q_LEVELS(MAX_N_INTF_A)      &
     & ,INTF_TR_LEVELS(MAX_N_INTF_A)  ,INTFWIDTHA(MAX_N_INTF_A)         &
     & ,Intf_ExtHalo_NS(Max_N_Intf_A) ,Intf_ExtHalo_EW(Max_N_Intf_A)    &
     & ,LBC_ND(Max_N_Intf_A)                                            &
     & ,A_INTF_START_HR(MAX_N_INTF_A) ,A_INTF_FREQ_HR(MAX_N_INTF_A)     &
     & ,A_INTF_FREQ_MN(MAX_N_INTF_A)  ,A_INTF_FREQ_SC(MAX_N_INTF_A)     &
     & ,A_INTF_END_HR(MAX_N_INTF_A)   ,LEN_INTFA_P(MAX_N_INTF_A)        &
     & ,LEN_INTFA_U(MAX_N_INTF_A)     ,LEN_INTFA_DATA(MAX_N_INTF_A)     &
     & ,LNEWBND(MAX_N_INTF_A)         ,INTF_VERT_INTERP(MAX_N_INTF_A)   &
     & ,INTF_PACK(MAX_N_INTF_A)       ,LBC_STREAM_A(MAX_N_INTF_A)       &
     & ,LBC_UNIT_NO_A(MAX_N_INTF_A)   ,LBC_FIRST_R_RHO(MAX_N_INTF_A)    &
     & ,LBC_BL_LEVELS(MAX_N_INTF_A)   ,LBC_Z_TOP_MODEL(MAX_N_INTF_A)    &
     & ,INTF_VERTLEVS(MAX_N_INTF_A)   ,LBC_Q_MIN                        &
     & ,INTF_L_VAR_LBC                ,INTF_HORZGRID                    &
     & ,LAMBDA_INTF_P                 ,LAMBDA_INTF_U                    &
     & ,PHI_INTF_P                    ,PHI_INTF_V
!---------------------------------------------------------------------

!  Namelist for atmos interface constants
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

!     Local variables
      INTEGER JINTF  !  loop index

!     Read in INTFCNSTA namelist to get output grids for
!     generating boundary data.

      REWIND 5
        LBC_ND (:) = 1
      READ (5,INTFCNSTA)
      REWIND 5

      IF (N_INTF_A >  0) THEN

!       Boundary data to be generated in this run.

        MAX_INTF_MODEL_LEVELS = 0
        MAX_LBCROW_LENGTH     = 0
        MAX_LBCROWS           = 0
        Do JINTF=1,N_INTF_A

          MAX_INTF_MODEL_LEVELS =                                       &
     &        MAX ( MAX_INTF_MODEL_LEVELS , INTF_P_LEVELS(JINTF)+1 )
          MAX_LBCROW_LENGTH =                                           &
     &        MAX ( MAX_LBCROW_LENGTH , INTF_ROW_LENGTH(JINTF)+1 )
          MAX_LBCROWS = MAX ( MAX_LBCROWS , INTF_P_ROWS(JINTF)+1 )  
        ENDDO

!       Check >= 1 to avoid zero dynamic allocation
        MAX_INTF_MODEL_LEVELS = MAX ( MAX_INTF_MODEL_LEVELS , 1)
        MAX_LBCROW_LENGTH     = MAX ( MAX_LBCROW_LENGTH , 1)
        MAX_LBCROWS           = MAX ( MAX_LBCROWS , 1)
        
        TOT_LEN_INTFA_P = 0
        TOT_LEN_INTFA_U = 0
        Do JINTF=1,N_INTF_A

          If (lbc_nd(jintf) == 0) Then  !  Old LBCs.

!         Calculate lengths for interface area JINTF

          LEN_INTFA_P(JINTF) = ( INTF_ROW_LENGTH(JINTF) +               &
     &    INTF_P_ROWS(JINTF) - 2*INTFWIDTHA(JINTF) )                    &
     &    * 2 * INTFWIDTHA(JINTF)
          LEN_INTFA_U(JINTF) = LEN_INTFA_P(JINTF) - 4*INTFWIDTHA(JINTF)

!         Add on to total length

          TOT_LEN_INTFA_P = TOT_LEN_INTFA_P + LEN_INTFA_P(JINTF)
          TOT_LEN_INTFA_U = TOT_LEN_INTFA_U + LEN_INTFA_U(JINTF)

          endif

        ENDDO


!       U_FIELD_INTFA Dimensions COEFF3 & COEFF4 in TYPINFA

        U_FIELD_INTFA = U_FIELD

        if (tot_len_intfa_p == 0 .and. tot_len_intfa_u == 0) Then
          tot_len_intfa_p = 1
          tot_len_intfa_u = 1
          u_field_intfa   = 1
        endif

        write (6,*) ' '
        write (6,*) ' Data lengths calculated in DERV_INTF_A.'
        do jintf=1,n_intf_a
        write (6,*) ' Area no ',jintf,                                  &
     &              ' lbc_nd ',lbc_nd(jintf),                           &
     &              ' len_intfa_p ',len_intfa_p(jintf),                 &
     &              ' len_intfa_u ',len_intfa_u(jintf)
        enddo
        write (6,*) ' n_intf_a ',n_intf_a
        write (6,*) ' tot_len_intfa_p ',tot_len_intfa_p
        write (6,*) ' tot_len_intfa_u ',tot_len_intfa_u
        write (6,*) ' max_intf_model_levels ',max_intf_model_levels
        write (6,*) ' max_lbcrow_length ',max_lbcrow_length
        write (6,*) ' max_lbcrows ',max_lbcrows
        write (6,*) ' u_field_intfa ',u_field_intfa

      ELSE

!       No boundary conditions to be generated.
!       Initialise to prevent zero length dynamic allocation.

        write (6,*) ' n_intf_a ',n_intf_a

        N_INTF_A = 1
        TOT_LEN_INTFA_P = 1
        TOT_LEN_INTFA_U = 1
        MAX_INTF_MODEL_LEVELS = 1
        MAX_LBCROW_LENGTH = 1
        MAX_LBCROWS = 1
        U_FIELD_INTFA = 1

      write (6,*) ' n_intf_a ',n_intf_a
      write (6,*) ' tot_len_intfa_p ',tot_len_intfa_p
      write (6,*) ' tot_len_intfa_u ',tot_len_intfa_u
      write (6,*) ' max_intf_model_levels ',max_intf_model_levels
      write (6,*) ' max_lbcrow_length ',max_lbcrow_length
      write (6,*) ' max_lbcrows ',max_lbcrows
      write (6,*) ' u_field_intfa ',u_field_intfa

      ENDIF

      RETURN
      END SUBROUTINE DERV_INTF_A
