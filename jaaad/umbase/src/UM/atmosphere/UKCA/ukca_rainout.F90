#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Calculates in-cloud aerosol wet deposition (nucleation-scavenging).
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Current Code Owner:       Graham Mann
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
! Subroutine Interface:
      SUBROUTINE UKCA_RAINOUT(NBOX,ND,MD,FCONV_CONV,CRAIN,DRAIN,        &
        CRAIN_UP,DRAIN_UP,T,DTC,BUD_AER_MAS,INUCSCAV)
!---------------------------------------------------------------------------
!
!     Calculates in-cloud aerosl wet deposition (nucleation-scavenging)
!
!     Includes large- (dynamic) & small- scale (convective) precip.
!
!     Include conversion rates FCONV_CONV and FCONV_DYN which represent
!     fraction of condensate this is converted to rain in 6 hours.
!
!     Currently takes FCONV_CONV as input with FCONV_DYN set constant
!
!     Inputs
!     ------
!     NBOX       : Number of grid boxes
!     ND         : Initial no. concentration of aer mode (ptcls/cm3)
!     MD         : Initial avg cpt mass conc of aer mode (molcules/ptcl)
!     FCONV_CONV : Fraction of box condensate converted to rain in
!                                                   6 hours (convective)
!     CRAIN      : Rain rate for conv. precip. in box (kgm^-2s^-1)
!     DRAIN      : Rain rate for dynamic precip. in box (kgm^-2s^-1)
!     CRAIN_UP   : Rain rate for conv. precip. in box above (kgm^-2s^-1)
!     DRAIN_UP   : Rain rate for dyn. precip. in box above (kgm^-2s^-1)
!     T          : Air temperature (K)
!     DTC        : Chemistry timestep (s)
!     INUCSCAV   : Switch for scheme for removal by nucl scav
!
!     Outputs
!     -------
!     ND           : Updated no. concentration in each mode (ptcls/cm3)
!     BUD_AER_MASS : Aerosol mass budgets
!
!     Local Variables
!     ---------------
!     SWITCH       : Rain created in box? (0=none,1=conv,2=dyn,3=both=3)
!     RSCAV        : Scavenging parameters for each mode
!     FCONV_DYN    : Fraction of box condensate converted to rain in
!                                                   6 hours (dynamic)
!     TAU_CONV_DYN : e-folding timescale for conversion of
!                                        condensate to rain (dynamic)
!     TAU_CONV_CONV: e-folding timescale for conversion of
!                                        condensate to rain (convective)
!     FBOX_CONV    : Gridbox fraction over which convective rain occurs
!     TICE         : Temperature below which insoluble aerosol can act
!                    as ice nucleii and hence be removed
!     DELN         : Change in number conc. due to nucleation-scavenging
!
!     Inputted by module UKCA_MODE_SETUP
!     ----------------------------------
!     NMODES       : Number of possible aerosol modes
!     NCP          : Number of possible aerosol components
!     MODE         : Defines which modes are set
!     COMPONENT    : Defines which components are set in each mode
!     MODESOL      : Defines whether mode is soluble of not (integer)
!     NUM_EPS      : Value of NEWN below which don't recalculate MD or
!                                                    carry out process
!     CP_SU        : Index of component in which H2SO4  cpt is stored
!     CP_BC        : Index of component in which BC     cpt is stored
!     CP_OC        : Index of component in which 1st OC cpt is stored
!     CP_CL        : Index of component in which NaCl   cpt is stored
!     CP_SO        : Index of component in which 2nd OC cpt is stored
!
!     Inputted by module UKCA_SETUP_INDICES
!     -------------------------------------
!     Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES
!
      IMPLICIT NONE
!
! .. Subroutine interface
      INTEGER :: NBOX
      INTEGER :: INUCSCAV
      REAL    :: ND(NBOX,NMODES)
      REAL    :: MD(NBOX,NMODES,NCP)
      REAL    :: FCONV_CONV(NBOX)
      REAL    :: CRAIN(NBOX)
      REAL    :: DRAIN(NBOX)
      REAL    :: CRAIN_UP(NBOX)
      REAL    :: DRAIN_UP(NBOX)
      REAL    :: T
      REAL    :: DTC
      REAL    :: NDNEW
      REAL    :: BUD_AER_MAS(NBOX,NBUDAER)
!
! .. Local variables
      INTEGER :: JL
      INTEGER :: IMODE
      INTEGER :: ICP
      INTEGER :: SWITCH
      REAL    :: RSCAV(NMODES)
      REAL    :: DELN
      REAL    :: TAU_CONV_CONV
      REAL    :: TAU_CONV_DYN
      REAL, PARAMETER :: FCONV_DYN=0.9999
      REAL, PARAMETER :: FBOX_CONV=0.3
      REAL, PARAMETER :: TICE=258.0
!
! ..  Calculate timescale for conversion to dynamic rain assuming
! ..  a factor FCONV_DYN is converted to rain in 6 hours
      TAU_CONV_DYN=(-6.0*3600.0)/(LOG(1.0-FCONV_DYN))
!
      IF(INUCSCAV == 1) RSCAV=(/0.00,0.00,1.00,1.00,0.00,0.00,0.00/)
      IF(INUCSCAV == 2) RSCAV=(/0.10,0.25,0.85,0.99,0.20,0.40,0.40/)
!
      DO JL=1,NBOX
!
       SWITCH=0
       IF(CRAIN(JL) > CRAIN_UP(JL))THEN
            SWITCH=1
       ENDIF
       IF(DRAIN(JL) > DRAIN_UP(JL))THEN
        IF(SWITCH == 0) SWITCH=2
        IF(SWITCH == 1) SWITCH=3
       ENDIF
!
! .. Switch > 0 only if rain is FORMED in that level
! .. Convective rain = 1, Dynamic rain = 2
! .. No rain occurs at the top level of the atm
!
       IF(SWITCH > 0)THEN
!
        DO IMODE=1,NMODES
         IF(MODE(IMODE)) THEN
          IF(ND(JL,IMODE) > NUM_EPS(IMODE)) THEN

! .. Only apply to soluble modes or insoluble when T<TICE
           IF((MODESOL(IMODE) == 1).or.(T < TICE)) THEN
!
!----------------------------------------------------------------------
!
! This section does removal by small-scale (conv) precipitation
! .. Convective rain : ND -> ND*(1-FCONV_CONV) over 6 hours
!                      only apply over fraction FBOX_CONV
!
            IF(SWITCH == 1) THEN ! if convective rain
             IF(FCONV_CONV(JL) < 1.0) THEN
              TAU_CONV_CONV=(-6.0*3600.0)/(LOG(1.0-FCONV_CONV(JL)))
              DELN=FBOX_CONV*ND(JL,IMODE)*                              &
               (1.0-EXP(-DTC/TAU_CONV_CONV))*RSCAV(IMODE)
             ELSE
              DELN=FBOX_CONV*ND(JL,IMODE)
             ENDIF
            ENDIF ! if small-scale (convective) rain

!-----------------------------------------------------------------------
!
! This section does removal by large-scale (dyn.) precipitation
!
! .. Dynamic rain    : ND -> ND*(1-FCONV_DYN) over 6 hours
!                      apply over all of box
!
            IF((SWITCH == 2).OR.(SWITCH == 3)) THEN ! if dynamic rain
             DELN=ND(JL,IMODE)*                                         &
               (1.0-EXP(-DTC/TAU_CONV_DYN ))*RSCAV(IMODE)
            ENDIF ! if large-scale (dynamic) rain

!-----------------------------------------------------------------------

! .. calculate updated number concentration due to nucleation-scavening
            NDNEW=ND(JL,IMODE)-DELN

            IF(NDNEW > NUM_EPS(IMODE)) THEN

! .. update number concentration following nucleation-scavening
             ND(JL,IMODE)=NDNEW

!-----------------------------------------------------------------------
!
! .. This section stores removal of each cpt mass for budget calculations

             DO ICP=1,NCP
              IF(COMPONENT(IMODE,ICP)) THEN
               IF(ICP == CP_SU) THEN
                IF((IMODE == 1).AND.(NMASNUSCSUNUCSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSUNUCSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSUNUCSOL)+DELN*MD(JL,IMODE,ICP)
                IF((IMODE == 2).AND.(NMASNUSCSUAITSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSUAITSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSUAITSOL)+DELN*MD(JL,IMODE,ICP)
                IF((IMODE == 3).AND.(NMASNUSCSUACCSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSUACCSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSUACCSOL)+DELN*MD(JL,IMODE,ICP)
                IF((IMODE == 4).AND.(NMASNUSCSUCORSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSUCORSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSUCORSOL)+DELN*MD(JL,IMODE,ICP)
               ENDIF
               IF(ICP == CP_BC) THEN
                IF((IMODE == 2).AND.(NMASNUSCBCAITSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCBCAITSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCBCAITSOL)+DELN*MD(JL,IMODE,ICP)
                IF((IMODE == 3).AND.(NMASNUSCBCACCSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCBCACCSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCBCACCSOL)+DELN*MD(JL,IMODE,ICP)
                IF((IMODE == 4).AND.(NMASNUSCBCCORSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCBCCORSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCBCCORSOL)+DELN*MD(JL,IMODE,ICP)
                IF((IMODE == 5).AND.(NMASNUSCBCAITINS > 0))             &
       BUD_AER_MAS(JL,NMASNUSCBCAITINS)=                                &
       BUD_AER_MAS(JL,NMASNUSCBCAITINS)+DELN*MD(JL,IMODE,ICP)
               ENDIF
               IF(ICP == CP_OC) THEN
                IF((IMODE == 1).AND.(NMASNUSCOCNUCSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCOCNUCSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCOCNUCSOL)+DELN*MD(JL,IMODE,ICP)
                IF((IMODE == 2).AND.(NMASNUSCOCAITSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCOCAITSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCOCAITSOL)+DELN*MD(JL,IMODE,ICP)
                IF((IMODE == 3).AND.(NMASNUSCOCACCSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCOCACCSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCOCACCSOL)+DELN*MD(JL,IMODE,ICP)
                IF((IMODE == 4).AND.(NMASNUSCOCCORSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCOCCORSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCOCCORSOL)+DELN*MD(JL,IMODE,ICP)
                IF((IMODE == 5).AND.(NMASNUSCOCAITINS > 0))             &
       BUD_AER_MAS(JL,NMASNUSCOCAITINS)=                                &
       BUD_AER_MAS(JL,NMASNUSCOCAITINS)+DELN*MD(JL,IMODE,ICP)
               ENDIF
               IF(ICP == CP_CL) THEN
                IF((IMODE == 3).AND.(NMASNUSCSSACCSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSSACCSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSSACCSOL)+DELN*MD(JL,IMODE,ICP)
                IF((IMODE == 4).AND.(NMASNUSCSSCORSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSSCORSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSSCORSOL)+DELN*MD(JL,IMODE,ICP)
               ENDIF
               IF(ICP == CP_SO) THEN
                IF((IMODE == 1).AND.(NMASNUSCSONUCSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSONUCSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSONUCSOL)+DELN*MD(JL,IMODE,ICP)
                IF((IMODE == 2).AND.(NMASNUSCSOAITSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSOAITSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSOAITSOL)+DELN*MD(JL,IMODE,ICP)
                IF((IMODE == 3).AND.(NMASNUSCSOACCSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSOACCSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSOACCSOL)+DELN*MD(JL,IMODE,ICP)
                IF((IMODE == 4).AND.(NMASNUSCSOCORSOL > 0))             &
       BUD_AER_MAS(JL,NMASNUSCSOCORSOL)=                                &
       BUD_AER_MAS(JL,NMASNUSCSOCORSOL)+DELN*MD(JL,IMODE,ICP)
               ENDIF

              ENDIF ! if component(imode,icp)
             ENDDO ! loop over components

!----------------------------------------------------------------------

            ENDIF ! IF NDNEW>NUM_EPS

           ENDIF ! if mode is soluble or T<TICE

          ENDIF ! IF ND>NUM_EPS
         ENDIF ! IF MODE is switched on
        ENDDO ! Loop over modes

       ENDIF ! IF RAIN IS PRODUCED IN THIS LEVEL

      ENDDO ! Loop over gridboxes
!
      RETURN
      END SUBROUTINE UKCA_RAINOUT
#endif
