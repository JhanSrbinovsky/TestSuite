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
!    Carry out mode-merging algorithm where average mass for
!    mixed nucl,Aitken,accum modes exceeds mid-point mass of next
!    mode up.
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
      SUBROUTINE UKCA_REMODE(NBOX,ND,MD,MDT,DRYDP,WETDP,VERBOSE,        &
        IMERGE,BUD_AER_MAS,N_MERGE_1D)
!----------------------------------------------------------------------
!
! Purpose
! -------
! Carry out mode-merging algorithm where average mass for
! mixed nucl,Aitken,accum modes exceeds mid-point mass of next
! mode up. Transfer fraction of mode number and mass which is greater
! than this threshold. Calculate fraction using error function to
! evaluate integrals of log-normal functions for number and mass
! to threshold. This then gives the amounts to transfer to next mode.
! Use numerical recipes to evaluate error function.
! Re-calculates ND,MD,MDT according to amounts to transfer.
!
! Inputs
! ------
! NBOX     : Number of grid boxes
! ND       : Aerosol ptcl no. concentration (ptcls per cc)
! MD       : Component median aerosol mass (molecules per ptcl)
! MDT      : Total median aerosol mass (molecules per ptcl)
! DRYDP    : Median dry diameter for particles in size mode (m)
! WETDP    : Median wet diameter for particles in size mode (m)
! VERBOSE  : Switch for printing out test print statements
! IMERGE   : Switch to use mid-pts (=1), edges (2) or dynamic (=3)
!
! Outputs
! -------
! ND       : Modified number concentration for each mode
! MD       : Modified average cpt particle mass for each mode
! MDT      : Modified average total particle mass for each mode
! BUD_AER_MAS : Aerosol mass budget terms
! N_MERGE_1D  : # of merges (grown out of bounds) in each box, each mode
!
! Local variables
! ---------------
! LNRATN   : Log of ratio of threshold diameter to no. median diameter
! LNRATM   : Log of ratio of threshold diameter to mass median diameter
! ERFNUM   : LNRATN/(sqrt(2)*log(sigmag))
! ERFMAS   : LNRATM/(sqrt(2)*log(sigmag))
! FRAC_N   : Fraction of ptcl number which is within bounds
! FRAC_M   : Fraction of ptcl mass which is within bounds
! DELN     : Number concentration to transfer due to mode-merging
! DM       : Cpt mass concentration to transfer due to mode-merging
! DP       : Number median diameter of mode (m)
! DP2      : Volume median diameter of mode (m)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES   : Number of possible aerosol modes
! NCP      : Number of possible aerosol components
! MODE     : Logical variable denoting where mode is defined
! COMPONENT: Logical variable denoting where cpt is defined
! DDPMID   : Mid-point of size mode = exp(0.5*(lndp0+lndp1)) (m)
! MMID     : Ptcl mass with dp=dpmed_g=exp(0.5*(lndp0+lndp1)) (ptcl^-1)
! MFRAC_0  : Initial mass fraction to set when no particles.
! SIGMAG   : Geometric standard deviation for each mode
! DDPLIM0  : Lower limit for dry diameter in mode (m)
! DDPLIM1  : Upper limit for dry diameter in mode (m)
! NUM_EPS  : Value of NEWN below which don't recalculate MD
!                                            or carry out process
! CP_SU    : Index of component where SO4    cpt is stored
! CP_BC    : Index of component where BC     cpt is stored
! CP_OC    : Index of component where 1st OC cpt is stored
! CP_CL    : Index of component where NaCl   cpt is stored
! CP_SO    : Index of component where 2nd OC cpt is stored
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES
!
      IMPLICIT NONE
!
!     Subroutine interface
      INTEGER :: NBOX
      INTEGER :: VERBOSE
      INTEGER :: IMERGE
      INTEGER :: N_MERGE_1D(NBOX,NMODES)
      REAL    :: ND(NBOX,NMODES)
      REAL    :: MD(NBOX,NMODES,NCP)
      REAL    :: MDT(NBOX,NMODES)
      REAL    :: DRYDP(NBOX,NMODES)
      REAL    :: WETDP(NBOX,NMODES)
      REAL    :: BUD_AER_MAS(NBOX,NBUDAER)
!
!     Local variables
      INTEGER :: JL
      INTEGER :: IMODE
      INTEGER :: ICP
      INTEGER :: IIMODE
      REAL    :: DP
      REAL    :: DP_IP1
      REAL    :: DP_THRESH1
      REAL    :: DP_THRESH2
      REAL    :: LNRATN
      REAL    :: ERFNUM
      REAL    :: FRAC_N
      REAL    :: DELN
      REAL    :: LOG2SG
      REAL    :: LNRATM
      REAL    :: ERFMAS
      REAL    :: FRAC_M
      REAL    :: DM(NCP)
      REAL    :: DMT
      REAL    :: DP2
      REAL    :: NEWN
      REAL    :: NEWNP1
      REAL    :: ERF
      REAL    :: DERF

      DO IMODE=1,3 ! only do for soluble nucl, Aitken, accum
       IF(MODE(IMODE)) THEN
        DO JL=1,NBOX

         DP=DRYDP(JL,IMODE)
         DP_IP1=DRYDP(JL,IMODE+1)

! IMERGE=1 --> Apply as originally done in M7 (next mode mid-pt)
         IF(IMERGE == 1) THEN
          DP_THRESH1=DDPMID(IMODE+1)
          DP_THRESH2=DDPMID(IMODE+1)
         ENDIF

! IMERGE=2 --> Apply whenever DP is out of bounds
         IF(IMERGE == 2) THEN
          DP_THRESH1=DDPLIM0(IMODE+1)
          DP_THRESH2=DDPLIM0(IMODE+1)
         ENDIF

! IMERGE=3 --> Keep transferring up top-tail of mode (dynamic)
         IF(IMERGE == 3) THEN
          DP_THRESH1=SQRT(DP*DP_IP1)
          DP_THRESH2=SQRT(DP*DP_IP1)
         ENDIF

! .. 1st threshold determines if mode-merging to take place.
         IF((DP > DP_THRESH1).OR.(IMERGE == 3)) THEN

          IF(ND(JL,IMODE) > NUM_EPS(IMODE)) THEN

           N_MERGE_1D(JL,IMODE)=N_MERGE_1D(JL,IMODE)+1

! .. 2nd threshold determines fraction of number/mass
           LNRATN=LOG(DP_THRESH2/DP)
           ERFNUM=LNRATN/SQRT(2.0)/LOG(SIGMAG(IMODE))

! .. FRAC_N is fraction of number remaining in mode
           FRAC_N=0.5*(1.0+DERF(dble(ERFNUM)))

! .. limit DELN to be max = half # of ptcls
           IF(FRAC_N < 0.5) FRAC_N=0.5

           DELN=ND(JL,IMODE)*(1.0-FRAC_N)

           LOG2SG=LOG(SIGMAG(IMODE))*LOG(SIGMAG(IMODE))
           DP2=EXP(LOG(DP)+3.0*LOG2SG) ! volume median diameter

! .. 2nd threshold determines fraction of number/mass
           LNRATM=LOG(DP_THRESH2/DP2)
           ERFMAS=LNRATM/SQRT(2.0)/LOG(SIGMAG(IMODE))
           FRAC_M=0.5*(1.0+DERF(dble(ERFMAS)))

! .. limit DELM to be max 99.9%
           IF(FRAC_M < 0.001) FRAC_M=0.001

! .. calculate new number concs for mode and larger mode
           NEWN=ND(JL,IMODE)-DELN
           NEWNP1=ND(JL,IMODE+1)+DELN

!-----------------------------------------------------------------------
!
! .. This section updates mode from which merging is taking place

           IF(NEWN > NUM_EPS(IMODE)) THEN

            DO ICP=1,NCP
             IF(COMPONENT(IMODE,ICP)) THEN

! .. calculate cpt mass conc to transfer to next mode (use old no/mass)
              DM(ICP)=MD(JL,IMODE,ICP)*ND(JL,IMODE)*(1.0-FRAC_M)

              IF(IMODE == 1) THEN
               IF((ICP == CP_SU).AND.(NMASMERGSUINTR12 > 0))            &
                      BUD_AER_MAS(JL,NMASMERGSUINTR12)=                 &
                      BUD_AER_MAS(JL,NMASMERGSUINTR12)+DM(ICP)
               IF((ICP == CP_OC).AND.(NMASMERGOCINTR12 > 0))            &
                      BUD_AER_MAS(JL,NMASMERGOCINTR12)=                 &
                      BUD_AER_MAS(JL,NMASMERGOCINTR12)+DM(ICP)
               IF((ICP == CP_SO).AND.(NMASMERGSOINTR12 > 0))            &
                      BUD_AER_MAS(JL,NMASMERGSOINTR12)=                 &
                      BUD_AER_MAS(JL,NMASMERGSOINTR12)+DM(ICP)
              ENDIF

              IF(IMODE == 2) THEN
               IF((ICP == CP_SU).AND.(NMASMERGSUINTR23 > 0))            &
                      BUD_AER_MAS(JL,NMASMERGSUINTR23)=                 &
                      BUD_AER_MAS(JL,NMASMERGSUINTR23)+DM(ICP)
               IF((ICP == CP_BC).AND.(NMASMERGBCINTR23 > 0))            &
                      BUD_AER_MAS(JL,NMASMERGBCINTR23)=                 &
                      BUD_AER_MAS(JL,NMASMERGBCINTR23)+DM(ICP)
               IF((ICP == CP_OC).AND.(NMASMERGOCINTR23 > 0))            &
                      BUD_AER_MAS(JL,NMASMERGOCINTR23)=                 &
                      BUD_AER_MAS(JL,NMASMERGOCINTR23)+DM(ICP)
               IF((ICP == CP_SO).AND.(NMASMERGSOINTR23 > 0))            &
                      BUD_AER_MAS(JL,NMASMERGSOINTR23)=                 &
                      BUD_AER_MAS(JL,NMASMERGSOINTR23)+DM(ICP)
              ENDIF

              IF(IMODE == 3) THEN
               IF((ICP == CP_SU).AND.(NMASMERGSUINTR34 > 0))            &
                      BUD_AER_MAS(JL,NMASMERGSUINTR34)=                 &
                      BUD_AER_MAS(JL,NMASMERGSUINTR34)+DM(ICP)
               IF((ICP == CP_BC).AND.(NMASMERGBCINTR34 > 0))            &
                      BUD_AER_MAS(JL,NMASMERGBCINTR34)=                 &
                      BUD_AER_MAS(JL,NMASMERGBCINTR34)+DM(ICP)
               IF((ICP == CP_OC).AND.(NMASMERGOCINTR34 > 0))            &
                      BUD_AER_MAS(JL,NMASMERGOCINTR34)=                 &
                      BUD_AER_MAS(JL,NMASMERGOCINTR34)+DM(ICP)
               IF((ICP == CP_CL).AND.(NMASMERGSSINTR34 > 0))            &
                      BUD_AER_MAS(JL,NMASMERGSSINTR34)=                 &
                      BUD_AER_MAS(JL,NMASMERGSSINTR34)+DM(ICP)
               IF((ICP == CP_SO).AND.(NMASMERGSOINTR34 > 0))            &
                      BUD_AER_MAS(JL,NMASMERGSOINTR34)=                 &
                      BUD_AER_MAS(JL,NMASMERGSOINTR34)+DM(ICP)
              ENDIF

             ELSE
              DM(ICP)=0.0
             ENDIF
            ENDDO
!
! .. first remove mass to be transferred from mode IMODE
            MDT(JL,IMODE)=0.0
            DO ICP=1,NCP
             IF(COMPONENT(IMODE,ICP)) THEN
              MD(JL,IMODE,ICP)=                                         &
                  (ND(JL,IMODE)*MD(JL,IMODE,ICP)-DM(ICP))/NEWN
              MDT(JL,IMODE)=MDT(JL,IMODE)+MD(JL,IMODE,ICP)
             ELSE
              MD(JL,IMODE,ICP)=0.0
             ENDIF ! COMPONENT(IMODE,ICP)
            ENDDO

! .. now set new number to mode IMODE
            ND(JL,IMODE)=NEWN ! set particle number to new value
!
!-----------------------------------------------------------------------
!
! .. This section updates mode IMODE+1

! .. Update MD for this mode
            MDT(JL,IMODE+1)=0.0
            DO ICP=1,NCP
             IF(COMPONENT(IMODE+1,ICP)) THEN
              MD(JL,IMODE+1,ICP)=                                       &
                  (ND(JL,IMODE+1)*MD(JL,IMODE+1,ICP)+DM(ICP))/NEWNP1
              MDT(JL,IMODE+1)=MDT(JL,IMODE+1)+MD(JL,IMODE+1,ICP)
             ELSE
              MD(JL,IMODE+1,ICP)=0.0
             ENDIF ! COMPONENT(IMODE+1,ICP)
            ENDDO

! .. now set new number to mode IMODE+1
            ND(JL,IMODE+1)=NEWNP1

           ENDIF ! IF NEWNP1>0
!
!----------------------------------------------------------------------

          ELSE

           DO ICP=1,NCP
            IF(COMPONENT(IMODE,ICP)) THEN
             MD(JL,IMODE,ICP)=MMID(IMODE)*MFRAC_0(IMODE,ICP)
            ENDIF
           ENDDO
           MDT(JL,IMODE)=MMID(IMODE)

          ENDIF ! if significant number of particles in lower mode

         ENDIF ! if DP outside limits (mode-merge criterion)

        ENDDO ! JL=1,NBOX

       ENDIF ! IF MODE(IMODE)
      ENDDO ! IMODE=1,3 (only merge for soluble nuc/Ait/acc)

      RETURN
      END SUBROUTINE UKCA_REMODE
#endif
