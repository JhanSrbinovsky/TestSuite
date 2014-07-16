!    COMDECK COCTROWS
!    ----------------
! History:
! Version  Date     Comment
! -------  ----     -------
!   6.2  23/11/05   Add extra arrays etc for vertical velocity
!                   fix when using free surface.  C.Harris
      REAL                                                              &
     & TTSEC,SWLDEG,FKMP_GLOBAL(IMT,JMT_GLOBAL)                         &
     &,ZU(IMT,JMT),ZV(IMT,JMT)                                          &
     &,ZUENG(IMT,8,JMT),ZVENG(IMT,8,JMT)                                &
     &,ZCONU(IMT_ZVRT,JMT_ZVRT,N_ZVRT)                                  &
                                       ! contributions to barotropic u
     &,ZCONV(IMT_ZVRT,JMT_ZVRT,N_ZVRT)                                  &
                                       ! and v tendencies
     &,SWZVRT(IMT_ZVRT,JMT_ZVRT,N_ZVRT)                                 &
                                       ! vorticity diagnostics
     &,P(IMT_STREAM,0:JMT_STREAM+1),PB(IMT_STREAM,0:JMT_STREAM+1)
      REAL :: PTD(IMT,JMT)
      REAL :: PTDB(IMT,JMT)                                             &
     &,UBT(IMT_FSF,JMTM1_FSF),VBT(IMT_FSF,JMTM1_FSF)                    &
     &,UBTBBC(IMT_FSF,JMTM1_FSF),VBTBBC(IMT_FSF,JMTM1_FSF)              &
     &,MLD(IMT_IPD_MIX,JMT_IPD_MIX)                                     &
     &,ZTD(IMT_STREAM,JMT_STREAM)
      REAL :: ETA(IMT_FSF,JMT_FSF)
      REAL :: ETAB(IMT_FSF,JMT_FSF)
      REAL :: PMEPASS(IMT,JMT)
      REAL :: PMEPASSB(IMT,JMT)
      REAL :: DUSP(IMT,JMT)
      INTEGER                                                           &
     & ITT
