
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE CALC_FSAT----------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE CALC_FSAT(L_GAMTOT,SOIL_PTS,SOIL_INDEX,                &
     &   NPNTS,TI_MEAN,TI_SIG,WUTOT,TOP_CRIT,GAMTOT,FSAT,FWETL)
!
! Description:
!     Calculates the surface saturation fraction.
!
!
! Documentation :
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.5    03/02/03   New Deck         Nic Gedney
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
      IMPLICIT NONE

      INTEGER                                                           &
     & NPNTS                                                            &
                      !IN No. of land points.
     &,SOIL_PTS                                                         &
     &,SOIL_INDEX(SOIL_PTS)

      REAL                                                              &
     & TI_MEAN(NPNTS)                                                   &
                      !IN Gridbox mean topographic index.
     &,TI_SIG(NPNTS)                                                    &
                      !IN Standard deviation in topographic index.
     &,GAMTOT(NPNTS)                                                    &
                      !IN Integrated complete Gamma function.
!                     !   Unless L_GAMTOT=TRUE, then INOUT
     &,TOP_CRIT(NPNTS)                                                  &
                      !IN Critical topographic index required
!                     !   to calculate the surface saturation fraction.
     &,WUTOT(NPNTS)   !IN UNFROZEN to TOTAL fraction at ZW.

      LOGICAL                                                           &
     & L_GAMTOT        !IN TRUE if need to calculate GAMTOT

      REAL                                                              &
     & FSAT(NPNTS)                                                      &
                      !INOUT Surface saturation fraction.
     &,FWETL(NPNTS)   !INOUT Wetland fraction.

      REAL                                                              &
     & TI_MAX_USE(NPNTS)                                                &
!                      !WORK TI_MAX dependent on ti_sig
     &,DTI_USE(NPNTS)                                                   &
                       !WORK increment which dependent on ti_sig
     &,ALF                                                              &
                       !WORK Parameter in incomplete Gamma function.
     &,ALF_KSAT                                                         &
                       !WORK Parameter in incomplete Gamma function
!                      !     including horizontal ksat variability.
     &,CUM                                                              &
                       !WORK Integrated incomplete Gamma fn.
     &,TI_SC                                                            &
                       !WORK Incremented Topographic index.
!                      !     Scaled by ti_mean
     &,DTI_SC                                                           &
                       !WORK Scale increment in TI
     &,TICR_SC                                                          &
                       !WORK Critical topographic index. (Any value
!                      !     above this gives surface saturation).
!                      !     Scaled by TI_MEAN
     &,TICR_SC_W                                                        &
                       !WORK As above, but for wetland fraction calc.
     &,FBOX_S                                                           &
                       !WORK Fraction of box for integration.
     &,FBOX_W          !WORK Fraction of box for integration.

      INTEGER                                                           &
     & I,J                                                              &
                       ! Loop counter for points.
     &,NTI                                                              &
                       ! Loop counter for topographic index.
     &,MTI                                                              &
                       ! Max loop count for topographic index.
     &,MTI_W           ! Max loop count for topographic index.

! Comdecks
! C_TOPOG start
! 5.5 17/02/03    Required for large-scale hydrology L_TOP code.
!
! Topographic index increment:
      REAL,PARAMETER :: DTI = 4.0
! Maximum topographic index considered:
      REAL,PARAMETER :: TI_MAX = 10.0
! Maximum allowed water table depth (m):
      REAL,PARAMETER :: ZW_MAX = 5.5
! Standard deviation of LOG(Ksat(0)):
      REAL,PARAMETER :: SIGMA_LOGK = 0.0
! Parameter to remove very high water tables
! from the calculated wetland fraction:
      REAL,PARAMETER :: TI_WETL = 1.5
!
! C_TOPOG end

!----------------------------------------------------------------------
! Set up maximum topographic index to integrate to and increment
! factor. These are topography dependent to maximise accuracy and
! minimise runtime:
!----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        TI_MAX_USE(I)=TI_MEAN(I)+5.*TI_SIG(I)
        DTI_USE(I)=DTI/(TI_SIG(I)*TI_SIG(I))
      ENDDO

!----------------------------------------------------------------------
! Calculate the total integral under the Gamma function. Carried
! out via the reconfiguration:
!----------------------------------------------------------------------
      IF(L_GAMTOT)THEN

        DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)


          MTI=NINT(TI_MAX_USE(I)/DTI_USE(I))
          DTI_SC=DTI_USE(I)/TI_MEAN(I)

          ALF=(TI_MEAN(I)/TI_SIG(I))**2
          ALF_KSAT=1./(1./ALF+(SIGMA_LOGK/TI_MEAN(I))**2)

          DO NTI=1,MTI
            TI_SC=(NTI-0.5)*DTI_SC
            GAMTOT(I)=GAMTOT(I)                                         &
     &        +TI_SC**(ALF_KSAT-1.0)*EXP(-ALF_KSAT*TI_SC)
          ENDDO

        ENDDO
      ELSE

!----------------------------------------------------------------------
! Calculate the integral under the incomplete Gamma function for the
! saturation surface fraction:
!----------------------------------------------------------------------
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)

          IF(TOP_CRIT(I) <  TI_MAX_USE(I))THEN

            TICR_SC=TOP_CRIT(I)/TI_MEAN(I)+1.0
            TICR_SC_W=TICR_SC+TI_WETL/TI_MEAN(I)

            IF(TICR_SC*TI_MEAN(I) <  TI_MAX_USE(I))THEN

              ALF=(TI_MEAN(I)/TI_SIG(I))**2
              ALF_KSAT=1./(1./ALF+(SIGMA_LOGK/TI_MEAN(I))**2)

              DTI_SC=DTI_USE(I)/TI_MEAN(I)
              MTI=INT(TICR_SC/DTI_SC)

              CUM=0.0
              DO NTI=1,MTI
                TI_SC=(NTI-0.5)*DTI_SC
                CUM=CUM+TI_SC**(ALF_KSAT-1.0)                           &
     &            *EXP(-ALF_KSAT*TI_SC)
              ENDDO

! Include the fractional increment:
              FBOX_S=(TICR_SC-MTI*DTI_SC)/DTI_SC
              TI_SC=(MTI+0.5*FBOX_S)*DTI_SC
              CUM=CUM+TI_SC**(ALF_KSAT-1.0)                             &
     &          *EXP(-ALF_KSAT*TI_SC)*FBOX_s

              FSAT(I)=1.0-CUM/GAMTOT(I)


!----------------------------------------------------------------------
! Calculate the integral under the incomplete Gamma function for the
! wetland fraction:
!----------------------------------------------------------------------
              IF(TICR_SC_W*TI_MEAN(I) <  TI_MAX_USE(I))THEN

                MTI_W=INT(TICR_SC_W/DTI_SC)

! Include the fractional increment:
                IF(MTI_W == MTI)THEN
                  FBOX_W=(TICR_SC_W-MTI*DTI_SC)/DTI_SC
                ELSE
                  FBOX_W=1.0
                ENDIF

                TI_SC=(MTI+0.5*FBOX_W)*DTI_SC
                CUM=CUM+TI_SC**(ALF_KSAT-1.0)                           &
     &            *EXP(-ALF_KSAT*TI_SC)*(FBOX_W-FBOX_S)

                DO NTI=MTI+2,MTI_W
                  TI_SC=(NTI-0.5)*DTI_SC
                    CUM=CUM+TI_SC**(ALF_KSAT-1.0)                       &
     &            *EXP(-ALF_KSAT*TI_SC)
                ENDDO

! Include the fractional increment:
                FBOX_W=(TICR_SC_W-MTI_W*DTI_SC)/DTI_SC
                TI_SC=(MTI_W+0.5*FBOX_W)*DTI_SC
                CUM=CUM+TI_SC**(ALF_KSAT-1.0)                           &
     &          *EXP(-ALF_KSAT*TI_SC)*FBOX_W

                FWETL(I)=-1.0+CUM/GAMTOT(I)+FSAT(I)

              ENDIF

              IF(FWETL(I) <  0.0)FWETL(I)=0.0
              IF(FSAT(I) <  0.0)FSAT(I)=0.0
              IF(FWETL(I) >  FSAT(I))FWETL(I)=FSAT(I)

!----------------------------------------------------------------------
! Assume that in partially frozen water, flow is inhibited, therefore
! the critical flow for wetland area is no longer valid:
!----------------------------------------------------------------------
              FWETL(I)=FWETL(I)*WUTOT(I)+FSAT(I)*(1.0-WUTOT(I))

            ENDIF
          ENDIF

        ENDDO

      ENDIF

      RETURN
      END SUBROUTINE CALC_FSAT

!-----------------------------------------------------------------------
