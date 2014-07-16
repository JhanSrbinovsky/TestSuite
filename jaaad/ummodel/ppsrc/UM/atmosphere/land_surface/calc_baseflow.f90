
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE CALC_BASEFLOW------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE CALC_BASEFLOW(                                         &
     & SOIL_PTS,SOIL_INDEX,NPNTS,NSHYD                                  &
     &,ZDEPTH,KSZ                                                       &
     &,B,FEXP,TI_MEAN,ZW,STHF,STHU                                      &
     &,WUTOT,TOP_CRIT,QBASE,QBASE_L                                     &
     & )





      IMPLICIT NONE
!
! Description:
!     Calculates methane emissions from wetland area.
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
!


      INTEGER                                                           &
     & NPNTS                                                            &
                           ! IN Number of gridpoints.
     &,NSHYD                                                            &
                           ! IN Number of soil moisture levels.
     &,SOIL_PTS            ! IN Number of soil points.

!   Array arguments with intent(IN) :
      INTEGER                                                           &
     & SOIL_INDEX(NPNTS)   ! IN Array of soil points.

      REAL                                                              &
     & B(NPNTS)                                                         &
                           ! IN Clapp-Hornberger exponent.
     &,FEXP(NPNTS)                                                      &
                           ! IN Decay factor in Sat. Conductivity
!                          !    in water table layer.
     &,TI_MEAN(NPNTS)                                                   &
                           ! IN Mean topographic index.
     &,ZW(NPNTS)                                                        &
                            ! IN   Water table depth (m).
     &,KSZ(NPNTS,0:NSHYD)                                               &
                            ! IN Saturated hydraulic condy
!                           !      for each layer (kg/m2/s).
     &,KSFZ(NPNTS,NSHYD+1)                                              &
                            ! WORK Function of sat. hydraulic
!                           !      conductivity, frozen soil
!                           !      and mean topographic index.
     &,STHF(NPNTS,NSHYD)                                                &
                           ! INOUT Frozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.
     &,STHU(NPNTS,NSHYD)   ! INOUT Unfrozen soil moisture content of
!                          !       each layer as a fraction of
!                          !       saturation.


      REAL                                                              &
     & QBASE(NPNTS)                                                     &
                            ! OUT Total base flow (kg/m2/s).
     &,QBASE_L(NPNTS,NSHYD+1)                                           &
!                           ! OUT Base flow from each layer (kg/m2/s).
     &,TOP_CRIT(NPNTS)                                                  &
                            ! OUT Critical topographic index required
!                           !     to calc surface saturation frac.
     &,WUTOT(NPNTS)         ! OUT UNFROZEN to TOTAL fraction at ZW.

      REAL                                                              &
     & QBASE_MAX(NPNTS)                                                 &
                            ! WORK Max possible base flow (kg/m2/s).
     &,QBASE_MAX_L(NPNTS,NSHYD+1)                                       &
!                           ! WORK Max possible base flow
!                           !     from each level (kg/m2/s).
     &,ZDEPTH(0:NSHYD)                                                  &
                            ! WORK Lower soil layer boundary depth (m).
     &,QBASE_MIN(NPNTS)     ! WORK Residual base flow at zw_max

! Local scalars:
      INTEGER                                                           &
     & I,J                                                              &
                            ! WORK Loop counters.
     &,N                                                                &
                            ! WORK Tile loop counter.                   
     &,ERRORSTATUS          ! WORK Error status for ereport.

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
!---Soil layer thicknesses (m)
!---6 layers => CABLE else revert to MOSES
   REAL,PARAMETER:: DZSOIL(6) =(/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/)


!-----------------------------------------------------------------------
! Initialise TOP_CRIT to maximum value.
! Initialise baseflow components to zero.
!-----------------------------------------------------------------------

      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        TOP_CRIT(I)=TI_MAX
        QBASE_MAX(I)=0.0
        DO N=1,NSHYD
          QBASE_MAX_L(I,N)=0.0
        ENDDO
      ENDDO

        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)

!-----------------------------------------------------------------------
! Calculate layer dependent variable with is dependent on
! effective saturated conductivity:
!-----------------------------------------------------------------------
          DO N=1,NSHYD
            IF(STHF(I,N) <  1.0)THEN
              KSFZ(I,N)=0.5*(KSZ(I,N-1)+KSZ(I,N))                       &
     &          *(1.0-STHF(I,N))**(2.0*B(I)+3.0)*EXP(-TI_MEAN(I))
            ELSE
              KSFZ(I,N)=0.0
            ENDIF
          ENDDO
          IF(STHF(I,NSHYD) <  1.0)THEN
            KSFZ(I,NSHYD+1)=KSZ(I,NSHYD)/FEXP(I)                        &
     &        *(1.0-STHF(I,NSHYD))**(2.0*B(I)+3.0)                      &
     &        *EXP(-TI_MEAN(I))
          ELSE
            KSFZ(I,NSHYD+1)=0.0
          ENDIF


          IF(EXP(-FEXP(I)*(ZW_MAX-ZDEPTH(NSHYD))) >  0.01)THEN
            WRITE(6,*)'maximum water table depth is too low!!!'
            WRITE(6,*)'at soil point',i,fexp(i)
            WRITE(6,*)exp(-fexp(i)*(zw_max-zdepth(NSHYD)))
            WRITE(6,*)'try zw_max>',-log(0.01)/fexp(i)+zdepth(NSHYD)
            ERRORSTATUS=1000
            CALL EREPORT('CALC_BASEFLOW', ERRORSTATUS,                  &
     &        'ERROR ZW_MAX is TOO SMALL')
          ENDIF









!-----------------------------------------------------------------------
! Calculate base flow between maximum allowed water table depth and
! "infinity":
!-----------------------------------------------------------------------
          QBASE_MIN(I)=KSFZ(I,NSHYD+1)                                  &
     &      *EXP(-FEXP(I)*(ZW_MAX-ZDEPTH(NSHYD)))

!-----------------------------------------------------------------------
! Calculate maximum possible and actual base flow for each layer:
!-----------------------------------------------------------------------

          DO N=1,NSHYD

            QBASE_MAX_L(I,N)=KSFZ(I,N)*(ZDEPTH(N)-ZDEPTH(N-1))

            IF(ZW(I) <= ZDEPTH(N-1))                                    &
     &        QBASE_L(I,N)=QBASE_MAX_L(I,N)

            IF(ZW(I) <  ZDEPTH(N).AND.ZW(I) >  ZDEPTH(N-1))THEN
              QBASE_L(I,N)=KSFZ(I,N)*(ZDEPTH(N)-ZW(I))
              IF(STHU(I,N)+STHF(I,N) >  0.0)                            &
     &          WUTOT(I)=STHU(I,N)/(STHU(I,N)+STHF(I,N))
            ENDIF
            IF(N == 1.AND.ZW(I) <  ZDEPTH(N))THEN
              QBASE_L(I,N)=KSFZ(I,N)*(ZDEPTH(N)-ZW(I))
              IF(STHU(I,N)+STHF(I,N) >  0.0)                            &
     &          WUTOT(I)=STHU(I,N)/(STHU(I,N)+STHF(I,N))
            ENDIF

          ENDDO

          QBASE_MAX_L(I,NSHYD+1)=KSFZ(I,NSHYD+1)                        &
     &      -QBASE_MIN(I)

          IF(ZW(I) <= ZDEPTH(NSHYD))THEN
            QBASE_L(I,NSHYD+1)=QBASE_MAX_L(I,NSHYD+1)
          ELSE
            QBASE_L(I,NSHYD+1)=KSFZ(I,NSHYD+1)                          &
     &        *EXP(-FEXP(I)*(ZW(I)-ZDEPTH(NSHYD)))                      &
     &        -QBASE_MIN(I)
            IF(STHU(I,NSHYD)+STHF(I,NSHYD) >  0.0)                      &
     &        WUTOT(I)=STHU(I,NSHYD)/(STHU(I,NSHYD)+STHF(I,NSHYD))
          ENDIF

!-----------------------------------------------------------------------
! Calculate total possible and actual base flow:
!-----------------------------------------------------------------------
          DO N=1,NSHYD+1
            QBASE_L(I,N)=MAX(0.0,QBASE_L(I,N))
            QBASE(I)=QBASE(I)+QBASE_L(I,N)
            QBASE_MAX(I)=QBASE_MAX(I)+QBASE_MAX_L(I,N)
          ENDDO

!-----------------------------------------------------------------------
! Calculate critical topographic index.
!-----------------------------------------------------------------------
          IF(QBASE(I) >  QBASE_MAX(I))QBASE(I)=QBASE_MAX(I)
!Check that QBASE_MAX(I)/QBASE(I) will not underflow.
          IF(QBASE_MAX(I) >  EPSILON(QBASE_MAX(I)).AND.                 &
     &      QBASE(I) >  QBASE_MAX(I)*(EPSILON(QBASE(I))))TOP_CRIT(I)    &
     &      =LOG(QBASE_MAX(I)/QBASE(I))

        ENDDO

      RETURN
      END SUBROUTINE CALC_BASEFLOW
