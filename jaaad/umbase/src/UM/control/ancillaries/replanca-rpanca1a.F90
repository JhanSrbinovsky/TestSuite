#if defined(C82_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine REPLANCA ---------------------------------------------
!LL
!LL Purpose:  Updates ancillary fields as requested in FIELDCODE array.
!LL   Tests whether update is required for each field, allowing for
!LL   dependencies between fields. Uses LOOKUP array to find data for
!LL   appropriate time, reads a record and checks for current data
!LL   type. Reads second record if time interpolation required. Updates
!LL   the field. Under DEF RECON, the interface to the routine is
!LL   modified for use in the reconfiguration rather than the model.
!LL   Under DEF CAL360 the 360 day rather than the Gregorian calender
!LL   is used.
!LL
!LL Level 2 control routine for CRAY YMP
!LL
!LL C.Wilson    <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.1  22/02/93  Changes to allow updating of SLAB ref SST and ice
!LL                 ancillary fields (items 178,179) from SST/ice files.
!LL                 Correct bug if SST updated but not ice fraction.
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL   3.2  15/04/93  Remove misleading warning messages if no time
!LL                  interpolation of SST
!LL   3.3  08/02/94  Modify calls to TIME2SEC/SEC2TIME to output/input
!LL                  elapsed times in days & secs, for portability. TCJ
!LL  3.3  22/11/93  Source term and aerosol ancillary fields added.
!LL   3.3  17/11/93 Initializing Integer UPDATE_MONTHS (N.Farnon)
!LL   3.3  08/12/93  Extra argument for READFLDS. D. Robinson
!LL   3.4  17/06/94  DEF CAL360 replaced by LOGICAL LCAL360
!LL                  Argument LCAL360 passed to SEC2TIME, TIME2SEC
!LL                                                   S.J.Swarbrick
!LL  3.4  20/07/94  Improve time interpolation by using reference time
!LL                  for ancillary updating.   R.T.H.Barnes.
!LL  3.4  05/09/94  Add murk and user ancillary fields.  R.T.H.Barnes.
!LL   3.4  18/05/94  Allow prognostic slabtemp under sea ice.J Thomson
!LL   3.4  31/8/94   More error trapping.        (William Ingram)
!LL  4.0  06/09/95  Only print time interpolation diagnostics when it
!LL                 is really done.  RTHBarnes.
!LL   4.0  08/09/95   Allow time interpolation of ozone fields and
!LL                   cater for zonal/full fields. D. Robinson
!LL   4.0  29/11/95   Set land points to zero for surface currents
!LL                   and Heat convergence fields. D. Robinson
!     4.1  18/06/96   Changes to cope with changes in STASH addressing
!                     Author D.M. Goddard.
!LL   4.1  22/05/96   Replace list of ancillary fields with call to
!LL                   comdeck CANCLSTA. D. Robinson.
!LL   4.2  08/11/96   Initialise fields to ensure haloes contain data
!LL                   for time interpolation for MPP runs. D. Robinson
!     4.1  16/12/96   Check ancillary files for non-constant polar rows,
!                     in reconfiguration only. Correct if LPOLARCHK=T
!                     Author D.M. Goddard
!LL   4.4  07/07/97   Alter SST and ice updating for AMIPII runs
!LL                   R A Stratton
!     4.4  13/11/97   Ancilary fields 72 - 89 no longer used for
!                     for user defined ancillaries. Code altered to
!                     ensure correct treatment of these fields.
!                     Author D.M. Goddard
!LL   4.4  25/07/97   (Reconfiguration only). Prevent failure when
!LL                   non-constant polar values for ancillary files are
!LL                   corrected.                         R. Rawlins
!LL   4.5  22/04/98   Add control of new NH3, soot aerosol emission
!LL                   ancillary fields. Plus minor message changes.
!LL                   R.Rawlins
!     4.5  22/10/98   Increase number of user ancillary fields by
!                     deleting existing four fields 68 - 72 and
!                     adding twenty to end 90 - 109
!                     Author D.M. Goddard
!LL   4.5  22/01/98   Correct level of second field read in for time
!LL                   interpolation.  D. Robinson
!LL   5.2  18/10/00   Account for U_FIELD and V_FIELD being different
!LL                   sizes when reading currents into D1 which are
!LL                   used, for example, by the slab model sea ice
!LL                   code.                               K.D.Williams
!LL   5.2  27/09/00   Correct setting of Artic ice depth for SN grid.
!LL   5.3  04/10/01   Removed land masking for ocn currents K.Williams
!     5.3     06/01   Add coastal tiling. Nic Gedney
!     5.4  27/08/02   Removed code which used slab sea ice fractional
!                     time ancil when slab SSTs were updated K.Williams
!LL   5.4  04/09/02   Allow time interpolation of effective vegetation
!LL                   parameter ancillaries.  M. Best.
!     5.5  05/02/03   Add control of biomass smoke emissions and mineral
!                     dust parent soil properties ancillary fields.
!                                                          P Davison
!     6.1  07/04/04   Add control for seawater DMS concentration.
!                                                            A. Jones
!     6.1  08/11/04   Add check for River Routing fields. R.Sharp
!     6.2  22/08/05   Fix operators and remove RECON. P.Selwood.
!     6.2  10/03/06   Allow for a single updating of soil moisture.
!                     Send level info back up. Clive Jones
!LL
!LL Programing standard : UMDP no 3, version no 2, dated 07/09/90
!LL
!LL Logical component covered : C71
!LL
!LL System task : C7
!LL
!LL   External Documentation: UMDP no C7
!LL
!LLEND-------------------------------------------------------------

       SUBROUTINE REPLANCA(I_YEAR,I_MONTH,I_DAY,I_HOUR,                 &
     &                     I_MINUTE,I_SECOND,I_DAY_NUMBER,              &
     &                     ANCIL_REFTIME,OFFSET_STEPS,                  &
     &                     P_FIELD,P_ROWS,U_FIELD,V_FIELD,D1,LAND,      &
     &                     A_STEP,LAND_FIELD,STEPS_PER_HR,              &
     &                     FLAND_CTILE,                                 &
     &                     TSTAR_LAND_CTILE,TSTAR_SEA_CTILE,            &
     &                     TSTAR_SICE_CTILE,                            &
     &                     ICE_FRACTION,TSTAR,TSTAR_ANOM,               &
     &                     SM_LEVELS,DZ_SOIL,SMC_UPDATED,               &
     &                     NS_SPACE,FIRST_LAT,                          &
     &                     LEN1_LOOKUP,LEN_FIXHD,LEN_INTHD,             &
     &                     LEN_REALHD,LEN_D1,FIXHD,INTHD,REALHD,        &
     &                     LOOKUP,RLOOKUP,FTNANCIL,LOOKUP_START,        &
     &                     NDATASETS,NLOOKUPS,                          &
#include "argppx.h"
     &                     ICODE,CMESSAGE,LCAL360)        ! Intent Out

      IMPLICIT NONE

      LOGICAL LCAL360

      INTEGER                                                           &
     &       I_YEAR,                                                    &
                                ! Curent Model Time
     &       I_MONTH,                                                   &
                                !   "      "     "
     &       I_DAY,                                                     &
                                !   "      "     "
     &       I_HOUR,                                                    &
                                !   "      "     "
     &       I_MINUTE,                                                  &
                                !   "      "     "
     &       I_SECOND,                                                  &
                                !   "      "     "
     &       I_DAY_NUMBER,                                              &
     &       ANCIL_REFTIME(6),                                          &
                                ! Reference time for ancillary updating
     &       OFFSET_STEPS,                                              &
                                ! Offset in timesteps of ref. from basis


     &       A_STEP,LAND_FIELD,STEPS_PER_HR,                            &


     &       P_FIELD,                                                   &
                                ! Size of horizontal fields
     &       P_ROWS,                                                    &
                                !
     &       U_FIELD,                                                   &
                                !   "  "      "         "
     &       V_FIELD,                                                   &
                                !   "  "      "         "
     &       NDATASETS,                                                 &
                                ! Number of ancillary datasets
     &       NLOOKUPS,                                                  &
                                ! Number of lookup tables
     &       LEN_D1             ! Size of primary data array

      INTEGER                                                           &
     &       LEN1_LOOKUP,                                               &
                                ! First dimension of lookup table
     &       LEN_FIXHD,                                                 &
                                ! Length of headers in data sets
     &       LEN_INTHD,                                                 &
                                !
     &       LEN_REALHD,                                                &
                         !
     &       FIXHD(LEN_FIXHD,NDATASETS),                                &
                                          ! Data set headers
     &       INTHD(LEN_INTHD,NDATASETS),                                &
                                          !
     &       LOOKUP(LEN1_LOOKUP,NLOOKUPS),                              &
                                          ! Data set lookup tables
     &       FTNANCIL(NDATASETS),                                       &
                                          ! FTN numbers of data sets
     &       LOOKUP_START(NDATASETS),                                   &
                                          ! Start of lookup tables
!                                         ! referring to each data set.
     &       SM_LEVELS          ! number of soil levels

      REAL                                                              &
     &       D1(LEN_D1),                                                &
                                !INOUT  Primary data array used to hold
!                               !       all fields except TSTAR and
!                               !       ICE_FRACTION
     &       ICE_FRACTION(P_FIELD),                                     &
                                    !INOUT  Ice frac of sea part of grid
!                                   !       box, updated if requested
     &       FLAND_CTILE(LAND_FIELD),                                   &
!                                   !IN  Fractional land on land pts.
     &       FLAND_G(P_FIELD),                                          &
                                 !WORK Frac land over all points.
     &       TSTAR(P_FIELD),                                            &
                                 !INOUT  TSTAR:updated if requested
     &       TSTAR_LAND_CTILE(P_FIELD),                                 &
!                                !INOUT  as above, but for land.
     &       TSTAR_SEA_CTILE(P_FIELD),                                  &
!                                !INOUT  as above, but for open sea.
     &       TSTAR_SICE_CTILE(P_FIELD),                                 &
!                                !INOUT  as above, but for sea-ice.
     &       TSTAR_ANOM(P_FIELD),                                       &
                                 !INOUT  SST anomaly,formed in recon;
                                 !       added if requested in model run
     &       REALHD(LEN_REALHD,NDATASETS),                              &
     &       RLOOKUP(LEN1_LOOKUP,NLOOKUPS)                              &
     &       ,NS_SPACE                                                  &
                               ! NS latitude spacing
     &       ,FIRST_LAT                                                 &
                               ! latitude of first gridpoint
     &       ,DZ_SOIL(SM_LEVELS) !OUT soil thicknesses

      LOGICAL                                                           &
     &       LAND(P_FIELD),                                             &
                                 ! WORK LAND mask
     &       SEA(P_FIELD),                                              &
                                 ! WORK SEA mask
     &       LTSTAR_SICE,                                               &
                                 ! IN TRUE if TSTAR_SICE has been read i
!                                ! from input dump.
!                                ! If FALSE set to TSTAR_SEA.
     &       SMC_UPDATED         ! OUT T if smc updated

      INTEGER                                                           &
     &       ICODE                                                      &
                       ! Return code
     &      ,IOUNIT       !OUT I/O unit passed out in RECON mode

      CHARACTER*(80)                                                    &
     &       CMESSAGE  ! Error message
!*
! Comdecks:------------------------------------------------------------
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "cancila.h"
#include "clookadd.h"
#include "cntlatm.h"
#include "cphyscon.h"
#if defined(MPP)
#include "parvars.h"
#endif

!*L   Subroutines called;

      EXTERNAL                                                          &
     &       TIME2SEC,                                                  &
     &       READFLDS,                                                  &
     &       T_INT,                                                     &

     &       TO_LAND_POINTS,                                            &
     &       SEC2TIME,TIME_DF,                                          &

     &       T_INT_C

!*
!*L   Local real arrays

      REAL                                                              &
     &       ANCIL1(P_FIELD),                                           &
                                ! Buffers to hold values of ancillary
!                               ! data for time interpolation.
     &       ANCIL2(P_FIELD),                                           &
                                !
     &       ANCIL_DATA(P_FIELD),                                       &
                                 ! Field of ancillary data held prior
!                               ! to selective updating.
     &       SNOW_CHANGE(P_FIELD),                                      &
                                  ! Fractional time of change of
!                               ! snow cover
     &       ICE_EXTENT(P_FIELD,2),                                     &
                                   ! Fractional time of change
!                               ! of ice cover
     &       PRES_VALUE(P_FIELD)                                        &
                                 ! Prescribed value of data when
!                               ! controlling field is zero.
     &,      NO_ICE_EXTENT(P_FIELD)                                     &
                                   ! Indicator for no sea ice
!                               ! =0 if ice cover
     &,       TSTAR_LAND(P_FIELD)                                       &
                                 !Temporary store for land surface temp.
     &,       TSTAR_SEA(P_FIELD)                                        &
                                 !as above, but for open sea.
     &,       TSTAR_SICE(P_FIELD)                                       &
                                 !as above, but for sea-ice.
     &,       TSTAR_SSI(P_FIELD) !as above, but for sea mean.
!*
!     Local variables

      INTEGER                                                           &
     &       I,                                                         &
                                !
     &       I1,                                                        &
                                !
     &       I2,                                                        &
                                !
     &       I3,                                                        &
     &       ID,                                                        &
                                !
     &       IM,                                                        &
                                !
     &       IY,                                                        &
                                !
     &       K,                                                         &
     &       L,                                                         &
                                    ! Land index
     &       FIELD,                                                     &
                                ! Current field number.
     &       FILE               !

      INTEGER                                                           &
     &       INTERVAL,                                                  &
                                ! Interval between data times
     &       STEP,                                                      &
                                ! Number of data times skipped.
     &       MONTHS,                                                    &
                                ! Used in calculation of position
!                               ! of data required.
     &       HOURS,                                                     &
                                !
     &       PERIOD,                                                    &
                                ! Period of periodic data
     &       START_MONTH,                                               &
                                !
     &       LEVEL,                                                     &
                                !
     &       NFTIN,                                                     &
                                ! Current FTN number for ancillary field
     &       ANCIL_REF_DAYS,                                            &
                                ! Ancil.reference time in whole days
     &       ANCIL_REF_SECS,                                            &
                                ! Ancil.reference time in extra seconds
     &       DAY,SEC,                                                   &
                                ! Times relative to reference time
     &       DAY1,SEC1,                                                 &
                                ! Times relative to reference time
     &       INCR_SEC,                                                  &
                                ! Increment in sec
     &       LEN                                                        &
     &      ,IEND                                                       &
     &      ,II,ROW_LENGTH,J
      INTEGER                                                           &
     &       I_YEAR1,                                                   &
                                 ! Copy of Curent Model Time year
     &       I_MONTH1,                                                  &
                                 !   "      "     "          month
     &       I_DAY1,                                                    &
                                 !   "      "     "          day
     &       I_HOUR1             !   "      "     "          hour

      INTEGER                                                           &
     &       UPDATE_MONTHS      ! update frequency (months) if Gregorian
      LOGICAL                                                           &
     &       LGREG_MONTHLY      ! True for Gregorian monthly updating
!
! *IF -DEF,CAL360
!
      INTEGER                                                           &
     &       I_YEAR_BASIS,                                              &
                                      ! Basis Model Time
     &       I_MONTH_BASIS,                                             &
                                      !   "     "     "
     &       I_DAY_BASIS,                                               &
                                      !   "     "     "
     &       I_HOUR_BASIS,                                              &
                                      !   "     "     "
     &       I_MINUTE_BASIS,                                            &
                                      !   "     "     "
     &       I_SECOND_BASIS,                                            &
                                      !   "     "     "
     &       I_DAY_NUMBER_BASIS
!
! *ENDIF
!

      INTEGER                                                           &
     &       I_YEAR_REF,                                                &
                                      ! Reference Time
     &       I_MONTH_REF,                                               &
                                      !    "       "
     &       I_DAY_REF,                                                 &
                                      !    "       "
     &       I_HOUR_REF,                                                &
                                      !    "       "
     &       I_MINUTE_REF,                                              &
                                      !    "       "
     &       I_SECOND_REF             !    "       "


      LOGICAL                                                           &
     &       LINTERPOLATE,                                              &
                                ! Indicates whether time
!                               ! interpolation needed.
     &       LT_INT_C,                                                  &
                                ! Indicates use of controlled time
!                               ! interpolation
     &       LMISMATCH,                                                 &
                                ! Used in header checks
     &       LICE_FRACTION,                                             &
                                !
     &       LSNOW_DEPTH,                                               &
                                !
     &       SINGLE_TIME,                                               &
                                ! Indicates that only one time is
!                               ! available in data set
     &       PERIODIC,                                                  &
                                ! Data set is periodic
     &       REGULAR                                                    &
                                ! Interval between data times in
!                               ! dataset is regular in model timesteps.
     &       ,LICE_DEPTH

      REAL                                                              &
     &       ZERO,                                                      &
                              !
     &       TIME1,                                                     &
                              ! Times if data used in time interpolation
     &       TIME2,                                                     &
                              !
     &       TIME                                                       &
                        !Target time for time interpolation
     &       ,LAT_P          ! latitude of point
     LOGICAL LAMIP2X ! Local flag

!L    Internal structure

!L  List of Atmosphere & Slab Ancillary fields.
#include "canclsta.h"

!L  1.  Initialisation for atmosphere
      ICODE=0
      IOUNIT=0
      SMC_UPDATED=.FALSE.
      UPDATE_MONTHS=0
      INCR_SEC = 0

!     Set up surface temperatures:

       IF(L_CTILE)THEN
         DO I=1,P_FIELD
            TSTAR_LAND(I)=TSTAR_LAND_CTILE(I)
            TSTAR_SEA(I)=TSTAR_SEA_CTILE(I)
            TSTAR_SICE(I)=TSTAR_SICE_CTILE(I)
            IF(ICE_FRACTION(I) <= 0.0)THEN
              TSTAR_SSI(I)=TSTAR_SEA(I)
            ELSE
              TSTAR_SSI(I)=ICE_FRACTION(I)*TSTAR_SICE(I)                &
     &          +(1.0-ICE_FRACTION(I))*TSTAR_SEA(I)
            ENDIF
         ENDDO
       ELSE
         DO I=1,P_FIELD
            TSTAR_LAND(I)=TSTAR(I)
            TSTAR_SSI(I)=TSTAR(I)
         ENDDO
       ENDIF


!     Initialise ANCIL1/2. Includes Halos for MPP runs.
      DO I=1,P_FIELD
        ANCIL1(I)=0.0
        ANCIL2(I)=0.0
      ENDDO
!L  1.1 Set logical UPDATE for each ancillary field independently

      DO FIELD=1,NANCIL_FIELDS


        UPDATE(FIELD)=.FALSE.
        IF(STEPS(FIELD) /= 0) THEN
!         UPDATE(FIELD)=MOD(A_STEP,STEPS(FIELD)) == 0
          UPDATE(FIELD)=(MOD(A_STEP+OFFSET_STEPS,STEPS(FIELD)) == 0     &
     &                   .OR.A_STEP == 0)                               &
     &                    .AND.FIELDCODE(1,FIELD) >  0                  &
     &                     .AND.D1_ANCILADD(FIELD) >  1
        END IF

!L  1.05 Copy ancillary updating reference time to local variables
      I_YEAR_REF   = ANCIL_REFTIME(1)
      I_MONTH_REF  = ANCIL_REFTIME(2)
      I_DAY_REF    = ANCIL_REFTIME(3)
      I_HOUR_REF   = ANCIL_REFTIME(4)
      I_MINUTE_REF = ANCIL_REFTIME(5)
      I_SECOND_REF = ANCIL_REFTIME(6)
!L       and convert to reference days & secs
! DEPENDS ON: time2sec
            CALL TIME2SEC(I_YEAR_REF,I_MONTH_REF,I_DAY_REF,             &
     &                    I_HOUR_REF,I_MINUTE_REF,I_SECOND_REF,         &
     &                    0,0,ANCIL_REF_DAYS,ANCIL_REF_SECS,LCAL360)

!
      IF (.NOT. LCAL360) THEN

!L  1.11 Set logical UPDATE for Gregorian calender updates at monthly
!L       or yearly intervals. NB STEPS value set to 1 day in INANCILA
        IF(FIELDCODE(1,FIELD) == 1.OR.FIELDCODE(1,FIELD) == 2) THEN
          MONTHS=I_MONTH+I_YEAR*12-(I_MONTH_REF+I_YEAR_REF*12)
          UPDATE_MONTHS= FIELDCODE(2,FIELD)*                            &
     &     ((3-FIELDCODE(1,FIELD))/2 *12+ 1-(3-FIELDCODE(1,FIELD))/2)
          UPDATE(FIELD)=MOD(MONTHS,UPDATE_MONTHS) == 0.AND.I_DAY == 1
        END IF
      END IF !  (.NOT.LCAL360)
!


      END DO

!L 1.2 Allow for dependencies between fields
! Sea surface temperature must be updated when sea ice is updated

      UPDATE(28)=UPDATE(27).OR.UPDATE(28)

! Both surface current components must be updated together

      UPDATE(30)=UPDATE(30).OR.UPDATE(31)
      UPDATE(31)=UPDATE(30)

!L Select method of time interpolation for SST. The interpolation
!L allows for sea ice if ice data is available at the same times
!L as the temperature data. Otherwise linear interpolation is used.

      LT_INT_C=.TRUE.

      IF(UPDATE(28)) THEN
      IF(FIXHD(10,FILEANCIL(27)) == 0) LT_INT_C=.FALSE.
        IF(LT_INT_C) THEN
        DO I=21,41
          IF(FIXHD(I,FILEANCIL(27)) /= FIXHD(I,                         &
     &      FILEANCIL(28))) THEN
            LT_INT_C=.FALSE.
            WRITE(6,*)' WARNING:controlled time interpolation for SST', &
     &      ' not available: Mismatch in SST and SEA-ICE ancillary data'&
     &     ,' times in FIXED HEADER'
            WRITE(6,*)' position=',I,' SEA-ICE=',FIXHD(I,FILEANCIL(27))
            WRITE(6,*)' position=',I,' SST=',FIXHD(I,FILEANCIL(28))
          END IF
        END DO
        ENDIF
      END IF


! Read in fractional land field
! Set up global fractional land field
         IF(L_CTILE)THEN
           L=0
           DO I=1,P_FIELD
             FLAND_G(I)=0.0
             IF(LAND(I))THEN
               L=L+1
               FLAND_G(I)=FLAND_CTILE(L)
            ENDIF
           ENDDO
         ELSE
           DO I=1,P_FIELD
             IF(LAND(I))THEN
               FLAND_G(I)=1.0
             ELSE
               FLAND_G(I)=0.0
            ENDIF
           ENDDO
         ENDIF
!

        DO I=1,P_FIELD
          SEA(I)=.FALSE.
          IF(FLAND_G(I) <  1.0)SEA(I)=.TRUE.
        ENDDO

!L Loop over ancillary fields(atmosphere)

      DO FIELD=1,NANCIL_FIELDS
       ! Turn this off for all but SST and sea-ice
        LAMIP2X = LAMIPII .and. (FIELD == 27 .or. FIELD == 28 .or. FIELD == 29)

       ! Turn this off for all but SST and sea-ice
        LAMIP2X = LAMIPII .and. (FIELD == 27 .or. FIELD == 28 .or. FIELD == 29)

        LICE_DEPTH=field == 29  ! required for LAMIPII

      IF (UPDATE(FIELD)) THEN  ! (1st level IF)
        FILE=FILEANCIL(FIELD)
        NFTIN=FTNANCIL(FILE)

       IF(LICE_DEPTH.AND.LAMIP2X) THEN

! Uses ice fraction set earlier in field loop.
! WARNING this will fail if the order of ancillary fields is ever
! changed so that ice-depth preceeds ice fraction
! Note : For complete sea ice cover
!        Arctic ice depth    = 2m
!        Antarctic ice depth = 1m
! For ice concentrations less than 1. ice depth is 1 or 2 times conc.
! This results in similar values to those from runs using ancillary
! files containing ice depths set to 1 or 2m.

          ROW_LENGTH=P_FIELD/P_ROWS
          DO I=1,P_ROWS
! work out latitude in radians
#if !defined(MPP)
            LAT_P=FIRST_LAT+NS_SPACE*(I-1)
#else
            LAT_P=FIRST_LAT+NS_SPACE*(I+datastart(2)-Offy-1)
#endif
            DO J=1,ROW_LENGTH
              II=J+(I-1)*ROW_LENGTH
              ANCIL_DATA(II)=0.0
              IF (ICE_FRACTION(II) >  0.0) THEN
                IF (LAT_P >  0.0) THEN   ! Arctic ice depth
                  ANCIL_DATA(II)=2.*ICE_FRACTION(II)
                ELSE                     ! Antarctic ice depth
                  ANCIL_DATA(II)=1.*ICE_FRACTION(II)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
!L     Sea ice thickness
!L       Update over all sea points (all sea ice points are the only
!L       ones strictly required, but this cannot be determined easily)

          DO I=1,P_FIELD
            IF(SEA(I)) THEN
              D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
            END IF
          END DO
       ELSE
!     Update required for field

        WRITE(6,*)'REPLANCA: UPDATE REQUIRED FOR FIELD',FIELD

          IF ( FIXHD(10,FILE)  <   0 .OR. FIXHD(10,FILE)  >   2 ) THEN
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: Error in fixed header(10) of ancillary&
     & file                           '
            RETURN
          ENDIF

!L    Check whether more than one data time available in data set

        SINGLE_TIME=FIXHD(10,FILE) == 0

!L    Set default values for time interpolation

        LINTERPOLATE=.TRUE.
        IF(SINGLE_TIME) THEN
          LINTERPOLATE=.FALSE.
        END IF

        IF (FIELD >  9 .AND. FIELD <  19) THEN
          LINTERPOLATE=.FALSE.
        END IF

!L 2.1 Find position of input record

!L    Default settings of search parameters if only one time present

        IF(SINGLE_TIME) THEN
          STEP=0
        ELSE


          LGREG_MONTHLY=.FALSE.
!
      IF (.NOT. LCAL360) THEN
          IF(FIELDCODE(1,FIELD) == 1.OR.FIELDCODE(1,FIELD) == 2) THEN
            LGREG_MONTHLY=.TRUE.
            UPDATE_MONTHS= FIELDCODE(2,FIELD)*                          &
     &      ((3-FIELDCODE(1,FIELD))/2 *12+ 1-(3-FIELDCODE(1,FIELD))/2)
          END IF
      END IF
!


          PERIODIC=FIXHD(10,FILE) == 2
          REGULAR=.TRUE.

!
      IF (.NOT. LCAL360) THEN
          REGULAR=FIXHD(35,FILE) == 0.AND.FIXHD(36,FILE) == 0
! i.e. data at intervals of days/hours & non-periodic
          IF(PERIODIC) REGULAR=REGULAR.AND.FIXHD(37,FILE) == 0
! i.e. data at intervals of hours & periodic
      END IF
!

!         Error checking on time information.

          IF ( FIXHD(35,FILE)  <   0 .OR.                               &
     &         FIXHD(36,FILE)  <   0 .OR. FIXHD(36,FILE)  >   12 .OR.   &
     & REGULAR .AND. ( FIXHD(37,FILE)  <   0 .OR. FIXHD(37,FILE)  >   31&
     &  .OR. FIXHD(38,FILE)  <   0 .OR. FIXHD(38,FILE)  >   24 ) ) THEN
!           FIXHD(39-40) are not used by REPLANCA.
!           FIXHD(35-37) have already been used if not CAL360.
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: Error in validity time interval given &
     &in ancillary file (FIXHD(35-38))'
            RETURN
          ENDIF

          IF ( FIXHD(21,FILE)  <   0 .AND. .NOT. PERIODIC               &
     &  .OR. .NOT. ( REGULAR .AND. PERIODIC ) .AND.                     &
!    !  If it is REGULAR & PERIODIC more detailed check is applied below
     &     ( FIXHD(22,FILE)  <   0 .OR. FIXHD(22,FILE)  >   12 .OR.     &
     &       FIXHD(23,FILE)  <   0 .OR. FIXHD(23,FILE)  >   31 .OR.     &
     &       FIXHD(24,FILE)  <   0 .OR. FIXHD(24,FILE)  >   24 .OR.     &
     &       FIXHD(25,FILE)  <   0 .OR. FIXHD(25,FILE)  >   60 .OR.     &
     &       FIXHD(26,FILE)  <   0 .OR. FIXHD(26,FILE)  >   60 ) ) THEN
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: Error in first validity time given in &
     & ancillary file (FIXHD(21-26))  '
            RETURN
          ENDIF

          IF(.NOT.PERIODIC) THEN

!L            If data taken from full time series of input data.

! DEPENDS ON: time2sec
            CALL TIME2SEC(I_YEAR,I_MONTH,I_DAY,I_HOUR                   &
     &                    ,I_MINUTE,I_SECOND                            &
     &                    ,ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC        &
     &                    ,LCAL360)


!L Adjust time to middle of updating interval

            IF(.NOT.LGREG_MONTHLY) THEN
              SEC=SEC+STEPS(FIELD)*1800/STEPS_PER_HR

!  If start-up, adjust for offset of reference time from initial time,
!  & update with values for half a period before first standard update.
              IF (A_STEP == 0) THEN
                DAY1 = DAY
                SEC1 = SEC
             INCR_SEC=-3600*MOD(OFFSET_STEPS,STEPS(FIELD))/STEPS_PER_HR
! DEPENDS ON: time_df
                CALL TIME_DF(DAY1,SEC1,0,INCR_SEC,DAY,SEC)
              END IF

            ELSE
              IM=MOD(I_MONTH+UPDATE_MONTHS-1,12) + 1
              IY=I_YEAR+(I_MONTH+UPDATE_MONTHS-1)/12
! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,I_DAY,I_HOUR                          &
     &                    ,I_MINUTE,I_SECOND                            &
     &                    ,ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY1,SEC1      &
     &                    ,LCAL360)
              IF (MOD(DAY+DAY1,2) == 0) THEN
                DAY=(DAY+DAY1)/2
                SEC=(SEC+SEC1)/2
              ELSE
                DAY=(DAY+DAY1-1)/2
                SEC=(SEC+SEC1+86400)/2
              ENDIF
!  If start-up, adjust for offset of reference time from initial time,
!  & update with values for half a period before first standard update.
              IF (A_STEP == 0) THEN
                DAY1 = DAY
                SEC1 = SEC
             INCR_SEC=-3600*MOD(OFFSET_STEPS,STEPS(FIELD))/STEPS_PER_HR
! DEPENDS ON: time_df
                CALL TIME_DF(DAY1,SEC1,0,INCR_SEC,DAY,SEC)
              END IF
            ENDIF


            IF(REGULAR) THEN
!L 2.1.1  Standard cases:360 day calender;
!L 2.1.1  or Gregorian calendar with
!L        interval between data times in days or hours
!L        updating interval may be regular in model timesteps,
!L        or (LGREG_MONTHLY=T) irregular in model timesteps,

              HOURS=SEC/3600+DAY*24
!L FInd time(in hours) of first ancillary data on file
! DEPENDS ON: time2sec
              CALL TIME2SEC(FIXHD(21,FILE),FIXHD(22,FILE),              &
     &                   FIXHD(23,FILE),FIXHD(24,FILE),                 &
     &                   FIXHD(25,FILE),FIXHD(26,FILE),                 &
     &                   ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,         &
     &                   LCAL360)
              HOURS=HOURS-SEC/3600-DAY*24

              IF(HOURS <  0) THEN
                ICODE=400+FIELD
           CMESSAGE='REPLANCA: Current time precedes start time of data'
                RETURN
              END IF

!L FInd interval(in hours) between ancillary data on file
              INTERVAL=FIXHD(35,FILE)*8640+FIXHD(36,FILE)*720+          &
     &                FIXHD(37,FILE)*24+FIXHD(38,FILE)

! Do not interpolate in time if data time exactly matches model time

              IF(MOD(HOURS,INTERVAL) == 0) THEN
                LINTERPOLATE=.FALSE.
              END IF

              STEP=HOURS/INTERVAL
              TIME=REAL(HOURS)
              TIME1=STEP*INTERVAL
              TIME2=(STEP+1)*INTERVAL

            ELSE

!L 2.1.2 Gregorian calender;ancillary data interval is in months or
!L       years,which is irregular in model timesteps.
!L original code is inaccurate for this section - corrected code under
!L LAMIPII makes use of dates in lookup headers
!L For a real calendar year the mid-point of each month is different
!L in terms of its hour and day. The old inaccurate method assumes
!L the hour and day are taken from the fixhd values. These are only
!L usually correct for the first month on the ancillary file.


!L Adjust YMD time to middle of updating interval

              I_YEAR1=I_YEAR
              I_MONTH1=I_MONTH
              I_DAY1=I_DAY
              I_HOUR1=I_HOUR
! DEPENDS ON: sec2time
              CALL SEC2TIME(DAY,SEC,ANCIL_REF_DAYS,ANCIL_REF_SECS,      &
     &                     I_YEAR,I_MONTH,I_DAY,                        &
     &                     I_HOUR,I_MINUTE,I_SECOND,I_DAY_NUMBER,       &
     &                     LCAL360)


!L FInd interval(in months) between ancillary data on file
              INTERVAL=FIXHD(35,FILE)*12+FIXHD(36,FILE)
              MONTHS=I_YEAR*12+I_MONTH
              START_MONTH=FIXHD(21,FILE)*12+FIXHD(22,FILE)
              MONTHS=MONTHS-START_MONTH
!  Check for time within month
           IF (LAMIP2X) THEN   ! corrected code uses pp header
              STEP=MONTHS/INTERVAL
              I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*STEP
              I1=I2+LOOKUP_START(FILE)-1
! Check against day and hour of actual lookup header not first field
              IF((I_DAY*24+I_HOUR) <                                    &
     &           (LOOKUP(3,I1)*24+LOOKUP(4,I1))) THEN
                MONTHS=MONTHS-1
              END IF
           ELSE              ! old less accurate code uses FIXHD
              IF((I_DAY*24+I_HOUR) <                                    &
     &           (FIXHD(23,FILE)*24+FIXHD(24,FILE))) THEN
                MONTHS=MONTHS-1
              END IF
           ENDIF ! LAMIP2X

              IF(MONTHS <  0) THEN
                ICODE=400+FIELD
           CMESSAGE='REPLANCA: Current time precedes start time of data'
                RETURN
              END IF


!L Adjust YMD time back to start of updating interval

              I_YEAR=I_YEAR1
              I_MONTH=I_MONTH1
              I_DAY=I_DAY1
              I_HOUR=I_HOUR1



              STEP=MONTHS/INTERVAL

           IF (LAMIP2X) THEN       ! corrected code
              TIME=REAL(SEC)/3600+REAL(DAY*24)
! correct calculation of dates uses lookup table dates not fixhd date
              I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*STEP
              I1=I2+LOOKUP_START(FILE)-1
              I_YEAR1=lookup(1,i1)
              I_MONTH1=lookup(2,i1)
              I_DAY1=lookup(3,i1)
              I_HOUR1=lookup(4,i1)
! DEPENDS ON: time2sec
              CALL TIME2SEC(I_YEAR1,I_MONTH1,I_DAY1,I_HOUR1,            &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,              &
     &              LCAL360)
              TIME1=REAL(SEC)/3600+REAL(DAY*24)
! I1+1 correct pointer to next field as only one field in ancil file
              I_YEAR1=lookup(1,i1+1)
              I_MONTH1=lookup(2,i1+1)
              I_DAY1=lookup(3,i1+1)
              I_HOUR1=lookup(4,i1+1)
! DEPENDS ON: time2sec
              CALL TIME2SEC(I_YEAR1,I_MONTH1,I_DAY1,I_HOUR1,            &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,              &
     &              LCAL360)
              TIME2=REAL(SEC)/3600+REAL(DAY*24)

           ELSE   ! LAMIP2X test - old inaccurate code using FIXHD
! NB INTERVAL may be > 1 month
              MONTHS=STEP*INTERVAL
! Calculate data times for time interpolation
              TIME=REAL(SEC)/3600+REAL(DAY*24)
              IM=MOD(FIXHD(22,FILE)+MONTHS-1,12)+1
              IY=FIXHD(21,FILE)+(MONTHS+FIXHD(22,FILE)-1)/12
! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),        &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,              &
     &              LCAL360)
              TIME1=REAL(SEC)/3600+REAL(DAY*24)
              IM=MOD(FIXHD(22,FILE)+MONTHS+INTERVAL-1,12)+1
              IY=FIXHD(21,FILE)+(MONTHS+INTERVAL+FIXHD(22,FILE)-1)/12
! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),        &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,              &
     &              LCAL360)
              TIME2=REAL(SEC)/3600+REAL(DAY*24)
           ENDIF     ! end LAMIP2X test

! Do not interpolate in time if data time exactly matches model time

              IF(TIME == TIME1) THEN
                LINTERPOLATE=.FALSE.
              END IF

            ENDIF ! End of REGULAR/not REGULAR

          ELSE  ! PERIODIC data

!L 2.2   If data is taken from ancillary periodic data.

! DEPENDS ON: time2sec
            CALL TIME2SEC(I_YEAR,I_MONTH,I_DAY,I_HOUR,                  &
     &                     I_MINUTE,I_SECOND,                           &
     &                     ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,       &
     &                     LCAL360)


!L Adjust time to middle of updating interval

            IF(.NOT.LGREG_MONTHLY) THEN
              SEC=SEC+STEPS(FIELD)*1800/STEPS_PER_HR

!  If start-up, adjust for offset of reference time from initial time,
!  & update with values for half a period before first standard update.
              IF (A_STEP == 0) THEN
                DAY1 = DAY
                SEC1 = SEC
             INCR_SEC=-3600*MOD(OFFSET_STEPS,STEPS(FIELD))/STEPS_PER_HR
! DEPENDS ON: time_df
                CALL TIME_DF(DAY1,SEC1,0,INCR_SEC,DAY,SEC)
              END IF

            ELSE
              IM=MOD(I_MONTH+UPDATE_MONTHS-1,12) + 1
              IY=I_YEAR+(I_MONTH+UPDATE_MONTHS-1)/12
! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,I_DAY,I_HOUR                          &
     &                    ,I_MINUTE,I_SECOND                            &
     &                    ,ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY1,SEC1      &
     &                    ,LCAL360)
              IF (MOD(DAY+DAY1,2) == 0) THEN
                DAY=(DAY+DAY1)/2
                SEC=(SEC+SEC1)/2
              ELSE
                DAY=(DAY+DAY1-1)/2
                SEC=(SEC+SEC1+86400)/2
              ENDIF
!  If start-up, adjust for offset of reference time from initial time,
!  & update with values for half a period before first standard update.
              IF (A_STEP == 0) THEN
                DAY1 = DAY
                SEC1 = SEC
             INCR_SEC=-3600*MOD(OFFSET_STEPS,STEPS(FIELD))/STEPS_PER_HR
! DEPENDS ON: time_df
                CALL TIME_DF(DAY1,SEC1,0,INCR_SEC,DAY,SEC)
              END IF
            ENDIF


!L Adjust YMD time to middle of updating interval

            I_YEAR1=I_YEAR
            I_MONTH1=I_MONTH
            I_DAY1=I_DAY
            I_HOUR1=I_HOUR
! DEPENDS ON: sec2time
            CALL SEC2TIME(DAY,SEC,ANCIL_REF_DAYS,ANCIL_REF_SECS,        &
     &                     I_YEAR,I_MONTH,I_DAY,                        &
     &                     I_HOUR,I_MINUTE,I_SECOND,I_DAY_NUMBER,       &
     &                     LCAL360)



            IF (REGULAR) THEN
!L 2.2.1 Standard cases:1) 360 day calender, with allowed periods of
!L       1 day, 1 month or 1 year;
!L
!L       2) Gregorian calender with update in hours,and period of
!L       data 1 day.
!L
!L       For both updating interval and number of
!L       data times to be skipped in data set calculated in hours.

              HOURS=SEC/3600+DAY*24
              INTERVAL=FIXHD(35,FILE)*8640+FIXHD(36,FILE)*720+          &
     &                FIXHD(37,FILE)*24+FIXHD(38,FILE)

              PERIOD=INTHD(3,FILE)*INTERVAL

!L   Do not allow non-standard periods
      IF (LCAL360) THEN
              IF(PERIOD /= 8640.AND.PERIOD /= 720.AND.PERIOD /= 24)THEN
                ICODE=600+FIELD
           CMESSAGE='REPLANCA: Non-standard period for periodic data'
                RETURN
              ENDIF
      ELSE
              IF(PERIOD /= 24)THEN
                ICODE=600+FIELD
           CMESSAGE='REPLANCA: Non-standard period for periodic data'
                RETURN
              ENDIF
      ENDIF
              IF(PERIOD == 24)THEN
! Ancillary data interval in hour(s), period is 1 day

                IY=I_YEAR
                IM=I_MONTH
                ID=I_DAY
                IF(I_HOUR <  FIXHD(24,FILE)) HOURS=HOURS+24

              ELSE IF(PERIOD == 720)THEN
! Ancillary data interval in day(s) or hours , period is 1 month

                IY=I_YEAR
                IM=I_MONTH
                ID=FIXHD(23,FILE)
                IF((I_DAY*24+I_HOUR) <                                  &
     &             (FIXHD(23,FILE)*24+FIXHD(24,FILE)))                  &
     &           HOURS=HOURS+720

              ELSE IF(PERIOD == 8640)THEN
! Ancillary data interval in month(s)or days or hours, period is 1 year

                IY=I_YEAR
                IM=FIXHD(22,FILE)
                ID=FIXHD(23,FILE)
                IF((I_MONTH*720+I_DAY*24+I_HOUR) <                      &
     &          (FIXHD(22,FILE)*720+FIXHD(23,FILE)*24+FIXHD(24,FILE)))  &
     &           HOURS=HOURS+8640

              END IF

! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,ID,FIXHD(24,FILE),                    &
     &                     FIXHD(25,FILE),FIXHD(26,FILE),               &
     &                     ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,       &
     &                     LCAL360)
              HOURS=HOURS-SEC/3600-DAY*24

! Do not interpolate in time if data time exactly matches model time

              IF(MOD(HOURS,INTERVAL) == 0) THEN
                LINTERPOLATE=.FALSE.
              END IF
              STEP=HOURS/INTERVAL
              TIME=REAL(HOURS)
              TIME1=STEP*INTERVAL
              TIME2=(STEP+1)*INTERVAL

            ELSE  ! non regular case

!L 2.2.2 Gregorian calender,and data interval is in months,
!L       period is 1 year
!L       Updating interval and number of data times to be skipped
!L       calculated in months.

              TIME=REAL(SEC)/3600+REAL(DAY*24)
              INTERVAL=FIXHD(36,FILE)+FIXHD(35,FILE)*12
              PERIOD=INTHD(3,FILE)*INTERVAL
              IF(PERIOD /= 12)THEN
                ICODE=600+FIELD
           CMESSAGE='REPLANCA: Non-standard period for periodic data'
                RETURN
              ENDIF
!  Difference between date now (month) & first date ancil file (month)
              MONTHS=I_MONTH-FIXHD(22,FILE)
#ifdef ACCESS
              ! required to avoid tripping -check uninit when
              ! writing out hours below
              HOURS=INT(TIME)
#endif

           IF (LAMIP2X) THEN ! correct code to use lookup header dates
! Correctly use day and hour from lookup header not fixhd which
! contains values for first field on ancillary file only.
             step=months/INTERVAL
             I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*step
             I1=I2+LOOKUP_START(FILE)-1
!  Check for time within month - using ppheader information
             IF((I_DAY*24+I_HOUR) <  (lookup(3,i1)*24+lookup(4,i1))) THEN
                 MONTHS=MONTHS-1
             END IF
             IF(MONTHS <  0) THEN
                MONTHS=MONTHS+12
             END IF
! recalculate STEP
             STEP=MONTHS/INTERVAL
! NB INTERVAL may be > 1 month
             MONTHS=STEP*INTERVAL
             IY=I_YEAR
             IM=MOD(FIXHD(22,FILE)+MONTHS-1,12)+1
             IF(IM >  I_MONTH) IY=IY-1
             I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*STEP
             I1=I2+LOOKUP_START(FILE)-1
! DEPENDS ON: time2sec
             CALL TIME2SEC(IY,IM,lookup(3,i1),lookup(4,i1),             &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,LCAL360)
             TIME1=REAL(SEC)/3600+REAL(DAY*24)
!  Calculate  TIME2 for second ancillary data time
!  set IY correctly for time interpolation calculations
             IY=I_YEAR
             IM=MOD(FIXHD(22,FILE)+MONTHS+INTERVAL-1,12)+1
             IF(IM <  I_MONTH) IY=IY+1
             I1=(IM-1)/INTERVAL
             I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*I1
             I1=I2+LOOKUP_START(FILE)-1
! DEPENDS ON: time2sec
             CALL TIME2SEC(IY,IM,lookup(3,i1),lookup(4,i1),             &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,LCAL360)
             TIME2=REAL(SEC)/3600+REAL(DAY*24)

           ELSE   ! original code inaccurate use of FIXHD dates
!  Check for time within month
              IF((I_DAY*24+I_HOUR) <                                    &
     &           (FIXHD(23,FILE)*24+FIXHD(24,FILE))) THEN
                MONTHS=MONTHS-1
              END IF
              IF(MONTHS <  0) THEN
                MONTHS=MONTHS+12
              END IF

              STEP=MONTHS/INTERVAL
! NB INTERVAL may be > 1 month
              MONTHS=STEP*INTERVAL
!  Calculate TIME1 for first ancillary data time
!  set IY correctly for time interpolation calculations
              IY=I_YEAR
              IM=MOD(FIXHD(22,FILE)+MONTHS-1,12)+1
              IF(IM >  I_MONTH) IY=IY-1
! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),        &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,              &
     &              LCAL360)
              TIME1=REAL(SEC)/3600+REAL(DAY*24)
!  Calculate  TIME2 for second ancillary data time
!  set IY correctly for time interpolation calculations
              IY=I_YEAR
              IM=MOD(FIXHD(22,FILE)+MONTHS+INTERVAL-1,12)+1
              IF(IM <  I_MONTH) IY=IY+1
! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),        &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,              &
     &              LCAL360)
              TIME2=REAL(SEC)/3600+REAL(DAY*24)
           ENDIF  ! end LAMIP2X test

! Do not interpolate in time if data time exactly matches model time

              IF(TIME == TIME1) THEN
                LINTERPOLATE=.FALSE.
              END IF

            ENDIF  ! regular/non-regular


!L Adjust YMD time back to start of updating interval

            I_YEAR=I_YEAR1
            I_MONTH=I_MONTH1
            I_DAY=I_DAY1
            I_HOUR=I_HOUR1


          ENDIF  ! non-periodic/periodic

        IF (LINTERPOLATE) THEN
        WRITE(6,*)' REPLANCA - time interpolation for field ',field
        WRITE(6,*)' time,time1,time2 ',time,time1,time2
        WRITE(6,*)' hours,int,period ',hours,interval,period
        END IF

        END IF ! singletime/non-singletime

!L 2.3   Check STASH Code

        I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*STEP

        I1=LOOKUP(ITEM_CODE,I2+LOOKUP_START(FILE)-1)

        LMISMATCH=.FALSE.
        WRITE(6,*)' Information used in checking ancillary data set:',  &
     &  ' position of lookup table in dataset:',I2
        WRITE(6,*)' Position of first lookup table referring to ',      &
     &  'data type ',NLOOKUP(FIELD)
        WRITE(6,*)' Interval between lookup tables referring to data ', &
     &  'type ', LOOKUP_STEP(FIELD),' Number of steps', STEP
        WRITE(6,*)' STASH code in dataset ',I1,                         &
     &  '  STASH code requested ',STASHANCIL(FIELD)
        WRITE(6,*)'''Start'' position of lookup tables for dataset ',   &
     &  'in overall lookup array ' ,LOOKUP_START(FILE)

        IF(I1 /= STASHANCIL(FIELD)) THEN
        WRITE(6,*)I1,STASHANCIL(FIELD),FIELD
          LMISMATCH=.TRUE.
        END IF

!L Error exit if checks fail

        IF(LMISMATCH) THEN
          ICODE=200+FIELD
         CMESSAGE='REPLANCA: PP HEADERS ON ANCILLARY FILE DO NOT MATCH'
         RETURN
        END IF

        IF(LINTERPOLATE.AND..NOT.SINGLE_TIME) THEN
!L Check time interpolation factors
          IF(TIME <  TIME1.OR.TIME >  TIME2) THEN
           WRITE(6,*)' Information used in interpolation/replacement:'
           WRITE(6,*)' Time of first data=', TIME1
           WRITE(6,*)' Validity Time for update=', TIME
           WRITE(6,*)' Time of second data=', TIME2

           ICODE=500+FIELD
           CMESSAGE='REPLANCA: TIME INTERPOLATION ERROR'
           RETURN
          END IF
        END IF

!L 3   Loop over levels of ancillary data for field I
!L Reset pointer for dataset


!L Includes loop over X and Y components of surface currents

         LICE_FRACTION=FIELD == 27
         LSNOW_DEPTH=FIELD == 9
         LICE_DEPTH=FIELD == 29

        DO 30 LEVEL=1,LEVELS(FIELD)

!L Do not go through loop for ice edge or snow edge

        IF((LICE_FRACTION.OR.LSNOW_DEPTH).AND.LEVEL == 2) THEN
          GOTO 30
        END IF

!L 3.1 Read data for single level of ancillary field.

        IF(.NOT.LICE_FRACTION) THEN
! AMIPII case ice depth field not read from ancillary file
         IF(.NOT.(LICE_DEPTH.and.LAMIP2X)) THEN
! DEPENDS ON: readflds
          CALL READFLDS(NFTIN,1,I2,LOOKUP(1,LOOKUP_START(FILE)),        &
     &                  LEN1_LOOKUP,ANCIL1,P_FIELD,FIXHD(1,FILE),       &
#include "argppx.h"
     &                  ICODE,CMESSAGE)

         ENDIF
          IF(ICODE /= 0)THEN
            ICODE=FIELD+100
            IOUNIT=NFTIN
           CMESSAGE='REPLANCA :I/O ERROR '
            RETURN
          END IF

        ELSE

!L If ice-fraction,read fractional time field as well
!L       UNLESS IT IS A SINGLE TIME FIELD
!L If snow-depth,read fractional time field as well only if time
!L interpolation required.

      IF(.NOT.SINGLE_TIME.and..NOT.LAMIP2X) THEN
         IF(LOOKUP(ITEM_CODE,I2+LOOKUP_START(FILE)) == 38) THEN
! DEPENDS ON: readflds
          CALL READFLDS(NFTIN,2,I2,LOOKUP(1,LOOKUP_START(FILE)),        &
     &                  LEN1_LOOKUP,ICE_EXTENT,P_FIELD,FIXHD(1,FILE),   &
#include "argppx.h"
     &                  ICODE,CMESSAGE)
          IF(ICODE /= 0)THEN
            ICODE=FIELD+100
            IOUNIT=NFTIN
           CMESSAGE='REPLANCA :I/O ERROR '
            RETURN
          END IF

         ELSE
           ICODE=FIELD+100
           IOUNIT=NFTIN
            CMESSAGE='REPLANCA :ICE CHANGE DATA MISSING'
            RETURN
         END IF
        ELSE    ! single time or LAMIP2X - ie no time change field
! DEPENDS ON: readflds
          CALL READFLDS(NFTIN,1,I2,LOOKUP(1,LOOKUP_START(FILE)),        &
     &                  LEN1_LOOKUP,ICE_EXTENT,P_FIELD,FIXHD(1,FILE),   &
#include "argppx.h"
     &                  ICODE,CMESSAGE)
          IF(ICODE /= 0)THEN
            ICODE=FIELD+100
            IOUNIT=NFTIN
            CMESSAGE='REPLANCA: I/O ERROR'
            RETURN
          ENDIF
        END IF
      ENDIF

        IF(LSNOW_DEPTH.AND.LINTERPOLATE) THEN
      IF(LOOKUP(ITEM_CODE,I2+LOOKUP_START(FILE)) == 27) THEN

! DEPENDS ON: readflds
           CALL READFLDS(NFTIN,1,I2+1,LOOKUP(1,LOOKUP_START(FILE)),     &
     &                   LEN1_LOOKUP,SNOW_CHANGE,P_FIELD,FIXHD(1,FILE), &
#include "argppx.h"
     &                   ICODE,CMESSAGE)
          IF(ICODE /= 0)THEN
             ICODE=FIELD+100
             IOUNIT=NFTIN
            CMESSAGE='REPLANCA :I/O ERROR '
             RETURN
           END IF

         ELSE
           ICODE=FIELD+100
           IOUNIT=NFTIN
           CMESSAGE='REPLANCA :SNOW CHANGE DATA MISSING'
           RETURN
         END IF
        END IF

!L If sea surface temperature or other ice fields, read ice fraction
!L and fractional time field if not already pressent and if required
!L by time interpolation.  Similar if SLAB ref SST or ice depth needed.

        IF(FIELD == 29.OR.(FIELD == 28.AND.LT_INT_C).OR.                &
     &     FIELD == 38)                                                 &
     &    THEN

         IF(.NOT.UPDATE(27)) THEN
          I3 = NLOOKUP(27) + LOOKUP_STEP(27)*STEP + LOOKUP_START(       &
     &       FILEANCIL(27))
          IF ( LOOKUP(ITEM_CODE,I3)  ==  38 ) THEN

! DEPENDS ON: readflds
            CALL READFLDS(FTNANCIL(FILEANCIL(27)),2,                    &
     &                    NLOOKUP(27)+LOOKUP_STEP(27)*STEP,             &
     &                    LOOKUP(1,LOOKUP_START(FILEANCIL(27))),        &
     &                    LEN1_LOOKUP,ICE_EXTENT,                       &
     &                    P_FIELD,FIXHD(1,FILEANCIL(27)),               &
#include "argppx.h"
     &                    ICODE,CMESSAGE)
          IF(ICODE /= 0)THEN
              ICODE=FIELD+100
              IOUNIT=NFTIN
             CMESSAGE='REPLANCA :I/O ERROR '
              RETURN
            END IF
          IF ( RLOOKUP(BMDI,I3-1)  /=  RMDI ) THEN
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: RMDI in lookup of ancillary file of ti&
     &mes of sea-ice chge not standard'
            RETURN
          ENDIF


          ELSE
            ICODE=FIELD+100
            IOUNIT=NFTIN
            CMESSAGE='REPLANCA :ICE FIELD DATA MISSING'
            RETURN
          END IF
         END IF
        END IF

!L 3.3 If time interpolation required, read second record

        IF(LINTERPOLATE) THEN

          I1=I2+ LOOKUP_STEP(FIELD)
          IF(I1 <= FIXHD(152,FILE)) THEN

! AMIP II and ice depth don't read in ice depth field
          IF (.NOT.(LAMIP2X.and.LICE_DEPTH)) THEN

! DEPENDS ON: readflds
            CALL READFLDS(NFTIN,1,I1,LOOKUP(1,LOOKUP_START(FILE)),      &
     &                    LEN1_LOOKUP,ANCIL2,P_FIELD,FIXHD(1,FILE),     &
#include "argppx.h"
     &                    ICODE,CMESSAGE)
          ENDIF
          IF(ICODE /= 0)THEN
              ICODE=FIELD+300
              IOUNIT=NFTIN
              CMESSAGE='REPLANCA :I/O ERROR '
              RETURN
            END IF

          ELSE !end of data on file

!L  If end of data has been reached go back to the start.If data is
!L  periodic.
!L  Otherwise cancel time interpolation

            IF(PERIODIC) THEN

              I1 = NLOOKUP(FIELD) + LEVEL - 1

! DEPENDS ON: readflds
              CALL READFLDS(NFTIN,1,I1,LOOKUP(1,LOOKUP_START(FILE)),    &
     &                      LEN1_LOOKUP,ANCIL2,P_FIELD,FIXHD(1,FILE),   &
#include "argppx.h"
     &                      ICODE,CMESSAGE)
          IF(ICODE /= 0)THEN
                ICODE=FIELD+300
                IOUNIT=NFTIN
               CMESSAGE='REPLANCA :I/O ERROR '
                RETURN
              END IF
            ELSE
              LINTERPOLATE=.FALSE.
            END IF
          END IF! End of position on file test

          ICODE=0
        END IF ! End LINTERPOLATE

!L 3.4 Perform time interpolation

        IF(LINTERPOLATE) THEN

          ZERO=0.0

!L Select appropriate time interpolation for each field
!  Snowdepth: set equal to zero if no snow cover

          IF(LSNOW_DEPTH) THEN
            DO I=1,P_FIELD
              PRES_VALUE(I)=ZERO
            END DO

! For the call to T_INT_C, need to know BMDI is OK for SNOW_CHANGE
!  which was read in from position I2+1.
          IF ( RLOOKUP(BMDI,LOOKUP_START(FILE)+I2)  /=  RMDI ) THEN
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: RMDI in lookup of ancillary file of ti&
     &mes of snow change non-standard '
            RETURN
          ENDIF

! DEPENDS ON: t_int_c
            CALL T_INT_C (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,         &
     &           TIME,P_FIELD,SNOW_CHANGE,ANCIL1,PRES_VALUE)

! Ice fraction: ice depth set equal to zero if no ice

          ELSE IF(FIELD == 27.OR.FIELD == 29.OR.FIELD == 38) THEN
            IF(FIELD == 27) THEN
! For the call to T_INT_C, need to know BMDI is OK for ICE_EXTENT(1,2)
!  which was read in from position I1+1
          IF(.NOT.LAMIP2X) THEN
          IF ( RLOOKUP(BMDI,LOOKUP_START(FILE)+I1)  /=  RMDI ) THEN
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: RMDI in lookup of ancillary file of ti&
     &mes of sea-ice chge non-standard'
            RETURN
          ENDIF
          ENDIF

             IF (LAMIP2X) THEN
! linear uncontrolled time interpolation
! DEPENDS ON: t_int
              CALL T_INT (ICE_EXTENT,TIME1,ANCIL2,TIME2,ANCIL_DATA,     &
     &             TIME,P_FIELD)

! For AMIP II strictly ice concentrations should range between
! 0.0 and 1.0 but because of assumptions on T* made by the boundary
! layer and radiation schemes ice concentrations are restricted to
! 0.3 to 1.0. This will allow SSTs in areas of less than 30% ice to
! be used rather than TFS=-1.8C.

              DO I=1,P_FIELD
                IF (ANCIL_DATA(I) <  0.3) ANCIL_DATA(I)=0.0
                IF (ANCIL_DATA(I) >  1.0) ANCIL_DATA(I)=1.0
              ENDDO

             ELSE       ! non AMIPII option
              DO I=1,P_FIELD
                PRES_VALUE(I)=0
              END DO

! DEPENDS ON: t_int_c
              CALL T_INT_C (ICE_EXTENT,TIME1,ANCIL2,TIME2,ANCIL_DATA,   &
     &             TIME,P_FIELD,ICE_EXTENT(1,2),ICE_EXTENT,PRES_VALUE)

             ENDIF     ! end AMIPII test

            ELSE IF (FIELD == 29.OR.FIELD == 38) THEN

              DO I=1,P_FIELD
                PRES_VALUE(I)=0
              END DO

! DEPENDS ON: t_int_c
              CALL T_INT_C (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,       &
     &             TIME,P_FIELD,ICE_EXTENT(1,2),ICE_EXTENT,PRES_VALUE)


            END IF


! Sea surface temperature, set equal to TFS if ice present

          ELSE IF (FIELD == 28.AND.LT_INT_C) THEN
           IF (LAMIP2X) THEN

! DEPENDS ON: t_int
            CALL T_INT (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,           &
     &              TIME,P_FIELD)
! remove any T below TFS
            DO I=1,P_FIELD
              IF (ANCIL_DATA(i) <  TFS)  ANCIL_DATA(I)=TFS
            ENDDO

           ELSE     ! non AMIPII option

            IF(.NOT.LTLEADS)THEN
            DO I=1,P_FIELD
                PRES_VALUE(I)=TFS

! Set no_ice_extent indicator for controlled SST interpolation
                IF(ICE_EXTENT(I,1) == 0) THEN
                  NO_ICE_EXTENT(I)=1.0
                ELSE
                  NO_ICE_EXTENT(I)=0.0
                ENDIF
            END DO

! DEPENDS ON: t_int_c
            CALL T_INT_C (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,         &
     &           TIME,P_FIELD,ICE_EXTENT(1,2),NO_ICE_EXTENT,PRES_VALUE)
            ELSE
! DEPENDS ON: t_int
            CALL T_INT (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,           &
     &              TIME,P_FIELD)
            ENDIF

           ENDIF   ! end AMIPII test
! Otherwise linear interpolation in time, unless missing data indicator
! present at either time.

          ELSE

! Time interpolation checks the data against the standard missing data
!   indicator - check that the field is labelled as using the same one.
!  (It is to have the right I1 here that I3 is used above.)
          IF ( RLOOKUP(BMDI,LOOKUP_START(FILE)+I1-1)  /=  RMDI .OR.     &
     &         RLOOKUP(BMDI,LOOKUP_START(FILE)+I2-1)  /=  RMDI ) THEN
            WRITE (6, *) 'LOOKUPS:',                                    &
     &         RLOOKUP(BMDI,LOOKUP_START(FILE)+I1-1),                   &
     &         RLOOKUP(BMDI,LOOKUP_START(FILE)+I2-1)
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: Missing data indicator in lookup of an&
     &cillary file is non-standard    '
            RETURN
          ENDIF

          LEN=P_FIELD
!L  Ozone, test for zonal mean or full field
          IF(FIELD == 7) THEN
            IF(LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) THEN
              LEN=P_ROWS
            END IF
!   Cariolle ozone, test for zonal mean or full field. 
!   Currently same test as for conventional ozone.
          ELSE IF (FIELD >= 178.AND.FIELD <= 185) THEN
            IF(LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) THEN
              LEN=P_ROWS
            END IF
          END IF

! DEPENDS ON: t_int
            CALL T_INT(ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,            &
     &                 TIME,LEN)

          END IF ! End Lsnow_depth

! If no interpolation, copy data into final array

        ELSE ! no interpolation
         IF(LICE_FRACTION) THEN
          IF (LAMIP2X) THEN
          DO I=1,P_FIELD

          ANCIL_DATA(I)=ICE_EXTENT(I,1)

! For AMIP II strictly ice concentrations should range between
! 0.0 and 1.0 but because of assumptions on T* made by the boundary
! layer and radiation schemes ice concentrations are restricted to
! 0.3 to 1.0. This will allow SSTs in areas of less than 30% ice to
! be used rather than TFS=-1.8C.

             IF (ANCIL_DATA(I) <  0.3) ANCIL_DATA(I)=0.0
             IF (ANCIL_DATA(I) >  1.0) ANCIL_DATA(I)=1.0

          ENDDO
          ELSE           ! non AMIP II option
            DO I=1,P_FIELD
             ANCIL_DATA(I)=ICE_EXTENT(I,1)
            ENDDO
          ENDIF           ! end of AMIPII test
         ELSE IF (LAMIP2X.AND.FIELD == 28) THEN
          DO I=1,P_FIELD
            ANCIL_DATA(I)=ANCIL1(I)
            IF (ANCIL_DATA(I) <  TFS) ANCIL_DATA(I)=TFS
          ENDDO
         ELSE
          DO I=1,P_FIELD
            ANCIL_DATA(I)=ANCIL1(I)
          END DO
         ENDIF
        END IF !End interpolate/no interpolate

!L 3.5 Updating action for each field at each level
!L     Fields replaced except that Sea Surface Temperature may be
!L     incremented. Take apropriate action for each field.

        IF(FIELD <= 2.OR.FIELD == 7.OR.FIELD == 39.OR.FIELD == 40       &
     &  .OR.FIELD == 41.OR.FIELD == 42.OR.FIELD == 43                   &
     &  .OR.FIELD == 44.OR.FIELD == 45                                  &
                                          ! multi-level murk
     &  .OR.(FIELD >= 48 .AND. FIELD <= 67 .AND. L_UKCA)                &
                                          ! single-level user ancillaries
     &  .OR.(FIELD >= 68.AND.FIELD <= 70)                               &
                                          !NH3,soot aerosol emissions
     &  .OR.(FIELD >= 72.AND.FIELD <= 77)                               &
                                          !Sulphur cycle
     &  .OR.FIELD == 78                                                 &
                                          !CO2 EMISSIONS
     &  .OR.FIELD == 82                                                 &
                                          !HADCM2 sulphate aerosol
     &  .OR.(FIELD >= 90.AND.FIELD <= 109)                              &
                                           !multi-level user ancillaries
     &  .OR.(FIELD >= 112.AND.FIELD <= 120)                             &
                                            !mineral dust fields
     &  .OR.(FIELD >= 121.AND.FIELD <= 122)                             &
                                            !Biomass emissions
     &  .OR.FIELD == 123                                                &
                                           !Seawater DMS concentration
     &  .OR.(FIELD >= 157.AND.FIELD <= 177)                             &
                                           !Aerosol climatologies
     &  .OR.(FIELD >=178 .AND.FIELD<=185)                               &
                                           !Cariolle ozone ancillaries
     &  .OR.(FIELD >= 186.AND.FIELD <= 187)                             &
                                           !OCFF emissions
     &  )THEN

!L 3.5.0 Updates at all points

          LEN=P_FIELD
!L  Ozone, test for zonal mean or full field
          IF(FIELD == 7) THEN
            IF(LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) THEN
              LEN=P_ROWS
            END IF

!   Cariolle ozone, test for zonal mean or full field. 
!   Currently same test as for conventional ozone.
          ELSE IF (FIELD >= 178.AND.FIELD <= 185) THEN
            IF(LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) THEN
              LEN=P_ROWS
            END IF
          END IF

          DO I=1,LEN
            D1(D1_ANCILADD(FIELD)+I-1+(LEVEL-1)*LEN)=ANCIL_DATA(I)
          END DO

!L 3.5.1 Updates over all land points

        ELSEIF((FIELD > 2.AND.FIELD < 7)                                &
     &   .OR.(FIELD > 7.AND.FIELD < 27)                                 &
     &   .OR.(FIELD == 32).OR.(FIELD >= 34.AND.FIELD <= 36)             &
     &   .OR.(FIELD >= 48 .AND. FIELD <= 67 .AND. .NOT. L_UKCA)         &
                                      ! single level user ancillaries
     &   .OR.(FIELD >= 46.AND.FIELD <= 47)                              &
                                                 !Orographic roughness
     &   .OR.(FIELD >= 155.AND.FIELD <= 156)                            &
                                      ! Orographic X & Y gradients
     &   .OR.(FIELD >= 79.AND.FIELD <= 81)                              &
                                                 !MOSES-I
     &   .OR.(FIELD >= 83.AND.FIELD <= 89)) THEN !MOSES-II



!L If not reconfiguration, set snowdepth values at all land points
!L Reset TSTAR to TM if snow cover present

          IF(LSNOW_DEPTH) THEN
            DO I=1,P_FIELD
              IF(LAND(I)) THEN
                D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
                IF(TSTAR_LAND(I) >  TM.AND.ANCIL_DATA(I) >  0.0) THEN
                  TSTAR_LAND(I)=TM
                END IF
              END IF
            END DO

!L Set all other fields , which are stored at land points only

          ELSE
!L If field is a single time soil moisture, only update if
!L model time matches data time and then deactivate to prevent
!L any further updating

            IF(FIELD == 36.AND.FIXHD(10,FILE) == 0) THEN
              IF(LOOKUP(LBYR,LOOKUP_START(FILE)) == I_YEAR.AND.         &
     &           LOOKUP(LBMON,LOOKUP_START(FILE)) == I_MONTH.AND.       &
     &           LOOKUP(LBDAT,LOOKUP_START(FILE)) == I_DAY.AND.         &
     &           LOOKUP(LBHR,LOOKUP_START(FILE)) == I_HOUR) THEN
                WRITE(6,*) 'Updating soil moisture at ',                &
     &                       I_YEAR,I_MONTH,I_DAY,I_HOUR,               &
     &          ' for level ',RLOOKUP(BLEV,LOOKUP_START(FILE)+LEVEL-1)
                DZ_SOIL(LEVEL)=RLOOKUP(BLEV,LOOKUP_START(FILE)+LEVEL-1)
! DEPENDS ON: to_land_points
                CALL TO_LAND_POINTS(ANCIL_DATA,D1(D1_ANCILADD(FIELD)+   &
     &                             (LEVEL-1)*LAND_FIELD),LAND,P_FIELD,I)
! Switch off to prevent further attempts to update
                FIELDCODE(1,FIELD)=0
                STEPS(FIELD)=0
! Set flag to indicate that soil moisture has been updated
                SMC_UPDATED=.TRUE.
              ELSE
                WRITE(6,*) 'Update of soil moisture skipped'
              END IF

            ELSE
! other fields
! DEPENDS ON: to_land_points
              CALL TO_LAND_POINTS(ANCIL_DATA,D1(D1_ANCILADD(FIELD)+     &
     &                           (LEVEL-1)*LAND_FIELD),LAND,P_FIELD,I)
            END IF

          END IF



!L 3.5.2 Ice fraction


        ELSE IF(FIELD == 27) THEN
          DO I=1,P_FIELD
            ICE_FRACTION(I)=0.
            IF (SEA(I)) THEN
              ICE_FRACTION(I)=ANCIL_DATA(I)
            END IF
          END DO

!L Reduce TSTAR to TFS where ice fraction greater than zero
! Required at present because radiation and boundary layer codes
! assume T* is TFS and ignore any value set in TSTAR.

          IF(.NOT.LTLEADS)THEN
            DO I=1,P_FIELD
              IF(ICE_FRACTION(I) >  0.0) THEN
                TSTAR_SSI(I)=AMIN1(TSTAR_SSI(I),TFS)
              ENDIF
            END DO
          ENDIF

!L 3.5.3 Sea surface temperatures for atmosphere, allow fields to be
!L       incremented rather than replaced

        ELSE IF (FIELD == 28) THEN


          DO I=1,P_FIELD
            IF(L_CTILE.OR.ICE_FRACTION(I) == 0.0)THEN
              IF (SEA(I)) THEN
                IF (L_SSTANOM) THEN
                TSTAR_SEA(I)=ANCIL_DATA(I)+TSTAR_ANOM(I)
                ELSE
                  TSTAR_SEA(I)=ANCIL_DATA(I)
                END IF
                IF(ICE_FRACTION(I) == 0.0)TSTAR_SSI(I)=TSTAR_SEA(I)
              END IF
            END IF
          END DO

!L 3.5.3.1 Reference SSTs for SLAB model

        ELSE IF (FIELD == 37) THEN

          DO I=1,P_FIELD
            IF (SEA(I)) THEN
              D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
            ELSE
              D1(D1_ANCILADD(FIELD)+I-1)=RMDI

            ENDIF
          END DO

!L 3.5.4 Sea ice thickness/Reference seaice thickness for SLAB
!L       Update over all sea points (all sea ice points are the only
!L       ones strictly required, but this cannot be determined easily)

        ELSE IF (FIELD == 29.OR.FIELD == 38) THEN

          DO I=1,P_FIELD
            IF(SEA(I)) THEN
              D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
            END IF
          END DO

!L 3.5.5 Surface currents

        ELSE IF (FIELD == 30) THEN
          DO I=1,U_FIELD
            D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
          END DO

        ELSE IF (FIELD == 31) THEN
          DO I=1,V_FIELD
            D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
          END DO
!L 3.5.6 Heat convergence (slab model)
!L       Update over all non-land points

        ELSE IF (FIELD == 33) THEN

          DO I=1,P_FIELD
            IF(SEA(I)) THEN
              D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
            ELSE
              D1(D1_ANCILADD(FIELD)+I-1) = 0.0
            END IF
          END DO

        ELSE IF (FIELD >= 124.AND.FIELD <= 126) THEN   ! Riv Routing
          ICODE=750
          CMESSAGE='REPLANCA: ERROR Trying to use Riv Route Ancils'
!         There is no code yet to support in model updateing of the
!         River Routing fields.
          RETURN

! Tracer Fluxes - kdcorbin, 05/10
        ELSE IF (FIELD .ge. 188 .and. FIELD .lt. 208) Then
            DO I=1,P_FIELD
              D1(D1_ANCILADD(FIELD)+I-1) = ANCIL_DATA(I)
            ENDDO

        ELSE

        WRITE(6,*)' REPLANCA: ERROR - FIELD ',FIELD,                    &
     &  ' omitted from update block'

        END IF !End tests on FIELD numbers

!L End loop over levels

      I2=I2+1

 30   CONTINUE

!L End loop over ancillary fields (atmosphere)
       ENDIF ! LAMIP2X and ice depth test

      END IF    ! End UPDATE(field) test     level 1 IF


      END DO


      IF(L_CTILE)THEN
        DO I=1,P_FIELD
          IF(SEA(I).AND.ICE_FRACTION(I) >  0.0)THEN
            IF(LTLEADS.OR.LAMIP2X)THEN

              TSTAR_SSI(I)=ICE_FRACTION(I)*TSTAR_SICE(I)                &
     &          +(1.-ICE_FRACTION(I))*TSTAR_SEA(I)

            ELSE

              TSTAR_SEA(I)=TFS
              TSTAR_SICE(I)=(TSTAR_SSI(I)                               &
     &          -(1.-ICE_FRACTION(I))*TSTAR_SEA(I))/ICE_FRACTION(I)

            ENDIF
          ENDIF
!
          TSTAR(I)=FLAND_G(I)*TSTAR_LAND(I)                             &
     &      +(1.-FLAND_G(I))*TSTAR_SSI(I)
        ENDDO
      ELSE
        DO I=1,P_FIELD
          IF(LAND(I))THEN
            TSTAR(I)=TSTAR_LAND(I)
          ELSE
            TSTAR(I)=TSTAR_SSI(I)
          ENDIF
        ENDDO
      ENDIF

!     Set up surface temperatures:
      IF(L_CTILE)THEN
        DO I=1,P_FIELD
          TSTAR_LAND_CTILE(I)=TSTAR_LAND(I)
          TSTAR_SEA_CTILE(I)=TSTAR_SEA(I)
          ! The use of TSTAR_SICE appears to cause problems in 
          ! some configurations (e.g. seasonal). Possibly because
          ! of the use of inconsistent ancillary fields.
          ! Hence this is commented out but retained for reference.
          ! [See also equivalent change in replanca-rcf_replanca.F90]
          ! TSTAR_SICE_CTILE(I)=TSTAR_SICE(I)
        ENDDO
      ENDIF

!
900   RETURN
      END SUBROUTINE REPLANCA

#endif
