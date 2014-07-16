#if defined(C82_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine INANCILA
!LL
!LL Purpose : Takes as input,the code defining the frequency of update
!LL           of ancillary fields as set by the user interface.
!LL           Converts them into a list of numbers of timesteps after
!LL           which each field must be updated, and calculates the
!LL           frequency with which this list must be interrogated.
!LL           Where the update interval is in months or years,
!LL           the check will be carried out each day. The physical
!LL           files required are also determined by input code,
!LL           and the headers and lookup tables are read into
!LL           the arguments FIXHD,INTHD,LOOKUP which are in
!LL           COMMON/ANCILHDA/ of calling routine INANCCTL.
!LL           Indexes for each possible ancillary field are set up in
!LL           COMMON/IXANCILA/
!LL
!LL Level 2 Control routine for CRAY YMP
!LL
!LL CW, DR      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.1  22/02/93  Changes to add 2 SLAB fields (STASH items 178,179)
!LL                  - to be updated from existing atmosphere files.
!LL   3.3  22/11/93  Add aerosol ancillary fields.  R T H Barnes.
!LL   3.3  21/12/93  Fix put in to prevent array 'out of bounds'
!LL                  problem in section 1.6. Problem to be investigated
!LL                  for 3.4 D. Robinson.
!LL   3.4  16/06/94  DEF CAL360 replaced by LOGICAL LCAL360
!LL                                                   S.J.Swarbrick
!LL  3.4  05/09/94  Add murk and user ancillary fields.  RTHBarnes.
!LL   3.4  22/06/94  Array 'out of bounds' problem solved. D. Robinson
!LL  3.4   11/10/94   Part of modset which sorts out some handling
!LL                   of unset data by recon_dump.
!LL                   Necessary to port model to a T3D.
!LL                   Author D.M. Goddard
!LL   3.5  24/03/95    Changed OPEN to FILE_OPEN  P.Burton
!LL  3.5  24/07/95  Check fields for updating have valid address. RTHB
!    4.0  01/09/95  Add diagnostic information to output about
!                   ozone ancillary fields and test correct ozone
!                   data provided.  D. Goddard & D. Robinson
!LL  4.0  10/10/95  Set LOOKUP(45) in ancillary files. D. Robinson.
!LL
!LL  4.0  29/09/95  Need extra rewind of namelist file. RTHBarnes.
!LL  4.0  05/08/95  Temporary solution to get round problem of
!LL                 no. of soil moisture levels being hard-wired
!LL                 to no. of deep soil temperature levels
!LL                 This causes a problem with introduction of
!LL                 Penman-Monteith BL code at 4.0 - use if test
!LL                 on number of deep soil temperature
!LL                 levels which is set to 4 for Penman-Monteith code
!LL                 (set to 3 for all other BL options)
!LL                 Permanent solution suggested for 4.1
!LL                 search on C**** for comments
!LL                 J Smith
!LL  4.0  06/01/96  SI array received for two internal models (atmos
!LL                 and slab) in argument list. Hardwire processing of
!LL                 slab ancillary field (item code 177) to use
!LL                 SI_SLAB. D. Robinson
!    4.1  03/05/96  Use READHEAD to read in ancillary file headers.
!                   D. Robinson
!    4.1  18/06/96  Changes to cope with changes in STASH addressing
!                   Author D.M. Goddard.
!LL  4.1  22/05/96  Call new CANC* comdecks. Use new arrays in
!LL                 CANCFLDA. Cater for new sulphur ancillary files.
!LL                 Remove hardwired fix for slab ancillary fields
!LL                 introduced at 4.0 D. Robinson.
!LL  4.4  28/07/97  Add LAMIPII to namelist for special updating of
!LL                 ice in AMIP II runs. R A Stratton
!LL  4.4  16/09/97  Set number of headers for multi-pseudo-level
!LL                 ancillary fields for surface and vegetation types.
!LL                                              Richard Betts
!LL  4.4  09/09/97  New namelist UPANCA for updating information.
!LL                 D. Robinson.
!LL  4.4  10/09/97  Check calendar indicator in Anc File. D Robinson.
!    4.5  22/10/98  Set LEVELS array for new user hulti-layer
!                   ancillary fields
!                   Author D.M Goddard
!LL  4.5  19/01/98  Remove SOIL_VARS and VEG_VARS. D. Robinson.
!LL  4.5  05/05/98  Improve error message for missing files. R. Rawlins
!LL  5.1  20/07/00  Remove potential OOB reference from STASHANCIL
!LL                                                     D.Robinson
!LL  5.2  09/03/01  Initialise FIELDCODE to zero. D.Robinson
!    5.3  23/01/02  Update vertical levels checking for ozone files.
!                   Dave Robinson
!LL  5.3  19/06/01  Added code for tropopause-based ozone. D.Tan
!    5.4  27/08/02  Update vertical levels checking for murk and
!                   mulit_level ancillary files. D Robinson.
!    6.1  18/08/04  Fix possible out of bounds error. S.Wilson
!    6.1  27/07/04  Correct dimensioning of LEVDEPC. D Robinson
!    6.2  22/08/05  Remove RECON def. P.Selwood
!LL
!LL System components covered : C710
!LL
!LL System task : C7
!LL
!LL Documentation : Unified Model Documentation Paper No C7
!LL                 Version No 4  dated 15/06/90
!LLEND
      SUBROUTINE INANCILA(LEN_FIXHD,LEN_INTHD,LEN_REALHD,               &
                                                           !Intent (In)
     &                    LEN1_LEVDEPC,LEN2_LEVDEPC,                    &
     &                    FIXHD,INTHD,REALHD,LOOKUP,                    &
     &                    A_FIXHD,A_REALHD,A_LEVDEPC,                   &
     &                    NDATASETS,NLOOKUPS,FTNANCIL,                  &
     &                    LOOKUP_START,LEN1_LOOKUP,ROW_LENGTH,          &
     &                    P_ROWS,U_ROWS,P_LEVELS,                       &
     &                    TR_LEVELS,ST_LEVELS,SM_LEVELS,                &
     &                    OZONE_LEVELS,tpps_ozone_levels,TITLE,         &
!kdcorbin, 08/10 - added variables to call
     &                    SI,NITEM,NSECT,N_MODEL,                       &
     &                    SI_ATMOS,SI_SLAB,SILEN,                       &
     &                    ANCILLARY_STEPS,STEPS_PER_HR,                 &
#include "argppx.h"
     &                    ICODE,CMESSAGE,LCAL360)         ! Intent (Out)

      USE CSENARIO_MOD
      IMPLICIT NONE

      LOGICAL LCAL360  ! Logical switch for 360-day calendar

      INTEGER                                                           &
     &        LEN_FIXHD,                                                &
                               ! Length of header blocks in ancillary
!                              ! data sets
     &        LEN_INTHD,                                                &
                               !
     &        LEN_REALHD,                                               &
                               !
     &        LEN1_LEVDEPC,                                             &
                               ! Dimension of LEVDEPC in model
     &      LEN2_LEVDEPC                                                &
     &     ,ANCILLARY_STEPS,                                            &
     &        STEPS_PER_HR


      INTEGER                                                           &
     &        NDATASETS,                                                &
                               ! No of physical files
     &        NLOOKUPS,                                                 &
                               ! No of lookups required(set by User I.)
     &                IOUNIT,                                           &
     &        FTNANCIL(NDATASETS),                                      &
                                   ! Fortran nos of physical files
     &        LOOKUP_START(NDATASETS),                                  &
                                      !start of each individual lookup
!                                     !in overall LOOKUP array
     &        LEN1_LOOKUP,                                              &
                               ! Length of PP header
     &        ROW_LENGTH,                                               &
                               ! Atmosphere model dimensions
     &        P_ROWS,                                                   &
                               ! No. of rows for pressure-type variables
     &        U_ROWS,                                                   &
                               ! No. of rows for wind-type variables
     &        P_LEVELS,                                                 &
                               ! No. of pressure levels
     &        TR_LEVELS,                                                &
                               ! No. of tracer levels
     &        FILE_LEVELS,                                              &
                               ! Number of levels of data in files
!                              ! contining multi-level data.
     &        ST_LEVELS,                                                &
                               ! No. of soil temperature levels
     &        SM_LEVELS,                                                &
                               ! No. of soil moisture levels
     &        OZONE_LEVELS,                                             &
                                 ! No. of ozone levels
     &        tpps_ozone_levels                                         &
                                 ! No of ozone levs in TppsOzon dataset

!      For atmos only runs SI_SLAB is a copy of SI_ATMOS
!      SI_SLAB is only used in SLAB runs.

     &       ,SILEN                                                     &
                                ! Length for SI_ATMOS/SLAB arrays
     &       ,SI_ATMOS(SILEN)                                           &
                                ! ) STASHin addresses of atmos and
     &       ,SI_SLAB(SILEN)    ! ) slab ancillary fields.

     !Added SI array - kdcorbin, 05/10
      INTEGER :: NITEM,NSECT,N_MODEL
      INTEGER :: SI(NITEM,0:NSECT,N_MODEL)

      !Added filename for ancillary data - kdcorbin, 05/10
      INTEGER, PARAMETER :: Max_Filename_Len=120
      CHARACTER*Max_Filename_Len :: AncFileName
      INTEGER :: len_anc_filename

      CHARACTER*80 TITLE(NDATASETS) ! Titles of each dataset

      Logical :: l_vert_mismatch    ! T : Vertical levels mismatch

      INTEGER                                                           &
     &        FIXHD(LEN_FIXHD,NDATASETS),                               &
                                         ! Overall Fixed header array
     &        A_FIXHD(LEN_FIXHD),                                       &
                                  ! Fixed header for Dump
     &        INTHD(LEN_INTHD,NDATASETS),                               &
                                         ! Overall Integer header array
     &        LOOKUP(LEN1_LOOKUP,NLOOKUPS),                             &
                                           !Overall Lookup array
     &        ICODE            ! Return code =0 Normal Exit
!                              !             >0 Error

      REAL                                                              &
     &      REALHD(LEN_REALHD,NDATASETS),                               &
                                         !
     &      A_REALHD(LEN_REALHD),                                       &
                                 !
     &      A_LEVDEPC(LEN1_LEVDEPC,LEN2_LEVDEPC),                       &
     &      LEVDEPC( (P_LEVELS+1)*4 ) ! Space to hold level dependent
                                      ! constants from Anc File

      CHARACTER*100                                                     &
     &        CMESSAGE         ! Out error message if I>0
      Character (Len=*), Parameter :: RoutineName = 'INANCILA'

      INTEGER :: isec, iind  !kdcorbin, 05/10

! Comdecks:----------------------------------------------------------
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "model.h"
#include "clookadd.h"
#include "cancila.h"
#include "nstypes.h"
#include "c_mdi.h"
#include "cenvir.h"
#include "parvars.h"
#include "cprintst.h"


! Comdecks for ancillary files/fields.
#include "cancflda.h"

!L External subroutines called:

      EXTERNAL                                                          &
     &        FILE_OPEN,                                                &
     &        READ_FLH, READHEAD, SETPOS

!L Namelist input

      NAMELIST/ANCILCTA/L_SSTANOM,LAMIPII

!     UPANCA Namelist
      INTEGER                                                           &
     &   ANC_REF_NO                                                     &
                          ! Ancil Ref. No : See comdeck CANCFLDA
     &  ,PERIOD                                                         &
                          ! Period of Updating Interval (Y/M/D/H)
     &  ,INTERVAL         ! Updating Interval

     !Added finput flag and filename to namelist info - kdcorbin, 04/10
      INTEGER :: FINPUT
      CHARACTER*Max_Filename_Len :: FNAME


      NAMELIST /UPANCA/ ANC_REF_NO,PERIOD,INTERVAL,FINPUT,FNAME

! Local Variables

      INTEGER                                                           &
     &        I,                                                        &
                               !
     &        ITEM,                                                     &
                               !
     &        J,                                                        &
                               !
     &        J1,                                                       &
                               !
     &        K,                                                        &
                               !
     &        LEN_IO,                                                   &
                               !
     &        LOOKUPS,                                                  &
                               !
     &        NFTIN,                                                    &
                               ! Current FTN number for ancillary data
     &        START_BLOCK                                               &
                               !
     &       ,STASH_CODE                                                &
                               ! Stash item code
     &       ,NREC_A,NREC_S                                             &
                               ! No of atmos & slab records
     &       ,STASH_ADDR                                                &
                               ! Stash address
     &       ,DUMMY                                                     &
                               !
     &       ,N_ANC_UPD        ! No of ancillaries to be updated
      DATA DUMMY /1/

      CHARACTER*8 CPERIOD      ! PERIOD in characters.
      LOGICAL                                                           &
     &        LFILE            !

!     SoilDepths : position of soil depths in level dependent constants
      Integer, Parameter :: SoilDepths = 4

      REAL P1,P2
      LOGICAL LNER
      LNER(P1,P2) = ((ABS(P1-P2))  >   (1.E-6*ABS(P1+P2)))

!L Internal Structure

      ICODE=0
      CMESSAGE=' '
      IOUNIT=0

!
!L  1.  Initialisation for atmosphere model

      DO I=1,NANCIL_FIELDS
        FILEANCIL(I) =ANCIL_FILE_NO(I)
        STASHANCIL(I)=ITEM_CODES_ANCIL(I)
      ENDDO
      FIELDCODE(:,:) = 0


! Set default values

      L_SSTANOM=.FALSE.
      LAMIPII=.FALSE.

!L  Read in control information from namelist

        REWIND 5
      READ(5,ANCILCTA)

!     Initialise FIELDCODE from Namelist UPANCA
      N_ANC_UPD = 0
      FINPUT = 0
      FNAME = ''
      DO I=1,NANCIL_FIELDS
        READ (5,UPANCA,ERR=101,END=101)
        FIELDCODE(1,ANC_REF_NO) = PERIOD
        FIELDCODE(2,ANC_REF_NO) = INTERVAL
        ANC_FILE_FINPUT(ANC_REF_NO) = FINPUT  !kdcorbin, 04/10
        ANC_FILE_FNAME(ANC_REF_NO) = FNAME  !kdcorbin, 04/10
        N_ANC_UPD = N_ANC_UPD+1
      ENDDO

 101  CONTINUE
      WRITE (6,*) ' '
      WRITE (6,*) N_ANC_UPD,' Atmos & Slab Ancillaries to be updated.'
      DO I=1,NANCIL_FIELDS
        IF (FIELDCODE(1,I) >  0) THEN
        IF (FIELDCODE(1,I) == 1) CPERIOD=' Years'
        IF (FIELDCODE(1,I) == 2) CPERIOD=' Months'
        IF (FIELDCODE(1,I) == 3) CPERIOD=' Days'
        IF (FIELDCODE(1,I) == 4) CPERIOD=' Hours'
        WRITE (6,*) 'Anc Ref No ',I,' Stash code ',ITEM_CODES_ANCIL(I), &
     &  ' Interval ',FIELDCODE(2,I),CPERIOD
        ENDIF
      ENDDO
      WRITE (6,*) ' '

! Check that ancillary field has valid address (>1) before proceding
!  to try and update it.  If not, switch off updating via FIELDCODE.
      DO I=1,NANCIL_FIELDS
      IF (STASHANCIL(I)  >   0) THEN
        if (model_codes_ancil(i) == slab_im) then
          stash_addr = si_slab(stashancil(i))
         else

          ! Including alternate sections from 0 to add 
          !   tracer fluxes - kdcorbin, 05/10
          !stash_addr = si_atmos(stashancil(i))
          isec=stashancil(i)/1000
          iind=stashancil(i)-isec*1000
          stash_addr = SI(iind,isec,1)

        endif
      ELSE
        stash_addr=0
      ENDIF
        IF (stash_addr  <=  1) THEN
          IF (FIELDCODE(1,I) >  0) THEN
           WRITE(6,*)' INANCILA: update requested for item ',i,         &
     &     ' STASHcode ',stashancil(i),' but prognostic address not set'
            WRITE(6,*)' FIELDCODE values reset to zeroes'
            FIELDCODE(1,I) = 0
            FIELDCODE(2,I) = 0
          END IF
        END IF
      END DO

!L  1.1 Set number of steps after which each ancillary field is updated
!       Zero is used for fields not to be updated

      DO I=1,NANCIL_FIELDS
        STEPS(I)=0
        IF (FIELDCODE(1,I) == 4)THEN
          STEPS(I)=FIELDCODE(2,I)*STEPS_PER_HR
        END IF
        IF (FIELDCODE(1,I) == 3) THEN
          STEPS(I)=FIELDCODE(2,I)*24*STEPS_PER_HR
        END IF

      IF (LCAL360) THEN
        IF (FIELDCODE(1,I) == 2) THEN
          STEPS(I)=FIELDCODE(2,I)*30*24*STEPS_PER_HR
        END IF
        IF (FIELDCODE(1,I) == 1) THEN
          STEPS(I)=FIELDCODE(2,I)*360*24*STEPS_PER_HR
        END IF
      ELSE
! Gregorian calender:
! If update interval is months or years, test each day. Further testing
! done in REPLANCA.

        IF (FIELDCODE(1,I) == 1.OR.FIELDCODE(1,I) == 2)THEN
         STEPS(I)=24*STEPS_PER_HR
        END IF
      END IF

      END DO

!L  1.2 Set master number of steps ANCILLARY_STEPS at which
!L      individual switches are tested.

!   Find first active field

      DO I=1,NANCIL_FIELDS
        IF (STEPS(I) >  0) THEN
          ANCILLARY_STEPS=STEPS(I)
          GOTO 121
        END IF
      END DO

! No above fields found

      ANCILLARY_STEPS=0

      GOTO 900
121   ITEM=I

!L      Set ANCILLARY_STEPS to lowest common denominater of
!L      frequencies for active fields

      DO I=ITEM+1,NANCIL_FIELDS
        IF (STEPS(I) <  ANCILLARY_STEPS                                 &
     &      .AND. STEPS(I) >  0) THEN
          IF (MOD(ANCILLARY_STEPS,STEPS(I)) == 0) THEN
            ANCILLARY_STEPS=STEPS(I)
          ELSE
            J1=STEPS(I)-1
            DO J=J1,1,-1
              IF ((MOD(ANCILLARY_STEPS,J) == 0).AND.                    &
     &           (MOD(STEPS(I),J) == 0)) THEN
                 GOTO 124
              ENDIF
            END DO
124         ANCILLARY_STEPS = J
          END IF
        END IF
      END DO

!L 1.2.4 Sea surface temperature must be updated when sea ice is update

      IF (STEPS(27) >  0.AND.STEPS(28) <= 0) THEN
         STEPS(28)=1
      END IF


!L 1.3 Set number of headers for each ancillary field

      DO I=1,NANCIL_FIELDS
        LEVELS(I)=1
!   Multilayer hydrology
        IF(I == 36)LEVELS(I)=SM_LEVELS
!   Multilayer aerosols
        IF(I >= 41.AND.I <= 43) LEVELS(I)=TR_LEVELS
!   Multilayer murk concentration and source
        IF(I >= 44.AND.I <= 45) LEVELS(I)=P_LEVELS
!   Multilayer user ancillaries
        IF(I >= 90.AND.I <= 109) LEVELS(I)=P_LEVELS
!   Multi-level ancillaries for sulphur cycle
        IF (I == 72) LEVELS(I) = P_LEVELS
        IF (I == 73) LEVELS(I) = P_LEVELS
        IF (I == 74) LEVELS(I) = P_LEVELS
        IF (I == 75) LEVELS(I) = P_LEVELS
        IF (I == 76) LEVELS(I) = P_LEVELS
        IF (I == 82) LEVELS(I) = NSULPAT
        IF (I == 83) LEVELS(I) = NTYPE
        IF (I == 84) LEVELS(I) = NPFT
        IF (I == 85) LEVELS(I) = NPFT
!   Multi-level ancillaries aerosol climatology
        IF(I >= 157.AND.I <= 177) LEVELS(I)=P_LEVELS
        IF(I >= 178.AND.I <= 185) LEVELS(I)=P_LEVELS
      END DO

      LEVELS(7)=OZONE_LEVELS
!! consider do a check for l_use_tpps_ozone and if set then
!! set levels(110)=tpps_ozone_levels
      LEVELS(10)=ST_LEVELS


!L 1.4 Read headers

      LOOKUPS=0

      DO I=1,NDATASETS

!  Initialise LOOKUP_START (=0 implies file I not required)
        LOOKUP_START(I)=0

!L Check whether each physical file is needed

        LFILE=.FALSE.
        DO 141 J=1,NANCIL_FIELDS


          IF (FILEANCIL(J) == I.AND.STEPS(J) >  0) THEN

            LFILE=.TRUE.

            !kdcorbin, 05/10 - added check for ancillary file name
            FINPUT=ANC_FILE_FINPUT(J)
            FNAME=ANC_FILE_FNAME(J)

          END IF
141     CONTINUE

        IF(LFILE) THEN

      ! Open the File
      ! Changed from opening with environmental variable to filename
      ! kdcorbin, 05/10
 
      NFTIN=FTNANCIL(I)

      if (ANC_FILE_FINPUT(NFTIN) == 1) Then
         AncFileName=ANC_FILE_FNAME(NFTIN)
      else
         CALL Fort_Get_Env(FT_ENVIRON(NFTIN),LEN_FT_ENVIR(NFTIN),  &
             AncFileName,Max_Filename_Len,icode)
      endif

        Write(6,*) ''
        Write(6,*) 'Field: ',NFTIN,FT_ENVIRON(NFTIN)
        Write(6,*) 'Opening Anc File: ',trim(AncFileName)

        len_anc_filename=len_trim( AncFileName)
! DEPENDS ON: file_open
        CALL FILE_OPEN(NFTIN,AncFileName,                      &
     &                 len_anc_filename,0,1,ICODE)

      !Original File Open:
      !! DEPENDS ON: file_open
      !  CALL FILE_OPEN(NFTIN,FT_ENVIRON(NFTIN),                         &
      !               LEN_FT_ENVIR(NFTIN),0,0,ICODE)

        IF(ICODE /= 0)THEN
          CMESSAGE='INANCLA: Error opening file'
          write(6,*) 'INANCILA: Error opening file on unit ',NFTIN,     &
     &               ' accessed from env.var.: ',FT_ENVIRON(NFTIN)
          RETURN
        ENDIF
! DEPENDS ON: setpos
        CALL SETPOS(NFTIN,0,ICODE)

!       Read in fixed header to get array dimensions
! DEPENDS ON: read_flh
        CALL READ_FLH(NFTIN,FIXHD(1,I),LEN_FIXHD,ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
          WRITE (6,*) ' Error in reading fixed header for file ',I
          GO TO 9999   !  Return
        ENDIF

!       Check for negative dimensions
        IF (FIXHD(101,I) <= 0) FIXHD(101,I)=1
        IF (FIXHD(106,I) <= 0) FIXHD(106,I)=1
        IF (FIXHD(111,I) <= 0) FIXHD(111,I)=1
        IF (FIXHD(112,I) <= 0) FIXHD(112,I)=1
        IF (FIXHD(151,I) <= 0) FIXHD(151,I)=1
        IF (FIXHD(152,I) <= 0) FIXHD(152,I)=1
        IF (FIXHD(161,I) <= 0) FIXHD(161,I)=1

! Set start position of boundary fields for file
        LOOKUP_START(I)=LOOKUPS+1

        IF (LOOKUPS+FIXHD(152,I) >  NLOOKUPS) THEN
          WRITE (6,*) 'No room in LOOKUP table for Ancillary File ',I
          CMESSAGE='INANCILA: Insufficient space for LOOKUP headers'
          ICODE=14
          GO TO 9999   !  Return
        END IF

! DEPENDS ON: setpos
        CALL SETPOS(NFTIN,0,ICODE)
        IF (ICODE >  0) THEN
          WRITE (6,*) ' ERROR in SETPOS called from INANCA1A'
          WRITE (6,*) ' SETPOS attempted with Unit No ',NFTIN
          CMESSAGE = 'INANCA1A : ERROR in SETPOS'
          GO TO 9999    !   Return
        ENDIF

! DEPENDS ON: readhead
        CALL READHEAD(NFTIN,                                            &
     &                FIXHD(1,I),LEN_FIXHD,                             &
     &                INTHD(1,I),FIXHD(101,I),                          &
     &                REALHD(1,I),FIXHD(106,I),                         &
     &                LEVDEPC,FIXHD(111,I),FIXHD(112,I),                &
     &                DUMMY,DUMMY,DUMMY,                                &
     &                DUMMY,DUMMY,DUMMY,                                &
     &                DUMMY,DUMMY,DUMMY,                                &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                LOOKUP(1,LOOKUPS+1),FIXHD(151,I),FIXHD(152,I),    &
     &                FIXHD(161,I),                                     &
#include "argppx.h"
     &                START_BLOCK,ICODE,CMESSAGE)

        IF (ICODE >  0) THEN
           WRITE(6,*) 'ERROR in READHEAD for Ancillary File ',I
           WRITE(6,*) 'Unit Number ',NFTIN
           GO TO 9999   !   Return
        ENDIF

!     Check calendar indicator
        IF ((     LCAL360 .and. FIXHD(8,I) /= 2) .or.                   &
     &      (.not.LCAL360 .and. FIXHD(8,I) /= 1) ) THEN
          ICODE=100+I
          CMESSAGE='INANCILA : Wrong calendar set in Ancillary File'
          WRITE (6,*) ' ******** Error in INANCILA ********'
          WRITE (6,*) ' Wrong calendar setting in Ancillary File ',I
          IF (LCAL360) THEN
            WRITE (6,*) ' Model run is set up for 360 day calendar.'
            WRITE (6,*) ' Ancillary File is for 365 day calendar.'
          ELSE
            WRITE (6,*) ' Model run is set up for 365 day calendar.'
            WRITE (6,*) ' Ancillary File is for 360 day calendar.'
          ENDIF
          WRITE (6,*) ' Rerun with correct ancillary file.'
          GO TO 9999   !  Return
        ENDIF

        FILE_LEVELS=1

        IF(I == 1) THEN
          FILE_LEVELS=OZONE_LEVELS
        ELSE IF(I == 2) THEN
          FILE_LEVELS=SM_LEVELS
! This is the maximum value that might be present on the ancillary
! file if it includes soil moisture in layers; otherwise only single
! level data is present and PR_FIXHD will not check value since
! FIXHD(110) will be zero
        ELSE IF(I == 3) THEN
            FILE_LEVELS=ST_LEVELS
        ELSE IF(I == 13) THEN   ! for multilevel aerosols
            FILE_LEVELS=TR_LEVELS
        ELSE IF(I == 14.or.I == 16) THEN   ! for murk and user ancil.
            FILE_LEVELS=P_LEVELS
        ELSE IF(I == 17.or.I == 18) THEN
!           multi-level sulphur cycle ancillary files.
            FILE_LEVELS=P_LEVELS
        else if (i == 25) then
! tropopause-based ozone file with tpps_ozone_levels
           file_levels=tpps_ozone_levels
        END IF


!L 1.4.2 Buffer in integer constants

           IF(FIXHD(100,I) >  0) THEN

! Check for error in file pointers

! Check validity of integer data and print out information
! All files except ozone should contain full fields

            IF(INTHD(6,I) /= ROW_LENGTH) THEN
! Ozone may contain zonal mean data
! also applies to tropopause-based ozone -- no 25.
              IF(.not.((I == 1) .or. (i  ==  25))                       &
     &        .OR. INTHD(6,I) /= 1) THEN
                ICODE=4
                CMESSAGE='INANCILA:integer header error'
                WRITE(6,*) ' INTHD(6) : ',INTHD(6,I),' ?'
                RETURN
              END IF
            END IF

            IF(INTHD(7,I) /= P_ROWS.AND.(I == 9.AND.INTHD               &
     &        (7,I) /= U_ROWS)) THEN
              ICODE=5
              CMESSAGE='INANCILA:integer header error'
              WRITE(6,*) ' INTHD(7) : ',INTHD(7,I),' ?'
              RETURN
            END IF

            IF (I == 1 .or. i  ==  25) THEN  ! Ozone or tpps-ozone file
              WRITE (6,*) ' '
              IF (INTHD(6,I) == 1)THEN
                WRITE (6,*) ' OZONE file contains zonal mean data for ',&
     &          INTHD(6,I),' points x ',INTHD(7,I),' rows'
              ELSEIF (INTHD(6,I) == ROW_LENGTH)THEN
                WRITE (6,*) ' OZONE file contains full fields for ',    &
     &          INTHD(6,I),' points x ',INTHD(7,I),' rows'
              ENDIF
! Check that correct ozone file has been provided.
              IF ((ZonAvOzone .and. i  ==  1) .or.                      &
     &          (ZonAvTppsOzone .and. i  ==  25)) THEN
!! Where is ZonAvOzone (and ZonAvTppsOzone) defined.
                IF (INTHD(6,I) /= 1) THEN
                  WRITE (6,*) ' Zonal Ozone Data is expected',          &
     &            ' for 1 point x ',P_ROWS,' rows'
                  ICODE = 51
                  CMESSAGE = 'INANCA1A : Wrong Ozone data provided.'
                  GO TO 9999   !  Return
                ENDIF
              ELSE
                IF (INTHD(6,I) /= ROW_LENGTH) THEN
                  WRITE (6,*) ' Ozone Data is expected for ',           &
     &            ROW_LENGTH,' points x ',P_ROWS,' rows.'
                  ICODE = 52
                  CMESSAGE = 'INANCA1A : Wrong Ozone data provided.'
                  GO TO 9999   !  Return
                ENDIF
              ENDIF
            ENDIF

          END IF

!L 1.4.3 Buffer in real constants

          IF(FIXHD(105,I) >  0) THEN

! Check validity of real header and print out information

           DO J=1,6
             IF(REALHD(J,I) >  (A_REALHD(J)+0.1).OR.                    &
     &         REALHD(J,I) <  (A_REALHD(J)-0.1))THEN
             IF(I /= 1.OR.(J /= 1.AND.J /= 4))THEN
               WRITE(6,*)(REALHD(K,I),K=1,6),(A_REALHD(K),K=1,6)
               ICODE=8
               CMESSAGE='INANCILA: REAL header Error.'
               RETURN
             END IF
             END IF
           END DO

         END IF

!L 1.4.4 Buffer in level dependent constants if required
!        Not retained in model after initial check

         IF(FIXHD(110,I) >  0) THEN

!L Only files 1 (Ozone), and 3 (Soil temperature)should contain multi
!L level data. File 2 (Soil moisture,snow depth,fractional snow time
!L and soil moisture in layers) may possibly also have multi level data.
!L FILES 13,14,16 (aerosols, murkiness, user ancil.) may also have
!L  multi level data.
!! Files 25 (TppsOzon) may also contain multi-level data
!! File 46 has multilevel data for cariolle ozone scheme

           If (I == 1   .or.                                            &
                                !  Ozone File
     &         I == 14  .or.                                            &
                                !  Murkiness File
     &         I == 16  .or.                                            &
                                !  User Ancillary File                  
     &         I == 46) Then    !  cariolle ozone File  

! Check that ancillary file is set up for correct vertical levels

            If (fixhd(111,I)-1 /= p_levels) Then
              icode=110
              write (CMESSAGE,*) ' Ancillary File set up for wrong',    &
     &        ' no of model levels. Anc ',fixhd(111,I)-1,               &
     &        ' Model ',p_levels
! DEPENDS ON: ereport
              Call Ereport ( RoutineName, icode, Cmessage )
            End If

            l_vert_mismatch = .false.

! Check eta_theta and eta_rho

            Do j=1,p_levels+1
              If (LNER( LEVDEPC(J), A_LEVDEPC(J,1) )) Then
                l_vert_mismatch = .true.
                exit
              End If
            End Do

            Do j=1,p_levels
              If (LNER( LEVDEPC(FIXHD(111,I)+J), A_LEVDEPC(J,2) )) Then
                l_vert_mismatch = .true.
                exit
              End If
            End Do

! Abort if there is a mis-match

            If (l_vert_mismatch) then
              write (6,*) 'Mismatch in vertical levels between model ', &
     &                    'and Ancillary File.'
              write (6,*) 'Anc File : ',title(i)
              write (6,*) 'Eta_Theta - Model'
              write (6,'(5F10.7)') (a_levdepc(k,1),k=1,p_levels+1)
              write (6,*) 'Eta_Theta - Anc File'
              write (6,'(5F10.7)') (levdepc(k),k=1,p_levels+1)
              write (6,*) 'Eta_Rho   - Model'
              write (6,'(5F10.7)') (a_levdepc(k,2),k=1,p_levels)
              write (6,*) 'Eta_Rho   - Anc File'
              write (6,'(5F10.7)') (levdepc(p_levels+1+k),k=1,p_levels)
                   ICODE=11
              Write (CMESSAGE,*) 'Mismatch in LEVDEPC ',                &
     &        'between model and Ancillary File.'
! DEPENDS ON: ereport
              Call Ereport ( RoutineName, icode, Cmessage )
            End If

           else if (i  ==  25) then !! tropopause-based ozone
             !! no checks to run....

           Else If (I == 2) Then  !  Soil Moisture File

             if (PrintStatus >= PrStatus_Diag .and. mype == 0 )then
               write (6,*)
               write (6,*) 'SoilDepths = ',SoilDepths
               write (6,*) 'SM_Levels  = ',sm_levels
               do j=1,sm_levels
                 write (6,*) 'model ',A_LEVDEPC(J,SoilDepths),          &
     &                       ' anc ',LEVDEPC(fixhd(111,I)*3+J)
               enddo
             endif

! Check Soil moisture levels

             Do J=1,SM_LEVELS
               If (LNER(LEVDEPC(fixhd(111,I)*3+J),                      &
     &                  A_LEVDEPC(J,SoilDepths))) Then
                 ICODE=12
                 CMESSAGE='INANCILA: error in LEVDEPC.'
                 RETURN
               End If
             End Do

           Else If (I == 3) Then  !  Deep Soil Temperature File

             If (PrintStatus >= PrStatus_Diag .and. mype == 0) Then
               write (6,*)
               write (6,*) 'SoilDepths = ',SoilDepths
               write (6,*) 'st_levels  = ',st_levels
               do j=1,st_levels
                 write (6,*) 'model ',A_LEVDEPC(J,SoilDepths),          &
     &                       ' anc ',LEVDEPC(fixhd(111,I)*3+J)
               End Do
             End If

! Check Deep Soil levels

             Do J=1,ST_LEVELS
               If (LNER(LEVDEPC(fixhd(111,I)*3+J),                      &
     &                  A_LEVDEPC(J,SoilDepths))) Then
                 ICODE=13
                 CMESSAGE='INANCILA: error in LEVDEPC.'
                 RETURN
               End If
             End Do

!L If aerosol file, check against model levels

           ELSE IF (I == 13) THEN

             DO J=1,TR_LEVELS
               DO J1=1,4
                 IF(LNER(LEVDEPC(J+(J1-1)*FIXHD(111,I)),A_LEVDEPC       &
     &                   (J,J1))) THEN
      WRITE(6,*)'Error in level dependent constants:Level=',J
                   WRITE(6,*)'Position=',J1
                   WRITE(6,*)'Value in model =',A_LEVDEPC(J,J1)
                   WRITE(6,*)'Value in ancillary data =',LEVDEPC(J+     &
     &                             (J1-1)*FIXHD(111,I))
                   ICODE=16
               CMESSAGE='INANCILA: error in LEVDEPC.'
                   RETURN
                 END IF
               END DO
             END DO

           END IF  !  If I

         END IF  !  If Fixhd(110,I) > 0

!L 1.4.5 Buffer in lookup table
! Set start position of boundary fields for file

         IF(FIXHD(150,I) >  0) THEN


           NREC_A = 0
           NREC_S = 0
           DO J = 1,FIXHD(152,I)
             IF (LOOKUP(MODEL_CODE,LOOKUPS+J)  ==  0 .or.               &
     &           LOOKUP(MODEL_CODE,LOOKUPS+J)  ==  imdi) THEN
               STASH_CODE = LOOKUP(ITEM_CODE,LOOKUPS+J)
               IF ((STASH_CODE >= 177 .and. STASH_CODE <= 179) .or.     &
     &             (STASH_CODE >= 210 .and. STASH_CODE <= 212)) THEN
                 LOOKUP(MODEL_CODE,LOOKUPS+J) = slab_im
                 NREC_S = NREC_S+1
               ELSE
                 LOOKUP(MODEL_CODE,LOOKUPS+J) = atmos_im
                 NREC_A = NREC_A+1
               END IF
             END IF
           END DO
           IF (NREC_A >  0) THEN
             WRITE (6,*) ' '
             WRITE (6,*) ' INANCA1A : submodel_id in ',NREC_A,          &
     &       ' records set to atmos_im in ancillary file ',I
           ENDIF
           IF (NREC_S >  0) THEN
             WRITE (6,*) ' '
             WRITE (6,*) ' INANCA1A : submodel_id in ',NREC_S,          &
     &       ' records set to slab_im in ancillary file ',I
           ENDIF

         END IF

         LOOKUPS=LOOKUPS+FIXHD(152,I)

       ELSE

!L  If file not required, zero fixed length header
         DO J=1,LEN_FIXHD
      FIXHD(J,I)=0
         END DO

         LOOKUP_START(I)=LOOKUPS+1
       END IF

      END DO

!L 1.5 Set positions in main data blocks

      DO I=1,NANCIL_FIELDS
        IF (STASHANCIL(I)  >   0) THEN
        IF (MODEL_CODES_ANCIL(I) == SLAB_IM) THEN
          D1_ANCILADD(I)=SI_SLAB(STASHANCIL(I))
        ELSE
          !Using si array directly to allow all sections - kdcorbin, 05/10
          !D1_ANCILADD(I)=SI_ATMOS(STASHANCIL(I))
          isec=stashancil(i)/1000
          iind=stashancil(i)-isec*1000
          D1_ANCILADD(I) = SI(iind,isec,1)
        ENDIF
        ELSE
          D1_ANCILADD(I)=0
        ENDIF
      ENDDO

!L 1.51 If a request is made to update a field, ensure that space for
!L     that field has been allocted in D1.

      DO I=1,NANCIL_FIELDS
        IF((FIELDCODE(1,I) >  0).AND.(D1_ANCILADD(I) <= 1)) THEN
          WRITE(6,*)' An address in D1 has not been set for ancillary   &
     & field number ',I
          ICODE=30
          CMESSAGE='INANCILA: updating for ancillary field is requested &
     & but no space has been allocated in D1'
          RETURN
        ENDIF
      END DO

!L 1.6 Set positions of data

      DO I=1,NANCIL_FIELDS
      NLOOKUP(I) =0
      LOOKUP_STEP(I)=0

! If LOOKUP_START=0 for file FILEANCIL(I), no fields required.
        IF   (STASHANCIL(I)  >   0) THEN
        IF (LOOKUP_START(FILEANCIL(I)) >  0) THEN

        DO J=LOOKUP_START(FILEANCIL(I)),LOOKUPS

          IF (LOOKUP(ITEM_CODE,J) == STASHANCIL(I)) THEN
            NLOOKUP(I)=J-LOOKUP_START(FILEANCIL(I))+1
            GOTO 161
          END IF

        END DO

! Find second occurence of data to set LOOKUP_STEP

161     LOOKUP_STEP(I)=0


        IF(J <  LOOKUPS) THEN

          DO J1=J+LEVELS(I),LOOKUPS
            IF (LOOKUP(ITEM_CODE,J1) == STASHANCIL(I)) THEN
              LOOKUP_STEP(I)=J1-NLOOKUP(I)-LOOKUP_START(FILEANCIL(I))+1
              GOTO 164
            END IF
          END DO
164      CONTINUE
        END IF

        END IF
        ENDIF

      END DO

!L SET LEVELS=2 FOR ICE FRACTION AND SNOW DEPTH, TO INDICATE PRESCENCE
!L fractional time fields

      LEVELS(27)=2
      LEVELS(9)=2

 900  CONTINUE
 9999 CONTINUE
      RETURN
      END SUBROUTINE INANCILA
#endif
