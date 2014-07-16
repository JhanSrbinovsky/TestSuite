
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Read in Ancillary Fields
!
!  Subroutine REPLANCA  - Read in Ancillary Fields
!
! Description:
!   Read in Ancillary Fields
!
! Method:
!   This routine reads in the required ancillary fields. Any time
!   interpolation is done as required. Both 360 and 365 day calendars
!   are catered for. The ANCILmaster files are used in the
!   reconfiguration.
!
! Current Code Owner: D. Robinson
!

       SUBROUTINE REPLANCA(I_YEAR,I_MONTH,I_DAY,I_HOUR,                 &
     &                     I_MINUTE,I_SECOND,                           &
     &                     P_FIELD,P_ROWS,U_FIELD,V_FIELD,RR_FIELD,     &
     &                     LAND_FIELD,                                  &
     &                     D1,LAND,                                     &
     &                     ICE_FRACTION,TSTAR,FLAND,                    &
     &                     TSTAR_LAND_CTILE,TSTAR_SEA_CTILE,            &
     &                     TSTAR_SICE_CTILE,                            &
     &                     TSTAR_ANOM,                                  &
     &                     NS_SPACE,FIRST_LAT,                          &
     &                     LEN1_LOOKUP,LEN_FIXHD,LEN_INTHD,             &
     &                     LEN_D1,FIXHD,INTHD,                          &
     &                     LOOKUP,RLOOKUP,LOOKUP_START,                 &
     &                     NLOOKUPS,                                    &
     &                     ICODE,CMESSAGE)

      Use ancil_mod, Only :                                             &
     &    Anc_File,       ancFiles,                                     &
     &    Anc_Record,     ancRecs,                                      &
     &    AncF_UnitNo,    ancil_add,                                    &
     &    StashAncil,     Levels,                                       &
     &    Lookup_Step,    NLookup

      Use Locate_Anc_mod, Only :                                        &
     &    Locate_Anc_Field,                                             &
     &    Locate_Anc_File

      Use Rcf_DecompTP_Mod, Only :                                      &
     &    Decomp_rcf_output

      Use Rcf_Recon_Mod, Only :                                         &
     &    lcal360,                                                      &
     &    lamipii,                                                      &
     &    l_sstanom

      Use Rcf_Parvars_Mod, Only :                                       &
     &    mype

      Use Rcf_PrintStatus_Mod, Only :                                   &
     &    PrintStatus, PrStatus_Diag

      Use Rcf_CntlAtm_Mod, Only :                                       &
     &    L_CTILE, LTLEADS, L_UKCA

      Use Csenario_Mod ! taken from r3540

      !kdcorbin, 05/10
      Use Rcf_Items_Mod, only :  &
          num_items,  &
          sctn_array,item_array,upaf_array

      Use Rcf_FortranIO_Mod, Only : &
          Max_Filename_Len

      IMPLICIT NONE

! Arguments
      Integer :: I_YEAR          ! Curent Model Time
      Integer :: I_MONTH         !   "      "     "
      Integer :: I_DAY           !   "      "     "
      Integer :: I_HOUR          !   "      "     "
      Integer :: I_MINUTE        !   "      "     "
      Integer :: I_SECOND        !   "      "     "

      Integer :: P_FIELD         ! Size of horizontal fields
      Integer :: P_ROWS          !
      Integer :: U_FIELD         !   "  "      "         "
      Integer :: V_FIELD         !   "  "      "         "
      Integer :: RR_FIELD        !   "  "      "         "
      Integer :: Land_Field      !   "  "      "         "

      Integer :: NLOOKUPS        ! Number of lookup tables
      Integer :: Len_FixHd       ! Length of Fixed Header
      Integer :: Len_IntHd       ! Length of Integer Header
      Integer :: Len1_Lookup     ! First dimension of Lookup
      Integer :: LEN_D1          ! Size of primary data array

      Integer FIXHD(LEN_FIXHD,ancFiles)    ! Anc fixed headers
      Integer INTHD(LEN_INTHD,ancFiles)    ! Anc Integer Headers
      Integer LOOKUP(LEN1_LOOKUP,NLOOKUPS) ! Anc Lookup Tables
      Integer LOOKUP_START(ancFiles)       ! Start of lookup tables for
                                           ! anc files opened.
      Integer ICODE                        ! Return Code
      Integer IPOS_111                     ! Land fraction
!                                          ! ancil position


      Real NS_Space              ! NS latitude spacing
      Real First_Lat             ! Latitude of first gridpoint
      Real D1(Len_D1)            ! INOUT Data array to hold fields
                                 !       except TStar and Ice Fraction
      Real Ice_Fraction(P_Field) ! INOUT Ice Fraction, updated if
                                 !       requested
      Real TStar (P_Field)       ! INOUT T Star, updated if requested
      Real TStar_Land_Ctile (P_Field)
!                                ! INOUT T*_land, updated if requested
      Real TStar_Sea_Ctile (P_Field)
!                                ! INOUT T*_sea, updated if requested
      Real TStar_Sice_Ctile (P_Field)
!                                ! INOUT T*_sice, updated if requested
      Real TStar_Anom(P_Field)   ! INOUT SST Anolomy, formed in Recon.
                                 !       Added in model run if requested

      Real RLookup(Len1_Lookup,NLookups) ! REAL copy of Lookup Table

      Logical Land (P_Field)     ! Land mask

      Character (Len=80) :: CMessage

! Local Variables

!                                  ! Buffers to hold values of ancillary
!                                  ! data for time interpolation.
!                                  ! Field of ancillary data held prior
!                                  ! to selective updating.
      Real, Dimension (:), Allocatable :: ANCIL1
      Real, Dimension (:), Allocatable :: ANCIL2
      Real, Dimension (:), Allocatable :: ANCIL_DATA

      Real SNOW_CHANGE(P_FIELD)   ! Fractional time of change of
                                  ! snow cover
      Real ICE_EXTENT(P_FIELD,2)  ! Fractional time of change
                                  ! of ice cover
      Real PRES_VALUE(P_FIELD)    ! Prescribed value of data when
                                  ! controlling field is zero.
      Real NO_ICE_EXTENT(P_FIELD) ! Indicator for no sea ice
                                  ! =0 if ice cover

      Real TStar_Land (P_Field)  ! T*_land, updated if requested
      Real TStar_Sea (P_Field)   ! T*_sea, updated if requested
      Real TStar_Sice (P_Field)  ! T*_sice, updated if requested
      Real TStar_Ssi (P_Field)   ! T*_ssi, updated if requested
      Real Flandg (P_Field)      ! Land fraction in gridbox
      Real Fland (Land_Field)    ! Land fraction in gridbox
!                                ! (on land-only points).
      Integer I,J,L              ! Loop indices

      Integer I1,I2,I3,II       ! Array indices
      Integer ID                !
      Integer IM                !
      Integer IY                !
      Integer FIELD             ! Current Ancil Ref Number.
      Integer FILE              ! Anc File number
      Integer IREC              ! Anc record number

      Integer INTERVAL          ! Interval between data times
      Integer STEP              ! Number of data times skipped.
      Integer MONTHS            ! )Used in calculation of position
      Integer HOURS             ! )of data required
      Integer PERIOD            ! Period of periodic data
      Integer START_MONTH       !
      Integer LEVEL             ! Loop index for levels
      Integer ANCIL_REF_DAYS    ! Ancil.reference time in whole days
      Integer ANCIL_REF_SECS    ! Ancil.reference time in extra seconds
      Integer DAY,SEC           ! Times relative to reference time
      Integer LEN
      Integer ROW_LENGTH
      Integer I_YEAR1          ! Copy of Curent Model Time year
      Integer I_MONTH1         !   "      "     "          month
      Integer I_DAY1           !   "      "     "          day
      Integer I_HOUR1          !   "      "     "          hour

      Integer len_anc_env_var  !  Length of env var
      Integer anc_file_open    !  Stores which ancillary file is open.
      Integer ipos_27          !  Pos of anc record for Sea Ice Frac
      Integer ipos_28          !  Pos of anc record for SST
      Integer ipos_30          !  Pos of anc record for surf u-current
      Integer ipos_31          !  Pos of anc record for surf v-current

      Integer FILEANCIL(ancRecs) ! Anc File No for anc record

      Integer FIELD_SIZE       !  Stores field size allowing use of
                               !  P_FIELD or RR_FIELD in calls.

      Logical LINTERPOLATE      ! Indicates whether time
                                ! interpolation needed.
      Logical LT_INT_C          ! Indicates use of controlled time
                                ! interpolation
      Logical LMISMATCH         ! Used in header checks
      Logical SINGLE_TIME       ! Indicates that only one time is
                                ! available in data set
      Logical PERIODIC          ! Data set is periodic
      Logical REGULAR           ! Interval between data times in
                                ! dataset is regular in model timesteps.
      Logical LICE_DEPTH        ! T : Field is Ice Depth
      Logical LICE_FRACTION     ! T : Field is Ice Fraction
      Logical LSNOW_DEPTH       ! T : Field is Snow depth

      Logical UPDATE (ancRecs)  ! T : Anc Field to be updated

      Logical Sea (P_Field)      ! Sea mask

      Real ZERO         !
      Real TIME1        !  Times if data used in time interpolation
      Real TIME2        !
      Real TIME         !  Target time for time interpolation
      Real LAT_P        !  Latitude of point

!kdcorbin, 05/10
Character (Len=Max_Filename_Len) :: AncFileName

! TM, TFS
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------
!------------------------ nstypes.h ----------------------------------
!jhan:further renovation of ths file may be necessary params are dependent on dataset
!jhan: ALSO nstypes_cable.h should be unecessary nsoil/soil is only difference
      !--- Number of non-vegetation surface types
      Integer, Parameter :: NNVG  = 4

      !--- Number of plant functional types.
      Integer, Parameter :: NPFT  = 13
      
      !--- Number of surface types.
      Integer, Parameter :: NTYPE =17 
      
      !--- Index of the surface type 'Soil'
      !Integer, Parameter :: SOIL  = 16 
      !dhb599, 20110615: change made as per Peter Vohralik, item 1:
      Integer, Parameter :: SOIL  = 14

!--- Land surface types :
!--- original veg. tiles 
!     1 - Broadleaf Tree
!     2 - Needleleaf Tree
!     3 - C3 Grass
!     4 - C4 Grass
!     5 - Shrub
!--- for testing these tiles are set = 1:5 
!     6 - Broadleaf Tree
!     7 - Needleleaf Tree
!     8 - C3 Grass
!     9 - C4 Grass
!    10 - Shrub
!--- for testing these tiles are set = 0
!    11 - 0 
!    11 - 0
!    11 - 0
!--- original non-veg tiles moved to these indices
!     14 - Urban
!     15 - Water
!     16 - Soil
!     17 - Ice


!*L------------------ COMDECK LOOKADD ----------------------------------
!LL
!LL Purpose : Contains information about the format
!LL           of the PP header
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   4.0  12/09/95   Change NPERIODS to LBUSER3, BRSVD1 to BULEV,
!LL                   BRSVD2 to BHULEV and definitions for BRLEV and
!LL                   BHRLEV. Corresponding changes made to STWORK1A
!LL                   and PPHEAD1A. (Andrew Brady)
!LL  4.0  12/10/95  Change item 45 from lbuser7 to model_code. RTHBarnes
!LL  5.1  17/04/00    Fixed/Free format. P.Selwood.
!LL  5.2  25/09/00    Add LBCC_xxxx variables for the compressed
!LL                   LBC LOOKUP array                  P.Burton
!LL
!LL Programming standard :
!LL
!LL Logical components covered : F092
!LL
!LL Project task :
!LL
!LL External documentation:
!LL
!LLEND -----------------------------------------------------------------
!

! Validity time
      Integer, Parameter :: LBYR   =1   ! Year
      Integer, Parameter :: LBMON  =2   ! Month
      Integer, Parameter :: LBDAT  =3   ! Day of month
      Integer, Parameter :: LBHR   =4   ! Hour
      Integer, Parameter :: LBMIN  =5   ! Minute
      Integer, Parameter :: LBDAY  =6   ! Day number

! Data time
      Integer, Parameter :: LBYRD  =7   ! Year
      Integer, Parameter :: LBMOND =8   ! Month
      Integer, Parameter :: LBDATD =9   ! Day of month
      Integer, Parameter :: LBHRD  =10  ! Hour
      Integer, Parameter :: LBMIND =11  ! Minute
      Integer, Parameter :: LBDAYD =12  ! Day number

      Integer, Parameter :: LBTIM  =13  ! Time indicator
      Integer, Parameter :: LBFT   =14  ! Forcast period (hours)
      Integer, Parameter :: LBLREC =15  ! Length of data record
      Integer, Parameter :: LBCODE =16  ! Grid type code
      Integer, Parameter :: LBHEM  =17  ! Hemisphere indicator
      Integer, Parameter :: LBROW  =18  ! Number of rows in grid
      Integer, Parameter :: LBNPT  =19  ! Number of points per row
      Integer, Parameter :: LBEXT  =20  ! Length of extra data
      Integer, Parameter :: LBPACK =21  ! Packing method indicator
      Integer, Parameter :: LBREL  =22  ! Header release number
      Integer, Parameter :: LBFC   =23  ! Field code
      Integer, Parameter :: LBCFC  =24  ! Second field code
      Integer, Parameter :: LBPROC =25  ! Processing code
      Integer, Parameter :: LBVC   =26  ! Vertical coordinate type
      Integer, Parameter :: LBRVC  =27  ! Coordinate type for reference
                                        ! level

      Integer, Parameter :: LBEXP  =28  ! Experiment number
      Integer, Parameter :: LBEGIN =29  ! Start record
      Integer, Parameter :: LBNREC =30  ! No of records-Direct access
                                        ! only
      Integer, Parameter :: LBPROJ =31  ! Met-O-8 projection number
      Integer, Parameter :: LBTYP  =32  ! Met-O-8 field type
      Integer, Parameter :: LBLEV  =33  ! Met-O-8 level code
      Integer, Parameter :: LBRSVD1=34  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD2=35  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD3=36  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD4=37  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBSRCE =38  ! =1111 to indicate following
                                        ! apply to UM
      Integer, Parameter :: DATA_TYPE =39  ! Indicator for real/int or
                                           ! timeseries
      Integer, Parameter :: NADDR  =40  ! Start address in DATA_REAL or
                                        ! DATA_INT
      Integer, Parameter :: LBUSER3=41  ! Free for user-defined function
      Integer, Parameter :: ITEM_CODE =42  !Stash code
      Integer, Parameter :: LBPLEV =43  ! Pseudo-level indicator (if
                                        ! defined)
      Integer, Parameter :: LBUSER6=44  ! Free for user-defined function
      Integer, Parameter :: MODEL_CODE =45 ! internal model identifier

      Integer, Parameter :: BULEV  =46  ! Upper level boundary
      Integer, Parameter :: BHULEV =47  ! Upper level boundary
      Integer, Parameter :: BRSVD3 =48  ! Reserved for future PP-package
                                        ! use
      Integer, Parameter :: BRSVD4 =49  ! Reserved for future PP-package
                                        ! use
      Integer, Parameter :: BDATUM =50  ! Datum value
      Integer, Parameter :: BACC   =51  ! (Packed fields) Packing
                                        ! accuracy
      Integer, Parameter :: BLEV   =52  ! Level
      Integer, Parameter :: BRLEV  =53  ! Lower level boundary
      Integer, Parameter :: BHLEV  =54  ! (Hybrid levels) A-level of
                                        ! value
      Integer, Parameter :: BHRLEV =55  ! Lower level boundary
      Integer, Parameter :: BPLAT  =56  ! Real latitude of 'pseudo'
                                        ! N Pole
      Integer, Parameter :: BPLON  =57  ! Real longitude of 'pseudo'
                                        ! N Pole
      Integer, Parameter :: BGOR   =58  ! Grid orientation
      Integer, Parameter :: BZY    =59  ! Zeroth latitude
      Integer, Parameter :: BDY    =60  ! Latitude interval
      Integer, Parameter :: BZX    =61  ! Zeroth longitude
      Integer, Parameter :: BDX    =62  ! Longitude interval
      Integer, Parameter :: BMDI   =63  ! Missing data indicator
      Integer, Parameter :: BMKS   =64  ! M,K,S scaling factor

      Integer, Parameter :: LBCC_LBYR   = 1  ! Year
      Integer, Parameter :: LBCC_LBMON  = 2  ! Month
      Integer, Parameter :: LBCC_LBDAT  = 3  ! Day of the month
      Integer, Parameter :: LBCC_LBHR   = 4  ! Hour
      Integer, Parameter :: LBCC_LBMIN  = 5  ! Minute
      Integer, Parameter :: LBCC_LBDAY  = 6  ! Day number
      Integer, Parameter :: LBCC_LBEGIN = 7  ! Start record
      Integer, Parameter :: LBCC_NADDR  = 8  ! Start address of DATA
! Mapping of MPP_LOOKUP; analogous to mapping in PP header

      Integer, Parameter :: P_NADDR=1    ! Address on local PE
      Integer, Parameter :: P_LBLREC=2   ! Local length of record

!*----------------------------------------------------------------------
! NADDR IS LOCATION IN PP-HEADER (LOOKUP) FOR START POSN OF VARIABLE
! ITEM_CODE is the location in PP header for a code defined as
!           (section number)*1000+item number
! DATA_TYPE is the location in the PP header defining data as REAL or
!           INTEGER.
! LBNPT is the location defining the number of points per row
!

      EXTERNAL                                                          &
     &       TIME2SEC,                                                  &
     &       Rcf_ReadFlds,                                              &
     &       T_INT,                                                     &
     &       T_INT_C

!   0. Set up local arrays for ANCIL data to required size.
      FIELD_SIZE=MAX(P_FIELD,RR_FIELD)

      If (PrintStatus >= PrStatus_Diag .and. mype == 0 ) Then
        write (6,*) "Rcf_ReplancA: chosen FIELD_SIZE=",FIELD_SIZE
      End If

      Allocate ( ANCIL1(FIELD_SIZE) )
      Allocate ( ANCIL2(FIELD_SIZE) )
      Allocate ( ANCIL_DATA(FIELD_SIZE) )

!L  1.  Initialisation for atmosphere

      ICODE=0
      CMESSAGE=' '
      anc_file_open = 0

!     Initialise ANCIL1/2. Includes Halos for mpp runs.
      ANCIL1(:)=0.0
      ANCIL2(:)=0.0

! Read in fractional land field first:
! use ANCIL_DATA as temporary storage

       If(L_CTILE)Then

         call locate_anc_field (111, ipos_111)

         If(anc_record(ipos_111) % anc_field_read == 1)Then !read ancil

           call locate_anc_file(anc_record(ipos_111) % anc_file_number, FILE)
           len_anc_env_var = len_trim ( anc_file(file) % anc_env_var )
! DEPENDS ON: file_open
           call file_open (AncF_UnitNo,                                 &
     &                      anc_file(file) % anc_env_var,               &
     &                      len_anc_env_var,0,0,icode)

! DEPENDS ON: rcf_readflds
           Call Rcf_ReadFlds                                            &
     &                 (AncF_UnitNo,1,NLOOKUP(ipos_111),                &
     &                  LOOKUP(1,LOOKUP_START(FILE)),                   &
     &                  LEN1_LOOKUP,ANCIL_DATA,P_FIELD,FIXHD(1,FILE),   &
     &                  ICODE,CMESSAGE)

           WRITE(6,*)'READ IN LAND FRACTION'

            Do I=1,P_FIELD
              D1(ANCIL_ADD(IPOS_111)+I-1)=ANCIL_DATA(I)
            End Do

            Do I=1,P_FIELD
              FLANDG(I)=0.0
              If(LAND(I))FLANDG(I)=ANCIL_DATA(I)
            End Do

          Else             ! Land frac already read in from input dump

            L=0
            Do I=1,P_FIELD
              FLANDG(I)=0.0
              If(LAND(I))Then
                L=L+1
                FLANDG(I)=FLAND(L)
              End If
            End Do
          End If

          Do I=1,P_FIELD
! If land or sea fraction is less than machine tolerance print warning
            If(LAND(I).AND.FLANDG(I) <  Epsilon(1.0))Then
              WRITE(6,*)'*****************WARNING********************'
              WRITE(6,*)'LAND FRACTION IS LESS THAN MACHINE TOLERANCE'
            End If
            If(.NOT.LAND(I).AND.1.0-FLANDG(I) <  Epsilon(1.0))Then
              WRITE(6,*)'*****************WARNING********************'
              WRITE(6,*)'SEA FRACTION IS LESS THAN MACHINE TOLERANCE'
            End If
    !
            If(FLANDG(I) <= 0.0.AND.LAND(I))Then
              WRITE(6,*)'*ERROR* a) LAND FRAC & LAND MASK INCONSISTENT'
              ICODE = 800
              CMESSAGE='REPLANCA:ERR:LAND FRAC & MASK ARE INCONSISTENT'
            End If
            If(FLANDG(I) >  0.0.AND..NOT.LAND(I))Then
              WRITE(6,*)'*ERROR* b) LAND FRAC & LAND MASK INCONSISTENT'
              ICODE = 801
              CMESSAGE='REPLANCA:ERR:LAND FRAC & MASK ARE INCONSISTENT'
            End If

          End Do

        Else                     ! Not coastal tiling:
          Do I=1,P_FIELD
            If(LAND(I))Then
              FLANDG(I)=1.0
            Else
              FLANDG(I)=0.0
            End If
          End Do
        End If                   ! End of coastal tiling loop


        Do I=1,P_FIELD
          If(FLANDG(I) <  1.0)Then
            SEA(I)=.TRUE.
          Else
            SEA(I)=.FALSE.
          End If
        End Do


!     Set up surface temperatures:
       If(L_CTILE)Then
         Do I=1,P_FIELD
            TSTAR_LAND(I)=TSTAR_LAND_CTILE(I)
            TSTAR_SEA(I)=TSTAR_SEA_CTILE(I)
            TSTAR_SICE(I)=TSTAR_SICE_CTILE(I)
            If(ICE_FRACTION(I) <= 0.0)Then
              TSTAR_SSI(I)=TSTAR_SEA(I)
            Else
              TSTAR_SSI(I)=ICE_FRACTION(I)*TSTAR_SICE(I)                &
     &          +(1.0-ICE_FRACTION(I))*TSTAR_SEA(I)
            End If
         End Do
       Else
         Do I=1,P_FIELD
            TSTAR_LAND(I)=TSTAR(I)
            TSTAR_SSI(I)=TSTAR(I)
         End Do
       End If

!L  1.1 Set logical UPDATE for each ancillary field independently

!  FILEANCIL not required. Remove later.

      Do i=1,ancRecs
        call locate_anc_file( anc_record(i) % anc_file_number, FILEANCIL(i) )
        UPDATE(i)    = anc_record(i) % anc_field_read > 0
      End Do

!  Initialise for valid time interpolation in reconfiguration mode.
      ANCIL_REF_DAYS = 0
      ANCIL_REF_SECS = 0

!L 1.2 Allow for dependencies between fields
! Sea surface temperature must be updated when sea ice is updated

      call locate_anc_field (27, ipos_27)
      call locate_anc_field (28, ipos_28)
      call locate_anc_field (30, ipos_30)
      call locate_anc_field (31, ipos_31)

      UPDATE(ipos_28) = UPDATE(ipos_27) .OR. UPDATE(ipos_28)

! Both surface current components must be updated together

      UPDATE(ipos_30) = UPDATE(ipos_30) .OR. UPDATE(ipos_31)
      UPDATE(ipos_31) = UPDATE(ipos_30)

!L Select method of time interpolation for SST. The interpolation
!L allows for sea ice if ice data is available at the same times
!L as the temperature data. Otherwise linear interpolation is used.

      LT_INT_C=.TRUE.

      If (UPDATE(ipos_28)) Then

      If(FIXHD(10,FILEANCIL(ipos_27)) == 0) LT_INT_C=.FALSE.
        If(LT_INT_C) Then
        Do I=21,41
          If(FIXHD(I,FILEANCIL(ipos_27)) /= FIXHD(I,                    &
     &      FILEANCIL(ipos_28))) Then
            LT_INT_C=.FALSE.
            WRITE(6,*)' WARNING:controlled time interpolation for SST', &
     &      ' not available: Mismatch in SST and SEA-ICE ancillary data'&
     &     ,' times in FIXED HEADER'
            WRITE(6,*)' position=',I,                                   &
     &                ' SEA-ICE=',FIXHD(I,FILEANCIL(ipos_27))
            WRITE(6,*)' position=',I,                                   &
     &                ' SST    =',FIXHD(I,FILEANCIL(ipos_28))
          End If
        End Do
        End If
      End If

!L Loop over ancillary fields(atmosphere)

      Do irec=1,ancRecs

        field = anc_record(irec) % ancil_ref_number

        LICE_DEPTH=field == 29  ! required for LAMIPII

        If (UPDATE(irec)) THEN  ! (1st level IF)

          call locate_anc_file(anc_record(irec) % anc_file_number, FILE)

! Should this be mype=0 only ? Failed when put in. Check out.

!        if (mype == 0) then   !  Open the Ancillary file

          If (file  /=  anc_file_open) Then ! File is not open

            !Open with environmental variables
            !Commented - kdcorbin, 05/10
            !len_anc_env_var = len_trim ( anc_file(file) % anc_env_var )

            !! DEPENDS ON: file_open
            !call file_open (AncF_UnitNo,                                &
            !                anc_file(file) % anc_env_var,               &
            !                len_anc_env_var,0,0,icode)

            !Get file name from user_prog_ancil_file - kdcorbin, 05/10
            Call Fort_Get_Env(anc_file(file)%anc_env_var,8,AncFileName, &
                    Max_Filename_Len,icode)

            If (icode /= 0) Then
               do j=1,num_items
                  if (anc_record(irec)%section_number == sctn_array(j) &
                     .and. anc_record(irec)%item_number == item_array(j)) Then
                     AncFileName = upaf_array(j)
                     icode=0
                   endif
               enddo
            endif

            len_anc_env_var = len_trim( AncFileName)
! DEPENDS ON: file_open
            call File_Open(AncF_UnitNo,AncFileName,  &
                 len_anc_env_var,0,1,icode)

            write(6,*) 'File Name: ',AncFileName(1:len_trim(AncFileName))

           If (icode /= 0) Then
              write (6,*) ' problem with file_open.'
              write (6,*) ' icode ',icode
              write (6,*) ' cmessage ',cmessage
              Go To 900
            End If

!           Record ancillary file opened
            anc_file_open = file

          End If ! If file not already opened

!        endif   !   if mype=0

       If(LICE_DEPTH.AND.LAMIPII) Then

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
          Do I=1,P_ROWS
! work out latitude in radians
            LAT_P=FIRST_LAT-NS_SPACE*(I-1)
            Do J=1,ROW_LENGTH
              II=J+(I-1)*ROW_LENGTH
              ANCIL_DATA(II)=0.0
              If (ICE_FRACTION(II) >  0.0) Then
                If (LAT_P >  0.0) THEN   ! Arctic ice depth
                  ANCIL_DATA(II)=2.*ICE_FRACTION(II)
                Else                     ! Antarctic ice depth
                  ANCIL_DATA(II)=1.*ICE_FRACTION(II)
                End If
              End If
            End Do
          End Do
!L     Sea ice thickness
!L       Update over all sea points (all sea ice points are the only
!L       ones strictly required, but this cannot be determined easily)

          Do I=1,P_FIELD
            If(SEA(I)) Then
              D1(ANCIL_ADD(irec)+I-1)=ANCIL_DATA(I)
            End If
          End Do
       Else
!     Update required for field

        If ( mype == 0 ) then
        !kdcorbin, 05/10 - added RCF to print statement
        WRITE(6,*)'REPLANCA-RCF: UPDATE REQUIRED FOR FIELD',FIELD,' : ',    &
     &  anc_record(irec) % anc_name
        End If

          If ( FIXHD(10,FILE)  <   0 .OR. FIXHD(10,FILE)  >   2 ) Then
        write (6,*) ' file ',file
        write (6,*) ' FIXHD(10,file) ',FIXHD(10,file)
            ICODE = 700 + file
            CMESSAGE = 'REPLANCA: Error in fixed header(10) of ancillary&
     & file                           '
            Go To 900
          End If

!L    Check whether more than one data time available in data set

        SINGLE_TIME=FIXHD(10,FILE) == 0

!L    Set default values for time interpolation

        LINTERPOLATE=.TRUE.

        If(SINGLE_TIME) Then
          LINTERPOLATE=.FALSE.
        End If

        If (FIELD >  9 .AND. FIELD <  19) Then
          LINTERPOLATE=.FALSE.
        End If

!L 2.1 Find position of input record

!L    Default settings of search parameters if only one time present

        If(SINGLE_TIME) Then
          STEP=0
        Else

          PERIODIC=FIXHD(10,FILE) == 2
          REGULAR=.TRUE.

      If (.NOT. LCAL360) Then
          REGULAR=FIXHD(35,FILE) == 0.AND.FIXHD(36,FILE) == 0
! i.e. data at intervals of days/hours & non-periodic
          If(PERIODIC) REGULAR=REGULAR.AND.FIXHD(37,FILE) == 0
! i.e. data at intervals of hours & periodic
      End If

!         Error checking on time information.

          If ( FIXHD(35,FILE)  <   0 .OR.                               &
     &         FIXHD(36,FILE)  <   0 .OR. FIXHD(36,FILE)  >   12 .OR.   &
     & REGULAR .AND. ( FIXHD(37,FILE)  <   0 .OR. FIXHD(37,FILE)  >   31&
     &  .OR. FIXHD(38,FILE)  <   0 .OR. FIXHD(38,FILE)  >   24 ) ) Then
!           FIXHD(39-40) are not used by REPLANCA.
!           FIXHD(35-37) have already been used if not CAL360.
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: Error in validity time interval given &
     &in ancillary file (FIXHD(35-38))'
            RETURN
          End If

          If ( FIXHD(21,FILE)  <   0 .AND. .NOT. PERIODIC               &
     &  .OR. .NOT. ( REGULAR .AND. PERIODIC ) .AND.                     &
!    !  If it is REGULAR & PERIODIC more detailed check is applied below
     &     ( FIXHD(22,FILE)  <   0 .OR. FIXHD(22,FILE)  >   12 .OR.     &
     &       FIXHD(23,FILE)  <   0 .OR. FIXHD(23,FILE)  >   31 .OR.     &
     &       FIXHD(24,FILE)  <   0 .OR. FIXHD(24,FILE)  >   24 .OR.     &
     &       FIXHD(25,FILE)  <   0 .OR. FIXHD(25,FILE)  >   60 .OR.     &
     &       FIXHD(26,FILE)  <   0 .OR. FIXHD(26,FILE)  >   60 ) ) Then
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: Error in first validity time given in &
     & ancillary file (FIXHD(21-26))  '
            RETURN
          End If

          If(.NOT.PERIODIC) Then

!L            If data taken from full time series of input data.

! DEPENDS ON: time2sec
            CALL TIME2SEC(I_YEAR,I_MONTH,I_DAY,I_HOUR                   &
     &                    ,I_MINUTE,I_SECOND                            &
     &                    ,ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC        &
     &                    ,LCAL360)


            If (REGULAR) Then

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

              If(HOURS <  0) Then
                ICODE=400+FIELD
           CMESSAGE='REPLANCA: Current time precedes start time of data'
                RETURN
              End If

!L FInd interval(in hours) between ancillary data on file
              INTERVAL=FIXHD(35,FILE)*8640+FIXHD(36,FILE)*720+          &
     &                FIXHD(37,FILE)*24+FIXHD(38,FILE)

! Do not interpolate in time if data time exactly matches model time

              If(MOD(HOURS,INTERVAL) == 0) Then
                LINTERPOLATE=.FALSE.
              End If

              STEP=HOURS/INTERVAL
              TIME=REAL(HOURS)
              TIME1=STEP*INTERVAL
              TIME2=(STEP+1)*INTERVAL

            Else

!L 2.1.2 Gregorian calender;ancillary data interval is in months or
!L       years,which is irregular in model timesteps.
!L original code is inaccurate for this section - corrected code under
!L LAMIPII makes use of dates in lookup headers
!L For a real calendar year the mid-point of each month is different
!L in terms of its hour and day. The old inaccurate method assumes
!L the hour and day are taken from the fixhd values. These are only
!L usually correct for the first month on the ancillary file.

!L FInd interval(in months) between ancillary data on file
              INTERVAL=FIXHD(35,FILE)*12+FIXHD(36,FILE)
              MONTHS=I_YEAR*12+I_MONTH
              START_MONTH=FIXHD(21,FILE)*12+FIXHD(22,FILE)
              MONTHS=MONTHS-START_MONTH
!  Check for time within month
           If (LAMIPII) THEN   ! corrected code uses pp header
              STEP=MONTHS/INTERVAL
!             I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*STEP
              I2=NLOOKUP(irec)+LOOKUP_STEP(irec)*STEP
              I1=I2+LOOKUP_START(FILE)-1
! Check against day and hour of actual lookup header not first field
              If((I_DAY*24+I_HOUR) <                                    &
     &           (LOOKUP(3,I1)*24+LOOKUP(4,I1))) Then
                MONTHS=MONTHS-1
              End If
           Else              ! old less accurate code uses FIXHD
              If((I_DAY*24+I_HOUR) <                                    &
     &           (FIXHD(23,FILE)*24+FIXHD(24,FILE))) Then
                MONTHS=MONTHS-1
              End If
           End If ! LAMIPII

              If(MONTHS <  0) Then
                ICODE=400+FIELD
           CMESSAGE='REPLANCA: Current time precedes start time of data'
                RETURN
              End If


              STEP=MONTHS/INTERVAL

           If (LAMIPII) THEN       ! corrected code
              TIME=REAL(SEC)/3600+REAL(DAY*24)
! correct calculation of dates uses lookup table dates not fixhd date
!             I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*STEP
              I2=NLOOKUP(Irec)+LOOKUP_STEP(irec)*STEP
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

           Else   ! LAMIPII test - old inaccurate code using FIXHD
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
           End If     ! end LAMIPII test

! Do not interpolate in time if data time exactly matches model time

              If(TIME == TIME1) Then
                LINTERPOLATE=.FALSE.
              End If

            End If ! End of REGULAR/not REGULAR

          Else  ! PERIODIC data

!L 2.2   If data is taken from ancillary periodic data.

! DEPENDS ON: time2sec
            CALL TIME2SEC(I_YEAR,I_MONTH,I_DAY,I_HOUR,                  &
     &                     I_MINUTE,I_SECOND,                           &
     &                     ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,       &
     &                     LCAL360)


            If (REGULAR) Then
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
      If (LCAL360) Then
              If(PERIOD /= 8640.AND.PERIOD /= 720.AND.PERIOD /= 24)Then
                ICODE=600+FIELD
           CMESSAGE='REPLANCA: Non-standard period for periodic data'
                RETURN
              End If
      Else
              If(PERIOD /= 24)Then
                ICODE=600+FIELD
           CMESSAGE='REPLANCA: Non-standard period for periodic data'
                RETURN
              End If
      End If
              If(PERIOD == 24)Then
! Ancillary data interval in hour(s), period is 1 day

                IY=I_YEAR
                IM=I_MONTH
                ID=I_DAY
                If(I_HOUR <  FIXHD(24,FILE)) HOURS=HOURS+24

              Else If(PERIOD == 720)Then
! Ancillary data interval in day(s) or hours , period is 1 month

                IY=I_YEAR
                IM=I_MONTH
                ID=FIXHD(23,FILE)
                If((I_DAY*24+I_HOUR) <                                  &
     &             (FIXHD(23,FILE)*24+FIXHD(24,FILE)))                  &
     &           HOURS=HOURS+720

              Else If(PERIOD == 8640)Then
! Ancillary data interval in month(s)or days or hours, period is 1 year

                IY=I_YEAR
                IM=FIXHD(22,FILE)
                ID=FIXHD(23,FILE)
                If((I_MONTH*720+I_DAY*24+I_HOUR) <                      &
     &          (FIXHD(22,FILE)*720+FIXHD(23,FILE)*24+FIXHD(24,FILE)))  &
     &           HOURS=HOURS+8640

              End If

! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,ID,FIXHD(24,FILE),                    &
     &                     FIXHD(25,FILE),FIXHD(26,FILE),               &
     &                     ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,       &
     &                     LCAL360)
              HOURS=HOURS-SEC/3600-DAY*24

! Do not interpolate in time if data time exactly matches model time

              If(MOD(HOURS,INTERVAL) == 0) Then
                LINTERPOLATE=.FALSE.
              End If
              STEP=HOURS/INTERVAL
              TIME=REAL(HOURS)
              TIME1=STEP*INTERVAL
              TIME2=(STEP+1)*INTERVAL

            Else  ! non regular case

!L 2.2.2 Gregorian calender,and data interval is in months,
!L       period is 1 year
!L       Updating interval and number of data times to be skipped
!L       calculated in months.

              TIME=REAL(SEC)/3600+REAL(DAY*24)
              INTERVAL=FIXHD(36,FILE)+FIXHD(35,FILE)*12
              PERIOD=INTHD(3,FILE)*INTERVAL
              If(PERIOD /= 12)Then
                ICODE=600+FIELD
           CMESSAGE='REPLANCA: Non-standard period for periodic data'
                RETURN
              End If
!  Difference between date now (month) & first date ancil file (month)
              MONTHS=I_MONTH-FIXHD(22,FILE)


           If (LAMIPII) THEN ! correct code to use lookup header dates
! Correctly use day and hour from lookup header not fixhd which
! contains values for first field on ancillary file only.
             step=months/INTERVAL
!            I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*step
             I2=NLOOKUP(irec)+LOOKUP_STEP(irec)*step
             I1=I2+LOOKUP_START(FILE)-1
!  Check for time within month - using ppheader information
          If((I_DAY*24+I_HOUR) <  (lookup(3,i1)*24+lookup(4,i1))) Then
                MONTHS=MONTHS-1
          End If
             If(MONTHS <  0) Then
                MONTHS=MONTHS+12
             End If
! recalculate STEP
              STEP=MONTHS/INTERVAL
! NB INTERVAL may be > 1 month
              MONTHS=STEP*INTERVAL
              IY=I_YEAR
              IM=MOD(FIXHD(22,FILE)+MONTHS-1,12)+1
              If(IM >  I_MONTH) IY=IY-1
!             I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*STEP
              I2=NLOOKUP(irec)+LOOKUP_STEP(irec)*STEP
              I1=I2+LOOKUP_START(FILE)-1
! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,lookup(3,i1),lookup(4,i1),            &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,LCAL360)
              TIME1=REAL(SEC)/3600+REAL(DAY*24)
!  Calculate  TIME2 for second ancillary data time
!  set IY correctly for time interpolation calculations
              IY=I_YEAR
              IM=MOD(FIXHD(22,FILE)+MONTHS+INTERVAL-1,12)+1
              If(IM <  I_MONTH) IY=IY+1
              I1=(IM-1)/INTERVAL
!             I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*I1
              I2=NLOOKUP(irec)+LOOKUP_STEP(irec)*I1
              I1=I2+LOOKUP_START(FILE)-1
! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,lookup(3,i1),lookup(4,i1),            &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,LCAL360)
              TIME2=REAL(SEC)/3600+REAL(DAY*24)

           Else   ! original code inaccurate use of FIXHD dates
!  Check for time within month
              If((I_DAY*24+I_HOUR) <                                    &
     &           (FIXHD(23,FILE)*24+FIXHD(24,FILE))) Then
                MONTHS=MONTHS-1
              End If
              If(MONTHS <  0) Then
                MONTHS=MONTHS+12
              End If

              STEP=MONTHS/INTERVAL
! NB INTERVAL may be > 1 month
              MONTHS=STEP*INTERVAL
!  Calculate TIME1 for first ancillary data time
!  set IY correctly for time interpolation calculations
              IY=I_YEAR
              IM=MOD(FIXHD(22,FILE)+MONTHS-1,12)+1
              If(IM >  I_MONTH) IY=IY-1
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
              If(IM <  I_MONTH) IY=IY+1
! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),        &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,              &
     &              LCAL360)
              TIME2=REAL(SEC)/3600+REAL(DAY*24)
           End If  ! end LAMIPII test

! Do not interpolate in time if data time exactly matches model time

              If(TIME == TIME1) Then
                LINTERPOLATE=.FALSE.
              End If

            End If  ! regular/non-regular


          End If  ! non-periodic/periodic


        End If ! singletime/non-singletime

!L 2.3   Check STASH Code

        I2=NLOOKUP(irec)+LOOKUP_STEP(irec)*STEP

        I1=LOOKUP(ITEM_CODE,I2+LOOKUP_START(FILE)-1)

        LMISMATCH=.FALSE.

      If (PrintStatus >= PrStatus_Diag .and. mype == 0 ) Then
        WRITE(6,*)' Information used in checking ancillary data set:',  &
     &  ' position of lookup table in dataset:',I2
        WRITE(6,*)' Position of first lookup table referring to ',      &
     &  'data type ',NLOOKUP(irec)
        WRITE(6,*)' Interval between lookup tables referring to data ', &
     &  'type ', LOOKUP_STEP(irec),' Number of steps', STEP
        WRITE(6,*)' STASH code in dataset ',I1,                         &
     &  '  STASH code requested ',STASHANCIL(irec)
        WRITE(6,*)'''Start'' position of lookup tables for dataset ',   &
     &  'in overall lookup array ' ,LOOKUP_START(FILE)
      End If

        If(I1 /= STASHANCIL(irec)) Then
        WRITE(6,*)I1,STASHANCIL(irec),irec
          LMISMATCH=.TRUE.
        End If

!L Error exit if checks fail

        If(LMISMATCH) Then
          ICODE=200+FIELD
         CMESSAGE='REPLANCA: PP HEADERS ON ANCILLARY FILE DO NOT MATCH'
         RETURN
        End If

        If(LINTERPOLATE.AND..NOT.SINGLE_TIME) Then
!L Check time interpolation factors
          If(TIME <  TIME1.OR.TIME >  TIME2) Then
           WRITE(6,*)' Information used in interpolation/replacement:'
           WRITE(6,*)' Time of first data=', TIME1
           WRITE(6,*)' Validity Time for update=', TIME
           WRITE(6,*)' Time of second data=', TIME2

           ICODE=500+FIELD
           CMESSAGE='REPLANCA: TIME INTERPOLATION ERROR'
           RETURN
          End If
        End If

!L 3   Loop over levels of ancillary data for field I
!L Reset pointer for dataset


!L Includes loop over X and Y components of surface currents

         LICE_FRACTION=FIELD == 27
         LSNOW_DEPTH=FIELD == 9
         LICE_DEPTH=FIELD == 29

        Do 30 LEVEL=1,LEVELS(IREC)

!L Do not go through loop for ice edge or snow edge

        If((LICE_FRACTION.OR.LSNOW_DEPTH).AND.LEVEL == 2) Then
          GOTO 30
        End If

!L 3.1 Read data for single level of ancillary field.

        If(.NOT.LICE_FRACTION) Then
! AMIPII case ice depth field not read from ancillary file
         If(.NOT.(LICE_DEPTH.and.LAMIPII)) Then

!          Check to see if field is one of the River Routing ones.
           If ( FIELD  ==  124 .OR.                                     &
     &          FIELD  ==  125 .OR.                                     &
     &          FIELD  ==  126 ) Then

             If (PrintStatus >= PrStatus_Diag .and. mype == 0 ) Then
               Write (6,*) 'Rcf_ReplancA: Resetting Ancil field size ', &
     &              'to RR_FIELD'
             End If

             FIELD_SIZE=RR_FIELD

           Else ! If not RR, then assume field size is PFIELD

             FIELD_SIZE=P_FIELD

           End If

! DEPENDS ON: rcf_readflds
          Call Rcf_ReadFlds                                             &
     &                 (AncF_UnitNo,1,I2,LOOKUP(1,LOOKUP_START(FILE)),  &
     &                  LEN1_LOOKUP,ANCIL1,FIELD_SIZE,FIXHD(1,FILE),    &
     &                  ICODE,CMESSAGE)

         End If

            If(ICODE >  0)Then
              ICODE=FIELD+100
              CMESSAGE='REPLANCA :I/O ERROR '
              Go To 900
            End If

        Else

!L If ice-fraction,read fractional time field as well
!L       UNLESS IT IS A SINGLE TIME FIELD
!L If snow-depth,read fractional time field as well only if time
!L interpolation required.

      If(.NOT.SINGLE_TIME.and..NOT.LAMIPII) Then
         If(LOOKUP(ITEM_CODE,I2+LOOKUP_START(FILE)) == 38) Then

! DEPENDS ON: rcf_readflds
          Call Rcf_ReadFlds                                             &
     &                 (AncF_UnitNo,2,I2,LOOKUP(1,LOOKUP_START(FILE)),  &
     &                  LEN1_LOOKUP,ICE_EXTENT,P_FIELD,FIXHD(1,FILE),   &
     &                  ICODE,CMESSAGE)

          If(ICODE >  0)Then
            ICODE=FIELD+100
            CMESSAGE='REPLANCA :I/O ERROR '
            Go To 900
          End If

         Else
           ICODE=FIELD+100
           CMESSAGE='REPLANCA :ICE CHANGE DATA MISSING'
           Go To 900
         End If
        Else    ! single time or LAMIPII - ie no time change field

! DEPENDS ON: rcf_readflds
          Call Rcf_ReadFlds                                             &
     &                 (AncF_UnitNo,1,I2,LOOKUP(1,LOOKUP_START(FILE)),  &
     &                  LEN1_LOOKUP,ICE_EXTENT,P_FIELD,FIXHD(1,FILE),   &
     &                  ICODE,CMESSAGE)

          If(ICODE >  0)Then
            ICODE=FIELD+100
            CMESSAGE='REPLANCA :I/O ERROR '
            Go To 900
          End If
         End If
      End If

        If(LSNOW_DEPTH.AND.LINTERPOLATE) Then
      If(LOOKUP(ITEM_CODE,I2+LOOKUP_START(FILE)) == 27) Then

! DEPENDS ON: rcf_readflds
          Call Rcf_ReadFlds                                             &
     &                  (AncF_UnitNo,1,I2+1,                            &
     &                   LOOKUP(1,LOOKUP_START(FILE)),                  &
     &                   LEN1_LOOKUP,SNOW_CHANGE,P_FIELD,FIXHD(1,FILE), &
     &                   ICODE,CMESSAGE)

          If(ICODE >  0)Then
            ICODE=FIELD+100
            CMESSAGE='REPLANCA :I/O ERROR '
            Go To 900
          End If

         Else
           ICODE=FIELD+100
           CMESSAGE='REPLANCA :SNOW CHANGE DATA MISSING'
           Go To 900
         End If
        End If

!L If sea surface temperature or other ice fields, read ice fraction
!L and fractional time field if not already pressent and if required
!L by time interpolation.  

        If(FIELD == 29.OR.(FIELD == 28.AND.LT_INT_C) )                  &
     &    Then

         If(.NOT.UPDATE(ipos_27)) Then
          I3 = NLOOKUP(ipos_27) + LOOKUP_STEP(ipos_27)*STEP             &
     &                          + LOOKUP_START(FILEANCIL(ipos_27))

          If ( LOOKUP(ITEM_CODE,I3)  ==  38 ) Then

! DEPENDS ON: rcf_readflds
          Call Rcf_ReadFlds                                             &
     &                   (AncF_UnitNo,2,                                &
     &                    NLOOKUP(ipos_27)+LOOKUP_STEP(ipos_27)*STEP,   &
     &                    LOOKUP(1,LOOKUP_START(FILEANCIL(ipos_27))),   &
     &                    LEN1_LOOKUP,ICE_EXTENT,                       &
     &                    P_FIELD,FIXHD(1,FILEANCIL(ipos_27)),          &
     &                    ICODE,CMESSAGE)

          If(ICODE /= 0)Then
            ICODE=FIELD+100
            CMESSAGE='REPLANCA :I/O ERROR '
            Go To 900
          End If
          If ( RLOOKUP(BMDI,I3-1)  /=  RMDI ) Then
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: RMDI in lookup of ancillary file of ti&
     &mes of sea-ice chge not standard'
            Go To 900
          End If


          Else
            ICODE=FIELD+100
            CMESSAGE='REPLANCA :ICE FIELD DATA MISSING'
            Go To 900
          End If
         End If
        End If

!L 3.3 If time interpolation required, read second record

        If(LINTERPOLATE) Then

          I1=I2+ LOOKUP_STEP(irec)
          If(I1 <= FIXHD(152,FILE)) Then

! AMIP II and ice depth don't read in ice depth field
          If (.NOT.(LAMIPII.and.LICE_DEPTH)) Then

!             Check assuption that Ancil is PFIELD size is safe.
              If ( FIELD  ==  124 .OR.                                  &
     &             FIELD  ==  125 .OR.                                  &
     &             FIELD  ==  126 ) Then
                Write (6,*) 'Rcf_ReplancA: May be about to use P_FIELD' &
     &                     ,' for a River Routing field'
                Write (6,*) 'Rcf_ReplancA: Field no.',Field
                ICODE=99
                CMessage='Rcf_ReplancA: Error in Ancil field size'
                Go To 900
              End If

! DEPENDS ON: rcf_readflds
          Call Rcf_ReadFlds                                             &
     &                   (AncF_UnitNo,1,I1,                             &
     &                    LOOKUP(1,LOOKUP_START(FILE)),                 &
     &                    LEN1_LOOKUP,ANCIL2,P_FIELD,FIXHD(1,FILE),     &
     &                    ICODE,CMESSAGE)

          End If

          If(ICODE /= 0)Then
            ICODE=FIELD+100
            CMESSAGE='REPLANCA :I/O ERROR '
            Go To 900
          End If

          Else !end of data on file

!L  If end of data has been reached go back to the start.If data is
!L  periodic.
!L  Otherwise cancel time interpolation

            If(PERIODIC) Then

              I1 = NLOOKUP(irec) + LEVEL - 1

!             Check assuption that Ancil is PFIELD size is safe.
              If ( FIELD  ==  124 .OR.                                  &
     &             FIELD  ==  125 .OR.                                  &
     &             FIELD  ==  126 ) Then
                Write (6,*) 'Rcf_ReplancA: May be about to use P_FIELD' &
     &                     ,' for a River Routing field'
                Write (6,*) 'Rcf_ReplancA: Field no.',Field
                ICODE=199
                CMessage='Rcf_ReplancA: Error in Ancil field size'
                Go To 900
              End If

! DEPENDS ON: rcf_readflds
          Call Rcf_ReadFlds                                             &
     &                     (AncF_UnitNo,1,I1,                           &
     &                      LOOKUP(1,LOOKUP_START(FILE)),               &
     &                      LEN1_LOOKUP,ANCIL2,P_FIELD,FIXHD(1,FILE),   &
     &                      ICODE,CMESSAGE)

          If(ICODE /= 0)Then
            ICODE=FIELD+100
            CMESSAGE='REPLANCA :I/O ERROR '
            Go To 900
          End If
            Else
              LINTERPOLATE=.FALSE.
            End If
          End If! End of position on file test

          ICODE=0
        End If ! End LINTERPOLATE

!L 3.4 Perform time interpolation

        If(LINTERPOLATE) Then

          ZERO=0.0

!L Select appropriate time interpolation for each field
!  Snowdepth: set equal to zero if no snow cover

          If(LSNOW_DEPTH) Then
            Do I=1,P_FIELD
              PRES_VALUE(I)=ZERO
            End Do

! For the call to T_INT_C, need to know BMDI is OK for SNOW_CHANGE
!  which was read in from position I2+1.
          If ( RLOOKUP(BMDI,LOOKUP_START(FILE)+I2)  /=  RMDI ) Then
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: RMDI in lookup of ancillary file of ti&
     &mes of snow change non-standard '
            Go To 900
          End If

! DEPENDS ON: t_int_c
            CALL T_INT_C (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,         &
     &           TIME,P_FIELD,SNOW_CHANGE,ANCIL1,PRES_VALUE)

! Ice fraction: ice depth set equal to zero if no ice

          Else If(FIELD == 27.OR.FIELD == 29) Then
            If(FIELD == 27) Then
! For the call to T_INT_C, need to know BMDI is OK for ICE_EXTENT(1,2)
!  which was read in from position I1+1
          If(.NOT.LAMIPII) Then
          If ( RLOOKUP(BMDI,LOOKUP_START(FILE)+I1)  /=  RMDI ) Then
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: RMDI in lookup of ancillary file of ti&
     &mes of sea-ice chge non-standard'
            Go To 900
          End If
          End If

             If (LAMIPII) Then
! linear uncontrolled time interpolation
! DEPENDS ON: t_int
              CALL T_INT (ICE_EXTENT,TIME1,ANCIL2,TIME2,ANCIL_DATA,     &
     &             TIME,P_FIELD)

! For AMIP II strictly ice concentrations should range between
! 0.0 and 1.0 but because of assumptions on T* made by the boundary
! layer and radiation schemes ice concentrations are restricted to
! 0.3 to 1.0. This will allow SSTs in areas of less than 30% ice to
! be used rather than TFS=-1.8C.

              Do I=1,P_FIELD
                If (ANCIL_DATA(I) <  0.3) ANCIL_DATA(I)=0.0
                If (ANCIL_DATA(I) >  1.0) ANCIL_DATA(I)=1.0
              End Do

             Else       ! non AMIPII option
              Do I=1,P_FIELD
                PRES_VALUE(I)=0
              End Do

! DEPENDS ON: t_int_c
              CALL T_INT_C (ICE_EXTENT,TIME1,ANCIL2,TIME2,ANCIL_DATA,   &
     &             TIME,P_FIELD,ICE_EXTENT(1,2),ICE_EXTENT,PRES_VALUE)

             End If     ! end AMIPII test

            Else If (FIELD == 29) Then

              Do I=1,P_FIELD
                PRES_VALUE(I)=0
              End Do

! DEPENDS ON: t_int_c
              CALL T_INT_C (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,       &
     &             TIME,P_FIELD,ICE_EXTENT(1,2),ICE_EXTENT,PRES_VALUE)


            End If


! Sea surface temperature, set equal to TFS if ice present

          Else If (FIELD == 28.AND.LT_INT_C) Then
           If (LAMIPII) Then

! DEPENDS ON: t_int
            CALL T_INT (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,           &
     &              TIME,P_FIELD)
! remove any T below TFS
            Do I=1,P_FIELD
              If (ANCIL_DATA(i) <  TFS)  ANCIL_DATA(I)=TFS
            End Do

           Else     ! non AMIPII option

            Do I=1,P_FIELD
                PRES_VALUE(I)=TFS

! Set no_ice_extent indicator for controlled SST interpolation
                If(ICE_EXTENT(I,1) == 0) Then
                  NO_ICE_EXTENT(I)=1.0
                Else
                  NO_ICE_EXTENT(I)=0.0
                End If
            End Do

! DEPENDS ON: t_int_c
            CALL T_INT_C (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,         &
     &           TIME,P_FIELD,ICE_EXTENT(1,2),NO_ICE_EXTENT,PRES_VALUE)

           End If   ! end AMIPII test
! Otherwise linear interpolation in time, unless missing data indicator
! present at either time.

          Else

! Time interpolation checks the data against the standard missing data
!   indicator - check that the field is labelled as using the same one.
!  (It is to have the right I1 here that I3 is used above.)
          If ( RLOOKUP(BMDI,LOOKUP_START(FILE)+I1-1)  /=  RMDI .OR.     &
     &         RLOOKUP(BMDI,LOOKUP_START(FILE)+I2-1)  /=  RMDI ) Then
            WRITE (6, *) 'LOOKUPS:',                                    &
     &         RLOOKUP(BMDI,LOOKUP_START(FILE)+I1-1),                   &
     &         RLOOKUP(BMDI,LOOKUP_START(FILE)+I2-1)
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: Missing data indicator in lookup of an&
     &cillary file is non-standard    '
            Go To 900
          End If

!         Check assuption that Ancil is PFIELD size is safe.
          If ( FIELD  ==  124 .OR.                                      &
     &         FIELD  ==  125 .OR.                                      &
     &         FIELD  ==  126 ) Then
            Write (6,*) 'Rcf_ReplancA: May be about to use P_FIELD'     &
     &                 ,' for a River Routing field'
            Write (6,*) 'Rcf_ReplancA: Field no.',Field
            ICODE=119
            CMessage='Rcf_ReplancA: Error in Ancil field size'
            Go To 900
          End If

          LEN=P_FIELD
!L  Ozone, test for zonal mean or full field
          If(FIELD == 7) Then
            If(LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) Then
              LEN=P_ROWS
            End If
!   Tropopause-based ozone, test for zonal mean or full field.
!   Currently the same test as for conventional ozone.
          Else If (FIELD == 110) Then
            If(LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) Then
              LEN=P_ROWS
            End If
!   Cariolle ozone, test for zonal mean or full field. 
!   Currently same test as for conventional ozone.
          Else if (FIELD >= 178.AND.FIELD <= 185) THEN
            IF(LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) THEN
              LEN=P_ROWS
            END IF
          End If

! DEPENDS ON: t_int
            CALL T_INT(ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,            &
     &                 TIME,LEN)

          End If ! End Lsnow_depth

! If no interpolation, copy data into final array

        Else ! no interpolation
         If(LICE_FRACTION) Then
          If (LAMIPII) Then
          Do I=1,P_FIELD

          ANCIL_DATA(I)=ICE_EXTENT(I,1)

! For AMIP II strictly ice concentrations should range between
! 0.0 and 1.0 but because of assumptions on T* made by the boundary
! layer and radiation schemes ice concentrations are restricted to
! 0.3 to 1.0. This will allow SSTs in areas of less than 30% ice to
! be used rather than TFS=-1.8C.

             If (ANCIL_DATA(I) <  0.3) ANCIL_DATA(I)=0.0
             If (ANCIL_DATA(I) >  1.0) ANCIL_DATA(I)=1.0

          End Do
          Else           ! non AMIP II option
            Do I=1,P_FIELD
             ANCIL_DATA(I)=ICE_EXTENT(I,1)
            End Do
          End If           ! end of AMIPII test
         Else If (LAMIPII.AND.FIELD == 28) Then
          Do I=1,P_FIELD
            ANCIL_DATA(I)=ANCIL1(I)
            If (ANCIL_DATA(I) <  TFS) ANCIL_DATA(I)=TFS
          End Do
         Else
!          Check to see if field is one of the River Routing ones.
           If ( FIELD  ==  124 .OR.                                     &
     &          FIELD  ==  125 .OR.                                     &
     &          FIELD  ==  126 ) Then

             If (PrintStatus >= PrStatus_Diag .and. mype == 0 ) Then
               Write (6,*) 'Rcf_ReplancA: Resetting Ancil field size ', &
     &              'to RR_FIELD'
             End If

             FIELD_SIZE=RR_FIELD

           Else ! If not RR, then assume field size is PFIELD

             FIELD_SIZE=P_FIELD

           End If

           Do I=1,FIELD_SIZE
             ANCIL_DATA(I)=ANCIL1(I)
           End Do

         End If
        End If !End interpolate/no interpolate

!L 3.5 Updating action for each field at each level
!L     Fields replaced except that Sea Surface Temperature may be
!L     incremented. Take appropriate action for each field.

        If(FIELD <= 2.OR.FIELD == 7.OR.FIELD == 39.OR.FIELD == 40       &
     & .OR.FIELD == 41.OR.FIELD == 42.OR.FIELD == 43                    &
     & .OR.FIELD == 44.OR.FIELD == 45                                   &
                                           ! multi-level murk
     & .OR.(FIELD >= 48 .AND. FIELD <= 67 .AND. L_UKCA)                 &
                                          ! single-level user ancillaries
     & .OR.(FIELD >= 68.AND.FIELD <= 70)                                &
                                           !NH3,soot aerosol emissions
     & .OR.(FIELD >= 72.AND.FIELD <= 77)                                &
                                           !Sulphur cycle
     & .OR.FIELD == 78                                                  &
                                           !CO2 EMISSIONS
     & .OR.FIELD == 82                                                  &
                                           !HADCM2 sulphate aerosol
     & .OR.(FIELD >= 90.AND.FIELD <= 109)                               &
                                           !multi-level user ancillaries
     & .OR.(FIELD >= 112.AND.FIELD <= 120)                              &
                                           !mineral dust fields
     & .OR.(FIELD >= 121.AND.FIELD <= 122)                              &
                                           !Biomass emissions
     & .OR.FIELD == 123                                                 &
                                           !Seawater DMS concentration
     & .OR.(FIELD >= 157.AND.FIELD <= 177)                              &
                                           !Aerosol climatologies
     & .OR.(FIELD >= 178.AND.FIELD <= 185)                              &
                                           !Cariolle ozone ancillaries
     & .OR.(FIELD >= 186.AND.FIELD <= 187)                              &
                                           !OCFF emissions
     &  )THEN

!L 3.5.0 Updates at all points

          LEN=P_FIELD
!L  Ozone, test for zonal mean or full field
          If(FIELD == 7) Then
            If(LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) Then
              LEN=P_ROWS
            End If
!   Tropopause-based ozone, test for zonal mean or full field.
!   Currently the same test as for conventional ozone.
          Else If (FIELD == 110) Then
            If(LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) Then
              LEN=P_ROWS
            END IF
!   Cariolle ozone, test for zonal mean or full field. 
!   Currently same test as for conventional ozone.
          ELSE IF (FIELD >= 178.AND.FIELD <= 185) THEN
            IF(LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) THEN
              LEN=P_ROWS
            END IF
          END IF

          Do I=1,LEN
            D1(ANCIL_ADD(irec)+I-1+(LEVEL-1)*LEN)=ANCIL_DATA(I)
          End Do

!L 3.5.1 Updates over all land points

        Else If((FIELD > 2.AND.FIELD < 7)                               &
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
     &   .OR.(FIELD >= 83.AND.FIELD <= 89)                              &
                                               !MOSES-II
     &   .OR.(FIELD == 111)) THEN                !COASTAL TILING

!  Set default value of Z0 over sea
          If (FIELD == 26) Then
            Do I=1,P_FIELD
              If(SEA(I)) Then
                D1(ANCIL_ADD(irec)+I-1+(LEVEL-1)*P_FIELD)=10.E-4
              End If
            End Do
          End If

          Do I=1,P_FIELD
            If (LAND(I)) Then
              D1(ANCIL_ADD(irec)+I-1+(LEVEL-1)*P_FIELD)=ANCIL_DATA(I)
            End If
          End Do

!L  Reset TSTAR to TM if snow cover present

          If(LSNOW_DEPTH) Then
            Do I=1,P_FIELD
              If (LAND(I).AND.ANCIL_DATA(I) >  0.0) Then
                If(TSTAR_LAND(I) >  TM) TSTAR_LAND(I)=TM
              End If
            End Do
          End If

!L 3.5.2 Ice fraction

        Else If(FIELD == 27) Then

          Do I=1,P_FIELD
            ICE_FRACTION(I)=0.
            If(SEA(I)) Then
              ICE_FRACTION(I)=ANCIL_DATA(I)
            End If
          End Do

!L Reduce TSTAR to TFS where ice fraction greater than zero
! Required at present because radiation and boundary layer codes
! assume T* is TFS and ignore any value set in TSTAR.

          If(.NOT.LTLEADS)Then
            Do I=1,P_FIELD
              If(ICE_FRACTION(I) >  0.0) Then
                TSTAR_SSI(I)=AMIN1(TSTAR_SSI(I),TFS)
              End If
            End Do
          End If

!L 3.5.3 Sea surface temperatures for atmosphere, allow fields to be
!L       incremented rather than replaced

        Else If (FIELD == 28) Then

          If(L_SSTANOM) Then
            Do I=1,P_FIELD
              TSTAR_ANOM(I)=0.0
            End Do
          End If

          Do I=1,P_FIELD
            If(L_CTILE.OR.ICE_FRACTION(I) == 0.0)Then
              If(L_SSTANOM) Then
                TSTAR_ANOM(I)=TSTAR_SEA(I)-ANCIL_DATA(I)
              Else
                TSTAR_SEA(I)=ANCIL_DATA(I)
              End If
                If (ICE_FRACTION(I) == 0.0) TSTAR_SSI(I)=TSTAR_SEA(I)
            End If
          End Do

!L 3.5.4 Sea ice thickness
!L       Update over all sea points (all sea ice points are the only
!L       ones strictly required, but this cannot be determined easily)

        Else If (FIELD == 29) Then

          Do I=1,P_FIELD
            If(SEA(I)) Then
              D1(ANCIL_ADD(irec)+I-1)=ANCIL_DATA(I)
            End If
          End Do

!L 3.5.5 Surface currents

        Else If (FIELD == 30) Then
          Do I=1,U_FIELD
            D1(ANCIL_ADD(irec)+I-1)=ANCIL_DATA(I)
          End Do

        Else If (FIELD == 31) Then
          Do I=1,V_FIELD
            D1(ANCIL_ADD(irec)+I-1)=ANCIL_DATA(I)
          End Do

!L 3.5.7 River Routing Ancils
!L       Update over all points

        Else If (FIELD == 124 .OR.                                      &
     &           FIELD == 125 .OR.                                      &
     &           FIELD == 126 ) Then

          Do I=1,RR_FIELD
            D1(ANCIL_ADD(irec)+I-1)=ANCIL_DATA(I)
          End Do

!Tracer Fluxes - kdcorbin, 04/10
        Else If (FIELD .ge. 188 .and. FIELD .lt. 208) Then
          Do I=1,P_FIELD
             D1(ANCIL_ADD(irec)+I-1) = ANCIL_DATA(I)
          End Do

        Else

        WRITE(6,*)' REPLANCA: ERROR - FIELD ',FIELD,                    &
     &  ' omitted from update block'

        End If !End tests on FIELD numbers

!L End loop over levels

      I2=I2+1

 30   CONTINUE

!L End loop over ancillary fields (atmosphere)
       End If ! LAMIPII and ice depth test

      End If    ! End UPDATE(irec) test     level 1 IF

      End Do

      If(L_CTILE)Then
        Do I=1,P_FIELD
          If(SEA(I).AND.ICE_FRACTION(I) >  0.0)Then

            If(LTLEADS.OR.LAMIPII)Then
              TSTAR_SICE(I)=AMIN1(TSTAR_SICE(I),TFS)
              TSTAR_SSI(I)=ICE_FRACTION(I)*TSTAR_SICE(I)                &
     &          +(1.-ICE_FRACTION(I))*TSTAR_SEA(I)
            Else
              TSTAR_SEA(I)=TFS
              TSTAR_SICE(I)=(TSTAR_SSI(I)                               &
     &          -(1.-ICE_FRACTION(I))*TSTAR_SEA(I))/ICE_FRACTION(I)
            End If

          End If
!
          TSTAR(I)=FLANDG(I)*TSTAR_LAND(I)                              &
     &      +(1.-FLANDG(I))*TSTAR_SSI(I)
        End Do
      Else
        Do I=1,P_FIELD
          If(LAND(I))Then
            TSTAR(I)=TSTAR_LAND(I)
          Else
            TSTAR(I)=TSTAR_SSI(I)
          End If
        End Do
      End If

!     Set up surface temperatures:
      If(L_CTILE)Then
        Do I=1,P_FIELD
          TSTAR_LAND_CTILE(I)=TSTAR_LAND(I)
          TSTAR_SEA_CTILE(I)=TSTAR_SEA(I)
          ! Ensure consistency with equivalent code in 
          ! replanca-rpanca1a.F90. Also helps to avoid 
          ! crazy values of TSTAR_SICE in reconfigured dumps.
          ! TSTAR_SICE_CTILE(I)=TSTAR_SICE(I)
        End Do
      End If

!     Deallocate the temporary storage arrays used for ancil data.
      DeAllocate ( ANCIL1 )
      DeAllocate ( ANCIL2 )
      DeAllocate ( ANCIL_DATA )

900   RETURN
      END SUBROUTINE REPLANCA
