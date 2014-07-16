
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Set up lists and space to store Header info for data created not read

Module Rcf_Grib_Spcl_Hdr_Mod

! SUBROUTINE Rcf_Grib_Spcl_Hdr
!
! Description: Routind which creates lists and list members to hold
!              the header info for fields which are created rather than
!              read from the original GRIB data.
!              -NB this routine should not be used to alter the headers
!               of fields which are altered (e.g. geopotential to orog)
!
! Method:
!
! Current Code Owner: Roddy Sharp
!
! History:
! Version   Date      Comment
! -------  --------   -------------------------
!  5.4     21/08/02   Original code. Roddy Sharp (frtz)
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
contains
Subroutine Rcf_Grib_Spcl_Hdr(Lists)

Use Rcf_GRIB_Block_Params_Mod, Only :      &
  List_Marker,          Grib_Record,       &
  p_Lvl_Type,           p_Lvl_Desc_1,      &
  Tb3_Pressure,         p_Orig_cntr,       &
  p_Param_ID,           EID_Temperature,   &
  Grb_Data_Real

Use Rcf_GRIB_Lookups_Mod, Only :    &
  grib_max_fields,           grib_Exner_field,           &
  GrbOrigECMWF

Use EReport_Mod, Only :     &
    EReport

Use Rcf_Parvars_mod, Only : &
    mype

Use Rcf_StashCodes_Mod, Only : &
   Stashcode_exner

Implicit None

! Global variables (#include statements etc):

! Subroutine arguments

!< Array  arguments with intent(InOut):>
Type (List_Marker), Intent(InOut) :: Lists(0:grib_max_fields)

! Comdecks
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
! contains LBLREC (amongst others)

! Local constants
Character (Len=*), Parameter     :: RoutineName='Rcf_GRIB_Spcl_Hdr'

! Local variables

Type (Grib_Record),Pointer       :: New_Record,Current

Character (Len=80)               :: Cmessage(2)   ! used for EReport
Integer                          :: ErrorStatus   ! used for EReport

Integer                          :: I
Integer                          :: Count,Criteria

!=======================================================================
!  Routine Code Start :
!=======================================================================

!=======================================================================
!  Loop across all lists, finding those which will require extra fields
!=======================================================================
Do I = 1, grib_max_fields

  If (Associated(Lists(I) % Begin) ) Then

  ! Lets only look at those ECMWF ones first
  If ( Lists(I) % Begin % Block_1 (p_Orig_cntr) == GrbOrigECMWF ) Then
    Select Case ( Lists(I) % Begin % Block_1 ( p_Param_ID ) )

      Case (EID_Temperature)
        ! Will need Exner for Height generation

        ! Set pointer to first entry in list
        Current => Lists(I) % Begin
        count =0

        ! Loop across members in list
        Do While (Associated(Current))
        count= count+1

        !Allocate a GRIB header to store info
        Write (6,*) 'About to allocate for Exner ', count
        Allocate(New_Record)

        New_Record % Block_1(:) = Current % Block_1(:)
        New_Record % Block_2(:) = Current % Block_2(:)
        New_Record % Block_3(:) = Current % Block_3(:)
        New_Record % Block_4(:) = Current % Block_4(:)
        New_Record % Block_R(:) = Current % Block_R(:)

        New_Record % Start_pos  = 0
        New_Record % Data_Type  = Grb_Data_Real
        New_Record % Num_Fp     = 0
        New_Record % Num_Vert   = 0
        New_Record % Num_Bitmap = 0
        New_Record % Num_Quasi  = 0

        ! Copy header info from Temp to Exner
        New_Record % StashCode = Stashcode_exner
        New_Record % Desc      = 'Exner field'

        New_Record % Block_1 (p_Orig_cntr) = (-1)
        New_Record % Block_1 (p_Param_ID)  = (1)

        ! Assign that record to it's list
        New_Record % Prev   => Lists(grib_Exner_field) % End
                               ! Point Prev pointer at end of
                               ! current list.(Null if first entry)
        If (Associated(Lists(grib_Exner_field) % End)) Then
                               ! If current end of list is a
                               ! valid record (Not first entry)
          Lists(grib_Exner_field) % End % Next  => New_Record
                               ! Point 'next' for previous entry
                               ! at current entry

        Else                         ! Else : must be 1st entry
          Lists(grib_Exner_field) % Begin  => New_Record
                               ! Point begining of List at New_Record
        End If

        Lists(grib_Exner_field) % End      => New_Record
                               ! Point End of List at (now complete)
                               ! New_Record Entry
        Nullify(New_Record % Next)   ! Ensure 'Next' is not associated

        Lists(grib_Exner_field) % LstCount =                          &
                                Lists(grib_Exner_field) % LstCount + 1
                               ! Add one to count of list size


      ! end do across list members
        Current => Current % Next
        End Do

    End Select

  End If
  End If

!=======================================================================
!  End Loop across all lists.
!=======================================================================
End Do  ! loop over all lists


Return

End Subroutine Rcf_Grib_Spcl_Hdr
End Module Rcf_Grib_Spcl_Hdr_Mod
